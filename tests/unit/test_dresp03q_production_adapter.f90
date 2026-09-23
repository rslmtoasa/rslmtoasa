!------------------------------------------------------------------------------
! DRESP-03Q -- production bcc-Fe adapter gate.
!
! The state is built through the ordinary lattice -> charge -> Hamiltonian ->
! reciprocal stack.  The finite-q fixture is then extracted from that live
! Hamiltonian and compared with the production reciprocal assembler at several
! k points.  The source Hamiltonian is intent(in) at the adapter boundary.
!------------------------------------------------------------------------------
program test_dresp03q_production_adapter
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use basis_mod, only: basis_init
   use control_mod, only: control
   use charge_mod, only: charge
   use hamiltonian_mod, only: hamiltonian
   use energy_mod, only: energy
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   use lr_kl_hessian_mod, only: lmto_live_hamiltonian_fixture, lmto_fixture_init, lmto_fixture_from_hamiltonian, &
      lmto_fixture_adapter_residual, assemble_lmto_hamiltonian, assemble_lmto_finite_q_torque, &
      assemble_lmto_finite_q_mixed_derivative, force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch
   use lr_rotation_response_mod, only: rotation_state, rotation_request, rotation_result, &
      prepare_rotation_response, evaluate_rotation_response, evaluate_rotation_response_oracle, &
      reduce_static_rotation_kernel, rotation_circular_unitary, rotation_axes
   use lr_kl_hessian_mod, only: finite_temperature_occupation
   use lr_lmto_turek_contour_mod, only: native_turek_contour_options, native_turek_contour_report, &
      native_turek_static_reference
   use math_mod, only: ang2au, i_unit, init_math_operators
   use precision_mod, only: rp
   use reciprocal_mod, only: reciprocal
   use timer_mod, only: g_timer, timer
   implicit none

   type(control) :: ctl
   type(lattice), target :: lat
   type(charge), target :: chg
   type(hamiltonian), target :: ham
   type(energy) :: ene
   type(reciprocal), target :: recip
   type(native_turek_contour_options) :: native_options
   type(native_turek_contour_report) :: native_report
   type(lmto_live_hamiltonian_fixture), target :: fixture
   complex(rp), allocatable :: h_production(:, :), h_fixture(:, :)
   real(rp), allocatable :: adapter_points(:, :), eigenvalues(:, :), endpoint_values(:, :), weights(:), shifted_points(:, :)
   complex(rp), allocatable :: eigenvectors(:, :, :), endpoint_vectors(:, :, :), torque_q(:, :, :, :), &
      torque_minus_q(:, :, :, :), mixed(:, :, :, :, :), jq_ud(:, :, :), jq_du(:, :, :), jq_sym(:, :, :), &
      native_delta(:, :, :), native_curvature(:, :, :)
   complex(rp) :: hessian(1,1), torque_torque(1,1), mixed_contact(1,1), complete(1,1)
   real(rp), parameter :: q_static(3,2)=reshape([0.0_rp,0.0_rp,0.0_rp, 0.125_rp,0.0_rp,0.0_rp],[3,2])
   real(rp), parameter :: rotation_axis(3)=[1.0_rp,0.0_rp,0.0_rp]
   real(rp) :: gamma_error, generic_error, fermi, kT, q0_hessian, finite_q_curvature, native_difference
   real(rp) :: max_error, max_relative, adapter_error, point_error, relative_error
   integer :: ik, i, iq, nk, nmat, gamma_index, generic_index

   call init_math_operators()
   call g_logger%init()
   g_timer = timer()

   ctl = control('input.nml')
   lat = lattice(ctl)
   call lat%build_data()
   call lat%bravais()
   call lat%structb(.true.)
   call lat%atomlist()
   call basis_init(lat%symbolic_atoms(1)%potential%lmax)

   chg = charge(lat)
   ham = hamiltonian(chg)
   do i = 1, lat%ntype
      call lat%symbolic_atoms(i)%build_pot()
   end do
   if (ctl%nsp == 2 .or. ctl%nsp == 4) call ham%build_lsham()
   call ham%build_bulkham()
   do i=1,lat%ntype
      call lat%symbolic_atoms(i)%predls(lat%wav*ang2au)
   end do

   if (.not. ham%hoh .or. ham%ccor_2c .or. ham%hubbard_u_general_check .or. ham%hubbard_v_check .or. &
       .not. allocated(ham%eeo) .or. .not. allocated(ham%enim)) then
      error stop 'DRESP-03Q production adapter gate requires clean active second-order HOH state'
   end if
   ene=energy(lat)
   recip = reciprocal(ham)
   recip%fermi_level=ene%fermi
   recip%reciprocal_mode = 'ham_only'
   if (trim(recip%kspace_ham_order) /= 'second') error stop 'DRESP-03Q reciprocal order is not second'
   call recip%generate_mp_mesh()
   call recip%build_kspace_hamiltonian()
   if (.not. allocated(recip%hk_bulk)) error stop 'DRESP-03Q accepted reciprocal H(k) cache is missing'
   call lmto_fixture_from_hamiltonian(ham, fixture)
   nmat = 2*fixture%norb*fixture%nsite
   allocate(h_production(nmat,nmat), h_fixture(nmat,nmat))
   gamma_index=0; generic_index=0; max_error=0.0_rp; max_relative=0.0_rp
   do ik=1,size(recip%k_points,2)
      if (maxval(abs(recip%k_points(:,ik)))<1.0e-12_rp) gamma_index=ik
      if (minval(abs(recip%k_points(:,ik)))>1.0e-3_rp) then
         if (maxval(abs(abs(recip%k_points(:,ik))-abs(recip%k_points(1,ik))))>1.0e-2_rp) generic_index=ik
      end if
      h_production=recip%hk_bulk(:,:,ik)
      call assemble_lmto_hamiltonian(fixture,recip%k_points(:,ik),h_fixture)
      point_error=maxval(abs(h_production-h_fixture))
      relative_error=sqrt(sum(abs(h_production-h_fixture)**2))/max(sqrt(sum(abs(h_production)**2)),tiny(1.0_rp))
      max_error=max(max_error,point_error); max_relative=max(max_relative,relative_error)
   end do
   if (gamma_index==0 .or. generic_index==0) error stop 'DRESP-03Q test mesh lacks Gamma or generic k point'
   allocate(adapter_points(3,2))
   adapter_points(:,1)=recip%k_points(:,gamma_index); adapter_points(:,2)=recip%k_points(:,generic_index)
   call lmto_fixture_adapter_residual(ham,fixture,adapter_points,adapter_error)
   call recip%build_hamiltonian_at_kpoint(adapter_points(:,1),h_production)
   call assemble_lmto_hamiltonian(fixture,adapter_points(:,1),h_fixture)
   gamma_error=maxval(abs(h_production-h_fixture))
   call recip%build_hamiltonian_at_kpoint(adapter_points(:,2),h_production)
   call assemble_lmto_hamiltonian(fixture,adapter_points(:,2),h_fixture)
   generic_error=maxval(abs(h_production-h_fixture))
   write (*, '(a,a,5(a,l1))') 'DRESP-03Q order/hoh/eeo/enim/ccor/soc = ', &
      trim(recip%kspace_ham_order), '/', ham%hoh, '/', allocated(ham%eeo), '/', allocated(ham%enim), '/', &
      ham%ccor_2c, '/', (ctl%nsp == 2 .or. ctl%nsp == 4)
   write (*, '(a,es12.4)') 'DRESP-03Q Gamma accepted-H absolute residual = ',gamma_error
   write (*, '(a,es12.4)') 'DRESP-03Q generic-k accepted-H absolute residual = ',generic_error
   write (*, '(a,es12.4)') 'DRESP-03Q reciprocal max absolute residual = ',max_error
   write (*, '(a,es12.4)') 'DRESP-03Q reciprocal max relative Frobenius residual = ',max_relative
   write (*, '(a,es12.4)') 'DRESP-03Q production real-space adapter residual = ', adapter_error
   if (adapter_error > 3.0e-12_rp .or. max_error > 3.0e-12_rp .or. max_relative>3.0e-12_rp .or. &
       gamma_error>3.0e-12_rp .or. generic_error>3.0e-12_rp) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'
   call check_material_rotation_adapter(fixture,q_static,rotation_axis)

   ! Both static references use the same reciprocal mesh, weights, EF and
   ! temperature. The native service evaluates P-S directly; it does not
   ! consume the finite-H eigenpairs or torque/contact matrices.
   nk=size(recip%k_points,2); fermi=recip%fermi_level
   kT=max(recip%temperature*6.3336814e-6_rp,1.0e-10_rp)
   write (*,'(a,i0,a,es14.6,a,es14.6)') 'DRESP-03Q second-order fixture metadata: nk=',nk,' EF=',fermi,' kT=',kT
   allocate(jq_ud(1,1,2),jq_du(1,1,2),jq_sym(1,1,2),native_delta(1,1,2),native_curvature(1,1,2))
   call native_turek_static_reference(lat,recip%k_points,recip%k_weights,q_static,fermi,kT,native_options, &
      jq_ud,jq_du,jq_sym,native_delta,native_curvature,native_report)
   write (*,'(a,i0,a,i0,a,l1)') 'native Turek points/poles/bounds = ',native_report%contour_points,'/', &
      native_report%native_spectral_poles,'/',native_report%native_bounds_verified
   write (*,'(a,2(es14.6,1x))') 'native Turek J_sym(q0), J_sym(qfinite) = ',real(jq_sym(1,1,1),rp),real(jq_sym(1,1,2),rp)
   write (*,'(a,2(es14.6,1x))') 'native Turek curvature(q0), curvature(qfinite) = ', &
      real(native_curvature(1,1,1),rp),real(native_curvature(1,1,2),rp)

   allocate(weights(nk),eigenvalues(nmat,nk),eigenvectors(nmat,nmat,nk),endpoint_values(nmat,nk), &
      endpoint_vectors(nmat,nmat,nk),shifted_points(3,nk),torque_q(nmat,nmat,1,nk), &
      torque_minus_q(nmat,nmat,1,nk),mixed(nmat,nmat,1,1,nk))
   weights=recip%k_weights
   call recip%calculate_eigenpairs_at_kpoints(recip%k_points,eigenvalues,eigenvectors)
   do iq=1,2
      if(iq==1) then
         endpoint_values=eigenvalues; endpoint_vectors=eigenvectors
      else
         do ik=1,nk
            shifted_points(:,ik)=recip%k_points(:,ik)+q_static(:,iq)
         end do
         call recip%calculate_eigenpairs_at_kpoints(shifted_points,endpoint_values,endpoint_vectors)
      end if
      do ik=1,nk
         call assemble_lmto_finite_q_torque(fixture,recip%k_points(:,ik),q_static(:,iq),1,rotation_axis,torque_q(:,:,1,ik))
         call assemble_lmto_finite_q_torque(fixture,recip%k_points(:,ik)+q_static(:,iq),-q_static(:,iq), &
            1,rotation_axis,torque_minus_q(:,:,1,ik))
         call assemble_lmto_finite_q_mixed_derivative(fixture,recip%k_points(:,ik),q_static(:,iq), &
            1,rotation_axis,1,rotation_axis,mixed(:,:,1,1,ik))
      end do
      call force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch(eigenvalues,eigenvectors,endpoint_values, &
         endpoint_vectors,fermi,kT,weights,torque_q,torque_minus_q,mixed,hessian,torque_torque,mixed_contact,complete)
      if(iq==1) then
         q0_hessian=abs(real(complete(1,1),rp))
         write (*,'(a,3(es14.6,1x))') 'second-order static Hessian q0 TT/contact/total = ', &
            real(torque_torque(1,1),rp),real(mixed_contact(1,1),rp),real(complete(1,1),rp)
         if(q0_hessian>2.0e-7_rp) error stop 'Q0_ROTATION_HESSIAN_OPEN'
      else
         finite_q_curvature=real(complete(1,1),rp)
         native_difference=finite_q_curvature-real(native_curvature(1,1,iq),rp)
         write (*,'(a,3(es14.6,1x))') 'second-order finite-q Hessian/native/difference = ', &
            finite_q_curvature,real(native_curvature(1,1,iq),rp),native_difference
         if(.not.ieee_is_finite(finite_q_curvature)) error stop 'second-order finite-q Hessian is non-finite'
      end if
   end do
   call check_rotation_response(fixture,recip)
   call fixture%clear()
   deallocate(h_production,h_fixture,adapter_points,weights,eigenvalues,eigenvectors,endpoint_values,endpoint_vectors, &
      shifted_points,torque_q,torque_minus_q,mixed,jq_ud,jq_du,jq_sym,native_delta,native_curvature)

contains

   subroutine check_rotation_response(base,recip)
      type(lmto_live_hamiltonian_fixture), target, intent(in) :: base
      type(reciprocal), target, intent(inout) :: recip
      type(rotation_state), target :: state0,stateq,statem
      type(rotation_request) :: request
      type(rotation_result) :: result0,resultq,resultm,plus,minus
      real(rp), parameter :: qset(3,2)=reshape([0.0_rp,0.0_rp,0.0_rp,0.125_rp,0.0_rp,0.0_rp],[3,2])
      real(rp), parameter :: eta_set(2)=[1.0e-4_rp,5.0e-5_rp]
      complex(rp), allocatable :: reduced(:, :),unitary(:, :),kernel_pm(:, :),bubble_o(:, :),contact_o(:, :),kernel_o(:, :)
      real(rp), allocatable :: values(:, :),endpoint_values(:, :),points(:, :),weights(:),axes(:, :)
      complex(rp), allocatable :: vectors(:, :, :),endpoint_vectors(:, :, :),torques_q(:, :, :, :), &
         torques_minus_q(:, :, :, :),mixed(:, :, :, :, :)
      complex(rp), allocatable :: rho(:, :),gx(:, :),gy(:, :),commutator(:, :)
      complex(rp) :: trace_value
      complex(rp), allocatable :: hessian(:, :),torque_torque(:, :),mixed_contact(:, :),complete(:, :)
      real(rp) :: q(3),q_cov(3),omega,eta,fermi,kT,weight_sum,weight,fn,m_band,berry_direct,slope_plus,slope_minus
      real(rp) :: static_residual,q0_residual,oracle_residual,covariance_residual,circular_residual,berry_residual,delta_omega
      real(rp) :: max_component,grid_residual,grid_best,grid_point(3)
      integer :: nk,nmat,ncoord,nband,norb,nsite,ik,ia,ib,n,m,site_a,site_b,axis_a,axis_b,up,dn,io,iq,iw,ie,jk

      nsite=base%nsite; norb=base%norb; nmat=2*norb*nsite; ncoord=2*nsite
      nk=size(recip%k_points,2); nband=nmat; fermi=recip%fermi_level
      kT=max(recip%temperature*6.3336814e-6_rp,1.0e-10_rp); weight_sum=sum(recip%k_weights)
      call prepare_rotation_response(base,recip,qset(:,1),state0)
      call prepare_rotation_response(base,recip,qset(:,2),stateq)

      allocate(reduced(ncoord,ncoord),unitary(ncoord,ncoord),kernel_pm(ncoord,ncoord))
      request%state=>state0; request%exact_static=.true.; request%omega=0.0_rp; request%eta=0.0_rp; request%want_inverse=.false.
      call evaluate_rotation_response(request,result0)
      call reduce_static_rotation_kernel(result0%kernel,reduced)
      q0_residual=maxval(abs(reduced))
      write(*,'(a,es12.4)') 'Rotation exact q0 Goldstone max residual = ',q0_residual
      if(q0_residual>2.0e-7_rp) error stop 'Q0_GOLDSTONE_OPEN'
      call rotation_circular_unitary(nsite,unitary)
      kernel_pm=matmul(unitary,matmul(reduced,conjg(transpose(unitary))))
      circular_residual=max(abs(kernel_pm(1,2)),abs(kernel_pm(2,1)))

      call independent_berry_oracle(base,recip,m_band,berry_direct)
      berry_residual=abs(berry_direct-result0%berry)/max(abs(berry_direct),1.0e-12_rp)
      write(*,'(a,4(es14.6,1x))') 'Rotation independent M_band / Berry / residuals = ',m_band,berry_direct, &
         berry_residual,state0%berry_residual
      if(berry_residual>2.0e-10_rp) error stop 'BERRY_COMMUTATOR_ORACLE_OPEN'
      delta_omega=1.0e-5_rp
      request%exact_static=.false.; request%eta=1.0e-9_rp; request%omega=delta_omega
      call evaluate_rotation_response(request,plus)
      request%omega=-delta_omega
      call evaluate_rotation_response(request,minus)
      slope_plus=real((plus%kernel_pm(1,1)-minus%kernel_pm(1,1))/(2.0_rp*delta_omega),rp)
      slope_minus=real((plus%kernel_pm(2,2)-minus%kernel_pm(2,2))/(2.0_rp*delta_omega),rp)
      berry_residual=max(abs(abs(slope_plus)-abs(berry_direct)),abs(abs(slope_minus)-abs(berry_direct)), &
         abs(slope_plus+slope_minus))/max(abs(berry_direct),1.0e-12_rp)
      circular_residual=max(circular_residual,abs(plus%kernel_pm(1,2)),abs(plus%kernel_pm(2,1)))
      write(*,'(a,3(es14.6,1x))') 'Rotation Berry slopes (+/-) / circular offdiag = ',slope_plus,slope_minus,circular_residual
      if(berry_residual>5.0e-2_rp) error stop 'BERRY_NORMALIZATION_OPEN'
      if(circular_residual>1.0e-7_rp) error stop 'CIRCULAR_CONVENTION_OPEN'

      ! Exact static reduction against the independent certified force-theorem
      ! Hessian at Gamma and one non-self-inverse finite q.
      do iq=1,2
         q=qset(:,iq)
         allocate(values(nband,nk),endpoint_values(nband,nk),vectors(nmat,nband,nk), &
            endpoint_vectors(nmat,nband,nk),weights(nk),torques_q(nmat,nmat,ncoord,nk), &
            torques_minus_q(nmat,nmat,ncoord,nk),mixed(nmat,nmat,ncoord,ncoord,nk), &
            hessian(ncoord,ncoord),torque_torque(ncoord,ncoord),mixed_contact(ncoord,ncoord),complete(ncoord,ncoord))
         values=recip%eigenvalues; vectors=recip%eigenvectors; weights=recip%k_weights
         if(maxval(abs(q))<1.0e-12_rp) then
            endpoint_values=values; endpoint_vectors=vectors
         else
            points=recip%k_points+spread(q,2,nk)
            call recip%calculate_eigenpairs_at_kpoints(points,endpoint_values,endpoint_vectors)
            deallocate(points)
         end if
         call rotation_axes(base%moments,axes)
         do ik=1,nk
            do ia=1,ncoord
               site_a=(ia+1)/2; axis_a=1+mod(ia-1,2)
               call assemble_lmto_finite_q_torque(base,recip%k_points(:,ik),q,site_a, &
                  axes(:,2*(site_a-1)+axis_a),torques_q(:,:,ia,ik))
               call assemble_lmto_finite_q_torque(base,recip%k_points(:,ik)+q,-q,site_a, &
                  axes(:,2*(site_a-1)+axis_a),torques_minus_q(:,:,ia,ik))
            end do
            do ia=1,ncoord
               site_a=(ia+1)/2; axis_a=1+mod(ia-1,2)
               do ib=1,ncoord
                  site_b=(ib+1)/2; axis_b=1+mod(ib-1,2)
                  call assemble_lmto_finite_q_mixed_derivative(base,recip%k_points(:,ik),q,site_a, &
                     axes(:,2*(site_a-1)+axis_a),site_b,axes(:,2*(site_b-1)+axis_b),mixed(:,:,ia,ib,ik))
               end do
            end do
         end do
         call force_theorem_finite_q_hessian_from_eigenbasis_metallic_batch(values,vectors,endpoint_values,endpoint_vectors, &
            fermi,kT,weights,torques_q,torques_minus_q,mixed,hessian,torque_torque,mixed_contact,complete)
         if(iq==1) then
            request%state=>state0; request%exact_static=.true.; request%omega=0.0_rp; request%eta=0.0_rp
            call evaluate_rotation_response(request,result0)
            call reduce_static_rotation_kernel(result0%kernel,reduced)
         else
            request%state=>stateq; request%exact_static=.true.; request%omega=0.0_rp; request%eta=0.0_rp
            call evaluate_rotation_response(request,resultq)
            call prepare_rotation_response(base,recip,-q,statem)
            request%state=>statem
            call evaluate_rotation_response(request,resultm)
            call reduce_static_rotation_kernel(resultq%kernel,reduced)
         end if
         static_residual=maxval(abs(reduced-complete))
         write(*,'(a,i0,a,es12.4)') 'Rotation static Hessian reduction q index ',iq,' max residual = ',static_residual
         if(iq==1 .and. static_residual>2.0e-7_rp) error stop 'Q0_STATIC_LIMIT_OPEN'
         if(iq==2 .and. static_residual>2.0e-8_rp) error stop 'FINITE_Q_STATIC_LIMIT_OPEN'
         if(iq==2) call statem%clear()
         deallocate(values,endpoint_values,vectors,endpoint_vectors,weights,axes,torques_q,torques_minus_q,mixed, &
            hessian,torque_torque,mixed_contact,complete)
      end do

      ! Literal direct band-loop oracle: q=0 and a generic non-self-inverse q,
      ! both frequency signs, and two causal regulators.
      allocate(bubble_o(ncoord,ncoord),contact_o(ncoord,ncoord),kernel_o(ncoord,ncoord))
      oracle_residual=0.0_rp
      do iq=1,2
         q=qset(:,iq)
         if(iq==1) then
            request%state=>state0
         else
            request%state=>stateq
         end if
         do ie=1,2
            eta=eta_set(ie)
            do iw=1,2
               omega=merge(-2.0e-4_rp,2.0e-4_rp,iw==1)
               request%exact_static=.false.; request%omega=omega; request%eta=eta; request%want_inverse=.false.
               call evaluate_rotation_response(request,plus)
               call evaluate_rotation_response_oracle(base,recip,q,omega,eta,bubble_o,contact_o,kernel_o)
               oracle_residual=max(oracle_residual,maxval(abs(plus%bubble-bubble_o)), &
                  maxval(abs(plus%contact-contact_o)),maxval(abs(plus%kernel-kernel_o)))
            end do
         end do
      end do
      write(*,'(a,es12.4)') 'Rotation direct band-loop oracle max absolute residual = ',oracle_residual
      if(oracle_residual>2.0e-10_rp) error stop 'DYNAMIC_BAND_ORACLE_OPEN'

      ! Use a genuine non-self-inverse reciprocal-mesh translation for the
      ! exact discrete q/-q/-omega covariance gate.  The direct oracle above
      ! still covers the generic noncommensurate q=1/8.
      q_cov=[1.0_rp/real(recip%nk_mesh(1),rp),0.0_rp,0.0_rp]
      call stateq%clear(); call statem%clear()
      grid_residual=0.0_rp
      do ik=1,nk
         grid_point=recip%k_points(:,ik)+q_cov
         grid_point=grid_point-floor(grid_point+0.5_rp)
         grid_best=huge(1.0_rp)
         do jk=1,nk
            grid_best=min(grid_best,maxval(abs(grid_point-recip%k_points(:,jk))))
         end do
         grid_residual=max(grid_residual,grid_best)
      end do
      write(*,'(a,es12.4)') 'Rotation covariance mesh-translation residual = ',grid_residual
      call prepare_rotation_response(base,recip,q_cov,stateq)
      call prepare_rotation_response(base,recip,-q_cov,statem)
      request%state=>stateq; request%omega=2.0e-4_rp; request%eta=5.0e-5_rp; request%exact_static=.false.
      call evaluate_rotation_response(request,plus)
      request%state=>statem; request%omega=-2.0e-4_rp
      call evaluate_rotation_response(request,minus)
      ! For the stored endpoint orientation the live reality relation keeps
      ! the same Cartesian labels: K_AB(q,w)=K_AB(-q,-w)^*.  The equivalent
      ! BA form in a transposed coordinate convention is not used here.
      covariance_residual=maxval(abs(plus%kernel-conjg(minus%kernel)))/ &
         max(1.0_rp,maxval(abs(plus%kernel)))
      write(*,'(a,3(es12.4,1x))') 'Rotation q/-q/-omega covariance q/residual = ',q_cov,covariance_residual
      if(grid_residual>1.0e-12_rp .or. covariance_residual>1.0e-8_rp) error stop 'Q_OMEGA_COVARIANCE_OPEN'
      call state0%clear(); call stateq%clear(); call statem%clear()
      deallocate(reduced,unitary,kernel_pm,bubble_o,contact_o,kernel_o)
   end subroutine check_rotation_response

   subroutine independent_berry_oracle(base,recip,m_band,berry)
      type(lmto_live_hamiltonian_fixture), intent(in) :: base
      type(reciprocal), intent(inout) :: recip
      real(rp), intent(out) :: m_band,berry
      real(rp), allocatable :: values(:, :)
      complex(rp), allocatable :: vectors(:, :, :),rho(:, :),gx(:, :),gy(:, :),commutator(:, :)
      complex(rp) :: trace_value
      real(rp) :: kT,weight,occ,spin_z,weight_sum
      integer :: nmat,nband,nk,nsite,norb,ik,n,site,io,up,dn,i,j
      nmat=2*base%norb*base%nsite; nband=nmat; nk=size(recip%k_points,2); nsite=base%nsite; norb=base%norb
      allocate(values(nband,nk),vectors(nmat,nband,nk),rho(nmat,nmat),gx(nmat,nmat),gy(nmat,nmat),commutator(nmat,nmat))
      call recip%calculate_eigenpairs_at_kpoints(recip%k_points,values,vectors)
      kT=max(recip%temperature*6.3336814e-6_rp,1.0e-10_rp); weight_sum=sum(recip%k_weights)
      rho=cmplx(0.0_rp,0.0_rp,rp); gx=rho; gy=rho
      do ik=1,nk
         weight=recip%k_weights(ik)/weight_sum
         do n=1,nband
            occ=finite_temperature_occupation(values(n,ik),recip%fermi_level,kT)
            do i=1,nmat
               do j=1,nmat
                  rho(i,j)=rho(i,j)+weight*occ*vectors(i,n,ik)*conjg(vectors(j,n,ik))
               end do
            end do
         end do
      end do
      do site=1,nsite
         do io=1,norb
            up=(site-1)*2*norb+io; dn=(site-1)*2*norb+norb+io
            gx(up,dn)=0.5_rp; gx(dn,up)=0.5_rp
            gy(up,dn)=-0.5_rp*i_unit; gy(dn,up)=0.5_rp*i_unit
         end do
      end do
      commutator=matmul(gx,gy)-matmul(gy,gx)
      trace_value=cmplx(0.0_rp,0.0_rp,rp); m_band=0.0_rp
      do i=1,nmat
         do j=1,nmat
            trace_value=trace_value+rho(i,j)*commutator(j,i)
         end do
      end do
      do site=1,nsite
         do io=1,norb
            up=(site-1)*2*norb+io; dn=(site-1)*2*norb+norb+io
            m_band=m_band+real(rho(up,up)-rho(dn,dn),rp)
         end do
      end do
      berry=real(-i_unit*trace_value,rp)
      deallocate(values,vectors,rho,gx,gy,commutator)
   end subroutine independent_berry_oracle

   subroutine check_material_rotation_adapter(base,qpoints,axis)
      type(lmto_live_hamiltonian_fixture), intent(in) :: base
      real(rp), intent(in) :: qpoints(:, :),axis(3)
      type(lmto_live_hamiltonian_fixture) :: super
      real(rp), parameter :: fd_steps(3)=[0.04_rp,0.01_rp,0.005_rp]
      complex(rp), allocatable :: hs(:, :)
      complex(rp), allocatable :: bplus(:, :), bminus(:, :), b0(:, :), analytic(:, :), fd(:, :)
      complex(rp), allocatable :: bpp(:, :), bpm(:, :), bmp(:, :), bmm(:, :), analytic_minus(:, :), analytic_pair(:, :)
      real(rp) :: torque_error(size(qpoints,2),size(fd_steps)), contact_error(size(qpoints,2),size(fd_steps))
      real(rp) :: mixed_xy_error(size(qpoints,2),size(fd_steps))
      real(rp) :: factor, value, axis_y(3)
      real(rp), parameter :: k_test(3)=[0.125_rp,0.0_rp,0.0_rp]
      integer :: ncell,nlocal,nsuper,iq,is

      ! q=(1/8,0,0) is represented by an eight-cell direct supercell.  The
      ! independent oracle rotates real-space moments, assembles H2, and then
      ! Fourier projects; it does not call either derivative routine.
      if (base%nsite/=1) error stop 'accepted Fe oracle currently requires one primitive sublattice'
      ncell=8; nlocal=2*base%norb*base%nsite; nsuper=ncell*nlocal
      call lift_material_fixture(base,ncell,super)
      axis_y=[0.0_rp,1.0_rp,0.0_rp]
      allocate(hs(nsuper,nsuper),bplus(nlocal,nlocal),bminus(nlocal,nlocal),b0(nlocal,nlocal), &
         analytic(nlocal,nlocal),fd(nlocal,nlocal),bpp(nlocal,nlocal),bpm(nlocal,nlocal), &
         bmp(nlocal,nlocal),bmm(nlocal,nlocal),analytic_minus(nlocal,nlocal),analytic_pair(nlocal,nlocal))
      do iq=1,size(qpoints,2)
         call assemble_lmto_finite_q_torque(base,k_test,qpoints(:,iq),1,axis,analytic)
         factor=1.0_rp
         if (maxval(abs(qpoints(:,iq)))>1.0e-12_rp) factor=2.0_rp
         do is=1,size(fd_steps)
            call material_supercell_block(super,base,qpoints(:,iq),axis,fd_steps(is),k_test+qpoints(:,iq),k_test,hs,bplus)
            call material_supercell_block(super,base,qpoints(:,iq),axis,-fd_steps(is),k_test+qpoints(:,iq),k_test,hs,bminus)
            fd=factor*(bplus-bminus)/(2.0_rp*fd_steps(is))
            torque_error(iq,is)=matrix_relative(fd,analytic)
         end do
         write(*,'(a,i0,a,3(es14.6,1x))') 'accepted-Fe second-order torque FD errors q index ',iq,' = ',torque_error(iq,:)

         call assemble_lmto_finite_q_mixed_derivative(base,k_test,qpoints(:,iq), &
            1,axis,1,axis,analytic)
         do is=1,size(fd_steps)
            value=fd_steps(is)
            call material_supercell_block(super,base,qpoints(:,iq),axis,2.0_rp*value,k_test,k_test,hs,bplus)
            call material_supercell_block(super,base,qpoints(:,iq),axis,0.0_rp,k_test,k_test,hs,b0)
            call material_supercell_block(super,base,qpoints(:,iq),axis,-2.0_rp*value,k_test,k_test,hs,bminus)
            fd=factor*(bplus-b0-b0+bminus)/(4.0_rp*value*value)
            contact_error(iq,is)=matrix_relative(fd,analytic)
         end do
         write(*,'(a,i0,a,3(es14.6,1x))') 'accepted-Fe second-order contact four-point errors q index ',iq,' = ',contact_error(iq,:)
         if (torque_error(iq,3)>3.0e-5_rp .or. torque_error(iq,2)>0.45_rp*torque_error(iq,1)) &
            error stop 'SECOND_ORDER_TORQUE_OPEN on accepted Fe fixture'
         if (contact_error(iq,3)>5.0e-5_rp .or. contact_error(iq,2)>0.45_rp*contact_error(iq,1)) &
            error stop 'SECOND_ORDER_CONTACT_OPEN on accepted Fe fixture'

         ! Independent four-point supercell oracle for the mixed transverse
         ! contact.  The exponential rotation coordinates generate the
         ! symmetrized x/y second moment in the accepted Hamiltonian.
         call assemble_lmto_finite_q_mixed_derivative(base,k_test,qpoints(:,iq), &
            1,axis,1,axis_y,analytic)
         if (maxval(abs(qpoints(:,iq)))>1.0e-12_rp) then
            call assemble_lmto_finite_q_mixed_derivative(base,k_test,-qpoints(:,iq), &
               1,axis,1,axis_y,analytic_minus)
            ! A real cosine perturbation contains both q and -q.  Its
            ! k-diagonal four-point derivative therefore checks the q-even
            ! combination of the two contact vertices.
            analytic_pair=0.5_rp*(analytic+analytic_minus)
         else
            analytic_pair=analytic
         end if
         do is=1,size(fd_steps)
            value=fd_steps(is)
            call material_supercell_block_xy(super,base,qpoints(:,iq),axis,value,axis_y,value,k_test,k_test,hs,bpp)
            call material_supercell_block_xy(super,base,qpoints(:,iq),axis,value,axis_y,-value,k_test,k_test,hs,bpm)
            call material_supercell_block_xy(super,base,qpoints(:,iq),axis,-value,axis_y,value,k_test,k_test,hs,bmp)
            call material_supercell_block_xy(super,base,qpoints(:,iq),axis,-value,axis_y,-value,k_test,k_test,hs,bmm)
            fd=factor*(bpp-bpm-bmp+bmm)/(4.0_rp*value*value)
            mixed_xy_error(iq,is)=sqrt(sum(abs(fd-analytic_pair)**2))/ &
               max(1.0_rp,sqrt(sum(abs(analytic_pair)**2)))
         end do
         write(*,'(a,i0,a,3(es14.6,1x))') 'accepted-Fe mixed xy q-even contact four-point errors q index ', &
            iq,' = ',mixed_xy_error(iq,:)
         if (mixed_xy_error(iq,3)>8.0e-5_rp .or. &
             (mixed_xy_error(iq,1)>1.0e-10_rp .and. mixed_xy_error(iq,2)>0.5_rp*mixed_xy_error(iq,1))) &
            error stop 'SECOND_ORDER_MIXED_XY_CONTACT_OPEN on accepted Fe fixture'
      end do
      call super%clear()
   end subroutine check_material_rotation_adapter

   subroutine lift_material_fixture(base,ncell,super)
      type(lmto_live_hamiltonian_fixture), intent(in) :: base
      integer, intent(in) :: ncell
      type(lmto_live_hamiltonian_fixture), intent(out) :: super
      integer :: cell,site,ib,dbond,bsource,btarget,target_cell
      call lmto_fixture_init(super,ncell*base%nsite,base%norb,ncell*base%nbond,base%hoh)
      super%include_enu=base%include_enu
      super%cartesian_to_spherical=base%cartesian_to_spherical
      do cell=0,ncell-1
         do site=1,base%nsite
            bsource=cell*base%nsite+site
            super%moments(:,bsource)=base%moments(:,site)
            super%site_position(:,bsource)=base%site_position(:,site)+[real(cell,rp),0.0_rp,0.0_rp]
            super%wx0(:,bsource)=base%wx0(:,site); super%wx1(:,bsource)=base%wx1(:,site)
            super%c0(:,bsource)=base%c0(:,site); super%c1(:,bsource)=base%c1(:,site)
            super%obar0(:,bsource)=base%obar0(:,site); super%obar1(:,bsource)=base%obar1(:,site)
            super%enu0(:,bsource)=base%enu0(:,site); super%enu1(:,bsource)=base%enu1(:,site)
         end do
         do ib=1,base%nbond
            dbond=cell*base%nbond+ib
            if (abs(base%bond_vector(1,ib)-real(nint(base%bond_vector(1,ib)),rp))>1.0e-8_rp) &
               error stop 'accepted Fe supercell oracle requires integer fractional-x bond translations'
            target_cell=modulo(cell+nint(base%bond_vector(1,ib)),ncell)
            super%bond_source(dbond)=cell*base%nsite+base%bond_source(ib)
            super%bond_target(dbond)=target_cell*base%nsite+base%bond_target(ib)
            super%bond_vector(:,dbond)=0.0_rp
            super%onsite(dbond)=base%onsite(ib); super%hhh(:,:,dbond)=base%hhh(:,:,ib)
         end do
      end do
   end subroutine lift_material_fixture

   subroutine material_supercell_block_xy(super,base,qpattern,axis_x,amplitude_x,axis_y,amplitude_y,krow,kcol,hs,block)
      type(lmto_live_hamiltonian_fixture), intent(inout) :: super
      type(lmto_live_hamiltonian_fixture), intent(in) :: base
      real(rp), intent(in) :: qpattern(3),axis_x(3),amplitude_x,axis_y(3),amplitude_y,krow(3),kcol(3)
      complex(rp), intent(out) :: hs(:, :),block(:, :)
      call set_material_rotation_xy(super,base,qpattern,axis_x,amplitude_x,axis_y,amplitude_y)
      call assemble_lmto_hamiltonian(super,[0.0_rp,0.0_rp,0.0_rp],hs)
      call supercell_fourier(hs,base%norb,base%nsite,krow,kcol,block)
   end subroutine material_supercell_block_xy

   subroutine material_supercell_block(super,base,qpattern,axis,amplitude,krow,kcol,hs,block)
      type(lmto_live_hamiltonian_fixture), intent(inout) :: super
      type(lmto_live_hamiltonian_fixture), intent(in) :: base
      real(rp), intent(in) :: qpattern(3),axis(3),amplitude,krow(3),kcol(3)
      complex(rp), intent(out) :: hs(:, :),block(:, :)
      call set_material_rotation(super,base,qpattern,axis,amplitude)
      call assemble_lmto_hamiltonian(super,[0.0_rp,0.0_rp,0.0_rp],hs)
      call supercell_fourier(hs,base%norb,base%nsite,krow,kcol,block)
   end subroutine material_supercell_block

   subroutine set_material_rotation(super,base,qpattern,axis,amplitude)
      type(lmto_live_hamiltonian_fixture), intent(inout) :: super
      type(lmto_live_hamiltonian_fixture), intent(in) :: base
      real(rp), intent(in) :: qpattern(3),axis(3),amplitude
      integer :: cell,site,super_site
      real(rp) :: angle
      do cell=0,size(super%moments,2)/base%nsite-1
         angle=amplitude*cos(2.0_rp*acos(-1.0_rp)*qpattern(1)*real(cell,rp))
         do site=1,base%nsite
            super_site=cell*base%nsite+site
            super%moments(:,super_site)=rotate_real_moment(base%moments(:,site),axis,angle)
         end do
      end do
   end subroutine set_material_rotation

   subroutine set_material_rotation_xy(super,base,qpattern,axis_x,amplitude_x,axis_y,amplitude_y)
      type(lmto_live_hamiltonian_fixture), intent(inout) :: super
      type(lmto_live_hamiltonian_fixture), intent(in) :: base
      real(rp), intent(in) :: qpattern(3),axis_x(3),amplitude_x,axis_y(3),amplitude_y
      integer :: cell,site,super_site
      real(rp) :: angle,rotation_vector(3)
      do cell=0,size(super%moments,2)/base%nsite-1
         angle=2.0_rp*acos(-1.0_rp)*qpattern(1)*real(cell,rp)
         rotation_vector=cos(angle)*(amplitude_x*axis_x+amplitude_y*axis_y)
         do site=1,base%nsite
            super_site=cell*base%nsite+site
            super%moments(:,super_site)=rotate_real_moment_vector(base%moments(:,site),rotation_vector)
         end do
      end do
   end subroutine set_material_rotation_xy

   function rotate_real_moment_vector(moment,rotation_vector) result(rotated)
      real(rp), intent(in) :: moment(3),rotation_vector(3)
      real(rp) :: rotated(3),angle,axis(3),cross1(3),cross2(3)
      angle=sqrt(dot_product(rotation_vector,rotation_vector))
      if (angle<=1.0e-12_rp) then
         cross1=[rotation_vector(2)*moment(3)-rotation_vector(3)*moment(2), &
            rotation_vector(3)*moment(1)-rotation_vector(1)*moment(3), &
            rotation_vector(1)*moment(2)-rotation_vector(2)*moment(1)]
         cross2=[rotation_vector(2)*cross1(3)-rotation_vector(3)*cross1(2), &
            rotation_vector(3)*cross1(1)-rotation_vector(1)*cross1(3), &
            rotation_vector(1)*cross1(2)-rotation_vector(2)*cross1(1)]
         rotated=moment+cross1+0.5_rp*cross2
      else
         axis=rotation_vector/angle
         rotated=rotate_real_moment(moment,axis,angle)
      end if
   end function rotate_real_moment_vector

   function rotate_real_moment(moment,axis,angle) result(rotated)
      real(rp), intent(in) :: moment(3),axis(3),angle
      real(rp) :: rotated(3),cross(3),dot
      cross=[axis(2)*moment(3)-axis(3)*moment(2),axis(3)*moment(1)-axis(1)*moment(3), &
         axis(1)*moment(2)-axis(2)*moment(1)]
      dot=sum(axis*moment)
      rotated=moment*cos(angle)+cross*sin(angle)+axis*dot*(1.0_rp-cos(angle))
   end function rotate_real_moment

   subroutine supercell_fourier(h,norb,nsite,krow,kcol,block)
      complex(rp), intent(in) :: h(:, :)
      integer, intent(in) :: norb,nsite
      real(rp), intent(in) :: krow(3),kcol(3)
      complex(rp), intent(out) :: block(:, :)
      complex(rp) :: phase
      real(rp) :: angle
      integer :: ncell,nlocal,row_first,row_last,col_first,col_last,row_cell,col_cell
      nlocal=2*norb*nsite; ncell=size(h,1)/nlocal
      if (size(h,1)/=size(h,2) .or. size(block,1)/=nlocal .or. size(block,2)/=nlocal .or. &
          ncell<1 .or. mod(size(h,1),nlocal)/=0) error stop 'accepted Fe supercell Fourier shape mismatch'
      block=cmplx(0.0_rp,0.0_rp,rp)
      do row_cell=0,ncell-1
         row_first=row_cell*nlocal+1; row_last=(row_cell+1)*nlocal
         do col_cell=0,ncell-1
            col_first=col_cell*nlocal+1; col_last=(col_cell+1)*nlocal
            angle=2.0_rp*acos(-1.0_rp)*(kcol(1)*real(col_cell,rp)-krow(1)*real(row_cell,rp))
            phase=cmplx(cos(angle),sin(angle),rp)
            ! Equivalent to exp(+i 2*pi*(kcol.R_col-krow.R_row)).
            block=block+phase*h(row_first:row_last,col_first:col_last)/real(ncell,rp)
         end do
      end do
   end subroutine supercell_fourier

   function matrix_relative(a,b) result(value)
      complex(rp), intent(in) :: a(:, :),b(:, :)
      real(rp) :: value
      value=sqrt(sum(abs(a-b)**2))/max(sqrt(sum(abs(b)**2)),tiny(1.0_rp))
   end function matrix_relative

end program test_dresp03q_production_adapter
