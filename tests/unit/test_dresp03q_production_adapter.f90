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
   use lr_lmto_turek_contour_mod, only: native_turek_contour_options, native_turek_contour_report, &
      native_turek_static_reference
   use math_mod, only: ang2au, init_math_operators
   use precision_mod, only: rp
   use reciprocal_mod, only: reciprocal
   use timer_mod, only: g_timer, timer
   implicit none

   type(control) :: ctl
   type(lattice), target :: lat
   type(charge), target :: chg
   type(hamiltonian), target :: ham
   type(energy) :: ene
   type(reciprocal) :: recip
   type(native_turek_contour_options) :: native_options
   type(native_turek_contour_report) :: native_report
   type(lmto_live_hamiltonian_fixture) :: fixture
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
   call fixture%clear()
   deallocate(h_production,h_fixture,adapter_points,weights,eigenvalues,eigenvectors,endpoint_values,endpoint_vectors, &
      shifted_points,torque_q,torque_minus_q,mixed,jq_ud,jq_du,jq_sym,native_delta,native_curvature)

contains

   subroutine check_material_rotation_adapter(base,qpoints,axis)
      type(lmto_live_hamiltonian_fixture), intent(in) :: base
      real(rp), intent(in) :: qpoints(:, :),axis(3)
      type(lmto_live_hamiltonian_fixture) :: super
      real(rp), parameter :: fd_steps(3)=[0.04_rp,0.01_rp,0.005_rp]
      complex(rp), allocatable :: hs(:, :)
      complex(rp), allocatable :: bplus(:, :), bminus(:, :), b0(:, :), analytic(:, :), fd(:, :)
      real(rp) :: torque_error(size(qpoints,2),size(fd_steps)), contact_error(size(qpoints,2),size(fd_steps))
      real(rp) :: factor, value
      real(rp), parameter :: k_test(3)=[0.125_rp,0.0_rp,0.0_rp]
      integer :: ncell,nlocal,nsuper,iq,is

      ! q=(1/8,0,0) is represented by an eight-cell direct supercell.  The
      ! independent oracle rotates real-space moments, assembles H2, and then
      ! Fourier projects; it does not call either derivative routine.
      if (base%nsite/=1) error stop 'accepted Fe oracle currently requires one primitive sublattice'
      ncell=8; nlocal=2*base%norb*base%nsite; nsuper=ncell*nlocal
      call lift_material_fixture(base,ncell,super)
      allocate(hs(nsuper,nsuper),bplus(nlocal,nlocal),bminus(nlocal,nlocal),b0(nlocal,nlocal), &
         analytic(nlocal,nlocal),fd(nlocal,nlocal))
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
