!------------------------------------------------------------------------------
! Native second-order local-rotation response.
!
! Coordinates are site-local rigid rotations about two orthonormal transverse
! axes.  For A=(site,axis), the accepted finite-q torque has rows at k+q and
! columns at k.  The retarded electronic term is kept ordered in A,B; only the
! exact static comparison applies the Hermitian/static reduction.
!------------------------------------------------------------------------------
module lr_rotation_response_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use precision_mod, only: rp
   use math_mod, only: i_unit
   use lr_kl_hessian_mod, only: lmto_live_hamiltonian_fixture, &
      assemble_lmto_finite_q_torque, assemble_lmto_finite_q_mixed_derivative, &
      finite_temperature_occupation, fermi_divided_difference
   use reciprocal_mod, only: reciprocal
   implicit none
   private

   real(rp), parameter :: kB_ry_per_k = 6.3336814e-6_rp
   real(rp), parameter :: q_tolerance = 1.0e-12_rp

   type, public :: rotation_state
      type(lmto_live_hamiltonian_fixture), pointer :: fixture => null()
      type(reciprocal), pointer :: reciprocal_state => null()
      real(rp) :: q(3) = 0.0_rp
      integer :: nsite = 0, ncoord = 0, nmat = 0, nbands = 0, nk = 0
      real(rp) :: fermi = 0.0_rp, kT = 0.0_rp, weight_sum = 0.0_rp
      real(rp) :: magnetization = 0.0_rp, berry = 0.0_rp, berry_residual = 0.0_rp
      real(rp) :: electron_count = 0.0_rp, electron_residual = 0.0_rp
      logical :: prepared = .false.
      real(rp), allocatable :: values(:, :), endpoint_values(:, :), weights(:)
      real(rp), allocatable :: occupations(:, :), endpoint_occupations(:, :)
      complex(rp), allocatable :: torque_band(:, :, :, :), partner_band(:, :, :, :)
      complex(rp), allocatable :: contact(:, :)
   contains
      procedure :: clear => rotation_state_clear
   end type rotation_state

   type, public :: rotation_request
      type(rotation_state), pointer :: state => null()
      real(rp) :: omega = 0.0_rp
      real(rp) :: eta = 0.0_rp
      logical :: exact_static = .false.
      logical :: want_inverse = .false.
   end type rotation_request

   type, public :: rotation_result
      real(rp) :: q(3) = 0.0_rp
      real(rp) :: omega = 0.0_rp, eta = 0.0_rp
      real(rp) :: berry = 0.0_rp, magnetization = 0.0_rp, berry_residual = 0.0_rp
      real(rp) :: electron_count = 0.0_rp, electron_residual = 0.0_rp
      logical :: inverse_available = .false.
      complex(rp), allocatable :: bubble(:, :), contact(:, :), kernel(:, :), kernel_pm(:, :)
      complex(rp), allocatable :: inverse_kernel(:, :), inverse_kernel_pm(:, :)
   end type rotation_result

   public :: prepare_rotation_response
   public :: evaluate_rotation_response
   public :: evaluate_rotation_response_oracle
   public :: reduce_static_rotation_kernel
   public :: rotation_axes
   public :: rotation_circular_unitary

   interface
      subroutine zgetrf(m, n, a, lda, ipiv, info)
         import :: rp
         integer, intent(in) :: m, n, lda
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         complex(rp), intent(inout) :: a(lda, *)
      end subroutine zgetrf
      subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
         import :: rp
         integer, intent(in) :: n, lda, lwork
         integer, intent(in) :: ipiv(*)
         integer, intent(out) :: info
         complex(rp), intent(inout) :: a(lda, *), work(*)
      end subroutine zgetri
   end interface

contains

   !> Freeze the accepted second-order state and pretransform both endpoints'
   !> torque vertices once.  The band-pair response then preserves Pi_AB and
   !> Pi_BA independently at every requested frequency.
   subroutine prepare_rotation_response(fixture, recip, q, state)
      type(lmto_live_hamiltonian_fixture), target, intent(in) :: fixture
      type(reciprocal), target, intent(inout) :: recip
      real(rp), intent(in) :: q(3)
      type(rotation_state), intent(inout) :: state
      real(rp), allocatable :: points(:, :), folded(:, :), axes(:, :)
      complex(rp), allocatable :: vectors(:, :, :), endpoint_vectors(:, :, :), vertex(:, :), mixed(:, :), band(:, :)
      integer :: ik, ia, ib, n, m, isite_a, isite_b, axis_a, axis_b, norb_site, i0, info
      real(rp) :: f, spin_z, norm_m, weight

      call state%clear()
      if (fixture%nsite < 1 .or. fixture%norb < 1 .or. .not. fixture%hoh .or. .not. fixture%include_enu) &
         error stop 'STATE_PROVENANCE_OPEN: response requires accepted H=B-QB+E_nu fixture'
      if (trim(recip%reciprocal_mode) /= 'ham_only' .or. trim(recip%kspace_ham_order) /= 'second') &
         error stop 'STATE_PROVENANCE_OPEN: response requires second-order ham_only reciprocal state'
      if (recip%include_so) error stop 'STATE_PROVENANCE_OPEN: native rotation response currently requires SOC off'
      if (.not. allocated(recip%k_points) .or. .not. allocated(recip%k_weights)) &
         error stop 'STATE_PROVENANCE_OPEN: reciprocal mesh is unavailable'
      if (fixture%nsite /= recip%lattice%nrec) error stop 'STATE_PROVENANCE_OPEN: fixture/state site counts differ'
      call recip%require_replicated_k_workset('prepare_rotation_response')

      state%fixture => fixture
      state%reciprocal_state => recip
      state%q = q
      state%nsite = fixture%nsite
      state%ncoord = 2*fixture%nsite
      state%nmat = 2*fixture%norb*fixture%nsite
      state%nk = size(recip%k_points,2)
      state%nbands = state%nmat
      state%fermi = recip%fermi_level
      state%kT = max(recip%temperature*kB_ry_per_k, 1.0e-10_rp)
      state%weight_sum = sum(recip%k_weights)
      if (state%nk < 1 .or. size(recip%k_weights) /= state%nk .or. state%weight_sum <= tiny(1.0_rp)) &
         error stop 'STATE_PROVENANCE_OPEN: invalid reciprocal k weights'
      if (abs(state%weight_sum-1.0_rp) > 1.0e-8_rp) &
         error stop 'STATE_PROVENANCE_OPEN: reciprocal k weights are not normalized'

      allocate(state%values(state%nbands,state%nk), state%endpoint_values(state%nbands,state%nk), &
         state%weights(state%nk), state%occupations(state%nbands,state%nk), &
         state%endpoint_occupations(state%nbands,state%nk), &
         state%torque_band(state%nbands,state%nbands,state%ncoord,state%nk), &
         state%partner_band(state%nbands,state%nbands,state%ncoord,state%nk), &
         state%contact(state%ncoord,state%ncoord))
      state%weights = recip%k_weights/state%weight_sum

      ! The normal endpoint comes from the accepted SCF eigensystem whenever
      ! it is present; the arbitrary-point API supplies every k+q endpoint.
      if (allocated(recip%eigenvalues) .and. allocated(recip%eigenvectors)) then
         if (all(shape(recip%eigenvalues) == [state%nbands,state%nk]) .and. &
             all(shape(recip%eigenvectors) == [state%nmat,state%nbands,state%nk])) then
            state%values = recip%eigenvalues
            allocate(vectors(state%nmat,state%nbands,state%nk))
            vectors = recip%eigenvectors
         else
            call recip%calculate_eigenpairs_at_kpoints(recip%k_points,state%values,vectors)
         end if
      else
         call recip%calculate_eigenpairs_at_kpoints(recip%k_points,state%values,vectors)
      end if
      if (maxval(abs(q)) <= q_tolerance) then
         state%endpoint_values = state%values
         allocate(endpoint_vectors(state%nmat,state%nbands,state%nk))
         endpoint_vectors = vectors
      else
         points = recip%k_points + spread(q,2,state%nk)
         call recip%calculate_eigenpairs_at_kpoints(points,state%endpoint_values,endpoint_vectors,folded)
      end if

      do ik = 1, state%nk
         do n = 1, state%nbands
            state%occupations(n,ik) = finite_temperature_occupation(state%values(n,ik),state%fermi,state%kT)
            state%endpoint_occupations(n,ik) = finite_temperature_occupation(state%endpoint_values(n,ik),state%fermi,state%kT)
         end do
      end do
      call rotation_axes(fixture%moments,axes)
      state%contact = cmplx(0.0_rp,0.0_rp,rp)
      state%magnetization = 0.0_rp
      state%electron_count = 0.0_rp
      norb_site = fixture%norb
      allocate(vertex(state%nmat,state%nmat),mixed(state%nmat,state%nmat),band(state%nbands,state%nbands))

      do ik = 1, state%nk
         weight = state%weights(ik)
         do n = 1, state%nbands
            f = state%occupations(n,ik)
            state%electron_count = state%electron_count + weight*f
            spin_z = 0.0_rp
            do isite_a = 1, state%nsite
               i0 = (isite_a-1)*2*norb_site
               spin_z = spin_z + sum(abs(vectors(i0+1:i0+norb_site,n,ik))**2) - &
                  sum(abs(vectors(i0+norb_site+1:i0+2*norb_site,n,ik))**2)
            end do
            state%magnetization = state%magnetization + weight*f*spin_z
         end do

         do ia = 1, state%ncoord
            isite_a = (ia+1)/2
            axis_a = 1+mod(ia-1,2)
            call assemble_lmto_finite_q_torque(fixture,recip%k_points(:,ik),q,isite_a,axes(:,2*(isite_a-1)+axis_a),vertex)
            band = matmul(conjg(transpose(endpoint_vectors(:,:,ik))),matmul(vertex,vectors(:,:,ik)))
            state%torque_band(:,:,ia,ik) = band
            call assemble_lmto_finite_q_torque(fixture,recip%k_points(:,ik)+q,-q,isite_a, &
               axes(:,2*(isite_a-1)+axis_a),vertex)
            band = matmul(conjg(transpose(vectors(:,:,ik))),matmul(vertex,endpoint_vectors(:,:,ik)))
            state%partner_band(:,:,ia,ik) = band
         end do

         do ia = 1, state%ncoord
            isite_a = (ia+1)/2
            axis_a = 1+mod(ia-1,2)
            do ib = 1, state%ncoord
               isite_b = (ib+1)/2
               axis_b = 1+mod(ib-1,2)
               call assemble_lmto_finite_q_mixed_derivative(fixture,recip%k_points(:,ik),q,isite_a, &
                  axes(:,2*(isite_a-1)+axis_a),isite_b,axes(:,2*(isite_b-1)+axis_b),mixed)
               band = matmul(conjg(transpose(vectors(:,:,ik))),matmul(mixed,vectors(:,:,ik)))
               do n = 1, state%nbands
                  state%contact(ia,ib) = state%contact(ia,ib) + weight*state%occupations(n,ik)*band(n,n)
               end do
            end do
         end do
      end do
      call global_spin_berry(vectors,state%occupations,state%weights,state%nsite,norb_site,state%berry)
      state%berry_residual=abs(state%berry-0.5_rp*state%magnetization)/ &
         max(abs(state%berry),abs(state%magnetization),1.0e-12_rp)
      state%electron_residual = state%electron_count-recip%total_electrons
      state%prepared = .true.
      deallocate(vectors,endpoint_vectors,axes,vertex,mixed,band)
      if (allocated(points)) deallocate(points)
      if (allocated(folded)) deallocate(folded)
   end subroutine prepare_rotation_response

   !> Evaluate the literal unsymmetrized retarded band sum.  At omega=eta=0,
   !> exact_static selects the finite-temperature divided-difference branch.
   subroutine evaluate_rotation_response(request,result)
      type(rotation_request), intent(in) :: request
      type(rotation_result), intent(out) :: result
      type(rotation_state), pointer :: state
      complex(rp), allocatable :: unitary(:, :)
      complex(rp) :: denominator, product
      real(rp) :: fn, fm, divided, energy_difference, weight
      integer :: ia, ib, ik, n, m, info, ncoord

      if (.not. associated(request%state)) error stop 'rotation response request has no prepared state'
      state => request%state
      if (.not. state%prepared) error stop 'rotation response state is not prepared'
      if (request%exact_static) then
         if (abs(request%omega) > q_tolerance .or. abs(request%eta) > q_tolerance) &
            error stop 'exact static response requires omega=eta=0'
      else if (request%eta <= 0.0_rp .or. .not.ieee_is_finite(request%eta)) then
         error stop 'retarded rotation response requires finite positive eta'
      end if
      ncoord = state%ncoord
      result%q = state%q
      result%omega = request%omega
      result%eta = request%eta
      result%berry = state%berry
      result%magnetization = state%magnetization
      result%berry_residual = state%berry_residual
      result%electron_count = state%electron_count
      result%electron_residual = state%electron_residual
      allocate(result%bubble(ncoord,ncoord),result%contact(ncoord,ncoord), &
         result%kernel(ncoord,ncoord),result%kernel_pm(ncoord,ncoord),unitary(ncoord,ncoord))
      result%bubble = cmplx(0.0_rp,0.0_rp,rp)
      do ik = 1, state%nk
         weight = state%weights(ik)
         do ia = 1, ncoord
            do ib = 1, ncoord
               do n = 1, state%nbands
                  fn = state%occupations(n,ik)
                  do m = 1, state%nbands
                     fm = state%endpoint_occupations(m,ik)
                     product = state%torque_band(m,n,ia,ik)*state%partner_band(n,m,ib,ik)
                     energy_difference = state%values(n,ik)-state%endpoint_values(m,ik)
                     if (request%exact_static) then
                        divided = fermi_divided_difference(state%values(n,ik),state%endpoint_values(m,ik), &
                           state%fermi,state%kT)
                        result%bubble(ia,ib) = result%bubble(ia,ib)+weight*divided*product
                     else
                        denominator = cmplx(request%omega+energy_difference,request%eta,rp)
                        result%bubble(ia,ib) = result%bubble(ia,ib)+weight*(fn-fm)*product/denominator
                     end if
                  end do
               end do
            end do
         end do
      end do
      result%contact = state%contact
      result%kernel = result%contact+result%bubble
      call rotation_circular_unitary(state%nsite,unitary)
      result%kernel_pm = matmul(unitary,matmul(result%kernel,conjg(transpose(unitary))))
      if (request%want_inverse) then
         ! The exact q=0 static Goldstone kernel is singular by construction.
         if (request%exact_static .and. maxval(abs(state%q)) <= q_tolerance) then
            result%inverse_available = .false.
         else
            allocate(result%inverse_kernel(ncoord,ncoord),result%inverse_kernel_pm(ncoord,ncoord))
            result%inverse_kernel = result%kernel
            call invert_complex_matrix(result%inverse_kernel,info)
            if (info == 0) then
               result%inverse_available = .true.
               result%inverse_kernel_pm = matmul(unitary, &
                  matmul(result%inverse_kernel,conjg(transpose(unitary))))
            else
               deallocate(result%inverse_kernel,result%inverse_kernel_pm)
               result%inverse_available = .false.
            end if
         end if
      end if
      deallocate(unitary)
   end subroutine evaluate_rotation_response

   !> Independent literal band-loop oracle.  It assembles orbital-space
   !> vertices and evaluates <m,k+q|T_A|n,k><n,k|T_B|m,k+q> directly; it does
   !> not call or consume the production accumulator or its transformed bands.
   subroutine evaluate_rotation_response_oracle(fixture,recip,q,omega,eta,bubble,contact,kernel)
      type(lmto_live_hamiltonian_fixture), intent(in) :: fixture
      type(reciprocal), intent(inout) :: recip
      real(rp), intent(in) :: q(3),omega,eta
      complex(rp), intent(out) :: bubble(:, :),contact(:, :),kernel(:, :)
      real(rp), allocatable :: values(:, :),endpoint_values(:, :),points(:, :),axes(:, :)
      complex(rp), allocatable :: vectors(:, :, :),endpoint_vectors(:, :, :),vq(:, :, :),vm(:, :, :),mixed(:, :)
      complex(rp) :: amp_a,amp_b,amp_c,denominator
      real(rp) :: weight_sum,weight,fermi,kT,fn,fm
      integer :: nk,nband,nmat,ncoord,nsite,ik,ia,ib,n,m,site_a,site_b,axis_a,axis_b
      integer :: norb_site,i0

      if (eta<=0.0_rp) error stop 'dynamic band oracle requires positive eta'
      if (.not.fixture%hoh .or. .not.fixture%include_enu .or. trim(recip%kspace_ham_order)/='second' .or. &
          trim(recip%reciprocal_mode)/='ham_only' .or. recip%include_so) &
         error stop 'STATE_PROVENANCE_OPEN: direct oracle requires second-order ham_only state with SOC off'
      call recip%require_replicated_k_workset('evaluate_rotation_response_oracle')
      nk=size(recip%k_points,2); nsite=fixture%nsite; ncoord=2*nsite; nmat=2*fixture%norb*nsite; nband=nmat
      if (any(shape(bubble)/=[ncoord,ncoord]) .or. any(shape(contact)/=[ncoord,ncoord]) .or. &
          any(shape(kernel)/=[ncoord,ncoord])) error stop 'evaluate_rotation_response_oracle: shape mismatch'
      fermi=recip%fermi_level; kT=max(recip%temperature*kB_ry_per_k,1.0e-10_rp); weight_sum=sum(recip%k_weights)
      allocate(values(nband,nk),endpoint_values(nband,nk),vectors(nmat,nband,nk), &
         endpoint_vectors(nmat,nband,nk),axes(3,ncoord),vq(nmat,nmat,ncoord),vm(nmat,nmat,ncoord),mixed(nmat,nmat))
      if (allocated(recip%eigenvalues) .and. allocated(recip%eigenvectors)) then
         if (all(shape(recip%eigenvalues)==[nband,nk]) .and. all(shape(recip%eigenvectors)==[nmat,nband,nk])) then
            values=recip%eigenvalues; vectors=recip%eigenvectors
         else
            call recip%calculate_eigenpairs_at_kpoints(recip%k_points,values,vectors)
         end if
      else
         call recip%calculate_eigenpairs_at_kpoints(recip%k_points,values,vectors)
      end if
      if (maxval(abs(q))<=q_tolerance) then
         endpoint_values=values; endpoint_vectors=vectors
      else
         points=recip%k_points+spread(q,2,nk)
         call recip%calculate_eigenpairs_at_kpoints(points,endpoint_values,endpoint_vectors)
         deallocate(points)
      end if
      call rotation_axes(fixture%moments,axes)
      bubble=cmplx(0.0_rp,0.0_rp,rp); contact=cmplx(0.0_rp,0.0_rp,rp)
      norb_site=fixture%norb
      do ik=1,nk
         weight=recip%k_weights(ik)/weight_sum
         do ia=1,ncoord
            site_a=(ia+1)/2; axis_a=1+mod(ia-1,2)
            call assemble_lmto_finite_q_torque(fixture,recip%k_points(:,ik),q,site_a, &
               axes(:,2*(site_a-1)+axis_a),vq(:,:,ia))
            call assemble_lmto_finite_q_torque(fixture,recip%k_points(:,ik)+q,-q,site_a, &
               axes(:,2*(site_a-1)+axis_a),vm(:,:,ia))
         end do
         do ia=1,ncoord
            site_a=(ia+1)/2; axis_a=1+mod(ia-1,2)
            do ib=1,ncoord
               site_b=(ib+1)/2; axis_b=1+mod(ib-1,2)
               call assemble_lmto_finite_q_mixed_derivative(fixture,recip%k_points(:,ik),q,site_a, &
                  axes(:,2*(site_a-1)+axis_a),site_b,axes(:,2*(site_b-1)+axis_b),mixed)
               do n=1,nband
                  fn=finite_temperature_occupation(values(n,ik),fermi,kT)
                  amp_c=dot_product(vectors(:,n,ik),matmul(mixed,vectors(:,n,ik)))
                  contact(ia,ib)=contact(ia,ib)+weight*fn*amp_c
                  do m=1,nband
                     fm=finite_temperature_occupation(endpoint_values(m,ik),fermi,kT)
                     denominator=cmplx(omega+values(n,ik)-endpoint_values(m,ik),eta,rp)
                     amp_a=dot_product(endpoint_vectors(:,m,ik),matmul(vq(:,:,ia),vectors(:,n,ik)))
                     amp_b=dot_product(vectors(:,n,ik),matmul(vm(:,:,ib),endpoint_vectors(:,m,ik)))
                     bubble(ia,ib)=bubble(ia,ib)+weight*(fn-fm)*amp_a*amp_b/denominator
                  end do
               end do
            end do
         end do
      end do
      kernel=contact+bubble
      deallocate(values,endpoint_values,vectors,endpoint_vectors,axes,vq,vm,mixed)
   end subroutine evaluate_rotation_response_oracle

   !> The stored torque orientation makes the exact static force-theorem
   !> Hessian the Cartesian-index symmetric part at the same q:
   !> H_AB(q)=0.5*(K_AB^R(q,0)+K_BA^R(q,0)).  The q/-q covariance is a
   !> separate identity; inserting a conjugation here retains a spurious
   !> antisymmetric-in-frequency component at finite q.
   pure subroutine reduce_static_rotation_kernel(kernel_q,reduced)
      complex(rp), intent(in) :: kernel_q(:, :)
      complex(rp), intent(out) :: reduced(:, :)
      integer :: a,b,n
      n = size(kernel_q,1)
      if (size(kernel_q,2)/=n .or. any(shape(reduced)/=[n,n])) &
         error stop 'reduce_static_rotation_kernel: shape mismatch'
      do a = 1,n
         do b = 1,n
            reduced(a,b)=0.5_rp*(kernel_q(a,b)+kernel_q(b,a))
         end do
      end do
   end subroutine reduce_static_rotation_kernel

   !> Evaluate -i Tr rho [Gx,Gy] directly in the accepted band basis.
   !> The spin generators are global spin-1/2 generators replicated over
   !> sites/orbitals; the result should independently equal M_band/2.
   subroutine global_spin_berry(vectors,occupations,weights,nsite,norb,berry)
      complex(rp), intent(in) :: vectors(:, :, :)
      real(rp), intent(in) :: occupations(:, :),weights(:)
      integer, intent(in) :: nsite,norb
      real(rp), intent(out) :: berry
      complex(rp), allocatable :: gx(:, :),gy(:, :),commutator(:, :)
      complex(rp) :: trace_value
      integer :: nmat,nband,nk,site,io,up,dn,ik,n
      nmat=2*nsite*norb; nband=size(vectors,2); nk=size(vectors,3)
      if (any(shape(vectors)/=[nmat,nband,nk]) .or. any(shape(occupations)/=[nband,nk]) .or. size(weights)/=nk) &
         error stop 'global_spin_berry: shape mismatch'
      allocate(gx(nmat,nmat),gy(nmat,nmat),commutator(nmat,nmat))
      gx=cmplx(0.0_rp,0.0_rp,rp); gy=gx
      do site=1,nsite
         do io=1,norb
            up=(site-1)*2*norb+io; dn=(site-1)*2*norb+norb+io
            gx(up,dn)=0.5_rp; gx(dn,up)=0.5_rp
            gy(up,dn)=-0.5_rp*i_unit; gy(dn,up)=0.5_rp*i_unit
         end do
      end do
      commutator=matmul(gx,gy)-matmul(gy,gx)
      berry=0.0_rp
      do ik=1,nk
         do n=1,nband
            trace_value=dot_product(vectors(:,n,ik),matmul(commutator,vectors(:,n,ik)))
            berry=berry+weights(ik)*occupations(n,ik)*real(-i_unit*trace_value,rp)
         end do
      end do
      deallocate(gx,gy,commutator)
   end subroutine global_spin_berry

   !> Construct local transverse rotation axes from each accepted unit moment.
   pure subroutine rotation_axes(moments,axes)
      real(rp), intent(in) :: moments(:, :)
      real(rp), allocatable, intent(out) :: axes(:, :)
      real(rp) :: m(3), reference(3), e1(3), e2(3), norm_m, dot_m
      integer :: site
      if (size(moments,1)/=3) error stop 'rotation_axes: moment array must have three Cartesian rows'
      allocate(axes(3,2*size(moments,2)))
      do site = 1,size(moments,2)
         m = moments(:,site)
         norm_m = sqrt(dot_product(m,m))
         if (norm_m<=tiny(1.0_rp)) error stop 'rotation_axes: zero site moment'
         m = m/norm_m
         if (abs(m(1))<0.8_rp) then
            reference=[1.0_rp,0.0_rp,0.0_rp]
         else
            reference=[0.0_rp,1.0_rp,0.0_rp]
         end if
         dot_m=dot_product(reference,m)
         e1=reference-dot_m*m
         e1=e1/sqrt(dot_product(e1,e1))
         e2=[m(2)*e1(3)-m(3)*e1(2),m(3)*e1(1)-m(1)*e1(3),m(1)*e1(2)-m(2)*e1(1)]
         axes(:,2*site-1)=e1
         axes(:,2*site)=e2
      end do
   end subroutine rotation_axes

   !> theta_+=(theta_x-i theta_y)/sqrt(2), theta_-=(theta_x+i theta_y)/sqrt(2).
   pure subroutine rotation_circular_unitary(nsite,unitary)
      integer, intent(in) :: nsite
      complex(rp), intent(out) :: unitary(:, :)
      complex(rp), parameter :: inv_sqrt2=cmplx(1.0_rp/sqrt(2.0_rp),0.0_rp,rp)
      integer :: site, i0
      if (any(shape(unitary)/=[2*nsite,2*nsite])) error stop 'rotation_circular_unitary: shape mismatch'
      unitary=cmplx(0.0_rp,0.0_rp,rp)
      do site=1,nsite
         i0=2*(site-1)
         unitary(i0+1,i0+1)=inv_sqrt2
         unitary(i0+1,i0+2)=-i_unit*inv_sqrt2
         unitary(i0+2,i0+1)=inv_sqrt2
         unitary(i0+2,i0+2)= i_unit*inv_sqrt2
      end do
   end subroutine rotation_circular_unitary

   subroutine invert_complex_matrix(matrix,info)
      complex(rp), intent(inout) :: matrix(:, :)
      integer, intent(out) :: info
      complex(rp), allocatable :: work(:)
      complex(rp) :: query(1)
      integer, allocatable :: pivots(:)
      integer :: n,lwork
      n=size(matrix,1)
      if (size(matrix,2)/=n) error stop 'invert_complex_matrix: matrix must be square'
      allocate(pivots(n))
      call zgetrf(n,n,matrix,n,pivots,info)
      if (info/=0) then
         deallocate(pivots)
         return
      end if
      call zgetri(n,matrix,n,pivots,query,-1,info)
      if (info/=0) then
         deallocate(pivots)
         return
      end if
      lwork=max(n,int(real(query(1),rp)))
      allocate(work(lwork))
      call zgetri(n,matrix,n,pivots,work,lwork,info)
      deallocate(pivots,work)
   end subroutine invert_complex_matrix

   subroutine rotation_state_clear(this)
      class(rotation_state), intent(inout) :: this
      if (allocated(this%values)) deallocate(this%values)
      if (allocated(this%endpoint_values)) deallocate(this%endpoint_values)
      if (allocated(this%weights)) deallocate(this%weights)
      if (allocated(this%occupations)) deallocate(this%occupations)
      if (allocated(this%endpoint_occupations)) deallocate(this%endpoint_occupations)
      if (allocated(this%torque_band)) deallocate(this%torque_band)
      if (allocated(this%partner_band)) deallocate(this%partner_band)
      if (allocated(this%contact)) deallocate(this%contact)
      nullify(this%fixture,this%reciprocal_state)
      this%q=0.0_rp; this%nsite=0; this%ncoord=0; this%nmat=0; this%nbands=0; this%nk=0
      this%fermi=0.0_rp; this%kT=0.0_rp; this%weight_sum=0.0_rp
      this%magnetization=0.0_rp; this%berry=0.0_rp
      this%electron_count=0.0_rp; this%electron_residual=0.0_rp
      this%prepared=.false.
   end subroutine rotation_state_clear

end module lr_rotation_response_mod
