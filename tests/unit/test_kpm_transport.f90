! TEST-15: lean production contracts for charge, spin, and orbital transport.
!
! The fixture is deliberately tiny: two sites, an sp basis, and reverse bonds.
! Current operators are built by the production Hamiltonian routines.  The
! moment check compares the production eigenpair moment kernel with a direct
! matrix Chebyshev recurrence for the same finite Hamiltonian; it does not
! implement a second Kubo-Bastin or conductivity solver.
program test_kpm_transport
   use precision_mod, only: rp
   use basis_mod, only: basis_init, nb, spin_off
   use control_mod, only: control
   use charge_mod, only: charge
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   use math_mod, only: L_z, S_z, hcpx, init_math_operators, i_unit
   use moment_kernel_mod, only: moment_onsite_block
   implicit none

   integer, parameter :: nsite = 2
   integer, parameter :: ll = 4
   real(rp), parameter :: tolerance = 2.0e-11_rp
   real(rp), parameter :: cheb_scale = 2.0_rp
   real(rp), parameter :: cheb_shift = 0.0_rp

   type(control), target :: ctl
   type(lattice), target :: lat
   type(charge), target :: charge_obj
   type(hamiltonian), target :: ham
   logical :: failed

   failed = .false.
   call g_logger%init()
   call basis_init(1)
   call init_math_operators()
   call setup_fixture()

   call test_charge_current()
   call test_spin_current()
   call test_orbital_current()
   call test_transport_moments()
   call test_s_only_orbital_symmetry()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   subroutine setup_fixture()
      integer :: i

      call ctl%restore_to_default()
      ctl%lmax = 1
      ctl%nsp = 2

      lat%ntype = 1
      lat%nrec = 1
      lat%kk = 3
      lat%nn_max = 3
      lat%alat = 1.0_rp
      allocate(lat%atlist(1), lat%iz(nsite + 1), lat%nn(1, 3), lat%cr(3, nsite + 1))
      lat%atlist = 1
      lat%iz = 1
      ! Site 1 sees site 2 at -x and site 3 at +x.  The two directed
      ! hopping blocks below are Hermitian reverses of one another.
      lat%nn(1, :) = [3, 2, 3]
      lat%cr = 0.0_rp
      lat%cr(1, 2) = -1.0_rp
      lat%cr(1, 3) = 1.0_rp

      charge_obj%lattice => lat
      ham%charge => charge_obj
      ham%lattice => lat
      ham%control => ctl
      ham%hoh = .false.
      ham%js_alpha = 'z'
      ham%jl_alpha = 'z'
      ham%v_alpha = [1.0_rp, 0.0_rp, 0.0_rp]
      ham%v_beta = ham%v_alpha

      allocate(ham%ee(nb, nb, 3, 1), ham%v_a(nb, nb, 3, 1), ham%v_b(nb, nb, 3, 1))
      allocate(ham%js_a(nb, nb, 3, 1), ham%jl_a(nb, nb, 3, 1))
      allocate(ham%jso_a(nb, nb, 3, 1), ham%jlo_a(nb, nb, 3, 1))
      allocate(ham%velocity_scale(1))
      ham%velocity_scale = 1.0_rp
      ham%ee = (0.0_rp, 0.0_rp)

      do i = 1, nb
         ham%ee(i, i, 1, 1) = cmplx(0.10_rp + 0.01_rp*real(i, rp), 0.0_rp, rp)
         ham%ee(i, i, 2, 1) = cmplx(0.03_rp + 0.002_rp*real(i, rp), 0.0_rp, rp)
         ham%ee(i, i, 3, 1) = conjg(ham%ee(i, i, 2, 1))
      end do
   end subroutine setup_fixture

   subroutine test_charge_current()
      complex(rp) :: v_forward(nb, nb)

      call ham%build_realspace_velocity_operators()
      v_forward = ham%v_a(:, :, 2, 1)
      call require(maxval(abs(ham%v_a(:, :, 3, 1) - conjg(transpose(v_forward)))) < tolerance, &
         'charge current reverse bond is Hermitian conjugate')
      call require(maxval(abs(ham%v_a(:, :, 2, 1) - ham%v_b(:, :, 2, 1))) < tolerance, &
         'equal charge-current directions share the same operator')
      call require(hermiticity_error(make_full_operator(v_forward)) < tolerance, &
         'charge current finite model is Hermitian')
   end subroutine test_charge_current

   subroutine test_spin_current()
      complex(rp) :: expected(nb, nb), v_forward(nb, nb)

      call ham%build_realspace_spin_operators()
      v_forward = ham%v_a(:, :, 2, 1)
      expected = 0.5_rp*(matmul(S_z, v_forward) + matmul(v_forward, S_z))
      call require(maxval(abs(ham%js_a(:, :, 2, 1) - expected)) < tolerance, &
         'spin current uses the production symmetrized spin-velocity product')
      call require(maxval(abs(ham%js_a(:, :, 3, 1) - conjg(transpose(ham%js_a(:, :, 2, 1))))) < tolerance, &
         'spin current reverse bond is Hermitian conjugate')
      call require(hermiticity_error(make_full_operator(ham%js_a(:, :, 2, 1))) < tolerance, &
         'spin current finite model is Hermitian')
   end subroutine test_spin_current

   subroutine test_orbital_current()
      complex(rp) :: expected(nb, nb), v_forward(nb, nb), l_operator(nb, nb)
      complex(rp) :: l_orbital(size(L_z, 1), size(L_z, 2))

      l_orbital = L_z
      call hcpx(l_orbital, 'cart2sph')
      l_operator = (0.0_rp, 0.0_rp)
      l_operator(1:norb_local(), 1:norb_local()) = l_orbital
      l_operator(norb_local() + 1:nb, norb_local() + 1:nb) = l_orbital

      call ham%build_realspace_orbital_velocity_operators()
      v_forward = ham%v_a(:, :, 2, 1)
      expected = 0.5_rp*(matmul(l_operator, v_forward) + matmul(v_forward, l_operator))
      call require(maxval(abs(ham%jl_a(:, :, 2, 1) - expected)) < tolerance, &
         'orbital current uses the production symmetrized orbital-velocity product')
      call require(maxval(abs(ham%jl_a(:, :, 3, 1) - conjg(transpose(ham%jl_a(:, :, 2, 1))))) < tolerance, &
         'orbital current reverse bond is Hermitian conjugate')
      call require(hermiticity_error(make_full_operator(ham%jl_a(:, :, 2, 1))) < tolerance, &
         'orbital current finite model is Hermitian')
      call require(maxval(abs(ham%jl_a(:, :, 2, 1))) > tolerance, &
         'orbital current is nonzero when p-orbital hopping is present')
   end subroutine test_orbital_current

   subroutine test_transport_moments()
      integer :: nmat
      real(rp), allocatable :: eigenvalues(:, :)
      complex(rp), allocatable :: eigenvectors(:, :, :), hk(:, :)
      complex(rp), allocatable :: va_k(:, :, :), vb_k(:, :, :)
      complex(rp) :: mu_exact(nb, nb, ll, ll)
      complex(rp), allocatable :: va(:, :), vb(:, :)
      character(len=*), parameter :: route_names(3) = ['charge  ', 'spin    ', 'orbital ']
      integer :: route

      nmat = nsite*nb
      allocate(eigenvalues(nmat, 1), eigenvectors(nmat, nmat, 1), hk(nmat, nmat))
      allocate(va_k(nmat, nmat, 1), vb_k(nmat, nmat, 1), va(nmat, nmat), vb(nmat, nmat))
      hk = make_hamiltonian()
      call diagonalize_hermitian(hk, eigenvalues(:, 1), eigenvectors(:, :, 1))
      vb = make_full_operator(ham%v_b(:, :, 2, 1))

      do route = 1, 3
         select case(route)
         case(1)
            va = make_full_operator(ham%v_a(:, :, 2, 1))
         case(2)
            va = make_full_operator(ham%js_a(:, :, 2, 1))
         case(3)
            va = make_full_operator(ham%jl_a(:, :, 2, 1))
         end select
         va_k(:, :, 1) = va
         vb_k(:, :, 1) = vb

         call moment_onsite_block(eigenvalues, eigenvectors, va_k, vb_k, cheb_scale, cheb_shift, &
            0, nb, ll, mu_exact)
         call compare_direct_chebyshev(trim(route_names(route)), hk, va, vb, mu_exact)
      end do
      deallocate(eigenvalues, eigenvectors, hk, va_k, vb_k, va, vb)
   end subroutine test_transport_moments

   subroutine compare_direct_chebyshev(name, hamiltonian_matrix, va, vb, exact_moments)
      character(len=*), intent(in) :: name
      complex(rp), intent(in) :: hamiltonian_matrix(:, :), va(:, :), vb(:, :)
      complex(rp), intent(in) :: exact_moments(:, :, :, :)
      complex(rp) :: tmat(size(hamiltonian_matrix, 1), size(hamiltonian_matrix, 1), ll)
      complex(rp) :: direct(size(exact_moments, 1), size(exact_moments, 2))
      real(rp) :: max_error
      integer :: n, m

      call matrix_chebyshev(hamiltonian_matrix, cheb_scale, cheb_shift, tmat)
      max_error = 0.0_rp
      do n = 1, ll
         do m = 1, ll
            direct = matmul(tmat(1:nb, :, m), matmul(va, matmul(tmat(:, :, n), vb(:, 1:nb))))
            max_error = max(max_error, maxval(abs(exact_moments(:, :, n, m) - direct)))
         end do
      end do
      write (*, '(a,a,es12.4)') 'Moment route ', trim(name)//' max_err = ', max_error
      call require(max_error < tolerance, trim(name)//' KPM/eigenpair moment agreement')
   end subroutine compare_direct_chebyshev

   subroutine test_s_only_orbital_symmetry()
      ham%ee = (0.0_rp, 0.0_rp)
      ham%ee(1, 1, 2, 1) = cmplx(0.08_rp, 0.0_rp, rp)
      ham%ee(1, 1, 3, 1) = cmplx(0.08_rp, 0.0_rp, rp)
      ham%ee(spin_off + 1, spin_off + 1, 2, 1) = cmplx(0.08_rp, 0.0_rp, rp)
      ham%ee(spin_off + 1, spin_off + 1, 3, 1) = cmplx(0.08_rp, 0.0_rp, rp)
      call ham%build_realspace_velocity_operators()
      call ham%build_realspace_orbital_velocity_operators()
      call require(maxval(abs(ham%jl_a(:, :, 2:3, 1))) < tolerance, &
         'orbital current vanishes for s-only hopping by angular-momentum symmetry')
   end subroutine test_s_only_orbital_symmetry

   function make_hamiltonian() result(hk)
      complex(rp) :: hk(nsite*nb, nsite*nb)

      hk = (0.0_rp, 0.0_rp)
      hk(1:nb, 1:nb) = ham%ee(:, :, 1, 1)
      hk(nb + 1:nsite*nb, nb + 1:nsite*nb) = ham%ee(:, :, 1, 1)
      hk(1:nb, nb + 1:nsite*nb) = ham%ee(:, :, 2, 1)
      hk(nb + 1:nsite*nb, 1:nb) = conjg(transpose(ham%ee(:, :, 2, 1)))
   end function make_hamiltonian

   function make_full_operator(block) result(operator)
      complex(rp), intent(in) :: block(:, :)
      complex(rp) :: operator(nsite*nb, nsite*nb)

      operator = (0.0_rp, 0.0_rp)
      operator(1:nb, nb + 1:nsite*nb) = block
      operator(nb + 1:nsite*nb, 1:nb) = conjg(transpose(block))
   end function make_full_operator

   subroutine matrix_chebyshev(matrix, a, b, tmat)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(in) :: a, b
      complex(rp), intent(out) :: tmat(:, :, :)
      complex(rp), allocatable :: scaled(:, :)
      integer :: p, i, n

      n = size(matrix, 1)
      allocate(scaled(n, n))
      scaled = matrix/cmplx(a, 0.0_rp, rp)
      do i = 1, n
         scaled(i, i) = scaled(i, i) - cmplx(b/a, 0.0_rp, rp)
      end do
      tmat = (0.0_rp, 0.0_rp)
      do i = 1, n
         tmat(i, i, 1) = (1.0_rp, 0.0_rp)
      end do
      if (size(tmat, 3) >= 2) tmat(:, :, 2) = scaled
      do p = 3, size(tmat, 3)
         tmat(:, :, p) = 2.0_rp*matmul(scaled, tmat(:, :, p - 1)) - tmat(:, :, p - 2)
      end do
      deallocate(scaled)
   end subroutine matrix_chebyshev

   subroutine diagonalize_hermitian(matrix, eigenvalues, eigenvectors)
      complex(rp), intent(in) :: matrix(:, :)
      real(rp), intent(out) :: eigenvalues(:)
      complex(rp), intent(out) :: eigenvectors(:, :)
      complex(rp) :: work_query(1)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      integer :: n, lwork, info
      external :: zheev

      n = size(matrix, 1)
      eigenvectors = matrix
      allocate(rwork(max(1, 3*n - 2)))
      call zheev('V', 'U', n, eigenvectors, n, eigenvalues, work_query, -1, rwork, info)
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      call zheev('V', 'U', n, eigenvectors, n, eigenvalues, work, lwork, rwork, info)
      call require(info == 0, 'tiny transport Hamiltonian diagonalization succeeds')
      deallocate(work, rwork)
   end subroutine diagonalize_hermitian

   real(rp) function hermiticity_error(matrix)
      complex(rp), intent(in) :: matrix(:, :)

      hermiticity_error = maxval(abs(matrix - conjg(transpose(matrix))))
   end function hermiticity_error

   integer function norb_local()
      norb_local = size(L_z, 1)
   end function norb_local

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description

      if (.not. condition) then
         write (*, '(a,a)') 'FAIL: ', trim(description)
         failed = .true.
      end if
   end subroutine require

end program test_kpm_transport
