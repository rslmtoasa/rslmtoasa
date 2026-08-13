! TEST-13: lean contracts for the production Liechtenstein Hubbard correction.
!
! These checks intentionally call hamiltonian%calculate_hubbard_u_potential_general
! directly.  The expected values below are only for the l=0 one-orbital case;
! the Hermiticity check exercises the production Coulomb-matrix path for l=1
! without introducing a second implementation of that algebra.
program test_lda_u_hamiltonian
   use precision_mod, only: rp
   use basis_mod, only: basis_init, nb, spin_off
   use control_mod, only: control
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   implicit none

   real(rp), parameter :: tolerance = 1.0e-11_rp
   logical :: failed

   type(hamiltonian), target :: ham
   type(lattice), target :: lat
   type(control), target :: ctl

   failed = .false.
   call g_logger%init()
   call basis_init(1)
   call setup_fixture()

   call test_zero_u_j_limit()
   call test_simple_diagonal_occupation()
   call test_collinear_hermiticity()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   subroutine setup_fixture()
      call ctl%restore_to_default()
      ctl%lmax = 1
      ctl%nsp = 2
      lat%ntype = 1
      allocate (lat%symbolic_atoms(1))
      call lat%symbolic_atoms(1)%restore_to_default()

      ham%lattice => lat
      ham%control => ctl
      allocate (ham%hubbard_u_pot(nb, nb, 1))
      ham%hubbard_u_pot = 0.0_rp
      ham%hubbard_u_potential_form = 'liechtenstein'
   end subroutine setup_fixture

   subroutine reset_fixture(u_s, j_s, u_p, j_p)
      real(rp), intent(in) :: u_s, j_s, u_p, j_p

      lat%symbolic_atoms(1)%potential%hubbard_u = 0.0_rp
      lat%symbolic_atoms(1)%potential%hubbard_j = 0.0_rp
      lat%symbolic_atoms(1)%potential%hubbard_u(1) = u_s
      lat%symbolic_atoms(1)%potential%hubbard_j(1) = j_s
      lat%symbolic_atoms(1)%potential%hubbard_u(2) = u_p
      lat%symbolic_atoms(1)%potential%hubbard_j(2) = j_p
      lat%symbolic_atoms(1)%potential%ldm = 0.0_rp
      ham%hubbard_u_pot = 0.0_rp
   end subroutine reset_fixture

   subroutine test_zero_u_j_limit()
      complex(rp) :: baseline(nb, nb), with_correction(nb, nb)
      integer :: i, j

      call reset_fixture(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
      lat%symbolic_atoms(1)%potential%ldm = 0.23_rp
      do i = 1, nb
         do j = 1, nb
            baseline(i, j) = cmplx(0.013_rp*real(i + 2*j, rp), &
               0.002_rp*real(i - j, rp), rp)
         end do
      end do

      call ham%calculate_hubbard_u_potential_general()
      with_correction = baseline + cmplx(ham%hubbard_u_pot(:, :, 1), 0.0_rp, rp)
      call require(maxval(abs(ham%hubbard_u_pot(:, :, 1))) < tolerance, &
         'U=J=0 produces zero Hubbard correction')
      call require(maxval(abs(with_correction - baseline)) < tolerance, &
         'U=J=0 leaves the Hamiltonian contribution unchanged')
   end subroutine test_zero_u_j_limit

   subroutine test_simple_diagonal_occupation()
      real(rp) :: expected_up, expected_down

      ! For one s orbital and J=0, the production expression reduces to
      ! V_up = U*(1/2 - n_up) and V_down = U*(1/2 - n_down).
      call reset_fixture(2.0_rp, 0.0_rp, 0.0_rp, 0.0_rp)
      lat%symbolic_atoms(1)%potential%ldm(1, 1, 1, 1) = 0.25_rp
      lat%symbolic_atoms(1)%potential%ldm(1, 2, 1, 1) = 0.75_rp

      ! Occupation outside the correlated shell must not leak into the s block
      ! or cause the routine to assume a d-sized Hamiltonian.
      lat%symbolic_atoms(1)%potential%ldm(2, 1, 1, 1) = 0.60_rp
      lat%symbolic_atoms(1)%potential%ldm(2, 2, 2, 2) = 0.30_rp

      call ham%calculate_hubbard_u_potential_general()
      expected_up = 2.0_rp*(0.5_rp - 0.25_rp)
      expected_down = 2.0_rp*(0.5_rp - 0.75_rp)
      call require(abs(ham%hubbard_u_pot(1, 1, 1) - expected_up) < tolerance, &
         'diagonal s-shell spin-up correction')
      call require(abs(ham%hubbard_u_pot(1 + spin_off, 1 + spin_off, 1) - expected_down) < tolerance, &
         'diagonal s-shell spin-down correction')
      call require(maxval(abs(ham%hubbard_u_pot(2:spin_off, 2:spin_off, 1))) < tolerance, &
         'uncorrelated p orbitals remain zero')
      call require(maxval(abs(ham%hubbard_u_pot(1:spin_off, spin_off + 1:nb, 1))) < tolerance, &
         'ordinary collinear path has no spin-mixing Hubbard block')
   end subroutine test_simple_diagonal_occupation

   subroutine test_collinear_hermiticity()
      real(rp) :: v_up(3, 3), v_down(3, 3)

      call reset_fixture(0.0_rp, 0.0_rp, 1.7_rp, 0.3_rp)
      call set_symmetric_p_occupation(1, [0.20_rp, 0.40_rp, 0.10_rp], &
         [0.07_rp, -0.03_rp, 0.05_rp])
      call set_symmetric_p_occupation(2, [0.30_rp, 0.15_rp, 0.25_rp], &
         [0.04_rp, -0.02_rp, 0.06_rp])

      call ham%calculate_hubbard_u_potential_general()
      v_up = ham%hubbard_u_pot(2:4, 2:4, 1)
      v_down = ham%hubbard_u_pot(2 + spin_off:4 + spin_off, 2 + spin_off:4 + spin_off, 1)
      call require(maxval(abs(v_up - transpose(v_up))) < tolerance, &
         'spin-up Hubbard potential is Hermitian')
      call require(maxval(abs(v_down - transpose(v_down))) < tolerance, &
         'spin-down Hubbard potential is Hermitian')
      call require(maxval(abs(ham%hubbard_u_pot(1:spin_off, spin_off + 1:nb, 1))) < tolerance, &
         'Hermitian collinear correction remains block diagonal in spin')
   end subroutine test_collinear_hermiticity

   subroutine set_symmetric_p_occupation(ispin, diagonal, off_diagonal)
      integer, intent(in) :: ispin
      real(rp), intent(in) :: diagonal(3), off_diagonal(3)

      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 1, 1) = diagonal(1)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 2, 2) = diagonal(2)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 3, 3) = diagonal(3)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 1, 2) = off_diagonal(1)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 2, 1) = off_diagonal(1)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 1, 3) = off_diagonal(2)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 3, 1) = off_diagonal(2)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 2, 3) = off_diagonal(3)
      lat%symbolic_atoms(1)%potential%ldm(2, ispin, 3, 2) = off_diagonal(3)
   end subroutine set_symmetric_p_occupation

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description

      if (.not. condition) then
         write (*, '(a,a)') 'FAIL: ', trim(description)
         failed = .true.
      end if
   end subroutine require

end program test_lda_u_hamiltonian
