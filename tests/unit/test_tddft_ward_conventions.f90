!------------------------------------------------------------------------------
! WR-00 -- pair-potential convention and Ward-identity oracle
!------------------------------------------------------------------------------
program test_tddft_ward_conventions
   use precision_mod, only: rp
   use math_mod, only: i_unit
   use response_components_mod, only: RESPONSE_MX, RESPONSE_MY, RESPONSE_PLUS, RESPONSE_MINUS
   use response_basis_mod, only: response_operator
   implicit none

   real(rp), parameter :: tol = 512.0_rp*epsilon(1.0_rp)
   logical :: failed

   failed = .false.
   call test_unhalved_circular_and_q_half()
   call test_cartesian_circular_dyson_equivalence()
   call test_uniform_scalar_kernel_reduction()
   call test_unequal_orbital_weighted_ward_identity()
   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_unhalved_circular_and_q_half()
      complex(rp), allocatable :: sigma_minus(:, :), q_minus(:, :)
      real(rp), parameter :: bxc = -0.60_rp, moment = 2.0_rp

      sigma_minus = response_operator(RESPONSE_MINUS, 1)
      q_minus = bxc/(2.0_rp*moment)*sigma_minus
      call check_complex('unhalved sigma-minus spin-flip element', sigma_minus(2, 1), cmplx(2.0_rp, 0.0_rp, rp))
      call check_complex('Q-minus contains circular half', q_minus(2, 1), cmplx(bxc/moment, 0.0_rp, rp))
      call check_true('wrong Q-minus without half is distinguishable', &
         abs(q_minus(2, 1) - cmplx(2.0_rp*bxc/moment, 0.0_rp, rp)) > 0.1_rp)
   end subroutine test_unhalved_circular_and_q_half

   subroutine test_cartesian_circular_dyson_equivalence()
      complex(rp), allocatable :: sx(:, :), sy(:, :), sp(:, :), sm(:, :)
      complex(rp) :: cartesian(2, 2), cartesian_dressed(2, 2), circular, circular_dressed
      real(rp), parameter :: b = 0.25_rp, moment = 1.0_rp, omega = 0.17_rp, eta = 0.013_rp
      real(rp) :: energies(2), occupations(2)

      ! H = -b sigma_z: the circular kernel is -b/(2M), while the Cartesian
      ! kernel is -b/M. The test fails if the circular half is reused there.
      energies = [-b, b]; occupations = [1.0_rp, 0.0_rp]
      sx = response_operator(RESPONSE_MX, 1)
      sy = response_operator(RESPONSE_MY, 1)
      sp = response_operator(RESPONSE_PLUS, 1)
      sm = response_operator(RESPONSE_MINUS, 1)
      cartesian(1, 1) = lehmann(energies, occupations, sx, sx, omega, eta)
      cartesian(1, 2) = lehmann(energies, occupations, sx, sy, omega, eta)
      cartesian(2, 1) = lehmann(energies, occupations, sy, sx, omega, eta)
      cartesian(2, 2) = lehmann(energies, occupations, sy, sy, omega, eta)
      circular = lehmann(energies, occupations, sp, sm, omega, eta)
      cartesian_dressed = dress_2x2(cartesian, -b/moment)
      circular_dressed = circular/(cmplx(1.0_rp, 0.0_rp, rp) - circular*cmplx(-b/(2.0_rp*moment), 0.0_rp, rp))
      call check_complex('bare Cartesian-to-circular response transform', reconstruct(cartesian), circular)
      call check_complex('Dyson Cartesian/circular equivalence', reconstruct(cartesian_dressed), circular_dressed)
      call check_true('Cartesian half-kernel negative control differs', &
         abs(reconstruct(dress_2x2(cartesian, -b/(2.0_rp*moment))) - circular_dressed) > 1.0e-4_rp)
   end subroutine test_cartesian_circular_dyson_equivalence

   subroutine test_uniform_scalar_kernel_reduction()
      real(rp), parameter :: b(2) = [0.40_rp, 0.40_rp], moment = 2.0_rp
      complex(rp) :: chi, xi_weighted, xi_scalar

      chi = circular_chi(b)
      xi_weighted = weighted_xi(b, moment)
      xi_scalar = chi*cmplx(-sum(b)/(2.0_rp*moment**2), 0.0_rp, rp)
      call check_complex('uniform pair-potential reduces to chi_KS K', xi_weighted, xi_scalar)
      call check_complex('uniform scalar kernel obeys Ward identity', xi_scalar, cmplx(1.0_rp, 0.0_rp, rp))
   end subroutine test_uniform_scalar_kernel_reduction

   subroutine test_unequal_orbital_weighted_ward_identity()
      real(rp), parameter :: b(2) = [0.20_rp, 0.80_rp], moment = 2.0_rp
      complex(rp) :: chi, xi_weighted, xi_site_scalar, xi_wrong_half
      real(rp) :: scalar_residual

      chi = circular_chi(b)
      xi_weighted = weighted_xi(b, moment)
      xi_site_scalar = chi*cmplx(-sum(b)/(2.0_rp*moment**2), 0.0_rp, rp)
      xi_wrong_half = 2.0_rp*xi_weighted
      scalar_residual = abs(xi_site_scalar - cmplx(1.0_rp, 0.0_rp, rp))
      call check_complex('exact weighted two-orbital Ward identity', xi_weighted, cmplx(1.0_rp, 0.0_rp, rp))
      call check_true('site-averaged scalar kernel has nonzero Goldstone residual', scalar_residual > 0.25_rp)
      call check_true('weighted and scalar Xi differ for unequal Bxc', abs(xi_weighted-xi_site_scalar) > 0.25_rp)
      call check_true('wrong circular half fails weighted Ward identity', abs(xi_wrong_half-cmplx(1.0_rp, 0.0_rp, rp)) > 0.25_rp)
   end subroutine test_unequal_orbital_weighted_ward_identity

   function circular_chi(b) result(chi)
      real(rp), intent(in) :: b(:)
      complex(rp) :: chi
      integer :: i

      chi = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, size(b)
         chi = chi + cmplx(-2.0_rp/b(i), 0.0_rp, rp)
      end do
   end function circular_chi

   function weighted_xi(b, moment) result(xi)
      real(rp), intent(in) :: b(:), moment
      complex(rp) :: xi
      integer :: i

      ! <down_i|B_i sigma_minus/(2M)|up_i>=-b_i/M because sigma_minus is unhalved.
      xi = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, size(b)
         xi = xi + cmplx((2.0_rp*(-b(i)/moment))/(-2.0_rp*b(i)), 0.0_rp, rp)
      end do
   end function weighted_xi

   function lehmann(energies, occupations, left, right, omega, eta) result(value)
      real(rp), intent(in) :: energies(:), occupations(:), omega, eta
      complex(rp), intent(in) :: left(:, :), right(:, :)
      complex(rp) :: value
      integer :: n, m

      value = cmplx(0.0_rp, 0.0_rp, rp)
      do n = 1, size(energies)
         do m = 1, size(energies)
            value = value + (occupations(n)-occupations(m))*left(n, m)*right(m, n)/ &
               cmplx(omega + energies(n)-energies(m), eta, rp)
         end do
      end do
   end function lehmann

   function dress_2x2(chi, cartesian_kernel) result(dressed)
      complex(rp), intent(in) :: chi(2, 2)
      real(rp), intent(in) :: cartesian_kernel
      complex(rp) :: dressed(2, 2), inverse(2, 2), a(2, 2), determinant

      a = -cmplx(cartesian_kernel, 0.0_rp, rp)*chi
      a(1, 1) = a(1, 1) + cmplx(1.0_rp, 0.0_rp, rp)
      a(2, 2) = a(2, 2) + cmplx(1.0_rp, 0.0_rp, rp)
      determinant = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
      inverse(1, 1) = a(2, 2)/determinant; inverse(1, 2) = -a(1, 2)/determinant
      inverse(2, 1) = -a(2, 1)/determinant; inverse(2, 2) = a(1, 1)/determinant
      dressed = matmul(inverse, chi)
   end function dress_2x2

   pure function reconstruct(chi) result(value)
      complex(rp), intent(in) :: chi(2, 2)
      complex(rp) :: value
      value = chi(1, 1) + chi(2, 2) + i_unit*chi(2, 1) - i_unit*chi(1, 2)
   end function reconstruct

   subroutine check_complex(label, actual, expected)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual, expected
      if (abs(actual-expected) > tol*max(1.0_rp, abs(expected))) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, abs(actual-expected)
         failed = .true.
      end if
   end subroutine check_complex

   subroutine check_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', label
         failed = .true.
      end if
   end subroutine check_true

end program test_tddft_ward_conventions
