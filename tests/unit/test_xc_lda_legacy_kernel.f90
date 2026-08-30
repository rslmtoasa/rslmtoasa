! Focused regression tests for the legacy spin-LDA kernel edge cases found by
! the XCR-01 fixed-density audit.
program test_xc_lda_legacy_kernel
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use xc_mod, only: xc
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   real(rp), parameter :: pi = 4.0_rp*atan(1.0_rp)
   real(rp), parameter :: radius = 1.0_rp
   real(rp), parameter :: zero_gradient(2) = 0.0_rp
   real(rp), parameter :: zero_laplacian(2) = 0.0_rp
   integer, parameter :: limit_selectors(5) = [1, 3, 4, 5, 11]
   real(rp), parameter :: limit_zeta(5) = [0.0_rp, 0.9_rp, 0.99_rp, 0.999999_rp, 1.0_rp]

   type(control) :: ctl
   logical :: failed
   integer :: iselector

   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 2
   failed = .false.

   ! Before the XCR-02 correction this point returned
   !   v_up = -3.1858480079535214 and v_down = -1.0223552813545744 Ry.
   ! Those values were a common charge-channel offset away from the exact
   ! derivative of the reported legacy energy.  The corrected values are
   ! retained here as a focused regression anchor.
   call check_bh_derivative(1, 0.5_rp, 0.99_rp, &
      -3.2013035412106694_rp, -1.0378108135640483_rp)

   do iselector = 1, size(limit_selectors)
      call check_energy_gradient(limit_selectors(iselector), 1.75_rp, 0.99_rp)
      call check_polarized_limits(limit_selectors(iselector))
   enddo

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   endif
   write (*, '(a)') 'RESULT: PASS'

contains

   pure function density_from_rs(rs) result(density)
      real(rp), intent(in) :: rs
      real(rp) :: density
      density = 3.0_rp/(4.0_rp*pi*rs**3)
   end function density_from_rs

   subroutine evaluate(functional, rho_down, rho_up, energy, v_down, v_up)
      type(xc), intent(in) :: functional
      real(rp), intent(in) :: rho_down, rho_up
      real(rp), intent(out) :: energy, v_down, v_up
      real(rp) :: rho_total

      rho_total = rho_down + rho_up
      call functional%XCPOT(rho_down, rho_up, rho_total, zero_gradient, zero_laplacian, radius, &
         v_down, v_up, energy)
   end subroutine evaluate

   subroutine check_bh_derivative(selector, rs, zeta, expected_v_up, expected_v_down)
      integer, intent(in) :: selector
      real(rp), intent(in) :: rs, zeta, expected_v_up, expected_v_down
      type(xc) :: functional
      real(rp) :: density, rho_down, rho_up, step
      real(rp) :: energy, energy_plus, energy_minus, v_down, v_up
      real(rp) :: v_down_fd, v_up_fd

      ctl%txc = selector
      functional = xc(ctl)
      density = density_from_rs(rs)
      rho_up = 0.5_rp*density*(1.0_rp + zeta)
      rho_down = 0.5_rp*density*(1.0_rp - zeta)
      step = 1.0e-6_rp*max(rho_down, rho_up)

      call evaluate(functional, rho_down, rho_up, energy, v_down, v_up)
      call evaluate(functional, rho_down + step, rho_up, energy_plus, v_down, v_up)
      call evaluate(functional, rho_down - step, rho_up, energy_minus, v_down, v_up)
      v_down_fd = ((rho_down + step + rho_up)*energy_plus - &
         (rho_down - step + rho_up)*energy_minus)/(2.0_rp*step)
      call evaluate(functional, rho_down, rho_up + step, energy_plus, v_down, v_up)
      call evaluate(functional, rho_down, rho_up - step, energy_minus, v_down, v_up)
      v_up_fd = ((rho_down + rho_up + step)*energy_plus - &
         (rho_down + rho_up - step)*energy_minus)/(2.0_rp*step)
      call evaluate(functional, rho_down, rho_up, energy, v_down, v_up)

      call require_close(v_up, expected_v_up, 'corrected BH spin-up potential', selector)
      call require_close(v_down, expected_v_down, 'corrected BH spin-down potential', selector)
      call require_close(v_up, v_up_fd, 'BH spin-up energy derivative', selector)
      call require_close(v_down, v_down_fd, 'BH spin-down energy derivative', selector)
   end subroutine check_bh_derivative

   subroutine check_energy_gradient(selector, rs, zeta)
      integer, intent(in) :: selector
      real(rp), intent(in) :: rs, zeta
      type(xc) :: functional
      real(rp) :: density, rho_down, rho_up, step
      real(rp) :: energy, energy_plus, energy_minus, v_down, v_up
      real(rp) :: v_down_fd, v_up_fd

      ctl%txc = selector
      functional = xc(ctl)
      density = density_from_rs(rs)
      rho_up = 0.5_rp*density*(1.0_rp + zeta)
      rho_down = 0.5_rp*density*(1.0_rp - zeta)
      step = 1.0e-6_rp*max(rho_down, rho_up)

      call evaluate(functional, rho_down + step, rho_up, energy_plus, v_down, v_up)
      call evaluate(functional, rho_down - step, rho_up, energy_minus, v_down, v_up)
      v_down_fd = ((rho_down + step + rho_up)*energy_plus - &
         (rho_down - step + rho_up)*energy_minus)/(2.0_rp*step)
      call evaluate(functional, rho_down, rho_up + step, energy_plus, v_down, v_up)
      call evaluate(functional, rho_down, rho_up - step, energy_minus, v_down, v_up)
      v_up_fd = ((rho_down + rho_up + step)*energy_plus - &
         (rho_down + rho_up - step)*energy_minus)/(2.0_rp*step)
      call evaluate(functional, rho_down, rho_up, energy, v_down, v_up)

      call require_close(v_up, v_up_fd, 'spin-up energy derivative', selector)
      call require_close(v_down, v_down_fd, 'spin-down energy derivative', selector)
   end subroutine check_energy_gradient

   subroutine check_polarized_limits(selector)
      integer, intent(in) :: selector
      type(xc) :: functional
      real(rp) :: density, rho_down, rho_up, energy, v_down, v_up
      integer :: izeta

      ctl%txc = selector
      functional = xc(ctl)
      density = density_from_rs(2.0_rp)
      do izeta = 1, size(limit_zeta)
         rho_up = 0.5_rp*density*(1.0_rp + limit_zeta(izeta))
         rho_down = 0.5_rp*density*(1.0_rp - limit_zeta(izeta))
         call evaluate(functional, rho_down, rho_up, energy, v_down, v_up)
         if (.not. ieee_is_finite(energy) .or. .not. ieee_is_finite(v_down) .or. &
            .not. ieee_is_finite(v_up)) then
            write (*, '(a,i0,a,f10.6)') 'FAIL: selector ', selector, &
               ' returned a non-finite polarized limit at zeta=', limit_zeta(izeta)
            failed = .true.
         endif
         if (limit_zeta(izeta) == 1.0_rp .and. abs(v_down) + abs(v_up) <= 1.0e-12_rp) then
            write (*, '(a,i0)') 'FAIL: selector discarded a finite-density fully polarized state: ', selector
            failed = .true.
         endif
      enddo
   end subroutine check_polarized_limits

   subroutine require_close(actual, expected, quantity, selector)
      real(rp), intent(in) :: actual, expected
      character(len=*), intent(in) :: quantity
      integer, intent(in) :: selector
      real(rp), parameter :: abs_tolerance = 5.0e-7_rp
      real(rp), parameter :: rel_tolerance = 5.0e-7_rp
      real(rp) :: scale

      scale = max(1.0_rp, abs(expected))
      if (abs(actual - expected) > abs_tolerance + rel_tolerance*scale) then
         write (*, '(a,i0,a,a,es16.8,a,es16.8)') 'FAIL: selector ', selector, ' ', trim(quantity)// &
            ' actual=', actual, ' expected=', expected
         failed = .true.
      endif
   end subroutine require_close

end program test_xc_lda_legacy_kernel
