!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
!> @brief Damping vector handed to the Chebyshev reconstruction (CPU and GPU plugin) for control%cheb_kernel.
!> @details green%cheb_damping(n, damp) must leave damp unallocated for the Jackson kernel (the reconstruction then
!>          builds the Jackson factor inline, bit-identical to the default), and must fill damp(i), i = 1..n, with
!>          g_(i-1) of the selected kernel otherwise. Closed forms written here: Lorentz
!>          g_m = sinh(lambda (1 - m/N))/sinh(lambda) (Weisse et al., Rev. Mod. Phys. 78, 275, Eq. 79), none g_m = 1.
!>          N = 604 is not a power of two; math_mod builds the Lorentz argument with default-kind reals, which
!>          measured 2.0e-8 against the closed form at N = 604, hence the tolerance 1.0e-7. A wrong N (n-1, n+1),
!>          a shifted index (g_(i) in place of g_(i-1)) or a stale vector kept from the previous call changes the
!>          values by 6e-3 or more at N = 604 and lambda = 3.5.
!------------------------------------------------------------------------------
program test_cheb_damping
   use control_mod, only: control
   use green_mod, only: green
   use precision_mod, only: rp
   implicit none

   integer, parameter :: n_mom = 604
   real(rp), parameter :: lambda = 3.5_rp, tolerance = 1.0e-7_rp
   type(control), target :: ctl
   type(green) :: green_obj
   real(rp), allocatable :: damp(:)
   real(rp) :: expected
   logical :: failed
   integer :: m

   failed = .false.
   call ctl%restore_to_default()
   green_obj%control => ctl

   ctl%cheb_kernel = 'lorentz'
   ctl%cheb_kernel_lambda = lambda
   call green_obj%cheb_damping(n_mom, damp)
   call require(allocated(damp), 'lorentz: damping is allocated')
   if (allocated(damp)) then
      call require(size(damp) == n_mom, 'lorentz: size equals the number of moments')
      if (size(damp) == n_mom) then
         do m = 0, n_mom - 1
            expected = sinh(lambda*(1.0_rp - real(m, rp)/real(n_mom, rp)))/sinh(lambda)
            if (abs(damp(m + 1) - expected) > tolerance) then
               write (*, '(a,i4,2(a,es14.6))') 'FAIL: lorentz order ', m, ' actual=', damp(m + 1), ' expected=', expected
               failed = .true.
               exit
            end if
         end do
      end if
   end if

   ! The same variable is reused: the Jackson call must deallocate the previous vector, not keep it.
   ctl%cheb_kernel = 'jackson'
   call green_obj%cheb_damping(n_mom, damp)
   call require(.not. allocated(damp), 'jackson: damping stays unallocated (inline Jackson factor)')

   ctl%cheb_kernel = 'none'
   call green_obj%cheb_damping(n_mom, damp)
   call require(allocated(damp), 'none: damping is allocated')
   if (allocated(damp)) then
      call require(size(damp) == n_mom, 'none: size equals the number of moments')
      call require(all(damp == 1.0_rp), 'none: every entry is exactly 1')
   end if

   ! Another N: the vector follows the requested size, not a cached one.
   ctl%cheb_kernel = 'lorentz'
   call green_obj%cheb_damping(n_mom/2, damp)
   if (allocated(damp)) then
      call require(size(damp) == n_mom/2, 'lorentz: size follows the requested N')
      expected = sinh(lambda*(1.0_rp - real(n_mom/2 - 1, rp)/real(n_mom/2, rp)))/sinh(lambda)
      call require(abs(damp(n_mom/2) - expected) <= tolerance, 'lorentz: last entry uses the requested N')
   else
      call require(.false., 'lorentz: damping is allocated for the second N')
   end if

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description

      if (.not. condition) then
         write (*, '(a,a)') 'FAIL: ', trim(description)
         failed = .true.
      end if
   end subroutine require

end program test_cheb_damping
