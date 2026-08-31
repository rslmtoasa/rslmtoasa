! XCF-05: exercise the legacy VXC0SP -> PBEGGA origin contract.
program test_vxc0sp_gga_origin
   use precision_mod, only: rp
   use control_mod, only: control
   use logger_mod, only: g_logger
   use self_mod, only: self, vxc0sp, vxc0sp_pbe_origin_derivatives
   use xc_mod, only: xc
   use xc_radial_mod, only: radgra
   use xc_response_kernel_mod, only: xc_response_radial_projection
   implicit none

   integer, parameter :: nr = 96
   real(rp), parameter :: mesh_a = log(1.0_rp + 5.0_rp/0.02_rp)/real(nr - 1, rp)
   real(rp), parameter :: mesh_b = 0.02_rp
   real(rp), parameter :: pi = 4.0_rp*atan(1.0_rp)

   type(control) :: ctl
   type(xc) :: functional
   type(self) :: production_path
   logical :: failed

   failed = .false.
   call g_logger%init()
   call ctl%restore_to_default()
   ctl%nsp = 2
   ctl%txc = 8
   functional = xc(ctl)

   call check_case(functional, production_path, 2, 0.60_rp, 0.04_rp, 0.001_rp, &
                   0.40_rp, 0.07_rp, 0.002_rp, 'spin-polarized')
   call check_case(functional, production_path, 1, 0.85_rp, 0.05_rp, 0.001_rp, &
                   0.0_rp, 0.0_rp, 0.0_rp, 'unpolarized total density')

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine check_case(functional, production_path, nsp, n0_up, a_up, b_up, n0_down, a_down, b_down, label)
      type(xc), intent(in) :: functional
      type(self), intent(inout) :: production_path
      integer, intent(in) :: nsp
      real(rp), intent(in) :: n0_up, a_up, b_up, n0_down, a_down, b_down
      character(len=*), intent(in) :: label
      real(rp), allocatable :: rofi(:), stored_rho(:, :), t_rho(:, :), d_rho(:, :), d2_rho(:, :), potential(:, :)
      real(rp), allocatable :: rho0(:), rhoeps(:), rhomu(:)
      type(xc_response_radial_projection) :: projection
      real(rp) :: expected_up, expected_down, expected_exc, wrong_up, wrong_down, wrong_exc
      real(rp) :: n(2), nd(2), ndd(2), wrong_ndd(2), origin_nd(2), origin_ndd(2), origin_rhopp(2), density
      real(rp) :: ob4pi, z
      integer :: ir, isp

      allocate (rofi(nr), stored_rho(nr, nsp), t_rho(nr, nsp), d_rho(nr, nsp), d2_rho(nr, nsp))
      allocate (potential(nr, nsp), rho0(2), rhoeps(nsp), rhomu(nsp))
      ob4pi = 1.0_rp/(4.0_rp*pi)
      z = 26.0_rp

      do ir = 1, nr
         rofi(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
         if (nsp == 1) then
            density = n0_up + a_up*rofi(ir)**2 + b_up*rofi(ir)**4
            stored_rho(ir, 1) = 4.0_rp*pi*rofi(ir)**2*density
         else
            stored_rho(ir, 1) = 4.0_rp*pi*rofi(ir)**2*(n0_up + a_up*rofi(ir)**2 + b_up*rofi(ir)**4)
            stored_rho(ir, 2) = 4.0_rp*pi*rofi(ir)**2*(n0_down + a_down*rofi(ir)**2 + b_down*rofi(ir)**4)
         end if
      end do

      do isp = 1, nsp
         t_rho(1, isp) = ob4pi*(stored_rho(2, isp)/rofi(2)**2*rofi(3) - &
                                stored_rho(3, isp)/rofi(3)**2*rofi(2))/(rofi(3) - rofi(2))
         do ir = 2, nr
            t_rho(ir, isp) = stored_rho(ir, isp)*ob4pi/rofi(ir)**2
         end do
         call radgra(mesh_a, mesh_b, nr, rofi, t_rho(:, isp), d_rho(:, isp))
         d_rho(1, isp) = 0.0_rp
         call radgra(mesh_a, mesh_b, nr, rofi, d_rho(:, isp), d2_rho(:, isp))
      end do

      potential = 0.0_rp
      call vxc0sp(production_path, functional, z, mesh_a, mesh_b, rofi, stored_rho, nr, potential, rho0, &
                  rhoeps, rhomu, nsp, 0.0_rp, projection)

      origin_rhopp = 0.0_rp
      origin_rhopp(1) = d2_rho(1, 1)
      if (nsp == 2) origin_rhopp(2) = d2_rho(1, 2)
      call vxc0sp_pbe_origin_derivatives(nsp, origin_rhopp, origin_nd, origin_ndd)
      call require_close(origin_nd(1), 0.0_rp, trim(label)//' origin dn/dr')
      call require_close(origin_nd(2), 0.0_rp, trim(label)//' origin dn/dr second slot')
      if (nsp == 1) then
         call require_close(origin_ndd(1), 0.5_rp*d2_rho(1, 1), trim(label)//' per-spin raw d2n/dr2')
         call require_close(origin_ndd(2), 0.5_rp*d2_rho(1, 1), trim(label)//' equal-spin raw d2n/dr2')
      else
         call require_close(origin_ndd(1), d2_rho(1, 2), trim(label)//' down raw d2n/dr2')
         call require_close(origin_ndd(2), d2_rho(1, 1), trim(label)//' up raw d2n/dr2')
      end if

      nd = origin_nd
      if (nsp == 1) then
         n = [0.5_rp*t_rho(1, 1), 0.5_rp*t_rho(1, 1)]
         ndd = origin_ndd
         wrong_ndd = 3.0_rp*ndd
      else
         ! XCPOT/PBEGGA consumes the historical down/up input order.
         n = [t_rho(1, 2), t_rho(1, 1)]
         ndd = origin_ndd
         wrong_ndd = 3.0_rp*ndd
      end if
      call functional%PBEGGA(1, 1, n, nd, ndd, 0.0_rp, expected_exc, expected_down, expected_up)
      call functional%PBEGGA(1, 1, n, nd, wrong_ndd, 0.0_rp, wrong_exc, wrong_down, wrong_up)

      call require_close(potential(1, 1), expected_up, trim(label)//' VXC0SP origin up channel')
      if (nsp == 2) then
         call require_close(potential(1, 2), expected_down, trim(label)//' VXC0SP origin down channel')
      else
         call require_close(potential(1, 1), expected_down, trim(label)//' VXC0SP origin equal-spin channel')
      end if
      if (nsp == 2) then
         call require(abs(potential(1, 1) - wrong_up) > 1.0e-8_rp, trim(label)//' rejects factor-three pre-scaling')
      end if
      call require(abs(d_rho(1, 1)) <= tiny(1.0_rp), trim(label)//' enforces zero first derivative at origin')
      call require(abs(d2_rho(1, 1)) > 1.0e-6_rp, trim(label)//' retains raw second derivative at origin')
   end subroutine check_case

   subroutine require_close(actual, expected, description)
      real(rp), intent(in) :: actual, expected
      character(len=*), intent(in) :: description
      if (abs(actual - expected) > 2.0e-12_rp*max(1.0_rp, abs(expected))) then
         write (*, '(a,2(1x,es18.10))') 'FAIL: '//trim(description)//' actual expected', actual, expected
         failed = .true.
      end if
   end subroutine require_close

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description
      if (.not. condition) then
         write (*, '(a)') 'FAIL: '//trim(description)
         failed = .true.
      end if
   end subroutine require

end program test_vxc0sp_gga_origin
