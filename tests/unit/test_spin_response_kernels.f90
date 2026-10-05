!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_spin_response_kernels
!
!> @brief spin_response C0 kernels against the developer's oracles A1-A3.
!> @details Usage: test_spin_response_kernels <oracle directory>. Inputs, expected
!>          values and tolerances are read from <dir>/a1_two_level.nml,
!>          a2_single_band.nml, a3_doubled_cell.nml and tolerances.nml; nothing is
!>          written there. Relative error: max|x - ref|/max|ref| on complex moduli
!>          for arrays, |x - ref|/|ref| for scalars.
!>
!>          Band layout of A1 and A2: band 2(mu-1)+s, s = 1 up, 2 down, one site,
!>          amplitude 1 at (mu, s) and 0 elsewhere. The oracle does not give the up
!>          bands at k+q; they are filled with the k values, which cannot matter
!>          because their T vanishes.
!------------------------------------------------------------------------------
program test_spin_response_kernels
   use precision_mod, only: rp
   use spin_response_mod, only: accumulate_chi0, u_juelich, solve_dyson, &
                                spectral_trace, pole_estimates
   use spin_response_mod, only: u_mills_kernel => u_mills
   implicit none

   real(rp) :: a1_rel, a2_chi_rel, a2_identity_rel, static_eta0_rel, a3_rel, &
               window_moment_rel, dyson_u0_rel, pole_peak_rel, pole_crossing_rel
   real(rp) :: electron_count_abs, eigen_contract_abs, moment_rel   ! read, used by the C1/C2 tests
   namelist /spin_response_tolerances/ a1_rel, a2_chi_rel, a2_identity_rel, static_eta0_rel, &
      a3_rel, window_moment_rel, dyson_u0_rel, pole_peak_rel, pole_crossing_rel, &
      electron_count_abs, eigen_contract_abs, moment_rel

   character(len=512) :: dir
   logical :: failed
   integer :: unit_tol

   failed = .false.
   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_spin_response_kernels <oracle directory>'
   open (newunit=unit_tol, file=trim(dir)//'/tolerances.nml', status='old', action='read')
   read (unit_tol, nml=spin_response_tolerances)
   close (unit_tol)

   write (*, '(a)') 'test | measured | tolerance'
   call test_a1()
   call test_a2()
   call test_a3()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   subroutine report(label, measured, tol)
      character(*), intent(in) :: label
      real(rp), intent(in) :: measured, tol

      logical :: ok

      ok = measured <= tol
      write (*, '(a,t58,es10.3,2x,es10.3,2x,a)') label, measured, tol, merge('ok  ', 'FAIL', ok)
      if (.not. ok) failed = .true.
   end subroutine report

   subroutine check(label, ok)
      character(*), intent(in) :: label
      logical, intent(in) :: ok

      write (*, '(a,t58,a)') label, merge('ok  ', 'FAIL', ok)
      if (.not. ok) failed = .true.
   end subroutine check

   function rel_arr(x, ref) result(r)
      complex(rp), intent(in) :: x(:), ref(:)
      real(rp) :: r

      r = maxval(abs(x - ref))/maxval(abs(ref))
   end function rel_arr

   function rel_scal(x, ref) result(r)
      complex(rp), intent(in) :: x, ref
      real(rp) :: r

      r = abs(x - ref)/abs(ref)
   end function rel_scal

   !> Integer value of `key = n` in an oracle file (needed to allocate before the namelist read).
   function scan_int(file, key) result(n)
      character(*), intent(in) :: file, key
      integer :: n

      character(len=512) :: line, rest
      integer :: u, ios, p

      open (newunit=u, file=file, status='old', action='read')
      do
         read (u, '(a)', iostat=ios) line
         if (ios /= 0) error stop 'scan_int: key not found: '//key
         line = adjustl(line)
         p = len_trim(key)
         if (line(1:p) /= key) cycle
         rest = adjustl(line(p + 1:))
         if (rest(1:1) /= '=') cycle
         read (rest(2:), *) n
         exit
      end do
      close (u)
   end function scan_int

   !> Energies and occupations of the (mu, s) bands, shape (nk, 2 norb).
   subroutine two_spin_bands(norb, e_up, e_dn, f_up, f_dn, e, f)
      integer, intent(in) :: norb
      real(rp), intent(in) :: e_up(:), e_dn(:), f_up(:), f_dn(:)
      real(rp), allocatable, intent(out) :: e(:, :), f(:, :)

      integer :: mu

      allocate (e(size(e_up), 2*norb), f(size(e_up), 2*norb))
      do mu = 1, norb
         e(:, 2*mu - 1) = e_up
         e(:, 2*mu) = e_dn
         f(:, 2*mu - 1) = f_up
         f(:, 2*mu) = f_dn
      end do
   end subroutine two_spin_bands

   !> Amplitude 1 at (band (mu, s), orbital mu, site 1, spin s), shape (nk, 2 norb, norb, 1, 2).
   subroutine band_amplitudes(nk, norb, a)
      integer, intent(in) :: nk, norb
      complex(rp), allocatable, intent(out) :: a(:, :, :, :, :)

      integer :: mu, s

      allocate (a(nk, 2*norb, norb, 1, 2))
      a = (0.0_rp, 0.0_rp)
      do mu = 1, norb
         do s = 1, 2
            a(:, 2*(mu - 1) + s, mu, 1, s) = (1.0_rp, 0.0_rp)
         end do
      end do
   end subroutine band_amplitudes

   !> Static value (omega = 0, eta = 0) of the one-site chi0 for the given bands.
   function static_chi0(kt, w, e_k, f_k, a, e_kq, f_kq) result(c)
      real(rp), intent(in) :: kt, w(:), e_k(:, :), f_k(:, :), e_kq(:, :), f_kq(:, :)
      complex(rp), intent(in) :: a(:, :, :, :, :)
      complex(rp) :: c

      complex(rp) :: chi0(1, 1, 1)

      chi0 = (0.0_rp, 0.0_rp)
      call accumulate_chi0(kt, w, e_k, f_k, a, e_kq, f_kq, a, [0.0_rp], 0.0_rp, chi0)
      c = chi0(1, 1, 1)
   end function static_chi0

   !> Real part of a static chi0, after requiring max|Im| <= a2_identity_rel * |Re|.
   function static_real(label, z) result(r)
      character(*), intent(in) :: label
      complex(rp), intent(in) :: z
      real(rp) :: r

      call check(label//' |Im| <= a2_identity_rel*|Re|', abs(aimag(z)) <= a2_identity_rel*abs(real(z)))
      r = real(z)
   end function static_real

   !> Trapezoid integral of y over x.
   function trapezoid(x, y) result(s)
      real(rp), intent(in) :: x(:), y(:)
      real(rp) :: s

      s = sum(0.5_rp*(y(1:size(y) - 1) + y(2:))*(x(2:) - x(1:size(x) - 1)))
   end function trapezoid

   !---------------------------------------------------------------------------
   subroutine test_a1()
      character(len=*), parameter :: name = '/a1_two_level.nml'
      integer :: norb, nomega, neta_static, iu, ip
      real(rp) :: delta, eta, e_up, e_dn, f_up, f_dn, u_test, mills_factor
      character(len=256) :: amplitude_rule
      real(rp), allocatable :: omega(:), eta_static(:)
      namelist /a1_inputs/ norb, delta, eta, amplitude_rule, e_up, e_dn, f_up, f_dn, u_test, &
         mills_factor, nomega, omega, neta_static, eta_static
      real(rp) :: moment, u_goldstone, chi0_static_eta0, pole, pole_weight, u_mills_variant, &
                  residual_mills_variant, gap_mills_variant
      real(rp), allocatable :: chi0_static_ladder_re(:), chi0_static_ladder_im(:), chi0_re(:), &
                               chi0_im(:), chi_re(:), chi_im(:), spectral(:)
      namelist /a1_expected/ moment, u_goldstone, chi0_static_eta0, chi0_static_ladder_re, &
         chi0_static_ladder_im, pole, pole_weight, u_mills_variant, residual_mills_variant, &
         gap_mills_variant, chi0_re, chi0_im, chi_re, chi_im, spectral

      real(rp), allocatable :: e(:, :), f(:, :), trl(:)
      complex(rp), allocatable :: a(:, :, :, :, :), chi0(:, :, :), chi(:, :, :)
      complex(rp) :: s1, s2
      real(rp) :: chi0s(1, 1), um(1), uj(1), res, peak, crossing, h, d, u_s
      real(rp), allocatable :: omega_s(:)
      complex(rp), allocatable :: chi0_s(:, :, :), chi_s(:, :, :)
      logical :: has_peak, has_crossing

      nomega = scan_int(trim(dir)//name, 'nomega')
      neta_static = scan_int(trim(dir)//name, 'neta_static')
      allocate (omega(nomega), eta_static(neta_static), chi0_static_ladder_re(neta_static), &
                chi0_static_ladder_im(neta_static), chi0_re(nomega), chi0_im(nomega), &
                chi_re(nomega), chi_im(nomega), spectral(nomega))
      open (newunit=iu, file=trim(dir)//name, status='old', action='read')
      read (iu, nml=a1_inputs)
      read (iu, nml=a1_expected)
      close (iu)

      call two_spin_bands(norb, [e_up], [e_dn], [f_up], [f_dn], e, f)
      call band_amplitudes(1, norb, a)
      allocate (chi0(1, 1, nomega), chi(1, 1, nomega))
      chi0 = (0.0_rp, 0.0_rp)
      call accumulate_chi0(1.0_rp, [1.0_rp], e, f, a, e, f, a, omega, eta, chi0)

      ! catches: wrong vertex factor, extra prefactor, sign of eta
      call report('A1 chi0 (601 omega)', rel_arr(chi0(1, 1, :), cmplx(chi0_re, chi0_im, rp)), a1_rel)

      ! catches: Dyson built as I - chi0 U or with the wrong side
      call solve_dyson(chi0, [u_test], chi)
      call report('A1 chi', rel_arr(chi(1, 1, :), cmplx(chi_re, chi_im, rp)), a1_rel)

      ! catches: sign or pi in the spectral trace
      trl = spectral_trace(chi)
      call report('A1 spectral', rel_arr(cmplx(trl, 0.0_rp, rp), cmplx(spectral, 0.0_rp, rp)), a1_rel)

      ! catches: finite eta in the static value; kt leaking into occupied/empty pairs
      s1 = static_chi0(1.0_rp, [1.0_rp], e, f, a, e, f)
      s2 = static_chi0(1.0e-3_rp, [1.0_rp], e, f, a, e, f)
      call check('A1 static: kt = 1 and kt = 1e-3 bitwise equal', s1 == s2)
      call report('A1 static chi0 at eta = 0', rel_scal(s1, cmplx(chi0_static_eta0, 0.0_rp, rp)), static_eta0_rel)

      ! catches: wrong sign or missing M in the Goldstone solve
      chi0s(1, 1) = static_real('A1 static', s1)
      call u_juelich(chi0s, [moment], uj)
      call report('A1 u_goldstone', rel_scal(cmplx(uj(1), 0.0_rp, rp), cmplx(u_goldstone, 0.0_rp, rp)), a1_rel)

      ! catches: wrong sign or normalisation in u_mills, residual not 1 + chi0(0,0) U
      um = u_mills_kernel([e_up], [e_up + mills_factor*delta], [moment])
      call report('A1 u_mills variant', rel_scal(cmplx(um(1), 0.0_rp, rp), cmplx(u_mills_variant, 0.0_rp, rp)), a1_rel)
      res = 1.0_rp + chi0s(1, 1)*um(1)
      call report('A1 Mills variant residual', rel_scal(cmplx(res, 0.0_rp, rp), cmplx(residual_mills_variant, 0.0_rp, rp)), a1_rel)
      call report('A1 Mills variant gap', rel_scal(cmplx(res*delta, 0.0_rp, rp), cmplx(gap_mills_variant, 0.0_rp, rp)), a1_rel)

      ! catches: bad parabola; Re part, sign or omega < 0 start in the crossing
      call pole_estimates(omega, trl, chi0(1, 1, :), u_test, peak, crossing, has_peak, has_crossing)
      call check('A1 pole: peak and crossing found', has_peak .and. has_crossing)
      call report('A1 pole from peak of tr L', abs(peak - pole)/abs(pole), pole_peak_rel)
      call report('A1 pole from Re(1 + chi0 U) crossing', abs(crossing - pole)/abs(pole), pole_crossing_rel)

      ! catches: a maximum on a grid end reported as a pole
      ip = minloc(abs(omega - pole), 1)
      call pole_estimates(omega(:ip - 3), trl(:ip - 3), chi0(1, 1, :ip - 3), u_test, peak, crossing, has_peak, has_crossing)
      call check('A1 boundary: maximum on the last point gives n/a', .not. has_peak)
      call pole_estimates(omega(ip + 3:), trl(ip + 3:), chi0(1, 1, ip + 3:), u_test, peak, crossing, has_peak, has_crossing)
      call check('A1 boundary: maximum on the first point gives n/a', .not. has_peak)

      ! catches: a crossing in the interval that straddles omega = 0 skipped. The grid is shifted by half
      ! a step so that 0 is not a grid point; the pole d*delta = h/4 lies in the straddling interval.
      ! eta is a tenth of the file's value so that its shift of the zero stays far below pole_crossing_rel.
      h = omega(2) - omega(1)
      d = 0.25_rp*h/delta
      u_s = (1.0_rp - d)*delta/moment
      omega_s = omega + 0.5_rp*h
      allocate (chi0_s(1, 1, nomega), chi_s(1, 1, nomega))
      chi0_s = (0.0_rp, 0.0_rp)
      call accumulate_chi0(1.0_rp, [1.0_rp], e, f, a, e, f, a, omega_s, 0.1_rp*eta, chi0_s)
      call solve_dyson(chi0_s, [u_s], chi_s)
      call pole_estimates(omega_s, spectral_trace(chi_s), chi0_s(1, 1, :), u_s, peak, crossing, has_peak, has_crossing)
      call check('A1 straddle: crossing found', has_crossing)
      call report('A1 straddle: crossing vs d*delta', abs(crossing - d*delta)/(d*delta), pole_crossing_rel)

      ! E3. catches: Dyson plumbing error at U = 0
      call solve_dyson(chi0, [0.0_rp], chi)
      call report('E3 A1 chi = chi0 at U = 0', rel_arr(chi(1, 1, :), chi0(1, 1, :)), dyson_u0_rel)
   end subroutine test_a1

   !---------------------------------------------------------------------------
   subroutine test_a2()
      character(len=*), parameter :: name = '/a2_single_band.nml'
      integer :: norb, nk, nomega, neta_static, nq, iu, iq
      real(rp) :: t, filling, u, kt, eta
      character(len=256) :: amplitude_rule, array_dims
      real(rp), allocatable :: omega(:), eta_static(:), k(:), w(:), e_up(:), e_dn(:), f_up(:), f_dn(:)
      namelist /a2_inputs/ norb, nk, t, filling, amplitude_rule, array_dims, u, kt, eta, nomega, omega, &
         neta_static, eta_static, k, w, e_up, e_dn, f_up, f_dn, nq
      real(rp) :: ef, moment, splitting, u_goldstone_exact, u_mills
      namelist /a2_ground_expected/ ef, moment, splitting, u_goldstone_exact, u_mills
      real(rp) :: q, chi0_static_eta0, magnon_eta0, continuum_edge, chi0_window_moment, chi0_full_moment
      real(rp), allocatable :: e_dn_kq(:), f_dn_kq(:), chi0_re(:), chi0_im(:), chi_re(:), chi_im(:), &
                               chi0_static_ladder_re(:), chi0_static_ladder_im(:)
      namelist /a2_q/ iq, q, e_dn_kq, f_dn_kq, chi0_re, chi0_im, chi_re, chi_im, chi0_static_eta0, &
         chi0_static_ladder_re, chi0_static_ladder_im, magnon_eta0, continuum_edge, chi0_window_moment, &
         chi0_full_moment

      real(rp), allocatable :: e(:, :), f(:, :), e_kq(:, :), f_kq(:, :), trl(:)
      complex(rp), allocatable :: a(:, :, :, :, :), chi0(:, :, :), chi(:, :, :)
      complex(rp) :: s
      real(rp) :: w_chi0, w_chi, w_static, w_window, chi0s(1, 1), uj(1), um(1)

      nk = scan_int(trim(dir)//name, 'nk')
      nomega = scan_int(trim(dir)//name, 'nomega')
      neta_static = scan_int(trim(dir)//name, 'neta_static')
      nq = scan_int(trim(dir)//name, 'nq')
      allocate (omega(nomega), eta_static(neta_static), k(nk), w(nk), e_up(nk), e_dn(nk), f_up(nk), f_dn(nk), &
                e_dn_kq(nk), f_dn_kq(nk), chi0_re(nomega), chi0_im(nomega), chi_re(nomega), chi_im(nomega), &
                chi0_static_ladder_re(neta_static), chi0_static_ladder_im(neta_static))
      open (newunit=iu, file=trim(dir)//name, status='old', action='read')
      read (iu, nml=a2_inputs)
      read (iu, nml=a2_ground_expected)

      call two_spin_bands(norb, e_up, e_dn, f_up, f_dn, e, f)
      call band_amplitudes(nk, norb, a)
      allocate (chi0(1, 1, nomega), chi(1, 1, nomega))
      w_chi0 = 0.0_rp; w_chi = 0.0_rp; w_static = 0.0_rp; w_window = 0.0_rp

      do iq = 1, nq
         read (iu, nml=a2_q)
         call two_spin_bands(norb, e_up, e_dn_kq, f_up, f_dn_kq, e_kq, f_kq)
         chi0 = (0.0_rp, 0.0_rp)
         call accumulate_chi0(kt, w, e, f, a, e_kq, f_kq, a, omega, eta, chi0)
         call solve_dyson(chi0, [u], chi)
         w_chi0 = max(w_chi0, rel_arr(chi0(1, 1, :), cmplx(chi0_re, chi0_im, rp)))
         w_chi = max(w_chi, rel_arr(chi(1, 1, :), cmplx(chi_re, chi_im, rp)))

         s = static_chi0(kt, w, e, f, a, e_kq, f_kq)
         w_static = max(w_static, rel_scal(s, cmplx(chi0_static_eta0, 0.0_rp, rp)))

         if (q == 0.0_rp) then
            ! catches: static value that is not -1/U at q = 0 (Goldstone)
            chi0s(1, 1) = static_real('A2 static q = 0', s)
            call u_juelich(chi0s, [moment], uj)
            call report('A2 Goldstone u_juelich (q = 0)', &
                        rel_scal(cmplx(uj(1), 0.0_rp, rp), cmplx(u_goldstone_exact, 0.0_rp, rp)), a2_identity_rel)
         end if

         trl = spectral_trace(chi0)
         w_window = max(w_window, abs(trapezoid(omega, trl) - chi0_window_moment)/abs(chi0_window_moment))
      end do
      close (iu)

      ! catches: f_n/f_m swapped, wrong k+q energies or occupations, weights, orbital sum
      call report('A2 chi0, worst of all q', w_chi0, a2_chi_rel)
      ! catches: Dyson error that a one-site A1 grid does not expose
      call report('A2 chi, worst of all q', w_chi, a2_chi_rel)
      ! catches: divided difference wrong or finite eta in the static value
      call report('A2 static chi0 at eta = 0, worst of all q', w_static, static_eta0_rel)
      ! catches: normalisation or sign error in chi0 seen through tr L (E1)
      call report('A2 E1 window moment, worst of all q', w_window, window_moment_rel)

      ! catches: wrong sign or normalisation in u_mills
      um = u_mills_kernel([0.0_rp], [splitting], [moment])
      call report('A2 u_mills', rel_scal(cmplx(um(1), 0.0_rp, rp), cmplx(u_mills, 0.0_rp, rp)), a2_identity_rel)
   end subroutine test_a2

   !---------------------------------------------------------------------------
   subroutine test_a3()
      character(len=*), parameter :: name = '/a3_doubled_cell.nml'
      integer :: norb, nsite, nband, nk, nomega, nq, iu, iq
      real(rp) :: u, kt, eta, ef
      character(len=256) :: array_dims, gauge
      real(rp), allocatable :: omega(:), k(:), w(:), e(:, :), f(:, :), a_re(:, :, :, :, :), a_im(:, :, :, :, :)
      namelist /a3_inputs/ norb, nsite, nband, nk, array_dims, gauge, u, kt, eta, ef, nomega, omega, k, w, &
         e, f, a_re, a_im, nq
      real(rp) :: q, g_super
      real(rp), allocatable :: e_kq(:, :), f_kq(:, :), a_kq_re(:, :, :, :, :), a_kq_im(:, :, :, :, :), &
                               chi0_prim_q_re(:), chi0_prim_q_im(:), chi0_prim_qg_re(:), chi0_prim_qg_im(:), &
                               tr_chi_re(:), tr_chi_im(:), tr_spectral(:), det_re(:), det_im(:)
      namelist /a3_q/ iq, q, g_super, e_kq, f_kq, a_kq_re, a_kq_im, chi0_prim_q_re, chi0_prim_q_im, &
         chi0_prim_qg_re, chi0_prim_qg_im, tr_chi_re, tr_chi_im, tr_spectral, det_re, det_im

      complex(rp), allocatable :: chi0(:, :, :), chi(:, :, :), det(:)
      real(rp) :: w_tr, w_det, w_spec, w_e3
      integer :: iw

      nk = scan_int(trim(dir)//name, 'nk')
      nband = scan_int(trim(dir)//name, 'nband')
      norb = scan_int(trim(dir)//name, 'norb')
      nsite = scan_int(trim(dir)//name, 'nsite')
      nomega = scan_int(trim(dir)//name, 'nomega')
      nq = scan_int(trim(dir)//name, 'nq')
      allocate (omega(nomega), k(nk), w(nk), e(nk, nband), f(nk, nband), a_re(nk, nband, norb, nsite, 2), &
                a_im(nk, nband, norb, nsite, 2), e_kq(nk, nband), f_kq(nk, nband), &
                a_kq_re(nk, nband, norb, nsite, 2), a_kq_im(nk, nband, norb, nsite, 2), &
                chi0_prim_q_re(nomega), chi0_prim_q_im(nomega), chi0_prim_qg_re(nomega), chi0_prim_qg_im(nomega), &
                tr_chi_re(nomega), tr_chi_im(nomega), tr_spectral(nomega), det_re(nomega), det_im(nomega))
      open (newunit=iu, file=trim(dir)//name, status='old', action='read')
      read (iu, nml=a3_inputs)
      allocate (chi0(nsite, nsite, nomega), chi(nsite, nsite, nomega), det(nomega))
      w_tr = 0.0_rp; w_det = 0.0_rp; w_spec = 0.0_rp; w_e3 = 0.0_rp

      do iq = 1, nq
         read (iu, nml=a3_q)
         chi0 = (0.0_rp, 0.0_rp)
         call accumulate_chi0(kt, w, e, f, cmplx(a_re, a_im, rp), e_kq, f_kq, cmplx(a_kq_re, a_kq_im, rp), &
                              omega, eta, chi0)
         call solve_dyson(chi0, [(u, iw=1, nsite)], chi)
         w_tr = max(w_tr, rel_arr(chi(1, 1, :) + chi(2, 2, :), cmplx(tr_chi_re, tr_chi_im, rp)))
         w_spec = max(w_spec, rel_arr(cmplx(spectral_trace(chi), 0.0_rp, rp), cmplx(tr_spectral, 0.0_rp, rp)))
         ! det(I + chi0 U) of the 2x2 cell matrix, U = u on both sites
         det = (1.0_rp + u*chi0(1, 1, :))*(1.0_rp + u*chi0(2, 2, :)) - u*chi0(1, 2, :)*u*chi0(2, 1, :)
         w_det = max(w_det, rel_arr(det, cmplx(det_re, det_im, rp)))

         call solve_dyson(chi0, [(0.0_rp, iw=1, nsite)], chi)
         w_e3 = max(w_e3, rel_arr(reshape(chi, [size(chi)]), reshape(chi0, [size(chi0)])))
      end do
      close (iu)

      ! catches: gauge-dependent vertex (dropped conj, |T|^2 per band), site ordering of T_i T_j^*
      call report('A3 tr chi, worst of all Q', w_tr, a3_rel)
      call report('A3 det(I + chi0 U), worst of all Q', w_det, a3_rel)
      ! catches: spectral trace that only handles one site
      call report('A3 tr L, worst of all Q', w_spec, a3_rel)
      ! catches: Dyson plumbing error at U = 0, 2x2 case
      call report('E3 A3 chi = chi0 at U = 0, worst of all Q', w_e3, dyson_u0_rel)
   end subroutine test_a3

end program test_spin_response_kernels
