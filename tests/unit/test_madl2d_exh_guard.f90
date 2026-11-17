!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_madl2d_exh_guard
!
!> @brief Pins the madl2d overflow-guard branch feeding the l=1 dipole matrix
!>        DSZ (B7.0 / B7 §2.5), dependency-free.
!>
!> @details madl2d evaluates, per reciprocal-lattice vector G, the pair
!>
!>            EXF = EXPP*ERFCP + EXPM*ERFCM
!>            EXH = EXPP*ERFCP - EXPM*ERFCM
!>
!>          with EXPP = exp(+DGI*QPPZ), EXPM = 1/EXPP,
!>          ERFCP = erfc(BETA + AQPPZ), ERFCM = erfc(BETA - AQPPZ),
!>          BETA = DGI/(2*lambda), AQPPZ = lambda*QPPZ.
!>
!>          For large BETAP = BETA + AQPPZ -- i.e. large interlayer separation
!>          QPPZ and/or large |G| -- ERFCP underflows to exactly zero while EXPP
!>          overflows. The guard avoids forming EXPP at all. Since EXPP*ERFCP is
!>          then identically zero, the guard must reproduce
!>
!>            EXF = EXPM*ERFCM,   EXH = -EXPM*ERFCM = -EXF.
!>
!>          The historical code wrote `exh = -exh`, reusing EXH from the PREVIOUS
!>          G iteration: a stale, possibly O(1) value substituted where an
!>          exponentially small one belongs (and undefined if the branch trips on
!>          the first G of a pair). EXH feeds only DMG -> SUM1G -> DMDL -> DSZ,
!>          so it corrupts the DIPOLE matrix specifically, precisely in the
!>          thick-slab / large-separation regime interface clusters live in.
!>
!>          Pins, all against an unguarded high-precision evaluation:
!>          1. Guard identity per G: in the ERFCP == 0 regime the guarded EXF/EXH
!>             equal the analytic limits (EXF = EXPM*ERFCM, EXH = -EXF), and the
!>             stale-EXH form does not.
!>          2. Continuity across the guard boundary: sweeping QPPZ through the
!>             threshold, EXH is smooth and strictly sign-definite (negative).
!>          3. DSZ accumulation: the full SUM1G -> DMDL reciprocal-space sum for a
!>             representative slab geometry, guarded vs unguarded quad-ish
!>             reference, agrees to tolerance; the stale form does not.
!>
!>          Exits non-zero (error stop) on any failure so ctest registers a fail.
!------------------------------------------------------------------------------
program test_madl2d_exh_guard
   use precision_mod, only: rp
   implicit none

   real(rp), parameter :: tol = 1.0e-12_rp
   logical :: failed

   failed = .false.

   call test_guard_identity()
   call test_guard_continuity()
   call test_dsz_accumulation()
   call test_source_has_fix()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !> The guarded branch as it now stands in madl2d.
   subroutine exf_exh_guarded(dgi, qppz, alamda, exf, exh)
      real(rp), intent(in) :: dgi, qppz, alamda
      real(rp), intent(out) :: exf, exh
      real(rp) :: beta, aqppz, betap, betam, erfcp, erfcm, expp, expm

      beta = dgi/(2.0_rp*alamda)
      aqppz = alamda*qppz
      betap = beta + aqppz
      betam = beta - aqppz
      erfcp = erfc(betap)
      erfcm = erfc(betam)

      if (erfcp == 0.0_rp) then
         expm = exp(-dgi*qppz)
         exf = expm*erfcm
         exh = -exf
      else
         expp = exp(dgi*qppz)
         expm = 1.0_rp/expp
         exf = expp*erfcp + expm*erfcm
         exh = expp*erfcp - expm*erfcm
      end if
   end subroutine exf_exh_guarded

   !> The historical buggy form, kept so the test can demonstrate it FAILS.
   !> exh_prev emulates the value left over from the previous G iteration.
   subroutine exf_exh_stale(dgi, qppz, alamda, exh_prev, exf, exh)
      real(rp), intent(in) :: dgi, qppz, alamda, exh_prev
      real(rp), intent(out) :: exf, exh
      real(rp) :: beta, aqppz, betap, betam, erfcp, erfcm, expp, expm

      beta = dgi/(2.0_rp*alamda)
      aqppz = alamda*qppz
      betap = beta + aqppz
      betam = beta - aqppz
      erfcp = erfc(betap)
      erfcm = erfc(betam)

      if (erfcp == 0.0_rp) then
         expm = exp(-dgi*qppz)
         exf = expm*erfcm
         exh = -exh_prev
      else
         expp = exp(dgi*qppz)
         expm = 1.0_rp/expp
         exf = expp*erfcp + expm*erfcm
         exh = expp*erfcp - expm*erfcm
      end if
   end subroutine exf_exh_stale

   !> Unguarded reference. EXPP*ERFCP is evaluated as a single exponentially
   !> scaled product exp(dgi*qppz - betap^2)*erfcx(betap), which never forms the
   !> overflowing factor on its own and so stays accurate where the naive
   !> expression cannot be evaluated at all.
   subroutine exf_exh_reference(dgi, qppz, alamda, exf, exh)
      real(rp), intent(in) :: dgi, qppz, alamda
      real(rp), intent(out) :: exf, exh
      real(rp) :: beta, aqppz, betap, betam, term_p, term_m, arg

      beta = dgi/(2.0_rp*alamda)
      aqppz = alamda*qppz
      betap = beta + aqppz
      betam = beta - aqppz

      ! EXPP*ERFCP = exp(dgi*qppz) * erfc(betap)
      !            = exp(dgi*qppz - betap^2) * erfcx(betap),  erfcx = exp(x^2)erfc(x)
      arg = dgi*qppz - betap*betap
      if (betap > 0.0_rp) then
         if (arg < -700.0_rp) then
            term_p = 0.0_rp
         else
            term_p = exp(arg)*erfcx_pos(betap)
         end if
      else
         term_p = exp(dgi*qppz)*erfc(betap)
      end if

      term_m = exp(-dgi*qppz)*erfc(betam)

      exf = term_p + term_m
      exh = term_p - term_m
   end subroutine exf_exh_reference

   !> erfcx(x) = exp(x^2) * erfc(x) for x > 0, via a continued-fraction /
   !> asymptotic split. Accurate to ~1e-14 relative over the range used here.
   pure function erfcx_pos(x) result(r)
      real(rp), intent(in) :: x
      real(rp) :: r
      real(rp), parameter :: one_over_sqrt_pi = 0.5641895835477563_rp
      real(rp) :: t, cf
      integer :: k

      if (x < 4.0_rp) then
         r = exp(x*x)*erfc(x)
      else
         ! Lentz-style continued fraction:
         ! erfcx(x) = (1/sqrt(pi)) / (x + 1/2/(x + 1/(x + 3/2/(x + ...))))
         cf = 0.0_rp
         do k = 60, 1, -1
            t = 0.5_rp*real(k, rp)
            cf = t/(x + cf)
         end do
         r = one_over_sqrt_pi/(x + cf)
      end if
   end function erfcx_pos

   !> 1. The guard reproduces the analytic limits; the stale form does not.
   subroutine test_guard_identity()
      real(rp) :: dgi, qppz, alamda, exf_g, exh_g, exf_r, exh_r, exf_s, exh_s
      real(rp) :: dg_list(3), qz_list(3)
      integer :: i, j
      logical :: any_tripped, stale_caught

      alamda = 4.0_rp   ! the value set in charge%wkinit
      dg_list = [0.5_rp, 2.0_rp, 6.0_rp]
      qz_list = [7.0_rp, 9.0_rp, 12.0_rp]   ! large separations: guard territory
      any_tripped = .false.
      stale_caught = .false.

      do i = 1, size(dg_list)
         do j = 1, size(qz_list)
            dgi = dg_list(i)
            qppz = qz_list(j)
            if (erfc(dgi/(2.0_rp*alamda) + alamda*qppz) /= 0.0_rp) cycle
            any_tripped = .true.

            call exf_exh_guarded(dgi, qppz, alamda, exf_g, exh_g)
            call exf_exh_reference(dgi, qppz, alamda, exf_r, exh_r)
            ! Stale value from a "previous iteration": O(1) and of the wrong sign.
            call exf_exh_stale(dgi, qppz, alamda, 0.37_rp, exf_s, exh_s)

            ! The guard must match the analytic limit EXH = -EXF exactly.
            if (exh_g /= -exf_g) then
               write (*, '(a,2f8.3)') 'FAIL guard EXH /= -EXF at dgi,qppz=', dgi, qppz
               failed = .true.
            end if
            if (.not. close_rel(exf_g, exf_r)) then
               write (*, '(a,2f8.3,2es16.7)') 'FAIL guard EXF vs ref at dgi,qppz=', dgi, qppz, exf_g, exf_r
               failed = .true.
            end if
            if (.not. close_rel(exh_g, exh_r)) then
               write (*, '(a,2f8.3,2es16.7)') 'FAIL guard EXH vs ref at dgi,qppz=', dgi, qppz, exh_g, exh_r
               failed = .true.
            end if
            ! The stale form must be demonstrably wrong -- otherwise this test
            ! would silently pass even with the bug reinstated.
            if (.not. close_rel(exh_s, exh_r)) stale_caught = .true.
         end do
      end do

      if (.not. any_tripped) then
         write (*, '(a)') 'FAIL guard branch never tripped -- test is vacuous'
         failed = .true.
      end if
      if (any_tripped .and. .not. stale_caught) then
         write (*, '(a)') 'FAIL stale-EXH form was not distinguishable from the reference'
         failed = .true.
      end if
      if (.not. failed) write (*, '(a)') 'ok  guard identity (EXH = -EXF) in the ERFCP==0 regime'
   end subroutine test_guard_identity

   !> 2. EXH is smooth and sign-definite across the guard boundary.
   subroutine test_guard_continuity()
      real(rp) :: alamda, dgi, qppz, exf, exh, exh_prev, ratio
      integer :: i
      logical :: crossed

      alamda = 4.0_rp
      dgi = 1.0_rp
      exh_prev = 0.0_rp
      crossed = .false.

      ! Sweep QPPZ through the ERFCP underflow threshold in fine steps.
      do i = 1, 400
         qppz = 5.0_rp + 0.02_rp*real(i, rp)
         call exf_exh_guarded(dgi, qppz, alamda, exf, exh)

         ! EXH = -EXPM*ERFCM < 0 throughout this regime.
         if (exh >= 0.0_rp) then
            write (*, '(a,f8.3,es16.7)') 'FAIL EXH not negative at qppz=', qppz, exh
            failed = .true.
            exit
         end if
         ! Monotone decay in magnitude, no jump: successive values differ by a
         ! bounded ratio. A stale-EXH reuse shows up here as an O(1) jump.
         if (i > 1) then
            ratio = abs(exh)/abs(exh_prev)
            if (ratio > 1.0_rp + 1.0e-9_rp) then
               write (*, '(a,f8.3,es16.7)') 'FAIL |EXH| increased at qppz=', qppz, ratio
               failed = .true.
               exit
            end if
         end if
         if (erfc(dgi/(2.0_rp*alamda) + alamda*qppz) == 0.0_rp) crossed = .true.
         exh_prev = exh
      end do

      if (.not. crossed) then
         write (*, '(a)') 'FAIL continuity sweep never crossed the guard threshold'
         failed = .true.
      end if
      if (.not. failed) write (*, '(a)') 'ok  EXH continuous and sign-definite across the guard boundary'
   end subroutine test_guard_continuity

   !> 3. The accumulated SUM1G -> DMDL (the DSZ driver) for a slab-like geometry.
   subroutine test_dsz_accumulation()
      real(rp) :: alamda, qppz, sum1g_g, sum1g_r, sum1g_s
      real(rp) :: dgi, phase, exf, exh, erfcp, erfcm, beta, aqppz
      real(rp) :: exh_prev_g, exh_prev_s, dmdl_g, dmdl_r, dmdl_s
      real(rp) :: facdk, ar2d, sws, pi
      integer :: ig
      logical :: tripped

      pi = acos(-1.0_rp)
      alamda = 4.0_rp
      ar2d = 7.0_rp          ! representative 2D cell area (a.u.^2)
      sws = 2.66_rp          ! representative Wigner-Seitz radius (a.u.)
      facdk = sws*sqrt(3.0_rp)*pi/ar2d

      ! A far layer pair, as found between a deep-bulk and a deep-vacuum row.
      qppz = 8.5_rp
      tripped = .false.

      ! SUM1G is seeded with (ERFCM - ERFCP) at G = 0, then accumulates -DMG.
      beta = 0.0_rp
      aqppz = alamda*qppz
      erfcp = erfc(aqppz)
      erfcm = 2.0_rp - erfcp
      sum1g_g = erfcm - erfcp
      sum1g_r = sum1g_g
      sum1g_s = sum1g_g

      exh_prev_g = 0.0_rp
      exh_prev_s = 0.0_rp

      ! A representative shell of reciprocal vectors.
      do ig = 1, 40
         dgi = 0.35_rp*real(ig, rp)
         phase = cos(0.3_rp*real(ig, rp))   ! stands in for the structure phase

         if (erfc(dgi/(2.0_rp*alamda) + alamda*qppz) == 0.0_rp) tripped = .true.

         call exf_exh_guarded(dgi, qppz, alamda, exf, exh)
         sum1g_g = sum1g_g - phase*exh
         exh_prev_g = exh

         call exf_exh_reference(dgi, qppz, alamda, exf, exh)
         sum1g_r = sum1g_r - phase*exh

         call exf_exh_stale(dgi, qppz, alamda, exh_prev_s, exf, exh)
         sum1g_s = sum1g_s - phase*exh
         exh_prev_s = exh
      end do

      dmdl_g = facdk*sum1g_g
      dmdl_r = facdk*sum1g_r
      dmdl_s = facdk*sum1g_s

      if (.not. tripped) then
         write (*, '(a)') 'FAIL accumulation test never entered the guard branch'
         failed = .true.
      end if
      if (.not. close_rel(dmdl_g, dmdl_r)) then
         write (*, '(a,2es16.7)') 'FAIL DMDL guarded vs reference: ', dmdl_g, dmdl_r
         failed = .true.
      end if
      if (close_rel(dmdl_s, dmdl_r)) then
         write (*, '(a)') 'FAIL stale-EXH DMDL indistinguishable from reference'
         failed = .true.
      end if
      if (.not. failed) then
         write (*, '(a,es12.4)') 'ok  DMDL (DSZ driver) matches reference; stale form deviates by ', &
            abs(dmdl_s - dmdl_r)
      end if
   end subroutine test_dsz_accumulation

   !> 4. Source-level guard. The three tests above pin the MATHEMATICS using local
   !> copies of the branch, but they cannot see charge.f90 itself: madl2d is a
   !> type-bound procedure needing a fully constructed charge object (lattice, 2D
   !> lattice vectors, symbolic atoms), which is out of reach for a
   !> dependency-free unit test. Without this check the suite would still pass if
   !> `exh = -exf` regressed to `exh = -exh` in the real kernel. So inspect the
   !> source directly. MADL2D_SOURCE is supplied by CMake as an absolute path.
   subroutine test_source_has_fix()
      character(len=*), parameter :: src = MADL2D_SOURCE
      integer :: unit, ios
      character(len=512) :: line
      character(len=512) :: squashed
      logical :: found_fix, found_bug, exists

      inquire (file=src, exist=exists)
      if (.not. exists) then
         write (*, '(a)') 'FAIL cannot open '//trim(src)//' for the source-level check'
         failed = .true.
         return
      end if

      found_fix = .false.
      found_bug = .false.
      open (newunit=unit, file=src, status='old', action='read')
      do
         read (unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         squashed = despace(lowercase(line))
         ! Ignore the commented-out prose describing the historical bug.
         if (index(squashed, '!') > 0) then
            if (index(squashed, '!') <= index(squashed, 'exh=-')) cycle
         end if
         if (index(squashed, 'exh=-exf') > 0) found_fix = .true.
         if (index(squashed, 'exh=-exh') > 0) found_bug = .true.
      end do
      close (unit)

      if (found_bug) then
         write (*, '(a)') 'FAIL charge.f90 still contains the stale `exh = -exh` overflow guard'
         failed = .true.
      end if
      if (.not. found_fix) then
         write (*, '(a)') 'FAIL charge.f90 does not contain the corrected `exh = -exf` assignment'
         failed = .true.
      end if
      if (.not. found_bug .and. found_fix) then
         write (*, '(a)') 'ok  charge.f90 carries the corrected `exh = -exf` overflow guard'
      end if
   end subroutine test_source_has_fix

   pure function lowercase(s) result(r)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: r
      integer :: i, c

      do i = 1, len(s)
         c = iachar(s(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) then
            r(i:i) = achar(c + 32)
         else
            r(i:i) = s(i:i)
         end if
      end do
   end function lowercase

   pure function despace(s) result(r)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: r
      integer :: i, n

      r = ''
      n = 0
      do i = 1, len(s)
         if (s(i:i) /= ' ' .and. s(i:i) /= achar(9)) then
            n = n + 1
            r(n:n) = s(i:i)
         end if
      end do
   end function despace

   pure function close_rel(a, b) result(ok)
      real(rp), intent(in) :: a, b
      logical :: ok
      real(rp) :: scale

      scale = max(abs(a), abs(b), 1.0e-300_rp)
      ok = abs(a - b) <= tol*scale
   end function close_rel

end program test_madl2d_exh_guard
