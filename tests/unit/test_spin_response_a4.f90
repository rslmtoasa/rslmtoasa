!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!> @brief spin_response accumulate_chi0 and solve_dyson (matrix kernel) on the developer oracle A4:
!>        self-consistent two-channel Stoner ferromagnet with an exchange matrix.
!> @details Usage: test_spin_response_a4 <oracle directory>. Bands, occupations, amplitudes (random band phases),
!>          expected chi0, chi and moments come from <dir>/a4_two_channel.nml, tolerances from <dir>/tolerances.nml.
!>          Matrix kernel of the oracle = the exchange matrix I entering (1 + chi0 I) chi = chi0, Delta = I m.
!>          Relative error: max|x - ref|/max|ref|, on complex values for arrays.
!------------------------------------------------------------------------------
program test_spin_response_a4
   use precision_mod, only: rp
   use spin_response_mod, only: accumulate_chi0, solve_dyson
   use spin_response_tolerances_mod
   implicit none

   integer, parameter :: nmax = 4096
   character(len=512) :: dir, array_dims, gauge
   integer :: nk, nband, norb, nkappa, nomega, nq, iq, iu, n5, i
   real(rp) :: ts, td, cs, cd, v, nel, kt, eta, ef, kernel(4), moment(2), delta(2), q
   real(rp) :: omega(nmax), k(nmax), w(nmax), e(nmax), f(nmax), a_re(nmax), a_im(nmax)
   namelist /a4_inputs/ nk, nband, norb, nkappa, ts, td, cs, cd, v, nel, kt, eta, kernel, ef, moment, delta, array_dims, &
      gauge, nomega, omega, k, w, e, f, a_re, a_im, nq
   real(rp) :: e_kq(nmax), f_kq(nmax), a_kq_re(nmax), a_kq_im(nmax), chi0_re(nmax), chi0_im(nmax), chi_re(nmax), chi_im(nmax)
   real(rp) :: chi0_static_re(4), chi0_static_im(4), eig_static(2), zeroth_moment(2), goldstone_residual
   real(rp) :: control_lambda(2), control_eigenvalue_along_m(2)
   namelist /a4_q/ iq, q, e_kq, f_kq, a_kq_re, a_kq_im, chi0_re, chi0_im, chi_re, chi_im, chi0_static_re, chi0_static_im, &
      eig_static, zeroth_moment, goldstone_residual, control_lambda, control_eigenvalue_along_m

   real(rp), allocatable :: e_k(:, :), f_k(:, :), e_q(:, :), f_q(:, :), kmat(:, :)
   complex(rp), allocatable :: a_k(:, :, :, :, :), a_q(:, :, :, :, :), chi0(:, :, :), chi(:, :, :), ref(:, :, :), s(:, :, :), hf(:, :, :)
   real(rp) :: ctrl(2), mom(2), chi0s(2, 2)
   logical :: failed

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_spin_response_a4 <oracle directory>'
   call read_tolerance_file(trim(dir))
   call require('a4_rel', a4_rel)
   call require('a4_identity_rel', a4_identity_rel)
   failed = .false.
   write (*, '(a)') 'test | measured | tolerance'

   open (newunit=iu, file=trim(dir)//'/a4_two_channel.nml', status='old', action='read')
   read (iu, nml=a4_inputs)
   n5 = nk*nband*norb*nkappa*2
   kmat = reshape(kernel(1:4), [2, 2])
   e_k = reshape(e(1:nk*nband), [nk, nband])
   f_k = reshape(f(1:nk*nband), [nk, nband])
   a_k = reshape(cmplx(a_re(1:n5), a_im(1:n5), rp), [nk, nband, norb, nkappa, 2])
   allocate (chi0(2, 2, nomega), chi(2, 2, nomega), s(2, 2, 1), hf(2, 2, 2))
   do iq = 1, nq
      read (iu, nml=a4_q)
      e_q = reshape(e_kq(1:nk*nband), [nk, nband])
      f_q = reshape(f_kq(1:nk*nband), [nk, nband])
      a_q = reshape(cmplx(a_kq_re(1:n5), a_kq_im(1:n5), rp), [nk, nband, norb, nkappa, 2])
      chi0 = (0.0_rp, 0.0_rp)
      call accumulate_chi0(kt, w(1:nk), e_k, f_k, a_k, e_q, f_q, a_q, omega(1:nomega), eta, chi0)
      call solve_dyson(chi0, kmat, chi)
      ref = reshape(cmplx(chi0_re(1:4*nomega), chi0_im(1:4*nomega), rp), [2, 2, nomega], order=[3, 1, 2])
      ! catches: wrong vertex (extra 2 or 4 pi, missing 1/2, conjugation), wrong k+q pairing, channel mixing
      call report('A4 chi0, q #'//char(48 + iq), maxval(abs(chi0 - ref))/maxval(abs(ref)), a4_rel)
      ref = reshape(cmplx(chi_re(1:4*nomega), chi_im(1:4*nomega), rp), [2, 2, nomega], order=[3, 1, 2])
      ! catches: Dyson with the kernel transposed or not applied as a matrix over channels
      call report('A4 chi, q #'//char(48 + iq), maxval(abs(chi - ref))/maxval(abs(ref)), a4_rel)
      s = (0.0_rp, 0.0_rp)
      call accumulate_chi0(kt, w(1:nk), e_k, f_k, a_k, e_q, f_q, a_q, [0.0_rp], 0.0_rp, s)
      ! catches: divided-difference rule wrong for the static value
      call report('A4 static chi0 at eta = 0, q #'//char(48 + iq), &
                  maxval(abs(s(:, :, 1) - reshape(cmplx(chi0_static_re, chi0_static_im, rp), [2, 2])))/ &
                  maxval(abs(s(:, :, 1))), a4_rel)
      ! zeroth moment per channel from the high-frequency limit: (w chi(w) - w chi(-w))/2 -> sum of residues, error (d/w)^2
      hf = (0.0_rp, 0.0_rp)
      call accumulate_chi0(kt, w(1:nk), e_k, f_k, a_k, e_q, f_q, a_q, [1.0e7_rp, -1.0e7_rp], eta, hf)
      mom = [(0.5_rp*real(1.0e7_rp*(hf(i, i, 1) - hf(i, i, 2)), rp), i=1, 2)]
      ! catches: normalisation of chi0 per channel (sum rule); expected is the oracle residue sum, which equals m
      call report('A4 zeroth moment per channel, q #'//char(48 + iq), maxval(abs(mom - zeroth_moment))/maxval(abs(moment)), a4_identity_rel)
      if (iq == 1) then
         chi0s = real(s(:, :, 1), rp)
         ! catches: static chi0 not satisfying chi0(0,0) Delta = -m
         call report('A4 Goldstone identity chi0 Delta + m', maxval(abs(matmul(chi0s, delta) + moment))/maxval(abs(moment)), a4_identity_rel)
         ! catches: zero eigenvalue along m lost with the exchange matrix as kernel
         call report('A4 (1 + chi0 I) m', maxval(abs(moment + matmul(chi0s, matmul(kmat, moment))))/maxval(abs(moment)), a4_identity_rel)
         do i = 1, 2
            ctrl = moment + control_lambda(i)*matmul(chi0s, matmul(kmat, moment)) - control_eigenvalue_along_m(i)*moment
            ! catches: eigenvalue along m not 1 - lambda for a scaled kernel
            call report('A4 control lambda ='//merge(' 2   ', ' 1/2 ', control_lambda(i) > 1.0_rp), maxval(abs(ctrl))/maxval(abs(moment)), a4_identity_rel)
         end do
      end if
   end do
   close (iu)
   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine report(label, measured, tol)
      character(*), intent(in) :: label
      real(rp), intent(in) :: measured, tol

      write (*, '(a,t58,es10.3,2x,es10.3,2x,a)') label, measured, tol, merge('ok  ', 'FAIL', measured <= tol)
      if (.not. measured <= tol) failed = .true.
   end subroutine report

end program test_spin_response_a4
