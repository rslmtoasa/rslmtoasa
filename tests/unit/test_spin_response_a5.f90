!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!> @brief spin_response kernels on the developer oracle A5: self-consistent two-sublattice antiferromagnet (af) and
!>        metallic ferrimagnet (ferri), global spin frame, site B moment negative.
!> @details Usage: test_spin_response_a5 <oracle directory>. Bands, occupations, amplitudes (the af amplitudes mix each
!>          degenerate up/down pair by a random unitary), expected chi0, chi, static chi0, Omega and exponents come from
!>          <dir>/a5_af.nml and <dir>/a5_ferri.nml, tolerances from <dir>/tolerances.nml. Kernel of the oracle = diag(u_model).
!>          Relative error: max|x - ref|/max|ref|, on complex values for arrays. The af q = 0 group has no omega_adiabatic
!>          (Omega is a Jordan block there), so it is not compared.
!------------------------------------------------------------------------------
program test_spin_response_a5
   use precision_mod, only: rp
   use spin_response_mod, only: accumulate_chi0, solve_dyson, u_juelich, adiabatic_frequencies
   use spin_response_tolerances_mod
   implicit none

   integer, parameter :: nmax = 4096
   character(len=512) :: dir, array_dims, gauge, frame
   character(len=16) :: set
   integer :: nk, nband, norb, nsite, nomega, nq, iq, iu, n5, iset
   real(rp) :: t, eps(2), u_model(2), nel, kt, eta, ef, moment(2), delta(2), q
   real(rp) :: omega(nmax), k(nmax), w(nmax), e(nmax), f(nmax), a_re(nmax), a_im(nmax)
   namelist /a5_inputs/ set, nk, nband, norb, nsite, t, eps, u_model, nel, kt, eta, ef, moment, delta, frame, array_dims, &
      gauge, nomega, omega, k, w, e, f, a_re, a_im, nq
   real(rp) :: e_kq(nmax), f_kq(nmax), a_kq_re(nmax), a_kq_im(nmax), chi0_re(nmax), chi0_im(nmax), chi_re(nmax), chi_im(nmax)
   real(rp) :: chi0_static_re(4), chi0_static_im(4), eig_static(2), zeroth_moment(2), goldstone_residual, u_goldstone(2)
   real(rp) :: omega_adiabatic(2)
   namelist /a5_q/ iq, q, e_kq, f_kq, a_kq_re, a_kq_im, chi0_re, chi0_im, chi_re, chi_im, chi0_static_re, chi0_static_im, &
      eig_static, zeroth_moment, goldstone_residual, u_goldstone, omega_adiabatic
   character(len=32) :: branch
   integer :: ladder_j(3)
   real(rp) :: omega_branch(3), exponent(2)
   namelist /a5_exponent/ branch, ladder_j, omega_branch, exponent

   real(rp), allocatable :: e_k(:, :), f_k(:, :), e_q(:, :), f_q(:, :), kmat(:, :)
   complex(rp), allocatable :: a_k(:, :, :, :, :), a_q(:, :, :, :, :), chi0(:, :, :), chi(:, :, :), ref(:, :, :), s(:, :, :), hf(:, :, :)
   real(rp) :: chi0s(2, 2), mom(2), uj(2), wad(2), imag_max, p(2), p_expected, br(3)
   character(len=3) :: tag
   logical :: failed

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_spin_response_a5 <oracle directory>'
   call read_tolerance_file(trim(dir))
   call require('a5_rel', a5_rel)
   call require('a5_identity_rel', a5_identity_rel)
   call require('a5_adiabatic_rel', a5_adiabatic_rel)
   call require('a5_exponent_abs', a5_exponent_abs)
   failed = .false.
   write (*, '(a)') 'test | measured | tolerance'

   do iset = 1, 2
      open (newunit=iu, file=trim(dir)//merge('/a5_af.nml   ', '/a5_ferri.nml', iset == 1), status='old', action='read')
      read (iu, nml=a5_inputs)
      tag = merge('af ', 'fe ', iset == 1)
      n5 = nk*nband*norb*nsite*2
      kmat = reshape([u_model(1), 0.0_rp, 0.0_rp, u_model(2)], [2, 2])
      e_k = reshape(e(1:nk*nband), [nk, nband])
      f_k = reshape(f(1:nk*nband), [nk, nband])
      a_k = reshape(cmplx(a_re(1:n5), a_im(1:n5), rp), [nk, nband, norb, nsite, 2])
      allocate (chi0(2, 2, nomega), chi(2, 2, nomega), s(2, 2, 1), hf(2, 2, 2))
      do iq = 1, nq
         read (iu, nml=a5_q)
         e_q = reshape(e_kq(1:nk*nband), [nk, nband])
         f_q = reshape(f_kq(1:nk*nband), [nk, nband])
         a_q = reshape(cmplx(a_kq_re(1:n5), a_kq_im(1:n5), rp), [nk, nband, norb, nsite, 2])
         chi0 = (0.0_rp, 0.0_rp)
         call accumulate_chi0(kt, w(1:nk), e_k, f_k, a_k, e_q, f_q, a_q, omega(1:nomega), eta, chi0)
         call solve_dyson(chi0, kmat, chi)
         ref = reshape(cmplx(chi0_re(1:4*nomega), chi0_im(1:4*nomega), rp), [2, 2, nomega], order=[3, 1, 2])
         ! catches: wrong vertex or k+q pairing, a dependence on the random up/down mixing of degenerate pairs (af),
         ! the vertex of the down sublattice taken in a local frame (cross-sublattice chi0 lost)
         call report('A5 '//tag//' chi0, q #'//char(48 + iq), maxval(abs(chi0 - ref))/maxval(abs(ref)), a5_rel)
         ref = reshape(cmplx(chi_re(1:4*nomega), chi_im(1:4*nomega), rp), [2, 2, nomega], order=[3, 1, 2])
         ! catches: Dyson with the signed-moment sites handled as a matrix wrongly
         call report('A5 '//tag//' chi, q #'//char(48 + iq), maxval(abs(chi - ref))/maxval(abs(ref)), a5_rel)
         s = (0.0_rp, 0.0_rp)
         call accumulate_chi0(kt, w(1:nk), e_k, f_k, a_k, e_q, f_q, a_q, [0.0_rp], 0.0_rp, s)
         ! catches: static value taken at a finite eta, or without the divided difference for degenerate pairs
         call report('A5 '//tag//' static chi0 at eta = 0, q #'//char(48 + iq), &
                     maxval(abs(s(:, :, 1) - reshape(cmplx(chi0_static_re, chi0_static_im, rp), [2, 2])))/maxval(abs(s(:, :, 1))), a5_rel)
         ! catches: the smallest eigenvalue of 1 + chi0(q, 0) U (stability) from a wrong static chi0 or U
         call report('A5 '//tag//' eig(1 + chi0 U), q #'//char(48 + iq), &
                     maxval(abs(sorted_eigenvalues(s(:, :, 1), u_model) - eig_static))/maxval(abs(eig_static)), a5_rel)
         hf = (0.0_rp, 0.0_rp)
         call accumulate_chi0(kt, w(1:nk), e_k, f_k, a_k, e_q, f_q, a_q, [1.0e7_rp, -1.0e7_rp], eta, hf)
         mom = [0.5_rp*real(1.0e7_rp*(hf(1, 1, 1) - hf(1, 1, 2)), rp), 0.5_rp*real(1.0e7_rp*(hf(2, 2, 1) - hf(2, 2, 2)), rp)]
         ! catches: normalisation of chi0 per site (zeroth moment = the signed moment)
         call report('A5 '//tag//' zeroth moment per site, q #'//char(48 + iq), &
                     maxval(abs(mom - zeroth_moment))/maxval(abs(moment)), a5_identity_rel)
         chi0s = real(s(:, :, 1), rp)
         if (iq == 1) then
            ! catches: static chi0 not satisfying chi0(0,0) Delta = -m for the signed m
            call report('A5 '//tag//' Goldstone identity chi0 Delta + m', maxval(abs(matmul(chi0s, delta) + moment))/maxval(abs(moment)), a5_identity_rel)
            ! catches: the zero eigenvalue lost along the signed moments
            call report('A5 '//tag//' (1 + chi0 U) m', maxval(abs(moment + matmul(chi0s, matmul(kmat, moment))))/maxval(abs(moment)), a5_identity_rel)
            call u_juelich(chi0s, moment, uj)
            ! catches: |M| instead of the signed M in the Goldstone solve, cross-site terms dropped
            call report('A5 '//tag//' U from u_juelich vs u_model', maxval(abs(uj - u_model))/maxval(u_model), a5_identity_rel)
         end if
         if (iset == 2 .or. iq > 1) then
            call adiabatic_frequencies(s(:, :, 1), kmat, moment, wad, imag_max)
            ! catches: unsigned N, chi0^-1 term missing or K sign wrong in Omega, cross-sublattice terms dropped
            call report('A5 '//tag//' omega_adiabatic, q #'//char(48 + iq), maxval(abs(wad - omega_adiabatic))/maxval(abs(omega_adiabatic)), a5_adiabatic_rel)
            call report('A5 '//tag//' max|Im Omega|/max|Re Omega|, q #'//char(48 + iq), imag_max/maxval(abs(wad)), a5_adiabatic_rel)
            if (iq > 1) br(iq - 1) = merge(wad(2), wad(minloc(abs(wad), 1)), iset == 1)
         end if
      end do
      read (iu, nml=a5_exponent)
      close (iu)
      p = log(abs(br(2:3)/br(1:2)))/log(2.0_rp)
      p_expected = merge(1.0_rp, 2.0_rp, iset == 1)
      write (*, '(a,3es13.5,a,2f9.5)') 'A5 '//tag//' branch at q, 2q, 4q (Ry):', br, ' ; exponents', p
      ! catches: a decoupled sublattice model (exponent 2 for the antiferromagnet) or a linear branch where the ferrimagnet
      ! has a quadratic one; the bound separates 1 from 2 by a factor of 20
      call report('A5 '//tag//' exponent minus expected, worst of two', maxval(abs(p - p_expected)), a5_exponent_abs)
      call report('A5 '//tag//' exponent vs the oracle exponent, worst of two', maxval(abs(p - exponent)), a5_adiabatic_rel)
      deallocate (chi0, chi, s, hf)
   end do
   if (failed) error stop 1
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine report(label, measured, tol)
      character(*), intent(in) :: label
      real(rp), intent(in) :: measured, tol

      write (*, '(a,t66,es10.3,2x,es10.3,2x,a)') label, measured, tol, merge('ok  ', 'FAIL', measured <= tol)
      if (.not. measured <= tol) failed = .true.
   end subroutine report

   !> Eigenvalues of I + chi0 diag(u), ascending: those of the Hermitian matrix I + u^(1/2) chi0 u^(1/2), 2 x 2 closed form.
   function sorted_eigenvalues(c, u) result(ev)
      complex(rp), intent(in) :: c(2, 2)
      real(rp), intent(in) :: u(2)
      real(rp) :: ev(2)

      real(rp) :: a, d, r

      a = 1.0_rp + u(1)*real(c(1, 1), rp)
      d = 1.0_rp + u(2)*real(c(2, 2), rp)
      r = sqrt(((a - d)/2.0_rp)**2 + u(1)*u(2)*abs(c(1, 2))**2)
      ev = [(a + d)/2.0_rp - r, (a + d)/2.0_rp + r]
   end function sorted_eigenvalues

end program test_spin_response_a5
