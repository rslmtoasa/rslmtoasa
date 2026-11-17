! UnitGbtKspaceHorth
! Relations C1 (eigenvalues of H^gamma are the roots of det(P^alpha(e) - S^alpha) = 0; Andersen & Jepsen, PRL 53, 2571 (1984)) and R2
! (theta = 0: spin up carries k - q/2, spin down k + q/2; Sandratskii 1998); oracle O3 (tests/gbt/oracles/o3_hgamma, amendment 07972fa).
! Chain under test: gbt_kspace_h1 with w = 1, c = 0 gives S^alpha_rot, gbt_kspace_horth gives H_rot, zheev its eigenvalues. Expected values come
! from the oracle file only. Bugs caught: wrong sign or side of Gamma, dele used as Delta, a spin/k swap in the rotated link, an anti-Hermitian H.
! usage: test_gbt_kspace_horth <oracle directory>
program test_gbt_kspace_horth
   use precision_mod, only: rp
   use gbt_kspace_mod, only: gbt_kspace_h1, gbt_kspace_horth
   use gbt_structure_mod, only: gbt_reference_frame
   implicit none

   integer, parameter :: nb = 15, norb = 9, kidx(5) = [0, 7, 19, 36, 53]
   real(rp), parameter :: pi = acos(-1.0_rp)
   character(len=:), allocatable :: txt
   character(len=512) :: dir
   real(rp) :: eig_identity_ry, omega_identity_ry, n_identity, model_eig_ry, model_ft_rel, scf_identity_ry, scf_moment_mub
   namelist /gbt_tolerances/ eig_identity_ry, omega_identity_ry, n_identity, model_eig_ry, model_ft_rel, scf_identity_ry, scf_moment_mub
   integer, parameter :: ones(nb) = 1
   character(len=17), parameter :: names(3) = [character(len=17) :: 'cc', 'dele', 'gamma_minus_alpha']
   real(rp) :: nd(3, nb), par(norb, 2, 3), kk(3, 5), antiherm, dev(3), ev(2*norb), flat(norb*norb*nb), flat2(norb*2)
   real(rp) :: want(2*norb), q(3), theta_deg(2), qq(3, 2), eg(norb, 2, 5)
   complex(rp) :: sblock(norb, norb, nb)
   integer :: u, ib, ik, iq, it, ia, ib2, ic, pairs(2, 3)
   character(len=64) :: key

   call get_command_argument(1, dir)
   if (len_trim(dir) == 0) error stop 'usage: test_gbt_kspace_horth <oracle directory>'
   open (newunit=u, file=trim(dir)//'/tolerances.nml', status='old', action='read')
   read (u, nml=gbt_tolerances)
   close (u)
   call load(trim(dir)//'/o3_hgamma/o3_fe_bcc.nml', txt)
   call get('bond_n_direct', flat(1:3*nb))
   nd = reshape(flat(1:3*nb), [3, nb])
   call get('bond_s_re', flat)
   sblock = reshape(cmplx(flat, 0.0_rp, rp), [norb, norb, nb])
   call get('bond_s_im', flat)
   sblock = sblock + reshape(cmplx(0.0_rp, flat, rp), [norb, norb, nb])
   do ic = 1, 3
      call get(trim(names(ic)), flat2)
      par(:, :, ic) = reshape(flat2, [norb, 2])
   end do
   do ik = 1, 5
      write (key, '(a,i0,a)') 'k', kidx(ik), '_direct'
      call get(trim(key), kk(:, ik))
      do ib = 1, 2
         write (key, '(a,i0,a,i0,a)') 'k', kidx(ik), '_spin', ib, '_h_gamma'
         call get(trim(key), eg(:, ib, ik))
      end do
   end do
   antiherm = 0.0_rp
   theta_deg = [5.0_rp, 30.0_rp]
   qq(:, 1) = [0.125_rp, -0.125_rp, 0.125_rp]
   qq(:, 2) = [0.25_rp, -0.25_rp, 0.25_rp]

   ! C1, collinear: theta = 0, q = 0 -> the eigenvalues of H^gamma_up(k) and H^gamma_dn(k)
   dev(1) = 0.0_rp
   do ik = 1, 5
      call rot_eigs(kk(:, ik), [0.0_rp, 0.0_rp, 0.0_rp], 0.0_rp, ev)
      dev(1) = max(dev(1), maxval(abs(ev - sorted([eg(:, 1, ik), eg(:, 2, ik)]))))
   end do
   ! C1, rotated: H_rot eigenvalues at theta = 5 and 30 deg, both q
   dev(2) = 0.0_rp
   do ik = 1, 5
      do iq = 1, 2
         do it = 1, 2
            write (key, '(a,i0,a,i0,a,i0,a)') 'k', kidx(ik), '_q', iq, '_theta', nint(theta_deg(it)), '_h_rot'
            call get(trim(key), want)
            call rot_eigs(kk(:, ik), qq(:, iq), theta_deg(it)*pi/180.0_rp, ev)
            dev(2) = max(dev(2), maxval(abs(ev - want)))
         end do
      end do
   end do
   ! R2: theta = 0 with k - q/2 = k_a and k + q/2 = k_b taken from the oracle k-points: up sees k_a, down sees k_b
   pairs = reshape([1, 2, 3, 4, 2, 5], [2, 3])
   dev(3) = 0.0_rp
   do ik = 1, 3
      ia = pairs(1, ik)
      ib2 = pairs(2, ik)
      q = kk(:, ib2) - kk(:, ia)
      call rot_eigs(0.5_rp*(kk(:, ia) + kk(:, ib2)), q, 0.0_rp, ev)
      dev(3) = max(dev(3), maxval(abs(ev - sorted([eg(:, 1, ia), eg(:, 2, ib2)]))))
   end do

   write (*, '(a,3(es10.2,1x),a,es9.1)') 'O3 max |eigenvalue - oracle| [Ry] (C1 collinear, C1 rotated, R2 theta=0): ', dev, ' tolerance ', eig_identity_ry
   write (*, '(a,es10.2)') 'max |H - H^dag| / max |H| before any symmetrisation (limit 1e-12): ', antiherm
   if (maxval(dev) >= eig_identity_ry .or. antiherm > 1.0e-12_rp) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine rot_eigs(k, qdir, theta, evals)
      real(rp), intent(in) :: k(3), qdir(3), theta
      real(rp), intent(out) :: evals(2*norb)
      complex(rp) :: srot(2*norb, 2*norb), h(2*norb, 2*norb), cmat(2*norb, 2*norb), frame(2, 2, 1), w(norb, 2, 1), c(norb, 2, 1), work(64)
      real(rp) :: rwork(6*norb), cc18(2*norb), dele18(2*norb), gam18(2*norb)
      integer :: info, i
      call gbt_reference_frame(theta, 0.0_rp, frame(:, :, 1))
      w = cmplx(1.0_rp, 0.0_rp, rp)
      c = cmplx(0.0_rp, 0.0_rp, rp)
      call gbt_kspace_h1(norb, 1, nb, ones, ones, nd, sblock, frame, w, c, k, qdir, srot)
      cc18 = [par(:, 1, 1), par(:, 2, 1)]
      dele18 = [par(:, 1, 2), par(:, 2, 2)]
      gam18 = [par(:, 1, 3), par(:, 2, 3)]
      cmat = (0.0_rp, 0.0_rp)
      do i = 1, 2*norb
         cmat(i, i) = cc18(i)
      end do
      call gbt_kspace_horth(2*norb, srot, cmat, dele18, gam18, h)
      antiherm = max(antiherm, maxval(abs(h - conjg(transpose(h))))/maxval(abs(h)))
      call zheev('N', 'U', 2*norb, h, 2*norb, evals, work, 64, rwork, info)
      if (info /= 0) error stop 'zheev failed'
   end subroutine rot_eigs

   function sorted(x) result(y)
      real(rp), intent(in) :: x(:)
      real(rp) :: y(size(x))
      integer :: i, j

      y = x
      do i = 2, size(y)
         j = i - 1
         do while (j >= 1)
            if (y(j) <= y(j + 1)) exit
            y(j:j + 1) = y(j + 1:j:-1)
            j = j - 1
         end do
      end do
   end function sorted

   subroutine load(path, text)
      character(len=*), intent(in) :: path
      character(len=:), allocatable, intent(out) :: text
      integer :: unit, sz

      inquire (file=path, size=sz)
      allocate (character(len=sz) :: text)
      open (newunit=unit, file=path, access='stream', form='unformatted', status='old', action='read')
      read (unit) text
      close (unit)
   end subroutine load

   !> Values of "  key = v1, v2, ..." (one line of the oracle file) in file order.
   subroutine get(key, x)
      character(len=*), intent(in) :: key
      real(rp), intent(out) :: x(:)
      integer :: i, j

      i = index(txt, new_line('a')//'  '//key//' =')
      if (i == 0) error stop 'oracle key missing'
      i = i + len(key) + 5
      j = i + index(txt(i:), new_line('a')) - 2
      read (txt(i:j), *) x
   end subroutine get

end program test_gbt_kspace_horth
