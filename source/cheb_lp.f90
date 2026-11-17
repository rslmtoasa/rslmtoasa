!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Cheb_lp
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Linear prediction (LP) of Chebyshev moments. Extends the known moments
!> mu(:,:,1:n_known) of one block to the full length of mu by a linear
!> recurrence x(n) = -sum_i a_i x(n-i) with one real coefficient vector shared
!> by the real and imaginary parts of all nb*nb elements. Hermiticity of the
!> block is therefore preserved by construction.
!> Refs: Ganahl et al., PRB 90, 045144; Wolf et al., PRB 90, 115124.
!------------------------------------------------------------------------------

module cheb_lp_mod

   use control_mod, only: control
   use logger_mod, only: g_logger
   use precision_mod, only: rp
   use string_mod, only: int2str, real2str
   implicit none

   private

   public :: cheb_lp, cheb_lp_extend

   !> Linear-prediction settings taken from the control namelist.
   type :: cheb_lp
      !> Extend the Chebyshev moments before reconstruction.
      logical :: enabled = .false.
      !> Extension factor: n_total = factor*(2*lld+2).
      integer :: factor = 2
   contains
      procedure :: n_total
      procedure :: apply
   end type cheb_lp

   interface cheb_lp
      procedure :: constructor
   end interface cheb_lp

   !> Order p = n_known/p_div and window start n_lo = n_known/nlo_div. The training data are all
   !> n_known moments (T = n_known) and p = T/2, the choice of Ganahl et al., PRB 90, 045144.
   integer, parameter :: p_div = 2
   integer, parameter :: nlo_div = 2
   !> Pseudo-inverse cutoff, relative to the largest eigenvalue of the normal matrix (PRB 90, 045144 uses 1e-5 to 1e-6).
   real(rp), parameter :: delta = 1.0e-6_rp
   !> Largest |mu(l,m,n) - conjg(mu(m,l,n))| over the known moments, relative to max|mu|, for which the block is treated as Hermitian.
   real(rp), parameter :: herm_tol = 1.0e-8_rp

contains

   function constructor(ctrl) result(obj)
      type(control), intent(in) :: ctrl
      type(cheb_lp) :: obj

      obj%enabled = ctrl%cheb_lp
      obj%factor = ctrl%cheb_lp_factor
   end function constructor

   !> @brief Length of the stored moment array: 2*lld+2, times the factor when enabled.
   pure function n_total(this, lld) result(n)
      class(cheb_lp), intent(in) :: this
      integer, intent(in) :: lld
      integer :: n

      n = 2*lld + 2
      if (this%enabled) n = n*this%factor
   end function n_total

   !> @brief Extend one moment block when enabled and log the fit diagnostics.
   !> @param[inout] mu Moments (nb,nb,n_total) of one atom or pair slot; 1:n_known in.
   !> @param[in] n_known Number of known moments (2*lld+2).
   !> @param[in] label Atom or pair identifier for the log line.
   subroutine apply(this, mu, n_known, label)
      class(cheb_lp), intent(in) :: this
      complex(rp), intent(inout) :: mu(:, :, :)
      integer, intent(in) :: n_known
      character(len=*), intent(in) :: label
      integer :: n_rank, n_unstable
      real(rp) :: lam_max, fit_res, imag_max, unstable_weight

      if (.not. this%enabled) return
      call cheb_lp_extend(mu, n_known, n_rank, n_unstable, lam_max, fit_res, imag_max, unstable_weight)
      !$omp critical (cheb_lp_log)
      call g_logger%info('Chebyshev LP '//label//': p='//int2str(n_known/p_div)//' rank='//int2str(n_rank)// &
                         ' unstable='//int2str(n_unstable)//' (weight '//real2str(unstable_weight, '(f6.3)')//') lam_max='//real2str(lam_max, '(f9.6)')// &
                         ' fit_res='//real2str(fit_res, '(es9.2)')//' imag='//real2str(imag_max, '(es9.2)'), &
                         __FILE__, __LINE__)
      !$omp end critical (cheb_lp_log)
   end subroutine apply

   !> @brief Extend the moments mu(:,:,n_known+1:) by linear prediction.
   !> @details Index i of mu is order i-1. Fit window: orders n_lo..n_known-1.
   !>          Roots with |lambda| >= 1 are dropped before extrapolating.
   !>          A Hermitian block (to herm_tol) is fitted on its nb*nb unique real channels, Re mu(l,m) for l <= m and
   !>          Im mu(l,m) for l < m, with weights 1 on the diagonal and 2 off it, which reproduces the sums over all 2*nb*nb
   !>          channels; the extension is stored in both triangles. Any other block is fitted on all 2*nb*nb channels.
   !> @param[inout] mu Moments (nb,nb,n_total); 1:n_known in, the rest out.
   !> @param[in] n_known Number of known moments.
   !> @param[out] n_rank Eigenvalues of the normal matrix kept by the cutoff.
   !> @param[out] n_unstable Roots with |lambda| >= 1 before stabilization.
   !> @param[out] lam_max Largest |lambda| before stabilization.
   !> @param[out] fit_res Relative residual of the fit on the window.
   !> @param[out] imag_max Largest imaginary part discarded in the extrapolation.
   !> @param[out] unstable_weight Share of the mode amplitude carried by the dropped roots (small: the drop is harmless).
   subroutine cheb_lp_extend(mu, n_known, n_rank, n_unstable, lam_max, fit_res, imag_max, unstable_weight)
      complex(rp), intent(inout) :: mu(:, :, :)
      integer, intent(in) :: n_known
      integer, intent(out) :: n_rank, n_unstable
      real(rp), intent(out) :: lam_max, fit_res, imag_max, unstable_weight
      integer :: nb, nb2, nch, p, nlo, n_ext, nsteps, nkeep, q, i, j, k, info, lwork
      integer, allocatable :: ipiv(:), cl(:), cm(:)
      logical, allocatable :: cim(:)
      logical :: herm
      real(rp) :: smax, ynorm, rnorm, mult
      real(rp), allocatable :: wch(:), xr(:, :), xw(:, :), tg(:, :), dd(:, :), gg(:, :), rr(:, :), rcopy(:, :), rv(:), ev(:), t(:), a(:)
      real(rp), allocatable :: comp(:, :), wr(:), wi(:), vr(:, :), work(:), amp(:), dre(:, :), bst(:, :), ext(:, :)
      logical, allocatable :: unst(:)
      complex(rp), allocatable :: vc(:, :), bm(:, :), lam(:), v1(:), cv(:), oc(:)
      complex(rp) :: lamk
      real(rp) :: dummy(1, 1)

      nb = size(mu, 1)
      nb2 = nb*nb
      n_ext = size(mu, 3)
      p = n_known/p_div
      nlo = n_known/nlo_div
      if (p < 2 .or. n_ext <= n_known .or. nlo < p) &
         call g_logger%fatal('cheb_lp_extend: need n_known >= 8 and size(mu,3) > n_known', __FILE__, __LINE__)

      herm = is_hermitian()
      nch = merge(nb2, 2*nb2, herm)
      allocate (cl(nch), cm(nch), cim(nch), wch(nch), xr(n_known, nch))
      call load_channels()

      ! Normal equations R a = -r over all (weighted) channels and the window, from the lag Gram matrix T(a,b) = sum_c w_c x_c(a) x_c(b)
      ! (positions a, b = 1..n_known). The window sums are differences along the diagonals of T:
      !   G(i,j) = sum_c w_c sum_{n=nlo+1}^{n_known} x_c(n-i) x_c(n-j) = D(n_known-i, n_known-j) - D(nlo-i, nlo-j),
      !   D(a,b) = T(a,b) + D(a-1,b-1), D = 0 for a = 0 or b = 0.
      ! Then R(i,j) = G(i,j), r(i) = G(i,0) (i,j = 1..p), |y|^2 = G(0,0), and the fit residual is |y|^2 + 2 a.r + a.R.a.
      allocate (tg(n_known, n_known), dd(0:n_known, 0:n_known), gg(0:p, 0:p))
      allocate (xw(n_known, nch))
      do i = 1, nch
         xw(:, i) = sqrt(wch(i))*xr(:, i)
      end do
      call dsyrk('U', 'N', n_known, nch, 1.0_rp, xw, n_known, 0.0_rp, tg, n_known)
      deallocate (xw)
      dd = 0.0_rp
      do j = 1, n_known
         do i = 1, j
            dd(i, j) = tg(i, j) + dd(i - 1, j - 1)
         end do
      end do
      do j = 0, p
         do i = j, p
            gg(i, j) = dd(n_known - i, n_known - j) - dd(nlo - i, nlo - j)
            gg(j, i) = gg(i, j)
         end do
      end do
      deallocate (tg, dd)
      allocate (rr(p, p), rcopy(p, p), rv(p), ev(p), t(p), a(p))
      rr = gg(1:p, 1:p)
      rcopy = rr
      rv = gg(1:p, 0)
      ynorm = gg(0, 0)
      deallocate (gg)
      if (ynorm == 0.0_rp) then
         mu(:, :, n_known + 1:) = (0.0_rp, 0.0_rp)
         n_rank = 0; n_unstable = 0; lam_max = 0.0_rp; fit_res = 0.0_rp; imag_max = 0.0_rp; unstable_weight = 0.0_rp
         return
      end if

      ! Pseudo-inverse solution via eigendecomposition of R.
      allocate (work(3*p))
      call dsyev('V', 'U', p, rr, p, ev, work, 3*p, info)
      if (info /= 0) call g_logger%fatal('cheb_lp_extend: dsyev failed', __FILE__, __LINE__)
      deallocate (work)
      smax = ev(p)
      n_rank = count(ev > delta*smax)
      call dgemv('T', p, p, 1.0_rp, rr, p, rv, 1, 0.0_rp, t, 1)
      where (ev > delta*smax)
         t = t/ev
      elsewhere
         t = 0.0_rp
      end where
      call dgemv('N', p, p, -1.0_rp, rr, p, t, 1, 0.0_rp, a, 1)

      call dgemv('N', p, p, 1.0_rp, rcopy, p, a, 1, 0.0_rp, t, 1)
      rnorm = ynorm + 2.0_rp*dot_product(a, rv) + dot_product(a, t)
      fit_res = sqrt(max(rnorm, 0.0_rp)/ynorm)

      ! Companion matrix, eigenpairs, stabilization (drop |lambda| >= 1).
      allocate (comp(p, p), wr(p), wi(p), vr(p, p), vc(p, p), lam(p))
      comp = 0.0_rp
      comp(1, :) = -a
      do i = 1, p - 1
         comp(i + 1, i) = 1.0_rp
      end do
      lwork = -1
      allocate (work(1))
      call dgeev('N', 'V', p, comp, p, wr, wi, dummy, 1, vr, p, work, lwork, info)
      lwork = int(work(1))
      deallocate (work)
      allocate (work(lwork))
      call dgeev('N', 'V', p, comp, p, wr, wi, dummy, 1, vr, p, work, lwork, info)
      if (info /= 0) call g_logger%fatal('cheb_lp_extend: dgeev failed', __FILE__, __LINE__)
      j = 1
      do while (j <= p)
         if (wi(j) == 0.0_rp) then
            vc(:, j) = cmplx(vr(:, j), 0.0_rp, rp)
            j = j + 1
         else
            vc(:, j) = cmplx(vr(:, j), vr(:, j + 1), rp)
            vc(:, j + 1) = conjg(vc(:, j))
            j = j + 2
         end if
      end do
      lam = cmplx(wr, wi, rp)
      lam_max = maxval(abs(lam))
      unst = abs(lam) >= 1.0_rp
      n_unstable = count(unst)
      where (unst) lam = (0.0_rp, 0.0_rp)

      ! x(n_known-1+k) = Re [V lam^k V^-1 s]_1, s = [x(n_known-1), ..., x(n_known-p)], i.e. Re sum_j v1_j lam_j^k bm_jc.
      ! Dropped roots have lam = 0. A conjugate pair is folded into its member with Im lam > 0 (weight 2), and
      ! Re(z b) = Re z Re b - Im z Im b makes the whole extension one real product over (modes) of steps x channels.
      allocate (v1(p), bm(p, nch), ipiv(p), cv(p), oc(nch))
      v1 = vc(1, :)
      do i = 1, p
         bm(i, :) = cmplx(xr(n_known - i + 1, :), 0.0_rp, rp)
      end do
      call zgesv(p, nch, vc, p, ipiv, bm, p, info)
      if (info /= 0) call g_logger%fatal('cheb_lp_extend: zgesv failed (singular eigenvector matrix)', __FILE__, __LINE__)
      allocate (amp(p))
      do i = 1, p
         amp(i) = abs(v1(i))*sum(wch*abs(bm(i, :)))
      end do
      unstable_weight = sum(amp, mask=unst)/sum(amp)
      nsteps = n_ext - n_known
      nkeep = count(lam /= (0.0_rp, 0.0_rp) .and. aimag(lam) >= 0.0_rp)
      allocate (dre(2*nkeep, nsteps), bst(2*nkeep, nch), ext(nch, nsteps))
      q = 0
      do j = 1, p
         if (lam(j) == (0.0_rp, 0.0_rp) .or. aimag(lam(j)) < 0.0_rp) cycle
         q = q + 1
         mult = merge(1.0_rp, 2.0_rp, aimag(lam(j)) == 0.0_rp)
         lamk = lam(j)
         do k = 1, nsteps
            dre(q, k) = mult*real(v1(j)*lamk, rp)
            dre(nkeep + q, k) = -mult*aimag(v1(j)*lamk)
            lamk = lamk*lam(j)
         end do
         bst(q, :) = real(bm(j, :), rp)
         bst(nkeep + q, :) = aimag(bm(j, :))
      end do
      ext = 0.0_rp
      if (nkeep > 0) call dgemm('T', 'N', nch, nsteps, 2*nkeep, 1.0_rp, bst, 2*nkeep, dre, 2*nkeep, 0.0_rp, ext, nch)
      call store_extension()
      ! imag_max: the imaginary part dropped by Re, evaluated for the last extrapolation step only (not the maximum over steps).
      cv = v1*lam**nsteps
      call zgemv('T', p, nch, (1.0_rp, 0.0_rp), bm, p, cv, 1, (0.0_rp, 0.0_rp), oc, 1)
      imag_max = maxval(abs(aimag(oc)))

   contains

      !> Hermitian to herm_tol: max |mu(l,m,n) - conjg(mu(m,l,n))| <= herm_tol * max|mu| over the known moments.
      logical function is_hermitian()
         real(rp) :: dev
         integer :: n
         dev = 0.0_rp
         do n = 1, n_known
            dev = max(dev, maxval(abs(mu(:, :, n) - conjg(transpose(mu(:, :, n))))))
         end do
         is_hermitian = dev <= herm_tol*maxval(abs(mu(:, :, 1:n_known)))
      end function is_hermitian

      !> Real channels c: Re or Im (cim) of mu(cl,cm,:) with weight wch. Hermitian: l <= m real, l < m imaginary, taken
      !> from the Hermitian part of mu. Otherwise every element, real part c = (m-1)*nb + l, imaginary part nb2 + that.
      subroutine load_channels()
         integer :: l, m, c, n
         complex(rp) :: z

         c = 0
         do m = 1, nb
            do l = 1, nb
               if (herm .and. l > m) cycle
               c = c + 1
               cl(c) = l; cm(c) = m; cim(c) = .false.
               wch(c) = merge(2.0_rp, 1.0_rp, herm .and. l < m)
               if (.not. herm) then
                  cl(nb2 + c) = l; cm(nb2 + c) = m; cim(nb2 + c) = .true.; wch(nb2 + c) = 1.0_rp
               end if
            end do
         end do
         if (herm) then
            do m = 2, nb
               do l = 1, m - 1
                  c = c + 1
                  cl(c) = l; cm(c) = m; cim(c) = .true.; wch(c) = 2.0_rp
               end do
            end do
         end if
         do c = 1, nch
            do n = 1, n_known
               z = mu(cl(c), cm(c), n)
               if (herm) z = 0.5_rp*(z + conjg(mu(cm(c), cl(c), n)))
               xr(n, c) = merge(aimag(z), real(z, rp), cim(c))
            end do
         end do
      end subroutine load_channels

      !> Write ext(:,k) to mu(:,:,n_known+k); a Hermitian block also gets its other triangle.
      subroutine store_extension()
         real(rp) :: er(nb, nb), ei(nb, nb)
         integer :: c, kk

         do kk = 1, nsteps
            er = 0.0_rp
            ei = 0.0_rp
            do c = 1, nch
               if (cim(c)) then
                  ei(cl(c), cm(c)) = ext(c, kk)
                  if (herm) ei(cm(c), cl(c)) = -ext(c, kk)
               else
                  er(cl(c), cm(c)) = ext(c, kk)
                  if (herm) er(cm(c), cl(c)) = ext(c, kk)
               end if
            end do
            mu(:, :, n_known + kk) = cmplx(er, ei, rp)
         end do
      end subroutine store_extension

   end subroutine cheb_lp_extend

end module cheb_lp_mod
