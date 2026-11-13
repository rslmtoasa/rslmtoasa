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
      integer :: nb, nb2, nch, p, nlo, n_ext, i, j, k, info, lwork
      integer, allocatable :: ipiv(:)
      real(rp) :: smax, ynorm, rnorm
      real(rp), allocatable :: xr(:, :), tg(:, :), dd(:, :), gg(:, :), rr(:, :), rcopy(:, :), rv(:), ev(:), t(:), a(:)
      real(rp), allocatable :: comp(:, :), wr(:), wi(:), vr(:, :), work(:), outc(:), amp(:)
      logical, allocatable :: unst(:)
      complex(rp), allocatable :: vc(:, :), bm(:, :), lam(:), v1(:), lamk(:), cv(:), oc(:)
      real(rp) :: dummy(1, 1)

      nb = size(mu, 1)
      nb2 = nb*nb
      nch = 2*nb2
      n_ext = size(mu, 3)
      p = n_known/p_div
      nlo = n_known/nlo_div
      if (p < 2 .or. n_ext <= n_known .or. nlo < p) &
         call g_logger%fatal('cheb_lp_extend: need n_known >= 8 and size(mu,3) > n_known', __FILE__, __LINE__)

      ! Real channels: c = (m-1)*nb + l for Re mu(l,m,:), plus nb2 for Im.
      allocate (xr(n_known, nch))
      call load_channels()

      ! Normal equations R a = -r over all channels and the window, from the lag Gram matrix T(a,b) = sum_c x_c(a) x_c(b)
      ! (positions a, b = 1..n_known). The window sums are differences along the diagonals of T:
      !   G(i,j) = sum_c sum_{n=nlo+1}^{n_known} x_c(n-i) x_c(n-j) = D(n_known-i, n_known-j) - D(nlo-i, nlo-j),
      !   D(a,b) = T(a,b) + D(a-1,b-1), D = 0 for a = 0 or b = 0.
      ! Then R(i,j) = G(i,j), r(i) = G(i,0) (i,j = 1..p), |y|^2 = G(0,0), and the fit residual is |y|^2 + 2 a.r + a.R.a.
      allocate (tg(n_known, n_known), dd(0:n_known, 0:n_known), gg(0:p, 0:p))
      call dsyrk('U', 'N', n_known, nch, 1.0_rp, xr, n_known, 0.0_rp, tg, n_known)
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

      ! x(n_known-1+k) = Re [V lam^k V^-1 s]_1, s = [x(n_known-1), ..., x(n_known-p)].
      allocate (v1(p), bm(p, nch), ipiv(p), lamk(p), cv(p), oc(nch), outc(nch))
      v1 = vc(1, :)
      do i = 1, p
         bm(i, :) = cmplx(xr(n_known - i + 1, :), 0.0_rp, rp)
      end do
      call zgesv(p, nch, vc, p, ipiv, bm, p, info)
      if (info /= 0) call g_logger%fatal('cheb_lp_extend: zgesv failed (singular eigenvector matrix)', __FILE__, __LINE__)
      allocate (amp(p))
      do i = 1, p
         amp(i) = abs(v1(i))*sum(abs(bm(i, :)))
      end do
      unstable_weight = sum(amp, mask=unst)/sum(amp)
      lamk = lam
      imag_max = 0.0_rp
      do k = 1, n_ext - n_known
         cv = v1*lamk
         call zgemv('T', p, nch, (1.0_rp, 0.0_rp), bm, p, cv, 1, (0.0_rp, 0.0_rp), oc, 1)
         outc = real(oc, rp)
         imag_max = max(imag_max, maxval(abs(aimag(oc))))
         mu(:, :, n_known + k) = reshape(cmplx(outc(1:nb2), outc(nb2 + 1:nch), rp), [nb, nb])
         lamk = lamk*lam
      end do

   contains

      subroutine load_channels()
         integer :: l, m
         do m = 1, nb
            do l = 1, nb
               xr(:, (m - 1)*nb + l) = real(mu(l, m, 1:n_known), rp)
               xr(:, nb2 + (m - 1)*nb + l) = aimag(mu(l, m, 1:n_known))
            end do
         end do
      end subroutine load_channels

   end subroutine cheb_lp_extend

end module cheb_lp_mod
