!------------------------------------------------------------------------------
! Optimised CPU block-Lanczos (Haydock) recursion kernels.
!
! Shadows the Chebyshev `fast` architecture: an optimised block-Lanczos engine
! producing the same A_n / B_n recursion coefficients as the legacy
! crecal_b/hop_b[_hoh] path. Two design choices give the speedup over legacy:
!
!  * Atom-major (nb, nb, kk) storage (like legacy) so the lightcone matvec reads
!    contiguous neighbour blocks.
!  * The per-step block contractions (A_n, B^2, orthonormalisation) are done as
!    SINGLE large contiguous GEMMs over the gathered active set (stacked
!    nb*irnum x nb buffers) instead of many small per-atom GEMMs.
!
! The lightcone (active set irlist(1:irnum)) is grown exactly as legacy
! (hop_b: one shell/step; hop_b_hoh: two, h-sweep then eeo-sweep), so results
! are bit-identical to a dense fast sweep (fp32) and match legacy to round-off
! (fp64). The nb x nb B-matrix square root (zheev) is ALWAYS fp64. Block Lanczos
! uses the RAW Hamiltonian (ee + lsham on-site, no (a,b) scaling).
!------------------------------------------------------------------------------
module haydock_fast_mod
   implicit none
   private
   public :: block_lanczos_fast

   integer, parameter :: sp = selected_real_kind(6, 37)
   integer, parameter :: rp = selected_real_kind(15, 307)

contains

   subroutine block_lanczos_fast(psi0, lld, a_b, b2_b, ee, hall, lsham, nn, iz, &
                                 kk, nb, nnmax, ntype, nmax, use_sp, &
                                 hoh, eeo, hallo, enim)
      integer, intent(in) :: lld, kk, nb, nnmax, ntype, nmax
      complex(rp), intent(in) :: psi0(nb, nb, kk)
      complex(rp), intent(out) :: a_b(nb, nb, lld), b2_b(nb, nb, lld)
      complex(rp), intent(in) :: ee(nb, nb, nnmax, ntype), hall(nb, nb, nnmax, *)
      complex(rp), intent(in) :: lsham(nb, nb, ntype)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      logical, intent(in) :: use_sp, hoh
      complex(rp), intent(in) :: eeo(nb, nb, nnmax, ntype), hallo(nb, nb, nnmax, *)
      complex(rp), intent(in) :: enim(nb, nb, ntype)

      if (use_sp) then
         call block_lanczos_sp(psi0, lld, a_b, b2_b, ee, hall, lsham, nn, iz, &
                               kk, nb, nnmax, ntype, nmax, hoh, eeo, hallo, enim)
      else
         call block_lanczos_dp(psi0, lld, a_b, b2_b, ee, hall, lsham, nn, iz, &
                               kk, nb, nnmax, ntype, nmax, hoh, eeo, hallo, enim)
      end if
   end subroutine block_lanczos_fast

   !> Hermitian square root of B^2 (nb x nb), always fp64.
   subroutine bsqrt_dp(b2, bmat, binv, nb)
      integer, intent(in) :: nb
      complex(rp), intent(in) :: b2(nb, nb)
      complex(rp), intent(out) :: bmat(nb, nb), binv(nb, nb)
      complex(rp), parameter :: cone = (1.0_rp, 0.0_rp), czero = (0.0_rp, 0.0_rp)
      complex(rp) :: u(nb, nb), lam(nb, nb), lami(nb, nb), dum(nb, nb)
      real(rp) :: ev(nb), rwork(3*nb - 2)
      complex(rp) :: zwork(2*nb - 1)
      integer :: i, info
      external :: zheev, zgemm
      u = b2
      call zheev('v', 'u', nb, u, nb, ev, zwork, 2*nb - 1, rwork, info)
      if (info /= 0) then
         bmat = czero; binv = czero; return
      end if
      lam = czero; lami = czero
      do i = 1, nb
         lam(i, i) = cmplx(sqrt(max(0.0_rp, ev(i))), 0.0_rp, kind=rp)
         if (abs(lam(i, i)) > 0.0_rp) then
            lami(i, i) = cone/lam(i, i)
         else
            lami(i, i) = czero
         end if
      end do
      call zgemm('N', 'N', nb, nb, nb, cone, u, nb, lam, nb, czero, dum, nb)
      call zgemm('N', 'C', nb, nb, nb, cone, dum, nb, u, nb, czero, bmat, nb)
      call zgemm('N', 'N', nb, nb, nb, cone, u, nb, lami, nb, czero, dum, nb)
      call zgemm('N', 'C', nb, nb, nb, cone, dum, nb, u, nb, czero, binv, nb)
   end subroutine bsqrt_dp

   !===========================================================================
   ! fp32 engine (atom-major blocks + gathered big-GEMM contractions)
   !===========================================================================
   subroutine block_lanczos_sp(psi0, lld, a_b, b2_b, ee, hall, lsham, nn, iz, &
                               kk, nb, nnmax, ntype, nmax, hoh, eeo, hallo, enim)
      integer, intent(in) :: lld, kk, nb, nnmax, ntype, nmax
      complex(rp), intent(in) :: psi0(nb, nb, kk)
      complex(rp), intent(out) :: a_b(nb, nb, lld), b2_b(nb, nb, lld)
      complex(rp), intent(in) :: ee(nb, nb, nnmax, ntype), hall(nb, nb, nnmax, *)
      complex(rp), intent(in) :: lsham(nb, nb, ntype)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      logical, intent(in) :: hoh
      complex(rp), intent(in) :: eeo(nb, nb, nnmax, ntype), hallo(nb, nb, nnmax, *)
      complex(rp), intent(in) :: enim(nb, nb, ntype)
      ! Folded fp32 operators (raw: no a,b scaling; lsham on shell-1 for he/ha)
      complex(sp), allocatable :: he(:, :, :, :), ha(:, :, :, :)
      complex(sp), allocatable :: oe(:, :, :, :), oa(:, :, :, :), ons(:, :, :)
      ! Atom-major Lanczos vectors and stacked active buffers
      complex(sp), allocatable :: psi(:, :, :), pmn(:, :, :), hpsi(:, :, :), wt(:, :, :)
      complex(sp), allocatable :: ps(:, :), hs(:, :), ms(:, :), psnew(:, :)
      complex(sp) :: an(nb, nb), sumb_sp(nb, nb), bmat_sp(nb, nb), binv_sp(nb, nb)
      complex(rp) :: sumb(nb, nb), bmat(nb, nb), binv(nb, nb)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp), cmone = (-1.0_sp, 0.0_sp)
      integer, allocatable :: izero(:), idum(:), irlist(:)
      integer :: irnum, k, l, t_, ll, n, m, lda
      external :: cgemm

      lda = nb*kk
      ! --- fp32 folded operators ---------------------------------------------
      allocate (he(nb, nb, nnmax, ntype))
      he = cmplx(real(ee, sp), aimag(ee), sp)
      do t_ = 1, ntype
         he(:, :, 1, t_) = he(:, :, 1, t_) + cmplx(real(lsham(:, :, t_), sp), aimag(lsham(:, :, t_)), sp)
      end do
      if (nmax > 0) then
         allocate (ha(nb, nb, nnmax, nmax))
         ha = cmplx(real(hall(:, :, :, 1:nmax), sp), aimag(hall(:, :, :, 1:nmax)), sp)
         do k = 1, nmax
            t_ = iz(k)
            ha(:, :, 1, k) = ha(:, :, 1, k) + cmplx(real(lsham(:, :, t_), sp), aimag(lsham(:, :, t_)), sp)
         end do
      else
         allocate (ha(nb, nb, nnmax, 1))
      end if
      if (hoh) then
         allocate (oe(nb, nb, nnmax, ntype), ons(nb, nb, ntype))
         oe = cmplx(real(eeo, sp), aimag(eeo), sp)
         ons = cmplx(real(enim, sp), aimag(enim), sp) + cmplx(real(lsham, sp), aimag(lsham), sp)
         if (nmax > 0) then
            allocate (oa(nb, nb, nnmax, nmax))
            oa = cmplx(real(hallo(:, :, :, 1:nmax), sp), aimag(hallo(:, :, :, 1:nmax)), sp)
         else
            allocate (oa(nb, nb, nnmax, 1))
         end if
      end if

      allocate (psi(nb, nb, kk), pmn(nb, nb, kk), hpsi(nb, nb, kk), wt(nb, nb, kk))
      allocate (ps(lda, nb), hs(lda, nb), ms(lda, nb), psnew(lda, nb))
      allocate (izero(kk), idum(kk), irlist(kk))

      psi = czero
      !$omp parallel do private(k) schedule(static)
      do k = 1, kk
         psi(:, :, k) = cmplx(real(psi0(:, :, k), sp), aimag(psi0(:, :, k)), sp)
      end do
      !$omp end parallel do

      izero = 0; irnum = 0
      do k = 1, kk
         if (any(psi0(:, :, k) /= (0.0_rp, 0.0_rp))) then
            izero(k) = 1; irnum = irnum + 1; irlist(irnum) = k
         end if
      end do

      pmn = czero
      sumb_sp = czero
      do l = 1, nb
         sumb_sp(l, l) = cone
      end do
      a_b = (0.0_rp, 0.0_rp); b2_b = (0.0_rp, 0.0_rp)

      do ll = 1, lld - 1
         call lc_matvec(hpsi)                  ! atom-major; grows active set
         ! gather active blocks (stacked nb*irnum x nb, contiguous columns)
         do n = 1, irnum
            ps(nb*(n - 1) + 1:nb*n, :) = psi(:, :, irlist(n))
            hs(nb*(n - 1) + 1:nb*n, :) = hpsi(:, :, irlist(n))
            ms(nb*(n - 1) + 1:nb*n, :) = pmn(:, :, irlist(n))
         end do
         m = nb*irnum
         ! A_n = Ps^H Hs
         call cgemm('C', 'N', nb, nb, m, cone, ps, lda, hs, lda, czero, an, nb)
         a_b(:, :, ll) = an
         b2_b(:, :, ll) = sumb_sp
         ! residual Ms = Hs - Ms - Ps*A_n
         ms(1:m, :) = hs(1:m, :) - ms(1:m, :)
         call cgemm('N', 'N', m, nb, nb, cmone, ps, lda, an, nb, cone, ms, lda)
         ! B^2 = Ms^H Ms
         call cgemm('C', 'N', nb, nb, m, cone, ms, lda, ms, lda, czero, sumb_sp, nb)
         sumb = sumb_sp
         call bsqrt_dp(sumb, bmat, binv, nb)
         bmat_sp = cmplx(real(bmat, sp), aimag(bmat), sp)
         binv_sp = cmplx(real(binv, sp), aimag(binv), sp)
         ! psi_{n+1} = Ms*Binv (-> psnew) ; pmn = Ps*B (-> hs reused)
         call cgemm('N', 'N', m, nb, nb, cone, ms, lda, binv_sp, nb, czero, psnew, lda)
         call cgemm('N', 'N', m, nb, nb, cone, ps, lda, bmat_sp, nb, czero, hs, lda)
         ! scatter back to atom-major
         do n = 1, irnum
            psi(:, :, irlist(n)) = psnew(nb*(n - 1) + 1:nb*n, :)
            pmn(:, :, irlist(n)) = hs(nb*(n - 1) + 1:nb*n, :)
         end do
      end do
      b2_b(:, :, lld) = sumb_sp

      deallocate (he, ha, psi, pmn, hpsi, wt, ps, hs, ms, psnew, izero, idum, irlist)
      if (allocated(oe)) deallocate (oe, oa, ons)

   contains

      !> Lightcone block matvec y = H psi (fp32, atom-major), growing the active
      !> set. Mirrors hop_b (non-hoh) / hop_b_hoh (h-sweep then eeo-sweep).
      subroutine lc_matvec(y)
         complex(sp), intent(out) :: y(nb, nb, kk)
         complex(sp) :: acc(nb, nb)
         integer :: kk_, s_, nbr, nr, t2, n
         y = czero
         ! --- sweep A: y = (ee+lsham) psi -----------------------------------
         !$omp parallel do private(kk_, s_, nbr, nr, t2, acc) schedule(dynamic, 32)
         do kk_ = 1, kk
            idum(kk_) = izero(kk_)
            t2 = iz(kk_); nr = nn(kk_, 1); acc = czero
            if (izero(kk_) /= 0) then
               if (kk_ <= nmax) then
                  call cgemm('N', 'N', nb, nb, nb, cone, ha(:, :, 1, kk_), nb, psi(:, :, kk_), nb, cone, acc, nb)
               else
                  call cgemm('N', 'N', nb, nb, nb, cone, he(:, :, 1, t2), nb, psi(:, :, kk_), nb, cone, acc, nb)
               end if
            end if
            do s_ = 2, nr
               nbr = nn(kk_, s_); if (nbr == 0) cycle
               if (izero(nbr) /= 0) then
                  if (kk_ <= nmax) then
                     call cgemm('N', 'N', nb, nb, nb, cone, ha(:, :, s_, kk_), nb, psi(:, :, nbr), nb, cone, acc, nb)
                  else
                     call cgemm('N', 'N', nb, nb, nb, cone, he(:, :, s_, t2), nb, psi(:, :, nbr), nb, cone, acc, nb)
                  end if
                  idum(kk_) = 1
               end if
            end do
            y(:, :, kk_) = acc
         end do
         !$omp end parallel do
         call rebuild_active()

         if (hoh) then
            wt = czero
            !$omp parallel do private(kk_, s_, nbr, nr, t2, acc) schedule(dynamic, 32)
            do kk_ = 1, kk
               idum(kk_) = izero(kk_)
               t2 = iz(kk_); nr = nn(kk_, 1); acc = czero
               if (izero(kk_) /= 0) then
                  if (kk_ <= nmax) then
                     call cgemm('N', 'N', nb, nb, nb, cone, oa(:, :, 1, kk_), nb, y(:, :, kk_), nb, cone, acc, nb)
                  else
                     call cgemm('N', 'N', nb, nb, nb, cone, oe(:, :, 1, t2), nb, y(:, :, kk_), nb, cone, acc, nb)
                  end if
               end if
               do s_ = 2, nr
                  nbr = nn(kk_, s_); if (nbr == 0) cycle
                  if (izero(nbr) /= 0) then
                     if (kk_ <= nmax) then
                        call cgemm('N', 'N', nb, nb, nb, cone, oa(:, :, s_, kk_), nb, y(:, :, nbr), nb, cone, acc, nb)
                     else
                        call cgemm('N', 'N', nb, nb, nb, cone, oe(:, :, s_, t2), nb, y(:, :, nbr), nb, cone, acc, nb)
                     end if
                     idum(kk_) = 1
                  end if
               end do
               wt(:, :, kk_) = acc
            end do
            !$omp end parallel do
            call rebuild_active()
            ! combine over active: y = h*psi - eeo*(h*psi) + (enim+lsham)*psi
            !$omp parallel do private(kk_, t2, acc, n) schedule(dynamic, 32)
            do n = 1, irnum
               kk_ = irlist(n); t2 = iz(kk_)
               acc = y(:, :, kk_) - wt(:, :, kk_)
               if (izero(kk_) /= 0) &
                  call cgemm('N', 'N', nb, nb, nb, cone, ons(:, :, t2), nb, psi(:, :, kk_), nb, cone, acc, nb)
               y(:, :, kk_) = acc
            end do
            !$omp end parallel do
         end if
      end subroutine lc_matvec

      subroutine rebuild_active()
         integer :: kk2
         irnum = 0
         do kk2 = 1, kk
            izero(kk2) = idum(kk2)
            if (izero(kk2) /= 0) then
               irnum = irnum + 1; irlist(irnum) = kk2
            end if
         end do
      end subroutine rebuild_active

   end subroutine block_lanczos_sp

   !===========================================================================
   ! fp64 engine (same structure, raw rp operators)
   !===========================================================================
   subroutine block_lanczos_dp(psi0, lld, a_b, b2_b, ee, hall, lsham, nn, iz, &
                               kk, nb, nnmax, ntype, nmax, hoh, eeo, hallo, enim)
      integer, intent(in) :: lld, kk, nb, nnmax, ntype, nmax
      complex(rp), intent(in) :: psi0(nb, nb, kk)
      complex(rp), intent(out) :: a_b(nb, nb, lld), b2_b(nb, nb, lld)
      complex(rp), intent(in) :: ee(nb, nb, nnmax, ntype), hall(nb, nb, nnmax, *)
      complex(rp), intent(in) :: lsham(nb, nb, ntype)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      logical, intent(in) :: hoh
      complex(rp), intent(in) :: eeo(nb, nb, nnmax, ntype), hallo(nb, nb, nnmax, *)
      complex(rp), intent(in) :: enim(nb, nb, ntype)
      complex(rp), allocatable :: psi(:, :, :), pmn(:, :, :), hpsi(:, :, :), wt(:, :, :)
      complex(rp), allocatable :: ps(:, :), hs(:, :), ms(:, :), psnew(:, :)
      complex(rp) :: an(nb, nb), sumb(nb, nb), bmat(nb, nb), binv(nb, nb)
      complex(rp), parameter :: cone = (1.0_rp, 0.0_rp), czero = (0.0_rp, 0.0_rp), cmone = (-1.0_rp, 0.0_rp)
      integer, allocatable :: izero(:), idum(:), irlist(:)
      integer :: irnum, k, l, ll, n, m, lda
      external :: zgemm

      lda = nb*kk
      allocate (psi(nb, nb, kk), pmn(nb, nb, kk), hpsi(nb, nb, kk), wt(nb, nb, kk))
      allocate (ps(lda, nb), hs(lda, nb), ms(lda, nb), psnew(lda, nb))
      allocate (izero(kk), idum(kk), irlist(kk))

      psi = czero
      !$omp parallel do private(k) schedule(static)
      do k = 1, kk
         psi(:, :, k) = psi0(:, :, k)
      end do
      !$omp end parallel do

      izero = 0; irnum = 0
      do k = 1, kk
         if (any(psi0(:, :, k) /= czero)) then
            izero(k) = 1; irnum = irnum + 1; irlist(irnum) = k
         end if
      end do

      pmn = czero
      sumb = czero
      do l = 1, nb
         sumb(l, l) = cone
      end do
      a_b = czero; b2_b = czero

      do ll = 1, lld - 1
         call lc_matvec_dp(hpsi)
         do n = 1, irnum
            ps(nb*(n - 1) + 1:nb*n, :) = psi(:, :, irlist(n))
            hs(nb*(n - 1) + 1:nb*n, :) = hpsi(:, :, irlist(n))
            ms(nb*(n - 1) + 1:nb*n, :) = pmn(:, :, irlist(n))
         end do
         m = nb*irnum
         call zgemm('C', 'N', nb, nb, m, cone, ps, lda, hs, lda, czero, an, nb)
         a_b(:, :, ll) = an
         b2_b(:, :, ll) = sumb
         ms(1:m, :) = hs(1:m, :) - ms(1:m, :)
         call zgemm('N', 'N', m, nb, nb, cmone, ps, lda, an, nb, cone, ms, lda)
         call zgemm('C', 'N', nb, nb, m, cone, ms, lda, ms, lda, czero, sumb, nb)
         call bsqrt_dp(sumb, bmat, binv, nb)
         call zgemm('N', 'N', m, nb, nb, cone, ms, lda, binv, nb, czero, psnew, lda)
         call zgemm('N', 'N', m, nb, nb, cone, ps, lda, bmat, nb, czero, hs, lda)
         do n = 1, irnum
            psi(:, :, irlist(n)) = psnew(nb*(n - 1) + 1:nb*n, :)
            pmn(:, :, irlist(n)) = hs(nb*(n - 1) + 1:nb*n, :)
         end do
      end do
      b2_b(:, :, lld) = sumb

      deallocate (psi, pmn, hpsi, wt, ps, hs, ms, psnew, izero, idum, irlist)

   contains

      subroutine lc_matvec_dp(y)
         complex(rp), intent(out) :: y(nb, nb, kk)
         complex(rp) :: acc(nb, nb), locham(nb, nb)
         integer :: kk_, s_, nbr, nr, t2, n
         y = czero
         !$omp parallel do private(kk_, s_, nbr, nr, t2, acc, locham) schedule(dynamic, 32)
         do kk_ = 1, kk
            idum(kk_) = izero(kk_)
            t2 = iz(kk_); nr = nn(kk_, 1); acc = czero
            if (izero(kk_) /= 0) then
               if (kk_ <= nmax) then
                  locham = hall(:, :, 1, kk_) + lsham(:, :, t2)
               else
                  locham = ee(:, :, 1, t2) + lsham(:, :, t2)
               end if
               call zgemm('N', 'N', nb, nb, nb, cone, locham, nb, psi(:, :, kk_), nb, cone, acc, nb)
            end if
            do s_ = 2, nr
               nbr = nn(kk_, s_); if (nbr == 0) cycle
               if (izero(nbr) /= 0) then
                  if (kk_ <= nmax) then
                     call zgemm('N', 'N', nb, nb, nb, cone, hall(:, :, s_, kk_), nb, psi(:, :, nbr), nb, cone, acc, nb)
                  else
                     call zgemm('N', 'N', nb, nb, nb, cone, ee(:, :, s_, t2), nb, psi(:, :, nbr), nb, cone, acc, nb)
                  end if
                  idum(kk_) = 1
               end if
            end do
            y(:, :, kk_) = acc
         end do
         !$omp end parallel do
         call rebuild_active_dp()

         if (hoh) then
            wt = czero
            !$omp parallel do private(kk_, s_, nbr, nr, t2, acc) schedule(dynamic, 32)
            do kk_ = 1, kk
               idum(kk_) = izero(kk_)
               t2 = iz(kk_); nr = nn(kk_, 1); acc = czero
               if (izero(kk_) /= 0) then
                  if (kk_ <= nmax) then
                     call zgemm('N', 'N', nb, nb, nb, cone, hallo(:, :, 1, kk_), nb, y(:, :, kk_), nb, cone, acc, nb)
                  else
                     call zgemm('N', 'N', nb, nb, nb, cone, eeo(:, :, 1, t2), nb, y(:, :, kk_), nb, cone, acc, nb)
                  end if
               end if
               do s_ = 2, nr
                  nbr = nn(kk_, s_); if (nbr == 0) cycle
                  if (izero(nbr) /= 0) then
                     if (kk_ <= nmax) then
                        call zgemm('N', 'N', nb, nb, nb, cone, hallo(:, :, s_, kk_), nb, y(:, :, nbr), nb, cone, acc, nb)
                     else
                        call zgemm('N', 'N', nb, nb, nb, cone, eeo(:, :, s_, t2), nb, y(:, :, nbr), nb, cone, acc, nb)
                     end if
                     idum(kk_) = 1
                  end if
               end do
               wt(:, :, kk_) = acc
            end do
            !$omp end parallel do
            call rebuild_active_dp()
            !$omp parallel do private(kk_, t2, acc, n) schedule(dynamic, 32)
            do n = 1, irnum
               kk_ = irlist(n); t2 = iz(kk_)
               acc = y(:, :, kk_) - wt(:, :, kk_)
               if (izero(kk_) /= 0) &
                  call zgemm('N', 'N', nb, nb, nb, cone, enim(:, :, t2) + lsham(:, :, t2), nb, psi(:, :, kk_), nb, cone, acc, nb)
               y(:, :, kk_) = acc
            end do
            !$omp end parallel do
         end if
      end subroutine lc_matvec_dp

      subroutine rebuild_active_dp()
         integer :: kk2
         irnum = 0
         do kk2 = 1, kk
            izero(kk2) = idum(kk2)
            if (izero(kk2) /= 0) then
               irnum = irnum + 1; irlist(irnum) = kk2
            end if
         end do
      end subroutine rebuild_active_dp

   end subroutine block_lanczos_dp

end module haydock_fast_mod
