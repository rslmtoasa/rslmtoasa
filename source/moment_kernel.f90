!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: moment_kernel
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Exact double-Chebyshev moment generator for the KPM conductivity pipeline
!> (milestone B5, task B5.1). Dependency-free (only `precision_mod`) so the
!> eigenbasis-identity known-answer test needs no LMTO machinery.
!>
!> The Kubo/KPM conductivity consumes the velocity-sandwiched DOUBLE Chebyshev
!> moment the recursion route builds in `recursion_transport::compute_moments_
!> stochastic`,
!>
!>   mu_nm(l1,l2) = <r| T_m(H~) v_a T_n(H~) v_b |r>,
!>
!> where H~ = (H - b)/a is the affinely rescaled Hamiltonian (the SAME window
!> `resolve_chebyshev_window` fixes for the recursion route and `gamma_nm`
!> uses), T_p is the Chebyshev polynomial of the first kind, v_a/v_b are the
!> velocity operators, and |r> is the reference block localized on one
!> unit-cell site (the "per_type" reference of the recursion route). Index n is
!> the INNER Chebyshev factor (between the two velocities), m the OUTER (leftmost)
!> factor -- exactly the `mu_nm_stochastic(:,:,n,m,:)` storage order.
!>
!> This kernel produces the SAME moment EXACTLY from the k-space eigenpairs. In
!> the eigenbasis of H(k) = U diag(eps) U^dagger,
!>
!>   T_p(H~(k)) = U diag(T_p(eps~)) U^dagger,   eps~ = (eps - b)/a,
!>
!> so, transforming the velocity operators to the eigenbasis
!> (v~_a = U^dagger v_a(k) U, likewise v~_b), the on-site block for the
!> unit-cell site whose orbitals occupy rows `ioff+1..ioff+nblk` is
!>
!>   mu_nm(l1,l2) = (1/N_k) sum_k [ U_s D_m v~_a D_n v~_b U_s^dagger ]_{l1,l2},
!>
!> with D_p = diag(T_p(eps~_alpha)) and U_s the site-block rows of U. The k-space
!> velocity operator v_a(k) = sum_R v_a(R) e^{i k.R} is the Fourier transform of
!> the SAME real-space velocity blocks the recursion route applies, so the two
!> routes evaluate one operator in two representations: they agree in the joint
!> large-cluster / dense-mesh limit, and the residual max|mu_exact - mu_rec| is
!> the direct KPM error bound (task B5.1 acceptance).
!>
!> Rescaling convention (do-not-get-this-wrong, B5.1): moments live on the scaled
!> spectrum; the caller passes the SAME a, b the recursion route/`gamma_nm` use
!> (a = (e_max - e_min)/(2 - 0.3), b = (e_max + e_min)/2 from the Chebyshev
!> window). Only H is scaled -- the velocity operators are applied unscaled,
!> matching `compute_moments_stochastic`.
!------------------------------------------------------------------------------
module moment_kernel_mod
   use precision_mod, only: rp
   implicit none
   private

   public :: moment_onsite_block
   public :: chebyshev_scaled

contains

   !> @brief Chebyshev polynomials T_0..T_{ll-1} of scaled scalars.
   !> @details T(:,1)=1, T(:,2)=x, T(:,p)=2 x T(:,p-1)-T(:,p-2), where the input
   !>          values are affinely scaled x = (e - b)/a. Column p holds T_{p-1},
   !>          matching the recursion route's 1-based moment index (m=1 -> T_0).
   !> @param[in]  e   Unscaled eigenvalues, shape (n).
   !> @param[in]  a   Scale of the affine map H -> (H-b)/a.
   !> @param[in]  b   Shift of the affine map.
   !> @param[in]  ll  Number of Chebyshev polynomials (cond_ll).
   !> @param[out] t   Chebyshev table, shape (n, ll); t(i,p) = T_{p-1}((e(i)-b)/a).
   pure subroutine chebyshev_scaled(e, a, b, ll, t)
      real(rp), intent(in) :: e(:)
      real(rp), intent(in) :: a, b
      integer, intent(in) :: ll
      real(rp), intent(out) :: t(:, :)
      integer :: n, p
      real(rp), allocatable :: x(:)

      n = size(e)
      allocate (x(n))
      x(:) = (e(:) - b)/a

      if (ll >= 1) t(:, 1) = 1.0_rp
      if (ll >= 2) t(:, 2) = x(:)
      do p = 3, ll
         t(:, p) = 2.0_rp*x(:)*t(:, p - 1) - t(:, p - 2)
      end do
      deallocate (x)
   end subroutine chebyshev_scaled

   !> @brief Accumulate the exact on-site double-Chebyshev moment block.
   !> @details Streams over all k-points, transforms the velocity operators into
   !>          the per-k eigenbasis, and forms the SAME on-site moment block the
   !>          recursion route builds (see the module header for the derivation
   !>          and the n/m ordering). `mu` is zeroed on entry, accumulated over k,
   !>          and normalized by 1/N_k on exit. The velocities enter unscaled;
   !>          only the eigenvalues are affinely rescaled by (a, b).
   !> @param[in]  eigenvalues  Real eigenvalues, shape (nmat, nk).
   !> @param[in]  eigenvectors Column eigenvectors U(:,alpha,ik), shape (nmat,nmat,nk).
   !> @param[in]  va_k         k-space velocity operator v_a(k), shape (nmat,nmat,nk).
   !> @param[in]  vb_k         k-space velocity operator v_b(k), shape (nmat,nmat,nk).
   !> @param[in]  a            Chebyshev window scale (same as recursion route).
   !> @param[in]  b            Chebyshev window shift (same as recursion route).
   !> @param[in]  ioff         Zero-based row offset of the reference site's block.
   !> @param[in]  nblk         Reference block size (nb = 2*norb).
   !> @param[in]  cond_ll      Number of Chebyshev moments per index.
   !> @param[out] mu           On-site moment block, shape (nblk,nblk,cond_ll,cond_ll);
   !>                          mu(:,:,n,m) = <site| T_m(H~) v_a T_n(H~) v_b |site>.
   subroutine moment_onsite_block(eigenvalues, eigenvectors, va_k, vb_k, a, b, &
                                  ioff, nblk, cond_ll, mu)
      real(rp), intent(in) :: eigenvalues(:, :)
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      complex(rp), intent(in) :: va_k(:, :, :)
      complex(rp), intent(in) :: vb_k(:, :, :)
      real(rp), intent(in) :: a, b
      integer, intent(in) :: ioff, nblk, cond_ll
      complex(rp), intent(out) :: mu(:, :, :, :)

      integer :: nmat, nk, ik, n, m, alpha
      real(rp), allocatable :: t(:, :)                 ! (nmat, cond_ll) T_{p-1}(eps~)
      complex(rp), allocatable :: u(:, :), udag(:, :)  ! (nmat, nmat)
      complex(rp), allocatable :: va_e(:, :), vb_e(:, :)
      complex(rp), allocatable :: us(:, :), bs(:, :)   ! (nblk, nmat), (nmat, nblk)
      complex(rp), allocatable :: utilde(:, :, :)      ! (nblk, nmat, cond_ll)
      complex(rp), allocatable :: btilde(:, :), wn(:, :)

      nmat = size(eigenvalues, 1)
      nk = size(eigenvalues, 2)

      mu = (0.0_rp, 0.0_rp)

      allocate (t(nmat, cond_ll))
      allocate (u(nmat, nmat), udag(nmat, nmat))
      allocate (va_e(nmat, nmat), vb_e(nmat, nmat))
      allocate (us(nblk, nmat), bs(nmat, nblk))
      allocate (utilde(nblk, nmat, cond_ll))
      allocate (btilde(nmat, nblk), wn(nmat, nblk))

      do ik = 1, nk
         call chebyshev_scaled(eigenvalues(:, ik), a, b, cond_ll, t)

         u(:, :) = eigenvectors(:, :, ik)
         udag(:, :) = conjg(transpose(u))

         ! Velocity operators in the per-k eigenbasis: v~ = U^dagger v(k) U.
         va_e(:, :) = matmul(udag, matmul(va_k(:, :, ik), u))
         vb_e(:, :) = matmul(udag, matmul(vb_k(:, :, ik), u))

         ! Reference-site block rows of U and B_s = v~_b U_s^dagger.
         us(:, :) = u(ioff + 1:ioff + nblk, :)
         bs(:, :) = matmul(vb_e, conjg(transpose(us)))

         ! U~_m[l1,alpha] = U_s[l1,alpha] * T_{m-1}(eps~_alpha) (column scaling).
         do m = 1, cond_ll
            do alpha = 1, nmat
               utilde(:, alpha, m) = us(:, alpha)*t(alpha, m)
            end do
         end do

         do n = 1, cond_ll
            ! B~_n[beta,l2] = T_{n-1}(eps~_beta) * B_s[beta,l2]; W_n = v~_a B~_n.
            do alpha = 1, nmat
               btilde(alpha, :) = t(alpha, n)*bs(alpha, :)
            end do
            wn(:, :) = matmul(va_e, btilde)
            do m = 1, cond_ll
               mu(:, :, n, m) = mu(:, :, n, m) + matmul(utilde(:, :, m), wn)
            end do
         end do
      end do

      mu = mu/real(nk, rp)

      deallocate (t, u, udag, va_e, vb_e, us, bs, utilde, btilde, wn)
   end subroutine moment_onsite_block

end module moment_kernel_mod
