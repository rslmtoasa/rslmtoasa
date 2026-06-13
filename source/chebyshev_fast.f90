!------------------------------------------------------------------------------
! Optimized CPU replacements for chebyshev_recur (recursion_mod) and
! chebyshev_green (bands_mod), applying the validated KPM tricks:
!
!  1. The window scaling (H-b)/a is folded ONCE into local copies of
!     ee/hall (lsham included), so the inner loop is psi2 = 2*Ht*psi1 - psi0
!     with no scale/shift passes.
!  2. psi is packed SITE-MAJOR as an (nb*kk, nb) matrix, so every block
!     moment sum_k psi^H phi is one GEMM ('C','N', nb, nb, nb*kk) -- and
!     mu(2ll+1) uses HERK (Hermitian, half flops). No per-atom scalar loop.
!  3. The matvec is one fused OpenMP sweep over atoms; the lightcone
!     (izero/idum) is dropped: zero psi contributes zero, results identical.
!  4. DOS/GF reconstruction is ONE GEMM per atom: g0_n = M_n * F, with the
!     Jackson kernel, x2 coefficients, -i*exp(-i*n*acos(x)) phase and
!     1/sqrt(a^2-(E-b)^2) prefactor folded into the transfer matrix F.
!     (Replaces the scalar exp/acos triple loop and the unused t_polynomial.)
!
! Drop-in: same mu_n ordering (2*lld+2 block moments) and same g0 output as
! the existing routines. hoh path not handled -- guard with
! if (.not. this%hamiltonian%hoh). Keeps low-precision mirrors for the
! Hamiltonian inputs and DOS-construction moments local to this module, while
! preserving double-precision interfaces to the rest of the code.
!------------------------------------------------------------------------------
module chebyshev_fast_mod
   implicit none
   private
   public :: cheb_moments_fast, cheb_green_fast
   integer, parameter :: sp = selected_real_kind(6, 37)
   integer, parameter :: rp = selected_real_kind(15, 307)
   real(rp), parameter :: pi = 3.14159265358979323846_rp

contains

   !> Block Chebyshev moments for one starting state.
   !> psi0  : (nb, nb, kk) starting block state (site or ij sign combos)
   !> ee    : (nb, nb, nnmax, ntype), hall: (nb, nb, nnmax, nmax) (or nmax=0)
   !> lsham : (nb, nb, ntype); nn: (kk, nnmax) with nn(k,1)=count; iz: (kk)
   !> mu    : (nb, nb, 2*lld+2) out, identical ordering to mu_n
   subroutine cheb_moments_fast(psi0, ee, hall, lsham, nn, iz, kk, nb, &
                                nnmax, ntype, nmax, lld, a, b, mu)
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax, lld
      complex(rp), intent(in) :: psi0(nb, nb, kk), ee(nb, nb, nnmax, ntype)
      complex(rp), intent(in) :: hall(nb, nb, nnmax, *), lsham(nb, nb, ntype)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      real(rp), intent(in) :: a, b
      complex(rp), intent(out) :: mu(nb, nb, 2*lld + 2)
      ! Locals
      complex(sp), allocatable :: hee(:, :, :, :), hha(:, :, :, :)
      complex(sp), allocatable :: p0(:, :), p1(:, :), p2(:, :)
      complex(sp) :: mu1_sp(nb, nb), mu2_sp(nb, nb), dum_sp(nb, nb)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp)
      integer :: k, l, c, ld, ll
      real(sp) :: inva, a_sp, b_sp

      ld = nb*kk
      a_sp = real(a, sp)
      b_sp = real(b, sp)
      inva = 1.0_sp/a_sp
      allocate (p0(ld, nb), p1(ld, nb), p2(ld, nb))

      ! --- 1. scaled operator copies: Ht = (H + lsham - b*I)/a -------------
      call prepare_scaled_hamiltonian_sp(ee, hall, lsham, iz, nb, nnmax, ntype, nmax, inva, b_sp, hee, hha)

      ! --- 2. pack psi0 site-major: p(l + nb*(k-1), c) ---------------------
      !$omp parallel do private(k, c, l) schedule(static)
      do k = 1, kk
         do c = 1, nb
            do l = 1, nb
               p0(l + nb*(k - 1), c) = cmplx(real(psi0(l, c, k), sp), aimag(psi0(l, c, k)), sp)
            end do
         end do
      end do
      !$omp end parallel do

      ! --- 3. moments: mu1 = p0^H p0, p1 = Ht p0, mu2 = p0^H p1 ------------
      call cherk_full(nb, ld, p0, mu1_sp)
      mu(:, :, 1) = mu1_sp
      call apply_step(p0, p0, p1, 1.0_sp, 0.0_sp)
      call cgemm('C', 'N', nb, nb, ld, cone, p0, ld, p1, ld, czero, mu2_sp, nb)
      mu(:, :, 2) = mu2_sp

      do ll = 1, lld
         call apply_step(p1, p0, p2, 2.0_sp, -1.0_sp)   ! p2 = 2 Ht p1 - p0
         call cherk_full(nb, ld, p1, dum_sp)            ! dum1 = p1^H p1
         mu(:, :, 2*ll + 1) = 2.0_sp*dum_sp - mu1_sp
         call cgemm('C', 'N', nb, nb, ld, cone, p2, ld, p1, ld, czero, dum_sp, nb)
         mu(:, :, 2*ll + 2) = 2.0_sp*dum_sp - mu2_sp
         call swapm(p0, p1); call swapm(p1, p2)         ! rotate registers
      end do
      deallocate (p0, p1, p2, hee)
      if (allocated(hha)) deallocate (hha)

   contains

      subroutine prepare_scaled_hamiltonian_sp(ee_in, hall_in, lsham_in, iz_in, nb_in, nnmax_in, ntype_in, nmax_in, inva_in, b_in, hee_out, hha_out)
         complex(rp), intent(in) :: ee_in(nb_in, nb_in, nnmax_in, ntype_in)
         complex(rp), intent(in) :: hall_in(nb_in, nb_in, nnmax_in, *)
         complex(rp), intent(in) :: lsham_in(nb_in, nb_in, ntype_in)
         integer, intent(in) :: iz_in(kk), nb_in, nnmax_in, ntype_in, nmax_in
         real(sp), intent(in) :: inva_in, b_in
         complex(sp), allocatable, intent(out) :: hee_out(:, :, :, :), hha_out(:, :, :, :)
         integer :: t_in, k_in, l_in

         allocate (hee_out(nb_in, nb_in, nnmax_in, ntype_in))
         hee_out = cmplx(real(ee_in, sp), aimag(ee_in), sp)*inva_in
         do t_in = 1, ntype_in
            hee_out(:, :, 1, t_in) = hee_out(:, :, 1, t_in) + cmplx(real(lsham_in(:, :, t_in), sp), aimag(lsham_in(:, :, t_in)), sp)*inva_in
            do l_in = 1, nb_in
               hee_out(l_in, l_in, 1, t_in) = hee_out(l_in, l_in, 1, t_in) - b_in*inva_in
            end do
         end do

         if (nmax_in > 0) then
            allocate (hha_out(nb_in, nb_in, nnmax_in, nmax_in))
            hha_out = cmplx(real(hall_in(:, :, :, 1:nmax_in), sp), aimag(hall_in(:, :, :, 1:nmax_in)), sp)*inva_in
            do k_in = 1, nmax_in
               t_in = iz_in(k_in)
               hha_out(:, :, 1, k_in) = hha_out(:, :, 1, k_in) + cmplx(real(lsham_in(:, :, t_in), sp), aimag(lsham_in(:, :, t_in)), sp)*inva_in
               do l_in = 1, nb_in
                  hha_out(l_in, l_in, 1, k_in) = hha_out(l_in, l_in, 1, k_in) - b_in*inva_in
               end do
            end do
         end if
      end subroutine prepare_scaled_hamiltonian_sp

      subroutine apply_step_manual(x1, x0, y, alpha, beta)
         complex(sp), intent(in)  :: x1(ld, nb), x0(ld, nb)
         complex(sp), intent(out) :: y(ld, nb)
         real(sp), intent(in)     :: alpha, beta
         complex(sp) :: acc(nb, nb), xm
         integer :: kk_, t_, s_, nbr, nr, r0, m, cc, l
         !$omp parallel do private(kk_,t_,s_,nbr,nr,r0,m,cc,l,xm,acc) schedule(dynamic,32)
         do kk_ = 1, kk
            acc = (0.0_sp, 0.0_sp)
            nr = nn(kk_, 1)
            t_ = iz(kk_)
            do s_ = 1, nr
               if (s_ == 1) then
                  nbr = kk_
               else
                  nbr = nn(kk_, s_); if (nbr == 0) cycle
               end if
               r0 = nb*(nbr - 1)
               if (kk_ <= nmax) then          ! acc += hha(:,:,s_,kk_) * x1_block
                  do cc = 1, nb
                     do m = 1, nb
                        xm = x1(r0 + m, cc)
                        do l = 1, nb          ! contiguous in l -> vectorizes
                           acc(l, cc) = acc(l, cc) + hha(l, m, s_, kk_)*xm
                        end do
                     end do
                  end do
               else                            ! acc += hee(:,:,s_,t_) * x1_block
                  do cc = 1, nb
                     do m = 1, nb
                        xm = x1(r0 + m, cc)
                        do l = 1, nb
                           acc(l, cc) = acc(l, cc) + hee(l, m, s_, t_)*xm
                        end do
                     end do
                  end do
               end if
            end do
            r0 = nb*(kk_ - 1)
            if (beta /= 0.0_sp) then
               y(r0+1:r0+nb, :) = alpha*acc + beta*x0(r0+1:r0+nb, :)
            else
               y(r0+1:r0+nb, :) = alpha*acc
            end if
         end do
         !$omp end parallel do
      end subroutine apply_step_manual

      !> y = alpha*(Ht x1) + beta*x0, one fused sweep (site-major arrays)
      subroutine apply_step(x1, x0, y, alpha, beta)
         complex(sp), intent(in) :: x1(ld, nb), x0(ld, nb)
         complex(sp), intent(out) :: y(ld, nb)
         real(sp), intent(in) :: alpha, beta
         complex(sp) :: acc(nb, nb)
         integer :: kk_, t_, s_, nbr, nr, r0
         !$omp parallel do private(kk_, t_, s_, nbr, nr, r0, acc) schedule(dynamic, 32)
         do kk_ = 1, kk
            acc = (0.0_sp, 0.0_sp)
            nr = nn(kk_, 1)
            do s_ = 1, nr
               if (s_ == 1) then
                  nbr = kk_
               else
                  nbr = nn(kk_, s_)
                  if (nbr == 0) cycle
               end if
               r0 = nb*(nbr - 1)
               if (kk_ <= nmax) then
                  ! acc += hha(:,:,s_,kk_) * x1_block(nbr)
                  call cgemm('N', 'N', nb, nb, nb, cone, hha(:, :, s_, kk_), nb, x1(r0 + 1, 1), ld, cone, acc, nb)
               else
                  t_ = iz(kk_)
                  call cgemm('N', 'N', nb, nb, nb, cone, hee(:, :, s_, t_), nb, x1(r0 + 1, 1), ld, cone, acc, nb)
               end if
            end do
            r0 = nb*(kk_ - 1)
            if (beta /= 0.0_sp) then
               y(r0 + 1:r0 + nb, :) = alpha*acc + beta*x0(r0 + 1:r0 + nb, :)
            else
               y(r0 + 1:r0 + nb, :) = alpha*acc
            end if
         end do
         !$omp end parallel do
      end subroutine apply_step

      !> Hermitian rank-k: C = A^H A via cherk, then fill the upper triangle
      subroutine cherk_full(n, k, Amat, C)
         integer, intent(in) :: n, k
         complex(sp), intent(in) :: Amat(k, n)
         complex(sp), intent(out) :: C(n, n)
         integer :: i, j
         C = (0.0_sp, 0.0_sp)
         call cherk('L', 'C', n, k, 1.0_sp, Amat, k, 0.0_sp, C, n)
         do j = 2, n
            do i = 1, j - 1
               C(i, j) = conjg(C(j, i))
            end do
         end do
      end subroutine cherk_full

      subroutine swapm(x, y)
         complex(sp), intent(inout) :: x(ld, nb), y(ld, nb)
         complex(sp) :: tmp
         integer :: i, j
         !$omp parallel do private(i, j, tmp) schedule(static)
         do j = 1, nb
            do i = 1, ld
               tmp = x(i, j); x(i, j) = y(i, j); y(i, j) = tmp
            end do
         end do
         !$omp end parallel do
      end subroutine swapm

   end subroutine cheb_moments_fast

   !> Green-function / DOS reconstruction, replacing the per-atom body of
   !> chebyshev_green: g0(:,:,ie,n) = sum_i mu(:,:,i,n)*F(i,ie). One GEMM
   !> per atom; Jackson kernel in the exact bands.f90 convention.
   !> mu : (nb, nb, n_mom, natoms) local slice; ene: (nv) inside the window
   !> g0 : (nb, nb, nv, natoms) out;  a, b from resolve_chebyshev_window
   subroutine cheb_green_fast(mu, nb, n_mom, natoms, ene, nv, a, b, g0)
      integer, intent(in) :: nb, n_mom, natoms, nv
      complex(rp), intent(in) :: mu(nb, nb, n_mom, natoms)
      real(rp), intent(in) :: ene(nv), a, b
      complex(rp), intent(out) :: g0(nb, nb, nv, natoms)
      complex(sp), allocatable :: F(:, :)
      complex(sp), allocatable :: mu_sp(:, :, :, :), g0_atom(:, :, :)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp)
      real(rp) :: th, thl, gj, cc, pref, ang, cot
      integer :: i, ie, n, bb

      allocate (F(n_mom, nv))
      allocate (mu_sp(nb, nb, n_mom, natoms), g0_atom(nb, nb, nv))
      call prepare_moments_sp(mu, nb, n_mom, natoms, mu_sp)
      cot = cos(pi/(n_mom + 1.0_rp))/sin(pi/(n_mom + 1.0_rp))
      !$omp parallel do private(ie, i, th, thl, gj, cc, pref, ang) schedule(static)
      do ie = 1, nv
         th = acos((ene(ie) - b)/a)
         pref = 1.0_rp/sqrt(a*a - (ene(ie) - b)**2)
         do i = 1, n_mom                       ! i-1 = Chebyshev order
            thl = pi*real(i - 1, rp)/(n_mom + 1.0_rp)
            gj = ((n_mom - (i - 1) + 1)*cos(thl) + sin(thl)*cot) &
                 /(n_mom + 1.0_rp)             ! Jackson, bands.f90 convention
            cc = merge(1.0_rp, 2.0_rp, i == 1)
            ang = real(i - 1, rp)*th           ! (-i) e^{-i ang}
            F(i, ie) = cmplx(gj*cc*pref, 0.0_rp, sp)*cmplx(-sin(ang), -cos(ang), sp)
         end do
      end do
      !$omp end parallel do

      bb = nb*nb
      do n = 1, natoms
         call cgemm('N', 'N', bb, nv, n_mom, cone, mu_sp(:, :, :, n), bb, F, n_mom, czero, g0_atom, bb)
         g0(:, :, :, n) = g0_atom
      end do
      deallocate (g0_atom, mu_sp, F)
   end subroutine cheb_green_fast

   subroutine prepare_moments_sp(mu_in, nb, n_mom, natoms, mu_out)
      integer, intent(in) :: nb, n_mom, natoms
      complex(rp), intent(in) :: mu_in(nb, nb, n_mom, natoms)
      complex(sp), intent(out) :: mu_out(nb, nb, n_mom, natoms)

      mu_out = cmplx(real(mu_in, sp), aimag(mu_in), sp)
   end subroutine prepare_moments_sp

end module chebyshev_fast_mod
