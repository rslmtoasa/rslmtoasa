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
   use iso_c_binding, only: c_ptr
   implicit none
   private
   public :: cheb_moments_fast, cheb_moments_fast_batched
   public :: cheb_moments_fast_mkl_batch, cheb_moments_fast_mkl_sparse
   public :: cheb_green_fast, cheb_fast_reset_cache
   integer, parameter :: sp = selected_real_kind(6, 37)
   integer, parameter :: rp = selected_real_kind(15, 307)
   real(rp), parameter :: pi = 3.14159265358979323846_rp
   complex(sp), allocatable, target, save :: hee_cache(:, :, :, :), hha_cache(:, :, :, :)
   complex(sp), allocatable, target, save :: oee_cache(:, :, :, :), oha_cache(:, :, :, :)
   complex(sp), allocatable, target, save :: bee_cache(:, :, :, :), bha_cache(:, :, :, :)
   complex(sp), allocatable, target, save :: hons_cache(:, :, :)
   complex(sp), allocatable, target, save :: hval_cache(:, :, :)
   complex(sp), allocatable, target, save :: work0_cache(:, :), work1_cache(:, :), work2_cache(:, :)
   complex(sp), allocatable, target, save :: workt_cache(:, :)
   complex(sp), allocatable, target, save :: block_products_cache(:, :, :)
   integer, allocatable, target, save :: hcol_cache(:), hrow_cache(:)
   type(c_ptr), allocatable, target, save :: mkl_a_ptr_cache(:), mkl_b_ptr_cache(:), mkl_c_ptr_cache(:)
   logical, save :: scaled_cache_valid = .false.
   logical, save :: bsr_cache_valid = .false.
   logical, save :: ortho_cache_valid = .false.
   integer, save :: scaled_cache_nb = 0, scaled_cache_nnmax = 0, scaled_cache_ntype = 0
   integer, save :: scaled_cache_nmax = -1, scaled_cache_kk = 0
   integer, save :: ortho_cache_nb = 0, ortho_cache_nnmax = 0, ortho_cache_ntype = 0
   integer, save :: ortho_cache_nmax = -1, ortho_cache_kk = 0
   real(sp), save :: ortho_cache_inva = -1.0_sp, ortho_cache_b = huge(1.0_sp)
   integer, save :: bsr_cache_nb = 0, bsr_cache_nnmax = 0, bsr_cache_ntype = 0
   integer, save :: bsr_cache_nmax = -1, bsr_cache_kk = 0
   real(sp), save :: scaled_cache_inva = -1.0_sp, scaled_cache_b = huge(1.0_sp)
   real(sp), save :: bsr_cache_inva = -1.0_sp, bsr_cache_b = huge(1.0_sp)
   integer, save :: work_cache_ld = 0, work_cache_nb = 0
   integer, save :: block_products_cache_nb = 0, block_products_cache_nblocks = 0
   integer, save :: mkl_ptr_cache_nblocks = 0

contains

   subroutine cheb_fast_reset_cache()
      if (allocated(hee_cache)) deallocate (hee_cache)
      if (allocated(hha_cache)) deallocate (hha_cache)
      if (allocated(oee_cache)) deallocate (oee_cache)
      if (allocated(oha_cache)) deallocate (oha_cache)
      if (allocated(bee_cache)) deallocate (bee_cache)
      if (allocated(bha_cache)) deallocate (bha_cache)
      if (allocated(hons_cache)) deallocate (hons_cache)
      if (allocated(hval_cache)) deallocate (hval_cache)
      if (allocated(work0_cache)) deallocate (work0_cache)
      if (allocated(work1_cache)) deallocate (work1_cache)
      if (allocated(work2_cache)) deallocate (work2_cache)
      if (allocated(workt_cache)) deallocate (workt_cache)
      if (allocated(block_products_cache)) deallocate (block_products_cache)
      if (allocated(hcol_cache)) deallocate (hcol_cache)
      if (allocated(hrow_cache)) deallocate (hrow_cache)
      if (allocated(mkl_a_ptr_cache)) deallocate (mkl_a_ptr_cache)
      if (allocated(mkl_b_ptr_cache)) deallocate (mkl_b_ptr_cache)
      if (allocated(mkl_c_ptr_cache)) deallocate (mkl_c_ptr_cache)
      scaled_cache_valid = .false.
      bsr_cache_valid = .false.
      ortho_cache_valid = .false.
      work_cache_ld = 0
      work_cache_nb = 0
      block_products_cache_nb = 0
      block_products_cache_nblocks = 0
      mkl_ptr_cache_nblocks = 0
   end subroutine cheb_fast_reset_cache

   subroutine ensure_work_buffers(ld_in, nb_in, w0, w1, w2)
      integer, intent(in) :: ld_in, nb_in
      complex(sp), pointer, intent(out) :: w0(:, :), w1(:, :), w2(:, :)

      if (.not. allocated(work0_cache) .or. work_cache_ld /= ld_in .or. work_cache_nb /= nb_in) then
         if (allocated(work0_cache)) deallocate (work0_cache)
         if (allocated(work1_cache)) deallocate (work1_cache)
         if (allocated(work2_cache)) deallocate (work2_cache)
         allocate (work0_cache(ld_in, nb_in), work1_cache(ld_in, nb_in), work2_cache(ld_in, nb_in))
         work_cache_ld = ld_in
         work_cache_nb = nb_in
      end if

      w0 => work0_cache
      w1 => work1_cache
      w2 => work2_cache
   end subroutine ensure_work_buffers

   !> Extra (ld, nb) scratch holding the sweep-A result t = (h/a)*x for the
   !> two-sweep hoh apply. Shares the work-buffer (ld, nb) shape; allocated only
   !> when hoh is active.
   subroutine ensure_hoh_buffer(ld_in, nb_in, wt)
      integer, intent(in) :: ld_in, nb_in
      complex(sp), pointer, contiguous, intent(out) :: wt(:, :)

      if (.not. allocated(workt_cache) .or. size(workt_cache, 1) /= ld_in .or. size(workt_cache, 2) /= nb_in) then
         if (allocated(workt_cache)) deallocate (workt_cache)
         allocate (workt_cache(ld_in, nb_in))
      end if
      wt => workt_cache
   end subroutine ensure_hoh_buffer

   subroutine ensure_block_products(nb_in, nblocks_in, block_products)
      integer, intent(in) :: nb_in, nblocks_in
      complex(sp), pointer, intent(out) :: block_products(:, :, :)

      if (.not. allocated(block_products_cache) &
          .or. block_products_cache_nb /= nb_in &
          .or. block_products_cache_nblocks /= nblocks_in) then
         if (allocated(block_products_cache)) deallocate (block_products_cache)
         allocate (block_products_cache(nb_in, nb_in, nblocks_in))
         block_products_cache_nb = nb_in
         block_products_cache_nblocks = nblocks_in
      end if

      block_products => block_products_cache
   end subroutine ensure_block_products

   subroutine ensure_mkl_batch_ptr_arrays(nblocks_in, a_ptr, b_ptr, c_ptrs)
      integer, intent(in) :: nblocks_in
      type(c_ptr), pointer, intent(out) :: a_ptr(:), b_ptr(:), c_ptrs(:)

      if (.not. allocated(mkl_a_ptr_cache) .or. mkl_ptr_cache_nblocks /= nblocks_in) then
         if (allocated(mkl_a_ptr_cache)) deallocate (mkl_a_ptr_cache)
         if (allocated(mkl_b_ptr_cache)) deallocate (mkl_b_ptr_cache)
         if (allocated(mkl_c_ptr_cache)) deallocate (mkl_c_ptr_cache)
         allocate (mkl_a_ptr_cache(nblocks_in), mkl_b_ptr_cache(nblocks_in), mkl_c_ptr_cache(nblocks_in))
         mkl_ptr_cache_nblocks = nblocks_in
      end if

      a_ptr => mkl_a_ptr_cache
      b_ptr => mkl_b_ptr_cache
      c_ptrs => mkl_c_ptr_cache
   end subroutine ensure_mkl_batch_ptr_arrays

   subroutine ensure_scaled_hamiltonian_sp(ee_in, hall_in, lsham_in, iz_in, kk_in, nb_in, nnmax_in, ntype_in, nmax_in, inva_in, b_in, hee_out, hha_out)
      complex(rp), intent(in) :: ee_in(nb_in, nb_in, nnmax_in, ntype_in)
      complex(rp), intent(in) :: hall_in(nb_in, nb_in, nnmax_in, *)
      complex(rp), intent(in) :: lsham_in(nb_in, nb_in, ntype_in)
      integer, intent(in) :: iz_in(kk_in), kk_in, nb_in, nnmax_in, ntype_in, nmax_in
      real(sp), intent(in) :: inva_in, b_in
      complex(sp), pointer, intent(out) :: hee_out(:, :, :, :), hha_out(:, :, :, :)
      integer :: t_in, k_in, l_in
      logical :: valid

      valid = scaled_cache_valid &
         .and. scaled_cache_nb == nb_in &
         .and. scaled_cache_nnmax == nnmax_in &
         .and. scaled_cache_ntype == ntype_in &
         .and. scaled_cache_nmax == nmax_in &
         .and. scaled_cache_kk == kk_in &
         .and. scaled_cache_inva == inva_in &
         .and. scaled_cache_b == b_in

      if (.not. valid) then
         if (allocated(hee_cache)) deallocate (hee_cache)
         if (allocated(hha_cache)) deallocate (hha_cache)

         allocate (hee_cache(nb_in, nb_in, nnmax_in, ntype_in))
         hee_cache = cmplx(real(ee_in, sp), aimag(ee_in), sp)*inva_in
         do t_in = 1, ntype_in
            hee_cache(:, :, 1, t_in) = hee_cache(:, :, 1, t_in) + &
               cmplx(real(lsham_in(:, :, t_in), sp), aimag(lsham_in(:, :, t_in)), sp)*inva_in
            do l_in = 1, nb_in
               hee_cache(l_in, l_in, 1, t_in) = hee_cache(l_in, l_in, 1, t_in) - b_in*inva_in
            end do
         end do

         if (nmax_in > 0) then
            allocate (hha_cache(nb_in, nb_in, nnmax_in, nmax_in))
            hha_cache = cmplx(real(hall_in(:, :, :, 1:nmax_in), sp), aimag(hall_in(:, :, :, 1:nmax_in)), sp)*inva_in
            do k_in = 1, nmax_in
               t_in = iz_in(k_in)
               hha_cache(:, :, 1, k_in) = hha_cache(:, :, 1, k_in) + &
                  cmplx(real(lsham_in(:, :, t_in), sp), aimag(lsham_in(:, :, t_in)), sp)*inva_in
               do l_in = 1, nb_in
                  hha_cache(l_in, l_in, 1, k_in) = hha_cache(l_in, l_in, 1, k_in) - b_in*inva_in
               end do
            end do
         end if

         scaled_cache_nb = nb_in
         scaled_cache_nnmax = nnmax_in
         scaled_cache_ntype = ntype_in
         scaled_cache_nmax = nmax_in
         scaled_cache_kk = kk_in
         scaled_cache_inva = inva_in
         scaled_cache_b = b_in
         scaled_cache_valid = .true.
      end if

      hee_out => hee_cache
      if (allocated(hha_cache)) then
         hha_out => hha_cache
      else
         nullify (hha_out)
      end if
   end subroutine ensure_scaled_hamiltonian_sp

   !> Build the fp32 operands for the two-sweep hoh apply, mirroring
   !> ham_hoh_vec_matmul (recursion.f90). The combine there is
   !>   y = ( h*x - eeo*(h*x) + (enim + lsham)*x - b*x ) / a
   !> so the eeo sweep must see the BARE h*x (no on-site enim/lsham/-bI fold).
   !> We therefore cache four fp32 operands, all consistent with one 1/a scaling:
   !>   bee/bha = (ee/hall)*inva                 -> sweep A: t = inva*h*x
   !>   oee/oha = eeo/hallo  (RAW, unscaled)     -> sweep B: eeo*t = inva*eeo*h*x
   !>   hons    = (enim + lsham - b*I)*inva      -> on-site, applied to x
   !> Cached independently of the standard scaled-Hamiltonian cache.
   subroutine ensure_scaled_ortho_sp(ee_in, hall_in, eeo_in, hallo_in, lsham_in, enim_in, &
                                     iz_in, kk_in, nb_in, nnmax_in, ntype_in, nmax_in, inva_in, b_in, &
                                     bee_out, bha_out, oee_out, oha_out, hons_out)
      complex(rp), intent(in) :: ee_in(nb_in, nb_in, nnmax_in, ntype_in)
      complex(rp), intent(in) :: hall_in(nb_in, nb_in, nnmax_in, *)
      complex(rp), intent(in) :: eeo_in(nb_in, nb_in, nnmax_in, ntype_in)
      complex(rp), intent(in) :: hallo_in(nb_in, nb_in, nnmax_in, *)
      complex(rp), intent(in) :: lsham_in(nb_in, nb_in, ntype_in)
      complex(rp), intent(in) :: enim_in(nb_in, nb_in, ntype_in)
      integer, intent(in) :: iz_in(kk_in), kk_in, nb_in, nnmax_in, ntype_in, nmax_in
      real(sp), intent(in) :: inva_in, b_in
      complex(sp), pointer, intent(out) :: bee_out(:, :, :, :), bha_out(:, :, :, :)
      complex(sp), pointer, intent(out) :: oee_out(:, :, :, :), oha_out(:, :, :, :)
      complex(sp), pointer, intent(out) :: hons_out(:, :, :)
      integer :: t_in, k_in, l_in
      logical :: valid

      valid = ortho_cache_valid &
         .and. ortho_cache_nb == nb_in &
         .and. ortho_cache_nnmax == nnmax_in &
         .and. ortho_cache_ntype == ntype_in &
         .and. ortho_cache_nmax == nmax_in &
         .and. ortho_cache_kk == kk_in &
         .and. ortho_cache_inva == inva_in &
         .and. ortho_cache_b == b_in

      if (.not. valid) then
         if (allocated(bee_cache)) deallocate (bee_cache)
         if (allocated(bha_cache)) deallocate (bha_cache)
         if (allocated(oee_cache)) deallocate (oee_cache)
         if (allocated(oha_cache)) deallocate (oha_cache)
         if (allocated(hons_cache)) deallocate (hons_cache)

         ! Sweep-A bare h, scaled by inva (no on-site fold).
         allocate (bee_cache(nb_in, nb_in, nnmax_in, ntype_in))
         bee_cache = cmplx(real(ee_in, sp), aimag(ee_in), sp)*inva_in

         ! Sweep-B eeo, raw.
         allocate (oee_cache(nb_in, nb_in, nnmax_in, ntype_in))
         oee_cache = cmplx(real(eeo_in, sp), aimag(eeo_in), sp)

         ! On-site (enim + lsham - b*I)*inva, per type.
         allocate (hons_cache(nb_in, nb_in, ntype_in))
         hons_cache = (cmplx(real(enim_in, sp), aimag(enim_in), sp) + &
                       cmplx(real(lsham_in, sp), aimag(lsham_in), sp))*inva_in
         do t_in = 1, ntype_in
            do l_in = 1, nb_in
               hons_cache(l_in, l_in, t_in) = hons_cache(l_in, l_in, t_in) - b_in*inva_in
            end do
         end do

         if (nmax_in > 0) then
            allocate (bha_cache(nb_in, nb_in, nnmax_in, nmax_in))
            bha_cache = cmplx(real(hall_in(:, :, :, 1:nmax_in), sp), aimag(hall_in(:, :, :, 1:nmax_in)), sp)*inva_in
            allocate (oha_cache(nb_in, nb_in, nnmax_in, nmax_in))
            oha_cache = cmplx(real(hallo_in(:, :, :, 1:nmax_in), sp), aimag(hallo_in(:, :, :, 1:nmax_in)), sp)
            ! Note: the impurity on-site enim/lsham fold is carried by hons via
            ! iz(k), identical to the bulk path (ham_hoh_vec_matmul applies
            ! lsham/enim on-site for both local and bulk regions).
         end if

         ortho_cache_nb = nb_in
         ortho_cache_nnmax = nnmax_in
         ortho_cache_ntype = ntype_in
         ortho_cache_nmax = nmax_in
         ortho_cache_kk = kk_in
         ortho_cache_inva = inva_in
         ortho_cache_b = b_in
         ortho_cache_valid = .true.
      end if

      bee_out => bee_cache
      oee_out => oee_cache
      hons_out => hons_cache
      if (allocated(bha_cache)) then
         bha_out => bha_cache
         oha_out => oha_cache
      else
         nullify (bha_out)
         nullify (oha_out)
      end if
   end subroutine ensure_scaled_ortho_sp

   subroutine ensure_scaled_bsr_sp(ee_in, hall_in, lsham_in, nn_in, iz_in, kk_in, nb_in, nnmax_in, ntype_in, nmax_in, inva_in, b_in, hval_out, hcol_out, hrow_out)
      complex(rp), intent(in) :: ee_in(nb_in, nb_in, nnmax_in, ntype_in)
      complex(rp), intent(in) :: hall_in(nb_in, nb_in, nnmax_in, *)
      complex(rp), intent(in) :: lsham_in(nb_in, nb_in, ntype_in)
      integer, intent(in) :: nn_in(kk_in, nnmax_in), iz_in(kk_in)
      integer, intent(in) :: kk_in, nb_in, nnmax_in, ntype_in, nmax_in
      real(sp), intent(in) :: inva_in, b_in
      complex(sp), pointer, intent(out) :: hval_out(:, :, :)
      integer, pointer, intent(out) :: hcol_out(:), hrow_out(:)
      integer :: atom, neighbor_idx, neighbor, block_col, block_idx, nblocks, num_neighbors, ih, l_in
      logical :: valid

      valid = bsr_cache_valid &
         .and. bsr_cache_nb == nb_in &
         .and. bsr_cache_nnmax == nnmax_in &
         .and. bsr_cache_ntype == ntype_in &
         .and. bsr_cache_nmax == nmax_in &
         .and. bsr_cache_kk == kk_in &
         .and. bsr_cache_inva == inva_in &
         .and. bsr_cache_b == b_in

      if (.not. valid) then
         if (allocated(hval_cache)) deallocate (hval_cache)
         if (allocated(hcol_cache)) deallocate (hcol_cache)
         if (allocated(hrow_cache)) deallocate (hrow_cache)

         nblocks = 0
         do atom = 1, kk_in
            num_neighbors = nn_in(atom, 1)
            do neighbor_idx = 1, num_neighbors
               if (neighbor_idx == 1 .or. nn_in(atom, neighbor_idx) /= 0) nblocks = nblocks + 1
            end do
         end do

         allocate (hval_cache(nb_in, nb_in, nblocks), hcol_cache(nblocks), hrow_cache(kk_in + 1))
         block_idx = 0
         hrow_cache(1) = 1

         do atom = 1, kk_in
            num_neighbors = nn_in(atom, 1)
            do neighbor_idx = 1, num_neighbors
               if (neighbor_idx == 1) then
                  block_col = atom
               else
                  neighbor = nn_in(atom, neighbor_idx)
                  if (neighbor == 0) cycle
                  block_col = neighbor
               end if

               block_idx = block_idx + 1
               if (nmax_in > 0 .and. atom <= nmax_in) then
                  hval_cache(:, :, block_idx) = cmplx(real(hall_in(:, :, neighbor_idx, atom), sp), &
                                                      aimag(hall_in(:, :, neighbor_idx, atom)), sp)*inva_in
                  ih = iz_in(atom)
               else
                  ih = iz_in(atom)
                  hval_cache(:, :, block_idx) = cmplx(real(ee_in(:, :, neighbor_idx, ih), sp), &
                                                      aimag(ee_in(:, :, neighbor_idx, ih)), sp)*inva_in
               end if

               if (neighbor_idx == 1) then
                  hval_cache(:, :, block_idx) = hval_cache(:, :, block_idx) + &
                     cmplx(real(lsham_in(:, :, ih), sp), aimag(lsham_in(:, :, ih)), sp)*inva_in
                  do l_in = 1, nb_in
                     hval_cache(l_in, l_in, block_idx) = hval_cache(l_in, l_in, block_idx) - b_in*inva_in
                  end do
               end if
               hcol_cache(block_idx) = block_col
            end do
            if (atom < kk_in) hrow_cache(atom + 1) = block_idx + 1
         end do
         hrow_cache(kk_in + 1) = nblocks + 1

         bsr_cache_nb = nb_in
         bsr_cache_nnmax = nnmax_in
         bsr_cache_ntype = ntype_in
         bsr_cache_nmax = nmax_in
         bsr_cache_kk = kk_in
         bsr_cache_inva = inva_in
         bsr_cache_b = b_in
         bsr_cache_valid = .true.
      end if

      hval_out => hval_cache
      hcol_out => hcol_cache
      hrow_out => hrow_cache
   end subroutine ensure_scaled_bsr_sp

   !> Block Chebyshev moments for one starting state.
   !> psi0  : (nb, nb, kk) starting block state (site or ij sign combos)
   !> ee    : (nb, nb, nnmax, ntype), hall: (nb, nb, nnmax, nmax) (or nmax=0)
   !> lsham : (nb, nb, ntype); nn: (kk, nnmax) with nn(k,1)=count; iz: (kk)
   !> mu    : (nb, nb, 2*lld+2) out, identical ordering to mu_n
   subroutine cheb_moments_fast(psi0, ee, hall, lsham, nn, iz, kk, nb, &
                                nnmax, ntype, nmax, lld, a, b, mu, &
                                hoh, eeo, hallo, enim)
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax, lld
      complex(rp), intent(in) :: psi0(nb, nb, kk), ee(nb, nb, nnmax, ntype)
      complex(rp), intent(in) :: hall(nb, nb, nnmax, *), lsham(nb, nb, ntype)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      real(rp), intent(in) :: a, b
      complex(rp), intent(out) :: mu(nb, nb, 2*lld + 2)
      logical, intent(in), optional :: hoh
      complex(rp), intent(in), optional :: eeo(nb, nb, nnmax, ntype)
      complex(rp), intent(in), optional :: hallo(nb, nb, nnmax, *)
      complex(rp), intent(in), optional :: enim(nb, nb, ntype)
      ! Locals
      complex(sp), pointer :: hee(:, :, :, :), hha(:, :, :, :)
      complex(sp), pointer :: bee(:, :, :, :), bha(:, :, :, :)
      complex(sp), pointer :: oee(:, :, :, :), oha(:, :, :, :)
      complex(sp), pointer :: hons(:, :, :)
      complex(sp), pointer :: w0(:, :), w1(:, :), w2(:, :)
      complex(sp), pointer, contiguous :: wt(:, :)
      complex(sp), pointer :: p0(:, :), p1(:, :), p2(:, :), ptmp(:, :)
      complex(sp) :: mu1_sp(nb, nb), mu2_sp(nb, nb), dum_sp(nb, nb)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp)
      integer :: k, l, c, ld, ll
      real(sp) :: inva, a_sp, b_sp
      logical :: do_hoh

      do_hoh = .false.
      if (present(hoh)) do_hoh = hoh

      ld = nb*kk
      a_sp = real(a, sp)
      b_sp = real(b, sp)
      inva = 1.0_sp/a_sp
      call ensure_work_buffers(ld, nb, w0, w1, w2)
      p0 => w0
      p1 => w1
      p2 => w2

      if (do_hoh) then
         ! --- two-sweep hoh operands (bare h, eeo, on-site) + extra buffer ---
         call ensure_hoh_buffer(ld, nb, wt)
         call ensure_scaled_ortho_sp(ee, hall, eeo, hallo, lsham, enim, iz, kk, nb, &
                                     nnmax, ntype, nmax, inva, b_sp, bee, bha, oee, oha, hons)
      else
         ! --- 1. scaled operator copies: Ht = (H + lsham - b*I)/a ----------
         call ensure_scaled_hamiltonian_sp(ee, hall, lsham, iz, kk, nb, nnmax, ntype, nmax, inva, b_sp, hee, hha)
      end if

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
      if (do_hoh) then
         call apply_step_hoh(p0, p0, p1, 1.0_sp, 0.0_sp, wt)
      else
         call apply_step(p0, p0, p1, 1.0_sp, 0.0_sp)
      end if
      call cgemm('C', 'N', nb, nb, ld, cone, p0, ld, p1, ld, czero, mu2_sp, nb)
      mu(:, :, 2) = mu2_sp

      do ll = 1, lld
         if (do_hoh) then
            call apply_step_hoh(p1, p0, p2, 2.0_sp, -1.0_sp, wt) ! p2 = 2 Ht p1 - p0
         else
            call apply_step(p1, p0, p2, 2.0_sp, -1.0_sp)   ! p2 = 2 Ht p1 - p0
         end if
         call cherk_full(nb, ld, p1, dum_sp)            ! dum1 = p1^H p1
         mu(:, :, 2*ll + 1) = 2.0_sp*dum_sp - mu1_sp
         call cgemm('C', 'N', nb, nb, ld, cone, p2, ld, p1, ld, czero, dum_sp, nb)
         mu(:, :, 2*ll + 2) = 2.0_sp*dum_sp - mu2_sp
         ptmp => p0
         p0 => p1
         p1 => p2
         p2 => ptmp
      end do

   contains

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

      !> Two-sweep hoh apply, mirroring ham_hoh_vec_matmul (recursion.f90):
      !>   sweep A:  t   = (h/a) * x1                  (bee/bha, into wt)
      !>   sweep B:  acc = t - eeo*t + hons*x1         (oee/oha raw, hons on-site)
      !>   combine:  y   = alpha*acc + beta*x0
      !> The 1/a scaling lives in bee/bha and hons; eeo is raw, so eeo*t already
      !> carries one factor of 1/a. The eeo sweep sees the BARE h*x (no on-site
      !> enim/lsham/-bI), exactly as in the legacy two-step recursion.
      subroutine apply_step_hoh(x1, x0, y, alpha, beta, t)
         complex(sp), intent(in) :: x1(ld, nb), x0(ld, nb)
         complex(sp), intent(out) :: y(ld, nb)
         real(sp), intent(in) :: alpha, beta
         complex(sp), intent(inout) :: t(ld, nb)   ! sweep-A scratch (= wt)
         complex(sp) :: acc(nb, nb)
         integer :: kk_, t_, s_, nbr, nr, r0

         ! --- sweep A: t = (h/a) * x1 ------------------------------------
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
                  call cgemm('N', 'N', nb, nb, nb, cone, bha(:, :, s_, kk_), nb, x1(r0 + 1, 1), ld, cone, acc, nb)
               else
                  t_ = iz(kk_)
                  call cgemm('N', 'N', nb, nb, nb, cone, bee(:, :, s_, t_), nb, x1(r0 + 1, 1), ld, cone, acc, nb)
               end if
            end do
            r0 = nb*(kk_ - 1)
            t(r0 + 1:r0 + nb, :) = acc
         end do
         !$omp end parallel do

         ! --- sweep B: acc = t - eeo*t + hons*x1 ; y = alpha*acc + beta*x0
         !$omp parallel do private(kk_, t_, s_, nbr, nr, r0, acc) schedule(dynamic, 32)
         do kk_ = 1, kk
            acc = (0.0_sp, 0.0_sp)
            nr = nn(kk_, 1)
            ! - eeo * t  (subtract the orthogonalisation contribution)
            do s_ = 1, nr
               if (s_ == 1) then
                  nbr = kk_
               else
                  nbr = nn(kk_, s_)
                  if (nbr == 0) cycle
               end if
               r0 = nb*(nbr - 1)
               if (kk_ <= nmax) then
                  call cgemm('N', 'N', nb, nb, nb, -cone, oha(:, :, s_, kk_), nb, t(r0 + 1, 1), ld, cone, acc, nb)
               else
                  t_ = iz(kk_)
                  call cgemm('N', 'N', nb, nb, nb, -cone, oee(:, :, s_, t_), nb, t(r0 + 1, 1), ld, cone, acc, nb)
               end if
            end do
            ! + hons * x1  (on-site enim + lsham - b*I, scaled by 1/a)
            t_ = iz(kk_)
            r0 = nb*(kk_ - 1)
            call cgemm('N', 'N', nb, nb, nb, cone, hons(:, :, t_), nb, x1(r0 + 1, 1), ld, cone, acc, nb)
            ! + t(kk_)  (the bare h*x/a contribution for this site)
            acc = acc + t(r0 + 1:r0 + nb, :)
            if (beta /= 0.0_sp) then
               y(r0 + 1:r0 + nb, :) = alpha*acc + beta*x0(r0 + 1:r0 + nb, :)
            else
               y(r0 + 1:r0 + nb, :) = alpha*acc
            end if
         end do
         !$omp end parallel do
      end subroutine apply_step_hoh

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

   !> Block Chebyshev moments using a scaled single-precision BSR operator.
   !> This is an optional tuning path for BLAS implementations that handle a
   !> smaller number of wider GEMMs better than many tiny per-neighbor GEMMs.
   !> It preserves the cheb_moments_fast interface and moment ordering.
   subroutine cheb_moments_fast_batched(psi0, ee, hall, lsham, nn, iz, kk, nb, &
                                        nnmax, ntype, nmax, lld, a, b, mu)
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax, lld
      complex(rp), intent(in) :: psi0(nb, nb, kk), ee(nb, nb, nnmax, ntype)
      complex(rp), intent(in) :: hall(nb, nb, nnmax, *), lsham(nb, nb, ntype)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      real(rp), intent(in) :: a, b
      complex(rp), intent(out) :: mu(nb, nb, 2*lld + 2)
      complex(sp), pointer :: hval(:, :, :)
      complex(sp), pointer :: w0(:, :), w1(:, :), w2(:, :)
      complex(sp), pointer :: p0(:, :), p1(:, :), p2(:, :), ptmp(:, :)
      complex(sp), pointer :: block_products(:, :, :)
      complex(sp) :: mu1_sp(nb, nb), mu2_sp(nb, nb), dum_sp(nb, nb)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp)
      integer, pointer :: hcol(:), hrow(:)
      integer :: k, l, c, ld, ll, nblocks
      real(sp) :: inva, a_sp, b_sp

      ld = nb*kk
      a_sp = real(a, sp)
      b_sp = real(b, sp)
      inva = 1.0_sp/a_sp
      call ensure_work_buffers(ld, nb, w0, w1, w2)
      p0 => w0
      p1 => w1
      p2 => w2

      call ensure_scaled_bsr_sp(ee, hall, lsham, nn, iz, kk, nb, nnmax, ntype, nmax, inva, b_sp, hval, hcol, hrow)
      nblocks = size(hcol)
      call ensure_block_products(nb, nblocks, block_products)

      !$omp parallel do private(k, c, l) schedule(static)
      do k = 1, kk
         do c = 1, nb
            do l = 1, nb
               p0(l + nb*(k - 1), c) = cmplx(real(psi0(l, c, k), sp), aimag(psi0(l, c, k)), sp)
            end do
         end do
      end do
      !$omp end parallel do

      call cherk_full_batched(nb, ld, p0, mu1_sp)
      mu(:, :, 1) = mu1_sp
      call apply_step_bsr(p0, p0, p1, 1.0_sp, 0.0_sp)
      call cgemm('C', 'N', nb, nb, ld, cone, p0, ld, p1, ld, czero, mu2_sp, nb)
      mu(:, :, 2) = mu2_sp

      do ll = 1, lld
         call apply_step_bsr(p1, p0, p2, 2.0_sp, -1.0_sp)
         call cherk_full_batched(nb, ld, p1, dum_sp)
         mu(:, :, 2*ll + 1) = 2.0_sp*dum_sp - mu1_sp
         call cgemm('C', 'N', nb, nb, ld, cone, p2, ld, p1, ld, czero, dum_sp, nb)
         mu(:, :, 2*ll + 2) = 2.0_sp*dum_sp - mu2_sp
         ptmp => p0
         p0 => p1
         p1 => p2
         p2 => ptmp
      end do


   contains

      !> y = alpha*(Ht x1) + beta*x0 using one packed GEMM per BSR row.
      subroutine apply_step_bsr(x1, x0, y, alpha, beta)
         complex(sp), intent(in) :: x1(ld, nb), x0(ld, nb)
         complex(sp), intent(out) :: y(ld, nb)
         real(sp), intent(in) :: alpha, beta
         complex(sp) :: acc(nb, nb)
         integer :: row, blk, col, r0

         !$omp parallel do private(blk, col, r0) schedule(static)
         do blk = 1, nblocks
            col = hcol(blk)
            r0 = nb*(col - 1)
            call cgemm('N', 'N', nb, nb, nb, cone, hval(:, :, blk), nb, &
                       x1(r0 + 1, 1), ld, czero, block_products(:, :, blk), nb)
         end do
         !$omp end parallel do

         !$omp parallel do private(row, blk, r0, acc) schedule(dynamic, 32)
         do row = 1, kk
            acc = (0.0_sp, 0.0_sp)
            do blk = hrow(row), hrow(row + 1) - 1
               acc = acc + block_products(:, :, blk)
            end do
            r0 = nb*(row - 1)
            if (beta /= 0.0_sp) then
               y(r0 + 1:r0 + nb, :) = alpha*acc + beta*x0(r0 + 1:r0 + nb, :)
            else
               y(r0 + 1:r0 + nb, :) = alpha*acc
            end if
         end do
         !$omp end parallel do
      end subroutine apply_step_bsr

      subroutine cherk_full_batched(n, kdim, Amat, C)
         integer, intent(in) :: n, kdim
         complex(sp), intent(in) :: Amat(kdim, n)
         complex(sp), intent(out) :: C(n, n)
         integer :: i, j
         C = (0.0_sp, 0.0_sp)
         call cherk('L', 'C', n, kdim, 1.0_sp, Amat, kdim, 0.0_sp, C, n)
         do j = 2, n
            do i = 1, j - 1
               C(i, j) = conjg(C(j, i))
            end do
         end do
      end subroutine cherk_full_batched

      subroutine swapm_batched(x, y)
         complex(sp), intent(inout) :: x(ld, nb), y(ld, nb)
         complex(sp) :: tmp
         integer :: i, j
         !$omp parallel do private(i, j, tmp) schedule(static)
         do j = 1, nb
            do i = 1, ld
               tmp = x(i, j)
               x(i, j) = y(i, j)
               y(i, j) = tmp
            end do
         end do
         !$omp end parallel do
      end subroutine swapm_batched

   end subroutine cheb_moments_fast_batched

   !> Optional Intel oneMKL cgemm_batch implementation.
   !> Compiled with -DUSE_MKL_BATCH and CMake option ENABLE_MKL_BATCH=ON.
   subroutine cheb_moments_fast_mkl_batch(psi0, ee, hall, lsham, nn, iz, kk, nb, &
                                          nnmax, ntype, nmax, lld, a, b, mu)
#ifdef USE_MKL_BATCH
      use iso_c_binding, only: c_int, c_float_complex, c_loc
#endif
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax, lld
      complex(rp), intent(in) :: psi0(nb, nb, kk), ee(nb, nb, nnmax, ntype)
      complex(rp), intent(in) :: hall(nb, nb, nnmax, *), lsham(nb, nb, ntype)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      real(rp), intent(in) :: a, b
      complex(rp), intent(out) :: mu(nb, nb, 2*lld + 2)

#ifndef USE_MKL_BATCH
      write (*, '(A)') 'ERROR: cheb_backend=mkl_batch requires CMake option ENABLE_MKL_BATCH=ON.'
      error stop
#else
      interface
         subroutine rslmto_mkl_cgemm_batch_nn(batch_count, m, n, k, alpha, a_array, lda, &
                                              b_array, ldb, beta, c_array, ldc, status) bind(C, name='rslmto_mkl_cgemm_batch_nn')
            import :: c_int, c_ptr, c_float_complex
            integer(c_int), intent(in) :: batch_count, m, n, k, lda, ldb, ldc
            complex(c_float_complex), intent(in) :: alpha, beta
            type(c_ptr), intent(in) :: a_array(*), b_array(*)
            type(c_ptr), intent(inout) :: c_array(*)
            integer(c_int), intent(out) :: status
         end subroutine rslmto_mkl_cgemm_batch_nn
      end interface

      complex(sp), pointer :: hval(:, :, :)
      complex(sp), pointer :: w0(:, :), w1(:, :), w2(:, :)
      complex(sp), pointer :: p0(:, :), p1(:, :), p2(:, :), ptmp(:, :)
      complex(sp), pointer :: block_products(:, :, :)
      complex(sp) :: mu1_sp(nb, nb), mu2_sp(nb, nb), dum_sp(nb, nb)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp)
      integer, pointer :: hcol(:), hrow(:)
      type(c_ptr), pointer :: a_ptr(:), b_ptr(:), c_ptrs(:)
      integer :: k, l, c, ld, ll, nblocks
      real(sp) :: inva, a_sp, b_sp

      ld = nb*kk
      a_sp = real(a, sp)
      b_sp = real(b, sp)
      inva = 1.0_sp/a_sp
      call ensure_work_buffers(ld, nb, w0, w1, w2)
      p0 => w0
      p1 => w1
      p2 => w2

      call ensure_scaled_bsr_sp(ee, hall, lsham, nn, iz, kk, nb, nnmax, ntype, nmax, inva, b_sp, hval, hcol, hrow)
      nblocks = size(hcol)
      call ensure_block_products(nb, nblocks, block_products)
      call ensure_mkl_batch_ptr_arrays(nblocks, a_ptr, b_ptr, c_ptrs)
      call init_mkl_batch_static_ptrs()

      !$omp parallel do private(k, c, l) schedule(static)
      do k = 1, kk
         do c = 1, nb
            do l = 1, nb
               p0(l + nb*(k - 1), c) = cmplx(real(psi0(l, c, k), sp), aimag(psi0(l, c, k)), sp)
            end do
         end do
      end do
      !$omp end parallel do

      call cherk_full_mkl_batch(nb, ld, p0, mu1_sp)
      mu(:, :, 1) = mu1_sp
      call apply_step_mkl_batch(p0, p0, p1, 1.0_sp, 0.0_sp)
      call cgemm('C', 'N', nb, nb, ld, cone, p0, ld, p1, ld, czero, mu2_sp, nb)
      mu(:, :, 2) = mu2_sp

      do ll = 1, lld
         call apply_step_mkl_batch(p1, p0, p2, 2.0_sp, -1.0_sp)
         call cherk_full_mkl_batch(nb, ld, p1, dum_sp)
         mu(:, :, 2*ll + 1) = 2.0_sp*dum_sp - mu1_sp
         call cgemm('C', 'N', nb, nb, ld, cone, p2, ld, p1, ld, czero, dum_sp, nb)
         mu(:, :, 2*ll + 2) = 2.0_sp*dum_sp - mu2_sp
         ptmp => p0
         p0 => p1
         p1 => p2
         p2 => ptmp
      end do


   contains

      subroutine init_mkl_batch_static_ptrs()
         integer :: block
         !$omp parallel do private(block) schedule(static)
         do block = 1, nblocks
            a_ptr(block) = c_loc(hval(1, 1, block))
            c_ptrs(block) = c_loc(block_products(1, 1, block))
         end do
         !$omp end parallel do
      end subroutine init_mkl_batch_static_ptrs

      subroutine apply_step_mkl_batch(x1, x0, y, alpha, beta)
         complex(sp), intent(in), target :: x1(ld, nb), x0(ld, nb)
         complex(sp), intent(out) :: y(ld, nb)
         real(sp), intent(in) :: alpha, beta
         complex(c_float_complex) :: alpha_c, beta_c
         complex(sp) :: acc(nb, nb)
         integer(c_int) :: batch_count_c, m_c, n_c, k_c, lda_c, ldb_c, ldc_c, status_c
         integer :: block, row, col, r0

         if (nblocks > huge(batch_count_c)) then
            write (*, '(A)') 'ERROR: mkl_batch backend has too many BSR blocks for the C wrapper.'
            error stop
         end if

         !$omp parallel do private(block, col, r0) schedule(static)
         do block = 1, nblocks
            col = hcol(block)
            r0 = nb*(col - 1)
            b_ptr(block) = c_loc(x1(r0 + 1, 1))
         end do
         !$omp end parallel do

         batch_count_c = int(nblocks, c_int)
         m_c = int(nb, c_int)
         n_c = int(nb, c_int)
         k_c = int(nb, c_int)
         lda_c = int(nb, c_int)
         ldb_c = int(ld, c_int)
         ldc_c = int(nb, c_int)
         alpha_c = cmplx(1.0_sp, 0.0_sp, c_float_complex)
         beta_c = cmplx(0.0_sp, 0.0_sp, c_float_complex)

         call rslmto_mkl_cgemm_batch_nn(batch_count_c, m_c, n_c, k_c, alpha_c, a_ptr, lda_c, &
                                        b_ptr, ldb_c, beta_c, c_ptrs, ldc_c, status_c)
         if (status_c /= 0_c_int) then
            write (*, '(A,I0)') 'ERROR: rslmto_mkl_cgemm_batch_nn failed with status ', status_c
            error stop
         end if

         !$omp parallel do private(row, block, r0, acc) schedule(dynamic, 32)
         do row = 1, kk
            acc = (0.0_sp, 0.0_sp)
            do block = hrow(row), hrow(row + 1) - 1
               acc = acc + block_products(:, :, block)
            end do
            r0 = nb*(row - 1)
            if (beta /= 0.0_sp) then
               y(r0 + 1:r0 + nb, :) = alpha*acc + beta*x0(r0 + 1:r0 + nb, :)
            else
               y(r0 + 1:r0 + nb, :) = alpha*acc
            end if
         end do
         !$omp end parallel do
      end subroutine apply_step_mkl_batch

      subroutine cherk_full_mkl_batch(n, kdim, Amat, C)
         integer, intent(in) :: n, kdim
         complex(sp), intent(in) :: Amat(kdim, n)
         complex(sp), intent(out) :: C(n, n)
         integer :: i, j
         C = (0.0_sp, 0.0_sp)
         call cherk('L', 'C', n, kdim, 1.0_sp, Amat, kdim, 0.0_sp, C, n)
         do j = 2, n
            do i = 1, j - 1
               C(i, j) = conjg(C(j, i))
            end do
         end do
      end subroutine cherk_full_mkl_batch

      subroutine swapm_mkl_batch(x, y)
         complex(sp), intent(inout) :: x(ld, nb), y(ld, nb)
         complex(sp) :: tmp
         integer :: i, j
         !$omp parallel do private(i, j, tmp) schedule(static)
         do j = 1, nb
            do i = 1, ld
               tmp = x(i, j)
               x(i, j) = y(i, j)
               y(i, j) = tmp
            end do
         end do
         !$omp end parallel do
      end subroutine swapm_mkl_batch
#endif
   end subroutine cheb_moments_fast_mkl_batch

   !> Optional Intel oneMKL Inspector-Executor Sparse BLAS implementation.
   !> Compiled with -DUSE_MKL_SPARSE.
   subroutine cheb_moments_fast_mkl_sparse(psi0, ee, hall, lsham, nn, iz, kk, nb, &
                                           nnmax, ntype, nmax, lld, a, b, mu)
#ifdef USE_MKL_SPARSE
      use mkl_spblas
#endif
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax, lld
      complex(rp), intent(in) :: psi0(nb, nb, kk), ee(nb, nb, nnmax, ntype)
      complex(rp), intent(in) :: hall(nb, nb, nnmax, *), lsham(nb, nb, ntype)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      real(rp), intent(in) :: a, b
      complex(rp), intent(out) :: mu(nb, nb, 2*lld + 2)

#ifndef USE_MKL_SPARSE
      write (*, '(A)') 'ERROR: cheb_backend=mkl_sparse requires CMake option ENABLE_MKL_SPARSE=ON.'
      error stop
#else
      type(sparse_matrix_t) :: mkl_A
      type(matrix_descr) :: descA
      complex(sp), pointer :: hval(:, :, :)
      complex(sp), pointer :: w0(:, :), w1(:, :), w2(:, :)
      complex(sp), pointer :: p0(:, :), p1(:, :), p2(:, :), ptmp(:, :)
      complex(sp) :: mu1_sp(nb, nb), mu2_sp(nb, nb), dum_sp(nb, nb)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp)
      integer, pointer :: hcol(:), hrow(:)
      integer :: status, k, l, c, ld, ll
      real(sp) :: inva, a_sp, b_sp

      ld = nb*kk
      a_sp = real(a, sp)
      b_sp = real(b, sp)
      inva = 1.0_sp/a_sp
      call ensure_work_buffers(ld, nb, w0, w1, w2)
      p0 => w0
      p1 => w1
      p2 => w2

      call ensure_scaled_bsr_sp(ee, hall, lsham, nn, iz, kk, nb, nnmax, ntype, nmax, inva, b_sp, hval, hcol, hrow)

      status = mkl_sparse_c_create_bsr(mkl_A, SPARSE_INDEX_BASE_ONE, SPARSE_LAYOUT_COLUMN_MAJOR, &
                                       kk, kk, nb, hrow(1:kk), hrow(2:kk + 1), hcol, hval)
      call check_mkl_sparse_status(status, 'mkl_sparse_c_create_bsr')

      descA%type = SPARSE_MATRIX_TYPE_GENERAL
      descA%mode = SPARSE_FILL_MODE_FULL
      descA%diag = SPARSE_DIAG_NON_UNIT

      status = mkl_sparse_set_mm_hint(mkl_A, SPARSE_OPERATION_NON_TRANSPOSE, descA, &
                                      SPARSE_LAYOUT_COLUMN_MAJOR, nb, max(2*lld + 1, 1))
      call check_mkl_sparse_status(status, 'mkl_sparse_set_mm_hint')

      status = mkl_sparse_optimize(mkl_A)
      call check_mkl_sparse_status(status, 'mkl_sparse_optimize')

      !$omp parallel do private(k, c, l) schedule(static)
      do k = 1, kk
         do c = 1, nb
            do l = 1, nb
               p0(l + nb*(k - 1), c) = cmplx(real(psi0(l, c, k), sp), aimag(psi0(l, c, k)), sp)
            end do
         end do
      end do
      !$omp end parallel do

      call cherk_full_mkl(nb, ld, p0, mu1_sp)
      mu(:, :, 1) = mu1_sp
      call apply_step_mkl(p0, p0, p1, 1.0_sp, 0.0_sp)
      call cgemm('C', 'N', nb, nb, ld, cone, p0, ld, p1, ld, czero, mu2_sp, nb)
      mu(:, :, 2) = mu2_sp

      do ll = 1, lld
         call apply_step_mkl(p1, p0, p2, 2.0_sp, -1.0_sp)
         call cherk_full_mkl(nb, ld, p1, dum_sp)
         mu(:, :, 2*ll + 1) = 2.0_sp*dum_sp - mu1_sp
         call cgemm('C', 'N', nb, nb, ld, cone, p2, ld, p1, ld, czero, dum_sp, nb)
         mu(:, :, 2*ll + 2) = 2.0_sp*dum_sp - mu2_sp
         ptmp => p0
         p0 => p1
         p1 => p2
         p2 => ptmp
      end do

      status = mkl_sparse_destroy(mkl_A)
      call check_mkl_sparse_status(status, 'mkl_sparse_destroy')

   contains

      subroutine apply_step_mkl(x1, x0, y, alpha, beta)
         complex(sp), intent(in) :: x1(ld, nb), x0(ld, nb)
         complex(sp), intent(out) :: y(ld, nb)
         real(sp), intent(in) :: alpha, beta
         complex(sp) :: alpha_c, beta_c
         integer :: status_mkl

         y = x0
         alpha_c = cmplx(alpha, 0.0_sp, sp)
         beta_c = cmplx(beta, 0.0_sp, sp)
         status_mkl = mkl_sparse_c_mm(SPARSE_OPERATION_NON_TRANSPOSE, alpha_c, mkl_A, descA, &
                                      SPARSE_LAYOUT_COLUMN_MAJOR, x1, nb, ld, beta_c, y, ld)
         call check_mkl_sparse_status(status_mkl, 'mkl_sparse_c_mm')
      end subroutine apply_step_mkl

      subroutine cherk_full_mkl(n, kdim, Amat, C)
         integer, intent(in) :: n, kdim
         complex(sp), intent(in) :: Amat(kdim, n)
         complex(sp), intent(out) :: C(n, n)
         integer :: i, j
         C = (0.0_sp, 0.0_sp)
         call cherk('L', 'C', n, kdim, 1.0_sp, Amat, kdim, 0.0_sp, C, n)
         do j = 2, n
            do i = 1, j - 1
               C(i, j) = conjg(C(j, i))
            end do
         end do
      end subroutine cherk_full_mkl

      subroutine swapm_mkl(x, y)
         complex(sp), intent(inout) :: x(ld, nb), y(ld, nb)
         complex(sp) :: tmp
         integer :: i, j
         !$omp parallel do private(i, j, tmp) schedule(static)
         do j = 1, nb
            do i = 1, ld
               tmp = x(i, j)
               x(i, j) = y(i, j)
               y(i, j) = tmp
            end do
         end do
         !$omp end parallel do
      end subroutine swapm_mkl

      subroutine check_mkl_sparse_status(status_in, routine)
         integer, intent(in) :: status_in
         character(len=*), intent(in) :: routine
         if (status_in /= SPARSE_STATUS_SUCCESS) then
            write (*, '(A,A,A,I0)') 'ERROR: ', trim(routine), ' failed with MKL sparse status ', status_in
            error stop
         end if
      end subroutine check_mkl_sparse_status
#endif
   end subroutine cheb_moments_fast_mkl_sparse

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
