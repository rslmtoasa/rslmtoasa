!------------------------------------------------------------------------------
! Optimized CPU replacements for chebyshev_recur (recursion_mod) and
! chebyshev_green (bands_mod), applying the validated KPM tricks:
!
!  1. The window scaling (H-b)/a is folded ONCE into local copies of
!     ee/hall (lsham included), so the inner loop is psi2 = 2*Ht*psi1 - psi0
!     with no scale/shift passes.
!  2. psi is packed SITE-MAJOR as an (nb*kk, nb) matrix, so every block
!     moment sum_k psi^H phi is one GEMM (’C’,’N’, nb, nb, nb*kk) -- and
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
   public :: cheb_moments_stochastic_fast, cheb_moments_orbital_fast
   integer, parameter :: sp = selected_real_kind(6, 37)
   integer, parameter :: rp = selected_real_kind(15, 307)
   real(rp), parameter :: pi = 3.14159265358979323846_rp

   !> @brief Cache-validity fingerprint for one cheb_cache_t buffer group.
   !> @details A fingerprint records the problem size and scaling parameters an
   !>          fp32 buffer group was built for. Each ensure_* routine compares
   !>          the current call’s parameters against its stored fingerprint. A
   !>          mismatch marks the group invalid and forces a rebuild. valid
   !>          starts false, so the first call always rebuilds.
   !> @note Not MPI-aware. The cache assumes a single rank owns and reuses it.
   type :: cheb_cache_fingerprint_t
      logical :: valid = .false.
      integer :: nb = 0
      integer :: nnmax = 0
      integer :: ntype = 0
      integer :: nmax = -1
      integer :: kk = 0
      integer :: ld = 0
      integer :: nblocks = 0
      real(sp) :: inva = -1.0_sp
      real(sp) :: b = huge(1.0_sp)
      logical :: hoh = .false.
   end type cheb_cache_fingerprint_t

   !> @brief Per-process cache of fp32 Chebyshev operator, work, and pointer
   !>        buffers.
   !> @details The cache holds fp32 mirrors of the Hamiltonian, ortho, BSR, and
   !>          velocity operator blocks (see the individual ensure_* routines
   !>          for what each group contains), plus shared work buffers and MKL
   !>          batch pointer arrays. Each group has its own
   !>          cheb_cache_fingerprint_t, so a change to only one group (for
   !>          example a new scaling window) does not invalidate the others.
   !>          cheb_moments_fast and its batched, MKL, hoh, stochastic, and
   !>          orbital variants read and populate this cache through the
   !>          ensure_* type-bound procedures; callers never touch the
   !>          pointers directly.
   !> @note One module-level instance (cheb_cache, below) is shared by every
   !>       recursion object in the process. It assumes a single MPI rank per
   !>       process. It is not safe to share across ranks, and it is not safe
   !>       to use from two recursion objects with different problem sizes at
   !>       the same time: each ensure_* call rebuilds the cache for the
   !>       CURRENT caller’s shape.
   type :: cheb_cache_t
      complex(sp), pointer :: hee_cache(:, :, :, :) => null(), hha_cache(:, :, :, :) => null()
      complex(sp), pointer :: oee_cache(:, :, :, :) => null(), oha_cache(:, :, :, :) => null()
      complex(sp), pointer :: bee_cache(:, :, :, :) => null(), bha_cache(:, :, :, :) => null()
      complex(sp), pointer :: hons_cache(:, :, :) => null()
      ! Raw fp32 mirrors for the conductivity/orbital fast kernels: velocity
      ! operators (fva/fvb), their hoh ortho (fvoa/fvob, on-site shell zeroed) and
      ! raw ee/hall (fvee/fvha) used by the velocity hoh two-sweep.
      complex(sp), pointer :: fva_cache(:, :, :, :) => null(), fvb_cache(:, :, :, :) => null()
      complex(sp), pointer :: fvoa_cache(:, :, :, :) => null(), fvob_cache(:, :, :, :) => null()
      complex(sp), pointer :: fvee_cache(:, :, :, :) => null(), fvha_cache(:, :, :, :) => null()
      complex(sp), pointer :: hval_cache(:, :, :) => null()
      complex(sp), pointer :: work0_cache(:, :) => null(), work1_cache(:, :) => null(), work2_cache(:, :) => null()
      complex(sp), pointer :: workt_cache(:, :) => null()
      complex(sp), pointer :: block_products_cache(:, :, :) => null()
      integer, pointer :: hcol_cache(:) => null(), hrow_cache(:) => null()
      type(c_ptr), pointer :: mkl_a_ptr_cache(:) => null(), mkl_b_ptr_cache(:) => null(), mkl_c_ptr_cache(:) => null()
      type(cheb_cache_fingerprint_t) :: scaled
      type(cheb_cache_fingerprint_t) :: bsr
      type(cheb_cache_fingerprint_t) :: ortho
      type(cheb_cache_fingerprint_t) :: vel
      type(cheb_cache_fingerprint_t) :: work
      type(cheb_cache_fingerprint_t) :: block_products
      type(cheb_cache_fingerprint_t) :: mkl_ptr
   contains
      procedure :: reset => cheb_cache_reset
      procedure :: ensure_work_buffers => cheb_cache_ensure_work_buffers
      procedure :: ensure_hoh_buffer => cheb_cache_ensure_hoh_buffer
      procedure :: ensure_block_products => cheb_cache_ensure_block_products
      procedure :: ensure_mkl_batch_ptr_arrays => cheb_cache_ensure_mkl_batch_ptr_arrays
      procedure :: ensure_scaled_hamiltonian_sp => cheb_cache_ensure_scaled_hamiltonian_sp
      procedure :: ensure_scaled_ortho_sp => cheb_cache_ensure_scaled_ortho_sp
      procedure :: ensure_scaled_bsr_sp => cheb_cache_ensure_scaled_bsr_sp
      procedure :: ensure_scaled_velocity_sp => cheb_cache_ensure_scaled_velocity_sp
   end type cheb_cache_t

   type(cheb_cache_t), save :: cheb_cache

contains

   subroutine cheb_fast_reset_cache()
      call cheb_cache%reset()
   end subroutine cheb_fast_reset_cache

   subroutine cheb_cache_reset(this)
      class(cheb_cache_t), intent(inout) :: this

      if (associated(this%hee_cache)) deallocate (this%hee_cache)
      if (associated(this%hha_cache)) deallocate (this%hha_cache)
      if (associated(this%oee_cache)) deallocate (this%oee_cache)
      if (associated(this%oha_cache)) deallocate (this%oha_cache)
      if (associated(this%bee_cache)) deallocate (this%bee_cache)
      if (associated(this%bha_cache)) deallocate (this%bha_cache)
      if (associated(this%hons_cache)) deallocate (this%hons_cache)
      if (associated(this%fva_cache)) deallocate (this%fva_cache)
      if (associated(this%fvb_cache)) deallocate (this%fvb_cache)
      if (associated(this%fvoa_cache)) deallocate (this%fvoa_cache)
      if (associated(this%fvob_cache)) deallocate (this%fvob_cache)
      if (associated(this%fvee_cache)) deallocate (this%fvee_cache)
      if (associated(this%fvha_cache)) deallocate (this%fvha_cache)
      if (associated(this%hval_cache)) deallocate (this%hval_cache)
      if (associated(this%work0_cache)) deallocate (this%work0_cache)
      if (associated(this%work1_cache)) deallocate (this%work1_cache)
      if (associated(this%work2_cache)) deallocate (this%work2_cache)
      if (associated(this%workt_cache)) deallocate (this%workt_cache)
      if (associated(this%block_products_cache)) deallocate (this%block_products_cache)
      if (associated(this%hcol_cache)) deallocate (this%hcol_cache)
      if (associated(this%hrow_cache)) deallocate (this%hrow_cache)
      if (associated(this%mkl_a_ptr_cache)) deallocate (this%mkl_a_ptr_cache)
      if (associated(this%mkl_b_ptr_cache)) deallocate (this%mkl_b_ptr_cache)
      if (associated(this%mkl_c_ptr_cache)) deallocate (this%mkl_c_ptr_cache)
      this%scaled = cheb_cache_fingerprint_t()
      this%bsr = cheb_cache_fingerprint_t()
      this%ortho = cheb_cache_fingerprint_t()
      this%vel = cheb_cache_fingerprint_t()
      this%work = cheb_cache_fingerprint_t()
      this%block_products = cheb_cache_fingerprint_t()
      this%mkl_ptr = cheb_cache_fingerprint_t()
   end subroutine cheb_cache_reset

   subroutine cheb_cache_ensure_work_buffers(this, ld_in, nb_in, w0, w1, w2)
      class(cheb_cache_t), intent(inout) :: this
      integer, intent(in) :: ld_in, nb_in
      complex(sp), pointer, intent(out) :: w0(:, :), w1(:, :), w2(:, :)

      if (.not. associated(this%work0_cache) .or. this%work%ld /= ld_in .or. this%work%nb /= nb_in) then
         if (associated(this%work0_cache)) deallocate (this%work0_cache)
         if (associated(this%work1_cache)) deallocate (this%work1_cache)
         if (associated(this%work2_cache)) deallocate (this%work2_cache)
         allocate (this%work0_cache(ld_in, nb_in), this%work1_cache(ld_in, nb_in), this%work2_cache(ld_in, nb_in))
         this%work%ld = ld_in
         this%work%nb = nb_in
      end if

      w0 => this%work0_cache
      w1 => this%work1_cache
      w2 => this%work2_cache
   end subroutine cheb_cache_ensure_work_buffers

   !> Extra (ld, nb) scratch holding the sweep-A result t = (h/a)*x for the
   !> two-sweep hoh apply. Shares the work-buffer (ld, nb) shape; allocated only
   !> when hoh is active.
   subroutine cheb_cache_ensure_hoh_buffer(this, ld_in, nb_in, wt)
      class(cheb_cache_t), intent(inout) :: this
      integer, intent(in) :: ld_in, nb_in
      complex(sp), pointer, contiguous, intent(out) :: wt(:, :)
      logical :: need

      ! NOTE: Fortran does not short-circuit .or., so size() must not be
      ! evaluated on an unallocated array. Test allocation status separately.
      need = .true.
      if (associated(this%workt_cache)) then
         need = (size(this%workt_cache, 1) /= ld_in .or. size(this%workt_cache, 2) /= nb_in)
      end if
      if (need) then
         if (associated(this%workt_cache)) deallocate (this%workt_cache)
         allocate (this%workt_cache(ld_in, nb_in))
      end if
      wt => this%workt_cache
   end subroutine cheb_cache_ensure_hoh_buffer

   subroutine cheb_cache_ensure_block_products(this, nb_in, nblocks_in, block_products)
      class(cheb_cache_t), intent(inout) :: this
      integer, intent(in) :: nb_in, nblocks_in
      complex(sp), pointer, intent(out) :: block_products(:, :, :)

      if (.not. associated(this%block_products_cache) &
          .or. this%block_products%nb /= nb_in &
          .or. this%block_products%nblocks /= nblocks_in) then
         if (associated(this%block_products_cache)) deallocate (this%block_products_cache)
         allocate (this%block_products_cache(nb_in, nb_in, nblocks_in))
         this%block_products%nb = nb_in
         this%block_products%nblocks = nblocks_in
      end if

      block_products => this%block_products_cache
   end subroutine cheb_cache_ensure_block_products

   subroutine cheb_cache_ensure_mkl_batch_ptr_arrays(this, nblocks_in, a_ptr, b_ptr, c_ptrs)
      class(cheb_cache_t), intent(inout) :: this
      integer, intent(in) :: nblocks_in
      type(c_ptr), pointer, intent(out) :: a_ptr(:), b_ptr(:), c_ptrs(:)

      if (.not. associated(this%mkl_a_ptr_cache) .or. this%mkl_ptr%nblocks /= nblocks_in) then
         if (associated(this%mkl_a_ptr_cache)) deallocate (this%mkl_a_ptr_cache)
         if (associated(this%mkl_b_ptr_cache)) deallocate (this%mkl_b_ptr_cache)
         if (associated(this%mkl_c_ptr_cache)) deallocate (this%mkl_c_ptr_cache)
         allocate (this%mkl_a_ptr_cache(nblocks_in), this%mkl_b_ptr_cache(nblocks_in), this%mkl_c_ptr_cache(nblocks_in))
         this%mkl_ptr%nblocks = nblocks_in
      end if

      a_ptr => this%mkl_a_ptr_cache
      b_ptr => this%mkl_b_ptr_cache
      c_ptrs => this%mkl_c_ptr_cache
   end subroutine cheb_cache_ensure_mkl_batch_ptr_arrays

   subroutine cheb_cache_ensure_scaled_hamiltonian_sp(this, ee_in, hall_in, lsham_in, iz_in, kk_in, nb_in, nnmax_in, ntype_in, nmax_in, inva_in, b_in, hee_out, hha_out)
      class(cheb_cache_t), intent(inout) :: this
      complex(rp), intent(in) :: ee_in(nb_in, nb_in, nnmax_in, ntype_in)
      complex(rp), intent(in) :: hall_in(nb_in, nb_in, nnmax_in, *)
      complex(rp), intent(in) :: lsham_in(nb_in, nb_in, ntype_in)
      integer, intent(in) :: iz_in(kk_in), kk_in, nb_in, nnmax_in, ntype_in, nmax_in
      real(sp), intent(in) :: inva_in, b_in
      complex(sp), pointer, intent(out) :: hee_out(:, :, :, :), hha_out(:, :, :, :)
      integer :: t_in, k_in, l_in
      logical :: valid

      valid = this%scaled%valid &
         .and. this%scaled%nb == nb_in &
         .and. this%scaled%nnmax == nnmax_in &
         .and. this%scaled%ntype == ntype_in &
         .and. this%scaled%nmax == nmax_in &
         .and. this%scaled%kk == kk_in &
         .and. this%scaled%inva == inva_in &
         .and. this%scaled%b == b_in

      if (.not. valid) then
         if (associated(this%hee_cache)) deallocate (this%hee_cache)
         if (associated(this%hha_cache)) deallocate (this%hha_cache)

         allocate (this%hee_cache(nb_in, nb_in, nnmax_in, ntype_in))
         this%hee_cache = cmplx(real(ee_in, sp), aimag(ee_in), sp)*inva_in
         do t_in = 1, ntype_in
            this%hee_cache(:, :, 1, t_in) = this%hee_cache(:, :, 1, t_in) + &
               cmplx(real(lsham_in(:, :, t_in), sp), aimag(lsham_in(:, :, t_in)), sp)*inva_in
            do l_in = 1, nb_in
               this%hee_cache(l_in, l_in, 1, t_in) = this%hee_cache(l_in, l_in, 1, t_in) - b_in*inva_in
            end do
         end do

         if (nmax_in > 0) then
            allocate (this%hha_cache(nb_in, nb_in, nnmax_in, nmax_in))
            this%hha_cache = cmplx(real(hall_in(:, :, :, 1:nmax_in), sp), aimag(hall_in(:, :, :, 1:nmax_in)), sp)*inva_in
            do k_in = 1, nmax_in
               t_in = iz_in(k_in)
               this%hha_cache(:, :, 1, k_in) = this%hha_cache(:, :, 1, k_in) + &
                  cmplx(real(lsham_in(:, :, t_in), sp), aimag(lsham_in(:, :, t_in)), sp)*inva_in
               do l_in = 1, nb_in
                  this%hha_cache(l_in, l_in, 1, k_in) = this%hha_cache(l_in, l_in, 1, k_in) - b_in*inva_in
               end do
            end do
         end if

         this%scaled%nb = nb_in
         this%scaled%nnmax = nnmax_in
         this%scaled%ntype = ntype_in
         this%scaled%nmax = nmax_in
         this%scaled%kk = kk_in
         this%scaled%inva = inva_in
         this%scaled%b = b_in
         this%scaled%valid = .true.
      end if

      hee_out => this%hee_cache
      if (associated(this%hha_cache)) then
         hha_out => this%hha_cache
      else
         nullify (hha_out)
      end if
   end subroutine cheb_cache_ensure_scaled_hamiltonian_sp

   !> Build the fp32 operands for the two-sweep hoh apply, mirroring
   !> ham_hoh_vec_matmul (recursion.f90). The combine there is
   !>   y = ( h*x - eeo*(h*x) + (enim + lsham)*x - b*x ) / a
   !> so the eeo sweep must see the BARE h*x (no on-site enim/lsham/-bI fold).
   !> We therefore cache four fp32 operands, all consistent with one 1/a scaling:
   !>   bee/bha = (ee/hall)*inva                 -> sweep A: t = inva*h*x
   !>   oee/oha = eeo/hallo  (RAW, unscaled)     -> sweep B: eeo*t = inva*eeo*h*x
   !>   hons    = (enim + lsham - b*I)*inva      -> on-site, applied to x
   !> Cached independently of the standard scaled-Hamiltonian cache.
   subroutine cheb_cache_ensure_scaled_ortho_sp(this, ee_in, hall_in, eeo_in, hallo_in, lsham_in, enim_in, &
                                                iz_in, kk_in, nb_in, nnmax_in, ntype_in, nmax_in, inva_in, b_in, &
                                                bee_out, bha_out, oee_out, oha_out, hons_out)
      class(cheb_cache_t), intent(inout) :: this
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

      valid = this%ortho%valid &
         .and. this%ortho%nb == nb_in &
         .and. this%ortho%nnmax == nnmax_in &
         .and. this%ortho%ntype == ntype_in &
         .and. this%ortho%nmax == nmax_in &
         .and. this%ortho%kk == kk_in &
         .and. this%ortho%inva == inva_in &
         .and. this%ortho%b == b_in

      if (.not. valid) then
         if (associated(this%bee_cache)) deallocate (this%bee_cache)
         if (associated(this%bha_cache)) deallocate (this%bha_cache)
         if (associated(this%oee_cache)) deallocate (this%oee_cache)
         if (associated(this%oha_cache)) deallocate (this%oha_cache)
         if (associated(this%hons_cache)) deallocate (this%hons_cache)

         ! Sweep-A bare h, scaled by inva (no on-site fold).
         allocate (this%bee_cache(nb_in, nb_in, nnmax_in, ntype_in))
         this%bee_cache = cmplx(real(ee_in, sp), aimag(ee_in), sp)*inva_in

         ! Sweep-B eeo, raw.
         allocate (this%oee_cache(nb_in, nb_in, nnmax_in, ntype_in))
         this%oee_cache = cmplx(real(eeo_in, sp), aimag(eeo_in), sp)

         ! On-site (enim + lsham - b*I)*inva, per type.
         allocate (this%hons_cache(nb_in, nb_in, ntype_in))
         this%hons_cache = (cmplx(real(enim_in, sp), aimag(enim_in), sp) + &
                       cmplx(real(lsham_in, sp), aimag(lsham_in), sp))*inva_in
         do t_in = 1, ntype_in
            do l_in = 1, nb_in
               this%hons_cache(l_in, l_in, t_in) = this%hons_cache(l_in, l_in, t_in) - b_in*inva_in
            end do
         end do

         if (nmax_in > 0) then
            allocate (this%bha_cache(nb_in, nb_in, nnmax_in, nmax_in))
            this%bha_cache = cmplx(real(hall_in(:, :, :, 1:nmax_in), sp), aimag(hall_in(:, :, :, 1:nmax_in)), sp)*inva_in
            allocate (this%oha_cache(nb_in, nb_in, nnmax_in, nmax_in))
            this%oha_cache = cmplx(real(hallo_in(:, :, :, 1:nmax_in), sp), aimag(hallo_in(:, :, :, 1:nmax_in)), sp)
            ! Note: the impurity on-site enim/lsham fold is carried by hons via
            ! iz(k), identical to the bulk path (ham_hoh_vec_matmul applies
            ! lsham/enim on-site for both local and bulk regions).
         end if

         this%ortho%nb = nb_in
         this%ortho%nnmax = nnmax_in
         this%ortho%ntype = ntype_in
         this%ortho%nmax = nmax_in
         this%ortho%kk = kk_in
         this%ortho%inva = inva_in
         this%ortho%b = b_in
         this%ortho%valid = .true.
      end if

      bee_out => this%bee_cache
      oee_out => this%oee_cache
      hons_out => this%hons_cache
      if (associated(this%bha_cache)) then
         bha_out => this%bha_cache
         oha_out => this%oha_cache
      else
         nullify (bha_out)
         nullify (oha_out)
      end if
   end subroutine cheb_cache_ensure_scaled_ortho_sp

   subroutine cheb_cache_ensure_scaled_bsr_sp(this, ee_in, hall_in, lsham_in, nn_in, iz_in, kk_in, nb_in, nnmax_in, ntype_in, nmax_in, inva_in, b_in, hval_out, hcol_out, hrow_out)
      class(cheb_cache_t), intent(inout) :: this
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

      valid = this%bsr%valid &
         .and. this%bsr%nb == nb_in &
         .and. this%bsr%nnmax == nnmax_in &
         .and. this%bsr%ntype == ntype_in &
         .and. this%bsr%nmax == nmax_in &
         .and. this%bsr%kk == kk_in &
         .and. this%bsr%inva == inva_in &
         .and. this%bsr%b == b_in

      if (.not. valid) then
         if (associated(this%hval_cache)) deallocate (this%hval_cache)
         if (associated(this%hcol_cache)) deallocate (this%hcol_cache)
         if (associated(this%hrow_cache)) deallocate (this%hrow_cache)

         nblocks = 0
         do atom = 1, kk_in
            num_neighbors = nn_in(atom, 1)
            do neighbor_idx = 1, num_neighbors
               if (neighbor_idx == 1 .or. nn_in(atom, neighbor_idx) /= 0) nblocks = nblocks + 1
            end do
         end do

         allocate (this%hval_cache(nb_in, nb_in, nblocks), this%hcol_cache(nblocks), this%hrow_cache(kk_in + 1))
         block_idx = 0
         this%hrow_cache(1) = 1

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
                  this%hval_cache(:, :, block_idx) = cmplx(real(hall_in(:, :, neighbor_idx, atom), sp), &
                                                      aimag(hall_in(:, :, neighbor_idx, atom)), sp)*inva_in
                  ih = iz_in(atom)
               else
                  ih = iz_in(atom)
                  this%hval_cache(:, :, block_idx) = cmplx(real(ee_in(:, :, neighbor_idx, ih), sp), &
                                                      aimag(ee_in(:, :, neighbor_idx, ih)), sp)*inva_in
               end if

               if (neighbor_idx == 1) then
                  this%hval_cache(:, :, block_idx) = this%hval_cache(:, :, block_idx) + &
                     cmplx(real(lsham_in(:, :, ih), sp), aimag(lsham_in(:, :, ih)), sp)*inva_in
                  do l_in = 1, nb_in
                     this%hval_cache(l_in, l_in, block_idx) = this%hval_cache(l_in, l_in, block_idx) - b_in*inva_in
                  end do
               end if
               this%hcol_cache(block_idx) = block_col
            end do
            if (atom < kk_in) this%hrow_cache(atom + 1) = block_idx + 1
         end do
         this%hrow_cache(kk_in + 1) = nblocks + 1

         this%bsr%nb = nb_in
         this%bsr%nnmax = nnmax_in
         this%bsr%ntype = ntype_in
         this%bsr%nmax = nmax_in
         this%bsr%kk = kk_in
         this%bsr%inva = inva_in
         this%bsr%b = b_in
         this%bsr%valid = .true.
      end if

      hval_out => this%hval_cache
      hcol_out => this%hcol_cache
      hrow_out => this%hrow_cache
   end subroutine cheb_cache_ensure_scaled_bsr_sp

   !> Raw fp32 mirrors of the velocity operators for the conductivity/orbital
   !> fast kernels. v_a/v_b are copied as-is; for hoh, vo_a/vo_b are copied with
   !> the on-site shell-1 block zeroed (velo_hoh_vec_matmul applies vo off-site
   !> only) and raw ee/hall are mirrored (the velocity hoh sweep uses the BARE,
   !> unscaled h). All raw (no 1/a, no -b*I, no lsham/enim fold).
   subroutine cheb_cache_ensure_scaled_velocity_sp(this, v_a, v_b, vo_a, vo_b, ee, hall, &
                                                   kk_in, nb_in, nnmax_in, ntype_in, nmax_in, hoh, &
                                                   fva, fvb, fvoa, fvob, fvee, fvha)
      class(cheb_cache_t), intent(inout) :: this
      integer, intent(in) :: kk_in, nb_in, nnmax_in, ntype_in, nmax_in
      complex(rp), intent(in) :: v_a(nb_in, nb_in, nnmax_in, ntype_in)
      complex(rp), intent(in) :: v_b(nb_in, nb_in, nnmax_in, ntype_in)
      complex(rp), intent(in) :: vo_a(nb_in, nb_in, nnmax_in, ntype_in)
      complex(rp), intent(in) :: vo_b(nb_in, nb_in, nnmax_in, ntype_in)
      complex(rp), intent(in) :: ee(nb_in, nb_in, nnmax_in, ntype_in)
      complex(rp), intent(in) :: hall(nb_in, nb_in, nnmax_in, *)
      logical, intent(in) :: hoh
      complex(sp), pointer, intent(out) :: fva(:, :, :, :), fvb(:, :, :, :)
      complex(sp), pointer, intent(out) :: fvoa(:, :, :, :), fvob(:, :, :, :)
      complex(sp), pointer, intent(out) :: fvee(:, :, :, :), fvha(:, :, :, :)
      logical :: valid
      integer :: t_in

      valid = this%vel%valid &
         .and. this%vel%nb == nb_in &
         .and. this%vel%nnmax == nnmax_in &
         .and. this%vel%ntype == ntype_in &
         .and. this%vel%nmax == nmax_in &
         .and. this%vel%kk == kk_in &
         .and. (this%vel%hoh .eqv. hoh)

      if (.not. valid) then
         if (associated(this%fva_cache)) deallocate (this%fva_cache)
         if (associated(this%fvb_cache)) deallocate (this%fvb_cache)
         if (associated(this%fvoa_cache)) deallocate (this%fvoa_cache)
         if (associated(this%fvob_cache)) deallocate (this%fvob_cache)
         if (associated(this%fvee_cache)) deallocate (this%fvee_cache)
         if (associated(this%fvha_cache)) deallocate (this%fvha_cache)

         allocate (this%fva_cache(nb_in, nb_in, nnmax_in, ntype_in))
         allocate (this%fvb_cache(nb_in, nb_in, nnmax_in, ntype_in))
         this%fva_cache = cmplx(real(v_a, sp), aimag(v_a), sp)
         this%fvb_cache = cmplx(real(v_b, sp), aimag(v_b), sp)

         if (hoh) then
            allocate (this%fvoa_cache(nb_in, nb_in, nnmax_in, ntype_in))
            allocate (this%fvob_cache(nb_in, nb_in, nnmax_in, ntype_in))
            this%fvoa_cache = cmplx(real(vo_a, sp), aimag(vo_a), sp)
            this%fvob_cache = cmplx(real(vo_b, sp), aimag(vo_b), sp)
            this%fvoa_cache(:, :, 1, :) = (0.0_sp, 0.0_sp)   ! on-site vo excluded
            this%fvob_cache(:, :, 1, :) = (0.0_sp, 0.0_sp)
            allocate (this%fvee_cache(nb_in, nb_in, nnmax_in, ntype_in))
            this%fvee_cache = cmplx(real(ee, sp), aimag(ee), sp)
            if (nmax_in > 0) then
               allocate (this%fvha_cache(nb_in, nb_in, nnmax_in, nmax_in))
               this%fvha_cache = cmplx(real(hall(:, :, :, 1:nmax_in), sp), aimag(hall(:, :, :, 1:nmax_in)), sp)
            end if
         end if

         this%vel%nb = nb_in
         this%vel%nnmax = nnmax_in
         this%vel%ntype = ntype_in
         this%vel%nmax = nmax_in
         this%vel%kk = kk_in
         this%vel%hoh = hoh
         this%vel%valid = .true.
      end if

      fva => this%fva_cache
      fvb => this%fvb_cache
      if (associated(this%fvoa_cache)) then
         fvoa => this%fvoa_cache
         fvob => this%fvob_cache
         fvee => this%fvee_cache
      else
         nullify (fvoa); nullify (fvob); nullify (fvee)
      end if
      if (associated(this%fvha_cache)) then
         fvha => this%fvha_cache
      else
         nullify (fvha)
      end if
   end subroutine cheb_cache_ensure_scaled_velocity_sp

   subroutine pack_psi0_sp(psi0, kk, nb, p0)
      integer, intent(in) :: kk, nb
      complex(rp), intent(in) :: psi0(nb, nb, kk)
      complex(sp), intent(out) :: p0(nb*kk, nb)
      integer :: k, l, c

      !$omp parallel do private(k, c, l) schedule(static)
      do k = 1, kk
         do c = 1, nb
            do l = 1, nb
               p0(l + nb*(k - 1), c) = cmplx(real(psi0(l, c, k), sp), aimag(psi0(l, c, k)), sp)
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine pack_psi0_sp

   !> Hermitian rank-k: C = A^H A via cherk, then fill the upper triangle.
   subroutine cherk_full_sp(n, kdim, Amat, C)
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
   end subroutine cherk_full_sp

   !> Shared per-shell fp32 block matvec (site-major):
   !>   y = alpha*(acc - bsc*x1)*inva + beta*x0,   acc = sum_neigh Op * x1_block
   !> Op is hee_op[iz(k)] (bulk) for k > nmax_op, hha_op[k] (impurity) otherwise.
   !> Pre-scaled operators use inva=1, bsc=0; raw operators carry scaling here.
   subroutine spmv_block_sp(hee_op, hha_op, nn, iz, kk, nb, nnmax, ntype, nmax_op, &
                            ld, x1, x0, y, alpha, beta, inva, bsc)
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax_op, ld
      complex(sp), intent(in) :: hee_op(nb, nb, nnmax, ntype)
      complex(sp), intent(in) :: hha_op(nb, nb, nnmax, *)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      complex(sp), intent(in) :: x1(ld, nb), x0(ld, nb)
      complex(sp), intent(out) :: y(ld, nb)
      real(sp), intent(in) :: alpha, beta, inva, bsc
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp)
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
            if (kk_ <= nmax_op) then
               call cgemm('N', 'N', nb, nb, nb, cone, hha_op(:, :, s_, kk_), nb, x1(r0 + 1, 1), ld, cone, acc, nb)
            else
               t_ = iz(kk_)
               call cgemm('N', 'N', nb, nb, nb, cone, hee_op(:, :, s_, t_), nb, x1(r0 + 1, 1), ld, cone, acc, nb)
            end if
         end do
         r0 = nb*(kk_ - 1)
         if (bsc /= 0.0_sp) acc = acc - bsc*x1(r0 + 1:r0 + nb, :)
         if (inva /= 1.0_sp) acc = acc*inva
         if (beta /= 0.0_sp) then
            y(r0 + 1:r0 + nb, :) = alpha*acc + beta*x0(r0 + 1:r0 + nb, :)
         else
            y(r0 + 1:r0 + nb, :) = alpha*acc
         end if
      end do
      !$omp end parallel do
   end subroutine spmv_block_sp

   !> In-place variant for the beta*x0 form when x0 and y are the same array.
   subroutine spmv_block_sp_inplace(hee_op, hha_op, nn, iz, kk, nb, nnmax, ntype, nmax_op, &
                                    ld, x1, y, alpha, beta, inva, bsc)
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax_op, ld
      complex(sp), intent(in) :: hee_op(nb, nb, nnmax, ntype)
      complex(sp), intent(in) :: hha_op(nb, nb, nnmax, *)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      complex(sp), intent(in) :: x1(ld, nb)
      complex(sp), intent(inout) :: y(ld, nb)
      real(sp), intent(in) :: alpha, beta, inva, bsc
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp)
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
            if (kk_ <= nmax_op) then
               call cgemm('N', 'N', nb, nb, nb, cone, hha_op(:, :, s_, kk_), nb, x1(r0 + 1, 1), ld, cone, acc, nb)
            else
               t_ = iz(kk_)
               call cgemm('N', 'N', nb, nb, nb, cone, hee_op(:, :, s_, t_), nb, x1(r0 + 1, 1), ld, cone, acc, nb)
            end if
         end do
         r0 = nb*(kk_ - 1)
         if (bsc /= 0.0_sp) acc = acc - bsc*x1(r0 + 1:r0 + nb, :)
         if (inva /= 1.0_sp) acc = acc*inva
         if (beta /= 0.0_sp) then
            y(r0 + 1:r0 + nb, :) = alpha*acc + beta*y(r0 + 1:r0 + nb, :)
         else
            y(r0 + 1:r0 + nb, :) = alpha*acc
         end if
      end do
      !$omp end parallel do
   end subroutine spmv_block_sp_inplace

   !> hoh combine for the H~ apply: y = alpha*(wt - e + hons*x1) + beta*x0,
   !> with wt = (h/a)x1 and e = eeo*wt (both precomputed); hons applied on-site.
   subroutine combine_hoh_sp(wt, e, hons, iz, kk, nb, ntype, ld, x1, x0, y, alpha, beta)
      integer, intent(in) :: kk, nb, ntype, ld
      complex(sp), intent(in) :: wt(ld, nb), e(ld, nb), x1(ld, nb), x0(ld, nb)
      complex(sp), intent(in) :: hons(nb, nb, ntype)
      integer, intent(in) :: iz(kk)
      complex(sp), intent(out) :: y(ld, nb)
      real(sp), intent(in) :: alpha, beta
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp)
      complex(sp) :: acc(nb, nb)
      integer :: kk_, t_, r0
      !$omp parallel do private(kk_, t_, r0, acc) schedule(dynamic, 32)
      do kk_ = 1, kk
         r0 = nb*(kk_ - 1)
         t_ = iz(kk_)
         acc = wt(r0 + 1:r0 + nb, :) - e(r0 + 1:r0 + nb, :)
         call cgemm('N', 'N', nb, nb, nb, cone, hons(:, :, t_), nb, x1(r0 + 1, 1), ld, cone, acc, nb)
         if (beta /= 0.0_sp) then
            y(r0 + 1:r0 + nb, :) = alpha*acc + beta*x0(r0 + 1:r0 + nb, :)
         else
            y(r0 + 1:r0 + nb, :) = alpha*acc
         end if
      end do
      !$omp end parallel do
   end subroutine combine_hoh_sp

   !> H operator sweep with null-impurity guard (impurity arrays are absent
   !> when nmax == 0): y = alpha*(Op x1 - bsc x1)*inva + beta*x0.
   subroutine hsweep_sp(opee, opha, nn, iz, kk, nb, nnmax, ntype, nmax, ld, x1, x0, y, al, be, inva, bsc)
      complex(sp), pointer, intent(in) :: opee(:, :, :, :), opha(:, :, :, :)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax, ld
      complex(sp), intent(in) :: x1(ld, nb), x0(ld, nb)
      complex(sp), intent(out) :: y(ld, nb)
      real(sp), intent(in) :: al, be, inva, bsc

      if (associated(opha)) then
         call spmv_block_sp(opee, opha, nn, iz, kk, nb, nnmax, ntype, nmax, ld, x1, x0, y, al, be, inva, bsc)
      else
         call spmv_block_sp(opee, opee, nn, iz, kk, nb, nnmax, ntype, 0, ld, x1, x0, y, al, be, inva, bsc)
      end if
   end subroutine hsweep_sp

   !> y = alpha*(H~ x1) + beta*x0, hoh-aware (two-sweep).
   subroutine happly_sp(do_hoh, hee, hha, bee, bha, oee, oha, hons, nn, iz, kk, nb, nnmax, ntype, nmax, &
                        ld, x1, x0, y, al, be, wt, etmp)
      logical, intent(in) :: do_hoh
      complex(sp), pointer, intent(in) :: hee(:, :, :, :), hha(:, :, :, :)
      complex(sp), pointer, intent(in) :: bee(:, :, :, :), bha(:, :, :, :)
      complex(sp), pointer, intent(in) :: oee(:, :, :, :), oha(:, :, :, :)
      complex(sp), pointer, intent(in) :: hons(:, :, :)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax, ld
      complex(sp), intent(in) :: x1(ld, nb), x0(ld, nb)
      complex(sp), intent(out) :: y(ld, nb)
      complex(sp), intent(inout) :: wt(ld, nb), etmp(ld, nb)
      real(sp), intent(in) :: al, be

      if (do_hoh) then
         call hsweep_sp(bee, bha, nn, iz, kk, nb, nnmax, ntype, nmax, ld, x1, x1, wt, 1.0_sp, 0.0_sp, 1.0_sp, 0.0_sp)
         call hsweep_sp(oee, oha, nn, iz, kk, nb, nnmax, ntype, nmax, ld, wt, wt, etmp, 1.0_sp, 0.0_sp, 1.0_sp, 0.0_sp)
         call combine_hoh_sp(wt, etmp, hons, iz, kk, nb, ntype, ld, x1, x0, y, al, be)
      else
         call hsweep_sp(hee, hha, nn, iz, kk, nb, nnmax, ntype, nmax, ld, x1, x0, y, al, be, 1.0_sp, 0.0_sp)
      end if
   end subroutine happly_sp

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
      complex(sp), pointer :: p0(:, :), p1(:, :), p2(:, :)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp)
      integer :: ld
      real(sp) :: inva, a_sp, b_sp
      logical :: do_hoh

      do_hoh = .false.
      if (present(hoh)) do_hoh = hoh

      ld = nb*kk
      a_sp = real(a, sp)
      b_sp = real(b, sp)
      inva = 1.0_sp/a_sp
      call cheb_cache%ensure_work_buffers(ld, nb, w0, w1, w2)
      p0 => w0
      p1 => w1
      p2 => w2

      if (do_hoh) then
         ! --- two-sweep hoh operands (bare h, eeo, on-site) + extra buffer ---
         call cheb_cache%ensure_hoh_buffer(ld, nb, wt)
         call cheb_cache%ensure_scaled_ortho_sp(ee, hall, eeo, hallo, lsham, enim, iz, kk, nb, &
                                                nnmax, ntype, nmax, inva, b_sp, bee, bha, oee, oha, hons)
      else
         ! --- 1. scaled operator copies: Ht = (H + lsham - b*I)/a ----------
         call cheb_cache%ensure_scaled_hamiltonian_sp(ee, hall, lsham, iz, kk, nb, nnmax, ntype, nmax, inva, b_sp, hee, hha)
      end if

      call pack_psi0_sp(psi0, kk, nb, p0)
      call run_three_term_selected()

   contains

      subroutine run_three_term_selected()
         complex(sp), pointer :: ptmp(:, :)
         complex(sp) :: mu1_sp(nb, nb), mu2_sp(nb, nb), dum_sp(nb, nb)
         integer :: ll

         call cherk_full_sp(nb, ld, p0, mu1_sp)
         mu(:, :, 1) = mu1_sp
         call apply_step_selected(p0, p0, p1, 1.0_sp, 0.0_sp, ld, nb)
         call cgemm('C', 'N', nb, nb, ld, cone, p0, ld, p1, ld, czero, mu2_sp, nb)
         mu(:, :, 2) = mu2_sp

         do ll = 1, lld
            call apply_step_selected(p1, p0, p2, 2.0_sp, -1.0_sp, ld, nb)
            call cherk_full_sp(nb, ld, p1, dum_sp)
            mu(:, :, 2*ll + 1) = 2.0_sp*dum_sp - mu1_sp
            call cgemm('C', 'N', nb, nb, ld, cone, p2, ld, p1, ld, czero, dum_sp, nb)
            mu(:, :, 2*ll + 2) = 2.0_sp*dum_sp - mu2_sp
            ptmp => p0
            p0 => p1
            p1 => p2
            p2 => ptmp
         end do
      end subroutine run_three_term_selected

      subroutine apply_step_selected(x1, x0, y, alpha, beta, ld_apply, nb_apply)
         integer, intent(in) :: ld_apply, nb_apply
         complex(sp), intent(in), target :: x1(ld_apply, nb_apply), x0(ld_apply, nb_apply)
         complex(sp), intent(out) :: y(ld_apply, nb_apply)
         real(sp), intent(in) :: alpha, beta

         if (do_hoh) then
            call apply_step_hoh(x1, x0, y, alpha, beta, wt)
         else
            call apply_step(x1, x0, y, alpha, beta, ld_apply, nb_apply)
         end if
      end subroutine apply_step_selected

      !> y = alpha*(Ht x1) + beta*x0, one fused sweep (site-major arrays)
      subroutine apply_step(x1, x0, y, alpha, beta, ld_apply, nb_apply)
         integer, intent(in) :: ld_apply, nb_apply
         complex(sp), intent(in), target :: x1(ld_apply, nb_apply), x0(ld_apply, nb_apply)
         complex(sp), intent(out) :: y(ld_apply, nb_apply)
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
      complex(sp), pointer :: p0(:, :), p1(:, :), p2(:, :)
      complex(sp), pointer :: block_products(:, :, :)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp)
      integer, pointer :: hcol(:), hrow(:)
      integer :: ld, nblocks
      real(sp) :: inva, a_sp, b_sp

      ld = nb*kk
      a_sp = real(a, sp)
      b_sp = real(b, sp)
      inva = 1.0_sp/a_sp
      call cheb_cache%ensure_work_buffers(ld, nb, w0, w1, w2)
      p0 => w0
      p1 => w1
      p2 => w2

      call cheb_cache%ensure_scaled_bsr_sp(ee, hall, lsham, nn, iz, kk, nb, nnmax, ntype, nmax, inva, b_sp, hval, hcol, hrow)
      nblocks = size(hcol)
      call cheb_cache%ensure_block_products(nb, nblocks, block_products)

      call pack_psi0_sp(psi0, kk, nb, p0)
      call run_three_term_bsr()

   contains

      subroutine run_three_term_bsr()
         complex(sp), pointer :: ptmp(:, :)
         complex(sp) :: mu1_sp(nb, nb), mu2_sp(nb, nb), dum_sp(nb, nb)
         integer :: ll

         call cherk_full_sp(nb, ld, p0, mu1_sp)
         mu(:, :, 1) = mu1_sp
         call apply_step_bsr(p0, p0, p1, 1.0_sp, 0.0_sp, ld, nb)
         call cgemm('C', 'N', nb, nb, ld, cone, p0, ld, p1, ld, czero, mu2_sp, nb)
         mu(:, :, 2) = mu2_sp

         do ll = 1, lld
            call apply_step_bsr(p1, p0, p2, 2.0_sp, -1.0_sp, ld, nb)
            call cherk_full_sp(nb, ld, p1, dum_sp)
            mu(:, :, 2*ll + 1) = 2.0_sp*dum_sp - mu1_sp
            call cgemm('C', 'N', nb, nb, ld, cone, p2, ld, p1, ld, czero, dum_sp, nb)
            mu(:, :, 2*ll + 2) = 2.0_sp*dum_sp - mu2_sp
            ptmp => p0
            p0 => p1
            p1 => p2
            p2 => ptmp
         end do
      end subroutine run_three_term_bsr

      !> y = alpha*(Ht x1) + beta*x0 using one packed GEMM per BSR row.
      subroutine apply_step_bsr(x1, x0, y, alpha, beta, ld_apply, nb_apply)
         integer, intent(in) :: ld_apply, nb_apply
         complex(sp), intent(in), target :: x1(ld_apply, nb_apply), x0(ld_apply, nb_apply)
         complex(sp), intent(out) :: y(ld_apply, nb_apply)
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

   end subroutine cheb_moments_fast_batched

   !> Optional Intel oneMKL cgemm_batch implementation.
   !> Compiled with -DUSE_MKL_BATCH and CMake option ENABLE_MKL_KERNELS=ON.
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
      write (*, '(A)') 'ERROR: cheb_backend=mkl_batch requires CMake option ENABLE_MKL_KERNELS=ON.'
      error stop
#else
      interface
         subroutine cblas_cgemm_batch(layout, transa, transb, m, n, k, alpha, a_array, lda, &
                                      b_array, ldb, beta, c_array, ldc, group_count, group_size) bind(C, name='cblas_cgemm_batch')
            import :: c_int, c_ptr, c_float_complex
            integer(c_int), value :: layout, group_count
            integer(c_int), intent(in) :: transa(*), transb(*), m(*), n(*), k(*), lda(*), ldb(*), ldc(*), group_size(*)
            complex(c_float_complex), intent(in) :: alpha(*), beta(*)
            type(c_ptr), intent(in) :: a_array(*), b_array(*)
            type(c_ptr), intent(inout) :: c_array(*)
         end subroutine cblas_cgemm_batch
      end interface

      complex(sp), pointer :: hval(:, :, :)
      complex(sp), pointer :: w0(:, :), w1(:, :), w2(:, :)
      complex(sp), pointer :: p0(:, :), p1(:, :), p2(:, :)
      complex(sp), pointer :: block_products(:, :, :)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp)
      integer, pointer :: hcol(:), hrow(:)
      type(c_ptr), pointer :: a_ptr(:), b_ptr(:), c_ptrs(:)
      integer :: ld, nblocks
      real(sp) :: inva, a_sp, b_sp

      ld = nb*kk
      a_sp = real(a, sp)
      b_sp = real(b, sp)
      inva = 1.0_sp/a_sp
      call cheb_cache%ensure_work_buffers(ld, nb, w0, w1, w2)
      p0 => w0
      p1 => w1
      p2 => w2

      call cheb_cache%ensure_scaled_bsr_sp(ee, hall, lsham, nn, iz, kk, nb, nnmax, ntype, nmax, inva, b_sp, hval, hcol, hrow)
      nblocks = size(hcol)
      call cheb_cache%ensure_block_products(nb, nblocks, block_products)
      call cheb_cache%ensure_mkl_batch_ptr_arrays(nblocks, a_ptr, b_ptr, c_ptrs)
      call init_mkl_batch_static_ptrs()

      call pack_psi0_sp(psi0, kk, nb, p0)
      call run_three_term_mkl_batch()

   contains

      subroutine run_three_term_mkl_batch()
         complex(sp), pointer :: ptmp(:, :)
         complex(sp) :: mu1_sp(nb, nb), mu2_sp(nb, nb), dum_sp(nb, nb)
         integer :: ll

         call cherk_full_sp(nb, ld, p0, mu1_sp)
         mu(:, :, 1) = mu1_sp
         call apply_step_mkl_batch(p0, p0, p1, 1.0_sp, 0.0_sp, ld, nb)
         call cgemm('C', 'N', nb, nb, ld, cone, p0, ld, p1, ld, czero, mu2_sp, nb)
         mu(:, :, 2) = mu2_sp

         do ll = 1, lld
            call apply_step_mkl_batch(p1, p0, p2, 2.0_sp, -1.0_sp, ld, nb)
            call cherk_full_sp(nb, ld, p1, dum_sp)
            mu(:, :, 2*ll + 1) = 2.0_sp*dum_sp - mu1_sp
            call cgemm('C', 'N', nb, nb, ld, cone, p2, ld, p1, ld, czero, dum_sp, nb)
            mu(:, :, 2*ll + 2) = 2.0_sp*dum_sp - mu2_sp
            ptmp => p0
            p0 => p1
            p1 => p2
            p2 => ptmp
         end do
      end subroutine run_three_term_mkl_batch

      subroutine init_mkl_batch_static_ptrs()
         integer :: block
         !$omp parallel do private(block) schedule(static)
         do block = 1, nblocks
            a_ptr(block) = c_loc(hval(1, 1, block))
            c_ptrs(block) = c_loc(block_products(1, 1, block))
         end do
         !$omp end parallel do
      end subroutine init_mkl_batch_static_ptrs

      subroutine apply_step_mkl_batch(x1, x0, y, alpha, beta, ld_apply, nb_apply)
         integer, intent(in) :: ld_apply, nb_apply
         complex(sp), intent(in), target :: x1(ld_apply, nb_apply), x0(ld_apply, nb_apply)
         complex(sp), intent(out) :: y(ld_apply, nb_apply)
         real(sp), intent(in) :: alpha, beta
         integer(c_int), parameter :: cblas_col_major = 102_c_int, cblas_no_trans = 111_c_int
         complex(c_float_complex) :: alpha_c(1), beta_c(1)
         complex(sp) :: acc(nb, nb)
         integer(c_int) :: group_count_c, group_size_c(1)
         integer(c_int) :: transa_c(1), transb_c(1), m_c(1), n_c(1), k_c(1), lda_c(1), ldb_c(1), ldc_c(1)
         integer :: block, row, col, r0

         !$omp parallel do private(block, col, r0) schedule(static)
         do block = 1, nblocks
            col = hcol(block)
            r0 = nb*(col - 1)
            b_ptr(block) = c_loc(x1(r0 + 1, 1))
         end do
         !$omp end parallel do

         group_count_c = 1_c_int
         group_size_c(1) = int(nblocks, c_int)
         transa_c(1) = cblas_no_trans
         transb_c(1) = cblas_no_trans
         m_c(1) = int(nb, c_int)
         n_c(1) = int(nb, c_int)
         k_c(1) = int(nb, c_int)
         lda_c(1) = int(nb, c_int)
         ldb_c(1) = int(ld, c_int)
         ldc_c(1) = int(nb, c_int)
         alpha_c(1) = cmplx(1.0_sp, 0.0_sp, c_float_complex)
         beta_c(1) = cmplx(0.0_sp, 0.0_sp, c_float_complex)

         call cblas_cgemm_batch(cblas_col_major, transa_c, transb_c, m_c, n_c, k_c, alpha_c, a_ptr, lda_c, &
                                b_ptr, ldb_c, beta_c, c_ptrs, ldc_c, group_count_c, group_size_c)

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

#endif
   end subroutine cheb_moments_fast_mkl_batch

   !> Optional Intel oneMKL Inspector-Executor Sparse BLAS implementation.
   !> Compiled with -DUSE_MKL_SPARSE and CMake option ENABLE_MKL_KERNELS=ON.
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
      write (*, '(A)') 'ERROR: cheb_backend=mkl_sparse requires CMake option ENABLE_MKL_KERNELS=ON.'
      error stop
#else
      type(sparse_matrix_t) :: mkl_A
      type(matrix_descr) :: descA
      complex(sp), pointer :: hval(:, :, :)
      complex(sp), pointer :: w0(:, :), w1(:, :), w2(:, :)
      complex(sp), pointer :: p0(:, :), p1(:, :), p2(:, :)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp)
      integer, pointer :: hcol(:), hrow(:)
      integer :: status, ld
      real(sp) :: inva, a_sp, b_sp

      ld = nb*kk
      a_sp = real(a, sp)
      b_sp = real(b, sp)
      inva = 1.0_sp/a_sp
      call cheb_cache%ensure_work_buffers(ld, nb, w0, w1, w2)
      p0 => w0
      p1 => w1
      p2 => w2

      call cheb_cache%ensure_scaled_bsr_sp(ee, hall, lsham, nn, iz, kk, nb, nnmax, ntype, nmax, inva, b_sp, hval, hcol, hrow)

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

      call pack_psi0_sp(psi0, kk, nb, p0)
      call run_three_term_mkl()

      status = mkl_sparse_destroy(mkl_A)
      call check_mkl_sparse_status(status, 'mkl_sparse_destroy')

   contains

      subroutine run_three_term_mkl()
         complex(sp), pointer :: ptmp(:, :)
         complex(sp) :: mu1_sp(nb, nb), mu2_sp(nb, nb), dum_sp(nb, nb)
         integer :: ll

         call cherk_full_sp(nb, ld, p0, mu1_sp)
         mu(:, :, 1) = mu1_sp
         call apply_step_mkl(p0, p0, p1, 1.0_sp, 0.0_sp, ld, nb)
         call cgemm('C', 'N', nb, nb, ld, cone, p0, ld, p1, ld, czero, mu2_sp, nb)
         mu(:, :, 2) = mu2_sp

         do ll = 1, lld
            call apply_step_mkl(p1, p0, p2, 2.0_sp, -1.0_sp, ld, nb)
            call cherk_full_sp(nb, ld, p1, dum_sp)
            mu(:, :, 2*ll + 1) = 2.0_sp*dum_sp - mu1_sp
            call cgemm('C', 'N', nb, nb, ld, cone, p2, ld, p1, ld, czero, dum_sp, nb)
            mu(:, :, 2*ll + 2) = 2.0_sp*dum_sp - mu2_sp
            ptmp => p0
            p0 => p1
            p1 => p2
            p2 => ptmp
         end do
      end subroutine run_three_term_mkl

      subroutine apply_step_mkl(x1, x0, y, alpha, beta, ld_apply, nb_apply)
         integer, intent(in) :: ld_apply, nb_apply
         complex(sp), intent(in), target :: x1(ld_apply, nb_apply), x0(ld_apply, nb_apply)
         complex(sp), intent(out) :: y(ld_apply, nb_apply)
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

   !> Fast-CPU stochastic conductivity moments for one reference vector:
   !>   mu_nm = sum_k left_m(k)^H (v_a T_{n-1}(H~) v_b |psiref>)(k),
   !>   left_m = T_{m-1}(H~)|psiref>.  Single precision; hoh-aware (two-sweep H~
   !>   and velocity v - vo*(h*.)). Mirrors compute_moments_stochastic.
   subroutine cheb_moments_stochastic_fast(psiref, ee, hall, lsham, nn, iz, kk, nb, &
                                           nnmax, ntype, nmax, cond_ll, a, b, mu_nm, &
                                           hoh, eeo, hallo, enim, v_a, v_b, vo_a, vo_b)
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax, cond_ll
      complex(rp), intent(in) :: psiref(nb, nb, kk), ee(nb, nb, nnmax, ntype)
      complex(rp), intent(in) :: hall(nb, nb, nnmax, *), lsham(nb, nb, ntype)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      real(rp), intent(in) :: a, b
      complex(rp), intent(out) :: mu_nm(nb, nb, cond_ll, cond_ll)
      logical, intent(in) :: hoh
      complex(rp), intent(in) :: eeo(nb, nb, nnmax, ntype), hallo(nb, nb, nnmax, *)
      complex(rp), intent(in) :: enim(nb, nb, ntype)
      complex(rp), intent(in) :: v_a(nb, nb, nnmax, ntype), v_b(nb, nb, nnmax, ntype)
      complex(rp), intent(in) :: vo_a(nb, nb, nnmax, ntype), vo_b(nb, nb, nnmax, ntype)
      ! Locals
      complex(sp), pointer :: hee(:, :, :, :), hha(:, :, :, :)
      complex(sp), pointer :: bee(:, :, :, :), bha(:, :, :, :), oee(:, :, :, :), oha(:, :, :, :), hons(:, :, :)
      complex(sp), pointer :: fva(:, :, :, :), fvb(:, :, :, :), fvoa(:, :, :, :), fvob(:, :, :, :)
      complex(sp), pointer :: fvee(:, :, :, :), fvha(:, :, :, :)
      complex(sp), allocatable :: leftst(:, :, :), p0(:, :), p1(:, :), p2(:, :), rr(:, :)
      complex(sp), allocatable :: wt(:, :), etmp(:, :), hwt(:, :)
      complex(sp) :: dum(nb, nb)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp)
      integer :: ld, k, c, l, n, m
      real(sp) :: inva, a_sp, b_sp
      logical :: do_hoh

      nullify (hee, hha, bee, bha, oee, oha, hons)
      nullify (fva, fvb, fvoa, fvob, fvee, fvha)
      do_hoh = hoh
      ld = nb*kk
      a_sp = real(a, sp); b_sp = real(b, sp); inva = 1.0_sp/a_sp

      if (do_hoh) then
         call cheb_cache%ensure_scaled_ortho_sp(ee, hall, eeo, hallo, lsham, enim, iz, kk, nb, &
            nnmax, ntype, nmax, inva, b_sp, bee, bha, oee, oha, hons)
      else
         call cheb_cache%ensure_scaled_hamiltonian_sp(ee, hall, lsham, iz, kk, nb, nnmax, ntype, nmax, inva, b_sp, hee, hha)
      end if
      call cheb_cache%ensure_scaled_velocity_sp(v_a, v_b, vo_a, vo_b, ee, hall, kk, nb, nnmax, ntype, nmax, &
         do_hoh, fva, fvb, fvoa, fvob, fvee, fvha)

      allocate (leftst(ld, nb, cond_ll), p0(ld, nb), p1(ld, nb), p2(ld, nb), rr(ld, nb))
      allocate (wt(ld, nb), etmp(ld, nb), hwt(ld, nb))

      !$omp parallel do private(k, c, l) schedule(static)
      do k = 1, kk
         do c = 1, nb
            do l = 1, nb
               p0(l + nb*(k - 1), c) = cmplx(real(psiref(l, c, k), sp), aimag(psiref(l, c, k)), sp)
            end do
         end do
      end do
      !$omp end parallel do

      ! --- left states L_m = T_{m-1}(H~)|psiref> -------------------------
      leftst(:, :, 1) = p0
      p1 = p0
      do m = 2, cond_ll
         if (m == 2) then
            call happly_sp(do_hoh, hee, hha, bee, bha, oee, oha, hons, nn, iz, kk, nb, nnmax, ntype, nmax, &
                           ld, p1, p1, p2, 1.0_sp, 0.0_sp, wt, etmp)
         else
            call happly_sp(do_hoh, hee, hha, bee, bha, oee, oha, hons, nn, iz, kk, nb, nnmax, ntype, nmax, &
                           ld, p1, p0, p2, 2.0_sp, -1.0_sp, wt, etmp)
         end if
         p0 = p1; p1 = p2
         leftst(:, :, m) = p1
      end do

      ! --- right recursion v_n = T_{n-1}(H~) (v_b|psiref>); R = v_a v_n ---
      call vapply(fvb, fvob, leftst(:, :, 1), p0)   ! p0 = v_b psiref
      do n = 1, cond_ll
         if (n == 1) then
            p1 = p0
         else if (n == 2) then
            p0 = p1
            call happly_sp(do_hoh, hee, hha, bee, bha, oee, oha, hons, nn, iz, kk, nb, nnmax, ntype, nmax, &
                           ld, p0, p0, p1, 1.0_sp, 0.0_sp, wt, etmp)
         else
            call happly_sp(do_hoh, hee, hha, bee, bha, oee, oha, hons, nn, iz, kk, nb, nnmax, ntype, nmax, &
                           ld, p1, p0, p2, 2.0_sp, -1.0_sp, wt, etmp)
            p0 = p1; p1 = p2
         end if
         call vapply(fva, fvoa, p1, rr)              ! rr = v_a v_n
         do m = 1, cond_ll
            call cgemm('C', 'N', nb, nb, ld, cone, leftst(1, 1, m), ld, rr, ld, czero, dum, nb)
            mu_nm(:, :, n, m) = dum
         end do
      end do

      deallocate (leftst, p0, p1, p2, rr, wt, etmp, hwt)

   contains

      !> Velocity apply (raw): y = v*x, hoh -> y = v*x - vo*(h_bare*x).
      subroutine vapply(fv, fvo, x, y)
         complex(sp), pointer, intent(in) :: fv(:, :, :, :), fvo(:, :, :, :)
         complex(sp), intent(in) :: x(ld, nb)
         complex(sp), intent(out) :: y(ld, nb)
         call spmv_block_sp(fv, fv, nn, iz, kk, nb, nnmax, ntype, 0, ld, x, x, y, 1.0_sp, 0.0_sp, 1.0_sp, 0.0_sp)
         if (do_hoh) then
            call hsweep_sp(fvee, fvha, nn, iz, kk, nb, nnmax, ntype, nmax, ld, x, x, hwt, 1.0_sp, 0.0_sp, 1.0_sp, 0.0_sp)
            call spmv_block_sp_inplace(fvo, fvo, nn, iz, kk, nb, nnmax, ntype, 0, ld, hwt, y, -1.0_sp, 1.0_sp, 1.0_sp, 0.0_sp)
         end if
      end subroutine vapply

   end subroutine cheb_moments_stochastic_fast

   !> Fast-CPU orbital moments for one reference vector:
   !>   mu_n = sum_k left(k)^H T_{n-1}(H~)|psiref>(k),  left fixed (host-built).
   !>   Single precision; hoh-aware. Mirrors chebyshev_orbital_mod.
   subroutine cheb_moments_orbital_fast(left, psiref, ee, hall, lsham, nn, iz, kk, nb, &
                                        nnmax, ntype, nmax, lld, a, b, mu, &
                                        hoh, eeo, hallo, enim)
      integer, intent(in) :: kk, nb, nnmax, ntype, nmax, lld
      complex(rp), intent(in) :: left(nb, nb, kk), psiref(nb, nb, kk)
      complex(rp), intent(in) :: ee(nb, nb, nnmax, ntype), hall(nb, nb, nnmax, *), lsham(nb, nb, ntype)
      integer, intent(in) :: nn(kk, nnmax), iz(kk)
      real(rp), intent(in) :: a, b
      complex(rp), intent(out) :: mu(nb, nb, lld)
      logical, intent(in) :: hoh
      complex(rp), intent(in) :: eeo(nb, nb, nnmax, ntype), hallo(nb, nb, nnmax, *)
      complex(rp), intent(in) :: enim(nb, nb, ntype)
      ! Locals
      complex(sp), pointer :: hee(:, :, :, :), hha(:, :, :, :)
      complex(sp), pointer :: bee(:, :, :, :), bha(:, :, :, :), oee(:, :, :, :), oha(:, :, :, :), hons(:, :, :)
      complex(sp), allocatable :: lp(:, :), p0(:, :), p1(:, :), p2(:, :), wt(:, :), etmp(:, :)
      complex(sp) :: dum(nb, nb)
      complex(sp), parameter :: cone = (1.0_sp, 0.0_sp), czero = (0.0_sp, 0.0_sp)
      integer :: ld, k, c, l, n
      real(sp) :: inva, a_sp, b_sp
      logical :: do_hoh

      nullify (hee, hha, bee, bha, oee, oha, hons)
      do_hoh = hoh
      ld = nb*kk
      a_sp = real(a, sp); b_sp = real(b, sp); inva = 1.0_sp/a_sp

      if (do_hoh) then
         call cheb_cache%ensure_scaled_ortho_sp(ee, hall, eeo, hallo, lsham, enim, iz, kk, nb, &
            nnmax, ntype, nmax, inva, b_sp, bee, bha, oee, oha, hons)
      else
         call cheb_cache%ensure_scaled_hamiltonian_sp(ee, hall, lsham, iz, kk, nb, nnmax, ntype, nmax, inva, b_sp, hee, hha)
      end if

      allocate (lp(ld, nb), p0(ld, nb), p1(ld, nb), p2(ld, nb), wt(ld, nb), etmp(ld, nb))

      !$omp parallel do private(k, c, l) schedule(static)
      do k = 1, kk
         do c = 1, nb
            do l = 1, nb
               lp(l + nb*(k - 1), c) = cmplx(real(left(l, c, k), sp), aimag(left(l, c, k)), sp)
               p1(l + nb*(k - 1), c) = cmplx(real(psiref(l, c, k), sp), aimag(psiref(l, c, k)), sp)
            end do
         end do
      end do
      !$omp end parallel do

      do n = 1, lld
         if (n == 1) then
            continue                            ! p1 = psiref
         else if (n == 2) then
            p0 = p1
            call happly_sp(do_hoh, hee, hha, bee, bha, oee, oha, hons, nn, iz, kk, nb, nnmax, ntype, nmax, &
                           ld, p0, p0, p1, 1.0_sp, 0.0_sp, wt, etmp)
         else
            call happly_sp(do_hoh, hee, hha, bee, bha, oee, oha, hons, nn, iz, kk, nb, nnmax, ntype, nmax, &
                           ld, p1, p0, p2, 2.0_sp, -1.0_sp, wt, etmp)
            p0 = p1; p1 = p2
         end if
         call cgemm('C', 'N', nb, nb, ld, cone, lp, ld, p1, ld, czero, dum, nb)
         mu(:, :, n) = dum
      end do

      deallocate (lp, p0, p1, p2, wt, etmp)
   end subroutine cheb_moments_orbital_fast

end module chebyshev_fast_mod
