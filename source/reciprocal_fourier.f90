submodule (reciprocal_mod) reciprocal_fourier
   implicit none

contains

   module subroutine calculate_structure_factors(this, k_vec, structure_factors)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: structure_factors  ! (nr_max, ntype)
      ! Local variables
      integer :: ntype, ineigh, ia, nr
      real(rp) :: k_dot_r
      real(rp), dimension(3) :: r_vec

      ! Loop over all atom types
      do ntype = 1, this%lattice%ntype
         ia = this%lattice%atlist(ntype)
         nr = this%lattice%nn(ia, 1)  ! Number of neighbors for this type

         ! Calculate structure factors exp(i*k·R) for each neighbor
         do ineigh = 1, min(nr, size(structure_factors, 1))
            if (ineigh == 1) then
               ! On-site term (R = 0)
               structure_factors(ineigh, ntype) = cmplx(1.0_rp, 0.0_rp, rp)
            else
               ! Off-site terms - get neighbor vector in fractional coordinates
               ! Use ham_vec_type_direct which is populated by build_neighbor_vectors
               if (allocated(this%ham_vec_type_direct)) then
                  r_vec(1:3) = this%ham_vec_type_direct(1:3, ineigh, ntype)
               else
                  call g_logger%error('calculate_structure_factors: ham_vec_type_direct not allocated!', __FILE__, __LINE__)
                  r_vec = 0.0_rp
               end if

               ! Phase factor: exp(i * 2π * k_frac · R_frac)
               ! k_vec is in fractional coordinates (dimensionless)
               ! r_vec is in fractional coordinates (dimensionless)
               k_dot_r = 2.0_rp * pi * dot_product(k_vec, r_vec)
               structure_factors(ineigh, ntype) = cmplx(cos(k_dot_r), sin(k_dot_r), rp)
            end if
         end do
      end do
   end subroutine calculate_structure_factors

   module subroutine fourier_transform_hamiltonian(this, k_vec, hk_result)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: hk_result  ! (n_orb*n_sites, n_orb*n_sites)

      ! First-order k-space Hamiltonian: H(k) = h(k) = Sum_R ee(R) exp(i*k·R).
      ! NOTE: This deliberately reproduces the historical first-order behaviour.
      !       The on-site E_nu (enim) and spin-orbit (lsham) terms are NOT added
      !       here; they are only included in the second-order path
      !       (fourier_transform_hamiltonian_second_order). See kspace_ham_order.
      call this%fourier_transform_array(this%hamiltonian%ee, k_vec, hk_result)
      if (this%hamiltonian%ccor_2c) then
         block
            complex(rp), allocatable :: hcck(:, :)
            allocate(hcck(size(hk_result, 1), size(hk_result, 2)))
            call this%fourier_transform_array(this%hamiltonian%eecc, k_vec, hcck)
            hk_result(:, :) = hk_result(:, :) + hcck(:, :)
            deallocate(hcck)
         end block
      end if
   end subroutine fourier_transform_hamiltonian

   module subroutine fourier_transform_array(this, array4d, k_vec, mk_result)
      class(reciprocal), intent(in) :: this
      complex(rp), dimension(:, :, :, :), intent(in) :: array4d  ! (nb,nb,neigh,ntype)
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: mk_result  ! (n_orb*n_sites, n_orb*n_sites)
      ! Local variables
      integer :: isite, jsite, ntype_i, ineigh, ia, ja, nr
      integer :: i_start, i_end, j_start, j_end
      integer :: n_orb, n_sites
      complex(rp), dimension(:, :), allocatable :: structure_factors  ! (nr_max, ntype)

      n_orb = nb
      n_sites = this%lattice%nrec

      allocate(structure_factors(this%lattice%nn_max, this%lattice%ntype))
      call this%calculate_structure_factors(k_vec, structure_factors)

      mk_result = cmplx(0.0_rp, 0.0_rp, rp)

      ! Loop over all sites in the unit cell.
      ! For each i_site -> j_site pair, sum over lattice vectors R.
      do isite = 1, n_sites
         ntype_i = this%lattice%ib(isite)  ! Type of site i
         ia = this%lattice%atlist(ntype_i) ! Cluster atom for this type
         nr = this%lattice%nn(ia, 1)       ! Number of neighbors

         i_start = (isite - 1) * n_orb + 1
         i_end = isite * n_orb

         do ineigh = 1, nr
            if (ineigh == 1) then
               jsite = isite                       ! On-site: R = 0
            else
               ja = this%lattice%nn(ia, ineigh)     ! Cluster atom index
               if (ja < 1 .or. ja > this%lattice%kk) cycle
               jsite = this%lattice%iz(ja)          ! Map to unit cell site
               if (jsite < 1 .or. jsite > n_sites) cycle
            end if

            j_start = (jsite - 1) * n_orb + 1
            j_end = jsite * n_orb

            mk_result(i_start:i_end, j_start:j_end) = &
               mk_result(i_start:i_end, j_start:j_end) + &
               array4d(:, :, ineigh, ntype_i) * structure_factors(ineigh, ntype_i)
         end do
      end do

      deallocate(structure_factors)
   end subroutine fourier_transform_array

   module subroutine fourier_transform_hamiltonian_second_order(this, k_vec, hk_result)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: hk_result  ! (n_orb*n_sites, n_orb*n_sites)
      ! Local variables
      integer :: isite, ntype_i, i_start, i_end
      integer :: n_orb, n_sites, ndim
      complex(rp), dimension(:, :), allocatable :: hk, eeok, hohk, hcck

      n_orb = nb
      n_sites = this%lattice%nrec
      ndim = n_orb * n_sites

      allocate(hk(ndim, ndim), eeok(ndim, ndim), hohk(ndim, ndim), hcck(ndim, ndim))

      ! h(k) = Sum_R ee(R) exp(i*k·R)
      call this%fourier_transform_array(this%hamiltonian%ee, k_vec, hk)
      ! eeo(k) = Sum_R eeo(R) exp(i*k·R)
      call this%fourier_transform_array(this%hamiltonian%eeo, k_vec, eeok)

      ! [hoh](k) = eeo(k) · h(k)
      ! GPU: this zgemm (and the two Bloch sums above) are the natural offload
      ! point; batch over k with a strided-batched cuBLAS zgemm, mirroring the
      ! existing rsrec_cuda plugin pattern.
      call zgemm('n', 'n', ndim, ndim, ndim, cone, eeok, ndim, hk, ndim, czero, hohk, ndim)

      ! H(k) = h(k) - [hoh](k)
      hk_result = hk - hohk
      if (this%hamiltonian%ccor_2c) then
         call this%fourier_transform_array(this%hamiltonian%eecc, k_vec, hcck)
         hk_result = hk_result + hcck
      end if

      ! Add on-site (R=0) terms: E_nu (enim) and L.S (lsham), per site/type,
      ! onto the block diagonal (no phase factor for R=0).
      do isite = 1, n_sites
         ntype_i = this%lattice%ib(isite)
         i_start = (isite - 1) * n_orb + 1
         i_end = isite * n_orb
         hk_result(i_start:i_end, i_start:i_end) = &
            hk_result(i_start:i_end, i_start:i_end) + &
            this%hamiltonian%enim(:, :, ntype_i) + &
            this%hamiltonian%lsham(:, :, ntype_i)
      end do

      deallocate(hk, eeok, hohk, hcck)
   end subroutine fourier_transform_hamiltonian_second_order

   module subroutine fourier_transform_overlap(this, k_vec, sk_result)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: sk_result
      integer :: isite, jsite, ntype_i, ineigh, ia, ja, nr, iorb
      integer :: i_start, i_end, j_start, j_end
      integer :: n_orb, n_sites
      complex(rp), dimension(:, :), allocatable :: structure_factors
      complex(rp), dimension(:, :), allocatable :: overlap_block

      n_orb = nb
      n_sites = this%lattice%nrec
      allocate(structure_factors(this%lattice%nn_max, this%lattice%ntype))
      allocate(overlap_block(n_orb, n_orb))
      call this%calculate_structure_factors(k_vec, structure_factors)

      sk_result = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, n_sites
         ntype_i = this%lattice%ib(isite)
         ia = this%lattice%atlist(ntype_i)
         nr = this%lattice%nn(ia, 1)
         i_start = (isite - 1) * n_orb + 1
         i_end = isite * n_orb
         do ineigh = 1, nr
            if (ineigh == 1) then
               jsite = isite
            else
               ja = this%lattice%nn(ia, ineigh)
               if (ja < 1 .or. ja > this%lattice%kk) cycle
               jsite = this%lattice%iz(ja)
               if (jsite < 1 .or. jsite > n_sites) cycle
            end if
            j_start = (jsite - 1) * n_orb + 1
            j_end = jsite * n_orb
            overlap_block(:, :) = this%hamiltonian%eeo(:, :, ineigh, ntype_i)
            if (ineigh == 1 .and. jsite == isite) then
               do iorb = 1, n_orb
                  overlap_block(iorb, iorb) = overlap_block(iorb, iorb) + cmplx(1.0_rp, 0.0_rp, rp)
               end do
            end if
            sk_result(i_start:i_end, j_start:j_end) = sk_result(i_start:i_end, j_start:j_end) + &
               overlap_block * structure_factors(ineigh, ntype_i)
         end do
      end do
      deallocate(overlap_block, structure_factors)
   end subroutine fourier_transform_overlap

   module subroutine build_kspace_hamiltonian(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ik, ik_global, nk, ntype
      character(len=200) :: debug_msg
      logical :: using_kpath, distribute_mesh
      logical :: use_second_order
      integer :: i, j

      ! Determine which k-point set to use
      using_kpath = .false.
      if (allocated(this%k_path)) then
         ! Use k-path for band structure
         nk = this%nk_path
         using_kpath = .true.
         call root_info('reciprocal%build_kspace_hamiltonian: Building H(k) for k-path', __FILE__, __LINE__)
      else if (allocated(this%k_points)) then
         ! Use k-mesh for DOS/SCF
         nk = this%nk_total
         call root_info('reciprocal%build_kspace_hamiltonian: Building H(k) for k-mesh', __FILE__, __LINE__)
      else
         call g_logger%error('reciprocal%build_kspace_hamiltonian: No k-points generated. ' // &
                           'Call generate_mp_mesh or generate k-path first.', __FILE__, __LINE__)
         return
      end if

      write(debug_msg, '(A,I0,A,I0,A)') 'build_kspace_hamiltonian: Building for ', nk, &
                                        ' k-points and ', this%lattice%ntype, ' atom types'
      call root_info(trim(debug_msg), __FILE__, __LINE__)

      distribute_mesh = (.not. using_kpath) .and. (trim(this%dos_method) == 'gaussian')
      call setup_k_mesh_distribution(this, nk, distribute_mesh)

      ! Build neighbor vectors for each atom type (required for multi-site H_k)
      call this%build_neighbor_vectors()
      !!! print *,'KSPACE neighbours'
      !!! do i=1, this%lattice%ntype
      !!!    print *,'Type ', i, ' nneigh ', this%lattice%nn(this%lattice%atlist(i),1)
      !!!    do j=2, this%lattice%nn(this%lattice%atlist(i),1)
      !!!       print '(a,2i4, a, 3f10.6)','  Neighbour ', this%lattice%nn(i, j), this%lattice%iz(this%lattice%nn(i,j)), ': ', this%ham_vec_type(1,j,i), &
      !!!                this%ham_vec_type(2,j,i), this%ham_vec_type(3,j,i)
      !!!    end do
      !!! end do
      !!! print *,'================================='

      ! Allocate k-space Hamiltonian for multi-site system
      ! Dimension: (n_orb * n_sites) x (n_orb * n_sites) x n_kpoints
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.hk_bulk', this%hk_bulk, &
                                [this%max_orbs*this%lattice%nrec, this%max_orbs*this%lattice%nrec, this%nk_local])
#else
   if (allocated(this%hk_bulk)) deallocate(this%hk_bulk)
   allocate(this%hk_bulk(this%max_orbs*this%lattice%nrec, this%max_orbs*this%lattice%nrec, this%nk_local))
#endif

      ! Decide first- vs second-order k-space Hamiltonian.
      ! Resolution of kspace_ham_order:
      !   'auto'   -> 'second' if hoh is enabled, else 'first' (default)
      !   'second' -> force second order
      !   'first'  -> force first order
      ! Second order requires the orthogonalization arrays (eeo/enim/lsham),
      ! which are only built when the real-space Hamiltonian was assembled with
      ! hoh enabled. Fall back to first order with a warning otherwise.
      select case (trim(this%kspace_ham_order))
      case ('second')
         use_second_order = .true.
      case ('first')
         use_second_order = .false.
      case default  ! 'auto'
         use_second_order = this%hamiltonian%hoh
         if (use_second_order) then
            call root_info('build_kspace_hamiltonian: kspace_ham_order=auto and hoh enabled ' // &
               '-> using second order', __FILE__, __LINE__)
         end if
      end select
      if (use_second_order) then
         if (.not. this%hamiltonian%hoh .or. .not. allocated(this%hamiltonian%eeo)) then
            call g_logger%warning('build_kspace_hamiltonian: second-order H(k) requires ' // &
               'hoh-enabled Hamiltonian (eeo not available). Falling back to first order.', __FILE__, __LINE__)
            use_second_order = .false.
         else
            call root_info('build_kspace_hamiltonian: using SECOND-order H(k) = E_nu + h - hoh + L.S', __FILE__, __LINE__)
         end if
      end if

      ! Parallelize over k-points (coarse-grained) when OpenMP is available.
#ifdef _OPENMP
      !$omp parallel do private(ik, ik_global) shared(this, using_kpath, use_second_order) default(none)
#endif
      do ik = 1, this%nk_local
         ik_global = local_k_index_to_global(this, ik)
         ! Fourier transform: builds full multi-site H(k)
         if (use_second_order) then
            if (using_kpath) then
               call this%fourier_transform_hamiltonian_second_order(this%k_path(:, ik), this%hk_bulk(:, :, ik))
            else
               call this%fourier_transform_hamiltonian_second_order(this%k_points(:, ik_global), this%hk_bulk(:, :, ik))
            end if
         else
            if (using_kpath) then
               call this%fourier_transform_hamiltonian(this%k_path(:, ik), this%hk_bulk(:, :, ik))
            else
               call this%fourier_transform_hamiltonian(this%k_points(:, ik_global), this%hk_bulk(:, :, ik))
            end if
         end if
      end do
#ifdef _OPENMP
      !$omp end parallel do
#endif

      if (this%hamiltonian%ccor_2c .and. this%hamiltonian%ccor_debug .and. this%nk_local > 0) then
         block
            complex(rp), allocatable :: hcck_diag(:, :)
            real(rp) :: max_eecc, max_hcck, max_hk
            integer :: diag_k_global
            character(len=300) :: ccor_msg

            allocate(hcck_diag(size(this%hk_bulk, 1), size(this%hk_bulk, 2)))
            diag_k_global = local_k_index_to_global(this, 1)
            if (using_kpath) then
               call this%fourier_transform_array(this%hamiltonian%eecc, this%k_path(:, 1), hcck_diag)
            else
               call this%fourier_transform_array(this%hamiltonian%eecc, this%k_points(:, diag_k_global), hcck_diag)
            end if
            max_eecc = maxval(abs(this%hamiltonian%eecc))
            max_hcck = maxval(abs(hcck_diag))
            max_hk = maxval(abs(this%hk_bulk(:, :, 1)))
            write(ccor_msg, '(a,i0,a,es12.4,a,es12.4,a,es12.4)') &
               'CCOR2C k-space diagnostic at local k=1/global k=', diag_k_global, &
               ': maxabs(eecc)=', max_eecc, ', maxabs(Hcc(k))=', max_hcck, &
               ', maxabs(H(k))=', max_hk
            call root_info(trim(ccor_msg), __FILE__, __LINE__)
            deallocate(hcck_diag)
         end block
      end if

      call root_info('reciprocal%build_kspace_hamiltonian: K-space Hamiltonian built successfully', __FILE__, __LINE__)
      if (trim(this%reciprocal_mode) == 'generalized_overlap_proxy') then
         call this%build_kspace_overlap()
      else if (trim(this%reciprocal_mode) == 'generalized_overlap_kanpur') then
         call g_logger%warning('reciprocal%build_kspace_hamiltonian: generalized_overlap_kanpur requested but not implemented yet. Using ham_only solve path.', __FILE__, __LINE__)
      end if
      
      ! Diagnostic: Check H(k) at Gamma point for multi-site systems
      if (this%lattice%nrec > 1 .and. this%nk_local > 0 .and. local_k_index_to_global(this, 1) == 1) then
         call this%check_multisite_hamiltonian_diagonal()
      end if
   end subroutine build_kspace_hamiltonian

   module subroutine build_kspace_overlap(this)
      class(reciprocal), intent(inout) :: this
      integer :: ik, ik_global, nk
      logical :: using_kpath

      using_kpath = .false.
      if (allocated(this%k_path)) then
         nk = this%nk_path
         using_kpath = .true.
      else if (allocated(this%k_points)) then
         nk = this%nk_total
      else
         call g_logger%error('build_kspace_overlap: No k-point set available.', __FILE__, __LINE__)
         return
      end if

      call setup_k_mesh_distribution(this, nk, (.not. using_kpath) .and. (trim(this%dos_method) == 'gaussian'))

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.sk_overlap', this%sk_overlap, [this%max_orbs*this%lattice%nrec, this%max_orbs*this%lattice%nrec, this%nk_local])
#else
      if (allocated(this%sk_overlap)) deallocate(this%sk_overlap)
      allocate(this%sk_overlap(this%max_orbs*this%lattice%nrec, this%max_orbs*this%lattice%nrec, this%nk_local))
#endif

#ifdef _OPENMP
      !$omp parallel do private(ik, ik_global) shared(this, using_kpath) default(none)
#endif
      do ik = 1, this%nk_local
         ik_global = local_k_index_to_global(this, ik)
         if (using_kpath) then
            call this%fourier_transform_overlap(this%k_path(:, ik), this%sk_overlap(:, :, ik))
         else
            call this%fourier_transform_overlap(this%k_points(:, ik_global), this%sk_overlap(:, :, ik))
         end if
      end do
#ifdef _OPENMP
      !$omp end parallel do
#endif
      call root_info('reciprocal%build_kspace_overlap: Built S(k) overlap proxy.', __FILE__, __LINE__)
   end subroutine build_kspace_overlap

module subroutine check_multisite_hamiltonian_diagonal(this)
   class(reciprocal), intent(in) :: this
   integer :: isite, iorb, ispin
   complex(rp) :: h_avg_site1_up, h_avg_site1_dn, h_avg_site2_up, h_avg_site2_dn
   integer :: n_sites, idx_up, idx_dn
   character(len=256) :: msg
   
   n_sites = this%lattice%nrec
   if (n_sites < 2) return
   
   ! Check on-site diagonal elements for first two sites at Gamma (ik=1)
   ! For spd basis: s, p, d for spin-up (indices 1-9), then spin-down (10-18)
   
   ! Site 1, spin-up d-orbital average (indices 5-9 in spin block)
   h_avg_site1_up = 0.0_rp
   do iorb = 5, 9
      idx_up = iorb  ! First 9 are spin-up
      h_avg_site1_up = h_avg_site1_up + this%hk_bulk(idx_up, idx_up, 1)
   end do
   h_avg_site1_up = h_avg_site1_up / 5.0_rp
   
   ! Site 1, spin-down d-orbital average (indices 14-18 in site block)
   h_avg_site1_dn = 0.0_rp
   do iorb = 5, 9
      idx_dn = iorb + 9  ! Next 9 are spin-down
      h_avg_site1_dn = h_avg_site1_dn + this%hk_bulk(idx_dn, idx_dn, 1)
   end do
   h_avg_site1_dn = h_avg_site1_dn / 5.0_rp
   
   ! Site 2, spin-up d-orbital average (indices 18+5 to 18+9)
   h_avg_site2_up = 0.0_rp
   do iorb = 5, 9
      idx_up = 18 + iorb  ! Site 2 block starts at index 19
      h_avg_site2_up = h_avg_site2_up + this%hk_bulk(idx_up, idx_up, 1)
   end do
   h_avg_site2_up = h_avg_site2_up / 5.0_rp
   
   ! Site 2, spin-down d-orbital average
   h_avg_site2_dn = 0.0_rp
   do iorb = 5, 9
      idx_dn = 18 + iorb + 9  ! Site 2, spin-down block
      h_avg_site2_dn = h_avg_site2_dn + this%hk_bulk(idx_dn, idx_dn, 1)
   end do
   h_avg_site2_dn = h_avg_site2_dn / 5.0_rp
   
   write(msg, '(A,2F12.6)') 'H(k=0) diagonal check - Site 1 d-orbital avg (up/dn): ', &
                             real(h_avg_site1_up), real(h_avg_site1_dn)
   call root_info(trim(msg), __FILE__, __LINE__)
   
   write(msg, '(A,2F12.6)') 'H(k=0) diagonal check - Site 2 d-orbital avg (up/dn): ', &
                             real(h_avg_site2_up), real(h_avg_site2_dn)
   call root_info(trim(msg), __FILE__, __LINE__)
   
   write(msg, '(A,F12.6)') 'H(k=0) diagonal difference (site2-site1) up: ', &
                            real(h_avg_site2_up - h_avg_site1_up)
   call root_info(trim(msg), __FILE__, __LINE__)
   
end subroutine check_multisite_hamiltonian_diagonal

module subroutine check_hamiltonian_hermiticity(this, ik)
   class(reciprocal), intent(in) :: this
   integer, intent(in) :: ik
   
   integer :: i, j, n
   real(rp) :: max_diff, diff
   complex(rp) :: h_ij, h_ji_conj
   
   if (.not. allocated(this%hk_total)) return
   
   n = size(this%hk_total, 1)
   max_diff = 0.0_rp
   
   do i = 1, n
      do j = 1, n
         h_ij = this%hk_total(i, j, ik)
         h_ji_conj = conjg(this%hk_total(j, i, ik))
         diff = abs(h_ij - h_ji_conj)
         max_diff = max(max_diff, diff)
      end do
   end do
   
   if (max_diff > 1.0e-8_rp) then
      call g_logger%warning('Hermiticity check: max violation = ' // &
                           trim(real2str(max_diff, '(ES12.4)')) // ' at k-point ' // &
                           trim(int2str(ik)), __FILE__, __LINE__)
   else
      call g_logger%info('Hermiticity check: H(k) is Hermitian (max diff = ' // &
                        trim(real2str(max_diff, '(ES12.4)')) // ')', __FILE__, __LINE__)
   end if
end subroutine check_hamiltonian_hermiticity

module subroutine print_hamiltonian_structure(this, ik)
   class(reciprocal), intent(in) :: this
   integer, intent(in) :: ik
   
   integer :: i, j, n_sites, block_size
   integer :: iblock, jblock, i_start, j_start
   real(rp) :: block_norm
   character(len=500) :: msg
   
   if (.not. allocated(this%hk_total)) return
   
   n_sites = this%lattice%nrec
   block_size = 18
   
   call g_logger%info('=== H(k) Block Structure at k-point ' // trim(int2str(ik)) // ' ===', &
                     __FILE__, __LINE__)
   
   ! Print block-wise norms to see coupling structure
   do iblock = 1, n_sites
      msg = 'Site ' // trim(int2str(iblock)) // ' couples to: '
      do jblock = 1, n_sites
         i_start = (iblock-1)*block_size + 1
         j_start = (jblock-1)*block_size + 1
         
         ! Calculate Frobenius norm of this block
         block_norm = 0.0_rp
         do i = 0, block_size-1
            do j = 0, block_size-1
               block_norm = block_norm + abs(this%hk_total(i_start+i, j_start+j, ik))**2
            end do
         end do
         block_norm = sqrt(block_norm)
         
         if (block_norm > 1.0e-6_rp) then
            msg = trim(msg) // ' site' // trim(int2str(jblock)) // &
                  '(' // trim(real2str(block_norm, '(F8.4)')) // ')'
         end if
      end do
      call g_logger%info(trim(msg), __FILE__, __LINE__)
   end do
   
   call g_logger%info('=== End H(k) Structure ===', __FILE__, __LINE__)
   
end subroutine print_hamiltonian_structure

   module function get_basis_type_from_size(this, ntype) result(basis_type)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: ntype
      character(len=10) :: basis_type

      if (.not. allocated(this%basis_size) .or. ntype > size(this%basis_size)) then
         basis_type = 'unknown'
         return
      end if

      select case (this%basis_size(ntype))
      case (4)
         basis_type = 'sp'
      case (9)
         basis_type = 'spd'
      case (16)
         basis_type = 'spdf'
      case default
         basis_type = 'custom'
      end select
   end function get_basis_type_from_size

end submodule reciprocal_fourier
