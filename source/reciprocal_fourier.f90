submodule (reciprocal_mod) reciprocal_fourier
   use lmto_magnetic_tangent_mod, only: lmto_hhmag_to_spinor, lmto_endpoint_tangent_record
   use lmto_pair_potential_mod, only: lmto_circular_pair_potential, lmto_circular_pair_potential_from_reverse, lmto_endpoint_phases, &
      lmto_pair_transition_metadata, lmto_transition_metadata
   implicit none

contains

   module subroutine reciprocal_workspace_ensure_capacity(this, nmat, tile_length, generalized, operator_generation, nnmax, ntype)
      class(reciprocal_workspace), intent(inout) :: this
      integer, intent(in) :: nmat, tile_length, operator_generation, nnmax, ntype
      logical, intent(in) :: generalized
      complex(rp) :: query(1)
      integer :: info, requested
      logical :: reusable

      if (nmat <= 0 .or. tile_length <= 0 .or. nnmax <= 0 .or. ntype <= 0) then
         call g_logger%fatal('reciprocal_workspace%ensure_capacity: invalid workspace dimensions.', __FILE__, __LINE__)
      end if
      requested = max(tile_length, 1)
      reusable = .false.
      if (allocated(this%phase)) reusable = this%nmat == nmat .and. this%tile_capacity >= requested .and. &
         (this%generalized .eqv. generalized) .and. size(this%phase, 1) == nnmax .and. size(this%phase, 2) == ntype
      if (reusable) then
         this%active_tile_length = tile_length
         this%generalized = generalized
         this%cached_operator_generation = operator_generation
         this%capacity_reuses = this%capacity_reuses + 1
         return
      end if

      call this%clear()
      this%nmat = nmat
      this%tile_capacity = requested
      this%active_tile_length = tile_length
      this%generalized = generalized
      this%cached_operator_generation = operator_generation
      this%storage_allocations = this%storage_allocations + 1
      allocate(this%h(nmat,nmat,requested), this%s(nmat,nmat,requested), this%eeo(nmat,nmat,requested), &
               this%hoh(nmat,nmat,requested), this%hcc(nmat,nmat,requested), this%phase(nnmax,ntype,requested), this%points(3,requested), &
               this%eigenvalue(nmat,requested), this%eigenvector(nmat,nmat,requested), this%info(requested), &
               this%lapack_rwork(max(1,3*nmat-2)))
      this%h = cmplx(0.0_rp,0.0_rp,rp); this%s = cmplx(0.0_rp,0.0_rp,rp)
      this%eeo = cmplx(0.0_rp,0.0_rp,rp); this%hoh = cmplx(0.0_rp,0.0_rp,rp)
      this%hcc = cmplx(0.0_rp,0.0_rp,rp); this%phase = cmplx(0.0_rp,0.0_rp,rp)
      if (generalized) then
         call zhegv(1, 'V', 'U', nmat, this%h(:,:,1), nmat, this%s(:,:,1), nmat, this%eigenvalue(:,1), &
                    query, -1, this%lapack_rwork, info)
      else
         call zheev('V', 'U', nmat, this%h(:,:,1), nmat, this%eigenvalue(:,1), query, -1, this%lapack_rwork, info)
      end if
      if (info /= 0) call g_logger%fatal('reciprocal_workspace%ensure_capacity: LAPACK workspace query failed.', __FILE__, __LINE__)
      this%lapack_workspace_queries = this%lapack_workspace_queries + 1
      this%lwork = max(1, int(real(query(1), rp)))
      allocate(this%lapack_work(this%lwork))
   end subroutine reciprocal_workspace_ensure_capacity

   module subroutine reciprocal_workspace_clear(this)
      class(reciprocal_workspace), intent(inout) :: this
      if (allocated(this%h)) deallocate(this%h)
      if (allocated(this%s)) deallocate(this%s)
      if (allocated(this%eeo)) deallocate(this%eeo)
      if (allocated(this%hoh)) deallocate(this%hoh)
      if (allocated(this%hcc)) deallocate(this%hcc)
      if (allocated(this%phase)) deallocate(this%phase)
      if (allocated(this%points)) deallocate(this%points)
      if (allocated(this%eigenvalue)) deallocate(this%eigenvalue)
      if (allocated(this%eigenvector)) deallocate(this%eigenvector)
      if (allocated(this%lapack_work)) deallocate(this%lapack_work)
      if (allocated(this%lapack_rwork)) deallocate(this%lapack_rwork)
      if (allocated(this%info)) deallocate(this%info)
      this%nmat = 0; this%tile_capacity = 0; this%active_tile_length = 0
      this%lwork = 0; this%cached_operator_generation = -1; this%generalized = .false.
   end subroutine reciprocal_workspace_clear

   module subroutine reciprocal_workspace_restore_to_default(this)
      class(reciprocal_workspace), intent(inout) :: this
      call this%clear()
      this%storage_allocations = 0
      this%lapack_workspace_queries = 0
      this%capacity_reuses = 0
   end subroutine reciprocal_workspace_restore_to_default

   module subroutine reciprocal_workspace_destructor(this)
      type(reciprocal_workspace) :: this
      call this%clear()
   end subroutine reciprocal_workspace_destructor

   module subroutine make_reciprocal_assembler(this, assembler)
      class(reciprocal), intent(inout) :: this
      type(reciprocal_assembler), intent(out) :: assembler
      call this%build_neighbor_vectors()
      assembler%hamiltonian => this%hamiltonian
      assembler%lattice => this%lattice
      assembler%control => this%control
      allocate(assembler%ham_vec_type_direct, source=this%ham_vec_type_direct)
      select case (trim(this%kspace_ham_order))
      case ('second'); assembler%use_second_order = .true.
      case ('first');  assembler%use_second_order = .false.
      case default;    assembler%use_second_order = this%hamiltonian%hoh
      end select
   end subroutine make_reciprocal_assembler

   ! Assemble a whole tile.  The loop is serial by design: a workspace is
   ! exclusive to its caller, avoiding nested OpenMP/BLAS oversubscription.
   ! Parallel callers should use one workspace per worker.
   module subroutine reciprocal_assembler_assemble_batch(this, k_points, workspace)
      class(reciprocal_assembler), intent(in) :: this
      real(rp), intent(in) :: k_points(:, :)
      class(reciprocal_workspace), intent(inout) :: workspace
      integer :: ik, isite, jsite, ntype_i, ineigh, ia, ja, nr, i_start, i_end, j_start, j_end, nmat
      real(rp) :: kfold(3), angle
      logical :: second_order, add_ccor

      if (size(k_points,1) /= 3 .or. .not. associated(this%hamiltonian) .or. .not. associated(this%lattice) .or. &
          .not. allocated(this%ham_vec_type_direct)) then
         call g_logger%fatal('reciprocal_assembler%assemble_batch: incomplete assembly state.', __FILE__, __LINE__)
      end if
      nmat = nb*this%lattice%nrec
      if (workspace%nmat /= nmat .or. workspace%tile_capacity < size(k_points,2)) then
         call workspace%ensure_capacity(nmat, size(k_points,2), workspace%generalized, this%hamiltonian%operator_generation, &
                                        this%lattice%nn_max, this%lattice%ntype)
      else
         workspace%active_tile_length = size(k_points,2)
         workspace%cached_operator_generation = this%hamiltonian%operator_generation
      end if
      second_order = this%use_second_order .and. this%hamiltonian%hoh .and. allocated(this%hamiltonian%eeo)
      add_ccor = this%hamiltonian%ccor_2c
      do ik = 1, size(k_points,2)
         kfold = k_points(:,ik) - floor(k_points(:,ik) + 0.5_rp)
         workspace%phase(:,:,ik) = cmplx(0.0_rp,0.0_rp,rp)
         do ntype_i = 1, this%lattice%ntype
            ia = this%lattice%atlist(ntype_i); nr = this%lattice%nn(ia,1)
            do ineigh = 1, nr
               if (ineigh == 1) then
                  workspace%phase(ineigh,ntype_i,ik) = cmplx(1.0_rp,0.0_rp,rp)
               else
                  angle = 2.0_rp*pi*dot_product(kfold, this%ham_vec_type_direct(:,ineigh,ntype_i))
                  workspace%phase(ineigh,ntype_i,ik) = cmplx(cos(angle),sin(angle),rp)
               end if
            end do
         end do
         workspace%h(:,:,ik) = cmplx(0.0_rp,0.0_rp,rp)
         workspace%eeo(:,:,ik) = cmplx(0.0_rp,0.0_rp,rp)
         workspace%hcc(:,:,ik) = cmplx(0.0_rp,0.0_rp,rp)
         do isite = 1, this%lattice%nrec
            ntype_i = this%lattice%ib(isite); ia = this%lattice%atlist(ntype_i); nr = this%lattice%nn(ia,1)
            i_start = (isite-1)*nb+1; i_end = isite*nb
            do ineigh = 1, nr
               if (ineigh == 1) then
                  jsite = isite
               else
                  ja = this%lattice%nn(ia,ineigh)
                  if (ja < 1 .or. ja > this%lattice%kk) cycle
                  jsite = this%lattice%iz(ja)
                  if (jsite < 1 .or. jsite > this%lattice%nrec) cycle
               end if
               j_start = (jsite-1)*nb+1; j_end = jsite*nb
               workspace%h(i_start:i_end,j_start:j_end,ik) = workspace%h(i_start:i_end,j_start:j_end,ik) + &
                  this%hamiltonian%ee(:,:,ineigh,ntype_i)*workspace%phase(ineigh,ntype_i,ik)
               if (second_order) workspace%eeo(i_start:i_end,j_start:j_end,ik) = workspace%eeo(i_start:i_end,j_start:j_end,ik) + &
                  this%hamiltonian%eeo(:,:,ineigh,ntype_i)*workspace%phase(ineigh,ntype_i,ik)
               if (add_ccor) workspace%hcc(i_start:i_end,j_start:j_end,ik) = workspace%hcc(i_start:i_end,j_start:j_end,ik) + &
                  this%hamiltonian%eecc(:,:,ineigh,ntype_i)*workspace%phase(ineigh,ntype_i,ik)
            end do
         end do
         if (second_order) then
            call zgemm('N','N',nmat,nmat,nmat,cone,workspace%eeo(:,:,ik),nmat,workspace%h(:,:,ik),nmat,czero,workspace%hoh(:,:,ik),nmat)
            workspace%h(:,:,ik) = workspace%h(:,:,ik) - workspace%hoh(:,:,ik)
            do isite = 1, this%lattice%nrec
               ntype_i = this%lattice%ib(isite); i_start=(isite-1)*nb+1; i_end=isite*nb
               workspace%h(i_start:i_end,i_start:i_end,ik) = workspace%h(i_start:i_end,i_start:i_end,ik) + &
                  this%hamiltonian%enim(:,:,ntype_i) + this%hamiltonian%lsham(:,:,ntype_i)
            end do
         end if
         if (add_ccor) workspace%h(:,:,ik) = workspace%h(:,:,ik) + workspace%hcc(:,:,ik)
      end do
   end subroutine reciprocal_assembler_assemble_batch

   module subroutine reciprocal_assembler_assemble_overlap_batch(this, k_points, workspace)
      class(reciprocal_assembler), intent(in) :: this
      real(rp), intent(in) :: k_points(:, :)
      class(reciprocal_workspace), intent(inout) :: workspace
      integer :: ik, isite, jsite, ntype_i, ineigh, ia, ja, nr, i_start, i_end, j_start, j_end, iorb
      if (workspace%active_tile_length /= size(k_points,2)) call this%assemble_batch(k_points, workspace)
      do ik = 1, size(k_points,2)
         workspace%s(:,:,ik) = cmplx(0.0_rp,0.0_rp,rp)
         do isite = 1, this%lattice%nrec
            ntype_i=this%lattice%ib(isite); ia=this%lattice%atlist(ntype_i); nr=this%lattice%nn(ia,1)
            i_start=(isite-1)*nb+1; i_end=isite*nb
            do ineigh=1,nr
               if (ineigh == 1) then; jsite=isite
               else
                  ja=this%lattice%nn(ia,ineigh); if (ja < 1 .or. ja > this%lattice%kk) cycle
                  jsite=this%lattice%iz(ja); if (jsite < 1 .or. jsite > this%lattice%nrec) cycle
               end if
               j_start=(jsite-1)*nb+1; j_end=jsite*nb
               workspace%s(i_start:i_end,j_start:j_end,ik) = workspace%s(i_start:i_end,j_start:j_end,ik) + &
                  this%hamiltonian%eeo(:,:,ineigh,ntype_i)*workspace%phase(ineigh,ntype_i,ik)
               if (ineigh == 1 .and. jsite == isite) then
                  do iorb=1,nb; workspace%s(i_start+iorb-1,i_start+iorb-1,ik) = workspace%s(i_start+iorb-1,i_start+iorb-1,ik) + 1.0_rp; end do
               end if
            end do
         end do
      end do
   end subroutine reciprocal_assembler_assemble_overlap_batch

   module subroutine reciprocal_assembler_assemble_one(this, k_point, hk_result, workspace)
      class(reciprocal_assembler), intent(in) :: this
      real(rp), intent(in) :: k_point(3)
      complex(rp), intent(out) :: hk_result(:, :)
      class(reciprocal_workspace), intent(inout) :: workspace
      real(rp) :: point(3,1)
      point(:,1) = k_point
      call this%assemble_batch(point, workspace)
      hk_result = workspace%h(:,:,1)
   end subroutine reciprocal_assembler_assemble_one

   !> @brief Calculate exp(i k.R) factors for each neighbor/type.
   !> @param[in] this Reciprocal object containing neighbor-vector tables.
   !> @param[in] k_vec k-point vector in reciprocal coordinates.
   !> @param[out] structure_factors Phase factors indexed by neighbor and atom type.
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

   !> @brief Fourier transform first-order real-space Hamiltonian blocks.
   !> @details Forms H(k)=sum_R h(R) exp(i k.R), preserving the historical
   !>          first-order path without onsite e_nu or spin-orbit terms.
   !> @param[in] this Reciprocal object containing real-space Hamiltonian data.
   !> @param[in] k_vec k-point vector.
   !> @param[out] hk_result Packed k-space Hamiltonian matrix.
   module subroutine fourier_transform_hamiltonian(this, k_vec, hk_result)
      class(reciprocal), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: hk_result  ! (n_orb*n_sites, n_orb*n_sites)
      type(reciprocal_assembler) :: assembler

      ! First-order k-space Hamiltonian: H(k) = h(k) = Sum_R ee(R) exp(i*k·R).
      ! NOTE: This deliberately reproduces the historical first-order behaviour.
      !       The on-site E_nu (enim) and spin-orbit (lsham) terms are NOT added
      !       here; they are only included in the second-order path
      !       (fourier_transform_hamiltonian_second_order). See kspace_ham_order.

      call this%make_reciprocal_assembler(assembler)
      assembler%use_second_order = .false.
      call this%workspace%ensure_capacity(size(hk_result,1), 1, .false., this%hamiltonian%operator_generation, &
                                          this%lattice%nn_max, this%lattice%ntype)
      call assembler%assemble_one(k_vec, hk_result, this%workspace)
   end subroutine fourier_transform_hamiltonian

   !> Build the finite-q LMTO pair potential in exactly the `ham_only` coefficient
   !> basis used by the reciprocal eigensolver.  The endpoint tangent is asked
   !> for before hxc/ee collapse; its directed-bond phase is the normal
   !> H(k)=sum_R H(R) exp(+i k.R) phase: left endpoints carry k and right
   !> endpoints carry k+q. The optional q=0 default preserves WR-02 callers.
   module subroutine build_lmto_pair_potential_at_kpoint(this, response_site, k_point, signed_moment, qminus, qplus, supported, reason, q_point, metadata)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: response_site
      real(rp), intent(in) :: k_point(3), signed_moment
      complex(rp), intent(out) :: qminus(:, :), qplus(:, :)
      logical, intent(out) :: supported
      character(len=*), intent(out), optional :: reason
      real(rp), intent(in), optional :: q_point(3)
      type(lmto_pair_transition_metadata), intent(out), optional :: metadata
      integer :: nmat, isite, jsite, ntype_i, ia, ja, it, jt, ineigh, nr, ni, idir
      integer :: i_start, i_end, j_start, j_end
      real(rp) :: vet(3)
      real(rp) :: ham_vec(3, this%lattice%nn_max)
      real(rp) :: hhh(norb, norb)
      complex(rp) :: left(norb, norb, 4, 3), right(norb, norb, 4, 3)
      complex(rp) :: left_cart(norb, norb, 4), right_cart(norb, norb, 4)
      complex(rp) :: left_spinor(2*norb, 2*norb), right_spinor(2*norb, 2*norb)
      complex(rp), allocatable :: dh_dx(:, :), dh_dy(:, :), dh_dx_reverse(:, :), dh_dy_reverse(:, :)
      complex(rp) :: left_phase, right_phase
      type(lmto_endpoint_tangent_record) :: record
      character(len=160) :: local_reason
      real(rp) :: q_direct(3)

      nmat = 2*norb*this%lattice%nrec
      q_direct = 0.0_rp; if (present(q_point)) q_direct = q_point
      if (present(metadata)) metadata = lmto_transition_metadata(k_point, q_direct)
      qminus = cmplx(0.0_rp, 0.0_rp, rp); qplus = cmplx(0.0_rp, 0.0_rp, rp)
      supported = .false.
      if (size(qminus,1) /= nmat .or. size(qminus,2) /= nmat .or. any(shape(qplus) /= shape(qminus))) then
         if (present(reason)) reason = 'pair-potential matrix shape does not match site-major ham_only basis'
         return
      end if
      if (response_site < 1 .or. response_site > this%lattice%nrec) then
         if (present(reason)) reason = 'invalid response-site identity'
         return
      end if
      if (trim(this%reciprocal_mode) /= 'ham_only') then
         if (present(reason)) reason = 'LMTO pair potential requires reciprocal_mode=ham_only'
         return
      end if

      call this%build_neighbor_vectors()
      allocate(dh_dx(nmat,nmat), dh_dy(nmat,nmat), dh_dx_reverse(nmat,nmat), dh_dy_reverse(nmat,nmat))
      dh_dx = cmplx(0.0_rp, 0.0_rp, rp); dh_dy = cmplx(0.0_rp, 0.0_rp, rp)
      dh_dx_reverse = cmplx(0.0_rp, 0.0_rp, rp); dh_dy_reverse = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, this%lattice%nrec
         ntype_i = this%lattice%ib(isite)
         ia = this%lattice%atlist(ntype_i); it = this%lattice%iz(ia)
         nr = this%lattice%nn(ia, 1)
         ham_vec = 0.0_rp
         ham_vec(:,1:nr) = this%ham_vec_type(:,1:nr,ntype_i)
         i_start = (isite-1)*2*norb + 1; i_end = isite*2*norb
         do ineigh = 1, nr
            if (ineigh == 1) then
               ja = ia; jsite = isite
            else
               ja = this%lattice%nn(ia, ineigh)
               if (ja < 1 .or. ja > this%lattice%kk) cycle
               jsite = this%lattice%iz(ja)
               if (jsite < 1 .or. jsite > this%lattice%nrec) cycle
            end if
            jt = this%lattice%iz(ja); vet = this%ham_vec_type(:,ineigh,ntype_i)
            call this%hamiltonian%hmfind(vet, nr, hhh, ineigh, ia, ineigh, ni, ham_vec(:,1:nr))
            if (ni == 0) cycle
            call this%hamiltonian%ham0m_nc_endpoint_tangents(ia, ja, it, jt, ineigh, vet, hhh, left, right, &
                                                               record, supported, local_reason)
            if (.not. supported) then
               if (present(reason)) reason = trim(local_reason)
               deallocate(dh_dx, dh_dy, dh_dx_reverse, dh_dy_reverse)
               return
            end if
            j_start = (jsite-1)*2*norb + 1; j_end = jsite*2*norb
            call lmto_endpoint_phases(k_point, q_direct, this%ham_vec_type_direct(:,ineigh,ntype_i), left_phase, right_phase)
            do idir = 1, 2
               if (isite == response_site) then
                  left_cart = left(:,:,:,idir)
                  call hcpx(left_cart(:,:,1), 'cart2sph'); call hcpx(left_cart(:,:,2), 'cart2sph')
                  call hcpx(left_cart(:,:,3), 'cart2sph'); call hcpx(left_cart(:,:,4), 'cart2sph')
                  call lmto_hhmag_to_spinor(left_cart, left_spinor)
                  if (idir == 1) dh_dx(i_start:i_end,j_start:j_end) = dh_dx(i_start:i_end,j_start:j_end) + left_spinor*left_phase
                  if (idir == 2) dh_dy(i_start:i_end,j_start:j_end) = dh_dy(i_start:i_end,j_start:j_end) + left_spinor*left_phase
                  ! Reverse transition K -> k: use the conjugate reverse-bond
                  ! contribution and its conjugate Bloch phase.  Keeping this
                  ! accumulator separate makes Q+ an independently assembled
                  ! reverse-q operator instead of a post-hoc Q- adjoint.
                  if (idir == 1) dh_dx_reverse(j_start:j_end,i_start:i_end) = dh_dx_reverse(j_start:j_end,i_start:i_end) + &
                     transpose(conjg(left_spinor*left_phase))
                  if (idir == 2) dh_dy_reverse(j_start:j_end,i_start:i_end) = dh_dy_reverse(j_start:j_end,i_start:i_end) + &
                     transpose(conjg(left_spinor*left_phase))
               end if
               if (jsite == response_site) then
                  right_cart = right(:,:,:,idir)
                  call hcpx(right_cart(:,:,1), 'cart2sph'); call hcpx(right_cart(:,:,2), 'cart2sph')
                  call hcpx(right_cart(:,:,3), 'cart2sph'); call hcpx(right_cart(:,:,4), 'cart2sph')
                  call lmto_hhmag_to_spinor(right_cart, right_spinor)
                  if (idir == 1) dh_dx(i_start:i_end,j_start:j_end) = dh_dx(i_start:i_end,j_start:j_end) + right_spinor*right_phase
                  if (idir == 2) dh_dy(i_start:i_end,j_start:j_end) = dh_dy(i_start:i_end,j_start:j_end) + right_spinor*right_phase
                  if (idir == 1) dh_dx_reverse(j_start:j_end,i_start:i_end) = dh_dx_reverse(j_start:j_end,i_start:i_end) + &
                     transpose(conjg(right_spinor*right_phase))
                  if (idir == 2) dh_dy_reverse(j_start:j_end,i_start:i_end) = dh_dy_reverse(j_start:j_end,i_start:i_end) + &
                     transpose(conjg(right_spinor*right_phase))
               end if
            end do
         end do
      end do
      call lmto_circular_pair_potential(dh_dx, dh_dy, signed_moment, qminus, qplus, supported, local_reason)
      if (supported) call lmto_circular_pair_potential_from_reverse(dh_dx_reverse, dh_dy_reverse, signed_moment, qplus, supported, local_reason)
      if (present(reason)) reason = trim(local_reason)
      deallocate(dh_dx, dh_dy, dh_dx_reverse, dh_dy_reverse)
   end subroutine build_lmto_pair_potential_at_kpoint

   !> @brief Fourier transform an arbitrary neighbor/type block array.
   !> @details Applies the reciprocal neighbor map to a (orbital, orbital,
   !>          neighbor, type) array and packs the result into site-major matrix form.
   !> @param[in] this Reciprocal object containing neighbor-vector tables.
   !> @param[in] array4d Real-space block array indexed by orbital, neighbor, and type.
   !> @param[in] k_vec k-point vector.
   !> @param[out] mk_result Packed k-space matrix.
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

   !> @brief Fourier transform the second-order ASA k-space Hamiltonian.
   !> @details Assembles onsite e_nu, first-order hopping, HOH, optional CCOR,
   !>          and spin-orbit contributions according to kspace_ham_order.
   !> @param[in] this Reciprocal object containing Hamiltonian and correction blocks.
   !> @param[in] k_vec k-point vector.
   !> @param[out] hk_result Packed second-order k-space Hamiltonian matrix.
   module subroutine fourier_transform_hamiltonian_second_order(this, k_vec, hk_result)
      class(reciprocal), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: hk_result  ! (n_orb*n_sites, n_orb*n_sites)
      type(reciprocal_assembler) :: assembler
      call this%make_reciprocal_assembler(assembler)
      assembler%use_second_order = .true.
      call this%workspace%ensure_capacity(size(hk_result,1), 1, .false., this%hamiltonian%operator_generation, &
                                          this%lattice%nn_max, this%lattice%ntype)
      call assembler%assemble_one(k_vec, hk_result, this%workspace)
   end subroutine fourier_transform_hamiltonian_second_order

   !> @brief Fourier transform overlap blocks into S(k).
   !> @details Builds the reciprocal-space overlap matrix used by generalized
   !>          eigenproblem modes and overlap diagnostics.
   !> @param[in] this Reciprocal object containing overlap and neighbor data.
   !> @param[in] k_vec k-point vector.
   !> @param[out] sk_result Packed k-space overlap matrix.
   module subroutine fourier_transform_overlap(this, k_vec, sk_result)
      class(reciprocal), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: sk_result
      type(reciprocal_assembler) :: assembler
      call this%make_reciprocal_assembler(assembler)
      call this%workspace%ensure_capacity(size(sk_result,1), 1, .true., this%hamiltonian%operator_generation, &
                                          this%lattice%nn_max, this%lattice%ntype)
      this%workspace%points(:,1) = k_vec
      call assembler%assemble_overlap_batch(this%workspace%points(:,1:1), this%workspace)
      sk_result = this%workspace%s(:,:,1)
   end subroutine fourier_transform_overlap

   !> Fold an arbitrary fractional reciprocal point into the first reciprocal
   !> cell.  The half-open convention [-1/2,1/2) gives one canonical
   !> representative while retaining the exact exp(i 2*pi*k.R) phase for every
   !> lattice translation R.
   module subroutine fold_kpoint(this, k_point, k_folded)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: k_point(3)
      real(rp), intent(out) :: k_folded(3)

      k_folded = k_point - floor(k_point + 0.5_rp)
   end subroutine fold_kpoint

   !> Build the normal-state H(k) at one arbitrary k point without consulting
   !> or changing the active reciprocal mesh.  This is the single assembly
   !> route shared by mesh bands/DOS and response-style arbitrary-k callers.
   module subroutine build_hamiltonian_at_kpoint(this, k_point, hk_result)
      class(reciprocal), intent(inout) :: this
      real(rp), intent(in) :: k_point(3)
      complex(rp), intent(out) :: hk_result(:, :)

      ! The neighbor table is a geometry cache, not a k-mesh cache.  It is the
      ! only object-local setup this service needs and never alters k_points,
      ! hk_bulk, eigenvalues, or DOS arrays.
      call this%build_neighbor_vectors()
      call assemble_hamiltonian_at_kpoint(this, k_point, hk_result)
   end subroutine build_hamiltonian_at_kpoint

   !> Return caller-owned eigenpairs at arbitrary fractional reciprocal points.
   !> Points supplied as k+q are folded only by reciprocal lattice vectors; no
   !> nearest-mesh mapping is ever performed.  Exact duplicate folded points
   !> within this call share one eigensolve, which is useful for repeated q
   !> work while keeping the service MPI-local and side-effect free for the
   !> standard bands/DOS mesh.
   module subroutine calculate_eigenpairs_at_kpoints(this, k_points, eigenvalues, eigenvectors, folded_k_points)
      class(reciprocal), intent(inout) :: this
      real(rp), intent(in) :: k_points(:, :)
      real(rp), allocatable, intent(out) :: eigenvalues(:, :)
      complex(rp), allocatable, intent(out) :: eigenvectors(:, :, :)
      real(rp), allocatable, intent(out), optional :: folded_k_points(:, :)

      integer :: nk, nmat, ik, iu, nunique, tile_first, tile_last, tile_length, slot
      integer, allocatable :: representative(:), unique_of_k(:)
      real(rp), allocatable :: folded(:, :)
      logical :: use_generalized
      type(reciprocal_execution_request) :: request
      type(reciprocal_execution_result) :: result

      if (size(k_points, 1) /= 3) then
         call g_logger%error('calculate_eigenpairs_at_kpoints: k_points must have shape (3,nk).', __FILE__, __LINE__)
         return
      end if
      if (.not. associated(this%lattice) .or. .not. associated(this%hamiltonian)) then
         call g_logger%error('calculate_eigenpairs_at_kpoints: reciprocal object is not initialized.', __FILE__, __LINE__)
         return
      end if
      call this%validate_nonzero_q_gbt('reciprocal%calculate_eigenpairs_at_kpoints')

      nk = size(k_points, 2)
      nmat = nb*this%lattice%nrec
      if (nmat <= 0) then
         call g_logger%error('calculate_eigenpairs_at_kpoints: invalid reciprocal basis size.', __FILE__, __LINE__)
         return
      end if

      allocate(eigenvalues(nmat, nk), eigenvectors(nmat, nmat, nk))
      allocate(folded(3, nk), representative(nk), unique_of_k(nk))
      do ik = 1, nk
         call this%fold_kpoint(k_points(:, ik), folded(:, ik))
      end do
      if (present(folded_k_points)) then
         allocate(folded_k_points(3, nk))
         folded_k_points = folded
      end if

      ! Only bit-identical canonical points are cached.  A tolerance-based
      ! merge would turn two distinct off-mesh Hamiltonians into an
      ! approximation, which is not acceptable for the response definition.
      nunique = 0
      do ik = 1, nk
         unique_of_k(ik) = 0
         do iu = 1, nunique
            if (all(folded(:, ik) == folded(:, representative(iu)))) then
               unique_of_k(ik) = iu
               exit
            end if
         end do
         if (unique_of_k(ik) == 0) then
            nunique = nunique + 1
            representative(nunique) = ik
            unique_of_k(ik) = nunique
         end if
      end do

      use_generalized = trim(this%reciprocal_mode) == 'generalized_overlap_proxy'
      ! Fold/deduplicate once, then process unique points in persistent CPU
      ! tiles.  The deferred call is once per tile, never per k-point or band.
      call this%make_execution_backend()
      do tile_first = 1, nunique, max(1, this%reciprocal_tile_size)
         tile_last = min(nunique, tile_first + max(1, this%reciprocal_tile_size) - 1)
         tile_length = tile_last - tile_first + 1
         request%assemble_hamiltonian = .true.
         request%assemble_overlap = use_generalized
         request%solve_eigensystem = .true.
         request%generalized = use_generalized
         request%request_eigenvectors = .true.
         request%request_assembled_hamiltonian = .false.
         request%request_assembled_overlap = .false.
         request%operator_generation = this%hamiltonian%operator_generation
         if (allocated(request%k_points)) deallocate(request%k_points)
         allocate(request%k_points(3,tile_length))
         do slot = 1, tile_length
            ik = representative(tile_first + slot - 1)
            request%k_points(:,slot) = folded(:,ik)
         end do
         call this%execution_backend%execute_batch(request, result)
         do slot = 1, tile_length
            ik = representative(tile_first + slot - 1)
            eigenvalues(:,ik) = result%eigenvalues(:,slot)
            eigenvectors(:,:,ik) = result%eigenvectors(:,:,slot)
         end do
      end do
      ! RF-04 exposed this workspace for allocation/cache diagnostics.  Keep
      ! that view coherent while the backend owns the active execution tile.
      select type (backend => this%execution_backend)
      type is (lapack_reciprocal_backend)
         this%workspace%tile_capacity = backend%workspace%tile_capacity
         this%workspace%cached_operator_generation = backend%workspace%cached_operator_generation
         this%workspace%storage_allocations = backend%workspace%storage_allocations
         this%workspace%lapack_workspace_queries = backend%workspace%lapack_workspace_queries
         this%workspace%capacity_reuses = backend%workspace%capacity_reuses
      end select
      do ik = 1, nk
         iu = unique_of_k(ik)
         if (representative(iu) /= ik) then
            eigenvalues(:, ik) = eigenvalues(:, representative(iu))
            eigenvectors(:, :, ik) = eigenvectors(:, :, representative(iu))
         end if
      end do

      deallocate(folded, representative, unique_of_k)
   end subroutine calculate_eigenpairs_at_kpoints

   !> Shared single-point Hamiltonian assembly.  Keeping this underneath both
   !> build_kspace_hamiltonian and the arbitrary-k API ensures that SOC,
   !> non-collinearity, Hubbard/CCOR blocks, and second-order terms follow the
   !> normal reciprocal calculation automatically.
   subroutine assemble_hamiltonian_at_kpoint(this, k_point, hk_result)
      class(reciprocal), intent(inout) :: this
      real(rp), intent(in) :: k_point(3)
      complex(rp), intent(out) :: hk_result(:, :)
      real(rp) :: k_folded(3)
      logical :: use_second_order

      call this%fold_kpoint(k_point, k_folded)
      select case (trim(this%kspace_ham_order))
      case ('second')
         use_second_order = .true.
      case ('first')
         use_second_order = .false.
      case default
         use_second_order = this%hamiltonian%hoh
      end select
      if (use_second_order .and. (.not. this%hamiltonian%hoh .or. .not. allocated(this%hamiltonian%eeo))) then
         use_second_order = .false.
      end if

      if (use_second_order) then
         call this%fourier_transform_hamiltonian_second_order(k_folded, hk_result)
      else
         call this%fourier_transform_hamiltonian(k_folded, hk_result)
      end if
   end subroutine assemble_hamiltonian_at_kpoint

   !> Diagonalize one exclusive workspace slice.  The caller prepared the
   !> workspace, including the mode-specific LAPACK work size, before entering
   !> its k-point tile loop.
   subroutine diagonalize_workspace_slot(this, workspace, slot, use_generalized)
      class(reciprocal), intent(in) :: this
      class(reciprocal_workspace), intent(inout) :: workspace
      integer, intent(in) :: slot
      logical, intent(in) :: use_generalized
      real(rp) :: max_herm, matrix_scale

      max_herm = maxval(abs(workspace%h(:,:,slot) - transpose(conjg(workspace%h(:,:,slot)))))
      matrix_scale = max(1.0_rp, maxval(abs(workspace%h(:,:,slot))))
      if (max_herm > 1.0e-10_rp*matrix_scale) then
         call g_logger%fatal('reciprocal workspace H(k) is non-Hermitian before eigensolution.', __FILE__, __LINE__)
      end if
      workspace%eigenvector(:,:,slot) = workspace%h(:,:,slot)
      if (use_generalized) then
         call this%check_overlap_properties(slot, workspace%s(:,:,slot))
         call zhegv(1, 'V', 'U', workspace%nmat, workspace%eigenvector(:,:,slot), workspace%nmat, workspace%s(:,:,slot), &
                    workspace%nmat, workspace%eigenvalue(:,slot), workspace%lapack_work, workspace%lwork, workspace%lapack_rwork, workspace%info(slot))
      else
         call zheev('V', 'U', workspace%nmat, workspace%eigenvector(:,:,slot), workspace%nmat, workspace%eigenvalue(:,slot), &
                    workspace%lapack_work, workspace%lwork, workspace%lapack_rwork, workspace%info(slot))
      end if
      if (workspace%info(slot) /= 0) call g_logger%fatal('reciprocal workspace LAPACK eigensolver failed.', __FILE__, __LINE__)
   end subroutine diagonalize_workspace_slot

   !> Solve one normal-state reciprocal eigenproblem using the same LAPACK
   !> convention as diagonalize_hamiltonian.  The optional proxy overlap is
   !> formed locally so the standard sk_overlap mesh cache remains untouched.
   subroutine diagonalize_single_kpoint(this, k_point, hk, use_generalized, eigenvals, eigenvecs)
      class(reciprocal), intent(inout) :: this
      real(rp), intent(in) :: k_point(3)
      complex(rp), intent(in) :: hk(:, :)
      logical, intent(in) :: use_generalized
      real(rp), intent(out) :: eigenvals(:)
      complex(rp), intent(out) :: eigenvecs(:, :)

      integer :: nmat, lwork, info
      complex(rp), allocatable :: h_copy(:, :), s_copy(:, :), work(:)
      real(rp), allocatable :: rwork(:)
      complex(rp) :: work_query(1)
      real(rp) :: max_herm, matrix_scale

      nmat = size(hk, 1)
      if (size(hk, 2) /= nmat .or. size(eigenvals) /= nmat .or. &
          size(eigenvecs, 1) /= nmat .or. size(eigenvecs, 2) /= nmat) then
         call g_logger%fatal('diagonalize_single_kpoint: inconsistent eigensystem dimensions.', __FILE__, __LINE__)
      end if
      max_herm = maxval(abs(hk - transpose(conjg(hk))))
      matrix_scale = max(1.0_rp, maxval(abs(hk)))
      if (max_herm > 1.0e-10_rp*matrix_scale) then
         call g_logger%fatal('diagonalize_single_kpoint: H(k) is non-Hermitian before eigensolution.', __FILE__, __LINE__)
      end if

      allocate(h_copy(nmat, nmat), rwork(max(1, 3*nmat - 2)))
      h_copy = hk
      if (use_generalized) then
         allocate(s_copy(nmat, nmat))
         call this%fourier_transform_overlap(k_point, s_copy)
         call this%check_overlap_properties(1, s_copy)
      end if

      if (use_generalized) then
         call zhegv(1, 'V', 'U', nmat, h_copy, nmat, s_copy, nmat, eigenvals, work_query, -1, rwork, info)
      else
         call zheev('V', 'U', nmat, h_copy, nmat, eigenvals, work_query, -1, rwork, info)
      end if
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      h_copy = hk
      if (use_generalized) then
         call this%fourier_transform_overlap(k_point, s_copy)
         call zhegv(1, 'V', 'U', nmat, h_copy, nmat, s_copy, nmat, eigenvals, work, lwork, rwork, info)
      else
         call zheev('V', 'U', nmat, h_copy, nmat, eigenvals, work, lwork, rwork, info)
      end if
      if (info /= 0) then
         call g_logger%fatal('diagonalize_single_kpoint: LAPACK eigensolver failed.', __FILE__, __LINE__)
      end if
      eigenvecs = h_copy
      if (allocated(s_copy)) deallocate(s_copy)
      deallocate(h_copy, rwork, work)
   end subroutine diagonalize_single_kpoint

   !> @brief Build H(k) for every active mesh or path k-point.
   !> @details Allocates hk_bulk/hk_total as needed, dispatches first- or
   !>          second-order Fourier transforms, and applies local k ownership.
   !> @param[inout] this Reciprocal object receiving k-space Hamiltonian arrays.
   module subroutine build_kspace_hamiltonian(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ik, ik_global, nk, ntype, tile_first, tile_last, tile_length
      character(len=200) :: debug_msg
      logical :: using_kpath, distribute_mesh
      integer :: i, j
      type(reciprocal_execution_request) :: request
      type(reciprocal_execution_result) :: result

      call this%validate_nonzero_q_gbt('reciprocal%build_kspace_hamiltonian')
      call this%force_full_bz_for_nonzero_q_gbt('reciprocal%build_kspace_hamiltonian')
      ! The Hamiltonian may have changed between SCF/frozen-magnon probes while
      ! the reciprocal object and k mesh are intentionally reused.  Invalidate
      ! every spectrum-derived object before rebuilding so no old eigenpairs,
      ! DOS moments, or canonical energy can survive the probe boundary.
      call this%invalidate_spectral_cache()

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

      ! Tile-level assembly crosses the execution backend boundary.  The CPU
      ! implementation owns the RF-04 scratch tile; a future device backend
      ! can keep this matrix resident for a following eigensolve.
      call this%make_execution_backend()
      do tile_first = 1, this%nk_local, max(1, this%reciprocal_tile_size)
         tile_last = min(this%nk_local, tile_first + max(1, this%reciprocal_tile_size) - 1)
         tile_length = tile_last - tile_first + 1
         request%assemble_hamiltonian = .true.
         request%assemble_overlap = .false.
         request%solve_eigensystem = .false.
         request%generalized = .false.
         request%request_eigenvectors = .false.
         request%request_assembled_hamiltonian = .true.
         request%request_assembled_overlap = .false.
         request%operator_generation = this%hamiltonian%operator_generation
         if (allocated(request%k_points)) deallocate(request%k_points)
         allocate(request%k_points(3,tile_length))
         if (using_kpath) then
            request%k_points = this%k_path(:,tile_first:tile_last)
         else
            request%k_points = this%k_workset%points(:,tile_first:tile_last)
         end if
         call this%execution_backend%execute_batch(request, result)
         this%hk_bulk(:,:,tile_first:tile_last) = result%hamiltonian
      end do

      ! H(k) now represents exactly this shared real-space operator generation.
      ! All later eigensystem/DOS/density entry points verify this identity.
      if (associated(this%hamiltonian)) then
         this%cached_operator_generation = this%hamiltonian%operator_generation
      end if

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
               call this%fourier_transform_array(this%hamiltonian%eecc, this%k_workset%points(:, 1), hcck_diag)
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

   !> @brief Build S(k) for every active mesh or path k-point.
   !> @param[inout] this Reciprocal object receiving sk_overlap arrays.
   module subroutine build_kspace_overlap(this)
      class(reciprocal), intent(inout) :: this
      integer :: ik, ik_global, nk, tile_first, tile_last, tile_length
      logical :: using_kpath
      type(reciprocal_assembler) :: assembler

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

      call this%make_reciprocal_assembler(assembler)
      do tile_first = 1, this%nk_local, max(1, this%reciprocal_tile_size)
         tile_last = min(this%nk_local, tile_first + max(1, this%reciprocal_tile_size) - 1)
         tile_length = tile_last - tile_first + 1
         call this%workspace%ensure_capacity(size(this%sk_overlap,1), tile_length, .true., this%hamiltonian%operator_generation, &
                                             this%lattice%nn_max, this%lattice%ntype)
         if (using_kpath) then
            this%workspace%points(:,1:tile_length) = this%k_path(:,tile_first:tile_last)
         else
            this%workspace%points(:,1:tile_length) = this%k_workset%points(:,tile_first:tile_last)
         end if
         call assembler%assemble_overlap_batch(this%workspace%points(:,1:tile_length), this%workspace)
         this%sk_overlap(:,:,tile_first:tile_last) = this%workspace%s(:,:,1:tile_length)
      end do
      call root_info('reciprocal%build_kspace_overlap: Built S(k) overlap proxy.', __FILE__, __LINE__)
   end subroutine build_kspace_overlap

!> @brief Diagnose onsite spin-diagonal blocks for multisite Hamiltonians.
!> @param[in] this Reciprocal object containing hk_total or hk_bulk data.
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

!> @brief Check Hermiticity of one k-space Hamiltonian matrix.
!> @param[in] this Reciprocal object containing k-space Hamiltonian data.
!> @param[in] ik k-point index to inspect.
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

!> @brief Print block-structure diagnostics for one k-space Hamiltonian.
!> @param[in] this Reciprocal object containing k-space Hamiltonian data.
!> @param[in] ik k-point index to inspect.
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

   !> @brief Return a basis label from the configured orbital count.
   !> @param[in] this Reciprocal object containing basis_size metadata.
   !> @param[in] ntype Atom type index.
   !> @return Basis label such as sp, spd, or spdf.
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
