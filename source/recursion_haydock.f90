submodule(recursion_mod) recursion_haydock

contains

   !> @brief Apply one block-Lanczos Hamiltonian hop with HOH correction.
   !> @details Forms H|Psi_ll> for the block recursion using the two-sweep
   !>          h - H O H + e_nu + l.s operator, accumulates the block alpha
   !>          coefficient, and updates the active light cone.
   !> @param[inout] this Recursion object; mutates hpsi/hohpsi/pmn_b/atemp_b.
   !> @param[in] ll Current recursion level whose alpha block is produced.
   !> @note Used by the legacy block recursion path for SCF and intersite workflows.
   module subroutine hop_b_hoh(this, ll)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, nr, nnmap, nlimplus1
      integer :: ll ! Recursion step
      integer :: ino ! Atom type
      integer :: ih ! Atom number in the clust
      integer, dimension(0:this%lattice%kk) :: idum
      complex(rp), dimension(nb, nb) :: summ
      complex(rp), dimension(nb, nb) :: locham

      idum(:) = 0
      this%hpsi(:, :, :) = (0.0d0, 0.0d0)
      this%hohpsi(:, :, :) = (0.0d0, 0.0d0)
      this%enupsi(:, :, :) = (0.0d0, 0.0d0)
      this%socpsi(:, :, :) = (0.0d0, 0.0d0)

      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(i, ino, nr, j, nnmap, locham) schedule(dynamic, 100)
         do i = 1, this%lattice%nmax ! Loop in the neighbouring
            idum(i) = this%izero(i)
            ino = this%lattice%iz(i)
            nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
            if (this%izero(i) /= 0) then
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, 1, i), nb, this%psi_b(:, :, i), nb, cone, this%hpsi(:, :, i), nb)
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%enim(:, :, ino), nb, this%psi_b(:, :, i), nb, cone, this%enupsi(:, :, i), nb)
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%lsham(:, :, ino), nb, this%psi_b(:, :, i), nb, cone, this%socpsi(:, :, i), nb)
               if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallcc(:, :, 1, i), nb, this%psi_b(:, :, i), nb, cone, this%enupsi(:, :, i), nb)
            end if
            if (nr >= 2) then
               do j = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, j)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, j, i), nb, this%psi_b(:, :, nnmap), nb, cone, this%hpsi(:, :, i), nb)
                        if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallcc(:, :, j, i), nb, this%psi_b(:, :, nnmap), nb, cone, this%enupsi(:, :, i), nb)
                        idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      !$omp parallel do default(shared) private(i, ih, nr, j, nnmap, locham) schedule(dynamic, 100)
      do i = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         idum(i) = this%izero(i)
         ih = this%lattice%iz(i) ! Atom type
         nr = this%lattice%nn(i, 1) ! Number of neighbours
         if (this%izero(i) /= 0) then
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, 1, ih), nb, this%psi_b(:, :, i), nb, cone, this%hpsi(:, :, i), nb)
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%enim(:, :, ih), nb, this%psi_b(:, :, i), nb, cone, this%enupsi(:, :, i), nb)
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%lsham(:, :, ih), nb, this%psi_b(:, :, i), nb, cone, this%socpsi(:, :, i), nb)
            if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eecc(:, :, 1, ih), nb, this%psi_b(:, :, i), nb, cone, this%enupsi(:, :, i), nb)
         end if
         if (nr >= 2) then
            do j = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(i, j)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, j, ih), nb, this%psi_b(:, :, nnmap), nb, cone, this%hpsi(:, :, i), nb)
                     if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eecc(:, :, j, ih), nb, this%psi_b(:, :, nnmap), nb, cone, this%enupsi(:, :, i), nb)
                     idum(i) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do
      !$omp end parallel do

      ! Mapping update for hoh calculation
      this%izero(:) = idum(:)

      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(i, ino, nr, j, nnmap, locham) schedule(dynamic, 100)
         do i = 1, this%lattice%nmax ! Loop in the neighbouring
            idum(i) = this%izero(i)
            ino = this%lattice%iz(i)
            nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
            if (this%izero(i) /= 0) then
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallo(:, :, 1, i), nb, this%hpsi(:, :, i), nb, cone, this%hohpsi(:, :, i), nb)
            end if
            if (nr >= 2) then
               do j = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, j)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallo(:, :, j, i), nb, this%hpsi(:, :, nnmap), nb, cone, this%hohpsi(:, :, i), nb)
                        idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      !$omp parallel do default(shared) private(i, ih, nr, j, nnmap, locham) schedule(dynamic, 100)
      do i = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         idum(i) = this%izero(i)
         ih = this%lattice%iz(i) ! Atom type
         nr = this%lattice%nn(i, 1) ! Number of neighbours
         if (this%izero(i) /= 0) then
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eeo(:, :, 1, ih), nb, this%hpsi(:, :, i), nb, cone, this%hohpsi(:, :, i), nb)
         end if
         if (nr >= 2) then
            do j = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(i, j)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eeo(:, :, j, ih), nb, this%hpsi(:, :, nnmap), nb, cone, this%hohpsi(:, :, i), nb)
                     idum(i) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do
      !$omp end parallel do

      summ = 0.0d0

      ! Redefines idum for all clust atoms
      this%irnum = 0
      do i = 1, this%lattice%kk
         this%izero(i) = idum(i)
         if (this%izero(i) /= 0) then
            this%irnum = this%irnum + 1
            this%irlist(this%irnum) = i
         end if
      end do

      !$omp parallel do default(shared) private(i) reduction(+:SUMM) schedule(dynamic, 100)
      do i = 1, this%irnum
         ! H = h - hoh + e_nu + l.s
         this%hpsi(:, :, this%irlist(i)) = this%hpsi(:, :, this%irlist(i)) - this%hohpsi(:, :, this%irlist(i)) + this%enupsi(:, :, this%irlist(i)) + this%socpsi(:, :, this%irlist(i))
         !
         this%pmn_b(:, :, this%irlist(i)) = this%hpsi(:, :, this%irlist(i)) - this%pmn_b(:, :, this%irlist(i))
         !
         call zgemm('c', 'n', nb, nb, nb, cone, this%psi_b(:, :, this%irlist(i)), nb, this%hpsi(:, :, this%irlist(i)), nb, cone, summ, nb)
      end do
      !$omp end parallel do

      this%atemp_b(:, :, ll) = summ
   end subroutine

   !> @brief Apply one block-Lanczos Hamiltonian hop without HOH correction.
   !> @details Forms H|Psi_ll> from local and bulk real-space hopping blocks,
   !>          accumulates the block alpha coefficient, and prepares pmn_b for
   !>          the orthogonalization step in crecal_b.
   !> @param[inout] this Recursion object; mutates hpsi/pmn_b/atemp_b and maps.
   !> @param[in] ll Current recursion level whose alpha block is produced.
   !> @note Includes ccor_2c contributions when enabled.
   module subroutine hop_b(this, ll)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, nr, nnmap, nlimplus1
      integer :: ll ! Recursion step
      integer :: ino ! Atom type
      integer :: ih ! Atom number in the clust
      integer, dimension(0:this%lattice%kk) :: idum
      complex(rp), dimension(nb, nb) :: summ
      complex(rp), dimension(nb, nb) :: locham

      idum(:) = 0
      this%hpsi(:, :, :) = (0.0d0, 0.0d0)

      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(i, ino, nr, j, nnmap, locham) schedule(dynamic, 100)
         do i = 1, this%lattice%nmax ! Loop in the neighbouring
            idum(i) = this%izero(i)
            ino = this%lattice%iz(i)
            nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
            if (this%izero(i) /= 0) then
               locham = this%hamiltonian%hall(:, :, 1, i) + this%hamiltonian%lsham(:, :, ino)
               if (this%hamiltonian%ccor_2c) locham = locham + this%hamiltonian%hallcc(:, :, 1, i)
               call zgemm('n', 'n', nb, nb, nb, cone, locham, nb, this%psi_b(:, :, i), nb, cone, this%hpsi(:, :, i), nb)
               !call zgemm(´n´,´n´, nb, nb,18,cone,this%hamiltonian%hall(:,:,1,i),18,this%psi_b(:,:,i),18,cone,this%hpsi(:,:,i),18)
               !call zgemm(´n´,´n´, nb, nb,18,cone,this%hamiltonian%lsham(:,:,ino),18,this%psi_b(:,:,i),18,cone,this%hpsi(:,:,I),18)
            end if
            if (nr >= 2) then
               do j = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, j)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        if (this%hamiltonian%ccor_2c) then
                           call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, j, i) + this%hamiltonian%hallcc(:, :, j, i), nb, this%psi_b(:, :, nnmap), nb, cone, this%hpsi(:, :, i), nb)
                        else
                           call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, j, i), nb, this%psi_b(:, :, nnmap), nb, cone, this%hpsi(:, :, i), nb)
                        end if
                        idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      !$omp parallel do default(shared) private(i, ih, nr, j, nnmap, locham) schedule(dynamic, 100)
      do i = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         idum(i) = this%izero(i)
         ih = this%lattice%iz(i) ! Atom type
         nr = this%lattice%nn(i, 1) ! Number of neighbours
         if (this%izero(i) /= 0) then
            locham = this%hamiltonian%ee(:, :, 1, ih) + this%hamiltonian%lsham(:, :, ih)
            if (this%hamiltonian%ccor_2c) locham = locham + this%hamiltonian%eecc(:, :, 1, ih)
            call zgemm('n', 'n', nb, nb, nb, cone, locham, nb, this%psi_b(:, :, i), nb, cone, this%hpsi(:, :, i), nb)
            !call zgemm(´n´,´n´, nb, nb,18,cone,this%hamiltonian%ee(:,:,1,ih),18,this%psi_b(:,:,i),18,cone,this%hpsi(:,:,i),18)
            !call zgemm(´n´,´n´, nb, nb,18,cone,this%hamiltonian%lsham(:,:,ih),18,this%psi_b(:,:,i),18,cone,this%hpsi(:,:,I),18)
         end if
         if (nr >= 2) then
            do j = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(i, j)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     if (this%hamiltonian%ccor_2c) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, j, ih) + this%hamiltonian%eecc(:, :, j, ih), nb, this%psi_b(:, :, nnmap), nb, cone, this%hpsi(:, :, i), nb)
                     else
                        call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, j, ih), nb, this%psi_b(:, :, nnmap), nb, cone, this%hpsi(:, :, i), nb)
                     end if
                     idum(i) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do
      !$omp end parallel do

      summ = 0.0d0
      ! Redefines idum for all clust atoms
      this%irnum = 0
      do i = 1, this%lattice%kk
         this%izero(i) = idum(i)
         if (this%izero(i) /= 0) then
            this%irnum = this%irnum + 1
            this%irlist(this%irnum) = i
         end if
      end do

      !$omp parallel do default(shared) private(i) reduction(+:SUMM) schedule(dynamic, 100)
      do i = 1, this%irnum
         !do i=1, this%lattice%kk
         this%pmn_b(:, :, this%irlist(i)) = this%hpsi(:, :, this%irlist(i)) - this%pmn_b(:, :, this%irlist(i))
         call zgemm('c', 'n', nb, nb, nb, cone, this%psi_b(:, :, this%irlist(i)), nb, this%hpsi(:, :, this%irlist(i)), nb, cone, summ, nb)
         !call zgemm(´c´,´n´, nb, nb,18,cone,PSI_B(1,1,irlist(I)),18,Hpsi(1,1,irlist(I)),18,cone,SUMM,18)
      end do
      !$omp end parallel do

      this%atemp_b(:, :, ll) = summ
   end subroutine hop_b

   !> @brief Generate block-Lanczos coefficients for the selected recursion atoms.
   !> @details Computes the block tridiagonal Haydock coefficients A_n and B_n
   !>          for each atom in lattice%irec. Real-space SCF and block Green's
   !>          functions consume the resulting a_b/b2_b and diagonal a/b2 views.
   !> @param[inout] this Recursion object; fills local slices of a_b, b2_b, a, and b2.
   !> @note Work is partitioned over MPI ranks and may use CUDA or haydock_fast.
   module subroutine recur_b(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, l, ll, kk, m
      integer :: llmax ! Recursion steps
      real(rp) :: rng

      ! Determine how many atoms each process should handle
      call get_mpi_variables(rank, this%lattice%nrec)

      llmax = this%lattice%control%lld
      if (gpu_plugin_ready(this, 'recur_b()', allow_hoh=.true.)) then
         call gpu_plugin_upload_hamiltonian(this)
         do i = start_atom, end_atom
            j = this%lattice%irec(i)
            call g_logger%info('CUDA block recursion on progress for atom '// &
               int2str(j), __FILE__, __LINE__)
            this%psi_b(:, :, :) = (0.0d0, 0.0d0)
            this%atemp_b(:, :, :) = (0.0d0, 0.0d0)
            this%b2temp_b(:, :, :) = (0.0d0, 0.0d0)
            do l = 1, nb
               this%psi_b(l, l, j) = (1.0d0, 0.0d0)
            end do
            call this%gpu_backend%block_lanczos(this%psi_b, llmax, &
               this%atemp_b, this%b2temp_b, &
               prec=merge(1, 0, trim(this%control%cheb_backend) == 'fast_dp'))
            do ll = 1, llmax
               do l = 1, nb
                  do m = 1, nb
                     this%a_b(l, m, ll, i - start_atom + 1) = this%atemp_b(l, m, ll)
                     this%b2_b(l, m, ll, i - start_atom + 1) = this%b2temp_b(l, m, ll)
                  end do
                  this%a(ll, l, i - start_atom + 1, 1) = real(this%atemp_b(l, l, ll))
                  this%b2(ll, l, i - start_atom + 1, 1) = real(this%b2temp_b(l, l, ll))
               end do
            end do
         end do
         return
      end if

      ! Fast CPU block-Lanczos kernels (haydock_fast). cheb_backend selects:
      !   'fast'    -> fp32 working precision, 'fast_dp' -> fp64. local_axis is
      ! not yet supported in the fast path -> fall through to legacy.
      if (this%hamiltonian%ccor_2c .and. this%hamiltonian%hoh .and. trim(this%control%cheb_backend) /= 'legacy') then
         call g_logger%warning('ccor_2c with hoh is not supported by fast block-Lanczos; falling back to legacy.', __FILE__, __LINE__)
      end if
      if (trim(this%control%cheb_backend) /= 'legacy' .and. .not. this%hamiltonian%local_axis .and. &
          .not. (this%hamiltonian%ccor_2c .and. this%hamiltonian%hoh)) then
         if (this%hamiltonian%ccor_2c) call ensure_ccor_operator_blocks(this)
         do i = start_atom, end_atom
            j = this%lattice%irec(i)
            call g_logger%info('Fast block recursion on progress for atom '//int2str(j), __FILE__, __LINE__)
            this%psi_b(:, :, :) = (0.0d0, 0.0d0)
            this%atemp_b(:, :, :) = (0.0d0, 0.0d0)
            this%b2temp_b(:, :, :) = (0.0d0, 0.0d0)
            do l = 1, nb
               this%psi_b(l, l, j) = (1.0d0, 0.0d0)
            end do
            if (this%hamiltonian%hoh) then
               call block_lanczos_fast(this%psi_b, llmax, this%atemp_b, this%b2temp_b, &
                  this%hamiltonian%ee, this%hamiltonian%hall, this%hamiltonian%lsham, &
                  this%lattice%nn, this%lattice%iz, this%lattice%kk, nb, size(this%lattice%nn, 2), &
                  this%lattice%ntype, this%lattice%nmax, &
                  trim(this%control%cheb_backend) == 'fast', .true., &
                  this%hamiltonian%eeo, this%hamiltonian%hallo, this%hamiltonian%enim)
            else
               call block_lanczos_fast_nohoh(this, llmax, this%hamiltonian%ccor_2c)
            end if
            do ll = 1, llmax
               do l = 1, nb
                  do m = 1, nb
                     this%a_b(l, m, ll, i - start_atom + 1) = this%atemp_b(l, m, ll)
                     this%b2_b(l, m, ll, i - start_atom + 1) = this%b2temp_b(l, m, ll)
                  end do
                  this%a(ll, l, i - start_atom + 1, 1) = real(this%atemp_b(l, l, ll))
                  this%b2(ll, l, i - start_atom + 1, 1) = real(this%b2temp_b(l, l, ll))
               end do
            end do
         end do
         return
      end if

      !do i=1, this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
      do i = start_atom, end_atom ! Loop on the number of atoms to be treat self-consistently
         j = this%lattice%irec(i) ! Atom number in the clust file
         ! Clear list of atoms in hopping region
         call g_logger%info('Block recursion on progress for atom '//int2str(j), __FILE__, __LINE__)
         this%izero(:) = 0
         ! Initializing wave functions
         this%psi_b(:, :, :) = (0.0d0, 0.0d0)
         this%pmn_b(:, :, :) = (0.0d0, 0.0d0)

         ! Rotate Hamiltonian if needed
         if (this%hamiltonian%local_axis) then
            call this%hamiltonian%rotate_to_local_axis(this%lattice%symbolic_atoms(i)%potential%mom)
         end if

         do l = 1, nb
            this%psi_b(l, l, j) = (1.0d0, 0.0d0)
            this%atemp_b(l, l, llmax) = (0.0d0, 0.0d0)
            this%b2temp_b(l, l, 1) = (1.0d0, 0.0d0)
         end do

         this%izero(j) = 1

         call this%crecal_b()

         do ll = 1, llmax
            do l = 1, nb
               do m = 1, nb
                  this%a_b(l, m, ll, i - start_atom + 1) = (this%atemp_b(l, m, ll))
                  this%b2_b(l, m, ll, i - start_atom + 1) = (this%b2temp_b(l, m, ll))
               end do
               this%a(ll, l, i - start_atom + 1, 1) = real(this%atemp_b(l, l, ll))
               this%b2(ll, l, i - start_atom + 1, 1) = real(this%b2temp_b(l, l, ll))
            end do
         end do
      end do ! End of the loop on the nrec
      ! For debug purposes
      do i = start_atom, end_atom
         do l = 1, nb  ! Loop on the orbital l
            write (1000*(1 + rank) + 122, *) 'orbital', l, 'atom', i
            do ll = 1, llmax ! Loop on the recursion steps
               write (1000*(1 + rank) + 122, '(2f12.8,4x,2f12.8)') this%a(ll, l, i - start_atom + 1, 1), this%b2(ll, l, i - start_atom + 1, 1)
               !write(1000*(1+rank)+122, ´(2f12.8,4x,2f12.8)´) this%a_b(l,l,ll,i-start_atom+1), this%b2_b(l,l,ll,i-start_atom+1)
          !!!write(1000*(1+rank)+122, *) this%a(ll, l, i,1-start_atom+1), this%b2(ll,l,i,1)
            end do
         end do
      end do
   end subroutine recur_b

   !> @brief Advance the legacy block-Lanczos recurrence for one starting block.
   !> @details Orthogonalizes H|Psi_n> against the current and previous block
   !>          states, diagonalizes the residual norm matrix to obtain B_{n+1},
   !>          and stores temporary A_n/B_n blocks for the caller.
   !> @param[inout] this Recursion object with psi_b initialized for one atom or pair.
   !> @note Uses LAPACK zheev and mutates psi_b, pmn_b, atemp_b, and b2temp_b.
   module subroutine crecal_b(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, info, lwork
      integer :: ll, llmax ! Recursion step
      integer :: nm1 ! LL-1
      integer :: nat ! Clust size
      integer :: hblocksize ! Hamiltonian size (18)
      complex(rp), dimension(nb, nb) :: sum_b, sum_a
      complex(rp), dimension(nb, nb) :: dum, b, b_i, u, lam, lam_i, u_ls, u_rs
!    complex(rp), dimension(nb, nb,this%lattice%kk) :: psi_t
      complex(rp), dimension(:, :, :), allocatable :: psi_t
      real(rp), dimension(nb) :: ev
      complex(rp), dimension(nb) :: zev
      real(rp), dimension(3*nb - 2) ::rwork
      complex(rp), dimension(3*nb - 2) ::zwork

      allocate (psi_t(nb, nb, this%lattice%kk))

      sum_b = this%b2temp_b(:, :, 1)

      hblocksize = nb
      nat = this%lattice%kk
      nm1 = this%lattice%control%lld - 1
      llmax = this%lattice%control%lld
      do ll = 1, nm1
         this%atemp_b(:, :, ll) = (0.0d0, 0.0d0)
         !call HOP to get H*w = A; psi=psi ; pmn=Hpsi - B_n-1|Psi_n-1>
         !PMN = H Wn - W_(n-1)B_(n-1) in Inoue syntax  (is zero for first recursion step)
         call g_timer%start('H|PSI_n>')
         if (this%hamiltonian%hoh) then
            call this%hop_b_hoh(ll) ! a_n=<PSI_n|H|PSI_n>
         else
            call this%hop_b(ll)     ! a_n=<PSI_n|H|PSI_n>
         end if
         call g_timer%stop('H|PSI_n>')

         !psi_t(:,:,:) = this%psi_b(:,:,:)
         !call zcopy(nb*nb*this%lattice%kk,this%psi_b,1,psi_t,1)
         call g_timer%start('H|Psi_n-A_n|Psi_n-B_n|Psi_n-1')
         do i = 1, this%irnum
            psi_t(:, :, this%irlist(i)) = this%psi_b(:, :, this%irlist(i))
         end do
         !psi_t(:,:,:) = this%psi_b(:,:,:)

         this%b2temp_b(:, :, ll) = sum_b(:, :)

         sum_b(:, :) = 0.0d0

         !$omp parallel do default(shared) private(i) reduction(+:SUM_b)
         do i = 1, this%irnum
            !do i=1, nat
            !  ( H|Psi_n>-A_n|Psi_n>-B_n-1|Psi_n-1> = PMN-A_n|Psi_n> )
         !! H Wn - Wn An  - Wn-1 Bn-1 in Inoue syntax
            call zgemm('n', 'n', nb, nb, nb, cmone, this%psi_b(:, :, this%irlist(i)), nb, this%atemp_b(:, :, ll), nb, cone, this%pmn_b(:, :, this%irlist(i)), nb)
            !call zgemm(´n´,´n´, nb, nb,18,cmone,this%psi_b(:,:,i),18,this%atemp_b(:,:,ll),18,cone,this%pmn_b(:,:,i), 18)
            ! B^2_n+1 = ( H|Psi_n-A_n|Psi_n-B_n|Psi_n-1 )´*( H|Psi_n-A_n|Psi_n-B_n|Psi_n-1 )
         !! (H Wn - An Wn  - Wn-1 Bn-1)*(H Wn - An Wn - Wn-1 Bn-1)´ -in Inoue syntax
            call zgemm('c', 'n', nb, nb, nb, cone, this%pmn_b(:, :, this%irlist(i)), nb, this%pmn_b(:, :, this%irlist(i)), nb, cone, sum_b, nb)
            !call zgemm(´c´,´n´, nb, nb,18,cone,this%pmn_b(:,:,i),18,this%pmn_b(:,:,i),18,cone,sum_b,18)
         end do
         !$omp end parallel do
         call g_timer%stop('H|Psi_n-A_n|Psi_n-B_n|Psi_n-1')
         call g_timer%start('B_n+1')
         ! Calculate B_n+1
         u(:, :) = sum_b(:, :)
         ! Replace sqrt with eigen solver to get lamda^2 and U
         ! get lamda^2 and U in lamda=U*BB´*U*
         call zheev('v', 'u', nb, u, nb, ev, dum, nb*nb, rwork, info)
         if (info /= 0) call g_logger%fatal('Diagonalization error', __FILE__, __LINE__)
         !
         !
         !

         lam = (0.0d0, 0.0d0)
         lam_i = (0.0d0, 0.0d0)
         do i = 1, nb
            lam(i, i) = cmplx(sqrt(max(0.0_rp, ev(i))), kind=kind(0.0d0))
            if (abs(lam(i, i)) > 0.0_rp) then
               lam_i(i, i) = 1.0d0/lam(i, i) ! 1/B_n
            else
               lam_i(i, i) = (0.0d0, 0.0d0)
            end if
         end do
         !
         !  Calc. U*lamda*U´= B
         call zgemm('n', 'n', nb, nb, nb, cone, u, nb, lam, nb, czero, dum, nb)
         call zgemm('n', 'c', nb, nb, nb, cone, dum, nb, u, nb, czero, b, nb)
         !  B^-1=U*lamda^-1*U´
         call zgemm('n', 'n', nb, nb, nb, cone, u, nb, lam_i, nb, czero, dum, nb)
         call zgemm('n', 'c', nb, nb, nb, cone, dum, nb, u, nb, czero, b_i, nb)
         call g_timer%stop('B_n+1')

         call g_timer%start('<PSI|B_n+1|PSI>')
         !$omp parallel do default(shared) private(i)
         do i = 1, this%irnum
            !do i=1, nat
            call zgemm('n', 'n', nb, nb, nb, cone, this%pmn_b(:, :, this%irlist(i)), nb, b_i, nb, czero, this%psi_b(:, :, this%irlist(i)), nb)
            call zgemm('n', 'n', nb, nb, nb, cone, psi_t(:, :, this%irlist(i)), nb, b, nb, czero, this%pmn_b(:, :, this%irlist(i)), nb)
         end do
         !$omp end parallel do
         call g_timer%stop('<PSI|B_n+1|PSI>')
      end do
      this%b2temp_b(:, :, llmax) = sum_b(:, :)
   end subroutine crecal_b

   !> @brief Replace stored block B^2 coefficients by their matrix square roots.
   !> @details Diagonalizes each Hermitian B2_b block and overwrites it with the
   !>          positive square-root matrix B. Block Green's-function reconstruction
   !>          calls this before continued-fraction evaluation.
   !> @param[inout] this Recursion object; overwrites local B2_b slices.
   !> @note Uses MPI partition sizes and LAPACK zheev; no collective communication occurs here.
   module subroutine zsqr(this)
      use mpi_mod
      use, intrinsic :: ieee_exceptions, only: ieee_divide_by_zero, &
                                               ieee_get_halting_mode, ieee_set_flag, ieee_set_halting_mode
      implicit none
      class(recursion), intent(inout) :: this
      !
      integer :: na, LDIM, I, J, K, L, N, M, info, n_glob
      real(rp), dimension(nb) :: ev
      real(rp), dimension(3*nb - 2) ::rwork
      complex(rp), dimension(nb*nb) ::zwork
      complex(rp), dimension(nb, nb) :: B2t, U, DUM, lam
      logical :: halt_divide_by_zero

      !
      if (this%lattice%njij == 0) then
         na = atoms_per_process !this%lattice%nrec
      else
         na = atoms_per_process*4 !4*this%lattice%njij
      end if
      !call get_mpi_variables(rank,na)
      lam = 0.0d0
      LDIM = nb
      ! do n_glob=start_atom, end_atom
      !  N = g2l_map(n_glob)
      do N = 1, na
         do L = 1, this%lattice%control%lld
            do I = 1, LDIM
               do J = 1, LDIM
                  U(J, I) = this%B2_b(J, I, L, N)
               end do
            end do
            ! replace sqrt with eigen solver to get lamda^2 and U
            ! get lamda^2 and U in lamda=U*BB´*U*
            ! oneMKL's zheev can evaluate an intermediate divide-by-zero in
            ! its scaling path; mask that external-library exception only.
            call ieee_get_halting_mode(ieee_divide_by_zero, halt_divide_by_zero)
            call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
            call zheev('V', 'U', LDIM, U, LDIM, ev, zwork, LDIM*LDIM, rwork, info)
            call ieee_set_flag(ieee_divide_by_zero, .false.)
            call ieee_set_halting_mode(ieee_divide_by_zero, halt_divide_by_zero)
            !
            if (info /= 0) print *, 'diag', info
            ! lam=sqrt(lamda^2) ; lam_i=1/sqrt(lamda^2)
            do I = 1, LDIM
               lam(I, I) = sqrt(max(0.0_rp, ev(i)))
            end do
            ! calc. U*lamda*U´= B
            call zgemm('n', 'n', LDIM, LDIM, LDIM, (1.0d0, 0.0d0), U, LDIM, lam, LDIM, (0.0d0, 0.0d0), DUM, LDIM)
            call zgemm('n', 'c', LDIM, LDIM, LDIM, (1.0d0, 0.0d0), DUM, LDIM, U, LDIM, (0.0d0, 0.0d0), this%B2_b(1, 1, L, N), LDIM)
         end do
      end do
   end subroutine zsqr

   !> @brief Estimate scalar terminator coefficients for continued fractions.
   !> @details Reduces finite Lanczos alpha/beta sequences to asymptotic constants
   !>          used by Green-function terminators. The result approximates the tail
   !>          of the recursion beyond the explicitly computed ll levels.
   !> @param[inout] this Recursion object; provides access to bpopt/emami helpers.
   !> @param[inout] alpha Recursion alpha coefficients, shape (np,ll).
   !> @param[in] beta Recursion beta coefficients, shape (np,ll).
   !> @param[inout] ll Number of available recursion levels.
   !> @param[in] np Number of scalar coefficient channels.
   !> @param[inout] nw Work-array size used by the legacy terminator routine.
   !> @param[out] alpha_inf Asymptotic alpha per channel.
   !> @param[out] beta_inf Asymptotic beta per channel.
   module subroutine get_cinf(this, alpha, beta, ll, np, nw, alpha_inf, beta_inf)
      implicit none
      class(recursion), intent(inout) :: this
      !.. Formal Arguments ..
      integer, intent(in) :: np
      integer, intent(inout) :: ll, nw
      real(rp), dimension(np, ll), intent(inout) :: alpha
      real(rp), dimension(np, ll), intent(in) :: beta
      real(rp), dimension(np), intent(out) :: alpha_inf
      real(rp), dimension(np), intent(out) :: beta_inf
      !
      !.. Local Scalars ..
      integer :: icode, iii, k, l, ll1, ineigh, nbp1, nl, nt, eidx, ne, nq, ifail
      real(rp) :: a1, a2, emax, emin, eps, err, e_shift, pi
      complex(rp) :: g_e
      !
      !.. Local Arrays ..
      integer, dimension(np) :: nb2
      integer, dimension(nw) :: jc
      integer, dimension(nw) :: iwk
      real(rp), dimension(ll) :: aa, am, bb, bm, sqbb
      real(rp), dimension(10) :: edge, weight, width
      real(rp), dimension(200, 10) :: bwk
      ! real(rp), dimension(np,ll) :: am2,bm2
      real(rp), dimension(np, 10) :: edge2, width2, weight2
      real(rp), dimension(nw, 2, 5) :: work
      real(rp), dimension(2) :: bandedges

      !
      !.. External Calls ..
      !
      !.. External Functions ..
      !
      !.. Intrinsic Functions ..
      intrinsic MAX, MIN
      !
      ! ... Executable Statements ...
      !
      !**************************************************************************
      icode = 0
      !*******************************************************
      eps = 1.0d-14
      err = 0.00001d0
      nbp1 = 2   !( Nband+1)
      !*************************************************************************
      do nl = 1, np
         do l = 1, ll
            aa(l) = alpha(nl, l)
            bb(l) = beta(nl, l)!**2
         end do
         ll1 = ll
         am = 0.0d0; bm = 0.0d0; edge = 0.0d0; width = 0.0d0; weight = 0.0d0
         call this%bpOPT(ll, AA, BB, LL - 1, alpha_inf(nl), beta_inf(nl), ifail)
         edge(1) = alpha_inf(nl) - 2.0d0*beta_inf(nl)
         width(1) = 4.0d0*beta_inf(nl)
      end do
   end subroutine get_cinf

   !> @brief Estimate block terminators for Haydock continued fractions.
   !> @details Converts complex block recursion coefficients to real scalar
   !>          channel sequences, calls get_cinf, sanitizes NaNs, and returns
   !>          diagonal/off-diagonal terminator matrices for block Green functions.
   !> @param[inout] this Recursion object; only helper state is used.
   !> @param[in] Acoef_b Block alpha coefficients, shape (ldim,ldim,ll,na).
   !> @param[in] B2coef_b Block beta/B coefficients, shape (ldim,ldim,ll,na).
   !> @param[in] na Number of atoms or intersite phase combinations.
   !> @param[in] ll Number of recursion levels.
   !> @param[in] ldim Block dimension, usually nb.
   !> @param[inout] nw Legacy work-array size.
   !> @param[out] a_inf Matrix terminator alpha values.
   !> @param[out] b_inf Matrix terminator beta values.
   !> @param[out] a_inf0 Site-averaged diagonal alpha terminator.
   !> @param[out] b_inf0 Site-averaged diagonal beta terminator.
   module subroutine get_terminf(this, Acoef_b, B2coef_b, na, ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)
      implicit none
      class(recursion), intent(inout) :: this
      integer, intent(in) :: na
      integer, intent(in) :: ll
      integer, intent(inout) :: nw
      integer, intent(in) :: ldim
      complex(rp), dimension(ldim, ldim, ll, na), intent(in) :: Acoef_b, B2coef_b
      real(rp), dimension(ldim, ldim, na), intent(out) :: a_inf, b_inf
      real(rp), dimension(na), intent(out) ::  a_inf0, b_inf0
      !
      integer :: n, i, j, ll_t, k
      complex(rp), dimension(ldim, ldim) :: MatIn, MatOut
      real(rp), dimension(ldim, ldim, ll) :: Acoef_r, B2coef_r
      real(rp) :: maxA, maxB, maxAinf, maxBinf, tmpval
      logical :: foundNaN_in, foundNaN_out
      !
      a_inf = 0.0d0; b_inf = 0.0d0
      do n = 1, na
         ll_t = ll
         do i = 1, ll
            MatIn = real(Acoef_b(:, :, i, n))
            Acoef_r(:, :, i) = MatIn
            MatIn = real(B2coef_b(:, :, i, n))
            B2coef_r(:, :, i) = MatIn
         end do
         ! Summarize input recursion coefficients (Acoef_r, B2coef_r)
         maxA = 0.0_rp
         maxB = 0.0_rp
         foundNaN_in = .false.
         do i = 1, ldim
            do j = 1, ldim
               do k = 1, ll
                  tmpval = abs(Acoef_r(i, j, k))
                  if (tmpval > maxA) maxA = tmpval
                  if (IsNaN(Acoef_r(i, j, k))) foundNaN_in = .true.
                  tmpval = abs(B2coef_r(i, j, k))
                  if (tmpval > maxB) maxB = tmpval
                  if (IsNaN(B2coef_r(i, j, k))) foundNaN_in = .true.
               end do
            end do
         end do
         ! write(*,*) 'DEBUG:get_terminf INPUT n=', n, ' maxA=', maxA, ' maxB=', maxB, ' NaN_in=', foundNaN_in
         call this%get_cinf(Acoef_r, B2coef_r, ll_t, ldim*ldim, nw, a_inf(:, :, n), b_inf(:, :, n))
         do j = 1, ldim
            do i = 1, ldim
               if (IsNaN(a_inf(i, j, n))) a_inf(i, j, n) = 0.0d0
               if (IsNaN(b_inf(i, j, n))) b_inf(i, j, n) = 0.0d0
            end do
            if (a_inf(j, j, n) == 0.0d0) a_inf(j, j, n) = 0.5d0
            if (b_inf(j, j, n) == 0.0d0) b_inf(j, j, n) = 0.5d0
         end do
         ! Summarize outputs a_inf / b_inf
         maxAinf = maxval(abs(a_inf(:, :, n)))
         maxBinf = maxval(abs(b_inf(:, :, n)))
         foundNaN_out = .false.
         do i = 1, ldim
            do j = 1, ldim
               if (IsNaN(a_inf(i, j, n)) .or. IsNaN(b_inf(i, j, n))) foundNaN_out = .true.
            end do
         end do
         ! write(*,*) 'DEBUG:get_terminf OUTPUT n=', n, ' a_inf0_avg=', a_inf0(n), ' maxAinf=', maxAinf, ' maxBinf=', maxBinf, ' NaN_out=', foundNaN_out
         a_inf0(n) = 0.0d0
         do i = 1, ldim
            a_inf0(n) = a_inf0(n) + a_inf(i, i, n)
         end do
         a_inf0(n) = a_inf0(n)/ldim
         b_inf(1, 1, n) = b_inf(1, 1, n)*1.01d0
         if (ldim >= 10) then
            b_inf(10, 10, n) = b_inf(10, 10, n)*1.01d0
         end if
         b_inf0(n) = 0.0d0
         do i = 1, ldim
            b_inf0(n) = b_inf0(n) + b_inf(i, i, n)
         end do
         b_inf0(n) = b_inf0(n)/ldim
      end do
   end subroutine get_terminf

   !> @brief Precompute light-cone support maps for Chebyshev recursion levels.
   !> @details Propagates the active-site mask one neighbor shell per level and
   !>          stores the result in izeroll(:,ll). Legacy Chebyshev routines use
   !>          this map to avoid visiting sites outside the reachable support.
   !> @param[inout] this Recursion object; fills izeroll.
   module subroutine create_ll_map(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, nr, nnmap, ll
      integer, dimension(0:this%lattice%kk) :: idumll

      idumll(:) = 0

      do ll = 1, this%lattice%control%lld
         idumll(:) = 0
         do i = 1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
            idumll(i) = this%izeroll(i, ll)
            nr = this%lattice%nn(i, 1) ! Number of neighbours
            if (nr >= 2) then
               do j = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, j)
                  if (nnmap /= 0) then
                     if (this%izeroll(nnmap, ll) /= 0) then
                        idumll(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do
         this%izeroll(:, ll + 1) = idumll(:)
      end do
   end subroutine create_ll_map

   !> @brief Apply one scalar Lanczos Hamiltonian hop.
   !> @details Forms H|psi_ll> over the current light cone, accumulates the scalar
   !>          alpha coefficient for each orbital channel, and prepares pmn for
   !>          the scalar Haydock recurrence.
   !> @param[inout] this Recursion object; mutates v, pmn, atemp, and izero.
   !> @param[in] ll Current recursion level whose alpha coefficient is produced.
   !> @note This legacy scalar path is used by control%recur='lanczos'.
   module subroutine hop(this, ll)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, nr, nnmap, nlimplus1
      integer :: ll ! Recursion step
      integer :: ino ! Atom type
      integer :: ih ! Atom number in the clust
      integer, dimension(0:this%lattice%kk) :: idum
      complex(rp), dimension(nb) :: dum
      real(rp) :: summ, start, finish

      summ = 0.0d0
      idum(:) = 0

      nlimplus1 = this%lattice%nmax + 1

      select case (this%lattice%control%nsp)

      case (1) ! Scalar relativistic case

         if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
            !$omp parallel do default(shared) private(i, nr, m, l, dum, j, nnmap) schedule(dynamic, 100)
            do i = 1, this%lattice%nmax ! Loop in the neighbouring
               idum(i) = this%izero(i)
               dum(:) = (0.0d0, 0.0d0)
               nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
               if (this%izero(i) /= 0) then
                  do m = 1, norb ! Loop on the orbital m
                     do l = 1, norb ! Loop on the orbital l
                                 dum(l) = dum(l) + this%hamiltonian%hall(l, m, 1, i)*this%psi(m, i)
                                 dum(l +spin_off) = dum(l +spin_off) + this%hamiltonian%hall(l +spin_off, m +spin_off, 1, i)*this%psi(m +spin_off, i)
                                 if (this%hamiltonian%ccor_2c) then
                                    dum(l) = dum(l) + this%hamiltonian%hallcc(l, m, 1, i)*this%psi(m, i)
                                    dum(l +spin_off) = dum(l +spin_off) + this%hamiltonian%hallcc(l +spin_off, m +spin_off, 1, i)*this%psi(m +spin_off, i)
                                 end if
                     end do ! End of loop on orbital m
                  end do ! End of loop on orbital l
               end if
               if (nr >= 2) then
                  do j = 2, nr ! Loop on the neighbouring
                     nnmap = this%lattice%nn(i, j)
                     if (nnmap /= 0) then
                        if (this%izero(nnmap) /= 0) then
                           do m = 1, norb ! Loop in the orbital m
                              do l = 1, norb ! Loop in the orbital l
                                 dum(l) = dum(l) + this%hamiltonian%hall(l, m, j, i)*this%psi(m, nnmap)
                                 dum(l +spin_off) = dum(l +spin_off) + this%hamiltonian%hall(l +spin_off, m +spin_off, j, i)*this%psi(m +spin_off, nnmap)
                                 if (this%hamiltonian%ccor_2c) then
                                    dum(l) = dum(l) + this%hamiltonian%hallcc(l, m, j, i)*this%psi(m, nnmap)
                                    dum(l +spin_off) = dum(l +spin_off) + this%hamiltonian%hallcc(l +spin_off, m +spin_off, j, i)*this%psi(m +spin_off, nnmap)
                                 end if
                              end do ! End of loop in orbital m
                           end do ! End of loop in orbital l
                           idum(i) = 1
                        end if
                     end if
                  end do ! End of loop in the neighbouring
               end if
               do l = 1, nb
                  this%v(l, i) = dum(l)
               end do
            end do ! End of loop in the neighbouring
            !$omp end parallel do
         end if ! End of local Hamiltonian loop

         !$omp parallel do default(shared) private(i, ih, nr, m, l, dum, j, nnmap) schedule(dynamic, 100)
         do i = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
            idum(i) = this%izero(i)
            ih = this%lattice%iz(i) ! Atom type
            dum(:) = (0.0d0, 0.0d0)
            nr = this%lattice%nn(i, 1) ! Number of neighbours
            if (this%izero(i) /= 0) then
               !write(125, *) ´i is ´, i
               do m = 1, norb ! Loop on the orbital m
                  do l = 1, norb ! Loop on the orbital l
                     dum(l) = dum(l) + this%hamiltonian%ee(l, m, 1, ih)*this%psi(m, i)
                     dum(l +spin_off) = dum(l +spin_off) + this%hamiltonian%ee(l +spin_off, m +spin_off, 1, ih)*this%psi(m +spin_off, i)
                     if (this%hamiltonian%ccor_2c) then
                        dum(l) = dum(l) + this%hamiltonian%eecc(l, m, 1, ih)*this%psi(m, i)
                        dum(l +spin_off) = dum(l +spin_off) + this%hamiltonian%eecc(l +spin_off, m +spin_off, 1, ih)*this%psi(m +spin_off, i)
                     end if
                  end do ! End of the loop on the orbital l
               end do ! End of loop on the orbital m
            end if
            if (nr >= 2) then
               do j = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, j)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        !write(125, *) i, j, this%lattice%nn(i, j)
                        do m = 1, norb ! Loop on the orbital m
                           do l = 1, norb ! Loop on orbital l
                              dum(l) = dum(l) + this%hamiltonian%ee(l, m, j, ih)*this%psi(m, nnmap)
                              dum(l +spin_off) = dum(l +spin_off) + this%hamiltonian%ee(l +spin_off, m +spin_off, j, ih)*this%psi(m +spin_off, nnmap)
                              if (this%hamiltonian%ccor_2c) then
                                 dum(l) = dum(l) + this%hamiltonian%eecc(l, m, j, ih)*this%psi(m, nnmap)
                                 dum(l +spin_off) = dum(l +spin_off) + this%hamiltonian%eecc(l +spin_off, m +spin_off, j, ih)*this%psi(m +spin_off, nnmap)
                              end if
                           end do ! End of loop on orbital l
                        end do ! End of loop on the orbital m
                        idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
            do l = 1, nb
               this%v(l, i) = dum(l)
            end do
         end do
         !$omp end parallel do

         ! Redefines idum for all clust atoms
         do i = 1, this%lattice%kk
            this%izero(i) = idum(i)
            do l = 1, nb
               dum(l) = this%v(l, i)
               summ = summ + real(dum(l)*conjg(this%psi(l, i)))
               this%pmn(l, i) = dum(l) + this%pmn(l, i)
            end do
         end do
         this%atemp(ll) = summ
      end select
   end subroutine hop

   !> @brief Advance the legacy scalar Lanczos recurrence for one orbital seed.
   !> @details Iteratively calls hop, subtracts the alpha projection, normalizes
   !>          the residual, and stores temporary scalar alpha/beta coefficients
   !>          for the caller.
   !> @param[inout] this Recursion object with psi initialized for one atom/orbital.
   !> @note Mutates psi, pmn, atemp, and b2temp.
   module subroutine crecal(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m
      integer :: ll, llmax ! Recursion step
      integer :: nm1 ! LL-1
      integer :: nat ! Clust size
      integer :: hblocksize ! Hamiltonian size (18)
      real(rp) :: s, summ, start, finish
      complex(rp) :: dum, ajc
      complex(rp), dimension(nb, this%lattice%kk) :: thpsi
      character(len=1) :: transa
      character, dimension(6) :: matdescra
      ! External functions
      complex(rp), external :: zdotc

      summ = this%b2temp(1)
      thpsi(:, :) = (0.0d0, 0.0d0)

      hblocksize = nb
      nat = this%lattice%kk
      nm1 = this%lattice%control%lld - 1
      llmax = this%lattice%control%lld
      do ll = 1, nm1
!     call g_timer%start(´hop´)
         call this%hop(ll)
!     call g_timer%stop(´hop´)
         this%b2temp(ll) = summ
         ajc = -this%atemp(ll)

         call zaxpy(nat*hblocksize, ajc, this%psi, 1, this%pmn, 1)

         summ = 0.0d0

         ! In case zdotc function lapack is not working
         do i = 1, nat
            do k = 1, nb
               summ = summ + real(conjg(this%pmn(k, i))*this%pmn(k, i))
            end do
         end do

         ! summ = real(zdotc(nat*hblocksize, this%pmn, 1, this%pmn, 1))

         s = 1.0d0/sqrt(MAX(summ, TINY(1.0_rp)))

         thpsi(:, :) = this%pmn(:, :)*s
         this%pmn(:, :) = this%psi(:, :)
         this%psi(:, :) = thpsi(:, :)

         s = sqrt(summ)

         this%pmn(:, :) = -this%pmn(:, :)*s
      end do

      this%b2temp(llmax) = summ
   end subroutine crecal

   !> @brief Generate scalar Lanczos coefficients for selected recursion atoms.
   !> @details Runs the Haydock scalar recursion for each orbital of every atom in
   !>          lattice%irec, filling a and b2 for the legacy lanczos Green-function
   !>          path used by real-space SCF and DOS.
   !> @param[inout] this Recursion object; fills local slices of a and b2.
   !> @note Work is MPI-partitioned and may dispatch to the scalar CUDA kernel for nsp=1.
   module subroutine recur(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, l, ll
      integer :: llmax ! Recursion steps
      integer :: i_loc

      llmax = this%lattice%control%lld

      if (gpu_plugin_ready(this, 'recur()', require_nsp1=.true.)) then
         call gpu_plugin_upload_hamiltonian(this)
         do i = start_atom, end_atom
            i_loc = g2l_map(i)
            j = this%lattice%irec(i)
            call g_logger%info('CUDA scalar Lanczos recursion on progress for atom '// &
               int2str(j), __FILE__, __LINE__)
            call this%gpu_backend%scalar_lanczos(j, llmax, this%a(:, :, i_loc, 1), &
               this%b2(:, :, i_loc, 1))
         end do
         return
      end if

      !do i=1, this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
      do i = start_atom, end_atom
         i_loc = g2l_map(i)
         j = this%lattice%irec(i) ! Atom number in the clust file
         do l = 1, nb ! Loop on the orbitals
            ! Clear list of atoms in hopping region

            this%izero(:) = 0
            ! Initializing wave functions
            this%psi(:, :) = (0.0d0, 0.0d0)
            this%pmn(:, :) = (0.0d0, 0.0d0)

            this%psi(l, j) = (1.0d0, 0.0d0)
            this%izero(j) = 1
            this%atemp(:) = 0.0d0; this%b2temp(:) = 0.0d0

            this%b2temp(1) = 1.0d0
            this%atemp(llmax) = 0.0d0
            !write(125, *) ´orbital ´, l
            call this%crecal()

            do ll = 1, llmax ! Loop on the recursion steps
               this%a(ll, l, i_loc, 1) = this%atemp(ll)
               this%b2(ll, l, i_loc, 1) = this%b2temp(ll)
            end do ! End of the loop on the recursion steps
         end do ! End of the loop on the orbitals
      end do ! End of the loop on the nrec

    !!!! For debug purposes
    !!!do i=1, this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
    !!!  do l = 1, nb  ! Loop on the orbital l
    !!!    write(123, *) ´orbital´, l
    !!!    do ll=1, llmax ! Loop on the recursion steps
    !!!      write(123, *) this%a(ll, l, i_loc, 1), this%b2(ll, l, i_loc, 1)
    !!!    end do
    !!!  end do
    !!!end do
   end subroutine

   !> @brief Optimize scalar terminator constants from a finite coefficient chain.
   !> @details Legacy Beer-Pettifor-style helper used by get_cinf to estimate the
   !>          asymptotic alpha and beta values that close a continued fraction.
   !> @param[inout] this Recursion object; no persistent state is modified.
   !> @param[in] ll Number of coefficients in the input arrays.
   !> @param[in] a Alpha coefficient sequence.
   !> @param[in] rb Beta coefficient sequence.
   !> @param[in] n Number of coefficients to use in the fit.
   !> @param[out] ainf Estimated asymptotic alpha.
   !> @param[out] rbinf Estimated asymptotic beta.
   !> @param[out] ifail Legacy status flag from the optimizer.
   module subroutine bpopt(this, ll, a, rb, n, ainf, rbinf, ifail)
      class(recursion), intent(inout) :: this
      ! Input
      integer, intent(in) :: n, ll
      real(rp), dimension(ll), intent(in) :: A, RB
      ! Output
      integer, intent(out) :: ifail
      real(rp), intent(out) :: ainf, rbinf
      ! Local variables
      integer :: I, JITER, NDIME
      real(rp) :: BM, BMAX, BMIN, EPS
      real(rp), dimension(ll) :: AZ, RBZ

      NDIME = ll
      IFAIL = 0
      EPS = 1.0d-05
      JITER = 0
      AINF = A(N)
      do
         JITER = JITER + 1
         AZ(1) = 0.5d0*(A(1) - AINF)
         do I = 2, N - 1
            AZ(I) = 0.5d0*(A(I) - AINF)
            RBZ(I) = 0.5d0*RB(I)
         end do
         AZ(N) = A(N) - AINF
         RBZ(N) = 1.0d0/SQRT(2.0d0)*RB(N)
         call this%EMAMI(NDIME, AZ, RBZ, N, BMAX, BMIN)
         BM = BMAX + BMIN
         BM = ABS(BM)
         AINF = AINF + (BMAX + BMIN)
         if (BM <= EPS) then
            exit
            !    elseif (JITER > 30) then
         elseif (JITER > 300) then
            IFAIL = 1
            exit
         end if
      end do
      RBINF = (BMAX - BMIN)/2.0d0
      ! write(700, *) AINF, RBINF
   end subroutine bpopt

   !> @brief Estimate band edges from scalar recursion coefficients.
   !> @details Legacy helper for terminator construction that scans a finite
   !>          continued-fraction coefficient set and returns approximate spectral
   !>          bounds.
   !> @param[inout] this Recursion object; no persistent state is modified.
   !> @param[in] nl Number of recursion levels in the coefficient arrays.
   !> @param[in] as Alpha coefficient sequence.
   !> @param[in] bs Beta coefficient sequence.
   !> @param[in] n Number of coefficients to use.
   !> @param[out] emax Estimated upper spectral edge.
   !> @param[out] emin Estimated lower spectral edge.
   module subroutine emami(this, nl, as, bs, n, emax, emin)
      class(recursion), intent(inout) :: this
      ! Input
      integer, intent(in) :: N, NL
      real(rp), dimension(NL), intent(in) :: AS, BS
      ! Output
      real(rp), intent(out) :: EMAX, EMIN
      ! Local Variables
      integer :: I, ISTOP, NUM
      real(rp) :: DELE, E, E1, E2, EMAX0, EMIN0, EPS, P, RELFEH, X1, X2
      real(rp), dimension(NL) :: A, B

      EMAX0 = -1.0d6
      EMIN0 = +1.0d6
      do I = 1, N
         A(I) = AS(I)
         B(I) = BS(I)
      end do
      B(1) = 0.0d0
      B(N + 1) = 0.0d0
      do I = 1, N
         X1 = A(I) + ABS(B(I)) + ABS(B(I + 1))
         X2 = A(I) - ABS(B(I)) - ABS(B(I + 1))
         if (EMAX0 <= X1) then
            EMAX0 = X1
         end if
         if (EMIN0 > X2) then
            EMIN0 = X2
         end if
      end do
      ! An identically zero scalar recursion channel has a zero-width spectrum.
      ! Its band edges are already known; avoid the relative 0/0 test below.
      if (EMAX0 == 0.0_rp .and. EMIN0 == 0.0_rp) then
         EMAX = 0.0_rp
         EMIN = 0.0_rp
         return
      end if
      RELFEH = 2.d0**(-39)
      EPS = 1.0d-6
      !C....CALCULATION OF EMAX
      ISTOP = 0
      EMAX = EMAX0
      EMIN = EMIN0
      do
         E = (EMAX + EMIN)/2.0d0
         ISTOP = ISTOP + 1
         if (ISTOP > 50) goto 1000
         NUM = 0
         P = A(1) - E
         if (P < 0.d0) then
            NUM = NUM + 1
         end if
         do I = 2, N
            if (P == 0.0) then
               P = (A(I) - E) - ABS(B(I))/RELFEH
               if (P < 0.d0) then
                  NUM = NUM + 1
               end if
            else
               P = (A(I) - E) - B(I)**2/P
               if (P < 0.d0) then
                  NUM = NUM + 1
               end if
            end if
         end do
         if (NUM == N) then
            EMAX = E
         end if
         if (NUM < N) then
            EMIN = E
         end if
         DELE = (EMAX - EMIN)/((EMAX + EMIN)/2.d0)
         DELE = ABS(DELE)
         if (DELE <= EPS) exit
      end do
      E1 = E
      !.....CALCULATION ON EMIN
      ISTOP = 0
      EMAX = E1
      EMIN = EMIN0
      do
         E = (EMAX + EMIN)/2.d0
         ISTOP = ISTOP + 1
         if (ISTOP > 50) goto 1000
         NUM = 0
         P = A(1) - E
         if (P < 0.d0) then
            NUM = NUM + 1
         end if
         do I = 2, N
            if (P == 0.d0) then
               P = (A(I) - E) - ABS(B(I))/RELFEH
               if (P < 0.d0) then
                  NUM = NUM + 1
               end if
            else
               P = (A(I) - E) - B(I)**2/P
               if (P < 0.d0) then
                  NUM = NUM + 1
               end if
            end if
         end do
         if (NUM == 0) then
            EMIN = E
         end if
         if (NUM > 0) then
            EMAX = E
         end if
         DELE = (EMAX - EMIN)/((EMAX + EMIN)/2.d0)
         DELE = ABS(DELE)
         if (DELE <= EPS) exit
      end do
      E2 = E
      !....
      EMAX = E1
      EMIN = E2
      return
1000  continue
      !1000 write (6, 10000)
      return
      !
      ! ... Format Declarations ...
      !
10000 format(" ", "NON-CONVERGE IN EMAMI")
   end subroutine emami

end submodule recursion_haydock
