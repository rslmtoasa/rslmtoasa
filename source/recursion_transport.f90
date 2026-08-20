submodule(recursion_mod) recursion_transport

   use iso_fortran_env, only: int64
   use kpm_profile_mod, only: g_kpm_profile

contains

   !> @brief Apply a real-space velocity-like operator to a Chebyshev block state.
   !> @details Multiplies psi_in by the supplied operator blocks v_op over the
   !>          current light cone. Conductivity and spin/orbital response moments
   !>          use this for v_a and v_b insertions when HOH is inactive.
   !> @param[inout] this Recursion object; updates idum to the propagated support.
   !> @param[in] c_or_n BLAS transpose selector for the operator blocks.
   !> @param[in] v_op Operator blocks in the real-space hopping layout.
   !> @param[in] psi_in Input block state, site-major as (nb,nb,nsite).
   !> @param[out] psi_out Output block state in the same layout.
   !> @note The local-Hamiltonian velocity path is intentionally not implemented.
   module subroutine velo_vec_matmul(this, c_or_n, v_op, psi_in, psi_out)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :, :), intent(in) :: v_op
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      character :: c_or_n 
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      ! Local variables
      integer :: ineigh, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(nb, nb) :: dum1, dum2, locham
      hblocksize = nb
      nat = this%lattice%kk

      psi_out(:, :, :) = (0.0d0, 0.0d0)

      ! NOT YET IMPLEMENTED THE VELOCITY OPERATOR FOR THE LOCAL HAMILTONIAN
      ! Write H*|phi_1> for the local Hamiltonian
      nlimplus1 = this%lattice%nmax + 1
      ! Write H*|phi_1> for the bulk Hamiltonian
      !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap)
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            !!!!! Probably useless because the onsite term is null. Here for consistency only
            call zgemm('n', 'n', nb, nb, nb, cone, v_op(:, :, 1, ih), nb, psi_in(:, :, k), nb, cone, psi_out(:, :, k), nb) 
            !!!!! Probably useless because the onsite term is null. Here for consistency only
         end if
         if (nr >= 2) then
            do ineigh = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm(c_or_n, 'n', nb, nb, nb, cone, v_op(:, :, ineigh, ih), nb, psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
      end do ! End loop in the clust
      !$omp end parallel do
   end subroutine velo_vec_matmul

   !> @brief Apply a velocity-like operator with the HOH two-sweep correction.
   !> @details Computes the response-operator action used by stochastic
   !>          conductivity moments when the orthogonalized Hamiltonian correction
   !>          is active, including the companion vo_op sweep through H psi.
   !> @param[inout] this Recursion object; uses psi1/psi2/hohpsi scratch arrays.
   !> @param[in] v_op Primary velocity or torque operator blocks.
   !> @param[in] vo_op Orthogonalization partner blocks for the HOH correction.
   !> @param[in] psi_in Input block state, site-major as (nb,nb,nsite).
   !> @param[out] psi_out Output block state in the same layout.
   !> @note Mutates izero/idum while propagating the active light cone.
   module subroutine velo_hoh_vec_matmul(this, v_op, vo_op, psi_in, psi_out)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      complex(rp), dimension(:, :, :, :), intent(in) :: v_op
      complex(rp), dimension(:, :, :, :), intent(in) :: vo_op
      ! Local variables
      integer :: ineigh, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(nb, nb) :: dum1, dum2, locham

      hblocksize = nb
      nat = this%lattice%kk

      this%hohpsi(:, :, :) = (0.0d0, 0.0d0)
      this%enupsi(:, :, :) = (0.0d0, 0.0d0)
      this%socpsi(:, :, :) = (0.0d0, 0.0d0)

      ! Write H*|phi_1> for the local Hamiltonian
      nlimplus1 = this%lattice%nmax + 1
      ! Write H*|phi_1> for the bulk Hamiltonian
      !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', nb, nb, nb, cone, v_op(:, :, 1, ih), nb, psi_in(:, :, k), nb, cone, this%psi2(:, :, k), nb)
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, 1, ih), nb, psi_in(:, :, k), nb, cone, this%psi1(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', nb, nb, nb, cone, v_op(:, :, ineigh, ih), nb, psi_in(:, :, nnmap), nb, cone, this%psi2(:, :, k), nb)
                  call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih), nb, psi_in(:, :, nnmap), nb, cone, this%psi1(:, :, k), nb)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
      end do ! End loop in the clust
      !$omp end parallel do

      ! Mapping update for hoh calculation
      this%izero(:) = this%idum(:)

      !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k) ! Atom type
         nr = this%lattice%nn(k, 1) ! Number of neighbours
         if (nr >= 2) then
            do ineigh = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', nb, nb, nb, cone, vo_op(:, :, ineigh, ih), nb, this%psi1(:, :, nnmap), nb, cone, this%hohpsi(:, :, k), nb)
                     this%idum(k) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do
      !$omp end parallel do

      ! H = h - hoh + e_nu + l.s
      psi_out(:, :, :) = this%psi2(:, :, :) - this%hohpsi(:, :, :) + this%enupsi(:, :, :) + this%socpsi(:, :, :)

      this%psi2(:, :, :) = (0.0_rp, 0.0_rp)
      this%psi1(:, :, :) = (0.0_rp, 0.0_rp)

   end subroutine velo_hoh_vec_matmul

   !> @brief Apply the scaled Hamiltonian with HOH correction to a block state.
   !> @details Forms ((H-b)/a) psi_in using the two-sweep h - H O H + e_nu + l.s
   !>          representation. Chebyshev stochastic transport uses this as the
   !>          recurrence operator when hamiltonian%hoh is enabled.
   !> @param[inout] this Recursion object; uses HOH scratch arrays and light-cone maps.
   !> @param[in] psi_in Input block state, site-major as (nb,nb,nsite).
   !> @param[out] psi_out Scaled Hamiltonian action in the same layout.
   !> @param[in] a Hamiltonian scaling half-width.
   !> @param[in] b Hamiltonian scaling center.
   !> @note Mutates izero/idum to the support reached by the Hamiltonian action.
   module subroutine ham_hoh_vec_matmul(this, psi_in, psi_out, a, b)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      real(rp), intent(in) :: a, b
      ! Local variables
      integer :: ineigh, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(nb, nb) :: dum1, dum2, locham

      hblocksize = nb
      nat = this%lattice%kk

      this%hohpsi(:, :, :) = (0.0d0, 0.0d0)
      this%enupsi(:, :, :) = (0.0d0, 0.0d0)
      this%socpsi(:, :, :) = (0.0d0, 0.0d0)

      ! Write H*|phi_1> for the local Hamiltonian
      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap,locham)
         do k = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1) ! Number of neighbours of atom i
            if (this%izero(k) /= 0) then
               locham = this%hamiltonian%hall(1:nb, 1:nb, 1, k)
               call zgemm('n', 'n', nb, nb, nb, cone, locham, nb, psi_in(:, :, k), nb, cone, this%psi2(:, :, k), nb)
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%lsham(:, :, ih), nb, psi_in(:, :, k), nb, cone, this%socpsi(:, :, k), nb)
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%enim(:, :, ih), nb, psi_in(:, :, k), nb, cone, this%enupsi(:, :, k), nb)
               if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallcc(:, :, 1, k), nb, psi_in(:, :, k), nb, cone, this%enupsi(:, :, k), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, ineigh, k), nb, psi_in(:, :, nnmap), nb, cone, this%psi2(:, :, k), nb)
                        if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallcc(:, :, ineigh, k), nb, psi_in(:, :, nnmap), nb, cone, this%enupsi(:, :, k), nb)
                        this%idum(k) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      ! Write H*|phi_1> for the bulk Hamiltonian
      !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            locham = this%hamiltonian%ee(1:nb, 1:nb, 1, ih)
            call zgemm('n', 'n', nb, nb, nb, cone, locham, nb, psi_in(:, :, k), nb, cone, this%psi2(:, :, k), nb)
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%enim(:, :, ih), nb, psi_in(:, :, k), nb, cone, this%enupsi(:, :, k), nb)
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%lsham(:, :, ih), nb, psi_in(:, :, k), nb, cone, this%socpsi(:, :, k), nb)
            if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eecc(:, :, 1, ih), nb, psi_in(:, :, k), nb, cone, this%enupsi(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih), nb, psi_in(:, :, nnmap), nb, cone, this%psi2(:, :, k), nb)
                  if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eecc(:, :, ineigh, ih), nb, psi_in(:, :, nnmap), nb, cone, this%enupsi(:, :, k), nb)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
      end do ! End loop in the clust
      !$omp end parallel do

      ! Mapping update for hoh calculation
      this%izero(:) = this%idum(:)

      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap,locham)
         do k = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1) ! Number of neighbours of atom i
            if (this%izero(k) /= 0) then
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallo(:, :, 1, k), nb, this%psi2(:, :, k), nb, cone, this%hohpsi(:, :, k), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallo(:, :, ineigh, k), nb, this%psi2(:, :, nnmap), nb, cone, this%hohpsi(:, :, k), nb)
                        this%idum(k) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k) ! Atom type
         nr = this%lattice%nn(k, 1) ! Number of neighbours
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eeo(:, :, 1, ih), nb, this%psi2(:, :, k), nb, cone, this%hohpsi(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eeo(:, :, ineigh, ih), nb, this%psi2(:, :, nnmap), nb, cone, this%hohpsi(:, :, k), nb)
                     this%idum(k) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do
      !$omp end parallel do

      ! H = h - hoh + e_nu + l.s
      psi_out(:, :, :) = this%psi2(:, :, :) - this%hohpsi(:, :, :) + this%enupsi(:, :, :) + this%socpsi(:, :, :)

      ! Do the scaling and shifting
      psi_out(:, :, :) = psi_out(:, :, :) - b*psi_in(:, :, :)
      psi_out(:, :, :) = psi_out(:, :, :)/a
     
      this%psi2(:, :, :) = (0.0_rp, 0.0_rp)
   end subroutine ham_hoh_vec_matmul

   !> @brief Apply the scaled no-HOH Hamiltonian to a block state.
   !> @details Forms ((H-b)/a) psi_in over the local and bulk real-space hopping
   !>          blocks. This is the legacy Chebyshev recurrence operator for
   !>          stochastic conductivity and related response moments.
   !> @param[inout] this Recursion object; updates idum to the propagated support.
   !> @param[in] psi_in Input block state, site-major as (nb,nb,nsite).
   !> @param[out] psi_out Scaled Hamiltonian action in the same layout.
   !> @param[in] a Hamiltonian scaling half-width.
   !> @param[in] b Hamiltonian scaling center.
   !> @note Includes ccor_2c blocks when hamiltonian%ccor_2c is set.
   module subroutine ham_vec_matmul(this, psi_in, psi_out, a, b)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      real(rp), intent(in) :: a, b
      ! Local variables
      integer :: ineigh, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(nb, nb) :: dum1, dum2, locham
      hblocksize = nb
      nat = this%lattice%kk

      psi_out(:, :, :) = (0.0d0, 0.0d0)

      ! Write H*|phi_1> for the local Hamiltonian
      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap,locham)
         do k = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1) ! Number of neighbours of atom i
            if (this%izero(k) /= 0) then
               locham = this%hamiltonian%hall(1:nb, 1:nb, 1, k) + this%hamiltonian%lsham(1:nb, 1:nb, ih)
               if (this%hamiltonian%ccor_2c) locham = locham + this%hamiltonian%hallcc(1:nb, 1:nb, 1, k)
               call zgemm('n', 'n', nb, nb, nb, cone, locham, nb, psi_in(:, :, k), nb, cone, psi_out(:, :, k), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        if (this%hamiltonian%ccor_2c) then
                           call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, ineigh, k) + this%hamiltonian%hallcc(:, :, ineigh, k), nb, psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
                        else
                           call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, ineigh, k), nb, psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
                        end if
                        this%idum(k) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
         !$omp end parallel do
      end if ! End of local Hamiltonian loop

      ! Write H*|phi_1> for the bulk Hamiltonian
      !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap,locham)
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            locham = this%hamiltonian%ee(1:nb, 1:nb, 1, ih) + this%hamiltonian%lsham(1:nb, 1:nb, ih)
            if (this%hamiltonian%ccor_2c) locham = locham + this%hamiltonian%eecc(1:nb, 1:nb, 1, ih)
            call zgemm('n', 'n', nb, nb, nb, cone, locham, nb, psi_in(:, :, k), nb, cone, psi_out(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  if (this%hamiltonian%ccor_2c) then
                     call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih) + this%hamiltonian%eecc(:, :, ineigh, ih), nb, psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
                  else
                     call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih), nb, psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
                  end if
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
      end do ! End loop in the clust
      !$omp end parallel do
      ! Do the scaling and shifting
      psi_out(:, :, :) = psi_out(:, :, :) - b*psi_in(:, :, :)
      psi_out(:, :, :) = psi_out(:, :, :)/a 
   end subroutine ham_vec_matmul

   !> @brief Compute two-index stochastic Chebyshev moments for transport.
   !> @details Builds mu_nm_stochastic = <r|T_m(H) v_a T_n(H) v_b|r> for charge,
   !>          spin, orbital, accumulation, and torque conductivity routes. The
   !>          starting states are either representative atom types or random
   !>          vectors, selected by control%cond_calctype.
   !> @param[inout] this Recursion object; allocates and fills mu_nm_stochastic.
   !> @note Builds the required real-space velocity/torque operators and may
   !>       dispatch to GPU or fast CPU stochastic kernels.
   module subroutine compute_moments_stochastic(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: ineigh, ih, i, j, k, nr, ll, m, n, l, hblocksize, nat, nnmap, loop_over, ie, lmax, ntype
      integer(int64) :: matrix_dimension, nnz, complex_bytes, integer_bytes, operator_h2d_bytes
      integer(int64) :: gpu_h2d_bytes, gpu_d2h_bytes
      complex(rp), dimension(nb, nb) :: dum, dum1, dum2
      complex(rp), dimension(:, :), allocatable :: S_op, L_op
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz
      complex(rp), dimension(:, :, :), allocatable :: psiref, w0, w1, w2, right_vec, v0, v1, v2, cn
      complex(rp), dimension(:, :, :, :), allocatable :: left_vec
      real(rp), dimension(this%en%channels_ldos + 10) :: w, wscale
      real(rp), dimension(this%control%cond_ll) :: kernel
      complex(rp), dimension(nb, nb, this%en%channels_ldos + 10) :: g0
      real(rp) :: a, b, rng, emin_win, emax_win
      real(rp) :: gpu_h2d_seconds, gpu_cheb_seconds, gpu_d2h_seconds
      complex(rp) :: exp_factor
      logical :: use_gpu
      character(len=48) :: precision_label

      lmax = lmax_basis
      hblocksize = nb
      nat = this%lattice%kk

      ! Memory allocation
#ifdef USE_SAFE_ALLOC
      select case(this%control%cond_calctype)
      case('per_type')
         call g_safe_alloc%allocate('recursion.mu_nm_stochastic', this%mu_nm_stochastic, (/2*(lmax + 1)**2, 2*(lmax + 1)**2, &
                                                                                       (this%lattice%control%cond_ll, &
                                                                               this%lattice%control%cond_ll, this%lattice%ntype/)
      case('random_vec')
         call g_safe_alloc%allocate('recursion.mu_nm_stochastic', this%mu_nm_stochastic, (/2*(lmax + 1)**2, 2*(lmax + 1)**2, &
                                                                                       (this%lattice%control%cond_ll, &
                                                                               this%lattice%control%cond_ll, this%control%random_vec_num/)
#else
      select case(this%control%cond_calctype)
      case('per_type')
         allocate (this%mu_nm_stochastic(2*(lmax + 1)**2, 2*(lmax + 1)**2, this%lattice%control%cond_ll, this%lattice%control%cond_ll,this%lattice%ntype))
      case('random_vec')
         allocate (this%mu_nm_stochastic(2*(lmax + 1)**2, 2*(lmax + 1)**2, this%lattice%control%cond_ll,this%lattice%control%cond_ll,this%control%random_vec_num))
      end select
#endif
      allocate(psiref(hblocksize, hblocksize, this%lattice%kk), left_vec(hblocksize, hblocksize, this%lattice%kk, this%control%cond_ll))
      allocate(w0(hblocksize, hblocksize, this%lattice%kk), w1(hblocksize, hblocksize, this%lattice%kk), right_vec(hblocksize, hblocksize, this%lattice%kk))
      allocate(w2(hblocksize, hblocksize, this%lattice%kk), v0(hblocksize, hblocksize, this%lattice%kk), v1(hblocksize, hblocksize, this%lattice%kk))
      allocate(v2(hblocksize, hblocksize, this%lattice%kk), S_op(hblocksize, hblocksize), L_op(hblocksize, hblocksize))

      ! General procedures
      call resolve_chebyshev_window(this, emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp

      ! Decide the backend once for the complete transport moment request.
      ! Repeating the readiness checks for every trace both obscures the
      ! profile and can produce misleading fallback diagnostics.
      use_gpu = gpu_plugin_ready(this, 'compute_moments_stochastic()', allow_hoh=.true.)
      if (use_gpu .and. trim(this%control%gpu_precision) == 'fp32' .and. &
          (trim(this%control%gpu_backend) == 'fft' .or. trim(this%control%gpu_backend) == 'conv')) then
         call g_logger%fatal('gpu_precision=fp32 is not supported by the structured CUDA stochastic transport backend; '// &
            'select gpu_backend=csr/bsr or gpu_precision=fp64.', __FILE__, __LINE__)
      end if

      select case(this%control%cond_calctype)
      case('per_type')
         loop_over = this%lattice%ntype
      case('random_vec')
         loop_over = this%control%random_vec_num
      end select

      matrix_dimension = int(nb, int64) * int(nat, int64)
      if (size(this%lattice%nn, 2) > 1) then
         nnz = int(count(this%lattice%nn(:, 2:) > 0), int64) * int(nb, int64) * int(nb, int64)
      else
         nnz = 0_int64
      end if
      complex_bytes = int(storage_size(this%hamiltonian%ee) / 8, int64)
      integer_bytes = int(storage_size(this%lattice%nn) / 8, int64)
      operator_h2d_bytes = 0_int64
      if (use_gpu) then
         operator_h2d_bytes = int(size(this%hamiltonian%ee), int64) * complex_bytes
         if (this%lattice%nmax > 0) operator_h2d_bytes = operator_h2d_bytes + &
            int(size(this%hamiltonian%hall), int64) * complex_bytes
         operator_h2d_bytes = operator_h2d_bytes + int(size(this%hamiltonian%lsham), int64) * complex_bytes + &
            int(size(this%lattice%nn), int64) * integer_bytes + int(size(this%lattice%iz), int64) * integer_bytes
         if (this%hamiltonian%hoh) then
            operator_h2d_bytes = operator_h2d_bytes + int(size(this%hamiltonian%eeo), int64) * complex_bytes + &
               int(size(this%hamiltonian%enim), int64) * complex_bytes
            if (this%lattice%nmax > 0) operator_h2d_bytes = operator_h2d_bytes + &
               int(size(this%hamiltonian%hallo), int64) * complex_bytes
         end if
         operator_h2d_bytes = operator_h2d_bytes + int(size(this%hamiltonian%v_a), int64) * complex_bytes + &
            int(size(this%hamiltonian%v_b), int64) * complex_bytes
         if (this%hamiltonian%hoh) operator_h2d_bytes = operator_h2d_bytes + &
            int(size(this%hamiltonian%vo_a), int64) * complex_bytes + int(size(this%hamiltonian%vo_b), int64) * complex_bytes
         if (trim(this%control%gpu_precision) == 'fp64') then
            precision_label = 'cuda_fp64'
         else
            precision_label = 'cuda_fp32_moments_fp64_host'
         end if
      else if (trim(this%control%cheb_backend) == 'legacy') then
         precision_label = 'cpu_fp64'
      else if (trim(this%control%cheb_backend) == 'fast_dp') then
         precision_label = 'cpu_fp64_fast'
      else
         precision_label = 'cpu_fp32_moments_fp64_host'
      end if
      call g_kpm_profile%configure(merge('cuda', 'cpu ', use_gpu), precision_label, &
         trim(this%control%cond_calctype), matrix_dimension, nnz, this%control%cond_ll, &
         this%control%lld, loop_over)

      call g_kpm_profile%start('T_operator')
      call this%hamiltonian%build_realspace_velocity_operators()

      ! Check the type of conductivity
      select case(this%control%cond_type)
      case ('charge')
         ! Redundant, but just testing
         this%hamiltonian%v_a(:, :, :, :) = this%hamiltonian%v_a(:, :, :, :)
         this%hamiltonian%vo_a(:, :, :, :) = this%hamiltonian%vo_a(:, :, :, :)
      case('spin')
         call this%hamiltonian%build_realspace_spin_operators()
         this%hamiltonian%v_a(:, :, :, :) = this%hamiltonian%js_a(:, :, :, :) 
         this%hamiltonian%vo_a(:, :, :, :) = this%hamiltonian%jso_a(:, :, :, :)
      case('orbital')
         call this%hamiltonian%build_realspace_orbital_velocity_operators()
         this%hamiltonian%v_a(:, :, :, :) = this%hamiltonian%jl_a(:, :, :, :)
         this%hamiltonian%vo_a(:, :, :, :) = this%hamiltonian%jlo_a(:, :, :, :)
      case('spin_accumulation')
         select case(this%hamiltonian%js_alpha)
         case('z')
            S_op = S_z
         case('x')
            S_op = S_x
         case('y')
            S_op = S_y
         end select
         this%hamiltonian%v_a(:, :, :, :) = (0.d0, 0.0d0)
         this%hamiltonian%vo_a(:, :, :, :) = (0.d0, 0.0d0)
         do ntype = 1, this%lattice%ntype
            this%hamiltonian%v_a(:, :, 1, ntype) = S_op(:, :)
         end do
      case('spin_torque')
         call this%hamiltonian%build_realspace_spin_torque_operators()
         this%hamiltonian%v_a(:, :, :, :) = this%hamiltonian%js_a(:, :, :, :)
         this%hamiltonian%vo_a(:, :, :, :) = this%hamiltonian%jso_a(:, :, :, :)
      case('orbital_torque')
         call this%hamiltonian%build_realspace_orbital_torque_operators()
         this%hamiltonian%v_a(:, :, :, :) = this%hamiltonian%jl_a(:, :, :, :)
         this%hamiltonian%vo_a(:, :, :, :) = this%hamiltonian%jlo_a(:, :, :, :)
      case('orbital_accumulation')
         !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
         mLx(:, :) = L_x(:, :)
         mLy(:, :) = L_y(:, :)
         mLz(:, :) = L_z(:, :)
   
         ! Transforming them into the spherical harmonics coordinates
         call hcpx(mLx, 'cart2sph')
         call hcpx(mLy, 'cart2sph')
         call hcpx(mLz, 'cart2sph') 
            
         ! Pick which orbital operator L_x, L_y, or L_z based on some user choice
         select case (this%hamiltonian%jl_alpha)   ! or whichever variable holds 'x','y','z'
         case ('x')
            L_op(1:norb, 1:norb) = mLx(:, :)
            L_op(norb+1:nb, norb+1:nb) = mLx(:, :)
         case ('y') 
            L_op(1:norb, 1:norb) = mLy(:, :)
            L_op(norb+1:nb, norb+1:nb) = mLy(:, :)
         case ('z')
            L_op(1:norb, 1:norb) = mLz(:, :)
            L_op(norb+1:nb, norb+1:nb) = mLz(:, :)
         end select
         this%hamiltonian%v_a(:, :, :, :) = (0.d0, 0.0d0)
         this%hamiltonian%vo_a(:, :, :, :) = (0.d0, 0.0d0)
         do ntype = 1, this%lattice%ntype
            this%hamiltonian%v_a(:, :, 1, ntype) = L_op(:, :)
         end do
      end select
      call g_kpm_profile%stop('T_operator')

      ! Check what kind of calculation
      do i = 1, loop_over
         call g_kpm_profile%start('T_trace_setup')
         call random_seed()

         ! Initializing wave functions
         v0(:, :, :) = (0.0d0, 0.0d0)
         v1(:, :, :) = (0.0d0, 0.0d0)
         v2(:, :, :) = (0.0d0, 0.0d0)
         psiref(:, :, :) = (0.0d0, 0.0d0)
         w0(:, :, :) = (0.0d0, 0.0d0)
         w1(:, :, :) = (0.0d0, 0.0d0)
         w2(:, :, :) = (0.0d0, 0.0d0)
         right_vec(:, :, :) = (0.0d0, 0.0d0)
         left_vec(:, :, :, :) = (0.0d0, 0.0d0)
         dum(:, :) = (0.0d0, 0.0d0)

         select case(this%control%cond_calctype)
         case('per_type')
            j =  this%lattice%atlist(i)
            call g_logger%info('Chebyshev moments being calculated taking atom type '//int2str(j), __FILE__, __LINE__)
            ! Initializing neighbording map
            this%izero(:) = 0
            this%izero(j) = 1
           ! Initializing psi
            do m = 1, nb
               psiref(m, m, j) = (1.0d0, 0.0d0)
            end do
         case('random_vec')
            call g_logger%info('Chebyshev moments being calculated for random vector '//int2str(i), __FILE__, __LINE__)
            this%izero(:) = 1
            ! Initialize random vector
            do k = 1, this%lattice%kk
               call random_number(rng)
               do m = 1, nb
                  psiref(m, m, k) =  (exp(2.0_rp * pi * i_unit * (rng))) !(2.0d0 * rng - 1.0d0)*sqrt(3.0d0) !(exp(2.0_rp * pi * i_unit * (rng)))
               end do
            end do
            ! Normalize the full matrix 
            psiref(:, :, :) = psiref(:, :, :) / sqrt(real(this%lattice%kk))
         end select
         call g_kpm_profile%stop('T_trace_setup')

         if (use_gpu) then
            if (i == 1) then
               call g_kpm_profile%start('T_H2D')
               call gpu_plugin_upload_hamiltonian(this)
               if (this%hamiltonian%hoh) then
                  call this%gpu_backend%set_velocity(this%hamiltonian%v_a, this%hamiltonian%v_b, &
                     vo_a=this%hamiltonian%vo_a, vo_b=this%hamiltonian%vo_b)
               else
                  call this%gpu_backend%set_velocity(this%hamiltonian%v_a, this%hamiltonian%v_b)
               end if
               call this%gpu_backend%set_precision(merge(1, 0, trim(this%control%gpu_precision) == 'fp64'))
               call g_kpm_profile%stop('T_H2D')
               call g_kpm_profile%add_bytes('H2D', operator_h2d_bytes)
            end if
            call this%gpu_backend%stochastic_moments(psiref, this%control%cond_ll, a, b, &
               this%mu_nm_stochastic(:, :, :, :, i))
            call this%gpu_backend%stochastic_profile(gpu_h2d_seconds, gpu_cheb_seconds, gpu_d2h_seconds, &
               gpu_h2d_bytes, gpu_d2h_bytes)
            call g_kpm_profile%add_seconds('T_H2D', gpu_h2d_seconds)
            call g_kpm_profile%add_seconds('T_cheb_moments', gpu_cheb_seconds)
            call g_kpm_profile%add_seconds('T_D2H', gpu_d2h_seconds)
            call g_kpm_profile%add_bytes('H2D', gpu_h2d_bytes)
            call g_kpm_profile%add_bytes('D2H', gpu_d2h_bytes)
            cycle
         end if

         ! Fast CPU stochastic kernel (single precision, hoh-aware). Falls back
         ! to the legacy loop below only when cheb_backend='legacy'.
         if (this%hamiltonian%ccor_2c .and. this%hamiltonian%hoh .and. trim(this%control%cheb_backend) /= 'legacy') then
            call g_logger%warning('ccor_2c with hoh is not supported by fast stochastic Chebyshev; falling back to legacy.', __FILE__, __LINE__)
         end if
         if (trim(this%control%cheb_backend) /= 'legacy' .and. &
             .not. (this%hamiltonian%ccor_2c .and. this%hamiltonian%hoh)) then
            if (this%hamiltonian%hoh) then
               call g_kpm_profile%start('T_cheb_moments')
               call cheb_moments_stochastic_fast(psiref, this%hamiltonian%ee, &
                  this%hamiltonian%hall, this%hamiltonian%lsham, this%lattice%nn, &
                  this%lattice%iz, this%lattice%kk, nb, size(this%lattice%nn, 2), &
                  this%lattice%ntype, this%lattice%nmax, this%control%cond_ll, a, b, &
                  this%mu_nm_stochastic(:, :, :, :, i), .true., this%hamiltonian%eeo, &
                  this%hamiltonian%hallo, this%hamiltonian%enim, this%hamiltonian%v_a, &
                  this%hamiltonian%v_b, this%hamiltonian%vo_a, this%hamiltonian%vo_b)
               call g_kpm_profile%stop('T_cheb_moments')
            else if (this%hamiltonian%ccor_2c) then
               call ensure_ccor_operator_blocks(this)
               call g_kpm_profile%start('T_cheb_moments')
               call cheb_moments_stochastic_fast(psiref, this%ee_ccor_work, &
                  this%hall_ccor_work, this%hamiltonian%lsham, this%lattice%nn, &
                  this%lattice%iz, this%lattice%kk, nb, size(this%lattice%nn, 2), &
                  this%lattice%ntype, this%lattice%nmax, this%control%cond_ll, a, b, &
                  this%mu_nm_stochastic(:, :, :, :, i), .false., this%ee_ccor_work, &
                  this%hall_ccor_work, this%hamiltonian%lsham, this%hamiltonian%v_a, &
                  this%hamiltonian%v_b, this%hamiltonian%v_a, this%hamiltonian%v_b)
               call g_kpm_profile%stop('T_cheb_moments')
            else
               call g_kpm_profile%start('T_cheb_moments')
               call cheb_moments_stochastic_fast(psiref, this%hamiltonian%ee, &
                  this%hamiltonian%hall, this%hamiltonian%lsham, this%lattice%nn, &
                  this%lattice%iz, this%lattice%kk, nb, size(this%lattice%nn, 2), &
                  this%lattice%ntype, this%lattice%nmax, this%control%cond_ll, a, b, &
                  this%mu_nm_stochastic(:, :, :, :, i), .false., this%hamiltonian%ee, &
                  this%hamiltonian%hall, this%hamiltonian%lsham, this%hamiltonian%v_a, &
                  this%hamiltonian%v_b, this%hamiltonian%v_a, this%hamiltonian%v_b)
               call g_kpm_profile%stop('T_cheb_moments')
            end if
            cycle
         end if

         ! Computing the left vector <r|Tm(H)
         call g_kpm_profile%start('T_cheb_moments')
         do m=1, this%control%cond_ll 
            if (m == 1) then
               w1(:, :, :) = psiref(:, :, :)
            else if (m == 2) then
               w0(:, :, :) = w1(:, :, :)
               if (this%hamiltonian%hoh) then
                  call this%ham_hoh_vec_matmul(w0, w1, a, b)
               else
                  call this%ham_vec_matmul(w0, w1, a, b)
               end if
               this%izero(:) = this%idum(:)
            else if (m > 2) then
               if (this%hamiltonian%hoh) then
                  call this%ham_hoh_vec_matmul(w1, w2, a, b)
               else
                  call this%ham_vec_matmul(w1, w2, a, b)
               end if
               this%izero(:) = this%idum(:)
               w2(:, :, :) = 2 * w2(:, :, :) - w0(:, :, :)
               w0(:, :, :) = w1(:, :, :)
               w1(:, :, :) = w2(:, :, :)
               w2(:, :, :) = (0.0d0, 0.0d0)
            end if
            left_vec(:, :, :, m) = w1(:, :, :) 
         end do

         ! Redifining neighboring map for the right vector v_a*Tn(H)*v_b|r>
         select case(this%control%cond_calctype)
         case('per_type')
            this%izero(:) = 0
            this%izero(j) = 1
         case('random_vec')
           ! this%izero(:) = 1
         end select

         ! Computing the right vector v_a*Tn(H)*v_b|r>
         ! Multiply with the velocity operator v_b
         if (this%hamiltonian%hoh) then
            call this%velo_hoh_vec_matmul(this%hamiltonian%v_b, this%hamiltonian%vo_b, psiref, v0)
         else
            call this%velo_vec_matmul('n',this%hamiltonian%v_b, psiref, v0)
         end if
         this%izero(:) = this%idum(:) 
         do n=1, this%control%cond_ll
            call show_progress(n, this%control%cond_ll)
            if (n == 1) then
               v1(:, :, :) = v0(:, :, :) 
            else if (n == 2) then
               v0(:, :, :) = v1(:, :, :)
               if (this%hamiltonian%hoh) then
                  call this%ham_hoh_vec_matmul(v0, v1, a, b)
               else
                  call this%ham_vec_matmul(v0, v1, a, b)
               end if
               this%izero(:) = this%idum(:)
            else if (n > 2) then
               if (this%hamiltonian%hoh) then
                  call this%ham_hoh_vec_matmul(v1, v2, a, b)
               else
                  call this%ham_vec_matmul(v1, v2, a, b)
               end if
               this%izero(:) = this%idum(:)
               v2(:, :, :) = 2 * v2(:, :, :) - v0(:, :, :)
               v0(:, :, :) = v1(:, :, :)
               v1(:, :, :) = v2(:, :, :)
               v2(:, :, :) = (0.0d0, 0.0d0)
            end if
            ! Multiply with the velocity operator v_a
            if (this%hamiltonian%hoh) then
               call this%velo_hoh_vec_matmul(this%hamiltonian%v_a, this%hamiltonian%vo_a, v1, right_vec)
            else
               call this%velo_vec_matmul('n',this%hamiltonian%v_a, v1, right_vec)
            end if
            this%izero(:) = this%idum(:)
            !$omp parallel do default(shared) private(m, dum, k) schedule(dynamic)
            do m=1, this%control%cond_ll
               dum(:, :) = (0.0d0, 0.0d0)
               do k=1, this%lattice%kk
                  call zgemm('c', 'n', nb, nb, nb, cone, left_vec(:, :, k, m), nb, right_vec(:, :, k), nb, cone, dum(:, :), nb)
               end do
               this%mu_nm_stochastic(:, :, n, m, i) = dum(:, :)
            end do
            !$omp end parallel do
         end do
      end do

      deallocate(psiref, w0, w1, w2, right_vec, v0, v1, v2, left_vec, S_op, L_op)

   end subroutine compute_moments_stochastic

   !> @brief Generate block-Lanczos coefficients for intersite atom pairs.
   !> @details Runs four phase-combination recursions per pair so downstream Green
   !>          function code can reconstruct G_ij for exchange and conductivity.
   !>          MPI ranks split lattice%ijpair, storing local results in a_b/b2_b.
   !> @param[inout] this Recursion object; fills pair-local block coefficients.
   !> @note May dispatch to CUDA or haydock_fast before falling back to legacy loops.
   module subroutine recur_b_ij(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, ij, j, l, ll, kk, m, reci
      integer :: llmax ! Recursion steps
      complex(rp) :: asign, bsign

      integer :: ij_loc

      llmax = this%lattice%control%lld
      if (gpu_plugin_ready(this, 'recur_b_ij()', allow_hoh=.true.)) then
         call gpu_plugin_upload_hamiltonian(this)
         do ij = start_atom, end_atom
            ij_loc = g2l_map(ij)
            i = this%lattice%ijpair(ij, 1)
            j = this%lattice%ijpair(ij, 2)
            call g_logger%info('CUDA block recursion on progress between atoms '// &
               int2str(i)//' and '//int2str(j), __FILE__, __LINE__)
            do reci = 1, 4
               this%psi_b(:, :, :) = (0.0d0, 0.0d0)
               this%atemp_b(:, :, :) = (0.0d0, 0.0d0)
               this%b2temp_b(:, :, :) = (0.0d0, 0.0d0)

               select case (reci)
               case (1)
                  asign = (1.0d0, 0.0d0)*one_over_sqrt_two
                  bsign = (1.0d0, 0.0d0)*one_over_sqrt_two
               case (2)
                  asign = (1.0d0, 0.0d0)*one_over_sqrt_two
                  bsign = (-1.0d0, 0.0d0)*one_over_sqrt_two
               case (3)
                  asign = (1.0d0, 0.0d0)*one_over_sqrt_two
                  bsign = (0.0d0, 1.0d0)*one_over_sqrt_two
               case (4)
                  asign = (1.0d0, 0.0d0)*one_over_sqrt_two
                  bsign = (0.0d0, -1.0d0)*one_over_sqrt_two
               end select

               if ((reci .eq. 1) .and. (i .eq. j)) then
                  asign = (1.0d0, 0.0d0)
                  bsign = (1.0d0, 0.0d0)
               else if (i .eq. j) then
                  cycle
               end if

               do l = 1, nb
                  this%psi_b(l, l, i) = asign
                  this%psi_b(l, l, j) = bsign
               end do

               call this%gpu_backend%block_lanczos(this%psi_b, llmax, &
                  this%atemp_b, this%b2temp_b, &
                  prec=merge(1, 0, trim(this%control%cheb_backend) == 'fast_dp'))

               do ll = 1, llmax
                  do l = 1, nb
                     do m = 1, nb
                        this%a_b(l, m, ll, ij_loc*4 - 4 + reci) = this%atemp_b(l, m, ll)
                        this%b2_b(l, m, ll, ij_loc*4 - 4 + reci) = this%b2temp_b(l, m, ll)
                     end do
                  end do
               end do
            end do
         end do
         return
      end if

      ! Fast CPU block-Lanczos kernels (haydock_fast) for the inter-site (ij)
      ! moments. cheb_backend: 'fast' -> fp32, 'fast_dp' -> fp64.
      if (this%hamiltonian%ccor_2c .and. this%hamiltonian%hoh .and. trim(this%control%cheb_backend) /= 'legacy') then
         call g_logger%warning('ccor_2c with hoh is not supported by fast block-Lanczos; falling back to legacy.', __FILE__, __LINE__)
      end if
      if (trim(this%control%cheb_backend) /= 'legacy' .and. .not. this%hamiltonian%local_axis .and. &
          .not. (this%hamiltonian%ccor_2c .and. this%hamiltonian%hoh)) then
         if (this%hamiltonian%ccor_2c) call ensure_ccor_operator_blocks(this)
         do ij = start_atom, end_atom
            ij_loc = g2l_map(ij)
            i = this%lattice%ijpair(ij, 1)
            j = this%lattice%ijpair(ij, 2)
            call g_logger%info('Fast block recursion on progress between atoms '// &
               int2str(i)//' and '//int2str(j), __FILE__, __LINE__)
            do reci = 1, 4
               this%psi_b(:, :, :) = (0.0d0, 0.0d0)
               this%atemp_b(:, :, :) = (0.0d0, 0.0d0)
               this%b2temp_b(:, :, :) = (0.0d0, 0.0d0)
               select case (reci)
               case (1)
                  asign = (1.0d0, 0.0d0)*one_over_sqrt_two
                  bsign = (1.0d0, 0.0d0)*one_over_sqrt_two
               case (2)
                  asign = (1.0d0, 0.0d0)*one_over_sqrt_two
                  bsign = (-1.0d0, 0.0d0)*one_over_sqrt_two
               case (3)
                  asign = (1.0d0, 0.0d0)*one_over_sqrt_two
                  bsign = (0.0d0, 1.0d0)*one_over_sqrt_two
               case (4)
                  asign = (1.0d0, 0.0d0)*one_over_sqrt_two
                  bsign = (0.0d0, -1.0d0)*one_over_sqrt_two
               end select
               if ((reci .eq. 1) .and. (i .eq. j)) then
                  asign = (1.0d0, 0.0d0)
                  bsign = (1.0d0, 0.0d0)
               else if (i .eq. j) then
                  cycle
               end if
               do l = 1, nb
                  this%psi_b(l, l, i) = asign
                  this%psi_b(l, l, j) = bsign
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
                        this%a_b(l, m, ll, ij_loc*4 - 4 + reci) = this%atemp_b(l, m, ll)
                        this%b2_b(l, m, ll, ij_loc*4 - 4 + reci) = this%b2temp_b(l, m, ll)
                     end do
                  end do
               end do
            end do
         end do
         return
      end if

      !do ij=1, this%lattice%njij ! Loop on the number of pair of atoms
      do ij = start_atom, end_atom
         ij_loc = g2l_map(ij)
         i = this%lattice%ijpair(ij, 1) ! Atom number in the clust file, atom i
         j = this%lattice%ijpair(ij, 2) ! Atom number in the clust file, atom j
         call g_logger%info('Block recursion on progress between atoms '//int2str(i)//' and '//int2str(j), __FILE__, __LINE__)
         do reci = 1, 4
            ! Clear list of atoms in hopping region
            this%izero(:) = 0
            ! Initializing wave functions
            this%psi_b(:, :, :) = (0.0d0, 0.0d0)
            this%pmn_b(:, :, :) = (0.0d0, 0.0d0)

            select case (reci)
            case (1)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (1.0d0, 0.0d0)*one_over_sqrt_two
               this%izero(i) = 1
               this%izero(j) = 1
            case (2)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (-1.0d0, 0.0d0)*one_over_sqrt_two
               this%izero(i) = 1
               this%izero(j) = 1
            case (3)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (0.0d0, 1.0d0)*one_over_sqrt_two
               this%izero(i) = 1
               this%izero(j) = 1
            case (4)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (0.0d0, -1.0d0)*one_over_sqrt_two
               this%izero(i) = 1
               this%izero(j) = 1
            end select

            if ((reci .eq. 1) .and. (i .eq. j)) then
               asign = (1.0d0, 0.0d0)
               bsign = (1.0d0, 0.0d0)
            else if (i .eq. j) then
               cycle
            end if

            do l = 1, nb
               this%psi_b(l, l, i) = asign
               this%psi_b(l, l, j) = bsign
               this%atemp_b(l, l, llmax) = (0.0d0, 0.0d0)
               this%b2temp_b(l, l, 1) = (1.0d0, 0.0d0)
            end do

            call this%crecal_b()

            do ll = 1, llmax
               do l = 1, nb
                  do m = 1, nb
                     this%a_b(l, m, ll, ij_loc*4 - 4 + reci) = (this%atemp_b(l, m, ll))
                     this%b2_b(l, m, ll, ij_loc*4 - 4 + reci) = (this%b2temp_b(l, m, ll))
                  end do
               end do
            end do
         end do ! End of the loop on the basis
      end do ! End of the loop on the atom pairs
      ! For debug purposes
      !do i=1, this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
      !  do l = 1, nb  ! Loop on the orbital l
      !    write(122, *) ´orbital´, l
      !    do ll=1, llmax ! Loop on the recursion steps
      !      write(122, *) this%a(ll, l, i,1), this%b2(ll,l,i,1)
      !    end do
      !  end do
      !end do
   end subroutine recur_b_ij

   !> @brief Generate Chebyshev moments for intersite atom pairs.
   !> @details Runs four phase-combination moment recursions per ij pair so
   !>          green%chebyshev_green_ij can reconstruct intersite Green functions
   !>          for exchange and conductivity workflows.
   !> @param[inout] this Recursion object; fills pair-local slices of mu_n.
   !> @note Work is MPI-partitioned and may dispatch to CUDA.
   module subroutine chebyshev_recur_ij(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      real(rp) :: a, b, emin_win, emax_win
      integer :: i, ij, j, l, ll, kk, m, reci
      integer :: llmax ! Recursion steps
      complex(rp) :: asign, bsign

      integer :: ij_loc

      call resolve_chebyshev_window(this, emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp
      call cheb_fast_reset_cache()

      llmax = this%lattice%control%lld
      if (gpu_plugin_ready(this, 'chebyshev_recur_ij()', allow_hoh=.true.)) then
         call gpu_plugin_upload_hamiltonian(this)
         do ij = start_atom, end_atom
            ij_loc = g2l_map(ij)
            i = this%lattice%ijpair(ij, 1)
            j = this%lattice%ijpair(ij, 2)
            call g_logger%info(int2str(rank)//': CUDA Chebyshev recursion on progress between atoms '// &
               int2str(i)//' and '//int2str(j), __FILE__, __LINE__)
            do reci = 1, 4
               this%psi0(:, :, :) = (0.0d0, 0.0d0)
               select case (reci)
               case (1)
                  asign = (1.0d0, 0.0d0)*one_over_sqrt_two
                  bsign = (1.0d0, 0.0d0)*one_over_sqrt_two
               case (2)
                  asign = (1.0d0, 0.0d0)*one_over_sqrt_two
                  bsign = (-1.0d0, 0.0d0)*one_over_sqrt_two
               case (3)
                  asign = (1.0d0, 0.0d0)*one_over_sqrt_two
                  bsign = (0.0d0, 1.0d0)*one_over_sqrt_two
               case (4)
                  asign = (1.0d0, 0.0d0)*one_over_sqrt_two
                  bsign = (0.0d0, -1.0d0)*one_over_sqrt_two
               end select
               ! On-site self-pair (i==j): both writes below target atom i, so a
               ! 1/sqrt2 start would leave norm^2=1/2 and halve G_ii (moments are
               ! quadratic in psi0). Use unit weight for reci=1 and skip the other
               ! phases -- mirrors the block path recur_b_ij.
               if ((reci .eq. 1) .and. (i .eq. j)) then
                  asign = (1.0d0, 0.0d0)
                  bsign = (1.0d0, 0.0d0)
               else if (i .eq. j) then
                  cycle
               end if
               do l = 1, nb
                  this%psi0(l, l, i) = asign
                  this%psi0(l, l, j) = bsign
               end do
               call this%gpu_backend%chebyshev_moments(this%psi0, llmax, a, b, &
                  this%mu_n(:, :, :, ij_loc*4 - 4 + reci))
               if (real(sum(this%mu_n(:, :, llmax + 2, ij_loc*4 - 4 + reci))) > 1000.d0) then
                  call g_logger%fatal('Chebyshev moments did not converge. Check energy limits energy_min and energy_max', __FILE__, __LINE__)
               end if
            end do
         end do
         return
      end if

      !do ij=1, this%lattice%njij ! Loop on the number of pair of atoms
      do ij = start_atom, end_atom
         ij_loc = g2l_map(ij)
         i = this%lattice%ijpair(ij, 1) ! Atom number in the clust file, atom i
         j = this%lattice%ijpair(ij, 2) ! Atom number in the clust file, atom j
!         call g_logger%info('Chebyshev recursion on progress between atoms '//int2str(i)//' and '//int2str(j), __FILE__, __LINE__)
         call g_logger%info(int2str(rank)//': Chebyshev recursion on progress between atoms '//int2str(i)//' and '//int2str(j), __FILE__, __LINE__)
         do reci = 1, 4
            ! Starting state |phi_0>: the four +/-, +/-i sign combinations on
            ! the two atoms i and j of the pair (RMP 78, 275). The same block
            ! also serves as the projection reference inside the drivers.
            this%psi0(:, :, :) = (0.0d0, 0.0d0)

            select case (reci)
            case (1)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (1.0d0, 0.0d0)*one_over_sqrt_two
            case (2)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (-1.0d0, 0.0d0)*one_over_sqrt_two
            case (3)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (0.0d0, 1.0d0)*one_over_sqrt_two
            case (4)
               asign = (1.0d0, 0.0d0)*one_over_sqrt_two
               bsign = (0.0d0, -1.0d0)*one_over_sqrt_two
            end select

            ! On-site self-pair (i==j): both writes below target atom i, so a
            ! 1/sqrt2 start would leave norm^2=1/2 and halve G_ii (moments are
            ! quadratic in psi0). Use unit weight for reci=1 and skip the other
            ! phases -- mirrors the block path recur_b_ij.
            if ((reci .eq. 1) .and. (i .eq. j)) then
               asign = (1.0d0, 0.0d0)
               bsign = (1.0d0, 0.0d0)
            else if (i .eq. j) then
               cycle
            end if

            do l = 1, nb
               this%psi0(l, l, i) = asign
               this%psi0(l, l, j) = bsign
            end do

            ! Dispatch to the selected CPU backend (fast kernels or, for hoh /
            ! cheb_backend='legacy', the wrapped legacy recursion).
            call cheb_moments_cpu(this, this%psi0, llmax, a, b, &
               this%mu_n(:, :, :, ij_loc*4 - 4 + reci))

            if (real(sum(this%mu_n(:, :, llmax + 2, ij_loc*4 - 4 + reci))) > 1000.d0) then
               call g_logger%fatal('Chebyshev moments did not converge. Check energy limits energy_min and energy_max', __FILE__, __LINE__)
            end if
         end do
      end do
   end subroutine chebyshev_recur_ij

   !> @brief Compute Chebyshev orbital-moment response data.
   !> @details Builds angular-momentum operator insertions and evaluates
   !>          Chebyshev moment/Green-function contractions for the
   !>          post_processing_orbital_modern workflow. The route is always
   !>          Chebyshev and loops over cluster atoms rather than only nrec atoms.
   !> @param[inout] this Recursion object; allocates temporary orbital moments
   !>                and uses mu_n-style Chebyshev work arrays.
   !> @note MPI ranks split the atom loop; the physically meaningful orbital data
   !>       are accumulated in the orbital-moment arrays built inside this routine.
   module subroutine chebyshev_orbital_mod(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: atom_neighbor, ih, i, j, k, nr, ll, m, n, l, hblocksize, ntype, nnmap, nlimplus1, llcheb, lmax, nv, ie, random, iz, ia
      complex(rp), dimension(nb, nb) :: dum, dum1, dum2
      complex(rp), dimension(nb, this%lattice%kk) :: v
      complex(rp), dimension(:, :, :, :), allocatable :: hcheb
      complex(rp), dimension(:, :, :), allocatable :: psiref, w0, w1, w2, right_vec, v0, v1, v2, left_vec, left_vec1, left_vec2
      complex(rp), dimension(:, :, :), allocatable :: mu_n_orb, mu_tmp
      complex(rp), dimension(nb, nb, this%en%channels_ldos + 10) :: g0
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz
      complex(rp), dimension(nb, nb) :: L_op
      complex(rp) :: exp_factor
      real(rp) :: lz
      real(rp), dimension(this%control%lld) :: kernel
      real(rp), dimension(this%en%channels_ldos + 10) :: w, wscale, lzi
      real(rp), dimension(3) :: rij, crossrij
      
      real(rp) :: a, b, start, finish, rng, emin_win, emax_win
      ! External functions
      complex(rp), external :: zdotc

      integer :: i_loc
      real(rp) :: maxg_tmp
      integer :: ie_tmp

      hblocksize = nb
      nlimplus1 = this%lattice%nmax + 1
      llcheb = (2*this%control%lld) + 2
      lmax = lmax_basis
      nv = this%en%channels_ldos + 10

      call resolve_chebyshev_window(this, emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp

      wscale(:) = (this%en%ene(:) - b)/a

      ! Calculating the Jackson Kernel
      call jackson_kernel((this%control%lld), kernel)

      allocate(psiref(hblocksize, hblocksize, this%lattice%kk), left_vec(hblocksize, hblocksize, this%lattice%kk), mu_n_orb(2*(lmax + 1)**2, 2*(lmax + 1)**2, this%control%lld))
      allocate(mu_tmp(nb, nb, this%control%lld))
      allocate(w0(hblocksize, hblocksize, this%lattice%kk), w1(hblocksize, hblocksize, this%lattice%kk), right_vec(hblocksize, hblocksize, this%lattice%kk))
      allocate(w2(hblocksize, hblocksize, this%lattice%kk), v0(hblocksize, hblocksize, this%lattice%kk), v1(hblocksize, hblocksize, this%lattice%kk))
      allocate(v2(hblocksize, hblocksize, this%lattice%kk), left_vec1(hblocksize, hblocksize, this%lattice%kk), left_vec2(hblocksize, hblocksize, this%lattice%kk))

      mLz(:, :) = 0.0d0; mLy(:, :) = 0.0d0; mLx(:, :) = 0.0d0; L_op(:, :) = 0.0d0
      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
      mLx(:, :) = L_x(:, :)
      mLy(:, :) = L_y(:, :)
      mLz(:, :) = L_z(:, :)

      ! Transforming them into the spherical harmonics coordinates
      call hcpx(mLx, 'cart2sph')
      call hcpx(mLy, 'cart2sph')
      call hcpx(mLz, 'cart2sph')

      L_op(1:norb, 1:norb) = mLz(:, :)
      L_op(norb+1:nb, norb+1:nb) = mLz(:, :)

      ll = 1  
      do random = 1, this%lattice%kk
         call random_seed()

         ! Initializing wave functions
         v0(:, :, :) = (0.0d0, 0.0d0)
         v1(:, :, :) = (0.0d0, 0.0d0)
         v2(:, :, :) = (0.0d0, 0.0d0)
         psiref(:, :, :) = (0.0d0, 0.0d0)
         w0(:, :, :) = (0.0d0, 0.0d0)
         w1(:, :, :) = (0.0d0, 0.0d0)
         w2(:, :, :) = (0.0d0, 0.0d0)
         right_vec(:, :, :) = (0.0d0, 0.0d0)
         left_vec(:, :, :) = (0.0d0, 0.0d0)
         left_vec1(:, :, :) = (0.0d0, 0.0d0)
         left_vec2(:, :, :) = (0.0d0, 0.0d0)
         dum(:, :) = (0.0d0, 0.0d0)
         dum1(:, :) = (0.0d0, 0.0d0)
         dum2(:, :) = (0.0d0, 0.0d0)
         !mu_n_orb(:, :, :) = (0.0d0, 0.0d0)
         !L_op(:, :) = (0.0d0, 0.0d0)
        
         this%izero(:) = 1
         !do k = 1, this%lattice%kk
         !   call random_number(rng)
         !   do m = 1, nb
         !      psiref(m, m, k) =  (exp(2.0_rp * pi * i_unit * (rng))) !(2.0d0 * rng - 1.0d0)*sqrt(3.0d0) !(exp(2.0_rp * pi * i_unit * (rng)))
         !   end do
         !end do
         !! Normalize the full matrix 
         !psiref(:, :, :) = psiref(:, :, :) / sqrt(real(this%lattice%kk))
         write(*,*) random
         do m = 1, nb
            psiref(m, m, random) = (1.0d0, 0.0d0)
         end do
!         do k = 1, this%lattice%kk
!            nr = this%lattice%nn(k, 1)     ! Number of neighbors for this atom type
!            ntype = this%lattice%iz(k)
!            call zgemm('n', 'n', nb, nb, 18, cone, L_op(:, :), 18, psiref(:, :, k), 18, cone, left_vec(:, :, k), 18)
!            ! Loop over neighbors
!            do m = 2, nr   ! Start from 2 to exclude the onsite term
!               atom_neighbor = this%lattice%nn(k, m)  ! Neighbor atom number
!               if (atom_neighbor /= 0) then
!                  L_op(:, :) = (0.0d0, 0.0d0)
!                  call zgemm('n', 'n', nb, nb, 18, cone, this%hamiltonian%ee(:, :, m, ntype), 18, psiref(:, :, atom_neighbor), 18, cone, left_vec(:, :, k), 18)
!                  crossrij(:) = cross_product(this%lattice%cr(:, k), this%lattice%cr(:, k) - this%lattice%cr(:, atom_neighbor))
!                  L_op(:, :) = (i_unit * crossrij(3) * this%hamiltonian%ee(:, :, m, ntype))
!                  call zgemm('n', 'n', nb, nb, 18, cone, L_op(:, :), 18, psiref(:, :, atom_neighbor), 18, cone, left_vec(:, :, k), 18)
!                  left_vec(:, :, k) = left_vec(:, :, k) * i_unit * crossrij(3)
!               end if
!            end do
!         end do
         !left_vec(:, :, :) =  left_vec(:, :, :) - b*psiref(:, :, :)
         !left_vec(:, :, :) =  left_vec(:, :, :)/a

         ! X*|r>
         do k = 1, this%lattice%kk
            left_vec1(:, :, k) = this%lattice%cr(1, k) * psiref(:, :, k) * this%lattice%alat
         end do

         ! H*X|r>
         call this%ham_vec_matmul(left_vec1, w0, a, b)
         
         ! Y*H*X|r>
         left_vec1(:, :, :) = (0.0d0, 0.0d0)
         do k = 1, this%lattice%kk
            left_vec1(:, :, k) = this%lattice%cr(2, k) * w0(:, :, k) * this%lattice%alat
         end do


         w0(:, :, :) = (0.0d0, 0.0d0)
         ! Y*|r>
         do k = 1, this%lattice%kk
            left_vec2(:, :, k) = this%lattice%cr(2, k) * psiref(:, :, k) * this%lattice%alat
         end do
         
         ! H*Y|r>
         call this%ham_vec_matmul(left_vec2, w0, a, b)

         ! X*H*Y|r>
         left_vec2(:, :, :) = (0.0d0, 0.0d0)
         do k = 1, this%lattice%kk
            left_vec2(:, :, k) = this%lattice%cr(1, k) * w0(:, :, k) * this%lattice%alat
         end do

         left_vec(:, :, :) = i_unit * (left_vec1(:, :, :) - left_vec2(:, :, :))
         write(*,*) sum(left_vec), sum(psiref), sum(left_vec1), sum(left_vec2)

         ! GPU path: offload the T_n(H~) right recursion + left^H projection.
         ! left_vec (orbital current) is built on the host above, exactly as in
         ! the legacy path; the GPU only evaluates mu_n = sum_k left^H T_n|r>.
         if (gpu_plugin_ready(this, 'chebyshev_orbital_mod()', allow_hoh=.true.)) then
            if (random == 1) call gpu_plugin_upload_hamiltonian(this)
            call this%gpu_backend%orbital_moments(left_vec, psiref, this%control%lld, a, b, mu_tmp)
            mu_n_orb(:, :, :) = mu_n_orb(:, :, :) + mu_tmp(:, :, :)
            cycle
         end if

         ! Fast CPU orbital kernel (single precision, hoh-aware).
         if (this%hamiltonian%ccor_2c .and. this%hamiltonian%hoh .and. trim(this%control%cheb_backend) /= 'legacy') then
            call g_logger%warning('ccor_2c with hoh is not supported by fast orbital Chebyshev; falling back to legacy.', __FILE__, __LINE__)
         end if
         if (trim(this%control%cheb_backend) /= 'legacy' .and. &
             .not. (this%hamiltonian%ccor_2c .and. this%hamiltonian%hoh)) then
            if (this%hamiltonian%hoh) then
               call cheb_moments_orbital_fast(left_vec, psiref, this%hamiltonian%ee, &
                  this%hamiltonian%hall, this%hamiltonian%lsham, this%lattice%nn, &
                  this%lattice%iz, this%lattice%kk, nb, size(this%lattice%nn, 2), &
                  this%lattice%ntype, this%lattice%nmax, this%control%lld, a, b, mu_tmp, &
                  .true., this%hamiltonian%eeo, this%hamiltonian%hallo, this%hamiltonian%enim)
            else if (this%hamiltonian%ccor_2c) then
               call ensure_ccor_operator_blocks(this)
               call cheb_moments_orbital_fast(left_vec, psiref, this%ee_ccor_work, &
                  this%hall_ccor_work, this%hamiltonian%lsham, this%lattice%nn, &
                  this%lattice%iz, this%lattice%kk, nb, size(this%lattice%nn, 2), &
                  this%lattice%ntype, this%lattice%nmax, this%control%lld, a, b, mu_tmp, &
                  .false., this%ee_ccor_work, this%hall_ccor_work, this%hamiltonian%lsham)
            else
               call cheb_moments_orbital_fast(left_vec, psiref, this%hamiltonian%ee, &
                  this%hamiltonian%hall, this%hamiltonian%lsham, this%lattice%nn, &
                  this%lattice%iz, this%lattice%kk, nb, size(this%lattice%nn, 2), &
                  this%lattice%ntype, this%lattice%nmax, this%control%lld, a, b, mu_tmp, &
                  .false., this%hamiltonian%ee, this%hamiltonian%hall, this%hamiltonian%lsham)
            end if
            mu_n_orb(:, :, :) = mu_n_orb(:, :, :) + mu_tmp(:, :, :)
            cycle
         end if

         do n=1, this%control%lld
            call show_progress(n, this%control%lld)
            if (n == 1) then
               v1(:, :, :) = psiref(:, :, :)
            else if (n == 2) then
               v0(:, :, :) = v1(:, :, :)
               if (this%hamiltonian%hoh) then
                  call this%ham_hoh_vec_matmul(v0, v1, a, b)
               else
                  call this%ham_vec_matmul(v0, v1, a, b)
               end if
               this%izero(:) = this%idum(:)
            else if (n > 2) then
               if (this%hamiltonian%hoh) then
                  call this%ham_hoh_vec_matmul(v1, v2, a, b)
               else
                  call this%ham_vec_matmul(v1, v2, a, b)
               end if
               this%izero(:) = this%idum(:)
               v2(:, :, :) = 2 * v2(:, :, :) - v0(:, :, :)
               v0(:, :, :) = v1(:, :, :)
               v1(:, :, :) = v2(:, :, :)
               v2(:, :, :) = (0.0d0, 0.0d0)
            end if
            right_vec(:, :, :) = v1(:, :, :)
            dum(:, :) = (0.0d0, 0.0d0)
            !$omp parallel do default(shared) private(k) reduction(+:dum)  schedule(dynamic)
            do k=1, this%lattice%kk
               call zgemm('c', 'n', nb, nb, nb, cone, left_vec(:, :, k), nb, right_vec(:, :, k), nb, cone, dum(:, :), nb)
            end do
            !$omp end parallel do
            mu_n_orb(:, :, n) = mu_n_orb(:, :, n) + dum(:, :)
         end do
      end do

      mu_n_orb(:, :, :) = mu_n_orb(:, :, :) / real(this%lattice%kk)

      do l = 1, nb
         do m = 1, nb
            mu_n_orb(l, m, :) = mu_n_orb(l, m, :)*kernel(:)
         end do
      end do

      mu_n_orb(:, :, 2:size(kernel)) = mu_n_orb(:, :, 2:size(kernel))*2.0_rp

      g0(:, :, :) = (0.0d0, 0.0d0)
      !$omp parallel do default(shared) private(ie, i, exp_factor, l,m)
      do ie = 1, this%en%channels_ldos + 10
         do i = 1, size(kernel)
            exp_factor = -i_unit*exp(-i_unit*(i - 1)*acos(wscale(ie)))
            do l = 1, nb
               do m = 1, nb
                  g0(l, m, ie) = g0(l, m, ie) + mu_n_orb(l, m, i)*aimag(exp_factor)
               end do
            end do
         end do
         do l = 1, nb
            do m = 1, nb
               g0(l, m, ie) = g0(l, m, ie)/((sqrt((a**2) - ((this%en%ene(ie) - b)**2))))
            end do
         end do
         call g_kpm_profile%stop('T_cheb_moments')
      end do
      !$omp end parallel do
      do ie = 1, this%en%channels_ldos + 10
         lzi(ie) = rtrace(g0(:, :, ie)) 
      end do
      ! Log a quick diagnostic: maximum absolute value in g0
      maxg_tmp = 0.0_rp
      do ie_tmp = 1, this%en%channels_ldos + 10
         maxg_tmp = max(maxg_tmp, maxval(abs(g0(:, :, ie_tmp))))
      end do
      call g_logger%info('chebyshev_orbital_mod: max|g0|='//fmt('e12.5', maxg_tmp), __FILE__, __LINE__)
      do ie = 1, this%en%channels_ldos + 10
         call simpson_f(lz, this%en%ene, this%en%ene(ie), this%en%nv1, lzi, .true., .false., 0.0d0)
         write(50, '(3es16.6)') this%en%ene(ie) - this%en%fermi, -(lz/pi), -(1/pi)*lzi(ie)
      end do
   end subroutine chebyshev_orbital_mod

end submodule recursion_transport
