submodule(recursion_mod) recursion_transport

   use iso_fortran_env, only: int64
   use kpm_profile_mod, only: g_kpm_profile

contains

   !> @brief Apply the scaled Hamiltonian block matrix of one site type to a block state (on-site only).
   !> @param[in] matrix On-site blocks indexed by atom type.
   !> @param[in] psi Input block state.
   !> @param[out] out Output block state.
   module subroutine onsite_vec_matmul(this, matrix, psi, out)
      class(recursion), intent(in) :: this
      complex(rp), contiguous, intent(in) :: matrix(:, :, :), psi(:, :, :)
      complex(rp), contiguous, intent(out) :: out(:, :, :)
      integer :: i, hblocksize, nrhs
      hblocksize = size(psi, 1)
      nrhs = size(psi, 2)
      out(:, :, :) = (0.0_rp, 0.0_rp)
      !$omp parallel do default(shared) private(i)
      do i = 1, size(psi, 3)
         if (this%izero(i) == 0) cycle
         call zgemm('n', 'n', hblocksize, nrhs, hblocksize, cone, matrix(:, :, this%lattice%iz(i)), &
                    hblocksize, psi(:, :, i), hblocksize, czero, out(:, :, i), hblocksize)
      end do
      !$omp end parallel do
   end subroutine onsite_vec_matmul

   !> @brief Apply the unscaled Hamiltonian (with or without H_SOC) to a block state for the Kubo operators.
   !> @details Uses the HOH action when hamiltonian%hoh is set. work holds the HOH scratch planes 1:3.
   module subroutine kubo_ham_vec_matmul(this, psi, out, work, include_soc, a, b)
      class(recursion), intent(inout) :: this
      complex(rp), contiguous, intent(in) :: psi(:, :, :)
      complex(rp), contiguous, intent(out) :: out(:, :, :)
      complex(rp), contiguous, intent(inout) :: work(:, :, :, :)
      logical, intent(in), optional :: include_soc
      real(rp), intent(in), optional :: a, b
      real(rp) :: scale, shift
      logical :: with_soc

      scale = 1.0_rp
      shift = 0.0_rp
      with_soc = .true.
      if (present(a)) scale = a
      if (present(b)) shift = b
      if (present(include_soc)) with_soc = include_soc
      if (this%hamiltonian%hoh) then
         call this%ham_hoh_vec_matmul(psi, out, scale, shift, work, with_soc)
      else
         call this%ham_vec_matmul(psi, out, scale, shift, with_soc)
      end if
   end subroutine kubo_ham_vec_matmul

   !> @brief Apply one Kubo response operator (charge or current, local, torque) to a block state.
   !> @details kind selects the action; q is the local 18x18 operator (polarization), direction the velocity direction,
   !>          weighted applies the bond weights velocity_scale (slot b). work: planes 1:3 H/V scratch, 4:6 operator products.
   module subroutine kubo_operator_vec_matmul(this, kind, q, direction, weighted, psi, out, work)
      class(recursion), intent(inout) :: this
      character(len=*), intent(in) :: kind
      complex(rp), contiguous, intent(in) :: q(:, :), psi(:, :, :)
      real(rp), intent(in) :: direction(3)
      logical, intent(in) :: weighted
      complex(rp), contiguous, intent(out) :: out(:, :, :)
      ! Caller-owned scratch: H/V use planes 1:3, operator products use 4:6.
      ! psi, out, and scratch must not overlap.
      complex(rp), contiguous, intent(inout) :: work(:, :, :, :)
      integer :: i, hblocksize, nrhs
      hblocksize = size(psi, 1)
      nrhs = size(psi, 2)
      this%idum(:) = this%izero(:)
      if (kind == 'charge') then
         if (this%hamiltonian%hoh) then
            call this%velo_hoh_vec_matmul(psi, out, work(:, :, :, 1:3), direction, weighted)
         else
            call this%velo_vec_matmul(psi, out, direction, weighted)
         end if
         return
      end if

      work(:, :, :, 4) = (0.0_rp, 0.0_rp)
      !$omp parallel do default(shared) private(i)
      do i = 1, size(psi, 3)
         if (this%izero(i) == 0) cycle
         call zgemm('n', 'n', hblocksize, nrhs, hblocksize, cone, q, hblocksize, &
                    psi(:, :, i), hblocksize, czero, work(:, :, i, 4), hblocksize)
      end do
      !$omp end parallel do
      if (kind == 'local') then
         out = work(:, :, :, 4)
         return
      end if
      select case(kind)
      case('current')
         ! Keep V*psi and V*(q*psi) separate until the final anticommutator.
         if (this%hamiltonian%hoh) then
            call this%velo_hoh_vec_matmul(psi, work(:, :, :, 5), work(:, :, :, 1:3), direction, weighted)
         else
            call this%velo_vec_matmul(psi, work(:, :, :, 5), direction, weighted)
         end if
         if (this%hamiltonian%hoh) then
            call this%velo_hoh_vec_matmul(work(:, :, :, 4), work(:, :, :, 6), work(:, :, :, 1:3), direction, weighted)
         else
            call this%velo_vec_matmul(work(:, :, :, 4), work(:, :, :, 6), direction, weighted)
         end if
      case('spin_torque')
         call this%kubo_ham_vec_matmul(psi, work(:, :, :, 5), work(:, :, :, 1:3), .false.)
         call this%kubo_ham_vec_matmul(work(:, :, :, 4), work(:, :, :, 6), work(:, :, :, 1:3), .false.)
      case('orbital_torque')
         call this%kubo_ham_vec_matmul(psi, work(:, :, :, 5), work(:, :, :, 1:3))
         call this%kubo_ham_vec_matmul(work(:, :, :, 4), work(:, :, :, 6), work(:, :, :, 1:3))
      case('soc_torque')
         call this%onsite_vec_matmul(this%hamiltonian%lsham, psi, work(:, :, :, 5))
         call this%onsite_vec_matmul(this%hamiltonian%lsham, work(:, :, :, 4), work(:, :, :, 6))
      case default
         call g_logger%fatal('Invalid Kubo response kind: '//trim(kind), __FILE__, __LINE__)
      end select
      ! Left product q*(V*psi) or q*(Hpart*psi), then add/subtract right product.
      out(:, :, :) = (0.0_rp, 0.0_rp)
      !$omp parallel do default(shared) private(i)
      do i = 1, size(psi, 3)
         if (this%idum(i) == 0) cycle
         call zgemm('n', 'n', hblocksize, nrhs, hblocksize, cone, q, hblocksize, &
                    work(:, :, i, 5), hblocksize, czero, out(:, :, i), hblocksize)
         if (kind == 'current') then
            out(:, :, i) = 0.5_rp*(out(:, :, i) + work(:, :, i, 6))
         else
            out(:, :, i) = cmplx(0.0_rp, -1.0_rp, rp)*(out(:, :, i) - work(:, :, i, 6))
         end if
      end do
      !$omp end parallel do
   end subroutine kubo_operator_vec_matmul

   !> @brief Apply the matrix-free velocity -i direction.d_ij H_ij to a block state.
   !> @details The physical bond vectors d_ij are cached in transport_bond (periodic images included). weighted applies velocity_scale.
   !> @param[inout] this Recursion object; updates idum to the propagated support.
   !> @param[in] psi_in Input block state, site-major as (nb,nrhs,nsite).
   !> @param[out] psi_out Output block state in the same layout.
   !> @param[in] direction Cartesian velocity direction.
   !> @param[in] weighted Apply the bond weights.
   module subroutine velo_vec_matmul(this, psi_in, psi_out, direction, weighted)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      real(rp), intent(in) :: direction(3)
      logical, intent(in) :: weighted
      ! Local variables
      integer :: ineigh, ih, k, nr, nrhs, nnmap, nlimplus1
      complex(rp) :: factor
      real(rp) :: distance

      if (.not. allocated(this%transport_bond)) call this%setup_kubo_bonds()
      nrhs = size(psi_in, 2)
      nlimplus1 = this%lattice%nmax + 1
      psi_out(:, :, :) = (0.0_rp, 0.0_rp)

      ! Local Hamiltonian: site-specific impurity-region blocks.
      if (this%lattice%nmax /= 0) then
         !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap, distance, factor)
         do k = 1, this%lattice%nmax
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1)
            ! The on-site velocity vanishes because r_k-r_k = 0.
            if (nr >= 2) then
               do ineigh = 2, nr
                  nnmap = this%lattice%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        distance = dot_product(direction, this%transport_bond(:, ineigh, k))
                        factor = cmplx(0.0_rp, -distance, rp)
                        if (weighted) then
                           factor = factor*max(this%hamiltonian%velocity_scale(ih), &
                                               this%hamiltonian%velocity_scale(this%lattice%iz(nnmap)))
                        end if
                        call zgemm('n', 'n', nb, nrhs, nb, factor, this%hamiltonian%hall(:, :, ineigh, k), nb, &
                                   psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
                        this%idum(k) = 1
                     end if
                  end if
               end do
            end if
         end do
         !$omp end parallel do
      end if

      ! Bulk Hamiltonian: blocks indexed by atom type.
      !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap, distance, factor)
      do k = nlimplus1, this%lattice%kk
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         ! The on-site velocity vanishes because r_k-r_k = 0.
         if (nr >= 2) then
            do ineigh = 2, nr
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     distance = dot_product(direction, this%transport_bond(:, ineigh, k))
                     factor = cmplx(0.0_rp, -distance, rp)
                     if (weighted) then
                        factor = factor*max(this%hamiltonian%velocity_scale(ih), &
                                            this%hamiltonian%velocity_scale(this%lattice%iz(nnmap)))
                     end if
                     call zgemm('n', 'n', nb, nrhs, nb, factor, this%hamiltonian%ee(:, :, ineigh, ih), nb, &
                                psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
                     this%idum(k) = 1
                  end if
               end if
            end do
         end if
      end do
      !$omp end parallel do
   end subroutine velo_vec_matmul

   !> @brief Apply the matrix-free velocity with the HOH corrections v - v O h - h O v to a block state.
   !> @param[inout] this Recursion object; mutates izero/idum while propagating the active light cone.
   !> @param[in] psi_in Input block state.
   !> @param[out] psi_out Output block state.
   !> @param[inout] work Scratch planes 1:3 (psi1, psi2, hohpsi).
   !> @param[in] direction Cartesian velocity direction.
   !> @param[in] weighted Apply the bond weights.
   module subroutine velo_hoh_vec_matmul(this, psi_in, psi_out, work, direction, weighted)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      complex(rp), target, contiguous, intent(inout) :: work(:, :, :, :)
      real(rp), intent(in) :: direction(3)
      logical, intent(in) :: weighted
      ! Local variables
      integer :: ineigh, ih, k, nr, nrhs, nnmap, nlimplus1
      integer :: input_map(size(psi_in, 3))
      complex(rp), dimension(nb, nb) :: locham
      complex(rp), pointer :: psi1(:, :, :), psi2(:, :, :), hohpsi(:, :, :)
      complex(rp) :: factor
      real(rp) :: distance

      ! Name existing storage, matching the original intermediate wavefunctions.
      psi1 => work(:, :, :, 1)
      psi2 => work(:, :, :, 2)
      hohpsi => work(:, :, :, 3)
      if (.not. allocated(this%transport_bond)) call this%setup_kubo_bonds()
      nrhs = size(psi_in, 2)
      nlimplus1 = this%lattice%nmax + 1
      input_map(:) = this%izero(:)
      psi1(:, :, :) = (0.0_rp, 0.0_rp)
      psi2(:, :, :) = (0.0_rp, 0.0_rp)
      hohpsi(:, :, :) = (0.0_rp, 0.0_rp)

      ! First pass: v*psi and h*psi, as in the original velocity routine.

      ! Local Hamiltonian: site-specific impurity-region blocks.
      if (this%lattice%nmax /= 0) then
         !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap, locham, distance, factor)
         do k = 1, this%lattice%nmax
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1)
            if (this%izero(k) /= 0) then
               call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%hall(:, :, 1, k), nb, &
                          psi_in(:, :, k), nb, cone, psi1(:, :, k), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr
                  nnmap = this%lattice%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%hall(:, :, ineigh, k), nb, &
                                   psi_in(:, :, nnmap), nb, cone, psi1(:, :, k), nb)
                        distance = dot_product(direction, this%transport_bond(:, ineigh, k))
                        factor = cmplx(0.0_rp, -distance, rp)
                        if (weighted) then
                           factor = factor*max(this%hamiltonian%velocity_scale(ih), &
                                               this%hamiltonian%velocity_scale(this%lattice%iz(nnmap)))
                        end if
                        call zgemm('n', 'n', nb, nrhs, nb, factor, this%hamiltonian%hall(:, :, ineigh, k), nb, &
                                   psi_in(:, :, nnmap), nb, cone, psi2(:, :, k), nb)
                        this%idum(k) = 1
                     end if
                  end if
               end do
            end if
         end do
         !$omp end parallel do
      end if

      ! Bulk Hamiltonian: blocks indexed by atom type.
      !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap, locham, distance, factor)
      do k = nlimplus1, this%lattice%kk
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%ee(:, :, 1, ih), nb, &
                       psi_in(:, :, k), nb, cone, psi1(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih), nb, &
                                psi_in(:, :, nnmap), nb, cone, psi1(:, :, k), nb)
                     distance = dot_product(direction, this%transport_bond(:, ineigh, k))
                     factor = cmplx(0.0_rp, -distance, rp)
                     if (weighted) then
                        factor = factor*max(this%hamiltonian%velocity_scale(ih), &
                                            this%hamiltonian%velocity_scale(this%lattice%iz(nnmap)))
                     end if
                     call zgemm('n', 'n', nb, nrhs, nb, factor, this%hamiltonian%ee(:, :, ineigh, ih), nb, &
                                psi_in(:, :, nnmap), nb, cone, psi2(:, :, k), nb)
                     this%idum(k) = 1
                  end if
               end if
            end do
         end if
      end do
      !$omp end parallel do

      ! Mapping update for both second-hop corrections.
      this%izero(:) = this%idum(:)

      ! Local Hamiltonian: site-specific impurity-region blocks.
      if (this%lattice%nmax /= 0) then
         !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap, locham, distance, factor)
         do k = 1, this%lattice%nmax
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1)
            if (this%izero(k) /= 0) then
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, 1, k), nb, &
                          this%hamiltonian%obarm(:, :, this%lattice%iz(k)), nb, czero, locham, nb)
               ! Missing old term: (h*o)*(v*psi), including the on-site h*o block.
               call zgemm('n', 'n', nb, nrhs, nb, cone, locham, nb, &
                          psi2(:, :, k), nb, cone, hohpsi(:, :, k), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr
                  nnmap = this%lattice%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, ineigh, k), nb, &
                                   this%hamiltonian%obarm(:, :, this%lattice%iz(nnmap)), nb, czero, locham, nb)
                        ! Missing old term: (h*o)*(v*psi), including the on-site h*o block.
                        call zgemm('n', 'n', nb, nrhs, nb, cone, locham, nb, &
                                   psi2(:, :, nnmap), nb, cone, hohpsi(:, :, k), nb)
                        ! Original term: (v*o)*(h*psi), with v_ij = factor*h_ij.
                        distance = dot_product(direction, this%transport_bond(:, ineigh, k))
                        factor = cmplx(0.0_rp, -distance, rp)
                        if (weighted) then
                           factor = factor*max(this%hamiltonian%velocity_scale(ih), &
                                               this%hamiltonian%velocity_scale(this%lattice%iz(nnmap)))
                        end if
                        call zgemm('n', 'n', nb, nrhs, nb, factor, locham, nb, &
                                   psi1(:, :, nnmap), nb, cone, hohpsi(:, :, k), nb)
                        this%idum(k) = 1
                     end if
                  end if
               end do
            end if
         end do
         !$omp end parallel do
      end if

      ! Bulk Hamiltonian: blocks indexed by atom type.
      !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap, locham, distance, factor)
      do k = nlimplus1, this%lattice%kk
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, 1, ih), nb, &
                       this%hamiltonian%obarm(:, :, this%lattice%iz(k)), nb, czero, locham, nb)
            ! Missing old term: (h*o)*(v*psi), including the on-site h*o block.
            call zgemm('n', 'n', nb, nrhs, nb, cone, locham, nb, &
                       psi2(:, :, k), nb, cone, hohpsi(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih), nb, &
                                this%hamiltonian%obarm(:, :, this%lattice%iz(nnmap)), nb, czero, locham, nb)
                     ! Missing old term: (h*o)*(v*psi), including the on-site h*o block.
                     call zgemm('n', 'n', nb, nrhs, nb, cone, locham, nb, &
                                psi2(:, :, nnmap), nb, cone, hohpsi(:, :, k), nb)
                     ! Original term: (v*o)*(h*psi), with v_ij = factor*h_ij.
                     distance = dot_product(direction, this%transport_bond(:, ineigh, k))
                     factor = cmplx(0.0_rp, -distance, rp)
                     if (weighted) then
                        factor = factor*max(this%hamiltonian%velocity_scale(ih), &
                                            this%hamiltonian%velocity_scale(this%lattice%iz(nnmap)))
                     end if
                     call zgemm('n', 'n', nb, nrhs, nb, factor, locham, nb, &
                                psi1(:, :, nnmap), nb, cone, hohpsi(:, :, k), nb)
                     this%idum(k) = 1
                  end if
               end if
            end do
         end if
      end do
      !$omp end parallel do

      ! Complete velocity: v*psi - v*o*h*psi - h*o*v*psi.
      psi_out(:, :, :) = psi2(:, :, :) - hohpsi(:, :, :)
      this%izero(:) = input_map(:)
      psi2(:, :, :) = (0.0_rp, 0.0_rp)
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
   !> @param[inout] work Optional scratch planes 1:2 (psi2, hohpsi); defaults to the recursion members.
   !> @param[in] include_soc Optional; false omits H_SOC (default true).
   !> @note Mutates izero/idum to the support reached by the Hamiltonian action.
   module subroutine ham_hoh_vec_matmul(this, psi_in, psi_out, a, b, work, include_soc)
      class(recursion), target, intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      real(rp), intent(in) :: a, b
      complex(rp), target, contiguous, optional, intent(inout) :: work(:, :, :, :)
      logical, optional, intent(in) :: include_soc
      ! Local variables
      integer :: ineigh, ih, k, nr, nrhs, nnmap, nlimplus1
      integer :: input_map(size(psi_in, 3))
      complex(rp), dimension(nb, nb) :: locham
      complex(rp), pointer :: psi2(:, :, :), hohpsi(:, :, :)
      logical :: with_soc

      ! Name existing storage; no allocation or copy is performed here.
      if (present(work)) then
         psi2 => work(:, :, :, 1)
         hohpsi => work(:, :, :, 2)
      else
         psi2 => this%psi2
         hohpsi => this%hohpsi
      end if
      nrhs = size(psi_in, 2)
      nlimplus1 = this%lattice%nmax + 1
      with_soc = .true.
      if (present(include_soc)) with_soc = include_soc
      input_map(:) = this%izero(:)
      psi2(:, :, :) = (0.0_rp, 0.0_rp)
      hohpsi(:, :, :) = (0.0_rp, 0.0_rp)
      psi_out(:, :, :) = (0.0_rp, 0.0_rp)

      ! First pass: h*psi and the on-site E_nu/SOC terms.

      ! Local Hamiltonian: site-specific impurity-region blocks.
      if (this%lattice%nmax /= 0) then
         !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap, locham)
         do k = 1, this%lattice%nmax
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1)
            if (this%izero(k) /= 0) then
               call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%hall(:, :, 1, k), nb, &
                          psi_in(:, :, k), nb, cone, psi2(:, :, k), nb)
               locham = this%hamiltonian%enim(:, :, ih)
               if (with_soc) locham = locham + this%hamiltonian%lsham(:, :, ih)
               call zgemm('n', 'n', nb, nrhs, nb, cone, locham, nb, &
                          psi_in(:, :, k), nb, cone, psi_out(:, :, k), nb)
               if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%hallcc(:, :, 1, k), nb, &
                                                         psi_in(:, :, k), nb, cone, psi_out(:, :, k), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr
                  nnmap = this%lattice%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%hall(:, :, ineigh, k), nb, &
                                   psi_in(:, :, nnmap), nb, cone, psi2(:, :, k), nb)
                        if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%hallcc(:, :, ineigh, k), nb, &
                                   psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
                        this%idum(k) = 1
                     end if
                  end if
               end do
            end if
         end do
         !$omp end parallel do
      end if

      ! Bulk Hamiltonian: blocks indexed by atom type.
      !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap, locham)
      do k = nlimplus1, this%lattice%kk
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%ee(:, :, 1, ih), nb, &
                       psi_in(:, :, k), nb, cone, psi2(:, :, k), nb)
            locham = this%hamiltonian%enim(:, :, ih)
            if (with_soc) locham = locham + this%hamiltonian%lsham(:, :, ih)
            call zgemm('n', 'n', nb, nrhs, nb, cone, locham, nb, &
                       psi_in(:, :, k), nb, cone, psi_out(:, :, k), nb)
            if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%eecc(:, :, 1, ih), nb, &
                                                      psi_in(:, :, k), nb, cone, psi_out(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih), nb, &
                                psi_in(:, :, nnmap), nb, cone, psi2(:, :, k), nb)
                     if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%eecc(:, :, ineigh, ih), nb, &
                                psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
                     this%idum(k) = 1
                  end if
               end if
            end do
         end if
      end do
      !$omp end parallel do

      ! Mapping update for the second hop in h*o*h.
      this%izero(:) = this%idum(:)

      ! Local Hamiltonian: site-specific impurity-region blocks.
      if (this%lattice%nmax /= 0) then
         !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap, locham)
         do k = 1, this%lattice%nmax
            this%idum(k) = this%izero(k)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k, 1)
            if (this%izero(k) /= 0) then
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, 1, k), nb, &
                          this%hamiltonian%obarm(:, :, this%lattice%iz(k)), nb, czero, locham, nb)
               call zgemm('n', 'n', nb, nrhs, nb, cone, locham, nb, &
                          psi2(:, :, k), nb, cone, hohpsi(:, :, k), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr
                  nnmap = this%lattice%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, ineigh, k), nb, &
                                   this%hamiltonian%obarm(:, :, this%lattice%iz(nnmap)), nb, czero, locham, nb)
                        call zgemm('n', 'n', nb, nrhs, nb, cone, locham, nb, &
                                   psi2(:, :, nnmap), nb, cone, hohpsi(:, :, k), nb)
                        this%idum(k) = 1
                     end if
                  end if
               end do
            end if
         end do
         !$omp end parallel do
      end if

      ! Bulk Hamiltonian: blocks indexed by atom type.
      !$omp parallel do default(shared) private(k, ih, nr, ineigh, nnmap, locham)
      do k = nlimplus1, this%lattice%kk
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, 1, ih), nb, &
                       this%hamiltonian%obarm(:, :, this%lattice%iz(k)), nb, czero, locham, nb)
            call zgemm('n', 'n', nb, nrhs, nb, cone, locham, nb, &
                       psi2(:, :, k), nb, cone, hohpsi(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih), nb, &
                                this%hamiltonian%obarm(:, :, this%lattice%iz(nnmap)), nb, czero, locham, nb)
                     call zgemm('n', 'n', nb, nrhs, nb, cone, locham, nb, &
                                psi2(:, :, nnmap), nb, cone, hohpsi(:, :, k), nb)
                     this%idum(k) = 1
                  end if
               end if
            end do
         end if
      end do
      !$omp end parallel do

      ! H = h - h*o*h + E_nu + H_SOC.
      psi_out(:, :, :) = psi2(:, :, :) - hohpsi(:, :, :) + psi_out(:, :, :)

      ! Do the scaling and shifting only after assembling the complete H.
      psi_out(:, :, :) = psi_out(:, :, :) - b*psi_in(:, :, :)
      psi_out(:, :, :) = psi_out(:, :, :)/a
      this%izero(:) = input_map(:)
      psi2(:, :, :) = (0.0_rp, 0.0_rp)
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
   !> @param[in] include_soc Optional; false omits H_SOC (default true).
   !> @note Includes ccor_2c blocks when hamiltonian%ccor_2c is set.
   module subroutine ham_vec_matmul(this, psi_in, psi_out, a, b, include_soc)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(:, :, :), intent(in) :: psi_in
      complex(rp), dimension(:, :, :), intent(out) :: psi_out
      real(rp), intent(in) :: a, b
      logical, intent(in), optional :: include_soc
      logical :: with_soc
      integer :: nrhs
      ! Local variables
      integer :: ineigh, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(nb, nb) :: dum1, dum2, locham
      hblocksize = nb
      nat = this%lattice%kk
      with_soc = .true.
      if (present(include_soc)) with_soc = include_soc
      nrhs = size(psi_in, 2)

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
               locham = this%hamiltonian%hall(1:nb, 1:nb, 1, k)
               if (with_soc) locham = locham + this%hamiltonian%lsham(1:nb, 1:nb, ih)
               if (this%hamiltonian%ccor_2c) locham = locham + this%hamiltonian%hallcc(1:nb, 1:nb, 1, k)
               call zgemm('n', 'n', nb, nrhs, nb, cone, locham, nb, psi_in(:, :, k), nb, cone, psi_out(:, :, k), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        if (this%hamiltonian%ccor_2c) then
                           call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%hall(:, :, ineigh, k) + this%hamiltonian%hallcc(:, :, ineigh, k), nb, psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
                        else
                           call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%hall(:, :, ineigh, k), nb, psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
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
            locham = this%hamiltonian%ee(1:nb, 1:nb, 1, ih)
            if (with_soc) locham = locham + this%hamiltonian%lsham(1:nb, 1:nb, ih)
            if (this%hamiltonian%ccor_2c) locham = locham + this%hamiltonian%eecc(1:nb, 1:nb, 1, ih)
            call zgemm('n', 'n', nb, nrhs, nb, cone, locham, nb, psi_in(:, :, k), nb, cone, psi_out(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap == 0) cycle
               if (this%izero(nnmap) /= 0) then
                  if (this%hamiltonian%ccor_2c) then
                     call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih) + this%hamiltonian%eecc(:, :, ineigh, ih), nb, psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
                  else
                     call zgemm('n', 'n', nb, nrhs, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih), nb, psi_in(:, :, nnmap), nb, cone, psi_out(:, :, k), nb)
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

   !> @brief Cache the physical bond vectors d_ij of every site and neighbor slot in transport_bond.
   !> @details Periodic systems use the wrapped minimum image (f_wrap_coord_diff); ambiguous images are fatal.
   module subroutine setup_kubo_bonds(this)
      class(recursion), intent(inout) :: this
      real(rp), allocatable :: coord(:, :)
      real(rp) :: bond(3)
      integer :: k, ineigh, j
      logical :: ambiguous

      if (allocated(this%transport_bond)) deallocate(this%transport_bond)
      allocate(this%transport_bond(3, size(this%lattice%nn, 2), this%lattice%kk))
      this%transport_bond = 0.0_rp
      if (this%lattice%pbc) then
         allocate(coord(3, this%lattice%kk))
         coord = this%lattice%cr(:, 1:this%lattice%kk)*this%lattice%alat
      end if
      do k = 1, this%lattice%kk
         do ineigh = 2, this%lattice%nn(k, 1)
            j = this%lattice%nn(k, ineigh)
            if (j == 0) cycle
            if (this%lattice%pbc) then
               call this%lattice%f_wrap_coord_diff(this%lattice%kk, coord, k, j, bond, ambiguous)
               if (ambiguous) then
                  call g_logger%fatal('Ambiguous periodic transport bond: enlarge the supercell', __FILE__, __LINE__)
               end if
               this%transport_bond(:, ineigh, k) = -bond
            else
               this%transport_bond(:, ineigh, k) = &
                  (this%lattice%cr(:, k)-this%lattice%cr(:, j))*this%lattice%alat
            end if
         end do
      end do
   end subroutine setup_kubo_bonds

   !> @brief Select the operators of the output (slot a) and input (slot b) of the linear response and their directions.
   module subroutine setup_kubo_operators(this, linear_out, linear_in)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      class(recursion), intent(inout) :: this
      character(len=*), intent(in) :: linear_out, linear_in
      real(rp) :: norm_a, norm_b
      call this%setup_kubo_bonds()
      norm_a = norm2(this%hamiltonian%v_alpha)
      norm_b = norm2(this%hamiltonian%v_beta)
      if (.not. ieee_is_finite(norm_a) .or. .not. ieee_is_finite(norm_b)) then
         call g_logger%fatal('Transport directions must be finite', __FILE__, __LINE__)
      end if
      if (norm_a <= tiny(1.0_rp) .or. norm_b <= tiny(1.0_rp)) then
         call g_logger%fatal('Transport directions must be nonzero', __FILE__, __LINE__)
      end if
      this%transport_dir_a = this%hamiltonian%v_alpha/norm_a
      this%transport_dir_b = this%hamiltonian%v_beta/norm_b
      call this%set_kubo_operator_slot(linear_in, 'b')
      call this%set_kubo_operator_slot(linear_out, 'a')
      ! All transport intermediates traverse the complete cluster. The legacy
      ! recursion mask is not used by this path and cannot truncate two-hop terms.
      this%idum = 1
   end subroutine setup_kubo_operators

   !> @brief Configure transport_kind_x and transport_q_x of one slot from the operator type and the polarization.
   module subroutine set_kubo_operator_slot(this, op_type, slot)
      class(recursion), intent(inout) :: this
      character(len=*), intent(in) :: op_type, slot
      character(len=10) :: pol
      character(len=16) :: kind
      complex(rp) :: q(nb,nb)
      select case(slot)
      case('a')
         pol = this%hamiltonian%pol_alpha
      case('b')
         pol = this%hamiltonian%pol_beta
      case default
         call g_logger%fatal('Invalid Kubo operator slot: '//trim(slot), __FILE__, __LINE__)
      end select
      q = (0.0_rp, 0.0_rp)
      kind = 'current'
      select case(trim(op_type))
      case('charge')
         kind = 'charge'
      case('spin', 'spin_accumulation', 'spin_torque', 'spin_soc_torque', 'soc_spin_torque')
         call select_spin_operator(pol, q)
      case('orbital', 'orbital_accumulation', 'orbital_torque')
         call select_orbital_operator(pol, q)
      case('quadrupole', 'oap', 'quadrupole_accumulation', 'oap_accumulation')
         call select_quadrupole_operator(pol, q)
      case('octupole', 'octupole_current', 'octupole_accumulation')
         call select_octupole_operator(pol, q)
      case default
         call g_logger%fatal('Unknown linear-response operator type: '//trim(op_type), __FILE__, __LINE__)
      end select
      select case(trim(op_type))
      case('spin_accumulation', 'orbital_accumulation', 'quadrupole_accumulation', &
           'oap_accumulation', 'octupole_accumulation')
         kind = 'local'
      case('spin_torque', 'orbital_torque')
         kind = trim(op_type)
      case('spin_soc_torque', 'soc_spin_torque')
         kind = 'soc_torque'
      end select
      select case(slot)
      case('a')
         this%transport_kind_a = kind
         this%transport_q_a = q
      case('b')
         this%transport_kind_b = kind
         this%transport_q_b = q
      end select
   end subroutine set_kubo_operator_slot

   !> @brief Apply the operator of slot a (linear_out) or b (linear_in) to a block state.
   module subroutine apply_kubo_operator(this, slot, psi, out, work)
      class(recursion), intent(inout) :: this
      character(len=*), intent(in) :: slot
      complex(rp), contiguous, intent(in) :: psi(:,:,:)
      complex(rp), contiguous, intent(out) :: out(:,:,:)
      complex(rp), contiguous, intent(inout) :: work(:, :, :, :)
      select case(slot)
      case('a')
         call this%kubo_operator_vec_matmul(this%transport_kind_a, this%transport_q_a, &
              this%transport_dir_a, .false., psi, out, work)
      case('b')
         call this%kubo_operator_vec_matmul(this%transport_kind_b, this%transport_q_b, &
              this%transport_dir_b, .true., psi, out, work)
      case default
         call g_logger%fatal('Invalid Kubo slot', __FILE__, __LINE__)
      end select
   end subroutine apply_kubo_operator

   !> @brief Compute two-index stochastic Chebyshev moments for transport.
   !> @details Builds mu_nm_stochastic = <r|T_m(H) v_a T_n(H) v_b|r> for charge,
   !>          spin, orbital, accumulation, and torque conductivity routes. The
   !>          starting states are either representative atom types or random
   !>          vectors, selected by control%cond_calctype.
   !> @param[inout] this Recursion object; allocates and fills mu_nm_stochastic.
   !> @note The operators are applied matrix-free (setup_kubo_operators); GPU requests are fatal until the kernels return.
   module subroutine compute_moments_stochastic(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, m, n, hblocksize, nat, loop_over, lmax
      integer :: seed_size, seed_i, seed_status, seed_ios, random_seed_value, seed_length
      integer, allocatable :: seed_values(:)
      integer(int64) :: matrix_dimension, nnz, host_moment_bytes
      complex(rp), dimension(nb, nb) :: dum
      complex(rp), dimension(:, :, :), allocatable :: psiref, w0, w1, w2, right_vec, v0, v1, v2
      complex(rp), dimension(:, :, :, :), allocatable :: left_vec, work
      real(rp) :: a, b, rng, emin_win, emax_win
      logical :: deterministic_random_seed
      character(len=64) :: random_seed_text

      lmax = lmax_basis
      hblocksize = nb
      nat = this%lattice%kk

      ! General procedures
      call resolve_chebyshev_window(this, emin_win, emax_win)
      ! Same scaling constant as calculate_gamma_nm, so that moments and reconstruction agree.
      a = (emax_win - emin_win)/(2 - 0.3)
      b = (emax_win + emin_win)/2.0_rp

      ! Transport applies the operators matrix-free (setup_kubo_operators), whatever
      ! cheb_backend says. The fast CPU and CUDA kernels return with the backend switch.
      if (gpu_plugin_ready(this, 'compute_moments_stochastic()', allow_hoh=.true.)) &
         call g_logger%fatal('GPU transport moments are not available with the Kubo operators yet', __FILE__, __LINE__)
      if (trim(this%control%cheb_backend) /= 'legacy') &
         call g_logger%info('Transport moments use the matrix-free Kubo operators for every cheb_backend', __FILE__, __LINE__)

      select case(this%control%cond_calctype)
      case('per_type')
         loop_over = this%lattice%ntype
      case('random_vec')
         loop_over = this%control%random_vec_num
      end select

      ! For paired CPU evidence, an optional benchmark environment seed
      ! initializes the intrinsic Fortran generator once per transport request.
      ! Without it, retain the historical random_seed() behavior.
      deterministic_random_seed = .false.
      random_seed_text = ''
      if (trim(this%control%cond_calctype) == 'random_vec') then
         call get_environment_variable('RSLMTO_KPM_RANDOM_SEED', random_seed_text, length=seed_length, status=seed_status)
         seed_length = min(seed_length, len(random_seed_text))
         if (seed_status == 0 .and. seed_length > 0) then
            read(random_seed_text(1:seed_length), *, iostat=seed_ios) random_seed_value
            if (seed_ios == 0) then
               call random_seed(size=seed_size)
               allocate(seed_values(seed_size))
               do seed_i = 1, seed_size
                  seed_values(seed_i) = random_seed_value + 104729*(seed_i - 1)
               end do
               call random_seed(put=seed_values)
               deallocate(seed_values)
               deterministic_random_seed = .true.
               call g_logger%info('KPM random-vector seed='//trim(random_seed_text(1:seed_length)), __FILE__, __LINE__)
            end if
         end if
      end if

#ifdef USE_SAFE_ALLOC
      select case(this%control%cond_calctype)
      case('per_type')
         call g_safe_alloc%allocate('recursion.mu_nm_stochastic', this%mu_nm_stochastic, &
            (/2*(lmax + 1)**2, 2*(lmax + 1)**2, this%control%cond_ll, this%control%cond_ll, this%lattice%ntype/))
      case('random_vec')
         call g_safe_alloc%allocate('recursion.mu_nm_stochastic', this%mu_nm_stochastic, &
            (/2*(lmax + 1)**2, 2*(lmax + 1)**2, this%control%cond_ll, this%control%cond_ll, this%control%random_vec_num/))
      end select
#else
      select case(this%control%cond_calctype)
      case('per_type')
         allocate(this%mu_nm_stochastic(2*(lmax + 1)**2, 2*(lmax + 1)**2, this%control%cond_ll, &
            this%control%cond_ll, this%lattice%ntype))
      case('random_vec')
         allocate(this%mu_nm_stochastic(2*(lmax + 1)**2, 2*(lmax + 1)**2, this%control%cond_ll, &
            this%control%cond_ll, this%control%random_vec_num))
      end select
#endif
      host_moment_bytes = int(size(this%mu_nm_stochastic), int64)*int(storage_size(this%mu_nm_stochastic)/8, int64)
      call g_kpm_profile%set_moment_memory(0_int64, host_moment_bytes)

      allocate(psiref(hblocksize, hblocksize, this%lattice%kk), &
         left_vec(hblocksize, hblocksize, this%lattice%kk, this%control%cond_ll))
      allocate(w0(hblocksize, hblocksize, this%lattice%kk), &
         w1(hblocksize, hblocksize, this%lattice%kk), right_vec(hblocksize, hblocksize, this%lattice%kk))
      allocate(w2(hblocksize, hblocksize, this%lattice%kk), &
         v0(hblocksize, hblocksize, this%lattice%kk), v1(hblocksize, hblocksize, this%lattice%kk))
      allocate(v2(hblocksize, hblocksize, this%lattice%kk), work(hblocksize, hblocksize, this%lattice%kk, 6))

      matrix_dimension = int(nb, int64)*int(nat, int64)
      if (size(this%lattice%nn, 2) > 1) then
         nnz = int(count(this%lattice%nn(:, 2:) > 0), int64)*int(nb, int64)*int(nb, int64)
      else
         nnz = 0_int64
      end if
      call g_kpm_profile%configure('cpu_legacy', 'fp64', 'cpu_blas', 'fp64', &
         trim(this%control%cond_calctype), matrix_dimension, nnz, this%control%cond_ll, &
         this%control%lld, loop_over)

      call g_kpm_profile%start('P_operator_setup')
      call this%setup_kubo_operators(this%control%linear_out, this%control%linear_in)
      call g_kpm_profile%stop('P_operator_setup')

      do i = 1, loop_over
         call g_kpm_profile%start('P_trace_setup')
         if (.not. deterministic_random_seed) call random_seed()

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
                  psiref(m, m, k) =  (exp(2.0_rp * pi * i_unit * (rng)))
               end do
            end do
            ! Normalize the full matrix
            psiref(:, :, :) = psiref(:, :, :) / sqrt(real(this%lattice%kk))
         end select
         call g_kpm_profile%stop('P_trace_setup')

         call g_kpm_profile%start('P_moments_total')

         ! Computing the left vector <r|Tm(H)
         do m=1, this%control%cond_ll
            if (m == 1) then
               w1(:, :, :) = psiref(:, :, :)
            else if (m == 2) then
               w0(:, :, :) = w1(:, :, :)
               call this%kubo_ham_vec_matmul(w0, w1, work, a=a, b=b)
               this%izero(:) = this%idum(:)
            else if (m > 2) then
               call this%kubo_ham_vec_matmul(w1, w2, work, a=a, b=b)
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
         ! Multiply with the input operator (slot b)
         call this%apply_kubo_operator('b', psiref, v0, work)
         this%izero(:) = this%idum(:)
         do n=1, this%control%cond_ll
            call show_progress(n, this%control%cond_ll)
            if (n == 1) then
               v1(:, :, :) = v0(:, :, :)
            else if (n == 2) then
               v0(:, :, :) = v1(:, :, :)
               call this%kubo_ham_vec_matmul(v0, v1, work, a=a, b=b)
               this%izero(:) = this%idum(:)
            else if (n > 2) then
               call this%kubo_ham_vec_matmul(v1, v2, work, a=a, b=b)
               this%izero(:) = this%idum(:)
               v2(:, :, :) = 2 * v2(:, :, :) - v0(:, :, :)
               v0(:, :, :) = v1(:, :, :)
               v1(:, :, :) = v2(:, :, :)
               v2(:, :, :) = (0.0d0, 0.0d0)
            end if
            ! Multiply with the output operator (slot a). Its output map belongs to
            ! right_vec; keep izero for the next Chebyshev step of v1.
            call this%apply_kubo_operator('a', v1, right_vec, work)
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
         call g_kpm_profile%note_trace_batch(1)
         call g_kpm_profile%stop('P_moments_total')
      end do

      deallocate(psiref, w0, w1, w2, right_vec, v0, v1, v2, left_vec, work)

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

      this%b2_b_is_sqrt = .false.

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
                  prec=this%gpu_precision_code(block=.true.))

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
      if (trim(this%control%cheb_backend) /= 'legacy' .and. .not. this%hamiltonian%local_axis) then
         ! With hoh, ccor_2c is an extra additive term in the kernel; the folded
         ! ee_ccor_work blocks are for the no-hoh route only.
         if (this%hamiltonian%ccor_2c .and. .not. this%hamiltonian%hoh) call ensure_ccor_operator_blocks(this)
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
               if (this%hamiltonian%hoh .and. this%hamiltonian%ccor_2c) then
                  call block_lanczos_fast(this%psi_b, llmax, this%atemp_b, this%b2temp_b, &
                     this%hamiltonian%ee, this%hamiltonian%hall, this%hamiltonian%lsham, &
                     this%lattice%nn, this%lattice%iz, this%lattice%kk, nb, size(this%lattice%nn, 2), &
                     this%lattice%ntype, this%lattice%nmax, &
                     trim(this%control%cheb_backend) == 'fast', .true., &
                     this%hamiltonian%eeo, this%hamiltonian%hallo, this%hamiltonian%enim, &
                     this%hamiltonian%eecc, this%hamiltonian%hallcc)
               else if (this%hamiltonian%hoh) then
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
         if (trim(this%control%cheb_backend) /= 'legacy') then
            if (this%hamiltonian%hoh .and. this%hamiltonian%ccor_2c) then
               call cheb_moments_orbital_fast(left_vec, psiref, this%hamiltonian%ee, &
                  this%hamiltonian%hall, this%hamiltonian%lsham, this%lattice%nn, &
                  this%lattice%iz, this%lattice%kk, nb, size(this%lattice%nn, 2), &
                  this%lattice%ntype, this%lattice%nmax, this%control%lld, a, b, mu_tmp, &
                  .true., this%hamiltonian%eeo, this%hamiltonian%hallo, this%hamiltonian%enim, &
                  this%hamiltonian%eecc, this%hamiltonian%hallcc)
            else if (this%hamiltonian%hoh) then
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
