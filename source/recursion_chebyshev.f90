submodule(recursion_mod) recursion_chebyshev

contains

   !> @brief Dispatch Chebyshev moment generation to the configured CPU backend.
   !> @details Builds moments mu_n = <psi0|T_n((H-b)/a)|psi0> for the real-space
   !>          Chebyshev SCF and DOS workflows. The selector honors the requested
   !>          fast, batched, MKL, or legacy backend, with fallbacks for HOH and
   !>          combined ccor_2c/HOH cases that are not implemented in every kernel.
   !> @param[inout] this Recursion object; may populate ccor operator work arrays.
   !> @param[in] psi0 Starting block state, stored site-major as (nb,nb,nsite).
   !> @param[in] lld Number of Chebyshev recursion steps requested.
   !> @param[in] a Hamiltonian scaling half-width.
   !> @param[in] b Hamiltonian scaling center.
   !> @param[inout] mu Moment tensor (nb,nb,nmoment) for one starting state.
   !> @note Logs the selected backend and may reuse fast-kernel caches.
   module subroutine cheb_moments_cpu(this, psi0, lld, a, b, mu)
      class(recursion), intent(inout) :: this
      complex(rp), intent(in) :: psi0(:, :, :)
      integer, intent(in) :: lld
      real(rp), intent(in) :: a, b
      complex(rp), intent(inout) :: mu(:, :, :)
      character(len=:), allocatable :: backend
      character(len=256) :: msg

      backend = trim(this%control%cheb_backend)

      ! 'fast_dp' is the fp64 selector for the block (Haydock) fast kernels; the
      ! Chebyshev fast kernel is single precision only, so treat it as 'fast'.
      if (backend == 'fast_dp') backend = 'fast'

      ! The orthogonalisation (hoh) term is now supported natively by the 'fast'
      ! kernel via a genuine two-sweep apply (h*x then eeo*(h*x)), mirroring the
      ! legacy recursion. The BSR-based kernels (batched / mkl_batch / mkl_sparse)
      ! do not yet carry the second sweep, so fall those back to legacy when hoh
      ! is active. 'legacy' stays available as an explicit choice.
      if (this%hamiltonian%hoh .and. backend /= 'fast' .and. backend /= 'legacy') then
         backend = 'legacy'
      end if
      if (this%hamiltonian%ccor_2c .and. this%hamiltonian%hoh .and. backend /= 'legacy') then
         call g_logger%warning('ccor_2c with hoh is not supported by fast Chebyshev backends; falling back to legacy.', __FILE__, __LINE__)
         backend = 'legacy'
      end if
      if (this%hamiltonian%ccor_2c .and. backend /= 'legacy') call ensure_ccor_operator_blocks(this)
      msg = 'Chebyshev moment driver backend: '//backend
      !print *, msg
      call g_logger%info(msg, __FILE__, __LINE__)
      if (this%hamiltonian%ccor_2c .and. backend /= 'legacy') then
         associate(ee_op => this%ee_ccor_work, hall_op => this%hall_ccor_work)
            call dispatch_cheb_backend(ee_op, hall_op)
         end associate
      else
         associate(ee_op => this%hamiltonian%ee, hall_op => this%hamiltonian%hall)
            call dispatch_cheb_backend(ee_op, hall_op)
         end associate
      end if

   contains

      subroutine dispatch_cheb_backend(ee_op, hall_op)
         complex(rp), intent(in), contiguous :: ee_op(:, :, :, :), hall_op(:, :, :, :)

         select case (backend)
         case ('fast')
            if (this%hamiltonian%hoh) then
               call cheb_moments_fast(psi0, this%hamiltonian%ee, this%hamiltonian%hall, &
                  this%hamiltonian%lsham, this%lattice%nn, this%lattice%iz, this%lattice%kk, &
                  nb, size(this%lattice%nn, 2), this%lattice%ntype, this%lattice%nmax, &
                  lld, a, b, mu, hoh=.true., eeo=this%hamiltonian%eeo, &
                  hallo=this%hamiltonian%hallo, enim=this%hamiltonian%enim)
            else
               call cheb_moments_fast(psi0, ee_op, hall_op, &
                  this%hamiltonian%lsham, this%lattice%nn, this%lattice%iz, this%lattice%kk, &
                  nb, size(this%lattice%nn, 2), this%lattice%ntype, this%lattice%nmax, &
                  lld, a, b, mu)
            end if
         case ('batched')
            call cheb_moments_fast_batched(psi0, ee_op, hall_op, &
               this%hamiltonian%lsham, this%lattice%nn, this%lattice%iz, this%lattice%kk, &
               nb, size(this%lattice%nn, 2), this%lattice%ntype, this%lattice%nmax, &
               lld, a, b, mu)
         case ('mkl_batch')
            call cheb_moments_fast_mkl_batch(psi0, ee_op, hall_op, &
               this%hamiltonian%lsham, this%lattice%nn, this%lattice%iz, this%lattice%kk, &
               nb, size(this%lattice%nn, 2), this%lattice%ntype, this%lattice%nmax, &
               lld, a, b, mu)
         case ('mkl_sparse')
            call cheb_moments_fast_mkl_sparse(psi0, ee_op, hall_op, &
               this%hamiltonian%lsham, this%lattice%nn, this%lattice%iz, this%lattice%kk, &
               nb, size(this%lattice%nn, 2), this%lattice%ntype, this%lattice%nmax, &
               lld, a, b, mu)
         case ('legacy')
            call cheb_moments_legacy(this, psi0, lld, a, b, mu)
         case default
            call g_logger%fatal('Unknown cheb_backend "'//backend//'"', __FILE__, __LINE__)
         end select
      end subroutine dispatch_cheb_backend
   end subroutine cheb_moments_cpu

   !> @brief Compute Chebyshev moments with the original light-cone recursion.
   !> @details Applies the three-term Chebyshev recurrence explicitly through
   !>          cheb_0th_mom, cheb_1st_mom[_hoh], and chebyshev_recur_ll[_hoh].
   !>          This path remains the reference implementation for cases not
   !>          covered by the optimized CPU or GPU backends.
   !> @param[inout] this Recursion object holding Hamiltonian, light-cone maps,
   !>                temporary wavefunctions, and output moment buffers.
   !> @param[in] psi0 Starting block state, stored site-major as (nb,nb,nsite).
   !> @param[in] lld Number of Chebyshev recursion steps requested.
   !> @param[in] a Hamiltonian scaling half-width.
   !> @param[in] b Hamiltonian scaling center.
   !> @param[inout] mu Moment tensor (nb,nb,nmoment) for one starting state.
   !> @note Rebuilds this%izero from the support of psi0 and mutates psi0/psi1/psi2.
   module subroutine cheb_moments_legacy(this, psi0, lld, a, b, mu)
      class(recursion), intent(inout) :: this
      complex(rp), intent(in) :: psi0(:, :, :)
      integer, intent(in) :: lld
      real(rp), intent(in) :: a, b
      complex(rp), intent(inout) :: mu(:, :, :)
      ! Local variables
      integer :: k, ll

      ! Rebuild the initial lightcone from the support of psi0. The leaf
      ! routines propagate it dynamically through this%izero/this%idum, so the
      ! izeroll/create_ll_map bookkeeping of the old inline drivers is not
      ! needed (those arrays are never read by the leaf routines).
      this%izero(:) = 0
      do k = 1, this%lattice%kk
         if (any(psi0(:, :, k) /= (0.0_rp, 0.0_rp))) this%izero(k) = 1
      end do

      ! Initialise the propagated states. psi0 also serves as the projection
      ! reference (diagonal double-trick), matching the original drivers.
      this%psi0(:, :, :) = psi0(:, :, :)
      this%psi1(:, :, :) = (0.0_rp, 0.0_rp)
      this%psi2(:, :, :) = (0.0_rp, 0.0_rp)

      ! 0th moment
      call this%cheb_0th_mom(psi0)
      mu(:, :, 1) = this%cheb_mom_temp(:, :)

      ! 1st moment
      if (this%hamiltonian%hoh) then
         call this%cheb_1st_mom_hoh(psi0, a, b)
      else
         call this%cheb_1st_mom(psi0, a, b)
      end if
      mu(:, :, 2) = this%cheb_mom_temp(:, :)

      this%izero(:) = this%idum(:)

      ! Higher moments via the double-trick recursion
      do ll = 1, lld
         if (this%hamiltonian%hoh) then
            call this%chebyshev_recur_ll_hoh(mu, ll, a, b)
         else
            call this%chebyshev_recur_ll(mu, ll, a, b)
         end if
      end do
   end subroutine cheb_moments_legacy

   !> @brief Return the energy window used to scale Chebyshev Hamiltonians.
   !> @details Uses Hamiltonian spectral bounds when they have been computed,
   !>          otherwise falls back to the namelist energy mesh. DOS and Green
   !>          reconstruction call this so that moment generation and evaluation
   !>          use the same affine map H -> (H-b)/a.
   !> @param[inout] this Recursion object with control, energy, and Hamiltonian bounds.
   !> @param[out] emin Lower edge of the scaling window.
   !> @param[out] emax Upper edge of the scaling window.
   !> @note Emits an informational log describing the chosen source of bounds.
   module subroutine resolve_chebyshev_window(this, emin, emax)
      class(recursion), intent(inout) :: this
      real(rp), intent(out) :: emin, emax
      character(len=256) :: msg

      emin = this%en%energy_min
      emax = this%en%energy_max

      if (this%control%recur /= 'chebyshev') return

      if (this%hamiltonian%bounds%e_max > this%hamiltonian%bounds%e_min) then
         emin = this%hamiltonian%bounds%e_min
         emax = this%hamiltonian%bounds%e_max
         msg = 'Chebyshev scaling window from Hamiltonian bounds: ['// &
               trim(fmt('f12.6', emin))//','//trim(fmt('f12.6', emax))//']'
         call g_logger%info(msg, __FILE__, __LINE__)
      else
         msg = 'Chebyshev Hamiltonian bounds from energy namelist window ['// &
               trim(fmt('f12.6', emin))//','//trim(fmt('f12.6', emax))//']'
         call g_logger%info(msg, __FILE__, __LINE__)
      end if
   end subroutine resolve_chebyshev_window

   !> @brief Accumulate the zeroth Chebyshev moment for one reference block.
   !> @details Computes <psiref|T_0|psiref> over the current light cone and stores
   !>          the result in cheb_mom_temp for the legacy Chebyshev driver.
   !> @param[inout] this Recursion object; writes cheb_mom_temp.
   !> @param[in] psiref Reference block state, site-major as (nb,nb,nsite).
   module subroutine cheb_0th_mom(this, psiref)
      class(recursion), intent(inout) :: this
      complex(rp), dimension(nb, nb, this%lattice%kk), intent(in) :: psiref
      complex(rp), dimension(nb, nb) :: dum
      integer :: nat, hblocksize, nlimplus1, k

      hblocksize = nb
      nat = this%lattice%kk
      nlimplus1 = this%lattice%nmax + 1

      this%cheb_mom_temp(:, :) = 0.0d0
      ! Write the 0th moment
      do k = 1, nat
         if (this%izero(k) /= 0) then
            call zgemm('c', 'n', nb, nb, nb, cone, psiref(:, :, k), nb, this%psi0(:, :, k), nb, cone, this%cheb_mom_temp, nb)
         end if
      end do
   end subroutine cheb_0th_mom

   !> @brief Accumulate the first Chebyshev moment without HOH correction.
   !> @details Applies the scaled Hamiltonian (H-b)/a once to the reference block
   !>          and stores <psiref|T_1|psiref> in cheb_mom_temp. This seeds the
   !>          legacy three-term Chebyshev recurrence.
   !> @param[inout] this Recursion object; mutates psi1, idum, and cheb_mom_temp.
   !> @param[in] psiref Reference block state, site-major as (nb,nb,nsite).
   !> @param[in] a Hamiltonian scaling half-width.
   !> @param[in] b Hamiltonian scaling center.
   module subroutine cheb_1st_mom(this, psiref, a, b)
      class(recursion), intent(inout) :: this
      real(rp), intent(in) :: a, b ! scale and shift
      complex(rp), dimension(nb, nb, this%lattice%kk), intent(in) :: psiref
      integer :: nat, hblocksize, nlimplus1, k, ih, nr, n, ineigh, nnmap, i

      hblocksize = nb
      nat = this%lattice%kk
      nlimplus1 = this%lattice%nmax + 1

      this%cheb_mom_temp(:, :) = 0.0d0

      ! Write |phi_1>=H|phi_0>
      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         do i = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(i) = this%izero(i)
            ih = this%lattice%iz(i)
            nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
            if (this%izero(i) /= 0) then
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, 1, i), nb, this%psi0(:, :, i), nb, cone, this%psi1(:, :, i), nb)
               if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallcc(:, :, 1, i), nb, this%psi0(:, :, i), nb, cone, this%psi1(:, :, i), nb)
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%lsham(:, :, ih), nb, this%psi0(:, :, i), nb, cone, this%psi1(:, :, i), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, ineigh, i), nb, this%psi0(:, :, nnmap), nb, cone, this%psi1(:, :, i), nb)
                        if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallcc(:, :, ineigh, i), nb, this%psi0(:, :, nnmap), nb, cone, this%psi1(:, :, i), nb)
                        this%idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
            ! Do the scaling and shifting
            this%psi1(:, :, i) = this%psi1(:, :, i) - b*this%psi0(:, :, i)
            this%psi1(:, :, i) = this%psi1(:, :, i)/a
         end do ! End of loop in the neighbouring
      end if ! End of local Hamiltonian loop

      ! Write |phi_1>=H|phi_0>
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(1, 1, 1, ih), nb, this%psi0(:, :, k), nb, cone, this%psi1(:, :, k), nb)
            if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eecc(:, :, 1, ih), nb, this%psi0(:, :, k), nb, cone, this%psi1(:, :, k), nb)
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%lsham(:, :, ih), nb, this%psi0(:, :, k), nb, cone, this%psi1(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(1, 1, ineigh, ih), nb, this%psi0(:, :, nnmap), nb, cone, this%psi1(:, :, k), nb)
                  if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eecc(:, :, ineigh, ih), nb, this%psi0(:, :, nnmap), nb, cone, this%psi1(:, :, k), nb)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
         ! Do the scaling and shifting
         this%psi1(:, :, k) = this%psi1(:, :, k) - b*this%psi0(:, :, k)
         this%psi1(:, :, k) = this%psi1(:, :, k)/a
      end do ! End loop in the clust

      ! Write the 1st moment
      do n = 1, this%lattice%kk
         if (this%izero(n) /= 0) then
            call zgemm('c', 'n', nb, nb, nb, cone, psiref(:, :, n), nb, this%psi1(:, :, n), nb, cone, this%cheb_mom_temp, nb)
         end if
      end do
   end subroutine cheb_1st_mom

   !> @brief Accumulate the first Chebyshev moment with HOH correction.
   !> @details Applies the scaled h - H O H + e_nu + l.s operator once to the
   !>          reference block and stores the T_1 projection in cheb_mom_temp.
   !> @param[inout] this Recursion object; uses psi1/psi2/hoh scratch and maps.
   !> @param[in] psiref Reference block state, site-major as (nb,nb,nsite).
   !> @param[in] a Hamiltonian scaling half-width.
   !> @param[in] b Hamiltonian scaling center.
   module subroutine cheb_1st_mom_hoh(this, psiref, a, b)
      class(recursion), intent(inout) :: this
      real(rp), intent(in) :: a, b ! scale and shift
      complex(rp), dimension(nb, nb, this%lattice%kk), intent(in) :: psiref
      integer :: nat, hblocksize, nlimplus1, k, ih, nr, n, ineigh, nnmap, i

      hblocksize = nb
      nat = this%lattice%kk
      nlimplus1 = this%lattice%nmax + 1

      this%cheb_mom_temp(:, :) = 0.0d0

      this%hohpsi(:, :, :) = (0.0d0, 0.0d0)
      this%enupsi(:, :, :) = (0.0d0, 0.0d0)
      this%socpsi(:, :, :) = (0.0d0, 0.0d0)

      ! Write |phi_1>=H|phi_0>
      nlimplus1 = this%lattice%nmax + 1
      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         do i = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(i) = this%izero(i)
            ih = this%lattice%iz(i)
            nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
            if (this%izero(i) /= 0) then
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, 1, i), nb, this%psi0(:, :, i), nb, cone, this%psi1(:, :, i), nb)
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%lsham(:, :, ih), nb, this%psi0(:, :, i), nb, cone, this%socpsi(:, :, i), nb)
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%enim(:, :, ih), nb, this%psi0(:, :, i), nb, cone, this%enupsi(:, :, i), nb)
               if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallcc(:, :, 1, i), nb, this%psi0(:, :, i), nb, cone, this%enupsi(:, :, i), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, ineigh, i), nb, this%psi0(:, :, nnmap), nb, cone, this%psi1(:, :, i), nb)
                        if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallcc(:, :, ineigh, i), nb, this%psi0(:, :, nnmap), nb, cone, this%enupsi(:, :, i), nb)
                        this%idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
      end if ! End of local Hamiltonian loop

      ! Write |phi_1>=H|phi_0>
      do k = nlimplus1, this%lattice%kk ! Loop in the clust
         this%idum(k) = this%izero(k)
         ih = this%lattice%iz(k)
         nr = this%lattice%nn(k, 1)
         if (this%izero(k) /= 0) then
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, 1, ih), nb, this%psi0(:, :, k), nb, cone, this%psi1(:, :, k), nb)
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%enim(:, :, ih), nb, this%psi0(:, :, k), nb, cone, this%enupsi(:, :, k), nb)
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%lsham(:, :, ih), nb, this%psi0(:, :, k), nb, cone, this%socpsi(:, :, k), nb)
            if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eecc(:, :, 1, ih), nb, this%psi0(:, :, k), nb, cone, this%enupsi(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih), nb, this%psi0(:, :, nnmap), nb, cone, this%psi1(:, :, k), nb)
                  if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eecc(:, :, ineigh, ih), nb, this%psi0(:, :, nnmap), nb, cone, this%enupsi(:, :, k), nb)
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
         ! Do the scaling and shifting
         !this%psi1(:, :, k) = this%psi1(:, :, k) - b*this%psi0(:, :, k)
         !this%psi1(:, :, k) = this%psi1(:, :, k)/a
      end do ! End loop in the clust

      ! Mapping update for hoh calculation
      this%izero(:) = this%idum(:)

      if (this%lattice%nmax /= 0) then ! In case of impurities using the local hamiltonian
         do i = 1, this%lattice%nmax ! Loop in the neighbouring
            this%idum(i) = this%izero(i)
            ih = this%lattice%iz(i)
            nr = this%lattice%nn(i, 1) ! Number of neighbours of atom i
            if (this%izero(i) /= 0) then
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallo(:, :, 1, i), nb, this%psi1(:, :, i), nb, cone, this%hohpsi(:, :, i), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(i, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallo(:, :, ineigh, i), nb, this%psi1(:, :, nnmap), nb, cone, this%hohpsi(:, :, i), nb)
                        this%idum(i) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
         end do ! End of loop in the neighbouring
      end if ! End of local Hamiltonian loop

      do i = nlimplus1, this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
         this%idum(i) = this%izero(i)
         ih = this%lattice%iz(i) ! Atom type
         nr = this%lattice%nn(i, 1) ! Number of neighbours
         if (this%izero(i) /= 0) then
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eeo(:, :, 1, ih), nb, this%psi1(:, :, i), nb, cone, this%hohpsi(:, :, i), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr ! Loop on the neighbouring
               nnmap = this%lattice%nn(i, ineigh)
               if (nnmap /= 0) then
                  if (this%izero(nnmap) /= 0) then
                     call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eeo(:, :, ineigh, ih), nb, this%psi1(:, :, nnmap), nb, cone, this%hohpsi(:, :, i), nb)
                     this%idum(i) = 1
                  end if
               end if
            end do ! End of loop in the neighbouring
         end if
      end do

      ! H = h - hoh + e_nu + l.s
      this%psi1(:, :, :) = this%psi1(:, :, :) - this%hohpsi(:, :, :) + this%enupsi(:, :, :) + this%socpsi(:, :, :)

      ! Do the scaling and shifting
      this%psi1(:, :, :) = this%psi1(:, :, :) - b*this%psi0(:, :, :)
      this%psi1(:, :, :) = this%psi1(:, :, :)/a

      ! Write the 1st moment
      do n = 1, this%lattice%kk
         if (this%izero(n) /= 0) then
            call zgemm('c', 'n', nb, nb, nb, cone, psiref(:, :, n), nb, this%psi1(:, :, n), nb, cone, this%cheb_mom_temp, nb)
         end if
      end do
   end subroutine cheb_1st_mom_hoh

   !> @brief Advance one no-HOH legacy Chebyshev recursion level.
   !> @details Applies the scaled Hamiltonian and the double-trick recurrence to
   !>          produce the next pair of moments for a single starting block.
   !> @param[inout] this Recursion object; mutates psi0/psi1/psi2 and light-cone maps.
   !> @param[inout] mu Moment tensor for the active starting block.
   !> @param[in] ll Current recursion level.
   !> @param[in] a Hamiltonian scaling half-width.
   !> @param[in] b Hamiltonian scaling center.
   module subroutine chebyshev_recur_ll(this, mu, ll, a, b)
      class(recursion), intent(inout) :: this
      complex(rp), intent(inout) :: mu(:, :, :)
      integer, intent(in) :: ll
      real(rp), intent(in) :: a, b
      ! Local variables
      integer :: ineigh, ih, j, k, nr, m, n, l, hblocksize, nat, nnmap, nlimplus1
      complex(rp), dimension(nb, nb) :: dum1, dum2, locham

      hblocksize = nb
      nat = this%lattice%kk

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
               call zgemm('n', 'n', nb, nb, nb, cone, locham, nb, this%psi1(:, :, k), nb, cone, this%psi2(:, :, k), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        if (this%hamiltonian%ccor_2c) then
                           call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, ineigh, k) + this%hamiltonian%hallcc(:, :, ineigh, k), nb, this%psi1(:, :, nnmap), nb, cone, this%psi2(:, :, k), nb)
                        else
                           call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, ineigh, k), nb, this%psi1(:, :, nnmap), nb, cone, this%psi2(:, :, k), nb)
                        end if
                        this%idum(k) = 1
                     end if
                  end if
               end do ! End of loop in the neighbouring
            end if
            ! Do the scaling and shifting
            this%psi2(:, :, k) = this%psi2(:, :, k) - b*this%psi1(:, :, k)
            this%psi2(:, :, k) = this%psi2(:, :, k)/a
            ! Write 2*H*|phi_1>
            this%psi2(:, :, k) = 2*this%psi2(:, :, k)
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
            call zgemm('n', 'n', nb, nb, nb, cone, locham, nb, this%psi1(:, :, k), nb, cone, this%psi2(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  if (this%hamiltonian%ccor_2c) then
                     call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih) + this%hamiltonian%eecc(:, :, ineigh, ih), nb, this%psi1(:, :, nnmap), nb, cone, this%psi2(:, :, k), nb)
                  else
                     call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih), nb, this%psi1(:, :, nnmap), nb, cone, this%psi2(:, :, k), nb)
                  end if
                  this%idum(k) = 1
               end if
            end do ! End of the loop in the neighbouring
         end if
         ! Do the scaling and shifting
         this%psi2(:, :, k) = this%psi2(:, :, k) - b*this%psi1(:, :, k)
         this%psi2(:, :, k) = this%psi2(:, :, k)/a
         ! Write 2*H*|phi_1>
         this%psi2(:, :, k) = 2*this%psi2(:, :, k)
      end do ! End loop in the clust
      !$omp end parallel do

      ! Write Moments
      dum1(:, :) = (0.0d0, 0.0d0)
      dum2(:, :) = (0.0d0, 0.0d0)

      ! Redefines idum for all clust atoms
      this%irnum = 0
      do n = 1, this%lattice%kk
         this%izero(n) = this%idum(n)
         if (this%izero(n) /= 0) then
            this%irnum = this%irnum + 1
            this%irlist(this%irnum) = n
         end if
      end do

      !$omp parallel do default(shared) private(n) reduction(+:dum1) reduction(+:dum2)
      do n = 1, this%irnum
         ! Write |phi_2>=2*H*|phi_1> - |phi_0>
         this%psi2(:, :, this%irlist(n)) = this%psi2(:, :, this%irlist(n)) - this%psi0(:, :, this%irlist(n))
         call zgemm('c', 'n', nb, nb, nb, cone, this%psi1(:, :, this%irlist(n)), nb, this%psi1(:, :, this%irlist(n)), nb, cone, dum1, nb)
         call zgemm('c', 'n', nb, nb, nb, cone, this%psi2(:, :, this%irlist(n)), nb, this%psi1(:, :, this%irlist(n)), nb, cone, dum2, nb)
         this%psi0(:, :, this%irlist(n)) = this%psi1(:, :, this%irlist(n))
         this%psi1(:, :, this%irlist(n)) = this%psi2(:, :, this%irlist(n))
         this%psi2(:, :, this%irlist(n)) = (0.0d0, 0.0d0)
      end do
      !$omp end parallel  do

      mu(:, :, 2*ll + 1) = 2.0_rp*dum1(:, :) - mu(:, :, 1)
      mu(:, :, 2*ll + 2) = 2.0_rp*dum2(:, :) - mu(:, :, 2)

      if (sum(real(mu(:, :, 2*ll + 2))) > 1000.d0) then
         call g_logger%fatal('Chebyshev moments did not converge. Check energy limits energy_min and energy_max', __FILE__, __LINE__)
      end if
   end subroutine chebyshev_recur_ll

   !> @brief Advance one HOH-aware legacy Chebyshev recursion level.
   !> @details Uses the two-sweep orthogonalized Hamiltonian action in the
   !>          Chebyshev three-term recurrence, producing the next moments for
   !>          a single starting block.
   !> @param[inout] this Recursion object; mutates psi0/psi1/psi2/hohpsi and maps.
   !> @param[inout] mu Moment tensor for the active starting block.
   !> @param[in] ll Current recursion level.
   !> @param[in] a Hamiltonian scaling half-width.
   !> @param[in] b Hamiltonian scaling center.
   module subroutine chebyshev_recur_ll_hoh(this, mu, ll, a, b)
      class(recursion), intent(inout) :: this
      complex(rp), intent(inout) :: mu(:, :, :)
      integer, intent(in) :: ll
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
               call zgemm('n', 'n', nb, nb, nb, cone, locham, nb, this%psi1(:, :, k), nb, cone, this%psi2(:, :, k), nb)
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%lsham(:, :, ih), nb, this%psi1(:, :, k), nb, cone, this%socpsi(:, :, k), nb)
               call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%enim(:, :, ih), nb, this%psi1(:, :, k), nb, cone, this%enupsi(:, :, k), nb)
               if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallcc(:, :, 1, k), nb, this%psi1(:, :, k), nb, cone, this%enupsi(:, :, k), nb)
            end if
            if (nr >= 2) then
               do ineigh = 2, nr ! Loop on the neighbouring
                  nnmap = this%lattice%nn(k, ineigh)
                  if (nnmap /= 0) then
                     if (this%izero(nnmap) /= 0) then
                        call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hall(:, :, ineigh, k), nb, this%psi1(:, :, nnmap), nb, cone, this%psi2(:, :, k), nb)
                        if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%hallcc(:, :, ineigh, k), nb, this%psi1(:, :, nnmap), nb, cone, this%enupsi(:, :, k), nb)
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
            call zgemm('n', 'n', nb, nb, nb, cone, locham, nb, this%psi1(:, :, k), nb, cone, this%psi2(:, :, k), nb)
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%enim(:, :, ih), nb, this%psi1(:, :, k), nb, cone, this%enupsi(:, :, k), nb)
            call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%lsham(:, :, ih), nb, this%psi1(:, :, k), nb, cone, this%socpsi(:, :, k), nb)
            if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eecc(:, :, 1, ih), nb, this%psi1(:, :, k), nb, cone, this%enupsi(:, :, k), nb)
         end if
         if (nr >= 2) then
            do ineigh = 2, nr ! Loop in the neighbouring
               nnmap = this%lattice%nn(k, ineigh)
               if (nnmap /= 0 .and. this%izero(nnmap) /= 0) then
                  call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%ee(:, :, ineigh, ih), nb, this%psi1(:, :, nnmap), nb, cone, this%psi2(:, :, k), nb)
                  if (this%hamiltonian%ccor_2c) call zgemm('n', 'n', nb, nb, nb, cone, this%hamiltonian%eecc(:, :, ineigh, ih), nb, this%psi1(:, :, nnmap), nb, cone, this%enupsi(:, :, k), nb)
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
      this%psi2(:, :, :) = this%psi2(:, :, :) - this%hohpsi(:, :, :) + this%enupsi(:, :, :) + this%socpsi(:, :, :)

      ! Do the scaling and shifting
      this%psi2(:, :, :) = this%psi2(:, :, :) - b*this%psi1(:, :, :)
      this%psi2(:, :, :) = this%psi2(:, :, :)/a
      ! Write 2*H*|phi_1>
      this%psi2(:, :, :) = 2*this%psi2(:, :, :)

      ! Write Moments
      dum1(:, :) = (0.0d0, 0.0d0)
      dum2(:, :) = (0.0d0, 0.0d0)

      ! Redefines idum for all clust atoms
      this%irnum = 0
      do n = 1, this%lattice%kk
         this%izero(n) = this%idum(n)
         if (this%izero(n) /= 0) then
            this%irnum = this%irnum + 1
            this%irlist(this%irnum) = n
         end if
      end do

      !$omp parallel do default(shared) private(n) reduction(+:dum1) reduction(+:dum2)
      do n = 1, this%irnum
         ! Write |phi_2>=2*H*|phi_1> - |phi_0>
         this%psi2(:, :, this%irlist(n)) = this%psi2(:, :, this%irlist(n)) - this%psi0(:, :, this%irlist(n))
         call zgemm('c', 'n', nb, nb, nb, cone, this%psi1(:, :, this%irlist(n)), nb, this%psi1(:, :, this%irlist(n)), nb, cone, dum1, nb)
         call zgemm('c', 'n', nb, nb, nb, cone, this%psi2(:, :, this%irlist(n)), nb, this%psi1(:, :, this%irlist(n)), nb, cone, dum2, nb)
         this%psi0(:, :, this%irlist(n)) = this%psi1(:, :, this%irlist(n))
         this%psi1(:, :, this%irlist(n)) = this%psi2(:, :, this%irlist(n))
         this%psi2(:, :, this%irlist(n)) = (0.0d0, 0.0d0)
      end do
      !$omp end parallel  do

      mu(:, :, 2*ll + 1) = 2.0_rp*dum1(:, :) - mu(:, :, 1)
      mu(:, :, 2*ll + 2) = 2.0_rp*dum2(:, :) - mu(:, :, 2)

      if (sum(real(mu(:, :, 2*ll + 2))) > 1000.d0) then
         call g_logger%fatal('Chebyshev moments did not converge. Check energy limits energy_min and energy_max', __FILE__, __LINE__)
      end if
   end subroutine chebyshev_recur_ll_hoh

   !> @brief Generate onsite Chebyshev moments for real-space SCF and DOS.
   !> @details For each recursion atom, initializes an orbital block state and
   !>          computes mu_n moments of T_n((H-b)/a). The Green/DOS layer later
   !>          applies kernels and reconstructs spectral quantities from mu_n.
   !> @param[inout] this Recursion object; fills local slices of mu_n.
   !> @note Work is MPI-partitioned and may dispatch to the CUDA plugin or fast CPU backends.
   module subroutine chebyshev_recur(this)
      use mpi_mod
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: ineigh, ih, i, j, k, nr, ll, m, n, l, hblocksize, nat, nnmap, nlimplus1, llcheb
      integer, dimension(0:this%lattice%kk) :: idumll
      !complex(rp) :: cone = (1.0d0, 0.0d0)
      complex(rp), dimension(nb, nb) :: dum, dum1, dum2
      complex(rp), dimension(nb, this%lattice%kk) :: v
      complex(rp), dimension(:, :, :, :), allocatable :: hcheb
      real(rp) :: a, b, start, finish, emin_win, emax_win
      ! External functions
      complex(rp), external :: zdotc

      integer :: i_loc

      hblocksize = nb
      nat = this%lattice%kk
      nlimplus1 = this%lattice%nmax + 1
      llcheb = (2*this%control%lld) + 2

      call resolve_chebyshev_window(this, emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp
      call cheb_fast_reset_cache()

      if (gpu_plugin_ready(this, 'chebyshev_recur()', allow_hoh=.true.)) then
         call gpu_plugin_upload_hamiltonian(this)
         do i = start_atom, end_atom
            i_loc = g2l_map(i)
            j = this%lattice%irec(i)
            call g_logger%info('CUDA Chebyshev recursion in progress for atom '// &
               int2str(j), __FILE__, __LINE__)
            this%psi0(:, :, :) = (0.0d0, 0.0d0)
            do l = 1, nb
               this%psi0(l, l, j) = (1.0d0, 0.0d0)
            end do
            call this%gpu_backend%chebyshev_moments(this%psi0, this%control%lld, a, b, &
               this%mu_n(:, :, :, i_loc))
            call g_logger%info('CUDA Chebyshev recursion done for atom '// &
               int2str(j), __FILE__, __LINE__)
         end do
         return
      end if

      do i = start_atom, end_atom ! Loop on the number of atoms to be treat self-consistently by this process
         i_loc = g2l_map(i)
         j = this%lattice%irec(i) ! Atom number in the clust file
         call g_logger%info('Chebyshev recursion on progress for atom '//int2str(j), __FILE__, __LINE__)

         ! Starting state |phi_0>: identity block on the recursion site j.
         this%psi0(:, :, :) = (0.0d0, 0.0d0)
         do l = 1, nb
            this%psi0(l, l, j) = (1.0d0, 0.0d0)
         end do

         ! Dispatch to the selected CPU backend (fast kernels or, for hoh /
         ! cheb_backend='legacy', the wrapped legacy recursion).
         call cheb_moments_cpu(this, this%psi0, this%control%lld, a, b, &
            this%mu_n(:, :, :, i_loc))
      end do ! End loop on the number of atoms to be treat self-consistently
   end subroutine chebyshev_recur

end submodule recursion_chebyshev
