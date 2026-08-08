!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! SUBMODULE: reciprocal_spin_density
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> The k-space producer for the shared rotating-frame density contract
!> (GBT completion blueprint WP7 / gate G7).
!>
!> This is the reciprocal-space twin of `bands%accumulate_spin_density_rs`:
!> it fills the SAME `spin_density` object, per site and per angular channel,
!> straight from the eigenvectors, k weights and Fermi-Dirac occupations --
!>
!>   rho_al(s,s') = sum_k w_k sum_n f(eps_nk) E_nk^p sum_{m in l}
!>                  c_{a l m s, n k} conjg(c_{a l m s', n k}),
!>
!> with p = 0, 1, 2 giving the three energy-moment orders the radial SCF
!> consumes. There is no energy grid, no broadening and, crucially, no
!> projection axis: the up/down channels are formed afterwards, in
!> `fill_band_moments_from_spin_density`, from the object's own explicit axis.
!>
!> Before WP7 the k-space route did the opposite: `project_dos_orbitals_*`
!> projected each eigenstate onto `potential%mom` (via `get_site_spin_axis`)
!> while building a broadened DOS, `calculate_band_moments` integrated that DOS
!> on the output grid, and the Cartesian moment was rebuilt separately
!> afterwards. Three different objects, one implicit axis, no way to compare the
!> result with the real-space route. The projected DOS remains -- as DOS output
!> -- but no longer defines the SCF density.
!>
!> Conventions: rho_{s s'} = psi_s conjg(psi_s'), so
!> m_x = 2 Re(conjg(u) d), m_y = 2 Im(conjg(u) d), m_z = |u|^2 - |d|^2, matching
!> both the pre-WP7 spinor projections and the real-space producer. See the
!> `spin_density_mod` header for the full contract.
!------------------------------------------------------------------------------
submodule(reciprocal_mod) reciprocal_spin_density
   use spin_density_mod, only: spin_density, sd_orders, sd_producer_kspace, &
                               sd_constrained_spiral
   implicit none

contains

   !> @brief Fill `this%rf_density` from eigenvectors, k weights and occupations.
   !> @details Requires an existing eigensystem (built and diagonalized by the
   !>          caller) and the canonical Fermi level. The object is sized
   !>          (n_sites, n_orb_types) and zeroed here; axes are NOT set -- the
   !>          caller states them through the SCF policy.
   module subroutine accumulate_spin_density_kspace(this)
      class(reciprocal), intent(inout) :: this

      integer :: ik, ik_global, ib, isite, iorb, io, iup, idn
      integer :: nbands, n_orb_site, n_orb_per_spin, site_orb_start
      integer :: iorder
      integer :: lstart(4), lend(4)
      real(rp) :: wk, occ, e, kT, farg, epow
      complex(rp) :: u, d
      complex(rp) :: blk(2, 2)
      real(rp), parameter :: kB_Ry_per_K = 6.3336814e-6_rp
#ifdef USE_MPI
      real(rp), allocatable :: rho_comm(:)
      integer :: nflat
#endif

      if (.not. allocated(this%eigenvectors) .or. .not. allocated(this%eigenvalues)) then
         call g_logger%fatal('reciprocal%accumulate_spin_density_kspace: no eigensystem; ' // &
            'build_kspace_hamiltonian/diagonalize_hamiltonian must run first.', __FILE__, __LINE__)
      end if

      this%n_sites = this%lattice%nrec
      this%n_orb_types = 4
      this%n_spin_components = 2
      nbands = size(this%eigenvalues, 1)

      if (mod(size(this%eigenvectors, 1), this%n_sites) /= 0) then
         call g_logger%fatal('reciprocal%accumulate_spin_density_kspace: eigenvector basis ' // &
            'size not divisible by the number of sites.', __FILE__, __LINE__)
      end if
      n_orb_site = size(this%eigenvectors, 1)/this%n_sites
      n_orb_per_spin = n_orb_site/this%n_spin_components

      this%rf_density = spin_density(this%n_sites, this%n_orb_types)
      this%rf_density%producer = sd_producer_kspace

      lstart = [1, 2, 5, 10]
      lend = [1, 4, 9, 16]
      kT = this%temperature*kB_Ry_per_K

      do ik = 1, size(this%eigenvalues, 2)
         ik_global = local_k_index_to_global(this, ik)
         wk = this%k_weights(ik_global)
         do ib = 1, nbands
            e = this%eigenvalues(ib, ik)
            if (kT > 1.0e-12_rp) then
               farg = (e - this%fermi_level)/kT
               if (farg > 50.0_rp) then
                  occ = 0.0_rp
               else if (farg < -50.0_rp) then
                  occ = 1.0_rp
               else
                  occ = 1.0_rp/(exp(farg) + 1.0_rp)
               end if
            else
               occ = merge(1.0_rp, 0.0_rp, e <= this%fermi_level)
            end if
            if (occ <= 1.0e-14_rp) cycle

            do isite = 1, this%n_sites
               site_orb_start = (isite - 1)*n_orb_site
               do iorb = 1, this%n_orb_types
                  if (lstart(iorb) > n_orb_per_spin) cycle
                  blk(:, :) = (0.0_rp, 0.0_rp)
                  do io = lstart(iorb), min(lend(iorb), n_orb_per_spin)
                     iup = site_orb_start + io
                     idn = site_orb_start + n_orb_per_spin + io
                     if (idn > size(this%eigenvectors, 1)) cycle
                     u = this%eigenvectors(iup, ib, ik)
                     d = this%eigenvectors(idn, ib, ik)
                     blk(1, 1) = blk(1, 1) + u*conjg(u)
                     blk(2, 2) = blk(2, 2) + d*conjg(d)
                     blk(1, 2) = blk(1, 2) + u*conjg(d)
                     blk(2, 1) = blk(2, 1) + d*conjg(u)
                  end do
                  epow = 1.0_rp
                  do iorder = 1, sd_orders
                     call this%rf_density%accumulate_block(isite, iorb, iorder, &
                                                           blk*cmplx(wk*occ*epow, 0.0_rp, rp))
                     epow = epow*e
                  end do
               end do
            end do
         end do
      end do

#ifdef USE_MPI
      if (this%k_mesh_distributed_active) then
         nflat = 2*size(this%rf_density%rho)
         allocate (rho_comm(nflat))
         rho_comm(1:nflat:2) = real(reshape(this%rf_density%rho, [size(this%rf_density%rho)]), rp)
         rho_comm(2:nflat:2) = aimag(reshape(this%rf_density%rho, [size(this%rf_density%rho)]))
         call MPI_ALLREDUCE(MPI_IN_PLACE, rho_comm, nflat, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         this%rf_density%rho = reshape(cmplx(rho_comm(1:nflat:2), rho_comm(2:nflat:2), rp), &
                                       shape(this%rf_density%rho))
         deallocate (rho_comm)
      end if
#endif
   end subroutine accumulate_spin_density_kspace

   !> @brief Fill `band_moments` by projecting the accumulated density.
   !> @details `band_moments(:,:,:,1..3)` keeps its established meaning --
   !>          occupation, band centre <E>, and rms width sqrt(<(E-<E>)^2>) --
   !>          but now comes from the shared contract, projected after
   !>          accumulation against the explicit per-site axis, instead of from
   !>          a `potential%mom`-projected DOS grid.
   !> @param[in] policy    SCF density policy governing the axis choice.
   !> @param[in] reference Per-site reference directions (3, n_sites).
   !> @param[out] axis_out Per-site axis actually used (3, n_sites).
   module subroutine fill_band_moments_from_spin_density(this, policy, reference, axis_out)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: policy
      real(rp), dimension(:, :), intent(in) :: reference
      real(rp), dimension(:, :), intent(out) :: axis_out

      integer :: isite, iorb, ispin
      real(rp) :: axis(3), ref(3), m_transverse(3), torque(3)
      real(rp) :: site_charge, m_long
      real(rp) :: q0, centre, q2
      logical :: ok
      character(len=256) :: message

      if (this%rf_density%nsite /= this%n_sites) then
         call g_logger%fatal('reciprocal%fill_band_moments_from_spin_density: the density ' // &
            'contract was not filled for this eigensystem.', __FILE__, __LINE__)
      end if

      this%rf_density%policy = trim(policy)

      do isite = 1, this%n_sites
         ref(:) = reference(:, isite)
         if (sqrt(sum(ref(:)**2)) <= 1.0e-12_rp) ref(:) = [0.0_rp, 0.0_rp, 1.0_rp]
         call this%rf_density%resolve_site_axis(isite, ref, axis, site_charge, m_long, &
                                                m_transverse, torque)
         call this%rf_density%set_axis(isite, axis)
         axis_out(:, isite) = axis(:)

         ! Same per-site policy report the real-space route emits, so the two
         ! routes can be diffed line for line.
         if (trim(this%rf_density%policy) == sd_constrained_spiral) then
            call g_logger%info('DENSITY_POLICY constrained_spiral atom'//fmt('i4', isite)// &
                               ' m_long='//fmt('f10.6', m_long)// &
                               ' |m_transverse|='//fmt('f10.6', sqrt(sum(m_transverse(:)**2)))// &
                               ' torque=('//fmt('f10.6', torque(1))//','// &
                               fmt('f10.6', torque(2))//','//fmt('f10.6', torque(3))//')', &
                               __FILE__, __LINE__)
         else
            call g_logger%info('DENSITY_POLICY relaxed_reference atom'//fmt('i4', isite)// &
                               ' |m|='//fmt('f10.6', m_long)// &
                               ' axis=('//fmt('f10.6', axis(1))//','// &
                               fmt('f10.6', axis(2))//','//fmt('f10.6', axis(3))//')', &
                               __FILE__, __LINE__)
         end if
      end do

      call this%rf_density%check_physicality(1.0e-6_rp, ok, message)
      if (.not. ok) then
         call g_logger%fatal('reciprocal%fill_band_moments_from_spin_density: rotating-frame ' // &
            'density failed a physicality assertion -- '//trim(message)//'.', __FILE__, __LINE__)
      end if

      if (allocated(this%band_moments)) deallocate (this%band_moments)
      allocate (this%band_moments(this%n_sites, this%n_orb_types, this%n_spin_components, 3))
      this%band_moments = 0.0_rp

      do isite = 1, this%n_sites
         do iorb = 1, this%n_orb_types
            do ispin = 1, this%n_spin_components
               call this%rf_density%radial_band_moments(isite, iorb, ispin, q0, centre, q2)
               this%band_moments(isite, iorb, ispin, 1) = q0
               this%band_moments(isite, iorb, ispin, 2) = centre
               if (abs(q0) > epsilon) then
                  this%band_moments(isite, iorb, ispin, 3) = sqrt(max(q2/q0, 0.0_rp))
               else
                  this%band_moments(isite, iorb, ispin, 3) = 0.0_rp
               end if
            end do
         end do
      end do

      call g_logger%info('reciprocal%fill_band_moments_from_spin_density: band moments ' // &
         'projected from the shared rotating-frame density (policy "'//trim(policy)// &
         '", electron count '//trim(real2str(this%rf_density%electron_count(), '(F12.6)'))//').', &
         __FILE__, __LINE__)
   end subroutine fill_band_moments_from_spin_density

end submodule reciprocal_spin_density
