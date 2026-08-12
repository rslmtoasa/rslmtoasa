!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! SUBMODULE: reciprocal_moments
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Exact k-space double-Chebyshev moment generator for the KPM conductivity
!> pipeline (milestone B5, task B5.1). A thin consumer of the B2 k-space engine
!> that fills the SAME `mu_nm_stochastic(:,:,n,m,ntype)` moment array the
!> recursion route fills in `recursion_transport::compute_moments_stochastic`,
!> so `conductivity::calculate_conductivity_tensor` runs UNCHANGED -- filling the
!> canonical moment array IS the route-agnosticism, exactly as `fill_green` does
!> for the Green's-function consumers.
!>
!> The moment is computed EXACTLY from the k-space eigenpairs (no stochastic
!> trace estimator, no Chebyshev truncation on the operator): in the eigenbasis
!> of H(k), T_p(H~(k)) = U diag(T_p(eps~)) U^dagger, so the velocity-sandwiched
!> double moment reduces to the eigenbasis contraction `moment_onsite_block`
!> performs (see `moment_kernel_mod`). The k-space velocity operator is the
!> Fourier transform of the SAME real-space velocity blocks the recursion route
!> applies (`build_realspace_velocity_operators` -> `fourier_transform_array`),
!> so the two routes evaluate one operator in two representations and the
!> residual max|mu_exact - mu_rec| is the direct KPM error bound (B5.1 acceptance).
!>
!> Scope / preconditions (B5.1):
!>   * cond_type = 'charge' and cond_calctype = 'per_type'. The other velocity
!>     flavors (spin/orbital/torque) and the random-vector stochastic reference
!>     are the recursion route's domain; the exact k-space generator resolves
!>     moments per atom type from the on-site cell block.
!>   * first-order k-space velocity v_a(k) = sum_R v_a(R) e^{i k.R} (the strict
!>     Lehmann / S=I regime backend E lives in); consistent with hoh=false.
!>   * full, unreduced, undistributed BZ mesh (each rank holds all k), inherited
!>     from backend E; k-parallel/symmetry reduction land with B4.
!>   * a, b are the recursion route's Chebyshev window (passed in) so the moments
!>     scale identically to the ones `gamma_nm` weights.
!------------------------------------------------------------------------------
submodule(reciprocal_mod) reciprocal_moments
   use moment_kernel_mod, only: moment_onsite_block
   implicit none

contains

   !> @brief Fill mu_nm_stochastic from the k-space eigenpairs (exact moments).
   !> @details Builds and diagonalizes H(k) on the full BZ mesh (if not already
   !>          done), forms the first-order k-space velocity operators from the
   !>          recursion route's real-space velocity blocks, and fills the on-site
   !>          double-Chebyshev moment block per atom type via
   !>          `moment_onsite_block`. The output array is allocated here with the
   !>          SAME shape `compute_moments_stochastic` uses, so the conductivity
   !>          consumer is untouched.
   !> @param[inout] this Reciprocal object (H(k) machinery, k-mesh, velocities).
   !> @param[out]   mu   Allocated (nb,nb,cond_ll,cond_ll,ntype) moment array.
   !> @param[in]    a    Chebyshev window scale (recursion route's resolve_chebyshev_window).
   !> @param[in]    b    Chebyshev window shift.
   module subroutine fill_moments(this, mu, a, b)
      class(reciprocal), intent(inout) :: this
      complex(rp), allocatable, intent(out) :: mu(:, :, :, :, :)
      real(rp), intent(in) :: a, b

      integer :: nk, ik, ikg, cond_ll, ntype, nmat, site, ioff
      complex(rp), allocatable :: va_k(:, :, :), vb_k(:, :, :)

      ! --- Scope guards (B5.1): charge / per_type only. --------------------------
      if (trim(this%lattice%control%cond_type) /= 'charge') then
         call g_logger%fatal('reciprocal%fill_moments: exact k-space moments are ' // &
            'implemented for cond_type="charge" only (got "' // &
            trim(this%lattice%control%cond_type) // '").', __FILE__, __LINE__)
      end if
      if (trim(this%lattice%control%cond_calctype) /= 'per_type') then
         call g_logger%fatal('reciprocal%fill_moments: exact k-space moments require ' // &
            'cond_calctype="per_type" (got "' // &
            trim(this%lattice%control%cond_calctype) // '").', __FILE__, __LINE__)
      end if

      ! --- Eigenpairs on the full BZ mesh (one Hermitian solve per k). -----------
      if (.not. allocated(this%eigenvectors) .or. .not. allocated(this%eigenvalues)) then
         call this%build_kspace_hamiltonian()
         call this%diagonalize_hamiltonian()
      end if
      if (.not. allocated(this%eigenvectors)) then
         call g_logger%fatal('reciprocal%fill_moments: eigenvectors unavailable after ' // &
            'diagonalization.', __FILE__, __LINE__)
      end if
      call this%require_replicated_k_workset('reciprocal%fill_moments')

      nk = size(this%eigenvalues, 2)
      if (nk /= this%nk_total) then
         call g_logger%fatal('reciprocal%fill_moments: exact moments require an ' // &
            'undistributed k-mesh (nk_local /= nk_total); k-parallel lands with B4.', &
            __FILE__, __LINE__)
      end if
      nmat = size(this%eigenvectors, 1)

      ! --- First-order k-space velocity operators v_{a,b}(k) = FT[v_{a,b}(R)]. ---
      ! The recursion route's real-space charge velocity blocks, Fourier-summed
      ! with the SAME neighbor map/phase as H(k) = FT[ee(R)].
      call this%hamiltonian%build_realspace_velocity_operators()
      allocate (va_k(nmat, nmat, nk), vb_k(nmat, nmat, nk))
      do ik = 1, nk
         ikg = local_k_index_to_global(this, ik)
         call this%fourier_transform_array(this%hamiltonian%v_a, this%k_points(:, ikg), va_k(:, :, ik))
         call this%fourier_transform_array(this%hamiltonian%v_b, this%k_points(:, ikg), vb_k(:, :, ik))
      end do

      ! --- Exact on-site double-Chebyshev moment per atom type. -----------------
      cond_ll = this%lattice%control%cond_ll
      ntype = this%lattice%ntype
      allocate (mu(nb, nb, cond_ll, cond_ll, ntype))

      block
         integer :: it
         do it = 1, ntype
            site = this%lattice%iz(this%lattice%atlist(it))
            ioff = (site - 1)*nb
            call moment_onsite_block(this%eigenvalues, this%eigenvectors, va_k, vb_k, &
                                     a, b, ioff, nb, cond_ll, mu(:, :, :, :, it))
         end do
      end block

      call g_logger%info('reciprocal%fill_moments: exact k-space Chebyshev moments ' // &
         'filled for ' // trim(int2str(ntype)) // ' atom type(s) on ' // &
         trim(int2str(nk)) // ' k-points.', __FILE__, __LINE__)

      deallocate (va_k, vb_k)
   end subroutine fill_moments

end submodule reciprocal_moments
