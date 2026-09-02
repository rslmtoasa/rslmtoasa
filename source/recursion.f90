 !------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Recursion
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module to handle recursion calculations
!------------------------------------------------------------------------------

module recursion_mod

   use hamiltonian_mod
   use chebyshev_fast_mod, only: cheb_moments_fast, cheb_moments_fast_batched, &
      cheb_moments_fast_mkl_batch, &
      cheb_moments_fast_mkl_sparse, cheb_fast_reset_cache, &
      cheb_moments_stochastic_fast, cheb_moments_orbital_fast
   use haydock_fast_mod, only: block_lanczos_fast
   use lattice_mod
   use control_mod
   use energy_mod
   use sparse_mod, only: sparse
   use precision_mod, only: rp
   use math_mod
   use string_mod
   use logger_mod, only: g_logger
   use rsrec_cuda_plugin_mod, only: rsrec_cuda_backend, get_gpu_context, &
      rsrec_cuda_plugin_compiled, decode_gpu_backend, gpu_backend_bsr, &
      gpu_backend_conv, gpu_backend_fft
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   use timer_mod, only: g_timer
   use basis_mod, only: nb, norb, spin_off, lmax_basis
   implicit none

   private

   !> Module´s main structure
   type, public :: recursion
      !> Hamiltonian
      class(hamiltonian), pointer :: hamiltonian
      !> Sparse algebra helper
      type(sparse), allocatable :: sparse_obj
      !> Lattice
      class(lattice), pointer :: lattice
      !> Energy
      class(energy), pointer :: en
      !> Control
      class(control), pointer :: control
      ! General variables

      !> Scalar recursion coefficients
      real(rp), dimension(:, :, :, :), allocatable :: a, b2
      real(rp), dimension(:), allocatable :: atemp, b2temp
      !> Block recursion coefficients
      complex(rp), dimension(:, :, :, :), allocatable :: a_b, b2_b
      logical :: b2_b_is_sqrt = .false.
      complex(rp), dimension(:, :, :), allocatable :: atemp_b, b2temp_b
      !> Atom list in hopping region
      integer, dimension(:), allocatable :: izero, idum, irlist, nzero, mzero, ndum, mdum
      integer :: irnum !< Atoms within recursion shells
      !> Atom list in hopping region as a function of recursion step
      integer, dimension(:, :), allocatable :: izeroll
      !> Wave functions for recursion hopping (Lanczos)
      complex(rp), dimension(:, :), allocatable :: psi, pmn
      !> Wave functions for block recursion hopping (Lanczos)
      complex(rp), dimension(:, :, :), allocatable :: psi_b, pmn_b, hpsi, hohpsi, enupsi, socpsi
      !> Wave functions for recursion hopping (Chebyshev)
      complex(rp), dimension(:, :, :), allocatable :: psi0, psi1, psi2
      !> Cached full first-order operator blocks for no-hoh fast/GPU CCOR paths
      complex(rp), dimension(:, :, :, :), allocatable :: ee_ccor_work, hall_ccor_work
      !> Chebyshev moments
      complex(rp), dimension(:, :, :, :), allocatable :: mu_n, mu_ng
      complex(rp), dimension(:, :, :, :, :), allocatable :: mu_nm_stochastic
      complex(rp), allocatable :: cheb_mom_temp(:,:)
      !> Variable to save H|psi>
      complex(rp), dimension(:, :), allocatable :: v
      !> Variable to save T(H)
      complex(rp), dimension(:, :, :, :, :), allocatable :: t_h
      type(rsrec_cuda_backend), pointer :: gpu_backend => null()
   contains
      procedure :: hop
      procedure :: crecal
      procedure :: recur
      procedure :: hop_b
      procedure :: hop_b_hoh
      procedure :: crecal_b
      procedure :: recur_b
      procedure :: recur_b_ij
      procedure :: zsqr
      procedure :: bpopt
      procedure :: get_terminf
      procedure :: get_cinf
      procedure :: emami
      procedure :: create_ll_map
      procedure :: cheb_0th_mom
      procedure :: cheb_1st_mom
      procedure :: cheb_1st_mom_hoh
      procedure :: resolve_chebyshev_window
      procedure :: chebyshev_recur_ij
      procedure :: chebyshev_recur_ll
      procedure :: chebyshev_recur_ll_hoh
      procedure :: chebyshev_orbital_mod
      procedure :: chebyshev_recur
      procedure :: compute_moments_stochastic
      procedure :: ham_vec_matmul
      procedure :: ham_hoh_vec_matmul
      procedure :: velo_vec_matmul
      procedure :: velo_hoh_vec_matmul
      procedure :: restore_to_default
      final :: destructor
   end type recursion

   interface recursion
      procedure :: constructor
   end interface recursion


   interface
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

   end subroutine resolve_chebyshev_window

   !> @brief Construct a recursion object from Hamiltonian, energy, and sparse helpers.
   !> @details The recursion object is the real-space spectral kernel used by SCF,
   !>          DOS, exchange, conductivity, and orbital-moment workflows. It stores
   !>          pointers to the upstream objects and allocates the coefficient and
   !>          moment arrays sized by the active lattice/control state.
   !> @param[in] hamiltonian_obj Hamiltonian object providing real-space blocks.
   !> @param[in] energy_obj Energy mesh used by DOS and Green reconstruction.
   !> @param[in] sparse_obj Sparse algebra helper for legacy kernels.
   !> @return Initialized recursion object.
   module function constructor(hamiltonian_obj, energy_obj, sparse_obj) result(obj)
      type(recursion) :: obj
      type(hamiltonian), target, intent(in) :: hamiltonian_obj
      type(energy), target, intent(in) :: energy_obj
      type(sparse), intent(in) :: sparse_obj

   end function constructor

   !> @brief Test whether the CUDA recursion plugin was requested and compiled in.
   !> @details Centralizes the build-time and namelist checks before GPU-specific
   !>          recursion paths attempt to acquire or upload a backend context.
   !> @param[in] this Recursion object whose control flags select plugin usage.
   !> @return True when GPU dispatch is allowed by both build and input state.
   module logical function gpu_plugin_enabled(this)
      class(recursion), intent(in) :: this

   end function gpu_plugin_enabled

   !> @brief Test whether the selected GPU backend expects periodic/k-space data.
   !> @details Separates real-space BSR/convolution recursion backends from GPU
   !>          backends intended for periodic representations so callers can avoid
   !>          uploading incompatible Hamiltonian layouts.
   !> @param[in] this Recursion object whose control flags encode the GPU backend.
   !> @return True for periodic GPU backend selectors.
   module logical function gpu_periodic_backend(this)
      class(recursion), intent(in) :: this
      integer :: backend

   end function gpu_periodic_backend

   !> @brief Decide whether a named recursion feature may run on the GPU plugin.
   !> @details Applies feature-level guards such as scalar-only spin support,
   !>          optional HOH support, backend compatibility, and plugin availability.
   !>          Callers use this before choosing CUDA paths for Lanczos, block
   !>          Lanczos, Chebyshev, stochastic, and intersite workflows.
   !> @param[in] this Recursion object with control/Hamiltonian state.
   !> @param[in] feature Human-readable feature name for logs and diagnostics.
   !> @param[in] require_nsp1 Optional guard for scalar-relativistic-only kernels.
   !> @param[in] allow_hoh Optional guard declaring whether HOH is supported.
   !> @return True when the GPU plugin can be used for the requested feature.
   module logical function gpu_plugin_ready(this, feature, require_nsp1, allow_hoh)
      class(recursion), intent(in) :: this
      character(len=*), intent(in) :: feature
      logical, intent(in), optional :: require_nsp1
      logical, intent(in), optional :: allow_hoh
      integer :: backend
      logical :: nsp1_only, hoh_ok

   end function gpu_plugin_ready

   !> @brief Ensure the CUDA backend has the current real-space Hamiltonian uploaded.
   !> @details Acquires the shared plugin context and transfers the Hamiltonian in
   !>          the layout required by the selected real-space GPU backend. The GPU
   !>          context owns its own fingerprinting, so repeated calls with unchanged
   !>          inputs are cheap.
   !> @param[inout] this Recursion object; sets this%gpu_backend when available.
   !> @note Fatal diagnostics are raised by the backend when a requested upload is unsupported.
   module subroutine gpu_plugin_upload_hamiltonian(this)
      class(recursion), intent(inout) :: this
      integer :: backend

   end subroutine gpu_plugin_upload_hamiltonian

   !> @brief Build cached first-order ccor_2c Hamiltonian blocks for fast kernels.
   !> @details Combines the ordinary and ccor_2c real-space operator blocks into
   !>          work arrays used by no-HOH Chebyshev, block-Lanczos, and transport
   !>          kernels that expect one effective hopping operator.
   !> @param[inout] this Recursion object; populates ee_ccor_work and hall_ccor_work.
   !> @note The cached arrays are reused across moment calls until the recursion object is reset.
   module subroutine ensure_ccor_operator_blocks(this)
      class(recursion), intent(inout) :: this

   end subroutine ensure_ccor_operator_blocks

   !> @brief Run the optimized no-HOH block-Lanczos kernel for one starting block.
   !> @details Wraps haydock_fast for the common block-recursion path where the
   !>          Hamiltonian is represented by a single hopping sweep. The caller has
   !>          initialized psi_b and copies atemp_b/b2temp_b into public coefficient
   !>          arrays after this routine returns.
   !> @param[inout] this Recursion object with current starting block in psi_b.
   !> @param[in] llmax Number of block-Lanczos steps.
   !> @param[in] use_ccor Use cached ccor_2c operator blocks instead of raw blocks.
   module subroutine block_lanczos_fast_nohoh(this, llmax, use_ccor)
      class(recursion), intent(inout) :: this
      integer, intent(in) :: llmax
      logical, intent(in) :: use_ccor

   end subroutine block_lanczos_fast_nohoh

   !> @brief Finalize a recursion object.
   !> @details Releases allocatable coefficient, wavefunction, moment, work, and
   !>          helper arrays owned by the object when it leaves scope.
   !> @param[inout] this Recursion object being finalized.
   module subroutine destructor(this)
      type(recursion) :: this
   end subroutine destructor

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
      complex(rp), dimension(nb, nb) :: dum, dum1, dum2
      complex(rp), dimension(:, :), allocatable :: S_op, L_op
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz
      complex(rp), dimension(:, :, :), allocatable :: psiref, w0, w1, w2, right_vec, v0, v1, v2, cn
      complex(rp), dimension(:, :, :, :), allocatable :: left_vec
      real(rp), dimension(this%en%channels_ldos + 10) :: w, wscale
      real(rp), dimension(this%control%cond_ll) :: kernel
      complex(rp), dimension(nb, nb, this%en%channels_ldos + 10) :: g0
      real(rp) :: a, b, rng, emin_win, emax_win
      complex(rp) :: exp_factor

   end subroutine compute_moments_stochastic

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

   end subroutine hop_b_hoh

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

   end subroutine hop_b

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

   end subroutine recur_b_ij

   !> @brief Generate block-Lanczos coefficients for the selected recursion atoms.
   !> @details Computes the block tridiagonal Haydock coefficients A_n and B_n
   !>          for each atom in lattice%irec. Real-space SCF and block Green’s
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

   end subroutine crecal_b

   !> @brief Replace stored block B^2 coefficients by their matrix square roots.
   !> @details Diagonalizes each Hermitian B2_b block and overwrites it with the
   !>          positive square-root matrix B. Block Green’s-function reconstruction
   !>          calls this before continued-fraction evaluation.
   !> @param[inout] this Recursion object; overwrites local B2_b slices.
   !> @note Uses MPI partition sizes and LAPACK zheev; no collective communication occurs here.
   module subroutine zsqr(this)
      use mpi_mod
      implicit none
      class(recursion), intent(inout) :: this
      !
      integer :: na, LDIM, I, J, K, L, N, M, info, n_glob
      real(rp), dimension(nb) :: ev
      real(rp), dimension(3*nb - 2) ::rwork
      complex(rp), dimension(nb*nb) ::zwork
      complex(rp), dimension(nb, nb) :: B2t, U, DUM, lam

      !
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
   end subroutine get_terminf

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

   end subroutine cheb_1st_mom_hoh

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

   end subroutine chebyshev_recur_ij

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

   end subroutine chebyshev_recur_ll_hoh

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

   end subroutine chebyshev_orbital_mod

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

   end subroutine chebyshev_recur

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

   end subroutine create_ll_map

   !> @brief Apply one scalar Lanczos Hamiltonian hop.
   !> @details Forms H|psi_ll> over the current light cone, accumulates the scalar
   !>          alpha coefficient for each orbital channel, and prepares pmn for
   !>          the scalar Haydock recurrence.
   !> @param[inout] this Recursion object; mutates v, pmn, atemp, and izero.
   !> @param[in] ll Current recursion level whose alpha coefficient is produced.
   !> @note This legacy scalar path is used by control%recur=’lanczos’.
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

   end subroutine recur

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

   end subroutine emami

   !> @brief Reset recursion-owned pointers, buffers, and backend state.
   !> @details Restores the object to a reusable baseline before construction or
   !>          teardown. A full reset also releases large coefficient, moment, and
   !>          wavefunction arrays so a new lattice/control shape can be installed.
   !> @param[inout] this Recursion object to reset.
   !> @param[in] full Optional flag requesting deallocation of all owned arrays.
   !> @note Clears the CUDA backend pointer association but does not own the plugin context.
   module subroutine restore_to_default(this, full)
      class(recursion) :: this
      logical, intent(in), optional :: full
      integer :: lmax

   end subroutine restore_to_default

   end interface

end module recursion_mod
