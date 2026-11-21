!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Reciprocal
!
!> @author
!> Generated following existing codebase patterns
!
! DESCRIPTION:
!> Module to handle reciprocal space calculations including k-point mesh
!> generation and Fourier transform of real-space Hamiltonians for
!> k-space band structure and DOS calculations
!>
!> This module provides functionality to:
!> - Generate Monkhorst-Pack k-point meshes
!> - Transform real-space Hamiltonians to k-space via Fourier transforms
!> - Handle spin-orbit coupling in k-space
!> - Support different basis sets (sp, spd, spdf)
!> - Build complete k-space Hamiltonians for band structure calculations
!>
!> Usage example:
!> ```fortran
!> type(reciprocal) :: recip
!> recip = reciprocal(hamiltonian_obj)
!> call recip%generate_mp_mesh()
!> call recip%build_kspace_hamiltonian()
!> if (include_so) call recip%build_kspace_hamiltonian_so()
!> call recip%build_total_hamiltonian()
!> ```
!------------------------------------------------------------------------------

module reciprocal_mod

   use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_double_complex, c_int, c_long_long, c_null_ptr, c_ptr
   use control_mod
   use lattice_mod
   use hamiltonian_mod
   use charge_mod
   use energy_mod, only: energy
   use green_mod, only: green
   use spectrum_bounds_mod, only: bounds, compute_spectrum_bounds
   use precision_mod, only: rp
   use math_mod
   use string_mod, only: int2str, real2str, fmt, lower
   use logger_mod, only: g_logger
   use timer_mod, only: g_timer
   use symmetry_mod, only: symmetry
   use basis_mod, only: lmax_basis, nb, norb, spin_off
   use spin_density_mod, only: spin_density
   use mpi_mod, only: rank, numprocs, ierr, get_mpi_range, g_parallel_context, get_cuda_device_override, &
                      chain_recv, chain_send_bcast
   use kpoint_workset_mod, only: kpoint_workset, make_kpoint_workset
#ifdef USE_MPI
   use mpi
#endif
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Reusable CPU storage for reciprocal-space tiles.  H, S, and all matrix
   !> temporaries use (nmat,nmat,nbatch): the first two dimensions are one
   !> Fortran column-major matrix and the third dimension is a contiguous,
   !> strided batch.  A workspace is deliberately owned by one reciprocal
   !> object/worker; callers must not share it across OpenMP workers.
   type, public :: reciprocal_workspace
      integer :: nmat = 0
      integer :: tile_capacity = 0
      integer :: active_tile_length = 0
      integer :: lwork = 0
      integer :: cached_operator_generation = -1
      !> Diagnostic counters: preparation events only.  A stable prepared
      !> tile must leave both counts unchanged inside its k-point loop.
      integer :: storage_allocations = 0
      integer :: lapack_workspace_queries = 0
      integer :: capacity_reuses = 0
      logical :: generalized = .false.
      complex(rp), allocatable :: h(:, :, :), s(:, :, :)
      complex(rp), allocatable :: eeo(:, :, :), hoh(:, :, :), hcc(:, :, :)
      complex(rp), allocatable :: phase(:, :, :)
      !> Private validation copy used by ZPOTRF before a generalized solve.
      !> It preserves S(k), which ZHEGV later consumes in place.
      complex(rp), allocatable :: overlap_cholesky(:, :)
      real(rp), allocatable :: points(:, :)
      real(rp), allocatable :: eigenvalue(:, :)
      complex(rp), allocatable :: eigenvector(:, :, :), lapack_work(:)
      real(rp), allocatable :: lapack_rwork(:)
      integer, allocatable :: info(:)
   contains
      procedure :: ensure_capacity => reciprocal_workspace_ensure_capacity
      procedure :: clear => reciprocal_workspace_clear
      procedure :: restore_to_default => reciprocal_workspace_restore_to_default
      final :: reciprocal_workspace_destructor
   end type reciprocal_workspace

   !> Non-owning view of the long-lived operator and geometry state required
   !> by Fourier assembly.  It owns no Hamiltonian matrices; phase and HOH
   !> scratch belong exclusively to reciprocal_workspace.
   type, public :: reciprocal_assembler
      class(hamiltonian), pointer :: hamiltonian => null()
      class(lattice), pointer :: lattice => null()
      class(control), pointer :: control => null()
      real(rp), allocatable :: ham_vec_type_direct(:, :, :)
      logical :: use_second_order = .false.
      logical :: use_orthogonal = .false.   !> kspace_ham_order = 'orthogonal': exact-ASA H^gamma (gbt_single_q only)
   contains
      procedure :: assemble_batch => reciprocal_assembler_assemble_batch
      procedure :: assemble_overlap_batch => reciprocal_assembler_assemble_overlap_batch
      procedure :: assemble_one => reciprocal_assembler_assemble_one
      procedure :: assemble_gbt => reciprocal_assembler_assemble_gbt
   end type reciprocal_assembler

   !> Residency values reported by reciprocal_execution_capabilities.  The
   !> first implementation is deliberately host-only; callers must not infer
   !> a device transfer from an implementation name.
   integer, parameter, public :: reciprocal_residency_host = 1
   integer, parameter, public :: reciprocal_residency_device = 2

   !> Backend feature and tile limits.  This is a typed contract rather than
   !> a collection of backend-name string tests in physics callers.
   type, public :: reciprocal_execution_capabilities
      logical :: standard_hermitian = .false.
      logical :: generalized_hermitian = .false.
      logical :: eigenvalues_only = .false.
      logical :: eigenvectors = .false.
      logical :: first_order_assembly = .false.
      logical :: second_order_assembly = .false.
      logical :: overlap = .false.
      integer :: residency = reciprocal_residency_host
   end type reciprocal_execution_capabilities

   !> One execution request owns its point list and any compatibility input
   !> matrices.  A resident backend is therefore free to retain a request until
   !> synchronize; the CPU adapter copies legacy mesh caches at this boundary.
   type, public :: reciprocal_execution_request
      logical :: assemble_hamiltonian = .true.
      logical :: assemble_overlap = .false.
      logical :: solve_eigensystem = .true.
      logical :: generalized = .false.
      logical :: request_eigenvectors = .true.
      logical :: request_assembled_hamiltonian = .false.
      logical :: request_assembled_overlap = .false.
      integer :: operator_generation = -1
      real(rp), allocatable :: k_points(:, :)
      complex(rp), allocatable :: input_hamiltonian(:, :, :)
      complex(rp), allocatable :: input_overlap(:, :, :)
   end type reciprocal_execution_request

   !> Results are caller-owned allocations.  Matrix allocation signals an
   !> explicit diagnostic request; eigenvectors are absent for values-only
   !> requests.  All slots describe exactly local_point_count points.
   type, public :: reciprocal_execution_result
      integer :: local_point_count = 0
      integer :: operator_generation = -1
      logical :: assembled_hamiltonian_valid = .false.
      logical :: assembled_overlap_valid = .false.
      logical :: eigenvalues_valid = .false.
      logical :: eigenvectors_valid = .false.
      integer :: resident_token = 0
      complex(rp), allocatable :: hamiltonian(:, :, :)
      complex(rp), allocatable :: overlap(:, :, :)
      real(rp), allocatable :: eigenvalues(:, :)
      complex(rp), allocatable :: eigenvectors(:, :, :)
   end type reciprocal_execution_result

   !> One host-side Lehmann request.  The execution backend consumes the
   !> established eigenpair representation and returns the canonical Green
   !> blocks; it does not own a second Green-function object or expose device
   !> pointers to physics code.
   type, public :: reciprocal_lehmann_request
      real(rp), allocatable :: eigenvalues(:, :)
      complex(rp), allocatable :: eigenvectors(:, :, :)
      real(rp), allocatable :: k_points(:, :)
      complex(rp), allocatable :: z_contour(:)
      real(rp), allocatable :: dr(:, :)
      integer, allocatable :: ioffset(:), joffset(:)
      integer :: nblk = 0
      ! ACC-11 resident requests may omit host eigenpairs; these dimensions
      ! identify the backend-owned device tile in that case.
      integer :: nmat = 0
      integer :: nk = 0
   end type reciprocal_lehmann_request

   type, public :: reciprocal_lehmann_result
      complex(rp), allocatable :: blocks(:, :, :, :)
      real(rp) :: h2d_seconds = 0.0_rp
      real(rp) :: contraction_seconds = 0.0_rp
      real(rp) :: d2h_seconds = 0.0_rp
      logical :: valid = .false.
   end type reciprocal_lehmann_result

   type, abstract, public :: reciprocal_execution_backend
   contains
      procedure(backend_initialize_if), deferred :: initialize
      procedure(backend_capabilities_if), deferred :: capabilities
      procedure(backend_prepare_if), deferred :: prepare_operator
      procedure(backend_execute_if), deferred :: execute_batch
      procedure(backend_execution_metrics_if), deferred :: execution_metrics
      procedure(backend_synchronize_if), deferred :: synchronize
      procedure(backend_release_if), deferred :: release
   end type reciprocal_execution_backend

   !> CPU implementation.  It composes the RF-04 assembler and exclusive
   !> workspace so an eventual device implementation can keep both assembly
   !> and eigensolution resident behind the same request/result boundary.
   type, extends(reciprocal_execution_backend), public :: lapack_reciprocal_backend
      type(reciprocal_assembler) :: assembler
      type(reciprocal_workspace) :: workspace
      !> Request-shape counters are backend-contract observability for the
      !> reciprocal tile tests.  They count execute_batch calls, not k points.
      integer :: execute_batch_requests = 0
      integer :: combined_assembly_solve_requests = 0
      integer :: assemble_only_requests = 0
      integer :: input_hamiltonian_solve_requests = 0
      real(rp) :: host_assembly_seconds = 0.0_rp
      logical :: initialized = .false.
   contains
      procedure :: initialize => lapack_backend_initialize
      procedure :: capabilities => lapack_backend_capabilities
      procedure :: prepare_operator => lapack_backend_prepare_operator
      procedure :: execute_batch => lapack_backend_execute_batch
      procedure :: execution_metrics => lapack_backend_execution_metrics
      procedure :: synchronize => lapack_backend_synchronize
      procedure :: release => lapack_backend_release
      final :: lapack_backend_destructor
   end type lapack_reciprocal_backend

   !> CUDA standard-Hermitian adapter.  The opaque context owns the CUDA
   !> stream, cuSOLVER handle, and reusable tile workspace.
   type, extends(reciprocal_execution_backend), public :: cuda_reciprocal_backend
      type(c_ptr) :: context = c_null_ptr
      integer :: device = -1
      integer :: prepared_operator_generation = -1
      integer :: execute_batch_requests = 0
      integer :: input_hamiltonian_solve_requests = 0
      integer :: operator_prepare_requests = 0
      integer :: operator_prepare_reuses = 0
      real(rp) :: h2d_seconds = 0.0_rp
      real(rp) :: gpu_solve_seconds = 0.0_rp
      real(rp) :: d2h_seconds = 0.0_rp
      real(rp) :: d2h_values_seconds = 0.0_rp
      real(rp) :: d2h_vectors_seconds = 0.0_rp
      real(rp) :: sync_seconds = 0.0_rp
      real(rp) :: host_staging_seconds = 0.0_rp
      real(rp) :: host_conversion_seconds = 0.0_rp
      real(rp) :: host_widen_seconds = 0.0_rp
      real(rp) :: total_reciprocal_seconds = 0.0_rp
      integer(c_long_long) :: h2d_bytes = 0_c_long_long
      integer(c_long_long) :: d2h_values_bytes = 0_c_long_long
      integer(c_long_long) :: d2h_vectors_bytes = 0_c_long_long
      integer(c_long_long) :: cuda_malloc_count = 0_c_long_long
      integer(c_long_long) :: cuda_free_count = 0_c_long_long
      integer(c_long_long) :: workspace_query_count = 0_c_long_long
      integer(c_long_long) :: workspace_reuse_count = 0_c_long_long
      integer(c_long_long) :: event_create_count = 0_c_long_long
      integer(c_long_long) :: event_destroy_count = 0_c_long_long
      integer(c_long_long) :: pinned_alloc_count = 0_c_long_long
      integer(c_long_long) :: pinned_free_count = 0_c_long_long
      integer :: pinned_host_active = 0
      integer :: timing_calls = 0
      integer :: resident_token = 0
      logical :: initialized = .false.
   contains
      procedure :: initialize => cuda_backend_initialize
      procedure :: capabilities => cuda_backend_capabilities
      procedure :: prepare_operator => cuda_backend_prepare_operator
      procedure :: execute_batch => cuda_backend_execute_batch
      procedure :: contract_lehmann => cuda_backend_contract_lehmann
      procedure :: contract_lehmann_resident => cuda_backend_contract_lehmann_resident
      procedure :: resident_eigensystem_matches => cuda_backend_resident_eigensystem_matches
      procedure :: execution_metrics => cuda_backend_execution_metrics
      procedure :: synchronize => cuda_backend_synchronize
      procedure :: release => cuda_backend_release
      final :: cuda_backend_destructor
   end type cuda_reciprocal_backend

   !> Module’s main procedure
   type, public :: reciprocal
      !> Hamiltonian
      class(hamiltonian), pointer :: hamiltonian
      !> Lattice
      class(lattice), pointer :: lattice
      !> Control
      class(control), pointer :: control

      ! K-point mesh variables
      !> Number of k-points in each direction (nx, ny, nz)
      integer, dimension(3) :: nk_mesh
      !> Dense mesh dimensions for post_processing=’fermi_surface’.
      !> Non-positive values fall back to the regular nk_mesh dimensions.
      integer :: fs_nk1, fs_nk2, fs_nk3
      !> Requested symmetry reduction for FS post-processing. The export
      !> itself is always promoted to the complete BZ mesh.
      logical :: fs_use_symmetry_reduction
      !> Authoritative k-point ownership.  points/weights in this object are
      !> local to the current rank; replicated worksets own the entire mesh.
      type(kpoint_workset) :: k_workset
      !> Total number of k-points.  Transitional compatibility view of
      !> k_workset%nk_global; do not write outside mesh construction.
      integer :: nk_total
      !> Transitional read-only full-list compatibility view.  New consumers
      !> must use k_workset%points and k_workset%weights.
      real(rp), dimension(:, :), allocatable :: k_points
      !> K-point weights for Brillouin zone integration
      real(rp), dimension(:), allocatable :: k_weights
      !> Local k-point ownership for distributed mesh workflows
      integer :: nk_local
      integer :: k_start
      integer :: k_end
      integer, dimension(:), allocatable :: k_l2g_map
      integer, dimension(:), allocatable :: k_g2l_map
      logical :: k_mesh_distributed_active
      !> Distribute a tetrahedron/Blochl mesh over MPI ranks (set by the k-space SCF cache only)
      logical :: distribute_k_mesh
      !> Include time-reversal symmetry in k-point generation
      logical :: use_time_reversal
      !> Offset for k-point mesh (for shifted grids)
      real(rp), dimension(3) :: k_offset

      ! Reciprocal lattice vectors
      !> Reciprocal lattice vectors (3x3 matrix)
      real(rp), dimension(3, 3) :: reciprocal_vectors
      !> Volume of reciprocal unit cell
      real(rp) :: reciprocal_volume

      ! K-space Hamiltonian
      !> Bulk Hamiltonian in k-space
      complex(rp), dimension(:, :, :), allocatable :: hk_bulk
      !> Spin-orbit Hamiltonian in k-space
      complex(rp), dimension(:, :, :), allocatable :: hk_so
      !> Total k-space Hamiltonian (bulk + SO)
      complex(rp), dimension(:, :, :), allocatable :: hk_total
      !> Overlap matrix in k-space
      complex(rp), dimension(:, :, :), allocatable :: sk_overlap
      !> Hamiltonian operator_generation used to build the current H(k) and
      !> every spectrum/DOS/density object derived from it.
      integer :: cached_operator_generation
      ! ACC-11 token for the last complete normal-mesh CUDA eigensystem.  The
      ! token is opaque to physics code and is invalidated by later backend
      ! solves or backend recreation.
      integer :: device_eigensystem_token = 0
      !> Reciprocal solver mode: ’ham_only’ or ’generalized_overlap_proxy’ (S = I + eeo; not a
      !> formal Kanpur overlap). ’generalized_overlap_kanpur’ is not implemented and stops the run.
      character(len=32) :: reciprocal_mode
      !> Production normal-mesh backend: ’lapack’ (default) or ’cuda’.
      character(len=16) :: reciprocal_backend = 'lapack'
      !> K-space Hamiltonian order:
      !>   ’first’  -> H(k) = h(k)            (first-order, current behaviour)
      !>   ’second’ -> H(k) = E_nu + h(k) - [hoh](k) + L.S   (second-order ASA)
      !> ’orthogonal’ -> exact-ASA H^gamma (gbt_kspace_horth).
      character(len=16) :: kspace_ham_order
      !> Kanpur-alignment diagnostics toggle
      logical :: kanpur_diagnostics
      !> Optional H(Gamma) bounds diagnostics
      logical :: gamma_bounds_diagnostics
      !> Experimental finite real-space HALL diagonalization
      logical :: hall_diag_experimental
      ! Conservative internal CPU tile size.  It is intentionally not a
      ! namelist option; tests may set it programmatically.
      integer :: reciprocal_tile_size
      type(reciprocal_workspace) :: workspace
      !> Persistent default execution context.  Legacy workspace remains a
      !> compatibility diagnostic view while this backend owns active tiles.
      class(reciprocal_execution_backend), allocatable :: execution_backend
      !> ACC-11: a Lehmann consumer requested a complete CUDA-resident normal
      !> mesh.  This is set by fill_green_lehmann before the first build; all
      !> other CUDA workflows retain their configured tile size.
      logical :: resident_lehmann_handoff_requested = .false.

      ! Band structure variables
      !> Maximum number of orbital channels per atom type
      integer :: max_orbs
      !> Basis set type indicator (sp=4, spd=9, spdf=16)
      integer, dimension(:), allocatable :: basis_size
      !> Include spin-orbit coupling
      logical :: include_so

      ! High-symmetry k-path variables
      !> K-path for band structure (high-symmetry points)
      real(rp), dimension(:, :), allocatable :: k_path
      !> Number of k-points along the path
      integer :: nk_path
      !> K-path labels (high-symmetry point names)
      character(len=10), dimension(:), allocatable :: k_labels
      !> K-path distances for plotting
      real(rp), dimension(:), allocatable :: k_distances
      
   ! Internal logging control
   !> When .true., suppress internal debug/info prints (structure factors, per-k debug)
   logical :: suppress_internal_logs

      ! Eigenvalue/eigenvector storage
      !> Eigenvalues for k-mesh [nbands, nk_total]
      real(rp), dimension(:, :), allocatable :: eigenvalues
      !> Eigenvalues for k-path [nbands, nk_path]
      real(rp), dimension(:, :), allocatable :: eigenvalues_path
      !> Eigenvectors for k-mesh [max_orb_channels, nbands, nk_total]
      complex(rp), dimension(:, :, :), allocatable :: eigenvectors
      !> Eigenvectors for k-path [max_orb_channels, nbands, nk_path]
      complex(rp), dimension(:, :, :), allocatable :: eigenvectors_path

      ! Symmetry analysis object
      type(symmetry) :: symmetry_analysis

      ! Magnetic symmetry group behind the reduced k mesh (see setup_magnetic_symmetry)
      !> Number of operations
      integer :: n_sym_ops = 0
      !> R(i,j) on fractional coordinates, f' = R f + t [3,3,n_sym_ops]
      integer, dimension(:, :, :), allocatable :: sym_rotations
      real(rp), dimension(:, :), allocatable :: sym_translations
      !> -1 for time-reversal-combined operations (k -> -R k), +1 otherwise
      integer, dimension(:), allocatable :: sym_spin_flips
      !> Site i is mapped onto site sym_atom_perm(i, op) by operation op [nrec, n_sym_ops]
      integer, dimension(:, :), allocatable :: sym_atom_perm

      ! Density of states variables
      !> Energy grid for DOS calculation [n_energy_points]
      real(rp), dimension(:), allocatable :: dos_energy_grid
      !> Number of energy points for DOS
      integer :: n_energy_points
      !> Energy range for DOS [min_energy, max_energy]
      real(rp), dimension(2) :: dos_energy_range
      !> Total DOS [n_energy_points]
      real(rp), dimension(:), allocatable :: total_dos
      !> Integrated total DOS / number of states [n_energy_points]
      real(rp), dimension(:), allocatable :: total_nos
      !> Projected DOS [n_sites, n_orb_types, n_spin, n_energy_points]
      real(rp), dimension(:, :, :, :), allocatable :: projected_dos
      !> Band moments [m0, m1, m2] for each projection
      real(rp), dimension(:, :, :, :), allocatable :: band_moments
      !> WP7 shared rotating-frame density contract, filled by the k-space
      !> producer `accumulate_spin_density_kspace`. `band_moments` is projected
      !> out of THIS object, after accumulation, against its own explicit axis.
      type(spin_density) :: rf_density
      !> Directional DOS x-component [n_energy_points]
      real(rp), dimension(:), allocatable :: dos_mx_tot
      !> Directional DOS y-component [n_energy_points]
      real(rp), dimension(:), allocatable :: dos_my_tot
      !> Directional DOS z-component [n_energy_points]
      real(rp), dimension(:), allocatable :: dos_mz_tot
      !> Projected directional DOS [n_sites, n_orb_types, n_spin, 3(x/y/z), n_energy_points]
      real(rp), dimension(:, :, :, :, :), allocatable :: projected_dos_moments
      !> DOS calculation method (’tetrahedron’ or ’gaussian’)
      character(len=20) :: dos_method
      !> Gaussian smearing parameter (in energy units)
      real(rp) :: gaussian_sigma
      !> Temperature for Fermi-Dirac distribution (in Kelvin)
      real(rp) :: temperature
      !> Fermi level for band moments integration
      real(rp) :: fermi_level
      !> Total number of valence electrons for Fermi level finding
      real(rp) :: total_electrons
      !> Projection-free electron count from the current eigenvalue occupations.
      real(rp) :: canonical_electron_count
      !> Projection-free occupied-eigenvalue band energy for the current spectrum.
      real(rp) :: canonical_band_energy
      !> Per-site spin axes used by the DOS projections, stored by prepare_density_of_states.
      real(rp), allocatable :: dos_site_axes(:, :)
      !> Global raw k-weight sum used by the canonical normalization contract.
      real(rp) :: canonical_weight_sum
      !> True only after the current eigensystem has been evaluated canonically.
      logical :: canonical_energy_valid
      !> Flag to automatically find Fermi level from the current eigensystem
      logical :: auto_find_fermi
      !> Number of sites for projections
      integer :: n_sites
      !> Number of orbital types (s, p, d, f)
      integer :: n_orb_types
      !> Number of spin components (RS-LMTO uses spin-polarized two-component blocks)
      integer :: n_spin_components

      ! Tetrahedron method variables
      !> Tetrahedron corners for each k-point [4, nk_total]
      integer, dimension(:, :), allocatable :: tetrahedra
      !> Tetrahedron volumes [nk_total]
      real(rp), dimension(:), allocatable :: tetrahedron_volumes
      !> Number of tetrahedra
      integer :: n_tetrahedra
      !> Irreducible tetrahedra and multiplicities of the dos_method='blochl' occupations, built once per mesh
      integer, dimension(:, :), allocatable :: occ_tetra
      integer, dimension(:), allocatable :: occ_tetra_mult
      !> nk_mesh followed by full_to_irred_k at the time occ_tetra was built; a different mesh or symmetry map rebuilds it
      integer, dimension(:), allocatable :: occ_tetra_key

      ! K-path (band structure) control variables
      !> Automatic k-path generation enabled
      logical :: auto_kpath
      !> Number of k-points per segment in band structure
      integer :: nk_per_segment
      !> Use symmetry reduction for k-mesh
      logical :: use_symmetry_reduction
      !> Enforce strict symmetry consistency checks
      logical :: strict_symmetry_checks
      !> Dump symmetry k-point mapping diagnostics to file
      logical :: dump_symmetry_kmap
      !> Full-mesh to irreducible-k mapping (size = nk1*nk2*nk3)
      integer, dimension(:), allocatable :: full_to_irred_k
      !> Irreducible-k representative indices in full mesh
      integer, dimension(:), allocatable :: irred_to_full_k

      ! Real-space neighbor vectors per atom type (for multi-site H_k)
      !> Neighbor vectors for each atom type [3, nn_max, ntype]
      !> These are the R vectors for Fourier transform H(k) = Σ_R H(R) e^(ik·R)
      !> Indexed properly for multi-site systems: ham_vec_type(coord, neighbor, atom_type)
      real(rp), dimension(:, :, :), allocatable :: ham_vec_type
      !> Neighbor vectors in fractional coordinates [3, nn_max, ntype]
      real(rp), dimension(:, :, :), allocatable :: ham_vec_type_direct

      ! k-space Green’s-function engine (milestone B2, reciprocal_green)
      !> Retarded broadening for the k-space Green contour z = E + i*green_eta (Ry)
      real(rp) :: green_eta
      !> Backend selector for fill_green: ’lehmann’ (E, Sigma=0) or ’dyson’ (D)
      character(len=16) :: green_backend

      ! Fermi-surface post-processing eigenpair export. The configured name is
      ! a base name; ’.bin’ and ’.meta’ are appended.
      logical :: write_eigenpair_projections
      !> Write independently diagonalized global collinear spin sectors to a
      !> separate ’<base>.spin.bin’ payload when H(k) and S(k) are spin-block
      !> diagonal over the complete exported mesh.
      logical :: write_spin_resolved_eigenpairs
      character(len=256) :: eigenpair_output_file

   contains
      procedure :: generate_mp_mesh
      procedure :: generate_reciprocal_vectors
      procedure :: build_kspace_hamiltonian
      procedure :: execute_normal_mesh_tiles
      procedure :: build_neighbor_vectors
      procedure :: calculate_structure_factors
      procedure :: fourier_transform_hamiltonian
      procedure :: fourier_transform_hamiltonian_second_order
      procedure :: fourier_transform_array
      procedure :: fourier_transform_overlap
      procedure :: fold_kpoint
      procedure :: build_hamiltonian_at_kpoint
      procedure :: make_reciprocal_assembler
      procedure :: make_execution_backend
      procedure :: calculate_eigenpairs_at_kpoints
      procedure :: require_replicated_k_workset
      procedure :: gather_k_columns
      procedure :: set_basis_sizes
      procedure :: check_multisite_hamiltonian_diagonal
      procedure :: build_kspace_overlap
      procedure :: diagonalize_hamiltonian
      procedure :: has_nonzero_q_gbt
      procedure :: validate_nonzero_q_gbt
      procedure :: force_full_bz_for_nonzero_q_gbt
      procedure :: print_kanpur_mapping
      procedure :: check_overlap_properties
      procedure :: run_gamma_bounds_diagnostics
      procedure :: diagonalize_hall_experimental
      procedure :: calculate_band_structure
      procedure :: calculate_density_of_states
      procedure :: prepare_density_of_states
      procedure :: calculate_dos_output
      procedure :: calculate_dos_tetrahedron
      procedure :: calculate_dos_gaussian
      procedure :: calculate_dos_blochl
      procedure :: setup_dos_energy_grid
      procedure :: setup_tetrahedra
      procedure :: get_kpoint_index
      procedure :: project_dos_orbitals
      procedure :: project_dos_orbitals_gaussian
      procedure :: project_dos_orbitals_tetrahedron
      procedure :: calculate_band_moments
      procedure :: accumulate_spin_density_kspace
      procedure :: fill_band_moments_from_spin_density
      procedure :: evaluate_eigenvalue_occupations
      procedure :: find_fermi_level_from_eigenvalues
      procedure :: calculate_canonical_band_energy
      procedure :: calculate_band_energy_from_total_dos
      procedure :: calculate_band_energy_from_moments
      procedure :: calculate_adaptive_sigma
      procedure :: calculate_ldm_from_projected_dos
      procedure :: calculate_spinor_dm
      procedure :: integrate_dos_up_to_energy
      procedure :: write_dos_to_file
      procedure :: restore_to_default
      procedure :: invalidate_spectral_cache
      procedure :: invalidate_if_operator_changed
      procedure :: build_from_file
      procedure :: generate_reduced_kpoint_mesh
      procedure :: setup_magnetic_symmetry
      procedure :: symmetrize_spin_density
      procedure :: symmetrize_site_moments
      procedure :: symmetrize_projected_dos
      procedure :: symmetrize_site_dm
      procedure :: validate_symmetry_kmap
      procedure :: write_symmetry_kmap_dump
      procedure :: full_mesh_operations
      procedure :: full_mesh_projection
      procedure :: build_irreducible_tetrahedra
      procedure :: blochl_occupations
      ! k-space Green’s-function engine (milestone B2)
      procedure :: fill_green
      procedure :: fill_moments
      procedure :: build_green_contour
      procedure :: calculate_bsf
      procedure :: write_kspace_eigenpairs
      final     :: destructor
   end type reciprocal

   interface reciprocal
      procedure :: constructor
   end interface

   public :: fermi_dirac_occupation
   public :: fermi_dirac_entropy_term, fermi_dirac_entropy
   public :: tetra_corner_weights
   public :: tetra_integrated_weights
   real(rp), parameter, public :: kB_Ry_per_K = 6.3336814e-6_rp

   interface
      !> Numerically stable Fermi-Dirac occupation at the effective kT used for
      !> EF, electron count, and EBAND.
      pure module function fermi_dirac_occupation(eigenvalue, fermi_level, kT) result(occupation)
         real(rp), intent(in) :: eigenvalue, fermi_level, kT
         real(rp) :: occupation
      end function fermi_dirac_occupation
      !> Entropy of one Fermi-Dirac state; exact empty/full states contribute zero.
      pure module function fermi_dirac_entropy_term(occupation) result(entropy)
         real(rp), intent(in) :: occupation
         real(rp) :: entropy
      end function fermi_dirac_entropy_term
      !> K-weighted, spin-inclusive entropy from a band eigenvalue table.
      pure module function fermi_dirac_entropy(eigenvalues, weights, fermi_level, kT) result(entropy)
         real(rp), intent(in) :: eigenvalues(:, :), weights(:), fermi_level, kT
         real(rp) :: entropy
      end function fermi_dirac_entropy
   end interface

#ifdef USE_CUDA_RECIPROCAL
   interface
      function rslmto_reciprocal_cuda_device_count(count) bind(C, name='rslmto_reciprocal_cuda_device_count')
         import c_int
         integer(c_int) :: count
         integer(c_int) :: rslmto_reciprocal_cuda_device_count
      end function rslmto_reciprocal_cuda_device_count

      function rslmto_reciprocal_cuda_create(device) bind(C, name='rslmto_reciprocal_cuda_create')
         import c_int, c_ptr
         integer(c_int), value :: device
         type(c_ptr) :: rslmto_reciprocal_cuda_create
      end function rslmto_reciprocal_cuda_create

      function rslmto_reciprocal_cuda_prepare_operator(context, operator_generation) &
         bind(C, name='rslmto_reciprocal_cuda_prepare_operator')
         import c_int, c_ptr
         type(c_ptr), value :: context
         integer(c_int), value :: operator_generation
         integer(c_int) :: rslmto_reciprocal_cuda_prepare_operator
      end function rslmto_reciprocal_cuda_prepare_operator

      function rslmto_reciprocal_cuda_solve_zheevd_batch(context, n, batch_size, hamiltonians, eigenvalues, eigenvectors, &
                                                         request_eigenvectors) &
         bind(C, name='rslmto_reciprocal_cuda_solve_zheevd_batch')
         import c_double, c_double_complex, c_int, c_ptr
         type(c_ptr), value :: context
         integer(c_int), value :: n, batch_size
         complex(c_double_complex), intent(in) :: hamiltonians(*)
         real(c_double), intent(out) :: eigenvalues(*)
         complex(c_double_complex), intent(out) :: eigenvectors(*)
         integer(c_int), value :: request_eigenvectors
         integer(c_int) :: rslmto_reciprocal_cuda_solve_zheevd_batch
      end function rslmto_reciprocal_cuda_solve_zheevd_batch

      function rslmto_reciprocal_cuda_get_resident_token(context) &
         bind(C, name='rslmto_reciprocal_cuda_get_resident_token')
         import c_int, c_ptr
         type(c_ptr), value :: context
         integer(c_int) :: rslmto_reciprocal_cuda_get_resident_token
      end function rslmto_reciprocal_cuda_get_resident_token

      function rslmto_reciprocal_cuda_resident_eigensystem_matches(context, n, batch_size, token) &
         bind(C, name='rslmto_reciprocal_cuda_resident_eigensystem_matches')
         import c_int, c_ptr
         type(c_ptr), value :: context
         integer(c_int), value :: n, batch_size, token
         integer(c_int) :: rslmto_reciprocal_cuda_resident_eigensystem_matches
      end function rslmto_reciprocal_cuda_resident_eigensystem_matches

      function rslmto_reciprocal_cuda_contract_lehmann(context, nmat, nk, ne, npair, nblk, eigenvalues, eigenvectors, &
                                                        k_points, z_contour, dr, ioffset, joffset, blocks, h2d_seconds, &
                                                        contraction_seconds, d2h_seconds) &
         bind(C, name='rslmto_reciprocal_cuda_contract_lehmann')
         import c_double, c_double_complex, c_int, c_ptr
         type(c_ptr), value :: context
         integer(c_int), value :: nmat, nk, ne, npair, nblk
         real(c_double), intent(in) :: eigenvalues(*), k_points(*), dr(*)
         complex(c_double_complex), intent(in) :: eigenvectors(*), z_contour(*)
         integer(c_int), intent(in) :: ioffset(*), joffset(*)
         complex(c_double_complex), intent(out) :: blocks(*)
         real(c_double), intent(out) :: h2d_seconds, contraction_seconds, d2h_seconds
         integer(c_int) :: rslmto_reciprocal_cuda_contract_lehmann
      end function rslmto_reciprocal_cuda_contract_lehmann

      function rslmto_reciprocal_cuda_contract_lehmann_resident(context, nmat, nk, ne, npair, nblk, resident_token, k_points, z_contour, &
                                                                 dr, ioffset, joffset, blocks, h2d_seconds, contraction_seconds, &
                                                                 d2h_seconds) &
         bind(C, name='rslmto_reciprocal_cuda_contract_lehmann_resident')
         import c_double, c_double_complex, c_int, c_ptr
         type(c_ptr), value :: context
         integer(c_int), value :: nmat, nk, ne, npair, nblk
         integer(c_int), value :: resident_token
         real(c_double), intent(in) :: k_points(*), dr(*)
         complex(c_double_complex), intent(in) :: z_contour(*)
         integer(c_int), intent(in) :: ioffset(*), joffset(*)
         complex(c_double_complex), intent(out) :: blocks(*)
         real(c_double), intent(out) :: h2d_seconds, contraction_seconds, d2h_seconds
         integer(c_int) :: rslmto_reciprocal_cuda_contract_lehmann_resident
      end function rslmto_reciprocal_cuda_contract_lehmann_resident

      subroutine rslmto_reciprocal_cuda_get_timings(context, h2d_seconds, solve_seconds, d2h_seconds, calls) &
         bind(C, name='rslmto_reciprocal_cuda_get_timings')
         import c_double, c_int, c_ptr
         type(c_ptr), value :: context
         real(c_double), intent(out) :: h2d_seconds, solve_seconds, d2h_seconds
         integer(c_int), intent(out) :: calls
      end subroutine rslmto_reciprocal_cuda_get_timings

      subroutine rslmto_reciprocal_cuda_get_detailed_timings(context, host_conversion_seconds, host_staging_seconds, h2d_seconds, &
                                                              solve_seconds, d2h_seconds, d2h_values_seconds, &
                                                              d2h_vectors_seconds, sync_seconds, host_widen_seconds, &
                                                              total_seconds, h2d_bytes, d2h_values_bytes, d2h_vectors_bytes, &
                                                              pinned_host_active, calls) &
         bind(C, name='rslmto_reciprocal_cuda_get_detailed_timings')
         import c_double, c_int, c_long_long, c_ptr
         type(c_ptr), value :: context
         real(c_double), intent(out) :: host_conversion_seconds, host_staging_seconds, h2d_seconds, solve_seconds
         real(c_double), intent(out) :: d2h_seconds, d2h_values_seconds, d2h_vectors_seconds, sync_seconds
         real(c_double), intent(out) :: host_widen_seconds, total_seconds
         integer(c_long_long), intent(out) :: h2d_bytes, d2h_values_bytes, d2h_vectors_bytes
         integer(c_int), intent(out) :: pinned_host_active
         integer(c_int), intent(out) :: calls
      end subroutine rslmto_reciprocal_cuda_get_detailed_timings

      subroutine rslmto_reciprocal_cuda_get_resource_counters(context, cuda_malloc_count, cuda_free_count, &
                                                               workspace_query_count, workspace_reuse_count, &
                                                               event_create_count, event_destroy_count, &
                                                               pinned_alloc_count, pinned_free_count) &
         bind(C, name='rslmto_reciprocal_cuda_get_resource_counters')
         import c_long_long, c_ptr
         type(c_ptr), value :: context
         integer(c_long_long), intent(out) :: cuda_malloc_count, cuda_free_count, workspace_query_count
         integer(c_long_long), intent(out) :: workspace_reuse_count, event_create_count, event_destroy_count
         integer(c_long_long), intent(out) :: pinned_alloc_count, pinned_free_count
      end subroutine rslmto_reciprocal_cuda_get_resource_counters

      function rslmto_reciprocal_cuda_synchronize(context) bind(C, name='rslmto_reciprocal_cuda_synchronize')
         import c_int, c_ptr
         type(c_ptr), value :: context
         integer(c_int) :: rslmto_reciprocal_cuda_synchronize
      end function rslmto_reciprocal_cuda_synchronize

      subroutine rslmto_reciprocal_cuda_destroy(context) bind(C, name='rslmto_reciprocal_cuda_destroy')
         import c_ptr
         type(c_ptr), value :: context
      end subroutine rslmto_reciprocal_cuda_destroy
   end interface
#endif


   interface
   subroutine backend_initialize_if(this, assembler)
      import reciprocal_execution_backend, reciprocal_assembler
      class(reciprocal_execution_backend), intent(inout) :: this
      type(reciprocal_assembler), intent(in) :: assembler
   end subroutine backend_initialize_if

   subroutine backend_capabilities_if(this, capabilities)
      import reciprocal_execution_backend, reciprocal_execution_capabilities
      class(reciprocal_execution_backend), intent(in) :: this
      type(reciprocal_execution_capabilities), intent(out) :: capabilities
   end subroutine backend_capabilities_if

   subroutine backend_prepare_if(this, operator_generation)
      import reciprocal_execution_backend
      class(reciprocal_execution_backend), intent(inout) :: this
      integer, intent(in) :: operator_generation
   end subroutine backend_prepare_if

   subroutine backend_execute_if(this, request, result)
      import reciprocal_execution_backend, reciprocal_execution_request, reciprocal_execution_result
      class(reciprocal_execution_backend), intent(inout) :: this
      type(reciprocal_execution_request), intent(in) :: request
      type(reciprocal_execution_result), intent(inout) :: result
   end subroutine backend_execute_if

   subroutine backend_execution_metrics_if(this, execute_requests, combined_requests, assemble_only, input_hamiltonian_solves)
      import reciprocal_execution_backend
      class(reciprocal_execution_backend), intent(in) :: this
      integer, intent(out) :: execute_requests, combined_requests, assemble_only, input_hamiltonian_solves
   end subroutine backend_execution_metrics_if

   subroutine backend_synchronize_if(this)
      import reciprocal_execution_backend
      class(reciprocal_execution_backend), intent(inout) :: this
   end subroutine backend_synchronize_if

   subroutine backend_release_if(this)
      import reciprocal_execution_backend
      class(reciprocal_execution_backend), intent(inout) :: this
   end subroutine backend_release_if

   module subroutine lapack_backend_initialize(this, assembler)
      class(lapack_reciprocal_backend), intent(inout) :: this
      type(reciprocal_assembler), intent(in) :: assembler
   end subroutine lapack_backend_initialize

   module subroutine lapack_backend_capabilities(this, capabilities)
      class(lapack_reciprocal_backend), intent(in) :: this
      type(reciprocal_execution_capabilities), intent(out) :: capabilities
   end subroutine lapack_backend_capabilities

   module subroutine lapack_backend_prepare_operator(this, operator_generation)
      class(lapack_reciprocal_backend), intent(inout) :: this
      integer, intent(in) :: operator_generation
   end subroutine lapack_backend_prepare_operator

   module subroutine lapack_backend_execute_batch(this, request, result)
      class(lapack_reciprocal_backend), intent(inout) :: this
      type(reciprocal_execution_request), intent(in) :: request
      type(reciprocal_execution_result), intent(inout) :: result
   end subroutine lapack_backend_execute_batch

   module subroutine lapack_backend_execution_metrics(this, execute_requests, combined_requests, assemble_only, input_hamiltonian_solves)
      class(lapack_reciprocal_backend), intent(in) :: this
      integer, intent(out) :: execute_requests, combined_requests, assemble_only, input_hamiltonian_solves
   end subroutine lapack_backend_execution_metrics

   module subroutine lapack_backend_synchronize(this)
      class(lapack_reciprocal_backend), intent(inout) :: this
   end subroutine lapack_backend_synchronize

   module subroutine lapack_backend_release(this)
      class(lapack_reciprocal_backend), intent(inout) :: this
   end subroutine lapack_backend_release

   module subroutine lapack_backend_destructor(this)
      type(lapack_reciprocal_backend), intent(inout) :: this
   end subroutine lapack_backend_destructor

   module subroutine cuda_backend_initialize(this, assembler)
      class(cuda_reciprocal_backend), intent(inout) :: this
      type(reciprocal_assembler), intent(in) :: assembler
   end subroutine cuda_backend_initialize

   module subroutine cuda_backend_capabilities(this, capabilities)
      class(cuda_reciprocal_backend), intent(in) :: this
      type(reciprocal_execution_capabilities), intent(out) :: capabilities
   end subroutine cuda_backend_capabilities

   module subroutine cuda_backend_prepare_operator(this, operator_generation)
      class(cuda_reciprocal_backend), intent(inout) :: this
      integer, intent(in) :: operator_generation
   end subroutine cuda_backend_prepare_operator

   module subroutine cuda_backend_execute_batch(this, request, result)
      class(cuda_reciprocal_backend), intent(inout) :: this
      type(reciprocal_execution_request), intent(in) :: request
      type(reciprocal_execution_result), intent(inout) :: result
   end subroutine cuda_backend_execute_batch

   module function cuda_backend_resident_eigensystem_matches(this, n, batch_size, token) result(matches)
      class(cuda_reciprocal_backend), intent(in) :: this
      integer, intent(in) :: n, batch_size, token
      logical :: matches
   end function cuda_backend_resident_eigensystem_matches

   module subroutine cuda_backend_contract_lehmann(this, request, result, status)
      class(cuda_reciprocal_backend), intent(inout) :: this
      type(reciprocal_lehmann_request), intent(in) :: request
      type(reciprocal_lehmann_result), intent(inout) :: result
      integer, intent(out) :: status
   end subroutine cuda_backend_contract_lehmann

   module subroutine cuda_backend_contract_lehmann_resident(this, request, result, status)
      class(cuda_reciprocal_backend), intent(inout) :: this
      type(reciprocal_lehmann_request), intent(in) :: request
      type(reciprocal_lehmann_result), intent(inout) :: result
      integer, intent(out) :: status
   end subroutine cuda_backend_contract_lehmann_resident

   module subroutine cuda_backend_execution_metrics(this, execute_requests, combined_requests, assemble_only, input_hamiltonian_solves)
      class(cuda_reciprocal_backend), intent(in) :: this
      integer, intent(out) :: execute_requests, combined_requests, assemble_only, input_hamiltonian_solves
   end subroutine cuda_backend_execution_metrics

   module subroutine cuda_backend_synchronize(this)
      class(cuda_reciprocal_backend), intent(inout) :: this
   end subroutine cuda_backend_synchronize

   module subroutine cuda_backend_release(this)
      class(cuda_reciprocal_backend), intent(inout) :: this
   end subroutine cuda_backend_release

   module subroutine cuda_backend_destructor(this)
      type(cuda_reciprocal_backend), intent(inout) :: this
   end subroutine cuda_backend_destructor

   module subroutine make_execution_backend(this, backend_name)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in), optional :: backend_name
   end subroutine make_execution_backend

   module subroutine reciprocal_workspace_ensure_capacity(this, nmat, tile_length, generalized, operator_generation, nnmax, ntype)
      class(reciprocal_workspace), intent(inout) :: this
      integer, intent(in) :: nmat, tile_length, operator_generation, nnmax, ntype
      logical, intent(in) :: generalized
   end subroutine reciprocal_workspace_ensure_capacity

   module subroutine reciprocal_workspace_clear(this)
      class(reciprocal_workspace), intent(inout) :: this
   end subroutine reciprocal_workspace_clear

   module subroutine reciprocal_workspace_restore_to_default(this)
      class(reciprocal_workspace), intent(inout) :: this
   end subroutine reciprocal_workspace_restore_to_default

   module subroutine reciprocal_workspace_destructor(this)
      type(reciprocal_workspace) :: this
   end subroutine reciprocal_workspace_destructor

   module subroutine reciprocal_assembler_assemble_batch(this, k_points, workspace)
      class(reciprocal_assembler), intent(in) :: this
      real(rp), intent(in) :: k_points(:, :)
      class(reciprocal_workspace), intent(inout) :: workspace
   end subroutine reciprocal_assembler_assemble_batch

   module subroutine reciprocal_assembler_assemble_overlap_batch(this, k_points, workspace)
      class(reciprocal_assembler), intent(in) :: this
      real(rp), intent(in) :: k_points(:, :)
      class(reciprocal_workspace), intent(inout) :: workspace
   end subroutine reciprocal_assembler_assemble_overlap_batch

   module subroutine reciprocal_assembler_assemble_one(this, k_point, hk_result, workspace)
      class(reciprocal_assembler), intent(in) :: this
      real(rp), intent(in) :: k_point(3)
      complex(rp), intent(out) :: hk_result(:, :)
      class(reciprocal_workspace), intent(inout) :: workspace
   end subroutine reciprocal_assembler_assemble_one

   module subroutine reciprocal_assembler_assemble_gbt(this, k_points, workspace)
      class(reciprocal_assembler), intent(in) :: this
      real(rp), intent(in) :: k_points(:, :)
      class(reciprocal_workspace), intent(inout) :: workspace
   end subroutine reciprocal_assembler_assemble_gbt

   module subroutine make_reciprocal_assembler(this, assembler)
      class(reciprocal), intent(inout) :: this
      type(reciprocal_assembler), intent(out) :: assembler
   end subroutine make_reciprocal_assembler

   !> @brief Print a root-rank informational message with source location.
   !> @param[in] message Message text.
   !> @param[in] file_name Source file name for diagnostics.
   !> @param[in] line_no Source line number for diagnostics.
   module subroutine root_info(message, file_name, line_no)
      character(len=*), intent(in) :: message, file_name
      integer, intent(in) :: line_no
   end subroutine root_info

   !> @brief Construct a reciprocal-space object from an initialized Hamiltonian.
   !> @details Wires reciprocal state to the Hamiltonian, lattice, and control
   !>          objects, restores defaults, and reads the &reciprocal namelist.
   !> @param[in] hamiltonian_obj Hamiltonian object providing real-space blocks and lattice state.
   !> @return Initialized reciprocal-space object.
   module function constructor(hamiltonian_obj) result(obj)
      type(reciprocal) :: obj
      type(hamiltonian), target, intent(in) :: hamiltonian_obj
   end function constructor

   !> @brief Finalize a reciprocal-space object.
   !> @details Releases k-point, reciprocal-vector, Hamiltonian, overlap, band,
   !>          DOS, tetrahedron, symmetry-map, and projection arrays.
   !> @param[inout] this Reciprocal object being finalized.
   module subroutine destructor(this)
      type(reciprocal) :: this
   end subroutine destructor

   !> @brief Reset reciprocal-space settings and owned storage to defaults.
   !> @details Restores k-mesh, solver-mode, DOS, band-path, symmetry, and
   !>          diagnostic options to baseline values and clears allocatables.
   !> @param[inout] this Reciprocal object to reset.
   module subroutine restore_to_default(this)
      class(reciprocal), intent(inout) :: this
   end subroutine restore_to_default

   !> @brief Read the &reciprocal namelist and install reciprocal-space options.
   !> @details Parses k-mesh, Fourier mode, band/DOS controls, symmetry options,
   !>          tetrahedron settings, and diagnostics from this%control%fname.
   !> @param[inout] this Reciprocal object whose input-facing options are populated.
   !> @note This is an input boundary and may raise fatal diagnostics for invalid options.
   module subroutine build_from_file(this)
      class(reciprocal), intent(inout) :: this
   end subroutine build_from_file

   !> @brief Generate reciprocal lattice vectors from the real-space lattice.
   !> @details Computes the reciprocal basis and reciprocal-cell volume from
   !>          lattice%a and lattice%alat for later k-point and Fourier work.
   !> @param[inout] this Reciprocal object receiving reciprocal_vectors and volume.
   module subroutine generate_reciprocal_vectors(this)
      class(reciprocal), intent(inout) :: this
   end subroutine generate_reciprocal_vectors

   !> @brief Generate the configured Monkhorst-Pack k-point mesh.
   !> @details Builds full or symmetry-reduced fractional k-points, weights, and
   !>          optional MPI ownership metadata according to reciprocal settings.
   !> @param[inout] this Reciprocal object receiving k_points/k_weights state.
   module subroutine generate_mp_mesh(this)
      class(reciprocal), intent(inout) :: this
   end subroutine generate_mp_mesh

   !> @brief Determine per-type basis sizes for reciprocal matrices.
   !> @details Maps each lattice atom type to sp/spd/spdf orbital counts used
   !>          when packing site-major k-space Hamiltonian and projection blocks.
   !> @param[inout] this Reciprocal object whose basis_size/max_orbs are updated.
   module subroutine set_basis_sizes(this)
      class(reciprocal), intent(inout) :: this
   end subroutine set_basis_sizes

   !> @brief Build real-space neighbor vectors used by Fourier transforms.
   !> @details Converts lattice neighbor geometry into Cartesian and direct
   !>          per-type R-vector tables aligned with Hamiltonian neighbor slots.
   !> @param[inout] this Reciprocal object receiving ham_vec_type arrays.
   module subroutine build_neighbor_vectors(this)
      class(reciprocal), intent(inout) :: this
   end subroutine build_neighbor_vectors

   !> @brief Calculate exp(i k.R) factors for each neighbor/type.
   !> @param[in] this Reciprocal object containing neighbor-vector tables.
   !> @param[in] k_vec k-point vector in reciprocal coordinates.
   !> @param[out] structure_factors Phase factors indexed by neighbor and atom type.
   module subroutine calculate_structure_factors(this, k_vec, structure_factors)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: structure_factors  ! (nr_max, ntype)
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
      ! First-order k-space Hamiltonian: H(k) = h(k) = Sum_R ee(R) exp(i*k·R).
      ! NOTE: This deliberately reproduces the historical first-order behaviour.
      !       The on-site E_nu (enim) and spin-orbit (lsham) terms are NOT added
      !       here; they are only included in the second-order path
      !       (fourier_transform_hamiltonian_second_order). See kspace_ham_order.
   end subroutine fourier_transform_hamiltonian

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
   end subroutine fourier_transform_overlap

   !> @brief Fold a fractional reciprocal-space point into [-1/2,1/2).
   !> @details Folding is component-wise in reciprocal lattice units, so k and
   !>          k+G are represented by exactly the same Bloch Hamiltonian.
   !> @param[in] this Reciprocal object (kept for a uniform type-bound API).
   !> @param[in] k_point Arbitrary fractional reciprocal-space point.
   !> @param[out] k_folded Equivalent point in the reciprocal primitive cell.
   module subroutine fold_kpoint(this, k_point, k_folded)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: k_point(3)
      real(rp), intent(out) :: k_folded(3)
   end subroutine fold_kpoint

   !> @brief Construct the normal-state reciprocal Hamiltonian at one arbitrary k point.
   !> @details This does not touch k_points, hk_bulk, or any bands/DOS cache.
   !>          It uses the same first/second-order selection and Fourier terms as
   !>          build_kspace_hamiltonian.
   module subroutine build_hamiltonian_at_kpoint(this, k_point, hk_result)
      class(reciprocal), intent(inout) :: this
      real(rp), intent(in) :: k_point(3)
      complex(rp), intent(out) :: hk_result(:, :)
   end subroutine build_hamiltonian_at_kpoint

   !> @brief Return caller-owned eigenpairs at arbitrary reciprocal-space points.
   !> @details Points are folded into the reciprocal primitive cell, repeated
   !>          folded points are solved once, and no standard mesh/bands/DOS
   !>          state is overwritten.  The routine is intentionally local to the
   !>          caller, making it suitable for an outer MPI-over-q decomposition.
   !> @param[in] this Reciprocal object with a completed normal-state Hamiltonian.
   !> @param[in] k_points Arbitrary fractional reciprocal points, shape (3,nk).
   !> @param[out] eigenvalues Caller-owned array, shape (nband,nk).
   !> @param[out] eigenvectors Caller-owned array, shape (nbasis,nband,nk).
   !> @param[out] folded_k_points Optional caller-owned folded points, shape (3,nk).
   module subroutine calculate_eigenpairs_at_kpoints(this, k_points, eigenvalues, eigenvectors, folded_k_points)
      class(reciprocal), intent(inout) :: this
      real(rp), intent(in) :: k_points(:, :)
      real(rp), allocatable, intent(out) :: eigenvalues(:, :)
      complex(rp), allocatable, intent(out) :: eigenvectors(:, :, :)
      real(rp), allocatable, intent(out), optional :: folded_k_points(:, :)
   end subroutine calculate_eigenpairs_at_kpoints

   !> Reject a distributed mesh at consumers whose mathematics still requires
   !> every rank to own every k point (tetrahedron, reciprocal Green, moments).
   module subroutine require_replicated_k_workset(this, consumer)
      class(reciprocal), intent(in) :: this
      character(len=*), intent(in) :: consumer
   end subroutine require_replicated_k_workset

   !> Collect per-k columns of a k-distributed array into the full mesh on every rank (plain copy if not distributed).
   !> @param[in] ld Number of values per k point.
   !> @param[in] local Local columns, ld x nk_local, any rank.
   !> @param[inout] full Full mesh, ld x nk_global, any rank.
   module subroutine gather_k_columns(this, ld, local, full)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: ld
      real(rp), intent(in) :: local(*)
      real(rp), intent(inout) :: full(*)
   end subroutine gather_k_columns

   !> @brief Configure MPI ownership for a k-point set.
   !> @param[inout] this Reciprocal object receiving local/global k maps.
   !> @param[in] nk_global Number of global k-points.
   !> @param[in] enable_distribution If true, distribute k-points over MPI ranks.
   module subroutine setup_k_mesh_distribution(this, nk_global, enable_distribution)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: nk_global
      logical, intent(in) :: enable_distribution
   end subroutine setup_k_mesh_distribution

   !> @brief Convert a local k-point index to a global k-point index.
   !> @param[in] this Reciprocal object containing distribution maps.
   !> @param[in] ik_local Local k-point index.
   !> @return Global k-point index.
   module integer function local_k_index_to_global(this, ik_local) result(ik_global)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: ik_local
   end function local_k_index_to_global

   !> @brief Build H(k) for every active mesh or path k-point.
   !> @details Allocates hk_bulk/hk_total as needed, dispatches first- or
   !>          second-order Fourier transforms, and applies local k ownership.
   !> @param[inout] this Reciprocal object receiving k-space Hamiltonian arrays.
   module subroutine build_kspace_hamiltonian(this)
      class(reciprocal), intent(inout) :: this
   end subroutine build_kspace_hamiltonian

   !> @brief Assemble and solve the active reciprocal mesh one tile at a time.
   !> @details The host compatibility caches are copied only from a completed
   !> execution result.  The typed request deliberately combines Fourier
   !> assembly and eigensolution so a resident backend never needs a host H(k)
   !> round trip between those operations.
   module subroutine execute_normal_mesh_tiles(this, using_kpath, generalized)
      class(reciprocal), intent(inout) :: this
      logical, intent(in) :: using_kpath, generalized
   end subroutine execute_normal_mesh_tiles

   !> @brief Build S(k) for every active mesh or path k-point.
   !> @param[inout] this Reciprocal object receiving sk_overlap arrays.
   module subroutine build_kspace_overlap(this)
      class(reciprocal), intent(inout) :: this
   end subroutine build_kspace_overlap

!> @brief Diagnose onsite spin-diagonal blocks for multisite Hamiltonians.
!> @param[in] this Reciprocal object containing hk_total or hk_bulk data.
module subroutine check_multisite_hamiltonian_diagonal(this)
   class(reciprocal), intent(in) :: this
end subroutine check_multisite_hamiltonian_diagonal

   !> @brief Diagonalize the active k-space Hamiltonians.
   !> @details Solves standard or generalized eigenproblems depending on
   !>          reciprocal_mode and stores eigenvalues/eigenvectors for bands and DOS.
   !> @param[inout] this Reciprocal object receiving eigenvalue/eigenvector arrays.
   module subroutine diagonalize_hamiltonian(this)
      class(reciprocal), intent(inout) :: this
   end subroutine diagonalize_hamiltonian

   !> @brief Return whether the legacy reciprocal GBT path has a finite q vector.
   !> @details WP0 treats this state specially: only the bare first-order
   !>          operator is permitted while the bond-gauge implementation is pending.
   module logical function has_nonzero_q_gbt(this)
      class(reciprocal), intent(in) :: this
   end function has_nonzero_q_gbt

   !> @brief Reject nonzero-q GBT combinations whose operator transformation is unverified.
   !> @param[in] context Caller name included in fatal diagnostics.
   module subroutine validate_nonzero_q_gbt(this, context)
      class(reciprocal), intent(in) :: this
      character(len=*), intent(in) :: context
   end subroutine validate_nonzero_q_gbt

   !> @brief Disable symmetry/time-reversal reduction and rebuild a full chemical BZ mesh for finite-q GBT.
   module subroutine force_full_bz_for_nonzero_q_gbt(this, context)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: context
   end subroutine force_full_bz_for_nonzero_q_gbt

   !> @brief Print mapping diagnostics for the Kanpur generalized-overlap path.
   !> @param[in] this Reciprocal object containing basis and overlap metadata.
   module subroutine print_kanpur_mapping(this)
      class(reciprocal), intent(in) :: this
   end subroutine print_kanpur_mapping

   !> @brief Check Hermiticity and basic diagnostics for an overlap matrix.
   !> @param[in] this Reciprocal object providing diagnostic context.
   !> @param[in] ik k-point index associated with s_k.
   !> @param[in] s_k Overlap matrix to inspect.
   module subroutine check_overlap_properties(this, ik, s_k)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: ik
      complex(rp), dimension(:, :), intent(in) :: s_k
   end subroutine check_overlap_properties

   !> @brief Run H(Gamma) spectral-bound diagnostics.
   !> @param[inout] this Reciprocal object used to build and bound the Gamma matrix.
   module subroutine run_gamma_bounds_diagnostics(this)
      class(reciprocal), intent(inout) :: this
   end subroutine run_gamma_bounds_diagnostics

   !> @brief Diagonalize the experimental real-space HALL matrix.
   !> @details Builds a finite local-cluster matrix from HALL blocks and prints
   !>          eigenvalue diagnostics for exploratory comparison only.
   !> @param[inout] this Reciprocal object providing Hamiltonian and lattice state.
   module subroutine diagonalize_hall_experimental(this)
      class(reciprocal), intent(inout) :: this
   end subroutine diagonalize_hall_experimental

   !> @brief Calculate and write a band structure along a crystal-specific path.
   !> @param[inout] this Reciprocal object used to generate path eigenvalues.
   !> @param[in] ham Hamiltonian object for k-space construction.
   !> @param[in] crystal_type Crystal/path selector.
   !> @param[in] npts_per_segment Number of interpolation points per path segment.
   !> @param[in] output_file Optional band-output file name.
   module subroutine calculate_band_structure(this, ham, crystal_type, npts_per_segment, output_file)
      class(reciprocal), intent(inout) :: this
      class(hamiltonian), intent(in) :: ham
      character(len=*), intent(in) :: crystal_type
      integer, intent(in) :: npts_per_segment
      character(len=*), intent(in), optional :: output_file
   end subroutine calculate_band_structure

   !> @brief Calculate and write a symmetry-derived band structure.
   !> @details Uses symmetry analysis or custom path settings to generate the
   !>          high-symmetry path before building and diagonalizing H(k).
   !> @param[inout] this Reciprocal object used to generate path eigenvalues.
   !> @param[in] ham Hamiltonian object for k-space construction.
   !> @param[in] npts_per_segment Optional number of interpolation points per segment.
   !> @param[in] output_file Optional band-output file name.
   module subroutine calculate_band_structure_auto(this, ham, npts_per_segment, output_file)
      class(reciprocal), intent(inout) :: this
      class(hamiltonian), intent(in) :: ham
      integer, intent(in), optional :: npts_per_segment
      character(len=*), intent(in), optional :: output_file
   end subroutine calculate_band_structure_auto

   !> @brief Generate a symmetry-reduced k-point mesh.
   !> @details Builds irreducible k-points, weights, and full/irreducible maps
   !>          using the configured symmetry and time-reversal settings.
   !> @param[inout] this Reciprocal object receiving reduced-mesh data.
   !> @param[in] mesh_dims Full mesh dimensions.
   !> @param[in] use_shift Optional flag selecting shifted-grid generation.
   module subroutine generate_reduced_kpoint_mesh(this, mesh_dims, use_shift)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: mesh_dims(3)
      logical, intent(in), optional :: use_shift
   end subroutine generate_reduced_kpoint_mesh

   !> @brief Build and store the symmetry group of the current site moments.
   !> @details The group is fixed from the moments at the time of the call; call again to rebuild it if the
   !>          moments rotate away from it.
   module subroutine setup_magnetic_symmetry(this)
      class(reciprocal), intent(inout) :: this
   end subroutine setup_magnetic_symmetry

   !> @brief Average this%rf_density over the symmetry group (no-op without a group).
   module subroutine symmetrize_spin_density(this)
      class(reciprocal), intent(inout) :: this
   end subroutine symmetrize_spin_density

   !> @brief Average one spin vector per site over the symmetry group (no-op without a group).
   module subroutine symmetrize_site_moments(this, v)
      class(reciprocal), intent(in) :: this
      real(rp), intent(inout) :: v(:, :)
   end subroutine symmetrize_site_moments

   !> @brief Average per-site density matrices over the symmetry group (no-op without a group).
   module subroutine symmetrize_site_dm(this, rho)
      class(reciprocal), intent(in) :: this
      complex(rp), intent(inout) :: rho(:, :, :)
   end subroutine symmetrize_site_dm

   !> @brief Calculate reciprocal-space density of states.
   !> @details Configures the energy grid and DOS method, ensures eigenvalues are
   !>          available, evaluates total/projected DOS, moments, and writes output.
   !> @param[inout] this Reciprocal object receiving DOS and moment arrays.
   !> @param[in] ham Hamiltonian object for k-space construction if needed.
   !> @param[in] n_energy_points Optional number of DOS grid points.
   !> @param[in] energy_range Optional energy window.
   !> @param[in] method Optional DOS method selector.
   !> @param[in] gaussian_sigma Optional Gaussian broadening.
   !> @param[in] temperature Optional Fermi-Dirac temperature.
   !> @param[in] fermi_level Optional fixed Fermi level.
   !> @param[in] total_electrons Optional electron count for automatic Fermi search.
   !> @param[in] auto_find_fermi Optional flag enabling Fermi-level search.
   !> @param[in] output_file Optional DOS output file.
   module subroutine calculate_density_of_states(this, ham, n_energy_points, energy_range, method, gaussian_sigma, temperature, fermi_level, total_electrons, auto_find_fermi, output_file)
      class(reciprocal), intent(inout) :: this
      class(hamiltonian), intent(in) :: ham
      integer, intent(in), optional :: n_energy_points
      real(rp), dimension(2), intent(in), optional :: energy_range
      character(len=*), intent(in), optional :: method
      real(rp), intent(in), optional :: gaussian_sigma
      real(rp), intent(in), optional :: temperature
      real(rp), intent(in), optional :: fermi_level
      real(rp), intent(in), optional :: total_electrons
      logical, intent(in), optional :: auto_find_fermi
      character(len=*), intent(in), optional :: output_file
   end subroutine calculate_density_of_states

   module subroutine prepare_density_of_states(this, ham, n_energy_points, energy_range, method, gaussian_sigma, temperature, fermi_level, total_electrons, auto_find_fermi)
      class(reciprocal), intent(inout) :: this
      class(hamiltonian), intent(in) :: ham
      integer, intent(in), optional :: n_energy_points
      real(rp), dimension(2), intent(in), optional :: energy_range
      character(len=*), intent(in), optional :: method
      real(rp), intent(in), optional :: gaussian_sigma
      real(rp), intent(in), optional :: temperature
      real(rp), intent(in), optional :: fermi_level
      real(rp), intent(in), optional :: total_electrons
      logical, intent(in), optional :: auto_find_fermi
   end subroutine prepare_density_of_states

   module subroutine calculate_dos_output(this, output_file, with_moments)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in), optional :: output_file
      logical, intent(in), optional :: with_moments
   end subroutine calculate_dos_output

   !> @brief Invalidate eigensystem, DOS, projection, and canonical-energy state.
   !> @details K-point geometry and weights are retained so a frozen-potential
   !>          probe can rebuild on the same mesh without reusing stale spectra.
   module subroutine invalidate_spectral_cache(this)
      class(reciprocal), intent(inout) :: this
   end subroutine invalidate_spectral_cache

   !> @brief Drop reciprocal caches when the shared real-space operator changed.
   !> @param[inout] this Reciprocal object owning operator-derived caches.
   !> @param[in] context_tag Caller label used in the invalidation diagnostic.
   !> @param[out] changed Optional flag reporting a generation mismatch.
   module subroutine invalidate_if_operator_changed(this, context_tag, changed)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: context_tag
      logical, intent(out), optional :: changed
   end subroutine invalidate_if_operator_changed

   !> @brief Validate full-to-irreducible k-point symmetry maps.
   !> @param[inout] this Reciprocal object containing symmetry maps and weights.
   !> @param[in] context_tag Diagnostic label for the caller/context.
   module subroutine validate_symmetry_kmap(this, context_tag)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: context_tag
   end subroutine validate_symmetry_kmap

   !> @brief Write symmetry k-point maps for diagnostics.
   !> @param[inout] this Reciprocal object containing symmetry maps.
   !> @param[in] filename Output file name.
   module subroutine write_symmetry_kmap_dump(this, filename)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: filename
   end subroutine write_symmetry_kmap_dump

   !> @brief Average the per-site, per-l projected DOS curves over the symmetry group (no-op without a group).
   !> @details The charge curve n is permuted over sites, the spin-vector curves projected_dos_moments are also acted on by
   !>          sym_spin_matrices. projected_dos (along dos_site_axes) and the totals dos_mx/my/mz_tot are rebuilt from them.
   module subroutine symmetrize_projected_dos(this)
      class(reciprocal), intent(inout) :: this
   end subroutine symmetrize_projected_dos

   !> @brief Operation that maps the irreducible point of each full-mesh point onto it, and the spin matrix of each operation.
   !> @details Searches the stored group, with time reversal for nsp = 1 (the sign is sym_spin_flips(op)); stops if none
   !>          maps the irreducible point onto the full-mesh point.
   !> @param[out] op Operation index of each full-mesh point [nk_full].
   !> @param[out] w Cartesian action of each stored operation on a spin vector [3, 3, n_sym_ops].
   module subroutine full_mesh_operations(this, op, w)
      class(reciprocal), intent(in) :: this
      integer, allocatable, intent(out) :: op(:)
      real(rp), allocatable, intent(out) :: w(:, :, :)
   end subroutine full_mesh_operations

   !> @brief Band projections of one full-mesh point from the irreducible mesh.
   !> @details The irreducible point of kf is rotated by the operation that maps it onto kf: the site is permuted, the spin
   !>          vector is acted on by sym_spin_matrices. An empty op means the mesh is not reduced and kf is irreducible.
   !> @param[in] proj_irr Projections (total, mx, my, mz; l; site; band; irreducible k).
   !> @param[in] op Operation index of each full-mesh point, from full_mesh_operations.
   !> @param[in] w Spin matrix of each operation, from full_mesh_operations.
   !> @param[in] kf Full-mesh point.
   !> @param[in] l Orbital type.
   !> @param[in] isite Site on the full mesh.
   !> @param[in] ib Band.
   !> @return (total, mx, my, mz) projection.
   module function full_mesh_projection(this, proj_irr, op, w, kf, l, isite, ib) result(p)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: proj_irr(:, :, :, :, :), w(:, :, :)
      integer, intent(in) :: op(:), kf, l, isite, ib
      real(rp) :: p(4)
   end function full_mesh_projection

   !> @brief Build irreducible tetrahedra and multiplicities from full-mesh cuts.
   !> @param[inout] this Reciprocal object containing symmetry maps and mesh dimensions.
   !> @param[out] tet_ir Irreducible tetrahedron corner indices.
   !> @param[out] tet_mult Multiplicity of each irreducible tetrahedron.
   !> @param[out] n_tet_ir Number of irreducible tetrahedra.
   module subroutine build_irreducible_tetrahedra(this, tet_ir, tet_mult, n_tet_ir)
      class(reciprocal), intent(in) :: this
      integer, allocatable, intent(out) :: tet_ir(:, :)
      integer, allocatable, intent(out) :: tet_mult(:)
      integer, intent(out) :: n_tet_ir
   end subroutine build_irreducible_tetrahedra

   !> @brief Bloechl-corrected tetrahedron occupations at the current Fermi level, occ(band, global k); dos_method='blochl' only.
   module subroutine blochl_occupations(this, occ)
      class(reciprocal), intent(inout) :: this
      real(rp), allocatable, intent(out) :: occ(:, :)
   end subroutine blochl_occupations

   !> @brief Select tetrahedron corner offsets for the current mesh geometry.
   !> @param[in] this Reciprocal object providing reciprocal lattice geometry.
   !> @param[out] tetra_cut Six tetrahedra per mesh cell, stored as corner offsets.
   module subroutine get_tetra_cut_offsets(this, tetra_cut)
      class(reciprocal), intent(in) :: this
      integer, intent(out) :: tetra_cut(3, 4, 6)
   end subroutine get_tetra_cut_offsets

   !> @brief Build the DOS energy grid from the configured range and size.
   !> @param[inout] this Reciprocal object receiving dos_energy_grid.
   module subroutine setup_dos_energy_grid(this)
      class(reciprocal), intent(inout) :: this
   end subroutine setup_dos_energy_grid

   !> @brief Calculate total DOS using the tetrahedron method.
   !> @details Accumulates DOS and integrated NOS over full or irreducible
   !>          tetrahedra, respecting MPI partitioning where active.
   !> @param[inout] this Reciprocal object receiving total_dos and total_nos.
   module subroutine calculate_dos_tetrahedron(this)
      class(reciprocal), intent(inout) :: this
   end subroutine calculate_dos_tetrahedron

   !> @brief Sort four real values in ascending order.
   !> @param[in] arr_in Four input values.
   !> @param[out] arr_out Sorted output values.
   module subroutine sort4(arr_in, arr_out)
      real(rp), dimension(4), intent(in) :: arr_in
      real(rp), dimension(4), intent(out) :: arr_out
   end subroutine sort4

   !> @brief Add one tetrahedron’s integrated-state contribution to a grid.
   !> @param[in] volwgt Tetrahedron volume/weight factor.
   !> @param[in] ecorn_in Tetrahedron corner energies.
   !> @param[in] emin Lower grid energy.
   !> @param[in] emax Upper grid energy.
   !> @param[inout] nos Integrated number-of-states grid.
   !> @param[in] npts Number of grid points.
   module subroutine tetra_add_nos(volwgt, ecorn_in, emin, emax, nos, npts)
      real(rp), intent(in) :: volwgt, ecorn_in(4), emin, emax
      integer, intent(in) :: npts
      real(rp), intent(inout) :: nos(npts)
   end subroutine tetra_add_nos

   !> @brief Add one tetrahedron’s DOS contribution to a grid.
   !> @param[in] volwgt Tetrahedron volume/weight factor.
   !> @param[in] ecorn_in Tetrahedron corner energies.
   !> @param[in] emin Lower grid energy.
   !> @param[in] emax Upper grid energy.
   !> @param[inout] dos DOS grid.
   !> @param[in] npts Number of grid points.
   module subroutine tetra_add_dos(volwgt, ecorn_in, emin, emax, dos, npts)
      real(rp), intent(in) :: volwgt, ecorn_in(4), emin, emax
      integer, intent(in) :: npts
      real(rp), intent(inout) :: dos(npts)
   end subroutine tetra_add_dos

   !> @brief Per-corner DOS weights of one linear tetrahedron and the Bloechl curvature correction at one energy.
   !> @details Per unit tetrahedron volume. w_lin is the energy derivative of the uncorrected weights of PRB 49, 16223
   !>          App. B; dw_corr is the energy derivative of the correction Eq. (22).
   !> @param[in] ecorn Corner energies in caller order.
   !> @param[in] energy Energy at which the weights are evaluated.
   !> @param[out] w_lin Linear-tetrahedron corner weights in caller order.
   !> @param[out] dw_corr Curvature-correction corner weights in caller order.
   module subroutine tetra_corner_weights(ecorn, energy, w_lin, dw_corr)
      real(rp), intent(in) :: ecorn(4), energy
      real(rp), intent(out) :: w_lin(4), dw_corr(4)
   end subroutine tetra_corner_weights

   !> @brief Integrated corner weights (PRB 49, 16223 App. B) and the Eq. (22) correction of one tetrahedron at one energy.
   !> @param[in] ecorn Corner energies in caller order.
   !> @param[in] energy Energy at which the weights are evaluated.
   !> @param[out] w_lin Linear-tetrahedron integration weights in caller order.
   !> @param[out] w_corr Curvature-correction weights in caller order.
   module subroutine tetra_integrated_weights(ecorn, energy, w_lin, w_corr)
      real(rp), intent(in) :: ecorn(4), energy
      real(rp), intent(out) :: w_lin(4), w_corr(4)
   end subroutine tetra_integrated_weights

!> @brief Calculate total DOS with Gaussian broadening.
!> @details Smears all eigenvalues over the DOS energy grid using the configured
!>          or adaptive sigma and k-point weights.
!> @param[inout] this Reciprocal object receiving total_dos.
module subroutine calculate_dos_gaussian(this)
   class(reciprocal), intent(inout) :: this
end subroutine calculate_dos_gaussian

   !> @brief Calculate total DOS with the Bloechl tetrahedron method.
   !> @param[inout] this Reciprocal object receiving total_dos.
   module subroutine calculate_dos_blochl(this)
      class(reciprocal), intent(inout) :: this
   end subroutine calculate_dos_blochl

   !> @brief Build full-mesh tetrahedron connectivity.
   !> @param[inout] this Reciprocal object receiving tetrahedra and volumes.
   module subroutine setup_tetrahedra(this)
      class(reciprocal), intent(inout) :: this
   end subroutine setup_tetrahedra

   !> @brief Convert 3D mesh coordinates to a periodic flat k-point index.
   !> @param[in] this Reciprocal object providing helper context.
   !> @param[in] i Mesh coordinate along axis 1.
   !> @param[in] j Mesh coordinate along axis 2.
   !> @param[in] k Mesh coordinate along axis 3.
   !> @param[in] nk1 Mesh size along axis 1.
   !> @param[in] nk2 Mesh size along axis 2.
   !> @param[in] nk3 Mesh size along axis 3.
   !> @return One-based periodic flat k-point index.
   module function get_kpoint_index(this, i, j, k, nk1, nk2, nk3) result(idx)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: i, j, k, nk1, nk2, nk3
      integer :: idx, ii, jj, kk
   end function get_kpoint_index

   !> @brief Estimate a Gaussian DOS broadening from k-point density.
   !> @param[in] this Reciprocal object containing mesh and reciprocal volume.
   !> @return Adaptive Gaussian sigma.
   module function calculate_adaptive_sigma(this) result(sigma)
      class(reciprocal), intent(in) :: this
      real(rp) :: sigma
   end function calculate_adaptive_sigma

   !> @brief Dispatch projected-DOS calculation for the active DOS method.
   !> @param[inout] this Reciprocal object receiving projected DOS arrays.
   module subroutine project_dos_orbitals(this)
      class(reciprocal), intent(inout) :: this
   end subroutine project_dos_orbitals

!> @brief Calculate orbital- and spin-projected DOS with Gaussian broadening.
!> @details Projects eigenvector weights onto site/orbital/spin channels and
!>          accumulates both scalar and directional spin DOS components.
!> @param[inout] this Reciprocal object receiving projected DOS arrays.
module subroutine project_dos_orbitals_gaussian(this)
   class(reciprocal), intent(inout) :: this
end subroutine project_dos_orbitals_gaussian

   !> @brief Return the local spin axis for a projected-DOS site.
   !> @param[in] this Reciprocal object containing lattice magnetic moments.
   !> @param[in] isite Site index in the reciprocal packing.
   !> @param[out] axis Normalized local spin axis.
   module subroutine get_site_spin_axis(this, isite, axis)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: isite
      real(rp), intent(out) :: axis(3)
   end subroutine get_site_spin_axis

   !> @brief Project one spinor eigenvector block onto scalar and vector weights.
   !> @param[in] this Reciprocal object containing eigenvectors.
   !> @param[in] ik k-point index.
   !> @param[in] ib Band index.
   !> @param[in] site_orb_start First orbital index for the site block.
   !> @param[in] n_orb_per_spin Number of orbital channels per spin.
   !> @param[in] i_first First local orbital index in the projection range.
   !> @param[in] i_last Last local orbital index in the projection range.
   !> @param[out] total_char Scalar projected character.
   !> @param[out] mx_char Projected x spin moment character.
   !> @param[out] my_char Projected y spin moment character.
   !> @param[out] mz_char Projected z spin moment character.
   module subroutine compute_spinor_block_projection(this, ik, ib, site_orb_start, n_orb_per_spin, i_first, i_last, &
                                              total_char, mx_char, my_char, mz_char)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: ik, ib, site_orb_start, n_orb_per_spin, i_first, i_last
      real(rp), intent(out) :: total_char, mx_char, my_char, mz_char
   end subroutine compute_spinor_block_projection

   !> @brief Calculate projected DOS with tetrahedron integration.
   !> @details Averages orbital/spinor projection weights over tetrahedron
   !>          corners and accumulates scalar plus directional projected DOS.
   !> @param[inout] this Reciprocal object receiving projected DOS arrays.
   module subroutine project_dos_orbitals_tetrahedron(this)
      class(reciprocal), intent(inout) :: this
   end subroutine project_dos_orbitals_tetrahedron

   !> @brief Integrate tabulated data with the trapezoidal rule.
   !> @param[in] x Grid coordinates.
   !> @param[in] y Values on the grid.
   !> @return Trapezoidal integral over the grid.
   module function trapezoidal_integral(x, y) result(integral)
      real(rp), dimension(:), intent(in) :: x, y
      real(rp) :: integral
   end function trapezoidal_integral

   !> @brief Linearly interpolate a tabulated grid value.
   !> @param[in] x Grid coordinates.
   !> @param[in] y Values on the grid.
   !> @param[in] n Number of grid points.
   !> @param[in] x0 Coordinate where y is requested.
   !> @return Interpolated value at x0.
   module function interpolate_grid_value(x, y, n, x0) result(y0)
      integer, intent(in) :: n
      real(rp), intent(in) :: x(n), y(n), x0
      real(rp) :: y0
   end function interpolate_grid_value

!> @brief Calculate zeroth, first, and second DOS band moments.
!> @details Integrates projected DOS channels with Fermi weighting to produce
!>          occupation, first-energy, and second-energy moments.
!> @param[inout] this Reciprocal object receiving band_moments.
module subroutine calculate_band_moments(this)
   class(reciprocal), intent(inout) :: this
end subroutine calculate_band_moments

!> @brief WP7 k-space producer for the shared rotating-frame density contract.
!> @param[inout] this Reciprocal object whose `rf_density` is filled.
module subroutine accumulate_spin_density_kspace(this)
   class(reciprocal), intent(inout) :: this
end subroutine accumulate_spin_density_kspace

!> @brief Project the accumulated density onto radial band moments.
!> @param[inout] this      Reciprocal object receiving band_moments.
!> @param[in]    policy    SCF density policy governing the axis choice.
!> @param[in]    reference Per-site reference directions (3, n_sites).
!> @param[out]   axis_out  Per-site axis actually used (3, n_sites).
module subroutine fill_band_moments_from_spin_density(this, policy, reference, axis_out)
   class(reciprocal), intent(inout) :: this
   character(len=*), intent(in) :: policy
   real(rp), dimension(:, :), intent(in) :: reference
   real(rp), dimension(:, :), intent(out) :: axis_out
end subroutine fill_band_moments_from_spin_density

!> @brief Evaluate normalized Fermi-Dirac electron count and eigenvalue EBAND.
!> @details Each rank accumulates only its owned eigenvalues.  The three raw
!>          sums (weight, occupation, energy) are reduced together exactly once
!>          when the mesh is distributed, then occupation and energy are divided
!>          by the explicit global weight sum.
   module subroutine evaluate_eigenvalue_occupations(this, fermi_level, electron_count, band_energy, weight_sum)
      class(reciprocal), intent(inout) :: this
      real(rp), intent(in) :: fermi_level
      real(rp), intent(out) :: electron_count, band_energy
      real(rp), intent(out), optional :: weight_sum
   end subroutine evaluate_eigenvalue_occupations

!> @brief Find EF directly from eigenvalues, k weights, and Fermi occupations.
   module function find_fermi_level_from_eigenvalues(this, total_electrons) result(fermi_level)
      class(reciprocal), intent(inout) :: this
      real(rp), intent(in) :: total_electrons
      real(rp) :: fermi_level
   end function find_fermi_level_from_eigenvalues

!> @brief Set/validate EF and return the canonical occupied-eigenvalue EBAND.
   module function calculate_canonical_band_energy(this, find_fermi, electron_count) result(eband)
      class(reciprocal), intent(inout) :: this
      logical, intent(in), optional :: find_fermi
      real(rp), intent(out), optional :: electron_count
      real(rp) :: eband
   end function calculate_canonical_band_energy

!> @brief Integrate DOS occupation up to an energy with Fermi weighting.
!> @param[in] this Reciprocal object containing DOS grid data.
!> @param[in] energy Trial Fermi energy.
!> @param[in] kT Thermal energy used in Fermi weighting.
!> @return Integrated electron count.
module function integrate_dos_up_to_energy(this, energy, kT) result(integral)
   class(reciprocal), intent(in) :: this
   real(rp), intent(in) :: energy, kT
   real(rp) :: integral
end function integrate_dos_up_to_energy

   !> @brief Write total and projected DOS data to output files.
   !> @param[in] this Reciprocal object containing DOS arrays.
   !> @param[in] filename Primary total-DOS output file.
   module subroutine write_dos_to_file(this, filename)
      class(reciprocal), intent(in) :: this
      character(len=*), intent(in) :: filename
   end subroutine write_dos_to_file

   !> @brief Calculate band energy from integrated DOS moments.
   !> @param[in] this Reciprocal object containing band_moments.
   !> @return Band energy accumulated over projected channels.
   module function calculate_band_energy_from_moments(this) result(eband)
      class(reciprocal), intent(in) :: this
      real(rp) :: eband
   end function calculate_band_energy_from_moments

   !> @brief Diagnostic integral of E*D_total(E)*f(E) on the configured DOS grid.
   module function calculate_band_energy_from_total_dos(this) result(eband)
      class(reciprocal), intent(in) :: this
      real(rp) :: eband
   end function calculate_band_energy_from_total_dos

   !> @brief Reconstruct lattice local density matrices from projected DOS moments.
   !> @param[inout] this Reciprocal object containing projected DOS/moments.
   !> @param[inout] lattice_obj Lattice object receiving local density matrices.
   module subroutine calculate_ldm_from_projected_dos(this, lattice_obj, fixed_liechtenstein, m1_out, m1_only)
      use lattice_mod
      class(reciprocal), intent(inout) :: this
      type(lattice), intent(inout) :: lattice_obj
      logical, intent(in), optional :: fixed_liechtenstein
      complex(rp), intent(out), optional :: m1_out(:, :, :, :, :)
      logical, intent(in), optional :: m1_only
   end subroutine calculate_ldm_from_projected_dos

   !> @brief Reconstruct lattice local density matrices directly from eigenvectors.
   !> @details Accumulates occupied spinor density matrices over k-points and
   !>          bands, then synchronizes the lattice result across MPI ranks.
   !> @param[inout] this Reciprocal object containing eigenvectors and occupations.
   !> @param[inout] lattice_obj Lattice object receiving local density matrices.
   module subroutine calculate_ldm_from_eigenvectors(this, lattice_obj, fixed_liechtenstein, m1_out, m1_only)
      use lattice_mod
      class(reciprocal), intent(inout) :: this
      type(lattice), intent(inout) :: lattice_obj
      logical, intent(in), optional :: fixed_liechtenstein
      complex(rp), intent(out), optional :: m1_out(:, :, :, :, :)
      logical, intent(in), optional :: m1_only
   end subroutine calculate_ldm_from_eigenvectors

   !> Full spinor d-shell P and M1 per site from the current eigenpairs (diagnostic, nsp>=2).
   module subroutine calculate_spinor_dm(this, lattice_obj)
      use lattice_mod
      class(reciprocal), intent(inout) :: this
      type(lattice), intent(inout) :: lattice_obj
   end subroutine calculate_spinor_dm

   ! k-space Green’s-function engine (milestone B2, submodule reciprocal_green)

   !> @brief Build the retarded complex-energy contour z = E + i*green_eta.
   !> @param[in]  this      Reciprocal object (supplies green_eta broadening).
   !> @param[in]  en        Energy object holding the prepared real-axis grid.
   !> @param[out] z_contour Allocated retarded contour, size = size(en%ene).
   module subroutine build_green_contour(this, en, z_contour)
      class(reciprocal), intent(in) :: this
      class(energy), intent(in) :: en
      complex(rp), allocatable, intent(out) :: z_contour(:)
   end subroutine build_green_contour

   !> @brief Fill the canonical `green` arrays from the k-space engine.
   !> @param[inout] this      Reciprocal object (k-mesh, H(k) machinery).
   !> @param[inout] green_obj Green object whose arrays are populated in place.
   module subroutine fill_green(this, green_obj)
      class(reciprocal), intent(inout) :: this
      type(green), intent(inout) :: green_obj
   end subroutine fill_green

   !> @brief Fill mu_nm_stochastic from the k-space eigenpairs (exact moments).
   !> @param[inout] this Reciprocal object (H(k) machinery, k-mesh, velocities).
   !> @param[out]   mu   Allocated (nb,nb,cond_ll,cond_ll,ntype) moment array.
   !> @param[in]    a    Chebyshev window scale (recursion route’s window).
   !> @param[in]    b    Chebyshev window shift.
   module subroutine fill_moments(this, mu, a, b)
      class(reciprocal), intent(inout) :: this
      complex(rp), allocatable, intent(out) :: mu(:, :, :, :, :)
      real(rp), intent(in) :: a, b
   end subroutine fill_moments

   !> @brief Compute and write the Bloch spectral function A(k,E) along the path.
   !> @param[inout] this        Reciprocal object (H(k) machinery, k-path, grids).
   !> @param[in]    output_file Optional base output filename (default ’bsf.dat’).
   module subroutine calculate_bsf(this, output_file)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in), optional :: output_file
   end subroutine calculate_bsf

   !> @brief Write a dense k-space eigensystem for external Fermi-surface use.
   !> @details The post-processing caller promotes any reduced mesh to the
   !>          complete Monkhorst-Pack mesh, refreshes a stale spectrum against
   !>          the current Hamiltonian generation, gathers distributed k-points
   !>          on rank zero, and writes a stream-binary payload plus a text
   !>          metadata/index sidecar. Eigenvectors are stored as separate real
   !>          and imaginary float64 arrays in Fortran column-major order.
   !> @param[inout] this Reciprocal object containing the eigensystem.
   !> @param[in] scf_iteration Optional source SCF iteration recorded in metadata.
   module subroutine write_kspace_eigenpairs(this, scf_iteration)
      class(reciprocal), intent(inout) :: this
      integer, intent(in), optional :: scf_iteration
   end subroutine write_kspace_eigenpairs

   end interface

end module reciprocal_mod
