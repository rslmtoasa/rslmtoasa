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

   use control_mod
   use lattice_mod
   use hamiltonian_mod
   use charge_mod
   use energy_mod, only: energy
   use green_mod, only: green
   use sigma_provider_mod, only: sigma_provider, sigma_zero
   use spectrum_bounds_mod, only: bounds, compute_spectrum_bounds
   use precision_mod, only: rp
   use math_mod
   use string_mod, only: int2str, real2str, fmt, lower
   use logger_mod, only: g_logger
   use timer_mod, only: g_timer
   use symmetry_mod, only: symmetry
   use basis_mod, only: nb, norb
   use mpi_mod, only: rank, numprocs, ierr, get_mpi_range
#ifdef USE_MPI
   use mpi
#endif
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Module's main procedure
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
      !> Total number of k-points
      integer :: nk_total
      !> K-point coordinates in reciprocal lattice units
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
      !> Reciprocal solver mode: 'ham_only', 'generalized_overlap_proxy',
      !> or 'generalized_overlap_kanpur'.
      character(len=32) :: reciprocal_mode
      !> K-space Hamiltonian order:
      !>   'first'  -> H(k) = h(k)            (first-order, current behaviour)
      !>   'second' -> H(k) = E_nu + h(k) - [hoh](k) + L.S   (second-order ASA)
      !> The second-order path is the seam where the combined correction (CCOR)
      !> and the "proper" k-space second order (generalized_overlap_kanpur) plug in.
      character(len=16) :: kspace_ham_order
      !> Kanpur-alignment diagnostics toggle
      logical :: kanpur_diagnostics
      !> Optional H(Gamma) bounds diagnostics
      logical :: gamma_bounds_diagnostics
      !> Experimental finite real-space HALL diagonalization
      logical :: hall_diag_experimental

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
      !> Directional DOS x-component [n_energy_points]
      real(rp), dimension(:), allocatable :: dos_mx_tot
      !> Directional DOS y-component [n_energy_points]
      real(rp), dimension(:), allocatable :: dos_my_tot
      !> Directional DOS z-component [n_energy_points]
      real(rp), dimension(:), allocatable :: dos_mz_tot
      !> Projected directional DOS [n_sites, n_orb_types, n_spin, 3(x/y/z), n_energy_points]
      real(rp), dimension(:, :, :, :, :), allocatable :: projected_dos_moments
      !> DOS calculation method ('tetrahedron' or 'gaussian')
      character(len=20) :: dos_method
      !> Gaussian smearing parameter (in energy units)
      real(rp) :: gaussian_sigma
      !> Temperature for Fermi-Dirac distribution (in Kelvin)
      real(rp) :: temperature
      !> Fermi level for band moments integration
      real(rp) :: fermi_level
      !> Total number of valence electrons for Fermi level finding
      real(rp) :: total_electrons
      !> Flag to automatically find Fermi level from DOS
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

      ! K-path (band structure) control variables
      !> Automatic k-path generation enabled
      logical :: auto_kpath
      !> Number of k-points per segment in band structure
      integer :: nk_per_segment
      !> Override space group number (0 = auto-detect)
      integer :: override_space_group
      !> Custom k-path specification string
      character(len=200) :: custom_kpath_spec
      !> Use symmetry reduction for k-mesh
      logical :: use_symmetry_reduction
      !> Enforce strict symmetry consistency checks
      logical :: strict_symmetry_checks
      !> Dump symmetry k-point mapping diagnostics to file
      logical :: dump_symmetry_kmap
      !> Tetra symmetry backend: 'irreducible_native' for scalar DOS or 'full_expand_ref'
      character(len=32) :: tetra_symmetry_mode
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

      ! k-space Green's-function engine (milestone B2, reciprocal_green)
      !> Retarded broadening for the k-space Green contour z = E + i*green_eta (Ry)
      real(rp) :: green_eta
      !> Backend selector for fill_green: 'lehmann' (E, Sigma=0) or 'dyson' (D)
      character(len=16) :: green_backend

   contains
      procedure :: generate_mp_mesh
      procedure :: generate_reciprocal_vectors
      procedure :: build_kspace_hamiltonian
      procedure :: build_neighbor_vectors
      procedure :: calculate_structure_factors
      procedure :: fourier_transform_hamiltonian
      procedure :: fourier_transform_gbt
      procedure :: fourier_transform_gbt_array
      procedure :: fourier_transform_hamiltonian_second_order
      procedure :: fourier_transform_array
      procedure :: fourier_transform_overlap
      procedure :: set_basis_sizes
      procedure :: get_basis_type_from_size
      procedure :: check_multisite_hamiltonian_diagonal
      procedure :: build_kspace_overlap
      procedure :: diagonalize_hamiltonian
      procedure :: print_kanpur_mapping
      procedure :: check_overlap_properties
      procedure :: run_gamma_bounds_diagnostics
      procedure :: diagonalize_hall_experimental
      procedure :: calculate_band_structure
      procedure :: calculate_density_of_states
      procedure :: calculate_dos_tetrahedron
      procedure :: calculate_dos_tetrahedron_with_symmetry
      procedure :: expand_eigenvalues_to_full_mesh
      procedure :: calculate_dos_gaussian
      procedure :: calculate_dos_blochl
      procedure :: setup_dos_energy_grid
      procedure :: setup_tetrahedra
      procedure :: tetrahedron_dos_contribution
      procedure :: blochl_dos_contribution
      procedure :: get_kpoint_index
      procedure :: project_dos_orbitals
      procedure :: project_dos_orbitals_gaussian
      procedure :: project_dos_orbitals_tetrahedron
      procedure :: calculate_band_moments
   procedure :: print_total_and_spin_dos
      procedure :: calculate_band_energy_from_moments
      procedure :: calculate_adaptive_sigma
      procedure :: find_fermi_level_from_dos
      procedure :: calculate_ldm_from_projected_dos
      procedure :: integrate_dos_up_to_energy
      procedure :: calculate_gaussian_weight_single
      procedure :: write_dos_to_file
      procedure :: restore_to_default
      procedure :: build_from_file
      procedure :: set_kpoint_mesh
      procedure :: generate_reduced_kpoint_mesh
      procedure :: validate_symmetry_kmap
      procedure :: write_symmetry_kmap_dump
      procedure :: ensure_tetra_symmetry_backend
      procedure :: ensure_full_mesh_for_spinor_integrations
      procedure :: build_irreducible_tetrahedra
      ! k-space Green's-function engine (milestone B2)
      procedure :: fill_green
      procedure :: build_green_contour
      final     :: destructor
   end type reciprocal

   interface reciprocal
      procedure :: constructor
   end interface


   interface
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

      ! Default k-point mesh settings
   end subroutine restore_to_default

   !> @brief Read the &reciprocal namelist and install reciprocal-space options.
   !> @details Parses k-mesh, Fourier mode, band/DOS controls, symmetry options,
   !>          tetrahedron settings, and diagnostics from this%control%fname.
   !> @param[inout] this Reciprocal object whose input-facing options are populated.
   !> @note This is an input boundary and may raise fatal diagnostics for invalid options.
   module subroutine build_from_file(this)
      class(reciprocal), intent(inout) :: this

      ! Reading process variables
      integer :: iostatus, funit

      ! Include namelist declarations for reciprocal module
   end subroutine build_from_file

   !> @brief Set Monkhorst-Pack mesh dimensions.
   !> @param[inout] this Reciprocal object whose nk_mesh is updated.
   !> @param[in] nk1 Number of k-points along reciprocal axis 1.
   !> @param[in] nk2 Number of k-points along reciprocal axis 2.
   !> @param[in] nk3 Number of k-points along reciprocal axis 3.
   module subroutine set_kpoint_mesh(this, nk1, nk2, nk3)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: nk1, nk2, nk3

   end subroutine set_kpoint_mesh

   !> @brief Generate reciprocal lattice vectors from the real-space lattice.
   !> @details Computes the reciprocal basis and reciprocal-cell volume from
   !>          lattice%a and lattice%alat for later k-point and Fourier work.
   !> @param[inout] this Reciprocal object receiving reciprocal_vectors and volume.
   module subroutine generate_reciprocal_vectors(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      real(rp), dimension(3, 3) :: real_vectors
      real(rp) :: det
      integer :: i

      ! Get real-space lattice vectors from lattice%a
   end subroutine generate_reciprocal_vectors

   !> @brief Generate the configured Monkhorst-Pack k-point mesh.
   !> @details Builds full or symmetry-reduced fractional k-points, weights, and
   !>          optional MPI ownership metadata according to reciprocal settings.
   !> @param[inout] this Reciprocal object receiving k_points/k_weights state.
   module subroutine generate_mp_mesh(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ik, ix, iy, iz, nk_irred
      real(rp) :: kx, ky, kz
      real(rp), parameter :: tol = 1.0e-8_rp
      integer :: shift(3), nfull
      real(rp), allocatable :: kpoints_full(:,:)

   end subroutine generate_mp_mesh

   !> @brief Determine per-type basis sizes for reciprocal matrices.
   !> @details Maps each lattice atom type to sp/spd/spdf orbital counts used
   !>          when packing site-major k-space Hamiltonian and projection blocks.
   !> @param[inout] this Reciprocal object whose basis_size/max_orbs are updated.
   module subroutine set_basis_sizes(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ntype
      character(len=10) :: basis_type

   end subroutine set_basis_sizes

   !> @brief Build real-space neighbor vectors used by Fourier transforms.
   !> @details Converts lattice neighbor geometry into Cartesian and direct
   !>          per-type R-vector tables aligned with Hamiltonian neighbor slots.
   !> @param[inout] this Reciprocal object receiving ham_vec_type arrays.
   module subroutine build_neighbor_vectors(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ntype, ia, nr, nn_max_loc, kk
      real(rp) :: r2
      real(rp), dimension(3, this%lattice%kk) :: cralat

   end subroutine build_neighbor_vectors

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
   end subroutine calculate_structure_factors

   !> @brief Fourier transform first-order real-space Hamiltonian blocks.
   !> @details Forms H(k)=sum_R h(R) exp(i k.R), preserving the historical
   !>          first-order path without onsite e_nu or spin-orbit terms.
   !> @param[in] this Reciprocal object containing real-space Hamiltonian data.
   !> @param[in] k_vec k-point vector.
   !> @param[out] hk_result Packed k-space Hamiltonian matrix.
   module subroutine fourier_transform_hamiltonian(this, k_vec, hk_result)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: hk_result  ! (n_orb*n_sites, n_orb*n_sites)

      ! First-order k-space Hamiltonian: H(k) = h(k) = Sum_R ee(R) exp(i*k·R).
      ! NOTE: This deliberately reproduces the historical first-order behaviour.
      !       The on-site E_nu (enim) and spin-orbit (lsham) terms are NOT added
      !       here; they are only included in the second-order path
      !       (fourier_transform_hamiltonian_second_order). See kspace_ham_order.
   end subroutine fourier_transform_hamiltonian

   !> @brief Generalized-Bloch-theorem spin-spiral k-space Hamiltonian.
   !> @param[in] this Reciprocal object with the collinear reference Hamiltonian.
   !> @param[in] k_vec k-point vector (fractional coordinates).
   !> @param[out] hk_result Packed GBT k-space Hamiltonian matrix.
   module subroutine fourier_transform_gbt(this, k_vec, hk_result)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: hk_result
   end subroutine fourier_transform_gbt

   !> @brief GBT Fourier transform for one spinor neighbor/type block array.
   !> @param[in] this Reciprocal object with neighbor-vector tables.
   !> @param[in] array4d Spinor block array indexed by orbital, neighbor, and type.
   !> @param[in] k_vec k-point vector (fractional coordinates).
   !> @param[out] mk_result Packed GBT k-space matrix.
   module subroutine fourier_transform_gbt_array(this, array4d, k_vec, mk_result)
      class(reciprocal), intent(in) :: this
      complex(rp), dimension(:, :, :, :), intent(in) :: array4d
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: mk_result
   end subroutine fourier_transform_gbt_array

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

   end subroutine fourier_transform_array

   !> @brief Fourier transform the second-order ASA k-space Hamiltonian.
   !> @details Assembles onsite e_nu, first-order hopping, HOH, optional CCOR,
   !>          and spin-orbit contributions according to kspace_ham_order.
   !> @param[in] this Reciprocal object containing Hamiltonian and correction blocks.
   !> @param[in] k_vec k-point vector.
   !> @param[out] hk_result Packed second-order k-space Hamiltonian matrix.
   module subroutine fourier_transform_hamiltonian_second_order(this, k_vec, hk_result)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: hk_result  ! (n_orb*n_sites, n_orb*n_sites)
      ! Local variables
      integer :: isite, ntype_i, i_start, i_end
      integer :: n_orb, n_sites, ndim
      complex(rp), dimension(:, :), allocatable :: hk, eeok, hohk, hcck

   end subroutine fourier_transform_hamiltonian_second_order

   !> @brief Fourier transform overlap blocks into S(k).
   !> @details Builds the reciprocal-space overlap matrix used by generalized
   !>          eigenproblem modes and overlap diagnostics.
   !> @param[in] this Reciprocal object containing overlap and neighbor data.
   !> @param[in] k_vec k-point vector.
   !> @param[out] sk_result Packed k-space overlap matrix.
   module subroutine fourier_transform_overlap(this, k_vec, sk_result)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: sk_result
      integer :: isite, jsite, ntype_i, ineigh, ia, ja, nr, iorb
      integer :: i_start, i_end, j_start, j_end
      integer :: n_orb, n_sites
      complex(rp), dimension(:, :), allocatable :: structure_factors
      complex(rp), dimension(:, :), allocatable :: overlap_block

   end subroutine fourier_transform_overlap

   !> @brief Configure MPI ownership for a k-point set.
   !> @param[inout] this Reciprocal object receiving local/global k maps.
   !> @param[in] nk_global Number of global k-points.
   !> @param[in] enable_distribution If true, distribute k-points over MPI ranks.
   module subroutine setup_k_mesh_distribution(this, nk_global, enable_distribution)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: nk_global
      logical, intent(in) :: enable_distribution
      integer :: local_count
      integer :: ik

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
      ! Local variables
      integer :: ik, ik_global, nk, ntype
      character(len=200) :: debug_msg
      logical :: using_kpath, distribute_mesh
      logical :: use_second_order
      integer :: i, j

      ! Determine which k-point set to use
   end subroutine build_kspace_hamiltonian

   !> @brief Build S(k) for every active mesh or path k-point.
   !> @param[inout] this Reciprocal object receiving sk_overlap arrays.
   module subroutine build_kspace_overlap(this)
      class(reciprocal), intent(inout) :: this
      integer :: ik, ik_global, nk
      logical :: using_kpath

   end subroutine build_kspace_overlap

!> @brief Diagnose onsite spin-diagonal blocks for multisite Hamiltonians.
!> @param[in] this Reciprocal object containing hk_total or hk_bulk data.
module subroutine check_multisite_hamiltonian_diagonal(this)
   class(reciprocal), intent(in) :: this
   integer :: isite, iorb, ispin
   complex(rp) :: h_avg_site1_up, h_avg_site1_dn, h_avg_site2_up, h_avg_site2_dn
   integer :: n_sites, idx_up, idx_dn
   character(len=256) :: msg
   
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
   
end subroutine print_hamiltonian_structure

   !> @brief Return a basis label from the configured orbital count.
   !> @param[in] this Reciprocal object containing basis_size metadata.
   !> @param[in] ntype Atom type index.
   !> @return Basis label such as sp, spd, or spdf.
   module function get_basis_type_from_size(this, ntype) result(basis_type)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: ntype
      character(len=10) :: basis_type

   end function get_basis_type_from_size

   !> @brief Diagonalize the active k-space Hamiltonians.
   !> @details Solves standard or generalized eigenproblems depending on
   !>          reciprocal_mode and stores eigenvalues/eigenvectors for bands and DOS.
   !> @param[inout] this Reciprocal object receiving eigenvalue/eigenvector arrays.
   module subroutine diagonalize_hamiltonian(this)
      class(reciprocal), intent(inout) :: this
      
      ! Local variables
      integer :: nk, ik, nmat, lwork, info, mode_fail_count
      complex(rp), dimension(:, :), allocatable :: h_k_copy
      complex(rp), dimension(:, :), allocatable :: s_k_copy
      real(rp), dimension(:), allocatable :: eigenvals
      complex(rp), dimension(:), allocatable :: work_complex
      real(rp), dimension(:), allocatable :: rwork
      character(len=100) :: info_msg
      logical :: use_generalized

      ! Check prerequisites
   end subroutine diagonalize_hamiltonian

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
      integer :: i, j, n
      real(rp) :: max_herm
   end subroutine check_overlap_properties

   !> @brief Run H(Gamma) spectral-bound diagnostics.
   !> @param[inout] this Reciprocal object used to build and bound the Gamma matrix.
   module subroutine run_gamma_bounds_diagnostics(this)
      class(reciprocal), intent(inout) :: this
      complex(rp), allocatable :: h_gamma(:, :)
      real(rp) :: egmin, egmax
      type(bounds) :: bnd

   end subroutine run_gamma_bounds_diagnostics

   !> @brief Diagonalize the experimental real-space HALL matrix.
   !> @details Builds a finite local-cluster matrix from HALL blocks and prints
   !>          eigenvalue diagnostics for exploratory comparison only.
   !> @param[inout] this Reciprocal object providing Hamiltonian and lattice state.
   module subroutine diagonalize_hall_experimental(this)
      class(reciprocal), intent(inout) :: this
      integer :: nsites, n_orb, n, i, jsite, isite, ineigh, ia, ja, nr, info, lwork
      integer :: i_start, i_end, j_start, j_end
      complex(rp), allocatable :: hall_mat(:, :), work(:)
      real(rp), allocatable :: evals(:), rwork(:)

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
      ! Local variables
      character(len=256) :: filename
      integer :: unit, i, j, nmat
      character(len=100) :: fmt_str

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
      
      ! Local variables
      integer :: npts
      character(len=256) :: filename
      integer :: unit, i, j, nmat
      character(len=100) :: fmt_str

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
      integer :: shift(3)
      integer :: num_ir_kpoints
      logical :: do_shift, effective_time_reversal
      real(rp), allocatable :: kpoints_frac(:,:), weights(:)
      integer :: i
      integer, allocatable :: full_to_irred(:), irred_to_full(:)

   end subroutine generate_reduced_kpoint_mesh

   !> @brief Write a complex matrix and its k-point label to a text file.
   !> @param[in] matrix Matrix to dump.
   !> @param[in] filename Output file name.
   !> @param[in] k_point k-point associated with the matrix.
   module subroutine dump_complex_matrix(matrix, filename, k_point)
      complex(rp), dimension(:, :), intent(in) :: matrix
      character(len=*), intent(in) :: filename
      real(rp), dimension(3), intent(in) :: k_point
      ! Local variables
      integer :: unit, i, j, n_rows, n_cols
      character(len=100) :: fmt_str

   end subroutine dump_complex_matrix

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

      ! Local variables
      character(len=100) :: filename

   end subroutine calculate_density_of_states

   !> @brief Validate full-to-irreducible k-point symmetry maps.
   !> @param[inout] this Reciprocal object containing symmetry maps and weights.
   !> @param[in] context_tag Diagnostic label for the caller/context.
   module subroutine validate_symmetry_kmap(this, context_tag)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: context_tag
      integer :: nk_full_expected, nk_irred, i, idx
      integer, allocatable :: counts(:)
      real(rp) :: wsum

   end subroutine validate_symmetry_kmap

   !> @brief Write symmetry k-point maps for diagnostics.
   !> @param[inout] this Reciprocal object containing symmetry maps.
   !> @param[in] filename Output file name.
   module subroutine write_symmetry_kmap_dump(this, filename)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer :: u, i, nk_full_expected

   end subroutine write_symmetry_kmap_dump

   !> @brief Ensure the configured tetrahedron symmetry backend is usable.
   !> @details Expands or preserves the k-mesh as needed for scalar tetrahedron
   !>          DOS versus spinor/projection integrations.
   !> @param[inout] this Reciprocal object whose mesh/backend state may be updated.
   module subroutine ensure_tetra_symmetry_backend(this)
      class(reciprocal), intent(inout) :: this
      integer :: nk_full_expected
      logical :: need_full

   end subroutine ensure_tetra_symmetry_backend

   !> @brief Ensure spinor integrations operate on a full k-mesh.
   !> @param[inout] this Reciprocal object whose eigenvalues/eigenvectors may be expanded.
   !> @param[in] context_tag Diagnostic label for the caller/context.
   module subroutine ensure_full_mesh_for_spinor_integrations(this, context_tag)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: context_tag
      integer :: nk_full_expected

   end subroutine ensure_full_mesh_for_spinor_integrations

   !> @brief Build irreducible tetrahedra and multiplicities from full-mesh cuts.
   !> @param[inout] this Reciprocal object containing symmetry maps and mesh dimensions.
   !> @param[out] tet_ir Irreducible tetrahedron corner indices.
   !> @param[out] tet_mult Multiplicity of each irreducible tetrahedron.
   !> @param[out] n_tet_ir Number of irreducible tetrahedra.
   module subroutine build_irreducible_tetrahedra(this, tet_ir, tet_mult, n_tet_ir)
      class(reciprocal), intent(inout) :: this
      integer, allocatable, intent(out) :: tet_ir(:, :)
      integer, allocatable, intent(out) :: tet_mult(:)
      integer, intent(out) :: n_tet_ir
      integer :: nk1, nk2, nk3, i, j, k, it, ic, idx
      integer :: n_tet_per_cube, tet_full_count, key_pos
      integer :: nk_full_expected
      integer, dimension(3, 4, 6) :: tetra_cut
      integer, dimension(4) :: key_sorted
      integer, allocatable :: keys(:, :), mult_tmp(:)
      logical :: found

   end subroutine build_irreducible_tetrahedra

   !> @brief Select tetrahedron corner offsets for the current mesh geometry.
   !> @param[in] this Reciprocal object providing reciprocal lattice geometry.
   !> @param[out] tetra_cut Six tetrahedra per mesh cell, stored as corner offsets.
   module subroutine get_tetra_cut_offsets(this, tetra_cut)
      class(reciprocal), intent(in) :: this
      integer, intent(out) :: tetra_cut(3, 4, 6)

      integer :: base_cut(3, 4, 6), test_cut(3, 4, 6)
      integer :: lx, ly, it, ic, jc, best_lx, best_ly
      real(rp) :: qvec(3, 3), p(3, 4), dp(3)
      real(rp) :: edge2, max_edge2, best_max_edge2

   end subroutine get_tetra_cut_offsets

   !> @brief Build the DOS energy grid from the configured range and size.
   !> @param[inout] this Reciprocal object receiving dos_energy_grid.
   module subroutine setup_dos_energy_grid(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i
      real(rp) :: energy_min, energy_max, delta_energy

   end subroutine setup_dos_energy_grid

   !> @brief Calculate total DOS using the tetrahedron method.
   !> @details Accumulates DOS and integrated NOS over full or irreducible
   !>          tetrahedra, respecting MPI partitioning where active.
   !> @param[inout] this Reciprocal object receiving total_dos and total_nos.
   module subroutine calculate_dos_tetrahedron(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i_tet, i_corner, i_band, nbands
      real(rp), dimension(4) :: e_corners
      integer, allocatable :: tet_ir(:, :), tet_mult(:)
      integer :: n_tet_ir
      real(rp) :: tet_weight, dos_integral, nos_integral
      real(rp) :: fermi_count, fermi_error
      real(rp), allocatable :: local_dos(:), local_nos(:)
      integer :: nk_full_expected, tet_start, tet_end, tet_count

   end subroutine calculate_dos_tetrahedron

   !> @brief Evaluate linear-tetrahedron DOS contribution at one energy.
   !> @param[in] this Reciprocal object providing method context.
   !> @param[in] energy Energy where the contribution is evaluated.
   !> @param[in] e_sorted Sorted tetrahedron corner energies.
   !> @return DOS contribution for one tetrahedron.
   module function tetrahedron_dos_contribution(this, energy, e_sorted) result(dos)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: energy
      real(rp), dimension(4), intent(in) :: e_sorted
      real(rp) :: dos

      ! Local variables
      real(rp) :: e1, e2, e3, e4, vol_factor
      real(rp), parameter :: eps = 1.0e-12_rp

      ! Tetrahedron volume factor (adjusted for correct normalization)
   end function tetrahedron_dos_contribution

   !> @brief Evaluate Bloechl-corrected tetrahedron DOS contribution.
   !> @param[in] this Reciprocal object providing method context.
   !> @param[in] energy Energy where the contribution is evaluated.
   !> @param[in] e_sorted Sorted tetrahedron corner energies.
   !> @return Corrected DOS contribution for one tetrahedron.
   module function blochl_dos_contribution(this, energy, e_sorted) result(dos)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: energy
      real(rp), dimension(4), intent(in) :: e_sorted
      real(rp) :: dos

      ! Local variables
      real(rp) :: e1, e2, e3, e4, C
      real(rp), parameter :: TOL = 1.0e-10_rp

   end function blochl_dos_contribution

   !> @brief Sort a real array while tracking original indices.
   !> @param[in] arr Values to sort.
   !> @param[out] sorted Sorted values.
   !> @param[out] indices Original indices corresponding to sorted values.
   module subroutine sort_real_array(arr, sorted, indices)
      real(rp), dimension(:), intent(in) :: arr
      real(rp), dimension(:), intent(out) :: sorted
      integer, dimension(:), intent(out) :: indices

      ! Local variables
      integer :: i, j, n, temp_idx
      real(rp) :: temp_val

   end subroutine sort_real_array

   !> @brief Sort four real values in ascending order.
   !> @param[in] arr_in Four input values.
   !> @param[out] arr_out Sorted output values.
   module subroutine sort4(arr_in, arr_out)
      real(rp), dimension(4), intent(in) :: arr_in
      real(rp), dimension(4), intent(out) :: arr_out
      real(rp) :: a1, a2, a3, a4, tmp

   end subroutine sort4

   !> @brief Add one tetrahedron's integrated-state contribution to a grid.
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

      integer :: i, i0(4), m
      real(rp) :: ecorn(4), de, e1, e2, e3, e4
      real(rp) :: c0, c1, c2, c3, x, x0
      real(rp), parameter :: tol = 1.0e-12_rp

   end subroutine tetra_add_nos

   !> @brief Add one tetrahedron's DOS contribution to a grid.
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

      integer :: i, i0(4), m
      real(rp) :: ecorn(4), de, e1, e2, e3, e4
      real(rp) :: c1, c2, c3, x
      real(rp), parameter :: tol = 1.0e-12_rp

   end subroutine tetra_add_dos

!> @brief Calculate total DOS with Gaussian broadening.
!> @details Smears all eigenvalues over the DOS energy grid using the configured
!>          or adaptive sigma and k-point weights.
!> @param[inout] this Reciprocal object receiving total_dos.
module subroutine calculate_dos_gaussian(this)
   class(reciprocal), intent(inout) :: this

   integer :: i_energy, i_k, i_band, i_k_global
   real(rp) :: energy, weight, gaussian_factor
   real(rp) :: sigma_squared, sigma_use
   real(rp) :: local_sum, dos_integral, norm_factor, kweight_sum
   integer :: nbands

   ! Debug: Check eigenvalue range
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

      ! Local variables
      integer :: nk1, nk2, nk3, n_tet_per_cube, i, j, k, tet_idx
      integer :: it, ic, idx
      integer, dimension(3, 4, 6) :: tetra_cut

   end subroutine setup_tetrahedra

   !> @brief Expand irreducible-k eigenvalues/eigenvectors to the full mesh.
   !> @param[inout] this Reciprocal object whose eigen arrays are expanded.
   module subroutine expand_eigenvalues_to_full_mesh(this)
      class(reciprocal), intent(inout) :: this
      integer :: ik_full, nk_full, nk_irred, nbands, ik_irred
      real(rp), dimension(:, :), allocatable :: eigenvalues_full
      complex(rp), dimension(:, :, :), allocatable :: eigenvectors_full
      integer :: nrow_evec, nband_evec, nk_evec
      
   end subroutine expand_eigenvalues_to_full_mesh

   !> @brief Calculate tetrahedron DOS using symmetry-aware mesh handling.
   !> @param[inout] this Reciprocal object receiving total DOS on the chosen backend.
   module subroutine calculate_dos_tetrahedron_with_symmetry(this)
      class(reciprocal), intent(inout) :: this
      
   end subroutine calculate_dos_tetrahedron_with_symmetry

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

      ! Apply periodic boundary conditions
   end function get_kpoint_index

   !> @brief Estimate a Gaussian DOS broadening from k-point density.
   !> @param[in] this Reciprocal object containing mesh and reciprocal volume.
   !> @return Adaptive Gaussian sigma.
   module function calculate_adaptive_sigma(this) result(sigma)
      class(reciprocal), intent(in) :: this
      real(rp) :: sigma
      real(rp) :: bz_volume, k_density, typical_spacing
      integer :: nk_total
      
      ! Calculate Brillouin zone k-point density
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
   integer :: ik, ik_global, ib, ie, iorb, i, isite
   integer :: n_orb_per_spin, orb_start, site_orb_start, n_orb_site
   integer :: lstart(4), lend(4)
   real(rp) :: weight, orbital_char, energy
   real(rp) :: gaussian_weight, sigma_squared, sigma_use
   complex(rp) :: psi_element
   real(rp) :: dos_integral, norm_factor
   integer :: nbands
   ! Additional locals for projected DOS integration/normalization
   real(rp) :: proj_integral, e_low, e_high
   integer :: iei
   ! Per-site orbital offsets for mixed atom types
   integer, dimension(:), allocatable :: site_orb_offset
   real(rp) :: mx_char, my_char, mz_char, local_char
   real(rp) :: axis(3)

end subroutine project_dos_orbitals_gaussian

   !> @brief Return the local spin axis for a projected-DOS site.
   !> @param[in] this Reciprocal object containing lattice magnetic moments.
   !> @param[in] isite Site index in the reciprocal packing.
   !> @param[out] axis Normalized local spin axis.
   module subroutine get_site_spin_axis(this, isite, axis)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: isite
      real(rp), intent(out) :: axis(3)
      integer :: atom_idx
      real(rp) :: axis_norm

   end subroutine get_site_spin_axis

   !> @brief Synchronize lattice local-density-matrix data across MPI ranks.
   !> @param[inout] lattice_obj Lattice object whose ldm-like storage is synchronized.
   module subroutine sync_lattice_ldm(lattice_obj)
      use lattice_mod
      type(lattice), intent(inout) :: lattice_obj
      real(rp), allocatable :: ldm_comm(:, :, :, :)
      integer :: max_flat_ldm, local_flat
      integer :: na_glob, plusbulk, lcount_ldm, l, ispin

   end subroutine sync_lattice_ldm

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

      integer :: i, idx_up, idx_dn
      real(rp) :: up_w, dn_w
      complex(rp) :: u, d, ud

   end subroutine compute_spinor_block_projection

   !> @brief Calculate projected DOS with tetrahedron integration.
   !> @details Averages orbital/spinor projection weights over tetrahedron
   !>          corners and accumulates scalar plus directional projected DOS.
   !> @param[inout] this Reciprocal object receiving projected DOS arrays.
   module subroutine project_dos_orbitals_tetrahedron(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i_energy, i_tet, i_corner, i_band, iorb, isite
      integer :: i_start, i_end
      integer :: n_orb_per_spin, orb_start, site_orb_start, ik, n_orb_site
      integer :: ie, nbands
      integer :: tet_start, tet_end, tet_count
      real(rp) :: energy, dos_contrib, orbital_char_avg, orbital_char, de
      real(rp) :: total_dos_integral, proj_dos_integral, norm_factor
      real(rp), dimension(4) :: e_corners, sorted_e, orbital_chars
      real(rp), dimension(4) :: mx_chars, my_chars, mz_chars
      integer :: lstart(4), lend(4)
      real(rp) :: e1, e2, e3, e4, x, C
      real(rp) :: tet_weight
      real(rp), parameter :: TOL = 1.0e-10_rp
      integer :: n_tet_ir
      ! Per-site orbital offsets for mixed atom types
      integer, dimension(:), allocatable :: site_orb_offset
      real(rp) :: mx_char_avg, my_char_avg, mz_char_avg, local_char
      real(rp) :: axis(3)
      real(rp), allocatable :: site_axes(:, :)
      real(rp), allocatable :: dos_line(:)
      real(rp), allocatable :: local_projected_dos(:, :, :, :)
      real(rp), allocatable :: local_projected_dos_moments(:, :, :, :, :)
      real(rp), allocatable :: local_dos_mx(:), local_dos_my(:), local_dos_mz(:)

   end subroutine project_dos_orbitals_tetrahedron

   !> @brief Integrate tabulated data with the trapezoidal rule.
   !> @param[in] x Grid coordinates.
   !> @param[in] y Values on the grid.
   !> @return Trapezoidal integral over the grid.
   module function trapezoidal_integral(x, y) result(integral)
      real(rp), dimension(:), intent(in) :: x, y
      real(rp) :: integral

      ! Local variables
      integer :: n, i
      real(rp) :: dx

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
      integer :: i
      real(rp) :: dx

   end function interpolate_grid_value

!> @brief Calculate zeroth, first, and second DOS band moments.
!> @details Integrates projected DOS channels with Fermi weighting to produce
!>          occupation, first-energy, and second-energy moments.
!> @param[inout] this Reciprocal object receiving band_moments.
module subroutine calculate_band_moments(this)
   class(reciprocal), intent(inout) :: this

   integer :: isite, iorb, ispin, ie, n_energy
   real(rp) :: energy, dos_value, fermi_weight
   real(rp) :: m0, m1, m2
   real(rp), dimension(:), allocatable :: integrand, fermi_dist
   real(rp) :: kT, fermi_arg
   real(rp) :: total_occupation, expected_electrons, total_occupation_alt
      real(rp), allocatable :: energy_grid_ry(:)
      real(rp), parameter :: eV_to_Ry = 0.073498618_rp

end subroutine calculate_band_moments

   !> @brief Print total and spin-resolved DOS occupation diagnostics.
   !> @param[in] this Reciprocal object containing DOS and projected DOS arrays.
   module subroutine print_total_and_spin_dos(this)
      class(reciprocal), intent(in) :: this

      integer :: ie, n_energy, isite, iorb, ispin
      real(rp), allocatable :: fermi_dist(:)
      real(rp), allocatable :: energy_grid(:)
      real(rp), allocatable :: integrand_up(:), integrand_dn(:)
      real(rp) :: kT, fermi_arg
      real(rp) :: electrons_up, electrons_dn, electrons_total
      real(rp), parameter :: kB_Ry_per_K = 6.3336814e-6_rp

   end subroutine print_total_and_spin_dos

!> @brief Find the Fermi level matching a target electron count.
!> @param[in] this Reciprocal object containing DOS/NOS information.
!> @param[in] total_electrons Target number of electrons.
!> @return Fermi level in the reciprocal DOS energy units.
module function find_fermi_level_from_dos(this, total_electrons) result(fermi_level)
   class(reciprocal), intent(in) :: this
   real(rp), intent(in) :: total_electrons
   real(rp) :: fermi_level

   integer :: ie, max_iter
   real(rp) :: energy, kT
   real(rp) :: e_min, e_max, e_mid, electrons_at_e
   real(rp) :: q1, q2, e1, e2, step
   real(rp), parameter :: eV_to_Ry = 0.073498618_rp
   real(rp), parameter :: kB_Ry_per_K = 6.3336814e-6_rp

end function find_fermi_level_from_dos

!> @brief Integrate DOS occupation up to an energy with Fermi weighting.
!> @param[in] this Reciprocal object containing DOS grid data.
!> @param[in] energy Trial Fermi energy.
!> @param[in] kT Thermal energy used in Fermi weighting.
!> @return Integrated electron count.
module function integrate_dos_up_to_energy(this, energy, kT) result(integral)
   class(reciprocal), intent(in) :: this
   real(rp), intent(in) :: energy, kT
   real(rp) :: integral

   integer :: ie
   real(rp) :: e, fermi_weight, fermi_weight_next, delta_e
   real(rp), parameter :: eV_to_Ry = 0.073498618_rp

end function integrate_dos_up_to_energy

   !> @brief Write total and projected DOS data to output files.
   !> @param[in] this Reciprocal object containing DOS arrays.
   !> @param[in] filename Primary total-DOS output file.
   module subroutine write_dos_to_file(this, filename)
      class(reciprocal), intent(in) :: this
      character(len=*), intent(in) :: filename

      ! Local variables
      integer :: unit, i_energy, isite, iorb, ispin
      character(len=256) :: proj_filename
      real(rp), parameter :: eV_to_Ry = 0.073498618_rp

   end subroutine write_dos_to_file

   !> @brief Calculate band energy from integrated DOS moments.
   !> @param[in] this Reciprocal object containing band_moments.
   !> @return Band energy accumulated over projected channels.
   module function calculate_band_energy_from_moments(this) result(eband)
      class(reciprocal), intent(in) :: this
      real(rp) :: eband
      integer :: isite, iorb, ispin
      
   end function calculate_band_energy_from_moments

   !> @brief Evaluate one normalized Gaussian smearing weight.
   !> @param[in] this Reciprocal object containing gaussian_sigma/adaptive settings.
   !> @param[in] grid_energy Energy-grid point.
   !> @param[in] eigenvalue Eigenvalue to smear.
   !> @return Gaussian weight at grid_energy.
   module function calculate_gaussian_weight_single(this, grid_energy, eigenvalue) result(weight)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: grid_energy, eigenvalue
      real(rp) :: weight

      ! Gaussian smearing: weight = exp(-(E_grid - E_eigen)²/(2σ²)) / (σ√(2π))
      real(rp) :: prefactor, exponent, delta_e

   end function calculate_gaussian_weight_single

   !> @brief Sort four tetrahedron eigenvalues and return their permutation.
   !> @param[in] e_in Input corner eigenvalues.
   !> @param[out] e_sorted Sorted corner eigenvalues.
   !> @param[out] sort_idx Original indices corresponding to sorted values.
   module pure subroutine sort_eigenvalues(e_in, e_sorted, sort_idx)
      real(rp), dimension(4), intent(in) :: e_in
      real(rp), dimension(4), intent(out) :: e_sorted
      integer, dimension(4), intent(out) :: sort_idx
      
      ! Local variables
      integer :: i, j, min_idx
      real(rp) :: temp_e
      integer :: temp_idx
      
      ! Initialize
   end subroutine sort_eigenvalues

   !> @brief Reconstruct lattice local density matrices from projected DOS moments.
   !> @param[inout] this Reciprocal object containing projected DOS/moments.
   !> @param[inout] lattice_obj Lattice object receiving local density matrices.
   module subroutine calculate_ldm_from_projected_dos(this, lattice_obj)
      use lattice_mod
      class(reciprocal), intent(inout) :: this
      type(lattice), intent(inout) :: lattice_obj
      
      ! Local variables
      integer :: isite
      
      ! Orbital indexing: 
      ! iorb = 1 (s, 1 orbital), iorb = 2 (p, 3 orbitals), iorb = 3 (d, 5 orbitals)
      
   end subroutine calculate_ldm_from_projected_dos

   !> @brief Reconstruct lattice local density matrices directly from eigenvectors.
   !> @details Accumulates occupied spinor density matrices over k-points and
   !>          bands, then synchronizes the lattice result across MPI ranks.
   !> @param[inout] this Reciprocal object containing eigenvectors and occupations.
   !> @param[inout] lattice_obj Lattice object receiving local density matrices.
   module subroutine calculate_ldm_from_eigenvectors(this, lattice_obj)
      use lattice_mod
      class(reciprocal), intent(inout) :: this
      type(lattice), intent(inout) :: lattice_obj
      
      ! Local variables
      integer :: ik, ik_global, ib, isite, alpha, beta, iorb, norb_l
      integer :: ispin, basis_offset, site_offset, n_orb_site, n_orb_per_spin
      real(rp) :: kweight, fermi_occ
      complex(rp), allocatable :: DM(:,:)
      complex(rp) :: psi_a, psi_b
      integer :: nbands, nsites, i
      real(rp) :: trace_tot, trace_up, trace_dn, max_imag
      integer :: idx_a, idx_b

   end subroutine calculate_ldm_from_eigenvectors

   ! k-space Green's-function engine (milestone B2, submodule reciprocal_green)

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
   !> @param[in]    sigma     Self-energy provider (sigma_zero for backend E).
   module subroutine fill_green(this, green_obj, sigma)
      class(reciprocal), intent(inout) :: this
      type(green), intent(inout) :: green_obj
      class(sigma_provider), intent(in) :: sigma
   end subroutine fill_green

   end interface

end module reciprocal_mod
