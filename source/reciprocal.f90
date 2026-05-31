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
   use spectrum_bounds_mod, only: bounds, compute_spectrum_bounds
   use precision_mod, only: rp
   use math_mod
   use string_mod, only: int2str, real2str, fmt, lower
   use logger_mod, only: g_logger
   use timer_mod, only: g_timer
   use symmetry_mod, only: symmetry
   use basis_mod, only: nb, norb
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
      !> Projected DOS [n_sites, n_orb_types, n_spin, n_energy_points]
      real(rp), dimension(:, :, :, :), allocatable :: projected_dos
      !> Band moments [m0, m1, m2] for each projection
      real(rp), dimension(:, :, :, :), allocatable :: band_moments
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
      !> Number of spin components (1 for non-spin-polarized, 2 for spin-polarized)
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

   contains
      procedure :: generate_mp_mesh
      procedure :: generate_reciprocal_vectors
      procedure :: build_kspace_hamiltonian
      procedure :: build_neighbor_vectors
      procedure :: calculate_structure_factors
      procedure :: fourier_transform_hamiltonian
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
      final     :: destructor
   end type reciprocal

   interface reciprocal
      procedure :: constructor
   end interface

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] hamiltonian_obj Variable that holds hamiltonian_mod properties
   !> @return type(reciprocal)
   !---------------------------------------------------------------------------
   function constructor(hamiltonian_obj) result(obj)
      type(reciprocal) :: obj
      type(hamiltonian), target, intent(in) :: hamiltonian_obj

      obj%hamiltonian => hamiltonian_obj
      obj%lattice => hamiltonian_obj%lattice
      obj%control => hamiltonian_obj%lattice%control

      call obj%restore_to_default()
      call obj%build_from_file()  ! Read parameters from input.nml
      call obj%generate_reciprocal_vectors()
      call obj%set_basis_sizes()
      call obj%symmetry_analysis%initialize(obj%lattice)
      
      ! Auto-calculate total electrons from valence if not set in input
      ! Similar to bands.f90: this%qqv = real(sum(this%symbolic_atom(1:this%lattice%nbulk_bulk)%element%valence))
      if (obj%total_electrons <= 1.0e-3_rp) then
         obj%total_electrons = real(sum(obj%lattice%symbolic_atoms(1:obj%lattice%nbulk_bulk)%element%valence), rp)
         call g_logger%info('reciprocal%constructor: Auto-calculated total_electrons = ' // trim(real2str(obj%total_electrons)) // ' from valence', __FILE__, __LINE__)
      end if
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(reciprocal) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%k_points)) call g_safe_alloc%deallocate('reciprocal.k_points', this%k_points)
      if (allocated(this%k_weights)) call g_safe_alloc%deallocate('reciprocal.k_weights', this%k_weights)
      if (allocated(this%hk_bulk)) call g_safe_alloc%deallocate('reciprocal.hk_bulk', this%hk_bulk)
      if (allocated(this%hk_so)) call g_safe_alloc%deallocate('reciprocal.hk_so', this%hk_so)
      if (allocated(this%hk_total)) call g_safe_alloc%deallocate('reciprocal.hk_total', this%hk_total)
      if (allocated(this%sk_overlap)) call g_safe_alloc%deallocate('reciprocal.sk_overlap', this%sk_overlap)
      if (allocated(this%basis_size)) call g_safe_alloc%deallocate('reciprocal.basis_size', this%basis_size)
      if (allocated(this%k_path)) call g_safe_alloc%deallocate('reciprocal.k_path', this%k_path)
      if (allocated(this%k_labels)) call g_safe_alloc%deallocate('reciprocal.k_labels', this%k_labels)
      if (allocated(this%k_distances)) call g_safe_alloc%deallocate('reciprocal.k_distances', this%k_distances)
      if (allocated(this%eigenvalues)) call g_safe_alloc%deallocate('reciprocal.eigenvalues', this%eigenvalues)
      if (allocated(this%eigenvalues_path)) call g_safe_alloc%deallocate('reciprocal.eigenvalues_path', this%eigenvalues_path)
      if (allocated(this%eigenvectors)) call g_safe_alloc%deallocate('reciprocal.eigenvectors', this%eigenvectors)
      if (allocated(this%eigenvectors_path)) call g_safe_alloc%deallocate('reciprocal.eigenvectors_path', this%eigenvectors_path)
      ! DOS arrays
      if (allocated(this%dos_energy_grid)) call g_safe_alloc%deallocate('reciprocal.dos_energy_grid', this%dos_energy_grid)
      if (allocated(this%total_dos)) call g_safe_alloc%deallocate('reciprocal.total_dos', this%total_dos)
      if (allocated(this%projected_dos)) call g_safe_alloc%deallocate('reciprocal.projected_dos', this%projected_dos)
      if (allocated(this%band_moments)) call g_safe_alloc%deallocate('reciprocal.band_moments', this%band_moments)
      if (allocated(this%tetrahedra)) call g_safe_alloc%deallocate('reciprocal.tetrahedra', this%tetrahedra)
      if (allocated(this%tetrahedron_volumes)) call g_safe_alloc%deallocate('reciprocal.tetrahedron_volumes', this%tetrahedron_volumes)
      if (allocated(this%ham_vec_type)) call g_safe_alloc%deallocate('reciprocal.ham_vec_type', this%ham_vec_type)
      if (allocated(this%ham_vec_type_direct)) call g_safe_alloc%deallocate('reciprocal.ham_vec_type_direct', this%ham_vec_type_direct)
#else
      if (allocated(this%k_points)) deallocate (this%k_points)
      if (allocated(this%k_weights)) deallocate (this%k_weights)
      if (allocated(this%hk_bulk)) deallocate (this%hk_bulk)
      if (allocated(this%hk_so)) deallocate (this%hk_so)
      if (allocated(this%hk_total)) deallocate (this%hk_total)
      if (allocated(this%sk_overlap)) deallocate (this%sk_overlap)
      if (allocated(this%basis_size)) deallocate (this%basis_size)
      if (allocated(this%k_path)) deallocate (this%k_path)
      if (allocated(this%k_labels)) deallocate (this%k_labels)
      if (allocated(this%k_distances)) deallocate (this%k_distances)
   if (allocated(this%eigenvalues)) deallocate(this%eigenvalues)
   if (allocated(this%eigenvectors)) deallocate(this%eigenvectors)
      ! DOS arrays
      if (allocated(this%dos_energy_grid)) deallocate (this%dos_energy_grid)
      if (allocated(this%total_dos)) deallocate (this%total_dos)
      if (allocated(this%projected_dos)) deallocate (this%projected_dos)
      if (allocated(this%band_moments)) deallocate (this%band_moments)
      if (allocated(this%tetrahedra)) deallocate (this%tetrahedra)
      if (allocated(this%tetrahedron_volumes)) deallocate (this%tetrahedron_volumes)
      if (allocated(this%ham_vec_type)) deallocate (this%ham_vec_type)
      if (allocated(this%ham_vec_type_direct)) deallocate (this%ham_vec_type_direct)
      if (allocated(this%full_to_irred_k)) deallocate(this%full_to_irred_k)
      if (allocated(this%irred_to_full_k)) deallocate(this%irred_to_full_k)
#endif
   end subroutine destructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Restore variables to default values
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      class(reciprocal), intent(inout) :: this

      ! Default k-point mesh settings
      this%nk_mesh = [8, 8, 8]  ! Default 8x8x8 mesh
      this%nk_total = 0
      this%use_time_reversal = .true.
      this%strict_symmetry_checks = .false.
      this%k_offset = [0.0_rp, 0.0_rp, 0.0_rp]  ! No shift by default
      this%include_so = .false.
      this%max_orbs = nb

   ! By default suppress internal verbose prints (can be enabled by user)
   this%suppress_internal_logs = .true.

      ! Initialize reciprocal lattice to zero
      this%reciprocal_vectors = 0.0_rp
      this%reciprocal_volume = 0.0_rp

      ! Default DOS settings
      ! All energy values in Ry (consistent with Hamiltonian)
      this%n_energy_points = 1000
      this%dos_energy_range = [-1.0_rp, 1.0_rp]  ! Default energy range in Ry
      this%dos_method = 'tetrahedron'  ! Default to tetrahedron method
      this%gaussian_sigma = 0.01_rp  ! Default Gaussian smearing in Ry
      this%temperature = 300.0_rp  ! Default temperature in Kelvin
      this%fermi_level = 0.0_rp  ! Default Fermi level
      this%total_electrons = 0.0_rp  ! 0 = auto-calculate from valence in constructor
      this%auto_find_fermi = .true.  ! Auto-find Fermi level from DOS (recommended default)
      this%reciprocal_mode = 'ham_only'
      this%kanpur_diagnostics = .true.
      this%gamma_bounds_diagnostics = .false.
      this%hall_diag_experimental = .false.
      this%n_sites = 0
      this%n_orb_types = 4  ! s, p, d, f
      this%n_spin_components = 1  ! Default to non-spin-polarized
      this%n_tetrahedra = 0

      ! Default k-path settings
      this%auto_kpath = .true.  ! Use automatic k-path generation by default
      this%nk_per_segment = 40  ! Default 40 k-points per segment
      this%override_space_group = 0  ! 0 = auto-detect
      this%custom_kpath_spec = ''  ! Empty = use automatic
      this%use_symmetry_reduction = .true.  ! Use symmetry reduction by default
   end subroutine restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read parameters from input file (namelist-based input)
   !---------------------------------------------------------------------------
   subroutine build_from_file(this)
      class(reciprocal), intent(inout) :: this

      ! Reading process variables
      integer :: iostatus, funit

      ! Include namelist declarations for reciprocal module
      include 'include_codes/namelists/reciprocal.f90'
      include 'include_codes/namelists/kpath.f90'

      ! Save previous values (from defaults or previous read)
      nk1 = this%nk_mesh(1)
      nk2 = this%nk_mesh(2)
      nk3 = this%nk_mesh(3)
      k_offset_x = this%k_offset(1)
      k_offset_y = this%k_offset(2)
      k_offset_z = this%k_offset(3)
      use_symmetry_reduction = this%use_symmetry_reduction
      use_time_reversal = this%use_time_reversal
      use_shift = .false.  ! Derived from k_offset
      n_energy_points = this%n_energy_points
      dos_energy_min = this%dos_energy_range(1)
      dos_energy_max = this%dos_energy_range(2)
      gaussian_sigma = this%gaussian_sigma
      temperature = this%temperature
      dos_method = this%dos_method
      auto_find_fermi = this%auto_find_fermi
      suppress_internal_logs = this%suppress_internal_logs
      reciprocal_mode = this%reciprocal_mode
      kanpur_diagnostics = this%kanpur_diagnostics
      gamma_bounds_diagnostics = this%gamma_bounds_diagnostics
      hall_diag_experimental = this%hall_diag_experimental
      
      ! K-path settings
      auto_kpath = this%auto_kpath
      nk_per_segment = this%nk_per_segment
      override_space_group = this%override_space_group
      custom_kpath_spec = this%custom_kpath_spec

      ! Read reciprocal namelist
      open (newunit=funit, file=this%control%fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(this%control%fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=reciprocal, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         ! Namelist not found or error - use defaults
         call g_logger%info('reciprocal namelist not found in input file, using defaults', __FILE__, __LINE__)
      end if
      
      ! Read kpath namelist
      rewind(funit)
      read (funit, nml=kpath, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         ! Namelist not found or error - use defaults
         call g_logger%info('kpath namelist not found in input file, using defaults', __FILE__, __LINE__)
      end if
      close (funit)

      ! Assign values back to type members
      this%nk_mesh = [nk1, nk2, nk3]
      this%k_offset = [k_offset_x, k_offset_y, k_offset_z]
      this%use_symmetry_reduction = use_symmetry_reduction
      this%use_time_reversal = use_time_reversal
      this%n_energy_points = n_energy_points
      this%dos_energy_range = [dos_energy_min, dos_energy_max]
      this%gaussian_sigma = gaussian_sigma
      this%temperature = temperature
      this%dos_method = dos_method
      this%auto_find_fermi = auto_find_fermi
      this%suppress_internal_logs = suppress_internal_logs
      this%reciprocal_mode = lower(trim(reciprocal_mode))
      this%kanpur_diagnostics = kanpur_diagnostics
      this%gamma_bounds_diagnostics = gamma_bounds_diagnostics
      this%hall_diag_experimental = hall_diag_experimental
      if (this%reciprocal_mode == 'generalized_overlap') then
         this%reciprocal_mode = 'generalized_overlap_proxy'
         call g_logger%warning("reciprocal_mode='generalized_overlap' is deprecated alias; using 'generalized_overlap_proxy'.", __FILE__, __LINE__)
      end if
      if (this%reciprocal_mode /= 'ham_only' .and. &
          this%reciprocal_mode /= 'generalized_overlap_proxy' .and. &
          this%reciprocal_mode /= 'generalized_overlap_kanpur') then
         call g_logger%warning("reciprocal_mode must be 'ham_only', 'generalized_overlap_proxy', or 'generalized_overlap_kanpur'. Falling back to ham_only.", __FILE__, __LINE__)
         this%reciprocal_mode = 'ham_only'
      end if

      ! K-path settings
      this%auto_kpath = auto_kpath
      this%nk_per_segment = nk_per_segment
      this%override_space_group = override_space_group
      this%custom_kpath_spec = custom_kpath_spec

      ! Log what was read
      call g_logger%info('reciprocal%build_from_file: Read k-mesh = ' // &
                        trim(int2str(nk1)) // ' x ' // trim(int2str(nk2)) // ' x ' // trim(int2str(nk3)), &
                        __FILE__, __LINE__)
      
      ! if (sum(abs(this%k_offset)) > 1.0e-8_rp) then
      !    call g_logger%info('reciprocal%build_from_file: k-offset = [' // &
      !                      trim(real2str(k_offset_x, '(F8.4)')) // ', ' // &
      !                      trim(real2str(k_offset_y, '(F8.4)')) // ', ' // &
      !                      trim(real2str(k_offset_z, '(F8.4)')) // ']', &
      !                      __FILE__, __LINE__)
      ! end if
      
      if (this%use_symmetry_reduction) then
         call g_logger%info('reciprocal%build_from_file: Symmetry reduction enabled', __FILE__, __LINE__)
      end if
      if (this%use_time_reversal) then
         call g_logger%info('reciprocal%build_from_file: use_time_reversal = true', __FILE__, __LINE__)
      else
         call g_logger%info('reciprocal%build_from_file: use_time_reversal = false', __FILE__, __LINE__)
      end if
      if (this%strict_symmetry_checks) then
         call g_logger%info('reciprocal%build_from_file: strict_symmetry_checks = true', __FILE__, __LINE__)
      else
         call g_logger%info('reciprocal%build_from_file: strict_symmetry_checks = false', __FILE__, __LINE__)
      end if
      
      if (this%auto_kpath) then
         call g_logger%info('reciprocal%build_from_file: Automatic k-path generation enabled', __FILE__, __LINE__)
      end if
      call g_logger%info('reciprocal%build_from_file: reciprocal_mode = ' // trim(this%reciprocal_mode), __FILE__, __LINE__)
   end subroutine build_from_file

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Set k-point mesh parameters for DOS calculations
   !---------------------------------------------------------------------------
   subroutine set_kpoint_mesh(this, nk1, nk2, nk3)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: nk1, nk2, nk3

      this%nk_mesh = [nk1, nk2, nk3]
      call g_logger%info('reciprocal%set_kpoint_mesh: Set k-point mesh to ' // &
         trim(int2str(nk1)) // 'x' // trim(int2str(nk2)) // 'x' // trim(int2str(nk3)), __FILE__, __LINE__)
   end subroutine set_kpoint_mesh

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Generate reciprocal lattice vectors from real space lattice
   !---------------------------------------------------------------------------
   subroutine generate_reciprocal_vectors(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      real(rp), dimension(3, 3) :: real_vectors
      real(rp) :: det
      integer :: i

      ! Get real-space lattice vectors from lattice%a
      real_vectors = this%lattice%a

      ! Calculate determinant (volume of unit cell)
      det = real_vectors(1, 1) * (real_vectors(2, 2) * real_vectors(3, 3) - real_vectors(2, 3) * real_vectors(3, 2)) &
          - real_vectors(1, 2) * (real_vectors(2, 1) * real_vectors(3, 3) - real_vectors(2, 3) * real_vectors(3, 1)) &
          + real_vectors(1, 3) * (real_vectors(2, 1) * real_vectors(3, 2) - real_vectors(2, 2) * real_vectors(3, 1))

      ! Calculate reciprocal lattice vectors: b_i = 2π * (a_j × a_k) / V
      ! b1 = 2π * (a2 × a3) / V
      this%reciprocal_vectors(1, 1) = two_pi * (real_vectors(2, 2) * real_vectors(3, 3) - real_vectors(2, 3) * real_vectors(3, 2)) / det
      this%reciprocal_vectors(2, 1) = two_pi * (real_vectors(2, 3) * real_vectors(3, 1) - real_vectors(2, 1) * real_vectors(3, 3)) / det
      this%reciprocal_vectors(3, 1) = two_pi * (real_vectors(2, 1) * real_vectors(3, 2) - real_vectors(2, 2) * real_vectors(3, 1)) / det

      ! b2 = 2π * (a3 × a1) / V
      this%reciprocal_vectors(1, 2) = two_pi * (real_vectors(3, 2) * real_vectors(1, 3) - real_vectors(3, 3) * real_vectors(1, 2)) / det
      this%reciprocal_vectors(2, 2) = two_pi * (real_vectors(3, 3) * real_vectors(1, 1) - real_vectors(3, 1) * real_vectors(1, 3)) / det
      this%reciprocal_vectors(3, 2) = two_pi * (real_vectors(3, 1) * real_vectors(1, 2) - real_vectors(3, 2) * real_vectors(1, 1)) / det

      ! b3 = 2π * (a1 × a2) / V
      this%reciprocal_vectors(1, 3) = two_pi * (real_vectors(1, 2) * real_vectors(2, 3) - real_vectors(1, 3) * real_vectors(2, 2)) / det
      this%reciprocal_vectors(2, 3) = two_pi * (real_vectors(1, 3) * real_vectors(2, 1) - real_vectors(1, 1) * real_vectors(2, 3)) / det
      this%reciprocal_vectors(3, 3) = two_pi * (real_vectors(1, 1) * real_vectors(2, 2) - real_vectors(1, 2) * real_vectors(2, 1)) / det

      this%reciprocal_volume = (two_pi)**3 / abs(det)

      call g_logger%info('reciprocal%generate_reciprocal_vectors: Reciprocal lattice vectors generated', __FILE__, __LINE__)
      
      ! Debug output - but use info level for now since debug is not enabled
      call g_logger%info('reciprocal%generate_reciprocal_vectors: Real cell volume = ' // real2str(det), __FILE__, __LINE__)
      call g_logger%info('reciprocal%generate_reciprocal_vectors: Reciprocal b1 = [' // &
         real2str(this%reciprocal_vectors(1, 1)) // ', ' // real2str(this%reciprocal_vectors(2, 1)) // ', ' // &
         real2str(this%reciprocal_vectors(3, 1)) // ']', __FILE__, __LINE__)
      call g_logger%info('reciprocal%generate_reciprocal_vectors: Reciprocal b2 = [' // &
         real2str(this%reciprocal_vectors(1, 2)) // ', ' // real2str(this%reciprocal_vectors(2, 2)) // ', ' // &
         real2str(this%reciprocal_vectors(3, 2)) // ']', __FILE__, __LINE__)
      call g_logger%info('reciprocal%generate_reciprocal_vectors: Reciprocal b3 = [' // &
         real2str(this%reciprocal_vectors(1, 3)) // ', ' // real2str(this%reciprocal_vectors(2, 3)) // ', ' // &
         real2str(this%reciprocal_vectors(3, 3)) // ']', __FILE__, __LINE__)
   end subroutine generate_reciprocal_vectors

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Generate Monkhorst-Pack k-point mesh
   !---------------------------------------------------------------------------
   subroutine generate_mp_mesh(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ik, ix, iy, iz, nk_irred
      real(rp) :: kx, ky, kz
      real(rp), parameter :: tol = 1.0e-8_rp
      integer :: shift(3), nfull
      real(rp), allocatable :: kpoints_full(:,:)

      this%nk_total = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)

      ! Allocate k-point arrays
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.k_points', this%k_points, [3, this%nk_total])
      call g_safe_alloc%allocate('reciprocal.k_weights', this%k_weights, [this%nk_total])
#else
       if (allocated(this%k_points)) deallocate(this%k_points)
       if (allocated(this%k_weights)) deallocate(this%k_weights)
      allocate(this%k_points(3, this%nk_total))
      allocate(this%k_weights(this%nk_total))
#endif
      if (allocated(this%full_to_irred_k)) deallocate(this%full_to_irred_k)
      if (allocated(this%irred_to_full_k)) deallocate(this%irred_to_full_k)
      allocate(this%full_to_irred_k(this%nk_total))
      allocate(this%irred_to_full_k(this%nk_total))

      shift = [0, 0, 0]
      if (abs(this%k_offset(1)) > tol) shift(1) = 1
      if (abs(this%k_offset(2)) > tol) shift(2) = 1
      if (abs(this%k_offset(3)) > tol) shift(3) = 1

#ifdef USE_SPGLIB
      if (this%symmetry_analysis%spglib%is_available()) then
         nfull = this%symmetry_analysis%spglib%get_full_kpoint_mesh_with_points(this%nk_mesh, shift, kpoints_full)
         if (nfull == this%nk_total) then
            this%k_points = kpoints_full
            this%k_weights = 1.0_rp / real(this%nk_total, rp)
            do ik = 1, this%nk_total
               this%full_to_irred_k(ik) = ik
               this%irred_to_full_k(ik) = ik
            end do
            deallocate(kpoints_full)
            call g_logger%info('reciprocal%generate_mp_mesh: Generated full mesh via spglib grid convention (' // &
                               trim(int2str(this%nk_total)) // ' k-points)', __FILE__, __LINE__)
            return
         end if
         if (allocated(kpoints_full)) deallocate(kpoints_full)
      end if
#endif

      ! Fallback: internal MP mesh formula
      ik = 0
      do iz = 1, this%nk_mesh(3)
         do iy = 1, this%nk_mesh(2)
            do ix = 1, this%nk_mesh(1)
               ik = ik + 1
               kx = (2.0_rp * ix - this%nk_mesh(1) - 1.0_rp) / (2.0_rp * this%nk_mesh(1)) + this%k_offset(1)
               ky = (2.0_rp * iy - this%nk_mesh(2) - 1.0_rp) / (2.0_rp * this%nk_mesh(2)) + this%k_offset(2)
               kz = (2.0_rp * iz - this%nk_mesh(3) - 1.0_rp) / (2.0_rp * this%nk_mesh(3)) + this%k_offset(3)
               this%k_points(1, ik) = kx
               this%k_points(2, ik) = ky
               this%k_points(3, ik) = kz
               this%k_weights(ik) = 1.0_rp / real(this%nk_total, rp)
               this%full_to_irred_k(ik) = ik
               this%irred_to_full_k(ik) = ik
            end do
         end do
      end do

      call g_logger%info('reciprocal%generate_mp_mesh: Generated Monkhorst-Pack mesh with ' // trim(int2str(this%nk_total)) // ' k-points', __FILE__, __LINE__)
   end subroutine generate_mp_mesh

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Set basis sizes for different atom types based on orbital configuration
   !---------------------------------------------------------------------------
   subroutine set_basis_sizes(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ntype
      character(len=10) :: basis_type

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.basis_size', this%basis_size, [this%lattice%ntype])
#else
    if (allocated(this%basis_size)) deallocate(this%basis_size)
      allocate(this%basis_size(this%lattice%ntype))
#endif

      ! Determine basis size for each atom type from active global basis.
      do ntype = 1, this%lattice%ntype
         this%basis_size(ntype) = norb
      end do

      this%max_orbs = nb

      call g_logger%info('reciprocal%set_basis_sizes: Basis sizes set: max_orb_channels = ' // trim(int2str(this%max_orbs)), __FILE__, __LINE__)
   end subroutine set_basis_sizes

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build neighbor vectors for each atom type using clusba
   !> This provides proper per-atom neighbor vectors for multi-site H_k generation
   !> Following the pattern from hamiltonian.f90 chbar_nc routine
   !---------------------------------------------------------------------------
   subroutine build_neighbor_vectors(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ntype, ia, nr, nn_max_loc, kk
      real(rp) :: r2
      real(rp), dimension(3, this%lattice%kk) :: cralat

      if (allocated(this%ham_vec_type) .and. allocated(this%ham_vec_type_direct)) then
         if (size(this%ham_vec_type, 1) == 3 .and. &
             size(this%ham_vec_type, 2) == this%lattice%nn_max .and. &
             size(this%ham_vec_type, 3) == this%lattice%ntype) then
            return
         end if
      end if

      call g_logger%info('reciprocal%build_neighbor_vectors: Building neighbor vectors for each atom type', __FILE__, __LINE__)

      r2 = this%lattice%r2
      kk = this%lattice%kk
      cralat(1:3, 1:kk) = this%lattice%cr(1:3, 1:kk) * this%lattice%alat

      ! Allocate storage for all atom types
      ! Use nn_max as maximum neighbors
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.ham_vec_type', this%ham_vec_type, &
                                 [3, this%lattice%nn_max, this%lattice%ntype])
      call g_safe_alloc%allocate('reciprocal.ham_vec_type_direct', this%ham_vec_type_direct, &
                                 [3, this%lattice%nn_max, this%lattice%ntype])
#else
      if (allocated(this%ham_vec_type)) deallocate(this%ham_vec_type)
      allocate(this%ham_vec_type(3, this%lattice%nn_max, this%lattice%ntype))
      if (allocated(this%ham_vec_type_direct)) deallocate(this%ham_vec_type_direct)
      allocate(this%ham_vec_type_direct(3, this%lattice%nn_max, this%lattice%ntype))
#endif

      this%ham_vec_type = 0.0_rp
      this%ham_vec_type_direct = 0.0_rp

      ! Build neighbor vectors for each atom type
      do ntype = 1, this%lattice%ntype
         ia = this%lattice%atlist(ntype)
         nr = this%lattice%nn(ia, 1)
         nn_max_loc = nr

         ! Use clusba directly - no need for lattice%sbarvec storage here
         ! We're building type-specific ham_vec_type arrays instead
         call this%lattice%clusba(r2, cralat, ia, kk, kk, nn_max_loc, &
                                  this%ham_vec_type(:, 1:nr, ntype))

         ! Update nr with actual number of neighbors found
         nr = nn_max_loc

         ! Convert to fractional coordinates for ham_vec_type_direct
         ! Note: ham_vec_type from clusba is already in absolute Cartesian units (cralat = cr * alat)
         do nn_max_loc = 1, nr
            if (this%lattice%a_cart_inv_ready) then
               ! a_cart_inv expects input in units of alat, so divide by alat
               this%ham_vec_type_direct(:, nn_max_loc, ntype) = &
                  matmul(this%lattice%a_cart_inv, this%ham_vec_type(:, nn_max_loc, ntype) ) !/ this%lattice%alat
            else
               ! Fallback: use inverse of lattice vectors directly.
               this%ham_vec_type_direct(:, nn_max_loc, ntype) = &
                  matmul(inverse_3x3(this%lattice%a), this%ham_vec_type(:, nn_max_loc, ntype) / this%lattice%alat)
            end if
         end do

         call g_logger%info('reciprocal%build_neighbor_vectors: Built ' // trim(int2str(nr)) // &
                          ' neighbor vectors for atom type ' // trim(int2str(ntype)), __FILE__, __LINE__)
         
      end do

      call g_logger%info('reciprocal%build_neighbor_vectors: Completed neighbor vector build for all types', __FILE__, __LINE__)

   end subroutine build_neighbor_vectors

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate structure factors for ALL atom types for Fourier transform
   !> Fills structure_factors(ineigh, ntype) for all neighbors and types
   !---------------------------------------------------------------------------
   subroutine calculate_structure_factors(this, k_vec, structure_factors)
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

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Fourier transform real-space Hamiltonian to k-space for multi-site system
   !> H(k) = Σ_{i,j sites} Σ_R H_ij(R) * exp(i*k·R_ij)
   !> where R_ij = R + r_j - r_i (lattice vector + basis positions)
   !---------------------------------------------------------------------------
   subroutine fourier_transform_hamiltonian(this, k_vec, hk_result)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: hk_result  ! (n_orb*n_sites, n_orb*n_sites)
      ! Local variables
      integer :: isite, jsite, ntype_i, ntype_j, ineigh, ia, ja, nr
      integer :: i_start, i_end, j_start, j_end
      integer :: n_orb, n_sites
      complex(rp), dimension(:, :), allocatable :: structure_factors  ! (nr_max, ntype)
      
      ! Get dimensions
      n_orb = nb
      n_sites = this%lattice%nrec
      
      ! Allocate structure factors for all types
      allocate(structure_factors(this%lattice%nn_max, this%lattice%ntype))
      
      ! Calculate structure factors for all atom types at this k-point
      call this%calculate_structure_factors(k_vec, structure_factors)
      
      ! Initialize result
      hk_result = cmplx(0.0_rp, 0.0_rp, rp)

      ! Loop over all sites in the unit cell
      ! For each i_site → j_site pair, sum over lattice vectors R
      do isite = 1, n_sites
         ntype_i = this%lattice%ib(isite)  ! Type of site i
         ia = this%lattice%atlist(ntype_i) ! Cluster atom for this type
         nr = this%lattice%nn(ia, 1)       ! Number of neighbors
         
         ! Orbital block indices for site i (row block)
         i_start = (isite - 1) * n_orb + 1
         i_end = isite * n_orb
         
         ! Loop over neighbors (which map to j_sites in different cells)
         do ineigh = 1, nr
            if (ineigh == 1) then
               ! On-site: i_site = j_site, R = 0
               jsite = isite
            else
               ! Off-site: determine which j_site this neighbor corresponds to
               ja = this%lattice%nn(ia, ineigh)     ! Cluster atom index
               if (ja < 1 .or. ja > this%lattice%kk) cycle
               jsite = this%lattice%iz(ja)          ! Map to unit cell site
               if (jsite < 1 .or. jsite > n_sites) cycle
            end if
            
            ! Orbital block indices for site j (column block)
            j_start = (jsite - 1) * n_orb + 1
            j_end = jsite * n_orb
            
            ! Add H_ij(R) * exp(i*k·R_ij) to the appropriate block
            ! Structure factor already contains the phase for this neighbor
            hk_result(i_start:i_end, j_start:j_end) = &
               hk_result(i_start:i_end, j_start:j_end) + &
               this%hamiltonian%ee(:, :, ineigh, ntype_i) * structure_factors(ineigh, ntype_i)
         end do
      end do

      deallocate(structure_factors)
      
   end subroutine fourier_transform_hamiltonian

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Fourier transform real-space overlap proxy to k-space.
   !> Current proxy uses S(R) = I + eeo(R) with on-site identity term.
   !---------------------------------------------------------------------------
   subroutine fourier_transform_overlap(this, k_vec, sk_result)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      complex(rp), dimension(:, :), intent(out) :: sk_result
      integer :: isite, jsite, ntype_i, ineigh, ia, ja, nr, iorb
      integer :: i_start, i_end, j_start, j_end
      integer :: n_orb, n_sites
      complex(rp), dimension(:, :), allocatable :: structure_factors
      complex(rp), dimension(:, :), allocatable :: overlap_block

      n_orb = nb
      n_sites = this%lattice%nrec
      allocate(structure_factors(this%lattice%nn_max, this%lattice%ntype))
      allocate(overlap_block(n_orb, n_orb))
      call this%calculate_structure_factors(k_vec, structure_factors)

      sk_result = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, n_sites
         ntype_i = this%lattice%ib(isite)
         ia = this%lattice%atlist(ntype_i)
         nr = this%lattice%nn(ia, 1)
         i_start = (isite - 1) * n_orb + 1
         i_end = isite * n_orb
         do ineigh = 1, nr
            if (ineigh == 1) then
               jsite = isite
            else
               ja = this%lattice%nn(ia, ineigh)
               if (ja < 1 .or. ja > this%lattice%kk) cycle
               jsite = this%lattice%iz(ja)
               if (jsite < 1 .or. jsite > n_sites) cycle
            end if
            j_start = (jsite - 1) * n_orb + 1
            j_end = jsite * n_orb
            overlap_block(:, :) = this%hamiltonian%eeo(:, :, ineigh, ntype_i)
            if (ineigh == 1 .and. jsite == isite) then
               do iorb = 1, n_orb
                  overlap_block(iorb, iorb) = overlap_block(iorb, iorb) + cmplx(1.0_rp, 0.0_rp, rp)
               end do
            end if
            sk_result(i_start:i_end, j_start:j_end) = sk_result(i_start:i_end, j_start:j_end) + &
               overlap_block * structure_factors(ineigh, ntype_i)
         end do
      end do
      deallocate(overlap_block, structure_factors)
   end subroutine fourier_transform_overlap

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build k-space Hamiltonian for all k-points (bulk contribution)
   !---------------------------------------------------------------------------
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build k-space Hamiltonian by Fourier transforming real-space H
   !>
   !> This routine builds H(k) for either:
   !>   - k-mesh (for DOS/SCF): uses k_points array
   !>   - k-path (for bands): uses k_path array
   !>
   !> Workflow:
   !>   1. Generate k-points (generate_kmesh or set up k-path)
   !>   2. Call build_kspace_hamiltonian() to compute hk_bulk
   !>   3. Call diagonalize_hamiltonian() to get eigenvalues/eigenvectors
   !---------------------------------------------------------------------------
   subroutine build_kspace_hamiltonian(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ik, nk, ntype
      character(len=200) :: debug_msg
      logical :: using_kpath
      integer :: i, j

      ! Determine which k-point set to use
      using_kpath = .false.
      if (allocated(this%k_path)) then
         ! Use k-path for band structure
         nk = this%nk_path
         using_kpath = .true.
         call g_logger%info('reciprocal%build_kspace_hamiltonian: Building H(k) for k-path', __FILE__, __LINE__)
      else if (allocated(this%k_points)) then
         ! Use k-mesh for DOS/SCF
         nk = this%nk_total
         call g_logger%info('reciprocal%build_kspace_hamiltonian: Building H(k) for k-mesh', __FILE__, __LINE__)
      else
         call g_logger%error('reciprocal%build_kspace_hamiltonian: No k-points generated. ' // &
                           'Call generate_mp_mesh or generate k-path first.', __FILE__, __LINE__)
         return
      end if

      write(debug_msg, '(A,I0,A,I0,A)') 'build_kspace_hamiltonian: Building for ', nk, &
                                        ' k-points and ', this%lattice%ntype, ' atom types'
      call g_logger%info(trim(debug_msg), __FILE__, __LINE__)

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
                                [this%max_orbs*this%lattice%nrec, this%max_orbs*this%lattice%nrec, nk])
#else
   if (allocated(this%hk_bulk)) deallocate(this%hk_bulk)
   allocate(this%hk_bulk(this%max_orbs*this%lattice%nrec, this%max_orbs*this%lattice%nrec, nk))
#endif

      ! Parallelize over k-points (coarse-grained) when OpenMP is available.
#ifdef _OPENMP
      !$omp parallel do private(ik) shared(this, nk, using_kpath) default(none)
#endif
      do ik = 1, nk
         ! Fourier transform: builds full multi-site H(k)
         if (using_kpath) then
            call this%fourier_transform_hamiltonian(this%k_path(:, ik), this%hk_bulk(:, :, ik))
         else
            call this%fourier_transform_hamiltonian(this%k_points(:, ik), this%hk_bulk(:, :, ik))
         end if
      end do
#ifdef _OPENMP
      !$omp end parallel do
#endif

      call g_logger%info('reciprocal%build_kspace_hamiltonian: K-space Hamiltonian built successfully', __FILE__, __LINE__)
      if (trim(this%reciprocal_mode) == 'generalized_overlap_proxy') then
         call this%build_kspace_overlap()
      else if (trim(this%reciprocal_mode) == 'generalized_overlap_kanpur') then
         call g_logger%warning('reciprocal%build_kspace_hamiltonian: generalized_overlap_kanpur requested but not implemented yet. Using ham_only solve path.', __FILE__, __LINE__)
      end if
      
      ! Diagnostic: Check H(k) at Gamma point for multi-site systems
      if (this%lattice%nrec > 1 .and. nk > 0) then
         call this%check_multisite_hamiltonian_diagonal()
      end if
   end subroutine build_kspace_hamiltonian

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build overlap proxy S(k) for all k-points.
   !---------------------------------------------------------------------------
   subroutine build_kspace_overlap(this)
      class(reciprocal), intent(inout) :: this
      integer :: ik, nk
      logical :: using_kpath

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

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.sk_overlap', this%sk_overlap, [this%max_orbs*this%lattice%nrec, this%max_orbs*this%lattice%nrec, nk])
#else
      if (allocated(this%sk_overlap)) deallocate(this%sk_overlap)
      allocate(this%sk_overlap(this%max_orbs*this%lattice%nrec, this%max_orbs*this%lattice%nrec, nk))
#endif

#ifdef _OPENMP
      !$omp parallel do private(ik) shared(this, nk, using_kpath) default(none)
#endif
      do ik = 1, nk
         if (using_kpath) then
            call this%fourier_transform_overlap(this%k_path(:, ik), this%sk_overlap(:, :, ik))
         else
            call this%fourier_transform_overlap(this%k_points(:, ik), this%sk_overlap(:, :, ik))
         end if
      end do
#ifdef _OPENMP
      !$omp end parallel do
#endif
      call g_logger%info('reciprocal%build_kspace_overlap: Built S(k) overlap proxy.', __FILE__, __LINE__)
   end subroutine build_kspace_overlap

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief
!> Diagnostic: Check diagonal blocks of multi-site H(k) at Gamma point
!---------------------------------------------------------------------------
subroutine check_multisite_hamiltonian_diagonal(this)
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
   call g_logger%info(trim(msg), __FILE__, __LINE__)
   
   write(msg, '(A,2F12.6)') 'H(k=0) diagonal check - Site 2 d-orbital avg (up/dn): ', &
                             real(h_avg_site2_up), real(h_avg_site2_dn)
   call g_logger%info(trim(msg), __FILE__, __LINE__)
   
   write(msg, '(A,F12.6)') 'H(k=0) diagonal difference (site2-site1) up: ', &
                            real(h_avg_site2_up - h_avg_site1_up)
   call g_logger%info(trim(msg), __FILE__, __LINE__)
   
end subroutine check_multisite_hamiltonian_diagonal

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief
!> Check Hermiticity of H(k) for debugging multi-site assembly
!---------------------------------------------------------------------------
subroutine check_hamiltonian_hermiticity(this, ik)
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

!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief
!> Print structure of H(k) to diagnose multi-site assembly
!---------------------------------------------------------------------------
subroutine print_hamiltonian_structure(this, ik)
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

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Get basis type string from basis size
   !---------------------------------------------------------------------------
   function get_basis_type_from_size(this, ntype) result(basis_type)
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



   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Diagonalize the pre-built k-space Hamiltonian hk_bulk at all k-points
   !> 
   !> This routine uses the Hamiltonian matrices that were previously computed
   !> by build_kspace_hamiltonian(). 
   !>
   !> Workflow:
   !>   1. Generate k-points (generate_kmesh or generate k-path)
   !>   2. Call build_kspace_hamiltonian() to compute hk_bulk
   !>   3. Call diagonalize_hamiltonian() to get eigenvalues/eigenvectors
   !>
   !> Parallelized with OpenMP over k-points for efficiency.
   !---------------------------------------------------------------------------
   subroutine diagonalize_hamiltonian(this)
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
      if (.not. allocated(this%hk_bulk)) then
         call g_logger%error('diagonalize_hamiltonian: hk_bulk not built - call build_kspace_hamiltonian first', &
                            __FILE__, __LINE__)
         return
      end if

      ! Get dimensions from hk_bulk
      nmat = size(this%hk_bulk, 1)
      nk = size(this%hk_bulk, 3)
      
      call g_logger%info('diagonalize_hamiltonian: Diagonalizing ' // trim(int2str(nk)) // ' k-points', __FILE__, __LINE__)
      call g_logger%info('diagonalize_hamiltonian: Matrix size = ' // &
                        trim(int2str(nmat)) // ' x ' // trim(int2str(nmat)), __FILE__, __LINE__)

      use_generalized = trim(this%reciprocal_mode) == 'generalized_overlap_proxy'
      if (use_generalized) then
         if (.not. allocated(this%sk_overlap)) call this%build_kspace_overlap()
         if (.not. allocated(this%sk_overlap)) then
            call g_logger%warning('diagonalize_hamiltonian: S(k) unavailable, falling back to ham_only.', __FILE__, __LINE__)
            use_generalized = .false.
         end if
      end if

      if (this%kanpur_diagnostics) call this%print_kanpur_mapping()

      ! Allocate eigenvalue and eigenvector storage
      if (allocated(this%eigenvalues)) deallocate(this%eigenvalues)
      if (allocated(this%eigenvectors)) deallocate(this%eigenvectors)
      allocate(this%eigenvalues(nmat, nk))
      allocate(this%eigenvectors(nmat, nmat, nk))
      mode_fail_count = 0

      ! Parallel diagonalization over k-points
      ! Each thread needs its own LAPACK workspace
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP& PRIVATE(ik, h_k_copy, s_k_copy, eigenvals, work_complex, rwork, lwork, info, info_msg) &
!$OMP& IF(nk > 10)
      
      ! Allocate thread-private work arrays
      allocate(h_k_copy(nmat, nmat))
      allocate(s_k_copy(nmat, nmat))
      allocate(eigenvals(nmat))
      allocate(rwork(3*nmat - 2))
      
      ! Query optimal LAPACK workspace size
      lwork = -1
      allocate(work_complex(1))
      if (use_generalized) then
         call zhegv(1, 'V', 'U', nmat, h_k_copy, nmat, s_k_copy, nmat, eigenvals, work_complex, lwork, rwork, info)
      else
         call zheev('V', 'U', nmat, h_k_copy, nmat, eigenvals, work_complex, lwork, rwork, info)
      end if
      lwork = int(real(work_complex(1)))
      deallocate(work_complex)
      allocate(work_complex(lwork))

      ! Loop over all k-points - use pre-computed hk_bulk
!$OMP DO SCHEDULE(DYNAMIC, 10)
      do ik = 1, nk
         ! Copy H(k) from pre-computed array
         h_k_copy = this%hk_bulk(:, :, ik)

         if (use_generalized) then
            s_k_copy = this%sk_overlap(:, :, ik)
            call this%check_overlap_properties(ik, s_k_copy)
            call zhegv(1, 'V', 'U', nmat, h_k_copy, nmat, s_k_copy, nmat, eigenvals, work_complex, lwork, rwork, info)
         else
            ! Diagonalize H(k) using LAPACK ZHEEV
            call zheev('V', 'U', nmat, h_k_copy, nmat, eigenvals, work_complex, lwork, rwork, info)
         end if
         
         if (info /= 0) then
            write(info_msg, '(A,I0,A,I0)') 'Diagonalization failed at k-point ', ik, ', info = ', info
            call g_logger%error('diagonalize_hamiltonian: ' // trim(info_msg), __FILE__, __LINE__)
!$OMP CRITICAL
            mode_fail_count = mode_fail_count + 1
!$OMP END CRITICAL
            cycle
         end if

         ! Store eigenvalues and eigenvectors
         this%eigenvalues(:, ik) = eigenvals
         this%eigenvectors(:, :, ik) = h_k_copy
      end do
!$OMP END DO
      
      ! Cleanup thread-private arrays
      deallocate(h_k_copy, s_k_copy, eigenvals, work_complex, rwork)
      
!$OMP END PARALLEL

      if (mode_fail_count > 0 .and. use_generalized) then
         call g_logger%warning('diagonalize_hamiltonian: generalized_overlap had failures; consider ham_only for robustness.', __FILE__, __LINE__)
      end if
      if (this%gamma_bounds_diagnostics) call this%run_gamma_bounds_diagnostics()
      if (this%hall_diag_experimental) call this%diagonalize_hall_experimental()
      call g_logger%info('diagonalize_hamiltonian: Completed successfully', __FILE__, __LINE__)
   end subroutine diagonalize_hamiltonian

   subroutine print_kanpur_mapping(this)
      class(reciprocal), intent(in) :: this
      call g_logger%info('Kanpur mapping: reciprocal_mode=' // trim(this%reciprocal_mode), __FILE__, __LINE__)
      if (trim(this%reciprocal_mode) == 'generalized_overlap_proxy') then
         call g_logger%info('Kanpur mapping: PROXY generalized solve H(k)c = E S_proxy(k)c, S_proxy from eeo + I.', __FILE__, __LINE__)
         call g_logger%warning('Kanpur mapping: PROXY mode is not a formal Kanpur LMTO overlap representation.', __FILE__, __LINE__)
      else if (trim(this%reciprocal_mode) == 'generalized_overlap_kanpur') then
         call g_logger%warning('Kanpur mapping: generalized_overlap_kanpur selected but not implemented; currently falling back to ham_only.', __FILE__, __LINE__)
      else
         call g_logger%info('Kanpur mapping: Hamiltonian-only mode (TB-like).', __FILE__, __LINE__)
      end if
      call g_logger%info('Kanpur mapping: non-orthogonality treatment is approximation-level diagnostic.', __FILE__, __LINE__)
   end subroutine print_kanpur_mapping

   subroutine check_overlap_properties(this, ik, s_k)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: ik
      complex(rp), dimension(:, :), intent(in) :: s_k
      integer :: i, j, n
      real(rp) :: max_herm
      max_herm = 0.0_rp
      n = size(s_k, 1)
      do i = 1, n
         do j = 1, n
            max_herm = max(max_herm, abs(s_k(i, j) - conjg(s_k(j, i))))
         end do
      end do
      if (ik == 1 .or. max_herm > 1.0e-6_rp) then
         call g_logger%info('S(k) hermiticity check ik=' // trim(int2str(ik)) // ' max_diff=' // &
            trim(real2str(max_herm, '(ES12.4)')), __FILE__, __LINE__)
      end if
   end subroutine check_overlap_properties

   subroutine run_gamma_bounds_diagnostics(this)
      class(reciprocal), intent(inout) :: this
      complex(rp), allocatable :: h_gamma(:, :)
      real(rp) :: egmin, egmax
      type(bounds) :: bnd

      if (.not. allocated(this%hk_bulk)) return
      allocate(h_gamma(size(this%hk_bulk, 1), size(this%hk_bulk, 2)))
      h_gamma = this%hk_bulk(:, :, 1)
      call compute_spectrum_bounds(h_gamma, bnd, method='both', verbose=.false.)
      call g_logger%info('Gamma bounds diagnostic: Gershgorin [' // &
         trim(real2str(bnd%e_min_gershgorin, '(F12.6)')) // ', ' // &
         trim(real2str(bnd%e_max_gershgorin, '(F12.6)')) // ']', __FILE__, __LINE__)
      call this%diagonalize_hall_experimental() ! reuse exact solver utility logs for finite matrix check
      if (allocated(this%eigenvalues)) then
         egmin = minval(this%eigenvalues(:, 1))
         egmax = maxval(this%eigenvalues(:, 1))
         call g_logger%info('Gamma bounds diagnostic: eig(H(Gamma))=[' // &
            trim(real2str(egmin, '(F12.6)')) // ', ' // trim(real2str(egmax, '(F12.6)')) // ']', __FILE__, __LINE__)
      end if
      deallocate(h_gamma)
   end subroutine run_gamma_bounds_diagnostics

   subroutine diagonalize_hall_experimental(this)
      class(reciprocal), intent(inout) :: this
      integer :: nsites, n_orb, n, i, jsite, isite, ineigh, ia, ja, nr, info, lwork
      integer :: i_start, i_end, j_start, j_end
      complex(rp), allocatable :: hall_mat(:, :), work(:)
      real(rp), allocatable :: evals(:), rwork(:)

      if (.not. this%hall_diag_experimental) return
      if (this%control%calctype /= 'I') then
         call g_logger%info('HALL experimental diagonalization skipped: calctype is not impurity.', __FILE__, __LINE__)
         return
      end if

      nsites = this%lattice%nmax
      n_orb = 18
      if (nsites <= 0) return
      n = nsites * n_orb
      allocate(hall_mat(n, n), evals(n), rwork(max(1, 3*n - 2)))
      hall_mat = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, nsites
         ia = isite
         nr = this%lattice%nn(ia, 1)
         i_start = (isite - 1) * n_orb + 1
         i_end = isite * n_orb
         do ineigh = 1, nr
            if (ineigh == 1) then
               jsite = isite
            else
               ja = this%lattice%nn(ia, ineigh)
               jsite = ja
               if (jsite < 1 .or. jsite > nsites) cycle
            end if
            j_start = (jsite - 1) * n_orb + 1
            j_end = jsite * n_orb
            hall_mat(i_start:i_end, j_start:j_end) = hall_mat(i_start:i_end, j_start:j_end) + this%hamiltonian%hall(:, :, ineigh, isite)
         end do
      end do
      allocate(work(1))
      call zheev('N', 'U', n, hall_mat, n, evals, work, -1, rwork, info)
      lwork = max(1, int(real(work(1))))
      deallocate(work)
      allocate(work(lwork))
      call zheev('N', 'U', n, hall_mat, n, evals, work, lwork, rwork, info)
      if (info == 0) then
         call g_logger%info('HALL experimental eig range: [' // trim(real2str(minval(evals), '(F12.6)')) // ', ' // &
            trim(real2str(maxval(evals), '(F12.6)')) // '] (diagnostic only)', __FILE__, __LINE__)
      else
         call g_logger%warning('HALL experimental diagonalization failed, info=' // trim(int2str(info)), __FILE__, __LINE__)
      end if
      deallocate(hall_mat, evals, rwork, work)
   end subroutine diagonalize_hall_experimental


   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate band structure along high-symmetry path
   !---------------------------------------------------------------------------
   subroutine calculate_band_structure(this, ham, crystal_type, npts_per_segment, output_file)
      class(reciprocal), intent(inout) :: this
      class(hamiltonian), intent(in) :: ham
      character(len=*), intent(in) :: crystal_type
      integer, intent(in) :: npts_per_segment
      character(len=*), intent(in), optional :: output_file
      ! Local variables
      character(len=256) :: filename
      integer :: unit, i, j, nmat
      character(len=100) :: fmt_str

      call g_logger%info('calculate_band_structure: Starting band structure calculation', __FILE__, __LINE__)

      ! Generate k-path using spglib-based canonical path generation
      call this%symmetry_analysis%generate_canonical_kpath(npts_per_segment)

      ! Copy k-path data from symmetry analysis to reciprocal object
      if (allocated(this%symmetry_analysis%k_path)) then
         ! Copy k-path arrays
         this%nk_path = this%symmetry_analysis%nk_path
         if (allocated(this%k_path)) deallocate(this%k_path)
         if (allocated(this%k_labels)) deallocate(this%k_labels)
         if (allocated(this%k_distances)) deallocate(this%k_distances)
         
         allocate(this%k_path(3, this%nk_path))
         allocate(this%k_labels(this%nk_path))
         allocate(this%k_distances(this%nk_path))
         
         this%k_path = this%symmetry_analysis%k_path
         this%k_labels = this%symmetry_analysis%k_labels
         this%k_distances = this%symmetry_analysis%k_distances
         
         call g_logger%info('calculate_band_structure: Copied k-path with ' // &
                           trim(int2str(this%nk_path)) // ' points from symmetry analysis', __FILE__, __LINE__)
      else
         call g_logger%error('calculate_band_structure: Symmetry analysis k_path not allocated!', __FILE__, __LINE__)
         return
      end if

      ! Debug logging for k-path
      call g_logger%info('calculate_band_structure: After k-path generation, nk_path = ' // &
                        trim(int2str(this%nk_path)), __FILE__, __LINE__)
      if (allocated(this%k_path)) then
         call g_logger%info('calculate_band_structure: k_path allocated with dimensions ' // &
                           trim(int2str(size(this%k_path,1))) // 'x' // &
                           trim(int2str(size(this%k_path,2))), __FILE__, __LINE__)
      else
         call g_logger%error('calculate_band_structure: k_path not allocated after k-path generation!', __FILE__, __LINE__)
      end if

      ! Build k-space Hamiltonian for the k-path
      call this%build_kspace_hamiltonian()

      ! Diagonalize Hamiltonian along k-path
      call this%diagonalize_hamiltonian()

      ! Set output filename
      filename = 'band_structure.dat'
      if (present(output_file)) filename = output_file

      ! Write band structure to file
      open(newunit=unit, file=trim(filename), status='replace', action='write')
      
      ! Write header
      write(unit, '(A)') '# Band structure calculation'
      write(unit, '(A)') '# Crystal type: ' // trim(crystal_type)
      write(unit, '(A,I0)') '# Number of k-points: ', this%nk_path
      nmat = size(this%eigenvalues, 1)
      write(unit, '(A,I0)') '# Number of bands: ', nmat
      write(unit, '(A)') '# Format: k_distance, eigenvalue_1, eigenvalue_2, ...'
      write(unit, '(A)') '#'

      ! Write data
      write(fmt_str, '(A,I0,A)') '(', nmat+1, '(ES16.8E3,1X))'
      do i = 1, this%nk_path
         write(unit, fmt_str) this%k_distances(i), (this%eigenvalues(j, i), j = 1, nmat)
      end do

      close(unit)

      call g_logger%info('calculate_band_structure: Band structure written to ' // trim(filename), __FILE__, __LINE__)

      ! Also write k-path information
      filename = 'kpath_info.dat'
      if (present(output_file)) then
         ! Replace extension with _kpath.dat
         i = index(output_file, '.', back=.true.)
         if (i > 0) then
            filename = output_file(1:i-1) // '_kpath.dat'
         else
            filename = trim(output_file) // '_kpath.dat'
         end if
      end if

      open(newunit=unit, file=trim(filename), status='replace', action='write')
      write(unit, '(A)') '# K-path information'
      write(unit, '(A)') '# Format: k_distance, kx, ky, kz, label'
      write(unit, '(A)') '#'
      
      do i = 1, this%nk_path
         if (trim(this%k_labels(i)) /= '') then
            write(unit, '(4(ES16.8E3,1X),A)') this%k_distances(i), this%k_path(:, i), ' # ' // trim(this%k_labels(i))
         else
            write(unit, '(4(ES16.8E3,1X))') this%k_distances(i), this%k_path(:, i)
         end if
      end do

      close(unit)

      call g_logger%info('calculate_band_structure: K-path information written to ' // trim(filename), __FILE__, __LINE__)
   end subroutine calculate_band_structure


   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate band structure with automatic k-path detection using spglib
   !> This is a simplified interface that auto-detects crystal structure
   !---------------------------------------------------------------------------
   subroutine calculate_band_structure_auto(this, ham, npts_per_segment, output_file)
      class(reciprocal), intent(inout) :: this
      class(hamiltonian), intent(in) :: ham
      integer, intent(in), optional :: npts_per_segment
      character(len=*), intent(in), optional :: output_file
      
      ! Local variables
      integer :: npts
      character(len=256) :: filename
      integer :: unit, i, j, nmat
      character(len=100) :: fmt_str

      npts = 40  ! Default
      if (present(npts_per_segment)) npts = npts_per_segment

      call g_logger%info('calculate_band_structure_auto: Starting automatic band structure calculation', __FILE__, __LINE__)

      ! Generate k-path automatically using symmetry analysis
      call this%symmetry_analysis%generate_automatic_kpath(npts)

      ! Copy k-path data from symmetry analysis to reciprocal object
      if (allocated(this%symmetry_analysis%k_path)) then
         this%nk_path = this%symmetry_analysis%nk_path
         
         if (allocated(this%k_path)) deallocate(this%k_path)
         if (allocated(this%k_labels)) deallocate(this%k_labels)
         if (allocated(this%k_distances)) deallocate(this%k_distances)
         
         allocate(this%k_path(3, this%nk_path))
         allocate(this%k_labels(this%nk_path))
         allocate(this%k_distances(this%nk_path))
         
         this%k_path = this%symmetry_analysis%k_path
         this%k_labels = this%symmetry_analysis%k_labels
         this%k_distances = this%symmetry_analysis%k_distances
         
         call g_logger%info('calculate_band_structure_auto: K-path has ' // &
                           trim(int2str(this%nk_path)) // ' points', __FILE__, __LINE__)
      else
         call g_logger%error('calculate_band_structure_auto: Failed to generate k-path', __FILE__, __LINE__)
         return
      end if

      ! Build k-space Hamiltonian for the k-path
      call this%build_kspace_hamiltonian()

      ! Diagonalize Hamiltonian along k-path
      call this%diagonalize_hamiltonian()

      ! Set output filename
      filename = 'band_structure.dat'
      if (present(output_file)) filename = output_file

      ! Write band structure to file
      open(newunit=unit, file=trim(filename), status='replace', action='write')
      
      ! Write header with auto-detected information
   write(unit, '(A)') '# Band structure calculation (automatic k-path)'
#ifdef USE_SPGLIB
   if (this%symmetry_analysis%spglib%is_available()) then
      write(unit, '(A)') '# Space group: ' // trim(this%symmetry_analysis%spglib%get_space_group_symbol()) // &
               ' (#' // trim(int2str(this%symmetry_analysis%spglib%get_space_group_number())) // ')'
      write(unit, '(A)') '# Crystal system: ' // trim(this%symmetry_analysis%spglib%get_crystal_system_name())
   end if
#endif
      write(unit, '(A,I0)') '# Number of k-points: ', this%nk_path
      nmat = size(this%eigenvalues_path, 1)
      write(unit, '(A,I0)') '# Number of bands: ', nmat
      write(unit, '(A)') '# Format: k_distance, eigenvalue_1, eigenvalue_2, ...'
      write(unit, '(A)') '#'

      ! Write data
      write(fmt_str, '(A,I0,A)') '(', nmat+1, '(ES16.8E3,1X))'
      do i = 1, this%nk_path
         write(unit, fmt_str) this%k_distances(i), (this%eigenvalues_path(j, i), j = 1, nmat)
      end do

      close(unit)

      call g_logger%info('calculate_band_structure_auto: Band structure written to ' // trim(filename), __FILE__, __LINE__)

      ! Write k-path information
      filename = 'kpath_info.dat'
      if (present(output_file)) then
         i = index(output_file, '.', back=.true.)
         if (i > 0) then
            filename = output_file(1:i-1) // '_kpath.dat'
         else
            filename = trim(output_file) // '_kpath.dat'
         end if
      end if

      open(newunit=unit, file=trim(filename), status='replace', action='write')
      write(unit, '(A)') '# K-path information'
      write(unit, '(A)') '# Format: k_distance, kx, ky, kz, label'
      write(unit, '(A)') '#'
      
      do i = 1, this%nk_path
         if (trim(this%k_labels(i)) /= '') then
            write(unit, '(4(ES16.8E3,1X),A)') this%k_distances(i), this%k_path(:, i), ' # ' // trim(this%k_labels(i))
         else
            write(unit, '(4(ES16.8E3,1X))') this%k_distances(i), this%k_path(:, i)
         end if
      end do

      close(unit)

      call g_logger%info('calculate_band_structure_auto: K-path information written to ' // trim(filename), __FILE__, __LINE__)
   end subroutine calculate_band_structure_auto





   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Generate symmetry-reduced k-point mesh using spglib
   !---------------------------------------------------------------------------
   subroutine generate_reduced_kpoint_mesh(this, mesh_dims, use_shift)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: mesh_dims(3)
      logical, intent(in), optional :: use_shift
      integer :: shift(3)
      integer :: num_ir_kpoints
      logical :: do_shift
      real(rp), allocatable :: kpoints_frac(:,:), weights(:)
      integer :: i
      integer, allocatable :: full_to_irred(:), irred_to_full(:)

      do_shift = .false.
      if (present(use_shift)) do_shift = use_shift

      if (do_shift) then
         shift = [1, 1, 1]  ! Offset by half a mesh spacing
      else
         shift = [0, 0, 0]  ! No offset
      end if

      ! Store mesh dimensions
      this%nk_mesh = mesh_dims

#ifdef USE_SPGLIB
   if (.not. this%symmetry_analysis%spglib%is_available()) then
      call g_logger%warning('generate_reduced_kpoint_mesh: spglib not available, using full mesh', __FILE__, __LINE__)
      call this%generate_mp_mesh()  ! Fall back to regular MP mesh
      return
   end if

   ! Get irreducible k-points and weights from spglib
   num_ir_kpoints = this%symmetry_analysis%spglib%get_reduced_kpoint_mesh_with_points( &
                           mesh_dims, shift, kpoints_frac, weights, this%use_time_reversal, &
                           full_to_irred, irred_to_full)
#else
   ! SPGLIB not enabled at compile time: fall back to full mesh
   call g_logger%warning('generate_reduced_kpoint_mesh: spglib support was not compiled in, using full mesh', __FILE__, __LINE__)
   call this%generate_mp_mesh()
   return
#endif

      if (num_ir_kpoints == 0) then
         call g_logger%error('generate_reduced_kpoint_mesh: Failed to get k-points from spglib', __FILE__, __LINE__)
         call this%generate_mp_mesh()  ! Fall back to regular MP mesh
         return
      end if

      ! Store k-point information
      this%nk_total = num_ir_kpoints

      ! Allocate k-point arrays
#ifdef USE_SAFE_ALLOC
      if (allocated(this%k_points)) call g_safe_alloc%deallocate('reciprocal.k_points', this%k_points)
      if (allocated(this%k_weights)) call g_safe_alloc%deallocate('reciprocal.k_weights', this%k_weights)
      call g_safe_alloc%allocate_real('reciprocal.k_points', this%k_points, (/3, this%nk_total/))
      call g_safe_alloc%allocate_real('reciprocal.k_weights', this%k_weights, (/this%nk_total/))
#else
      if (allocated(this%k_points)) deallocate(this%k_points)
      if (allocated(this%k_weights)) deallocate(this%k_weights)
      allocate(this%k_points(3, this%nk_total))
      allocate(this%k_weights(this%nk_total))
#endif

      ! Copy k-points and weights
      this%k_points = kpoints_frac
      this%k_weights = weights
      if (allocated(this%full_to_irred_k)) deallocate(this%full_to_irred_k)
      if (allocated(this%irred_to_full_k)) deallocate(this%irred_to_full_k)
      allocate(this%full_to_irred_k(size(full_to_irred)))
      allocate(this%irred_to_full_k(size(irred_to_full)))
      this%full_to_irred_k = full_to_irred
      this%irred_to_full_k = irred_to_full

      ! Clean up temporary arrays
      deallocate(kpoints_frac, weights, full_to_irred, irred_to_full)

      call g_logger%info('generate_reduced_kpoint_mesh: Generated ' // trim(int2str(this%nk_total)) // &
                        ' irreducible k-points from ' // trim(int2str(product(mesh_dims))) // ' total points', &
                        __FILE__, __LINE__)
      call g_logger%info('generate_reduced_kpoint_mesh: Reduction factor: ' // &
                        trim(real2str(real(product(mesh_dims), rp)/real(this%nk_total, rp), '(F6.2)')) // 'x', &
                        __FILE__, __LINE__)
      
      ! Verify weights sum to 1
      if (abs(sum(this%k_weights) - 1.0_rp) > 1.0e-6_rp) then
         if (this%strict_symmetry_checks) then
            call g_logger%fatal('generate_reduced_kpoint_mesh: K-point weights sum to ' // &
                               trim(real2str(sum(this%k_weights), '(F12.8)')) // ' (expected 1.0)', &
                               __FILE__, __LINE__)
         else
            call g_logger%warning('generate_reduced_kpoint_mesh: K-point weights sum to ' // &
                                 trim(real2str(sum(this%k_weights), '(F12.8)')) // ' (should be 1.0)', &
                                 __FILE__, __LINE__)
         end if
      end if
      if (allocated(this%full_to_irred_k)) then
         if (any(this%full_to_irred_k < 1) .or. any(this%full_to_irred_k > this%nk_total)) then
            if (this%strict_symmetry_checks) then
               call g_logger%fatal('generate_reduced_kpoint_mesh: Invalid full_to_irred mapping detected', __FILE__, __LINE__)
            else
               call g_logger%warning('generate_reduced_kpoint_mesh: Invalid full_to_irred mapping detected', __FILE__, __LINE__)
            end if
         end if
      end if
   end subroutine generate_reduced_kpoint_mesh

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Helper subroutine to dump a complex matrix to file for debugging
   !---------------------------------------------------------------------------
   subroutine dump_complex_matrix(matrix, filename, k_point)
      complex(rp), dimension(:, :), intent(in) :: matrix
      character(len=*), intent(in) :: filename
      real(rp), dimension(3), intent(in) :: k_point
      ! Local variables
      integer :: unit, i, j, n_rows, n_cols
      character(len=100) :: fmt_str

      n_rows = size(matrix, 1)
      n_cols = size(matrix, 2)

      open(newunit=unit, file=trim(filename), status='replace', action='write')

      ! Write header
      write(unit, '(A)') '# Complex matrix dump for debugging'
      write(unit, '(A,3F12.8)') '# k-point: ', k_point
      write(unit, '(A,I0,A,I0)') '# Matrix dimensions: ', n_rows, ' x ', n_cols
      write(unit, '(A)') '# Format: row, col, real_part, imag_part'
      write(unit, '(A)') '#'

      ! Write matrix elements
      do i = 1, n_rows
         do j = 1, n_cols
            write(unit, '(2I6,2ES20.12E3)') i, j, real(matrix(i,j)), aimag(matrix(i,j))
         end do
      end do

      close(unit)
   end subroutine dump_complex_matrix

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Main DOS calculation method
   !
   !> @param[in] ham Hamiltonian object
   !> @param[in] n_energy_points Number of energy points (optional)
   !> @param[in] energy_range Energy range [min, max] (optional)
   !> @param[in] method DOS method ('tetrahedron' or 'gaussian')
   !> @param[in] gaussian_sigma Gaussian smearing parameter (optional)
   !> @param[in] temperature Temperature for Fermi-Dirac distribution (optional)
   !> @param[in] fermi_level Fermi level for band moments integration (optional)
   !> @param[in] output_file Output filename (optional)
   !---------------------------------------------------------------------------
   subroutine calculate_density_of_states(this, ham, n_energy_points, energy_range, method, gaussian_sigma, temperature, fermi_level, total_electrons, auto_find_fermi, output_file)
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

      call g_logger%info('calculate_density_of_states: Starting DOS calculation', __FILE__, __LINE__)

      ! Set parameters from optional arguments
      if (present(n_energy_points)) this%n_energy_points = n_energy_points
      if (present(energy_range)) this%dos_energy_range = energy_range
      if (present(method)) this%dos_method = trim(method)
      if (present(gaussian_sigma)) this%gaussian_sigma = gaussian_sigma
      if (present(temperature)) this%temperature = temperature
      if (present(fermi_level)) this%fermi_level = fermi_level
      if (present(total_electrons)) this%total_electrons = total_electrons
      if (present(auto_find_fermi)) this%auto_find_fermi = auto_find_fermi

      ! Set output filename
      filename = 'density_of_states.dat'
      if (present(output_file)) filename = output_file

      ! Setup energy grid
      call this%setup_dos_energy_grid()

      ! Build k-space Hamiltonian and diagonalize if not already done
      if (.not. allocated(this%eigenvalues)) then
         call g_logger%info('calculate_density_of_states: Building and diagonalizing Hamiltonian on k-mesh', __FILE__, __LINE__)
         
         ! Build H(k) for all k-points in mesh
         call this%build_kspace_hamiltonian()
         
         ! Diagonalize to get eigenvalues
         call this%diagonalize_hamiltonian()
      end if

      ! IMPORTANT:
      ! Full-mesh expansion currently expands eigenvalues only. Eigenvectors are
      ! still stored on the irreducible mesh and are required for projected DOS,
      ! band moments, and spin-moment extraction in k-space SCF.
      ! Therefore, do not expand in Gaussian mode to avoid eigenvalue/eigenvector
      ! mesh mismatch in downstream moment/charge routines.
      if (this%use_symmetry_reduction) then
         if (trim(this%dos_method) == 'tetrahedron' .or. trim(this%dos_method) == 'blochl') then
            if (this%nk_total /= this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)) then
               call this%expand_eigenvalues_to_full_mesh()
            end if
         end if
      end if

      ! Calculate DOS based on method
      select case (trim(this%dos_method))
      case ('tetrahedron')
         call g_logger%info('calculate_density_of_states: Using tetrahedron method', __FILE__, __LINE__)
         call this%calculate_dos_tetrahedron()
      case ('blochl')
         call g_logger%info('calculate_density_of_states: Using Blöchl modified tetrahedron method', __FILE__, __LINE__)
         call this%calculate_dos_blochl()
      case ('gaussian')
         call g_logger%info('calculate_density_of_states: Using Gaussian smearing method', __FILE__, __LINE__)
         call this%calculate_dos_gaussian()
      case default
         call g_logger%error('calculate_density_of_states: Unknown DOS method: ' // trim(this%dos_method), __FILE__, __LINE__)
         return
      end select

      ! Calculate orbital projections and band moments
      call this%project_dos_orbitals()
      call this%calculate_band_moments()

      ! Write results to file
      call this%write_dos_to_file(filename)

      call g_logger%info('calculate_density_of_states: DOS calculation completed', __FILE__, __LINE__)
   end subroutine calculate_density_of_states

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Setup energy grid for DOS calculation (all energies in Ry)
   !---------------------------------------------------------------------------
   subroutine setup_dos_energy_grid(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i
      real(rp) :: energy_min, energy_max, delta_energy

      energy_min = this%dos_energy_range(1)
      energy_max = this%dos_energy_range(2)
      delta_energy = (energy_max - energy_min) / real(this%n_energy_points - 1, rp)

      ! Allocate energy grid
      if (allocated(this%dos_energy_grid)) deallocate(this%dos_energy_grid)
      allocate(this%dos_energy_grid(this%n_energy_points))

      ! Fill energy grid in Ry (consistent with Hamiltonian and eigenvalues)
      do i = 1, this%n_energy_points
         this%dos_energy_grid(i) = energy_min + real(i-1, rp) * delta_energy
      end do

      call g_logger%info('setup_dos_energy_grid: Created energy grid with ' // &
                        trim(int2str(this%n_energy_points)) // ' points from ' // &
                        trim(real2str(energy_min, '(F 8.5)')) // ' to ' // trim(real2str(energy_max, '(F 8.5)')) // ' Ry', &
                        __FILE__, __LINE__)
   end subroutine setup_dos_energy_grid

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate DOS using tetrahedron method
   !---------------------------------------------------------------------------
   subroutine calculate_dos_tetrahedron(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i_energy, i_tet, i_corner, i_band, nbands
      integer :: i_start, i_end
      real(rp) :: energy, dos_contrib, de
      real(rp) :: e1, e2, e3, e4, x
      real(rp) :: c1, c2, c3, dos_integral, norm_factor
      real(rp), dimension(4) :: e_corners
      real(rp), allocatable :: sorted_cache(:, :, :)
      real(rp), parameter :: eps = 1.0e-12_rp

      call g_logger%info('calculate_dos_tetrahedron: Calculating DOS using tetrahedron method', __FILE__, __LINE__)

      ! Setup tetrahedra if not already done
      if (.not. allocated(this%tetrahedra)) then
         call this%setup_tetrahedra()
      end if

      ! Allocate DOS arrays
      if (allocated(this%total_dos)) deallocate(this%total_dos)
      allocate(this%total_dos(this%n_energy_points))
      this%total_dos = 0.0_rp
      nbands = size(this%eigenvalues, 1)
      de = (this%dos_energy_grid(this%n_energy_points) - this%dos_energy_grid(1)) / &
           real(this%n_energy_points - 1, rp)
      if (de <= 0.0_rp) then
         call g_logger%fatal('calculate_dos_tetrahedron: Invalid DOS energy grid spacing', __FILE__, __LINE__)
      end if

      ! Precompute sorted tetra-corner energies once per (tetra, band).
      allocate(sorted_cache(4, nbands, this%n_tetrahedra))
      do i_tet = 1, this%n_tetrahedra
         do i_band = 1, nbands
            do i_corner = 1, 4
               e_corners(i_corner) = this%eigenvalues(i_band, this%tetrahedra(i_corner, i_tet))
            end do
            call sort4(e_corners, sorted_cache(:, i_band, i_tet))
         end do
      end do

      ! LMTO-style traversal: loop tetrahedron/band first, update only touched energy bins.
      do i_tet = 1, this%n_tetrahedra
         do i_band = 1, nbands
            e1 = sorted_cache(1, i_band, i_tet)
            e2 = sorted_cache(2, i_band, i_tet)
            e3 = sorted_cache(3, i_band, i_tet)
            e4 = sorted_cache(4, i_band, i_tet)

            if (e4 <= this%dos_energy_grid(1) .or. e1 >= this%dos_energy_grid(this%n_energy_points)) cycle
            if (abs(e2 - e1) < eps .or. abs(e3 - e1) < eps .or. abs(e4 - e1) < eps .or. &
                abs(e3 - e2) < eps .or. abs(e4 - e2) < eps .or. abs(e4 - e3) < eps) cycle

            i_start = max(1, int(floor((e1 - this%dos_energy_grid(1)) / de)) + 1)
            i_end = min(this%n_energy_points, int(ceiling((e4 - this%dos_energy_grid(1)) / de)) + 1)

            c3 = 3.0_rp / ((e2 - e1) * (e3 - e1) * (e4 - e1))
            c2 = 6.0_rp / ((e3 - e1) * (e4 - e1))
            c1 = c2 * (e2 - e1) * 0.5_rp

            do i_energy = i_start, i_end
               energy = this%dos_energy_grid(i_energy)
               if (energy < e1 .or. energy >= e4) cycle
               if (energy < e2) then
                  x = energy - e1
                  dos_contrib = c3 * x * x
               else if (energy < e3) then
                  x = energy - e2
                  dos_contrib = c1 + x * (c2 + x * (3.0_rp * (e1 + e2 - e3 - e4) / &
                                ((e3 - e1) * (e4 - e1) * (e3 - e2) * (e4 - e2))))
               else
                  x = energy - e4
                  dos_contrib = 3.0_rp * x * x / ((e4 - e3) * (e4 - e2) * (e4 - e1))
               end if
               this%total_dos(i_energy) = this%total_dos(i_energy) + dos_contrib * this%tetrahedron_volumes(i_tet)
            end do
         end do
      end do

      dos_integral = trapezoidal_integral(this%dos_energy_grid, this%total_dos)
      if (abs(dos_integral) > 1.0e-12_rp) then
         norm_factor = real(nbands, rp) / dos_integral
         this%total_dos = this%total_dos * norm_factor
         call g_logger%info('calculate_dos_tetrahedron: DOS normalization factor = ' // &
                           trim(real2str(norm_factor, '(ES12.4)')), __FILE__, __LINE__)
      else
         call g_logger%warning('calculate_dos_tetrahedron: DOS integral is ~0, skipping normalization', __FILE__, __LINE__)
      end if

      deallocate(sorted_cache)

      call g_logger%info('calculate_dos_tetrahedron: Tetrahedron DOS calculation completed', __FILE__, __LINE__)
   end subroutine calculate_dos_tetrahedron

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate DOS contribution from a single tetrahedron
   !> Uses linear tetrahedron method (Blöchl)
   !---------------------------------------------------------------------------
   function tetrahedron_dos_contribution(this, energy, e_sorted) result(dos)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: energy
      real(rp), dimension(4), intent(in) :: e_sorted
      real(rp) :: dos

      ! Local variables
      real(rp) :: e1, e2, e3, e4, vol_factor
      real(rp), parameter :: eps = 1.0e-12_rp

      ! Tetrahedron volume factor (adjusted for correct normalization)
      vol_factor = 0.4595_rp

      e1 = e_sorted(1)
      e2 = e_sorted(2)
      e3 = e_sorted(3)
      e4 = e_sorted(4)

      dos = 0.0_rp

      ! Handle different energy ranges
      if (energy < e1 - eps) then
         ! Energy below all eigenvalues
         dos = 0.0_rp
      else if (energy >= e1 - eps .and. energy < e2 - eps) then
         ! Energy in [e1, e2)
         dos = vol_factor * 3.0_rp * (energy - e1)**2 / ((e2 - e1) * (e3 - e1) * (e4 - e1))
      else if (energy >= e2 - eps .and. energy < e3 - eps) then
         ! Energy in [e2, e3)
         dos = vol_factor * (3.0_rp * (energy - e1)**2 - 6.0_rp * (energy - e1) * (energy - e2) + &
                           3.0_rp * (e3 - e1) * (energy - e2) * 2.0_rp / (e3 - e2) + &
                           3.0_rp * (energy - e1) * (energy - e2) * 2.0_rp / (e4 - e2)) / &
                           ((e3 - e1) * (e4 - e1))
      else if (energy >= e3 - eps .and. energy < e4 - eps) then
         ! Energy in [e3, e4)
         dos = vol_factor * 3.0_rp * (e4 - energy)**2 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
      else
         ! Energy above all eigenvalues
         dos = 0.0_rp
      end if

   end function tetrahedron_dos_contribution

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate DOS contribution from a single tetrahedron using Blöchl method
   !> Uses Blöchl improved tetrahedron method (PRB 49, 16223 (1994))
   !---------------------------------------------------------------------------
   function blochl_dos_contribution(this, energy, e_sorted) result(dos)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: energy
      real(rp), dimension(4), intent(in) :: e_sorted
      real(rp) :: dos

      ! Local variables
      real(rp) :: e1, e2, e3, e4, C
      real(rp), parameter :: TOL = 1.0e-10_rp

      e1 = e_sorted(1)
      e2 = e_sorted(2)
      e3 = e_sorted(3)
      e4 = e_sorted(4)

      dos = 0.0_rp

      ! Skip if degenerate (would cause division by zero)
      if (abs(e2-e1) < TOL .or. abs(e3-e1) < TOL .or. abs(e4-e1) < TOL .or. &
          abs(e4-e2) < TOL .or. abs(e4-e3) < TOL .or. abs(e3-e2) < TOL) then
         return
      end if

      ! Blöchl modified tetrahedron DOS contribution
      if (energy <= e1) then
         dos = 0.0_rp
      else if (energy >= e4) then
         dos = 0.0_rp
      else if (energy <= e2) then
         ! Region I: e1 < E <= e2
         ! Blöchl Eq. (23): D(E) = 3(E-e1)²/[(e4-e1)(e3-e1)(e2-e1)]
         dos = 3.0_rp * (energy - e1)**2 / ((e4 - e1) * (e3 - e1) * (e2 - e1))
      else if (energy <= e3) then
         ! Region II: e2 < E <= e3
         ! Blöchl Eq. (24): D(E) = 1/[(e4-e1)(e3-e1)] * 
         !   [3(e2-e1) + 6(E-e2) - 3(e3+e4-e1-e2)(E-e2)²/[(e3-e2)(e4-e2)]]
         C = 1.0_rp / ((e4 - e1) * (e3 - e1))
         dos = C * (3.0_rp * (e2 - e1) + 6.0_rp * (energy - e2) - &
                   3.0_rp * (e3 + e4 - e1 - e2) * (energy - e2)**2 / &
                   ((e3 - e2) * (e4 - e2)))
      else
         ! Region III: e3 < E < e4
         ! Blöchl Eq. (25): D(E) = 3(e4-E)²/[(e4-e1)(e4-e2)(e4-e3)]
         dos = 3.0_rp * (e4 - energy)**2 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
      end if

   end function blochl_dos_contribution

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Sort real array and return sorted values and indices
   !---------------------------------------------------------------------------
   subroutine sort_real_array(arr, sorted, indices)
      real(rp), dimension(:), intent(in) :: arr
      real(rp), dimension(:), intent(out) :: sorted
      integer, dimension(:), intent(out) :: indices

      ! Local variables
      integer :: i, j, n, temp_idx
      real(rp) :: temp_val

      n = size(arr)
      sorted = arr
      do i = 1, n
         indices(i) = i
      end do

      ! Simple bubble sort
      do i = 1, n-1
         do j = 1, n-i
            if (sorted(j) > sorted(j+1)) then
               temp_val = sorted(j)
               sorted(j) = sorted(j+1)
               sorted(j+1) = temp_val
               temp_idx = indices(j)
               indices(j) = indices(j+1)
               indices(j+1) = temp_idx
            end if
         end do
      end do
   end subroutine sort_real_array

   !---------------------------------------------------------------------------
   ! Small optimized sorter for 4 elements (inlined compare-swap network)
   ! This replaces bubble sort for the tetrahedron eigenvalue sorting (size=4)
   subroutine sort4(arr_in, arr_out)
      real(rp), dimension(4), intent(in) :: arr_in
      real(rp), dimension(4), intent(out) :: arr_out
      real(rp) :: a1, a2, a3, a4, tmp

      a1 = arr_in(1)
      a2 = arr_in(2)
      a3 = arr_in(3)
      a4 = arr_in(4)

      ! compare-swap sequence
      if (a1 > a2) then
         tmp = a1
         a1 = a2
         a2 = tmp
      end if
      if (a3 > a4) then
         tmp = a3
         a3 = a4
         a4 = tmp
      end if
      if (a1 > a3) then
         tmp = a1
         a1 = a3
         a3 = tmp
      end if
      if (a2 > a4) then
         tmp = a2
         a2 = a4
         a4 = tmp
      end if
      if (a2 > a3) then
         tmp = a2
         a2 = a3
         a3 = tmp
      end if

      arr_out(1) = a1
      arr_out(2) = a2
      arr_out(3) = a3
      arr_out(4) = a4
   end subroutine sort4

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate DOS using Gaussian smearing
   !---------------------------------------------------------------------------
subroutine calculate_dos_gaussian(this)
   class(reciprocal), intent(inout) :: this

   integer :: i_energy, i_k, i_band
   real(rp) :: energy, weight, gaussian_factor
   real(rp) :: sigma_squared, sigma_use
   real(rp) :: local_sum, dos_integral, norm_factor, kweight_sum
   integer :: nbands

   ! Debug: Check eigenvalue range
   if (.not. this%suppress_internal_logs) then
      call g_logger%info('calculate_dos_gaussian: Eigenvalue range = [' // &
         trim(real2str(minval(this%eigenvalues), '(F12.6)')) // ', ' // &
         trim(real2str(maxval(this%eigenvalues), '(F12.6)')) // '] Ry', __FILE__, __LINE__)
      call g_logger%info('calculate_dos_gaussian: DOS energy range = [' // &
         trim(real2str(this%dos_energy_range(1), '(F12.6)')) // ', ' // &
         trim(real2str(this%dos_energy_range(2), '(F12.6)')) // '] Ry', __FILE__, __LINE__)
      call g_logger%info('calculate_dos_gaussian: Number of k-points = ' // &
         trim(int2str(this%nk_total)), __FILE__, __LINE__)
      call g_logger%info('calculate_dos_gaussian: K-point weight = ' // &
         trim(real2str(this%k_weights(1), '(ES12.4)')), __FILE__, __LINE__)
   end if

   ! Determine sigma (already in Ry from input)
   if (this%gaussian_sigma < 0.001_rp) then
      sigma_use = this%calculate_adaptive_sigma()
      call g_logger%info('calculate_dos_gaussian: Using adaptive sigma = ' // &
                        trim(real2str(sigma_use, '(F8.5)')) // ' Ry', __FILE__, __LINE__)
   else
      sigma_use = this%gaussian_sigma
      call g_logger%info('calculate_dos_gaussian: Using input sigma = ' // &
                        trim(real2str(sigma_use, '(F8.5)')) // ' Ry', __FILE__, __LINE__)
   end if

   ! Allocate DOS arrays
   if (allocated(this%total_dos)) deallocate(this%total_dos)
   allocate(this%total_dos(this%n_energy_points))
   this%total_dos = 0.0_rp

   sigma_squared = sigma_use**2
   nbands = size(this%eigenvalues, 1)
   
   ! DEBUG: Check k-point weights sum
   kweight_sum = sum(this%k_weights)
   
   call g_logger%info('calculate_dos_gaussian: nbands = ' // trim(int2str(nbands)) // &
                     ', nk_total = ' // trim(int2str(this%nk_total)) // &
                     ', k_weights sum = ' // trim(real2str(kweight_sum, '(F12.8)')), __FILE__, __LINE__)
   
   ! DEBUG: Check eigenvalue array size
   call g_logger%info('calculate_dos_gaussian: eigenvalues array size = ' // &
                     trim(int2str(size(this%eigenvalues, 1))) // ' x ' // &
                     trim(int2str(size(this%eigenvalues, 2))), __FILE__, __LINE__)

   ! Calculate raw DOS
   do i_energy = 1, this%n_energy_points
      energy = this%dos_energy_grid(i_energy)  ! Already in Ry
      local_sum = 0.0_rp

!$OMP PARALLEL DO PRIVATE(i_k, i_band, weight, gaussian_factor) REDUCTION(+:local_sum) &
!$OMP& SCHEDULE(STATIC) IF(this%nk_total > 100)
      do i_k = 1, this%nk_total
         weight = this%k_weights(i_k)
         do i_band = 1, nbands
            gaussian_factor = exp(-((energy - this%eigenvalues(i_band, i_k))**2) / (2.0_rp * sigma_squared))
            gaussian_factor = gaussian_factor / (sigma_use * sqrt(2.0_rp * 3.141592653589793_rp))
            local_sum = local_sum + weight * gaussian_factor
         end do
      end do
!$OMP END PARALLEL DO
      this%total_dos(i_energy) = local_sum
   end do

   ! DEBUG: Check DOS integral WITHOUT normalization first
   dos_integral = 0.0_rp
   do i_energy = 1, this%n_energy_points - 1
      dos_integral = dos_integral + 0.5_rp * (this%total_dos(i_energy) + this%total_dos(i_energy+1)) * &
                                   (this%dos_energy_grid(i_energy+1) - this%dos_energy_grid(i_energy))
   end do
   
   call g_logger%info('calculate_dos_gaussian: Raw DOS (before norm) integrates to ' // &
                     trim(real2str(dos_integral, '(F12.6)')) // ' (should be ' // trim(int2str(nbands)) // ')', &
                     __FILE__, __LINE__)
   
   ! Normalize DOS to integrate to nbands (total number of states)
   if (abs(dos_integral) > 1.0e-10_rp) then
      norm_factor = 1.0_rp ! real(nbands, rp) / dos_integral
      this%total_dos = this%total_dos * norm_factor
      
      call g_logger%info('calculate_dos_gaussian: DOS normalized by factor ' // &
                        trim(real2str(norm_factor, '(F10.6)')), __FILE__, __LINE__)
      call g_logger%info('calculate_dos_gaussian: DOS integrates to ' // &
                        trim(real2str(dos_integral * norm_factor, '(F10.4)')) // &
                        ' (should be ' // trim(int2str(nbands)) // ')', __FILE__, __LINE__)
   else
      call g_logger%error('calculate_dos_gaussian: DOS integral is zero!', __FILE__, __LINE__)
   end if

   call g_logger%info('calculate_dos_gaussian: Gaussian DOS calculation completed', __FILE__, __LINE__)
end subroutine calculate_dos_gaussian

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate DOS using Blöchl modified tetrahedron method (PRB 49, 16223 (1994))
   !> Provides improved convergence over standard tetrahedron method
   !> All energies in Rydberg
   !---------------------------------------------------------------------------
   subroutine calculate_dos_blochl(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i_energy, i_tet, i_corner, i_band, i_start, i_end
      real(rp) :: energy, dos_contrib, de
      real(rp), dimension(4) :: e_corners
      integer :: nbands
      real(rp) :: e1, e2, e3, e4, E, C
      real(rp) :: dos_integral, norm_factor
      integer :: skipped_degenerate
      real(rp), parameter :: TOL = 1.0e-10_rp
      real(rp), allocatable :: sorted_cache(:, :, :)
      logical, allocatable :: degenerate_cache(:, :)
      real(rp), allocatable :: e1_cache(:, :), e2_cache(:, :), e3_cache(:, :), e4_cache(:, :)

      call g_logger%info('calculate_dos_blochl: Calculating DOS using Blöchl modified tetrahedron method', __FILE__, __LINE__)

      ! Setup tetrahedra if not already done
      if (.not. allocated(this%tetrahedra)) then
         call this%setup_tetrahedra()
      end if

      ! Debug: Check eigenvalue range
      if (.not. this%suppress_internal_logs) then
         call g_logger%info('calculate_dos_blochl: Eigenvalue range = [' // &
            trim(real2str(minval(this%eigenvalues), '(F12.6)')) // ', ' // &
            trim(real2str(maxval(this%eigenvalues), '(F12.6)')) // '] Ry', __FILE__, __LINE__)
         call g_logger%info('calculate_dos_blochl: DOS energy range = [' // &
            trim(real2str(this%dos_energy_range(1), '(F12.6)')) // ', ' // &
            trim(real2str(this%dos_energy_range(2), '(F12.6)')) // '] Ry', __FILE__, __LINE__)
         call g_logger%info('calculate_dos_blochl: Number of tetrahedra = ' // &
            trim(int2str(this%n_tetrahedra)), __FILE__, __LINE__)
         call g_logger%info('calculate_dos_blochl: Tetrahedron volume = ' // &
            trim(real2str(this%tetrahedron_volumes(1), '(ES12.4)')), __FILE__, __LINE__)
      end if

      ! Allocate DOS arrays
      if (allocated(this%total_dos)) deallocate(this%total_dos)
      allocate(this%total_dos(this%n_energy_points))
      this%total_dos = 0.0_rp
      de = (this%dos_energy_grid(this%n_energy_points) - this%dos_energy_grid(1)) / &
           real(this%n_energy_points - 1, rp)
      if (de <= 0.0_rp) then
         call g_logger%fatal('calculate_dos_blochl: Invalid DOS energy grid spacing', __FILE__, __LINE__)
      end if

      nbands = size(this%eigenvalues, 1)
      allocate(sorted_cache(4, nbands, this%n_tetrahedra))
      allocate(degenerate_cache(nbands, this%n_tetrahedra))
      allocate(e1_cache(nbands, this%n_tetrahedra), e2_cache(nbands, this%n_tetrahedra))
      allocate(e3_cache(nbands, this%n_tetrahedra), e4_cache(nbands, this%n_tetrahedra))

      ! Precompute sorted tetra-corner energies and degeneracy masks once.
      do i_tet = 1, this%n_tetrahedra
         do i_band = 1, nbands
            do i_corner = 1, 4
               e_corners(i_corner) = this%eigenvalues(i_band, this%tetrahedra(i_corner, i_tet))
            end do
            call sort4(e_corners, sorted_cache(:, i_band, i_tet))
            e1 = sorted_cache(1, i_band, i_tet)
            e2 = sorted_cache(2, i_band, i_tet)
            e3 = sorted_cache(3, i_band, i_tet)
            e4 = sorted_cache(4, i_band, i_tet)
            e1_cache(i_band, i_tet) = e1
            e2_cache(i_band, i_tet) = e2
            e3_cache(i_band, i_tet) = e3
            e4_cache(i_band, i_tet) = e4
            degenerate_cache(i_band, i_tet) = &
               (abs(e2-e1) < TOL .or. abs(e3-e1) < TOL .or. abs(e4-e1) < TOL .or. &
                abs(e4-e2) < TOL .or. abs(e4-e3) < TOL .or. abs(e3-e2) < TOL)
         end do
      end do

      skipped_degenerate = 0
      ! Bounded energy-window traversal (tetra/band outer loop).
      do i_tet = 1, this%n_tetrahedra
         do i_band = 1, nbands
            e1 = e1_cache(i_band, i_tet)
            e2 = e2_cache(i_band, i_tet)
            e3 = e3_cache(i_band, i_tet)
            e4 = e4_cache(i_band, i_tet)
            if (e4 <= this%dos_energy_grid(1) .or. e1 >= this%dos_energy_grid(this%n_energy_points)) cycle
            if (degenerate_cache(i_band, i_tet)) then
               skipped_degenerate = skipped_degenerate + 1
               cycle
            end if

            i_start = max(1, int(floor((e1 - this%dos_energy_grid(1)) / de)) + 1)
            i_end = min(this%n_energy_points, int(ceiling((e4 - this%dos_energy_grid(1)) / de)) + 1)

            do i_energy = i_start, i_end
               E = this%dos_energy_grid(i_energy)
               if (E <= e1 .or. E >= e4) cycle
               if (E <= e2) then
                  dos_contrib = 3.0_rp * (E - e1)**2 / ((e4 - e1) * (e3 - e1) * (e2 - e1))
               else if (E <= e3) then
                  C = 1.0_rp / ((e4 - e1) * (e3 - e1))
                  dos_contrib = C * (3.0_rp * (e2 - e1) + 6.0_rp * (E - e2) - &
                                    3.0_rp * (e3 + e4 - e1 - e2) * (E - e2)**2 / &
                                    ((e3 - e2) * (e4 - e2)))
               else
                  dos_contrib = 3.0_rp * (e4 - E)**2 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
               end if
               this%total_dos(i_energy) = this%total_dos(i_energy) + dos_contrib * this%tetrahedron_volumes(i_tet)
            end do
         end do
      end do

      if (skipped_degenerate > 0) then
         call g_logger%info('calculate_dos_blochl: Skipped ' // trim(int2str(skipped_degenerate)) // &
                           ' degenerate tetrahedron contributions', __FILE__, __LINE__)
      end if

      ! Normalize to get correct number of states
      dos_integral = 0.0_rp
      do i_energy = 1, this%n_energy_points - 1
         dos_integral = dos_integral + 0.5_rp * (this%total_dos(i_energy) + this%total_dos(i_energy+1)) * &
                                      (this%dos_energy_grid(i_energy+1) - this%dos_energy_grid(i_energy))
      end do

      if (abs(dos_integral) > TOL) then
         norm_factor = real(nbands, rp) / dos_integral
         this%total_dos = this%total_dos * norm_factor
         
         call g_logger%info('calculate_dos_blochl: DOS normalized by factor ' // &
                           trim(real2str(norm_factor, '(F10.6)')), __FILE__, __LINE__)
         call g_logger%info('calculate_dos_blochl: DOS integrates to ' // &
                           trim(real2str(dos_integral * norm_factor, '(F10.4)')) // &
                           ' (should be ' // trim(int2str(nbands)) // ')', __FILE__, __LINE__)
      else
         call g_logger%warning('calculate_dos_blochl: DOS integral is zero', __FILE__, __LINE__)
      end if
      deallocate(degenerate_cache, sorted_cache, e1_cache, e2_cache, e3_cache, e4_cache)

      call g_logger%info('calculate_dos_blochl: Blöchl DOS calculation completed', __FILE__, __LINE__)
   end subroutine calculate_dos_blochl


!---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Setup tetrahedra for tetrahedron DOS method
   !---------------------------------------------------------------------------
   subroutine setup_tetrahedra(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: nk1, nk2, nk3, n_tet_per_cube, i, j, k, tet_idx
      integer :: n1, n2, n3, idx
      integer, dimension(4, 6) :: tetrahedra_corners
      integer, dimension(4) :: corner_indices

      call g_logger%info('setup_tetrahedra: Setting up tetrahedra for Brillouin zone integration', __FILE__, __LINE__)

      ! Get mesh dimensions
      nk1 = this%nk_mesh(1)
      nk2 = this%nk_mesh(2)
      nk3 = this%nk_mesh(3)

      ! Number of tetrahedra per cube (6 for standard decomposition)
      n_tet_per_cube = 6

      ! Total number of tetrahedra
      this%n_tetrahedra = n_tet_per_cube * nk1 * nk2 * nk3

      ! Allocate tetrahedra array
      if (allocated(this%tetrahedra)) deallocate(this%tetrahedra)
      allocate(this%tetrahedra(4, this%n_tetrahedra))

      ! Allocate tetrahedron volumes (all equal for uniform mesh)
      if (allocated(this%tetrahedron_volumes)) deallocate(this%tetrahedron_volumes)
      allocate(this%tetrahedron_volumes(this%n_tetrahedra))

      ! Volume of each tetrahedron (1/6 of cube volume, times number of cubes)
      this%tetrahedron_volumes = 1.0_rp / (6.0_rp * real(nk1 * nk2 * nk3, rp))

      ! Standard tetrahedron decomposition for a cube
      ! Each tetrahedron is defined by 4 corner indices relative to cube
      tetrahedra_corners(:, 1) = [1, 2, 4, 5]  ! Tetrahedron 1
      tetrahedra_corners(:, 2) = [2, 3, 4, 5]  ! Tetrahedron 2
      tetrahedra_corners(:, 3) = [3, 4, 5, 6]  ! Tetrahedron 3
      tetrahedra_corners(:, 4) = [4, 5, 7, 8]  ! Tetrahedron 4
      tetrahedra_corners(:, 5) = [5, 6, 7, 8]  ! Tetrahedron 5
      tetrahedra_corners(:, 6) = [4, 5, 6, 7]  ! Tetrahedron 6

      ! Build tetrahedra
      tet_idx = 0
      do i = 1, nk1
         do j = 1, nk2
            do k = 1, nk3
               ! For each cube in the mesh
               do n1 = 1, n_tet_per_cube
                  tet_idx = tet_idx + 1

                  ! Get corner indices for this tetrahedron
                  corner_indices = tetrahedra_corners(:, n1)

                  ! Convert relative corner indices to absolute k-point indices
                  do n2 = 1, 4
                     n3 = corner_indices(n2)
                     ! Convert corner number to i,j,k offsets
                     select case (n3)
                     case (1)
                        idx = this%get_kpoint_index(i, j, k, nk1, nk2, nk3)
                     case (2)
                        idx = this%get_kpoint_index(i+1, j, k, nk1, nk2, nk3)
                     case (3)
                        idx = this%get_kpoint_index(i+1, j+1, k, nk1, nk2, nk3)
                     case (4)
                        idx = this%get_kpoint_index(i, j+1, k, nk1, nk2, nk3)
                     case (5)
                        idx = this%get_kpoint_index(i, j, k+1, nk1, nk2, nk3)
                     case (6)
                        idx = this%get_kpoint_index(i+1, j, k+1, nk1, nk2, nk3)
                     case (7)
                        idx = this%get_kpoint_index(i, j+1, k+1, nk1, nk2, nk3)
                     case (8)
                        idx = this%get_kpoint_index(i+1, j+1, k+1, nk1, nk2, nk3)
                     end select
                     this%tetrahedra(n2, tet_idx) = idx
                  end do
               end do
            end do
         end do
      end do

      call g_logger%info('setup_tetrahedra: Created ' // trim(int2str(this%n_tetrahedra)) // &
                        ' tetrahedra from ' // trim(int2str(nk1)) // 'x' // trim(int2str(nk2)) // &
                        'x' // trim(int2str(nk3)) // ' k-mesh', __FILE__, __LINE__)
   end subroutine setup_tetrahedra

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Expand eigenvalues from irreducible k-points to full mesh using symmetry
   !> This allows tetrahedron method to work with symmetry-reduced calculations
   !---------------------------------------------------------------------------
   subroutine expand_eigenvalues_to_full_mesh(this)
      class(reciprocal), intent(inout) :: this
      integer :: ik_full, nk_full, nk_irred, nbands, ik_irred
      real(rp), dimension(:, :), allocatable :: eigenvalues_full
      complex(rp), dimension(:, :, :), allocatable :: eigenvectors_full
      integer :: nrow_evec, nband_evec, nk_evec
      
      if (.not. this%use_symmetry_reduction) then
         call g_logger%info('expand_eigenvalues_to_full_mesh: No symmetry reduction, skipping', __FILE__, __LINE__)
         return
      end if
      
      nk_irred = this%nk_total  ! Current number of irreducible k-points
      nk_full = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)  ! Full mesh size
      nbands = size(this%eigenvalues, 1)
      
      call g_logger%info('expand_eigenvalues_to_full_mesh: Expanding ' // trim(int2str(nk_irred)) // &
                        ' irreducible k-points to ' // trim(int2str(nk_full)) // ' full mesh points', __FILE__, __LINE__)
      
      ! Allocate full mesh eigenvalues
      allocate(eigenvalues_full(nbands, nk_full))
      eigenvalues_full = 0.0_rp
      
      if (.not. allocated(this%full_to_irred_k)) then
         if (this%strict_symmetry_checks) then
            call g_logger%fatal('expand_eigenvalues_to_full_mesh: full_to_irred mapping not available', __FILE__, __LINE__)
         else
            call g_logger%warning('expand_eigenvalues_to_full_mesh: full_to_irred mapping not available; fallback copy', __FILE__, __LINE__)
            do ik_irred = 1, min(nk_irred, nk_full)
               eigenvalues_full(:, ik_irred) = this%eigenvalues(:, ik_irred)
            end do
         end if
      else
         do ik_full = 1, nk_full
            ik_irred = this%full_to_irred_k(ik_full)
            if (ik_irred >= 1 .and. ik_irred <= nk_irred) then
               eigenvalues_full(:, ik_full) = this%eigenvalues(:, ik_irred)
            else
               if (this%strict_symmetry_checks) then
                  call g_logger%fatal('expand_eigenvalues_to_full_mesh: invalid mapping index ' // trim(int2str(ik_irred)), __FILE__, __LINE__)
               end if
            end if
         end do
      end if
      
      ! Store expanded eigenvalues
      deallocate(this%eigenvalues)
      allocate(this%eigenvalues(nbands, nk_full))
      this%eigenvalues = eigenvalues_full

      ! If eigenvectors are available on irreducible mesh, expand them with
      ! the same full->irred mapping so projection/moment integrations remain
      ! consistent with expanded eigenvalues.
      if (allocated(this%eigenvectors)) then
         nrow_evec = size(this%eigenvectors, 1)
         nband_evec = size(this%eigenvectors, 2)
         nk_evec = size(this%eigenvectors, 3)
         if (nband_evec == nbands .and. nk_evec == nk_irred) then
            allocate(eigenvectors_full(nrow_evec, nband_evec, nk_full))
            eigenvectors_full = cmplx(0.0_rp, 0.0_rp, rp)
            if (allocated(this%full_to_irred_k)) then
               do ik_full = 1, nk_full
                  ik_irred = this%full_to_irred_k(ik_full)
                  if (ik_irred >= 1 .and. ik_irred <= nk_irred) then
                     eigenvectors_full(:, :, ik_full) = this%eigenvectors(:, :, ik_irred)
                  end if
               end do
            end if
            deallocate(this%eigenvectors)
            allocate(this%eigenvectors(nrow_evec, nband_evec, nk_full))
            this%eigenvectors = eigenvectors_full
            deallocate(eigenvectors_full)
         else
            call g_logger%warning('expand_eigenvalues_to_full_mesh: eigenvector dimensions do not match irreducible eigenvalue mesh; skipping eigenvector expansion', __FILE__, __LINE__)
         end if
      end if
      
      ! Update nk_total to full mesh
      this%nk_total = nk_full
      ! Keep k-points/weights/maps consistent with the expanded full mesh.
      call this%generate_mp_mesh()
      
      deallocate(eigenvalues_full)
      
      call g_logger%info('expand_eigenvalues_to_full_mesh: Expansion complete using explicit full_to_irred map', __FILE__, __LINE__)
   end subroutine expand_eigenvalues_to_full_mesh

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate DOS using tetrahedron method with optional symmetry support
   !---------------------------------------------------------------------------
   subroutine calculate_dos_tetrahedron_with_symmetry(this)
      class(reciprocal), intent(inout) :: this
      
      ! If using symmetry reduction, expand eigenvalues to full mesh first
      if (this%use_symmetry_reduction) then
         call g_logger%info('calculate_dos_tetrahedron_with_symmetry: Expanding eigenvalues for symmetry-reduced mesh', __FILE__, __LINE__)
         call this%expand_eigenvalues_to_full_mesh()
      end if
      
      ! Now call standard tetrahedron method
      call this%calculate_dos_tetrahedron()
   end subroutine calculate_dos_tetrahedron_with_symmetry

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Get k-point index from i,j,k coordinates (with periodic boundary conditions)
   !---------------------------------------------------------------------------
   function get_kpoint_index(this, i, j, k, nk1, nk2, nk3) result(idx)
      class(reciprocal), intent(in) :: this
      integer, intent(in) :: i, j, k, nk1, nk2, nk3
      integer :: idx, ii, jj, kk

      ! Apply periodic boundary conditions
      ii = mod(i-1, nk1) + 1
      jj = mod(j-1, nk2) + 1
      kk = mod(k-1, nk3) + 1

      ! Convert to 1D index (assuming k-points are stored as k1 varying fastest)
      idx = ii + (jj-1)*nk1 + (kk-1)*nk1*nk2
   end function get_kpoint_index

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate adaptive Gaussian sigma based on k-mesh density
   !> Uses rule: sigma ≈ α × ΔE_avg where ΔE_avg is typical band spacing
   !---------------------------------------------------------------------------
   function calculate_adaptive_sigma(this) result(sigma)
      class(reciprocal), intent(in) :: this
      real(rp) :: sigma
      real(rp) :: bz_volume, k_density, typical_spacing
      integer :: nk_total
      
      ! Calculate Brillouin zone k-point density
      nk_total = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)
      
      if (nk_total > 0) then
         ! BZ volume per k-point (in reciprocal space units)
         bz_volume = this%reciprocal_volume
         k_density = real(nk_total, rp) / bz_volume
         
         ! Typical energy spacing scales as 1/k_density^(1/3)
         ! For metallic systems: ΔE ∝ E_F / N_k^(1/3)
         ! Use heuristic: sigma = C / nk^(1/3) where C is tuned for accuracy
      ! Use a smaller scale factor (1.0 eV) for typical spacing to avoid huge sigma
      typical_spacing = 1.0_rp / (real(nk_total, rp)**(1.0_rp/3.0_rp))
         
         ! Adaptive sigma: smaller for denser meshes, larger for coarse meshes
         ! Factor of 0.5-1.0 gives good balance between accuracy and smoothness
         sigma = 0.7_rp * typical_spacing
         
         ! Clamp to reasonable range
         ! Clamp sigma to a conservative range (5 meV .. 50 meV) to avoid over-broadening
         sigma = max(0.005_rp, min(sigma, 0.05_rp))
      else
         sigma = 0.1_rp  ! Default fallback
      end if
      
   call g_logger%info('calculate_adaptive_sigma: Adaptive sigma = ' // trim(real2str(sigma, '(F8.5)')) // &
            ' Ry for ' // trim(int2str(nk_total)) // ' k-points', __FILE__, __LINE__)
   end function calculate_adaptive_sigma

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Project DOS onto orbitals, sites, and spin
   !---------------------------------------------------------------------------
   subroutine project_dos_orbitals(this)
      class(reciprocal), intent(inout) :: this

      call g_logger%info('project_dos_orbitals: Starting orbital projection calculation', __FILE__, __LINE__)

      ! Use tetrahedron or Gaussian method based on dos_method
      if (trim(this%dos_method) == 'tetrahedron' .or. trim(this%dos_method) == 'blochl') then
         call this%project_dos_orbitals_tetrahedron()
      else
         call this%project_dos_orbitals_gaussian()
      end if
   end subroutine project_dos_orbitals

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Project DOS onto orbitals, sites, and spin using Gaussian method
   !---------------------------------------------------------------------------
subroutine project_dos_orbitals_gaussian(this)
   class(reciprocal), intent(inout) :: this
   integer :: ik, ib, ie, iorb, ispin, i, isite
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

   call g_logger%info('project_dos_orbitals_gaussian: Starting projection', __FILE__, __LINE__)

   ! Determine sigma (same as in calculate_dos_gaussian, already in Ry)
   if (this%gaussian_sigma < 0.001_rp) then
      sigma_use = this%calculate_adaptive_sigma()
   else
      sigma_use = this%gaussian_sigma
   end if
   sigma_squared = sigma_use**2

   this%n_sites = this%lattice%nrec
   this%n_orb_types = 4
   this%n_spin_components = 2
   nbands = size(this%eigenvalues, 1)

   if (this%n_sites <= 0) then
      call g_logger%fatal('project_dos_orbitals_gaussian: invalid number of sites', __FILE__, __LINE__)
   end if
   if (size(this%eigenvectors, 1) <= 0) then
      call g_logger%fatal('project_dos_orbitals_gaussian: eigenvectors not available', __FILE__, __LINE__)
   end if
   if (mod(size(this%eigenvectors, 1), this%n_sites) /= 0) then
      call g_logger%fatal('project_dos_orbitals_gaussian: eigenvector basis size not divisible by number of sites', __FILE__, __LINE__)
   end if
   n_orb_site = size(this%eigenvectors, 1)/this%n_sites
   if (mod(n_orb_site, this%n_spin_components) /= 0) then
      call g_logger%fatal('project_dos_orbitals_gaussian: per-site basis size incompatible with spin components', __FILE__, __LINE__)
   end if
   n_orb_per_spin = n_orb_site/this%n_spin_components

   if (allocated(this%projected_dos)) deallocate(this%projected_dos)
   allocate(this%projected_dos(this%n_sites, this%n_orb_types, this%n_spin_components, this%n_energy_points))
   this%projected_dos = 0.0_rp

   allocate(site_orb_offset(this%n_sites + 1))
   do isite = 1, this%n_sites + 1
      site_orb_offset(isite) = (isite - 1) * n_orb_site
   end do
   lstart = [1, 2, 5, 10]
   lend = [1, 4, 9, 16]

   ! Diagnostic logging
   call g_logger%info('project_dos_orbitals_gaussian: n_sites = ' // trim(int2str(this%n_sites)) // &
                     ', nbands = ' // trim(int2str(nbands)) // &
                     ', max_orb_channels = ' // trim(int2str(this%max_orbs)) // &
                     ', eigenvector size = ' // trim(int2str(size(this%eigenvectors, 1))), __FILE__, __LINE__)
   ! call g_logger%info('project_dos_orbitals_gaussian: site_orb_offset = [' // &
   !                   trim(int2str(site_orb_offset(1))) // ', ' // &
   !                   trim(int2str(site_orb_offset(2))) // ', ' // &
   !                   trim(int2str(site_orb_offset(3))) // ', ' // &
   !                   trim(int2str(site_orb_offset(4))) // ', ' // &
   !                   trim(int2str(site_orb_offset(5))) // ']', __FILE__, __LINE__)

   ! Calculate projected DOS (raw) - same weights as total DOS
   do ie = 1, this%n_energy_points
      energy = this%dos_energy_grid(ie)  ! Already in Ry

      do ik = 1, this%nk_total
         do ib = 1, nbands  ! Loop over all bands, not just max_orb_channels
            ! Skip if eigenvalue is far from current energy
            if (abs(this%eigenvalues(ib, ik) - energy) > 5.0_rp * sigma_use) cycle

            ! Calculate Gaussian weight (already normalized)
            gaussian_weight = exp(-((energy - this%eigenvalues(ib, ik))**2) / (2.0_rp * sigma_squared))
            gaussian_weight = gaussian_weight / (sigma_use * sqrt(2.0_rp * 3.141592653589793_rp))
            
            if (abs(gaussian_weight) < 1.0e-10_rp) cycle

            ! Apply k-point weight
            weight = gaussian_weight * this%k_weights(ik)

            do isite = 1, this%n_sites
               ! Site-blocked layout: eigenvectors are [site1(18), site2(18), ...]
               ! where each site block contains [orb1_up...orb9_up, orb1_dn...orb9_dn]
               site_orb_start = site_orb_offset(isite)

               do ispin = 1, this%n_spin_components
                  orb_start = site_orb_start + (ispin - 1) * n_orb_per_spin
                  do iorb = 1, 4
                     orbital_char = 0.0_rp
                     if (lstart(iorb) <= n_orb_per_spin) then
                        do i = lstart(iorb), min(lend(iorb), n_orb_per_spin)
                           if (orb_start + i <= size(this%eigenvectors, 1)) then
                              psi_element = this%eigenvectors(orb_start + i, ib, ik)
                              orbital_char = orbital_char + real(conjg(psi_element) * psi_element, rp)
                           end if
                        end do
                     end if
                     this%projected_dos(isite, iorb, ispin, ie) = &
                        this%projected_dos(isite, iorb, ispin, ie) + orbital_char * weight
                  end do
               end do
            end do
         end do
      end do
   end do

   ! Normalize projected DOS so integrated sum over projections equals nbands
   call g_logger%info('project_dos_orbitals_gaussian: Projection completed (raw)', __FILE__, __LINE__)
   if (this%n_energy_points > 2) then
      proj_integral = 0.0_rp
      do iei = 1, this%n_energy_points - 1
         e_low = this%dos_energy_grid(iei)
         e_high = this%dos_energy_grid(iei+1)
         proj_integral = proj_integral + 0.5_rp * ( sum(this%projected_dos(:, :, :, iei)) + sum(this%projected_dos(:, :, :, iei+1)) ) * (e_high - e_low)
      end do

      if (abs(proj_integral) > 1.0e-12_rp) then
         norm_factor = 1.0_rp ! real(nbands, rp) / proj_integral
         this%projected_dos = this%projected_dos * norm_factor
         call g_logger%info('project_dos_orbitals_gaussian: Normalized projected DOS by factor ' // trim(real2str(norm_factor, '(F10.6)')), __FILE__, __LINE__)
      else
         call g_logger%warning('project_dos_orbitals_gaussian: projected DOS integral is zero, skipping normalization', __FILE__, __LINE__)
      end if

      ! Diagnostic: check mid-energy ratio after normalization
      ie = this%n_energy_points / 2
      if (abs(this%total_dos(ie)) > 1.0e-12_rp) then
         call g_logger%info('project_dos_orbitals_gaussian: At mid-energy, proj/total ratio (post-norm) = ' // &
            trim(real2str(sum(this%projected_dos(:, :, :, ie)) / this%total_dos(ie), '(F10.6)')), __FILE__, __LINE__)
      end if
   end if

   deallocate(site_orb_offset)
end subroutine project_dos_orbitals_gaussian
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Project DOS onto orbitals, sites, and spin using tetrahedron method
   !---------------------------------------------------------------------------
   subroutine project_dos_orbitals_tetrahedron(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i_energy, i_tet, i_corner, i_band, iorb, ispin, i, isite
      integer :: i_start, i_end
      integer :: n_orb_per_spin, orb_start, site_orb_start, ik, n_orb_site
      integer :: ie, nbands
      real(rp) :: energy, dos_contrib, orbital_char_avg, orbital_char, de
      real(rp) :: total_dos_integral, proj_dos_integral, norm_factor
      real(rp), dimension(4) :: e_corners, sorted_e, orbital_chars
      integer :: lstart(4), lend(4)
      complex(rp) :: psi_element
      real(rp) :: e1, e2, e3, e4, x, C
      real(rp), parameter :: TOL = 1.0e-10_rp
      ! Per-site orbital offsets for mixed atom types
      integer, dimension(:), allocatable :: site_orb_offset

      call g_logger%info('project_dos_orbitals_tetrahedron: Starting tetrahedron orbital projection calculation', __FILE__, __LINE__)

      ! Initialize dimensions - get number of sites from lattice
      this%n_sites = this%lattice%nrec
      this%n_orb_types = 4  ! s, p, d, f
      this%n_spin_components = 2  ! spin up/down
      if (this%n_sites <= 0) then
         call g_logger%fatal('project_dos_orbitals_tetrahedron: invalid number of sites', __FILE__, __LINE__)
      end if
      if (size(this%eigenvectors, 1) <= 0) then
         call g_logger%fatal('project_dos_orbitals_tetrahedron: eigenvectors not available', __FILE__, __LINE__)
      end if
      if (mod(size(this%eigenvectors, 1), this%n_sites) /= 0) then
         call g_logger%fatal('project_dos_orbitals_tetrahedron: eigenvector basis size not divisible by number of sites', __FILE__, __LINE__)
      end if
      n_orb_site = size(this%eigenvectors, 1)/this%n_sites
      if (mod(n_orb_site, this%n_spin_components) /= 0) then
         call g_logger%fatal('project_dos_orbitals_tetrahedron: per-site basis size incompatible with spin components', __FILE__, __LINE__)
      end if
      n_orb_per_spin = n_orb_site/this%n_spin_components
      lstart = [1, 2, 5, 10]
      lend = [1, 4, 9, 16]

      call g_logger%info('project_dos_orbitals_tetrahedron: Projecting onto ' // trim(int2str(this%n_sites)) // &
                        ' site(s)', __FILE__, __LINE__)

      ! Allocate projected DOS array
      if (allocated(this%projected_dos)) deallocate(this%projected_dos)
      allocate(this%projected_dos(this%n_sites, this%n_orb_types, this%n_spin_components, this%n_energy_points))
      this%projected_dos = 0.0_rp

      allocate(site_orb_offset(this%n_sites + 1))
      do isite = 1, this%n_sites + 1
         site_orb_offset(isite) = (isite - 1) * n_orb_site
      end do

      ! Setup tetrahedra if not already done
      if (.not. allocated(this%tetrahedra)) then
         call this%setup_tetrahedra()
      end if
      de = (this%dos_energy_grid(this%n_energy_points) - this%dos_energy_grid(1)) / &
           real(this%n_energy_points - 1, rp)
      if (de <= 0.0_rp) then
         call g_logger%fatal('project_dos_orbitals_tetrahedron: Invalid DOS energy grid spacing', __FILE__, __LINE__)
      end if

      ! Tetra/band-first traversal with bounded energy windows.
      do i_tet = 1, this%n_tetrahedra
         do i_band = 1, size(this%eigenvalues, 1)
            do i_corner = 1, 4
               ik = this%tetrahedra(i_corner, i_tet)
               e_corners(i_corner) = this%eigenvalues(i_band, ik)
            end do
            call sort4(e_corners, sorted_e)
            e1 = sorted_e(1)
            e2 = sorted_e(2)
            e3 = sorted_e(3)
            e4 = sorted_e(4)
            if (e4 <= this%dos_energy_grid(1) .or. e1 >= this%dos_energy_grid(this%n_energy_points)) cycle
            if (abs(e2 - e1) < TOL .or. abs(e3 - e1) < TOL .or. abs(e4 - e1) < TOL .or. &
                abs(e3 - e2) < TOL .or. abs(e4 - e2) < TOL .or. abs(e4 - e3) < TOL) cycle

            i_start = max(1, int(floor((e1 - this%dos_energy_grid(1)) / de)) + 1)
            i_end = min(this%n_energy_points, int(ceiling((e4 - this%dos_energy_grid(1)) / de)) + 1)

            do isite = 1, this%n_sites
               site_orb_start = site_orb_offset(isite)
               do ispin = 1, this%n_spin_components
                  orb_start = site_orb_start + (ispin - 1) * n_orb_per_spin
                  do iorb = 1, 4
                     orbital_char_avg = 0.0_rp
                     if (lstart(iorb) <= n_orb_per_spin) then
                        do i_corner = 1, 4
                           ik = this%tetrahedra(i_corner, i_tet)
                           orbital_char = 0.0_rp
                           do i = lstart(iorb), min(lend(iorb), n_orb_per_spin)
                              if (orb_start + i <= size(this%eigenvectors, 1)) then
                                 psi_element = this%eigenvectors(orb_start + i, i_band, ik)
                                 orbital_char = orbital_char + real(conjg(psi_element) * psi_element, rp)
                              end if
                           end do
                           orbital_chars(i_corner) = orbital_char
                        end do
                        orbital_char_avg = sum(orbital_chars) / 4.0_rp
                     end if
                     if (abs(orbital_char_avg) < 1.0e-14_rp) cycle

                     do i_energy = i_start, i_end
                        energy = this%dos_energy_grid(i_energy)
                        if (trim(this%dos_method) == 'blochl') then
                           if (energy <= e1 .or. energy >= e4) cycle
                           if (energy <= e2) then
                              dos_contrib = 3.0_rp * (energy - e1)**2 / ((e4 - e1) * (e3 - e1) * (e2 - e1))
                           else if (energy <= e3) then
                              C = 1.0_rp / ((e4 - e1) * (e3 - e1))
                              dos_contrib = C * (3.0_rp * (e2 - e1) + 6.0_rp * (energy - e2) - &
                                             3.0_rp * (e3 + e4 - e1 - e2) * (energy - e2)**2 / &
                                             ((e3 - e2) * (e4 - e2)))
                           else
                              dos_contrib = 3.0_rp * (e4 - energy)**2 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
                           end if
                        else
                           if (energy < e1 .or. energy >= e4) cycle
                           if (energy < e2) then
                              x = energy - e1
                              dos_contrib = 3.0_rp * x * x / ((e2 - e1) * (e3 - e1) * (e4 - e1))
                           else if (energy < e3) then
                              x = energy - e2
                              dos_contrib = 3.0_rp * (e2 - e1) / ((e3 - e1) * (e4 - e1)) + &
                                          x * (6.0_rp / ((e3 - e1) * (e4 - e1)) + &
                                          x * (3.0_rp * (e1 + e2 - e3 - e4) / &
                                          ((e3 - e1) * (e4 - e1) * (e3 - e2) * (e4 - e2))))
                           else
                              x = energy - e4
                              dos_contrib = 3.0_rp * x * x / ((e4 - e3) * (e4 - e2) * (e4 - e1))
                           end if
                        end if
                        dos_contrib = dos_contrib * this%tetrahedron_volumes(i_tet)
                        this%projected_dos(isite, iorb, ispin, i_energy) = this%projected_dos(isite, iorb, ispin, i_energy) + &
                                                                           orbital_char_avg * dos_contrib
                     end do
                  end do
               end do
            end do
         end do
      end do

      deallocate(site_orb_offset)

      ! Normalize projected DOS to match total DOS normalization
      ! When using Blöchl or standard tetrahedron method, the raw DOS needs normalization
      ! to integrate to the correct number of states. We must apply the same normalization
      ! to projected DOS for consistency.
      if (allocated(this%total_dos) .and. this%n_energy_points > 2) then
         nbands = size(this%eigenvalues, 1)
         total_dos_integral = 0.0_rp
         do ie = 1, this%n_energy_points - 1
            total_dos_integral = total_dos_integral + &
               0.5_rp * (this%total_dos(ie) + this%total_dos(ie+1)) * &
               (this%dos_energy_grid(ie+1) - this%dos_energy_grid(ie))
         end do
         
         ! The total_dos is already normalized, so total_dos_integral ≈ nbands
         ! Calculate integral of projected DOS (raw, unnormalized)
         proj_dos_integral = 0.0_rp
         do ie = 1, this%n_energy_points - 1
            ! Sum over all sites, orbitals, and spins
            proj_dos_integral = proj_dos_integral + &
               0.5_rp * (sum(this%projected_dos(:,:,:,ie)) + sum(this%projected_dos(:,:,:,ie+1))) * &
               (this%dos_energy_grid(ie+1) - this%dos_energy_grid(ie))
         end do
         
         ! Normalize projected DOS by the same factor
         if (abs(proj_dos_integral) > 1.0e-10_rp) then
            norm_factor = total_dos_integral / proj_dos_integral
            this%projected_dos = this%projected_dos * norm_factor
            call g_logger%info('project_dos_orbitals_tetrahedron: Normalized projected DOS by factor ' // &
                              trim(real2str(norm_factor, '(F10.6)')), __FILE__, __LINE__)
         end if
      end if

      call g_logger%info('project_dos_orbitals_tetrahedron: Tetrahedron orbital projection calculation completed', __FILE__, __LINE__)
   end subroutine project_dos_orbitals_tetrahedron

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate trapezoidal integral of y(x) over x grid
   !---------------------------------------------------------------------------
   function trapezoidal_integral(x, y) result(integral)
      real(rp), dimension(:), intent(in) :: x, y
      real(rp) :: integral

      ! Local variables
      integer :: n, i
      real(rp) :: dx

      n = size(x)
      if (n /= size(y)) then
         call g_logger%error('trapezoidal_integral: x and y arrays must have same size', __FILE__, __LINE__)
         integral = 0.0_rp
         return
      end if

      integral = 0.0_rp

      ! Trapezoidal rule: ∫ y dx ≈ Σ (y_i + y_{i+1}) * (x_{i+1} - x_i) / 2
      do i = 1, n-1
         dx = x(i+1) - x(i)
         integral = integral + 0.5_rp * (y(i) + y(i+1)) * dx
      end do
   end function trapezoidal_integral

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate band moments from projected DOS
   !---------------------------------------------------------------------------
subroutine calculate_band_moments(this)
   class(reciprocal), intent(inout) :: this

   integer :: isite, iorb, ispin, ie, n_energy
   real(rp) :: energy, dos_value, fermi_weight
   real(rp) :: m0, m1, m2
   real(rp), dimension(:), allocatable :: integrand, fermi_dist
   real(rp) :: kT, fermi_arg
   real(rp) :: total_occupation, expected_electrons, total_occupation_alt
      real(rp), allocatable :: energy_grid_ry(:)
      real(rp), parameter :: eV_to_Ry = 0.073498618_rp

   call g_logger%info('calculate_band_moments: Starting calculation', __FILE__, __LINE__)
   call g_logger%info('calculate_band_moments: DEBUG - this%auto_find_fermi = ' // &
                     merge('TRUE ', 'FALSE', this%auto_find_fermi) // &
                     ', this%total_electrons = ' // trim(real2str(this%total_electrons, '(F10.5)')), &
                     __FILE__, __LINE__)
   call g_logger%info('calculate_band_moments: DEBUG - Fermi level on entry = ' // &
                     trim(real2str(this%fermi_level, '(F10.6)')) // ' Ry', __FILE__, __LINE__)

   ! Auto-find Fermi level if requested
   if (this%auto_find_fermi .and. this%total_electrons > 0.0_rp) then
      this%fermi_level = this%find_fermi_level_from_dos(this%total_electrons)
      call g_logger%info('calculate_band_moments: Auto-found Fermi level = ' // &
                        trim(real2str(this%fermi_level, '(F 8.5)')) // ' Ry', __FILE__, __LINE__)
   else
      call g_logger%info('calculate_band_moments: Using pre-set Fermi level = ' // &
                        trim(real2str(this%fermi_level, '(F 8.5)')) // ' Ry (auto_find disabled)', &
                        __FILE__, __LINE__)
   end if

   if (allocated(this%band_moments)) deallocate(this%band_moments)
   allocate(this%band_moments(this%n_sites, this%n_orb_types, this%n_spin_components, 3))
   this%band_moments = 0.0_rp

   n_energy = this%n_energy_points
   allocate(integrand(n_energy))
   allocate(fermi_dist(n_energy))
   allocate(energy_grid_ry(n_energy))
   energy_grid_ry = this%dos_energy_grid  ! Already in Ry

   ! Boltzmann constant: kB = 6.3336814e-6 Ry/K
   kT = this%temperature * 6.3336814e-6_rp  ! Ry/K
   
   call g_logger%info('calculate_band_moments: kT = ' // trim(real2str(kT, '(ES12.5)')) // &
                     ' Ry at T = ' // trim(real2str(this%temperature, '(F8.2)')) // ' K', &
                     __FILE__, __LINE__)

   ! Pre-calculate Fermi-Dirac distribution
   do ie = 1, n_energy
      energy = energy_grid_ry(ie)
      fermi_arg = (energy - this%fermi_level) / kT
      
      if (fermi_arg > 50.0_rp) then
         fermi_dist(ie) = 0.0_rp
      else if (fermi_arg < -50.0_rp) then
         fermi_dist(ie) = 1.0_rp
      else
         fermi_dist(ie) = 1.0_rp / (exp(fermi_arg) + 1.0_rp)
      end if
   end do

   ! Calculate total occupation to verify normalization
   integrand = this%total_dos * fermi_dist
   total_occupation = trapezoidal_integral(energy_grid_ry, integrand)
   
   ! Also calculate using the same method as find_fermi_level_from_dos for comparison
   total_occupation_alt = this%integrate_dos_up_to_energy(this%fermi_level, kT)
   
   call g_logger%info('calculate_band_moments: Total occupation (method 1) = ' // &
                     trim(real2str(total_occupation, '(F10.5)')) // &
                     ', (method 2) = ' // trim(real2str(total_occupation_alt, '(F10.5)')) // &
                     ', expected = ' // trim(real2str(this%total_electrons, '(F10.5)')), &
                     __FILE__, __LINE__)

   ! Calculate moments for each projection
   do isite = 1, this%n_sites
      do iorb = 1, this%n_orb_types
         do ispin = 1, this%n_spin_components

            ! m0: occupation = ∫ DOS(E) * f(E) dE
            integrand = this%projected_dos(isite, iorb, ispin, :) * fermi_dist
            m0 = trapezoidal_integral(energy_grid_ry, integrand)

            ! m1: band center = ∫ E * DOS(E) * f(E) dE / m0
            if (abs(m0) > 1.0e-12_rp) then
               integrand = energy_grid_ry * this%projected_dos(isite, iorb, ispin, :) * fermi_dist
               m1 = trapezoidal_integral(energy_grid_ry, integrand) / m0
            else
               m1 = 0.0_rp
            end if

            ! m2: band width = sqrt(∫ (E - m1)² * DOS(E) * f(E) dE / m0)
            if (abs(m0) > 1.0e-12_rp) then
               ! Use energy grid in Ry for the variance integral (units must match projected_dos which
               ! was integrated/normalized in Ry). Previously we accidentally passed the eV grid here.
               integrand = (energy_grid_ry - m1)**2 * &
                          this%projected_dos(isite, iorb, ispin, :) * fermi_dist
               m2 = trapezoidal_integral(energy_grid_ry, integrand) / m0
               m2 = sqrt(max(m2, 0.0_rp))
            else
               m2 = 0.0_rp
            end if

            this%band_moments(isite, iorb, ispin, 1) = m0
            this%band_moments(isite, iorb, ispin, 2) = m1
            this%band_moments(isite, iorb, ispin, 3) = m2
         end do
      end do
   end do

   ! Calculate sum of all m0 moments (should equal total occupation)
   m0 = 0.0_rp
   do isite = 1, this%n_sites
      do iorb = 1, this%n_orb_types
         do ispin = 1, this%n_spin_components
            m0 = m0 + this%band_moments(isite, iorb, ispin, 1)
         end do
      end do
   end do
   
   call g_logger%info('calculate_band_moments: Sum of all m0 moments = ' // &
                     trim(real2str(m0, '(F10.5)')) // &
                     ', total occupation = ' // trim(real2str(total_occupation, '(F10.5)')), &
                     __FILE__, __LINE__)

   deallocate(integrand, fermi_dist)

   call g_logger%info('calculate_band_moments: Completed', __FILE__, __LINE__)
end subroutine calculate_band_moments

   !--------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print total and spin-polarized electron counts (DOS integrated) using
   !> the current DOS and Fermi level. Intended to be called each SCF
   !> iteration for diagnostics.
   !--------------------------------------------------------------------------
   subroutine print_total_and_spin_dos(this)
      class(reciprocal), intent(in) :: this

      integer :: ie, n_energy, isite, iorb, ispin
      real(rp), allocatable :: fermi_dist(:)
      real(rp), allocatable :: energy_grid(:)
      real(rp), allocatable :: integrand_up(:), integrand_dn(:)
      real(rp) :: kT, fermi_arg
      real(rp) :: electrons_up, electrons_dn, electrons_total
      real(rp), parameter :: kB_Ry_per_K = 6.3336814e-6_rp

      if (.not. allocated(this%total_dos) .or. .not. allocated(this%dos_energy_grid)) then
         call g_logger%warning('print_total_and_spin_dos: DOS not available to print', __FILE__, __LINE__)
         return
      end if

      n_energy = this%n_energy_points
      allocate(fermi_dist(n_energy))
      allocate(energy_grid(n_energy))
      energy_grid = this%dos_energy_grid

      ! Boltzmann factor (Ry)
      kT = this%temperature * kB_Ry_per_K

      ! Pre-calc fermi distribution on energy grid
      do ie = 1, n_energy
         fermi_arg = (energy_grid(ie) - this%fermi_level)
         if (kT > 1.0e-10_rp) then
            fermi_dist(ie) = 1.0_rp / (exp(fermi_arg / kT) + 1.0_rp)
         else
            if (energy_grid(ie) <= this%fermi_level) then
               fermi_dist(ie) = 1.0_rp
            else
               fermi_dist(ie) = 0.0_rp
            end if
         end if
      end do

      ! Initialize integrands
      allocate(integrand_up(n_energy))
      allocate(integrand_dn(n_energy))
      integrand_up = 0.0_rp
      integrand_dn = 0.0_rp

      if (allocated(this%projected_dos)) then
         ! Sum projected DOS across sites and orbitals for each spin
         do isite = 1, this%n_sites
            do iorb = 1, this%n_orb_types
               if (this%n_spin_components >= 1) then
                  integrand_up = integrand_up + this%projected_dos(isite, iorb, 1, :)
               end if
               if (this%n_spin_components >= 2) then
                  integrand_dn = integrand_dn + this%projected_dos(isite, iorb, 2, :)
               end if
            end do
         end do
      else
         ! No projected DOS available: fall back to total DOS split by spin if possible
         if (this%n_spin_components == 2) then
            ! If no projection but spin resolved, attempt to split total DOS equally
            integrand_up = 0.5_rp * this%total_dos
            integrand_dn = 0.5_rp * this%total_dos
         else
            integrand_up = this%total_dos
            integrand_dn = 0.0_rp
         end if
      end if

      ! Apply Fermi weighting and integrate
      integrand_up = integrand_up * fermi_dist
      integrand_dn = integrand_dn * fermi_dist

      electrons_up = trapezoidal_integral(energy_grid, integrand_up)
      electrons_dn = trapezoidal_integral(energy_grid, integrand_dn)
      electrons_total = electrons_up + electrons_dn

      ! Log results
      if (this%n_spin_components >= 2) then
         call g_logger%info('DOS summary: Total electrons (from DOS) = ' // trim(real2str(electrons_total, '(F10.6)')) // &
                           ', Up = ' // trim(real2str(electrons_up, '(F10.6)')) // &
                           ', Down = ' // trim(real2str(electrons_dn, '(F10.6)')) // &
                           ', M = ' // trim(real2str(electrons_up - electrons_dn, '(F10.6)')), __FILE__, __LINE__)
      else
         call g_logger%info('DOS summary: Total electrons (from DOS) = ' // trim(real2str(electrons_total, '(F10.6)')), __FILE__, __LINE__)
      end if

      deallocate(fermi_dist, energy_grid, integrand_up, integrand_dn)
   end subroutine print_total_and_spin_dos

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Find Fermi level from calculated DOS by integrating to find electron count
   !---------------------------------------------------------------------------
function find_fermi_level_from_dos(this, total_electrons) result(fermi_level)
   class(reciprocal), intent(in) :: this
   real(rp), intent(in) :: total_electrons
   real(rp) :: fermi_level

   integer :: ie, max_iter
   real(rp) :: energy, kT
   real(rp) :: e_min, e_max, e_mid, electrons_at_e
   real(rp), parameter :: eV_to_Ry = 0.073498618_rp
   real(rp), parameter :: kB_Ry_per_K = 6.3336814e-6_rp

   call g_logger%info('find_fermi_level_from_dos: Finding Fermi level for ' // &
                     trim(real2str(total_electrons, '(F 8.5)')) // ' electrons at T = ' // &
                     trim(real2str(this%temperature, '(F 8.5)')) // ' K', __FILE__, __LINE__)

   if (.not. allocated(this%total_dos)) then
      call g_logger%error('find_fermi_level_from_dos: Total DOS not calculated', __FILE__, __LINE__)
      fermi_level = 0.0_rp
      return
   end if

   ! Boltzmann constant in Ry/K
   kT = this%temperature * kB_Ry_per_K

   ! Energy range in Ry (dos_energy_grid is already in Ry, no conversion needed)
   e_min = this%dos_energy_grid(1)
   e_max = this%dos_energy_grid(this%n_energy_points)
   max_iter = 100

   ! Bisection method
   do ie = 1, max_iter
      e_mid = (e_min + e_max) / 2.0_rp
      electrons_at_e = this%integrate_dos_up_to_energy(e_mid, kT)

      ! ! DEBUG: Print first few and last iterations
      ! if (ie <= 5 .or. ie >= max_iter-2) then
      !    call g_logger%info('  Bisection iter ' // trim(int2str(ie)) // ': E=' // &
      !                      trim(real2str(e_mid, '(F10.6)')) // ' Ry, electrons=' // &
      !                      trim(real2str(electrons_at_e, '(F12.8)')) // ', target=' // &
      !                      trim(real2str(total_electrons, '(F12.8)')), __FILE__, __LINE__)
      ! end if

      if (abs(electrons_at_e - total_electrons) < 1.0e-6_rp) then
         fermi_level = e_mid
         call g_logger%info('  Bisection converged at iteration ' // trim(int2str(ie)), __FILE__, __LINE__)
         exit
      else if (electrons_at_e < total_electrons) then
         e_min = e_mid
      else
         e_max = e_mid
      end if
   end do

   if (ie >= max_iter) then
      call g_logger%warning('  Bisection did NOT converge after ' // trim(int2str(max_iter)) // ' iterations!', __FILE__, __LINE__)
   end if

   fermi_level = e_mid

   ! Final check
   electrons_at_e = this%integrate_dos_up_to_energy(fermi_level, kT)
   call g_logger%info('find_fermi_level_from_dos: Found Fermi level at ' // &
                     trim(real2str(fermi_level, '(F 8.5)')) // ' Ry (integrated ' // &
                     trim(real2str(electrons_at_e, '(F 8.5)')) // ' electrons)', __FILE__, __LINE__)
   
   ! DEBUG: Check total DOS integral
   call g_logger%info('find_fermi_level_from_dos: DEBUG - Checking DOS integral...', __FILE__, __LINE__)
   call g_logger%info('  Total electrons requested: ' // trim(real2str(total_electrons, '(F10.6)')), __FILE__, __LINE__)
   call g_logger%info('  Total electrons found: ' // trim(real2str(electrons_at_e, '(F10.6)')), __FILE__, __LINE__)
   call g_logger%info('  Error: ' // trim(real2str(electrons_at_e - total_electrons, '(F10.6)')) // &
                     ' (' // trim(real2str(100.0_rp*(electrons_at_e - total_electrons)/total_electrons, '(F6.3)')) // '%)', &
                     __FILE__, __LINE__)
end function find_fermi_level_from_dos


!---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Integrate DOS up to given energy with Fermi-Dirac weighting
   !---------------------------------------------------------------------------
function integrate_dos_up_to_energy(this, energy, kT) result(integral)
   class(reciprocal), intent(in) :: this
   real(rp), intent(in) :: energy, kT
   real(rp) :: integral

   integer :: ie
   real(rp) :: e, fermi_weight, delta_e
   real(rp), parameter :: eV_to_Ry = 0.073498618_rp

   integral = 0.0_rp

   do ie = 1, this%n_energy_points - 1
      ! dos_energy_grid is already in Ry, no conversion needed
      e = this%dos_energy_grid(ie)
      delta_e = this%dos_energy_grid(ie+1) - this%dos_energy_grid(ie)

      ! Fermi-Dirac weight at current energy
      if (kT > 1.0e-10_rp) then
         fermi_weight = 1.0_rp / (exp((e - energy) / kT) + 1.0_rp)
      else
         ! T=0 limit
         if (e <= energy) then
            fermi_weight = 1.0_rp
         else
            fermi_weight = 0.0_rp
         end if
      end if

      ! Trapezoidal integration
      integral = integral + 0.5_rp * delta_e * (this%total_dos(ie) * fermi_weight + &
                                               this%total_dos(ie+1) * (1.0_rp / (exp((this%dos_energy_grid(ie+1) - energy) / kT) + 1.0_rp)))
   end do
end function integrate_dos_up_to_energy

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Write DOS results to file
   !---------------------------------------------------------------------------
   subroutine write_dos_to_file(this, filename)
      class(reciprocal), intent(in) :: this
      character(len=*), intent(in) :: filename

      ! Local variables
      integer :: unit, i_energy, isite, iorb, ispin
      character(len=256) :: proj_filename
      real(rp), parameter :: eV_to_Ry = 0.073498618_rp

      call g_logger%info('write_dos_to_file: Writing DOS to ' // trim(filename), __FILE__, __LINE__)

      ! Write total DOS
      open(newunit=unit, file=trim(filename), status='replace', action='write')

      ! Write header
      write(unit, '(A)') '# Density of States'
      write(unit, '(A,A)') '# Method: ', trim(this%dos_method)
      if (trim(this%dos_method) == 'gaussian') then
         write(unit, '(A,F8.5,A)') '# Gaussian sigma: ', this%gaussian_sigma, ' Ry'
      end if
      write(unit, '(A,I0)') '# Energy points: ', this%n_energy_points
      write(unit, '(A,2F10.6,A)') '# Energy range: ', this%dos_energy_range(1), &
                                   this%dos_energy_range(2), ' Ry'
      write(unit, '(A)') '# Energy (Ry)    Total DOS'

      ! Write DOS data (energy grid already in Ry)
      do i_energy = 1, this%n_energy_points
         write(unit, '(2F12.6)') this%dos_energy_grid(i_energy) - this%fermi_level, this%total_dos(i_energy)
      end do

      close(unit)

      ! Write projected DOS if available
      if (allocated(this%projected_dos)) then
         proj_filename = 'projected_dos.dat'
         call g_logger%info('write_dos_to_file: Writing projected DOS to ' // trim(proj_filename), __FILE__, __LINE__)

         open(newunit=unit, file=trim(proj_filename), status='replace', action='write')

         ! Write header
         write(unit, '(A)') '# Projected Density of States'
         write(unit, '(A,A)') '# Method: ', trim(this%dos_method)
         if (trim(this%dos_method) == 'gaussian') then
            write(unit, '(A,F8.5,A)') '# Gaussian sigma: ', this%gaussian_sigma, ' Ry'
         end if
         write(unit, '(A,I0)') '# Energy points: ', this%n_energy_points
         write(unit, '(A,2F10.6,A)') '# Energy range: ', this%dos_energy_range(1), &
                                      this%dos_energy_range(2), ' Ry'
         write(unit, '(A)') '# Columns: Energy(Ry), s_up, p_up, d_up, f_up, s_down, p_down, d_down, f_down'

         ! Write projected DOS data (energy grid already in Ry)
         do i_energy = 1, this%n_energy_points
            write(unit, '(9F12.6)') this%dos_energy_grid(i_energy) - this%fermi_level, &
                                  (this%projected_dos(1, iorb, 1, i_energy), iorb=1,4), &  ! spin up: s,p,d,f
                                  (this%projected_dos(1, iorb, 2, i_energy), iorb=1,4)     ! spin down: s,p,d,f
         end do

         close(unit)
         call g_logger%info('write_dos_to_file: Projected DOS written to file', __FILE__, __LINE__)
      end if

      ! Write band moments if available
      if (allocated(this%band_moments)) then
         proj_filename = 'band_moments.dat'
         call g_logger%info('write_dos_to_file: Writing band moments to ' // trim(proj_filename), __FILE__, __LINE__)

         open(newunit=unit, file=trim(proj_filename), status='replace', action='write')

         ! Write header
         write(unit, '(A)') '# Band Moments'
         write(unit, '(A,A)') '# Method: ', trim(this%dos_method)
         if (trim(this%dos_method) == 'gaussian') then
            write(unit, '(A,F8.5,A)') '# Gaussian sigma: ', this%gaussian_sigma, ' Ry'
         end if
         write(unit, '(A,F8.4,A)') '# Temperature: ', this%temperature, ' K'
         write(unit, '(A,F12.6,A)') '# Fermi level: ', this%fermi_level, ' Ry'
         write(unit, '(A)') '# Format: site, orbital_type, spin, m0, m1, m2'
         write(unit, '(A)') '# Orbital types: 1=s, 2=p, 3=d, 4=f'
         write(unit, '(A)') '# Spin: 1=up, 2=down'
         write(unit, '(A)') '# m0 = integrated DOS up to Fermi level with Fermi-Dirac weighting'
         write(unit, '(A)') '# m1 = center of gravity (occupied states)'
         write(unit, '(A)') '# m2 = width of occupied states'

         ! Write band moments data
         do isite = 1, this%n_sites
            do iorb = 1, this%n_orb_types
               do ispin = 1, this%n_spin_components
                  write(unit, '(3I3,3F12.6)') isite, iorb, ispin, &
                                             this%band_moments(isite, iorb, ispin, 1), &  ! m0
                                             this%band_moments(isite, iorb, ispin, 2), &  ! m1
                                             this%band_moments(isite, iorb, ispin, 3)     ! m2
               end do
            end do
         end do

         close(unit)
         call g_logger%info('write_dos_to_file: Band moments written to file', __FILE__, __LINE__)
      end if

      call g_logger%info('write_dos_to_file: DOS written to file', __FILE__, __LINE__)
   end subroutine write_dos_to_file

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate band energy from first moment of DOS
   !> Band energy = sum over sites, orbitals, spins of (m0 * m1)
   !---------------------------------------------------------------------------
   function calculate_band_energy_from_moments(this) result(eband)
      class(reciprocal), intent(in) :: this
      real(rp) :: eband
      integer :: isite, iorb, ispin
      
      eband = 0.0_rp
      
      if (.not. allocated(this%band_moments)) return
      
      ! Band energy = ∫ E * DOS(E) * f(E) dE
      !             = sum over all orbitals of (m0 * m1)
      ! where m0 = occupation, m1 = band center
      do isite = 1, this%n_sites
         do iorb = 1, this%n_orb_types
            do ispin = 1, this%n_spin_components
               eband = eband + this%band_moments(isite, iorb, ispin, 1) * &
                              this%band_moments(isite, iorb, ispin, 2)
            end do
         end do
      end do
      
   end function calculate_band_energy_from_moments

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate Gaussian weight for a single energy pair
   !---------------------------------------------------------------------------
   function calculate_gaussian_weight_single(this, grid_energy, eigenvalue) result(weight)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: grid_energy, eigenvalue
      real(rp) :: weight

      ! Gaussian smearing: weight = exp(-(E_grid - E_eigen)²/(2σ²)) / (σ√(2π))
      real(rp) :: prefactor, exponent, delta_e

      delta_e = grid_energy - eigenvalue
      prefactor = 1.0_rp / (this%gaussian_sigma * sqrt(2.0_rp * 3.141592653589793_rp))
      exponent = -0.5_rp * (delta_e / this%gaussian_sigma)**2

      if (exponent > -20.0_rp) then  ! Avoid underflow
         weight = prefactor * exp(exponent)
      else
         weight = 0.0_rp
      end if
   end function calculate_gaussian_weight_single

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Sort eigenvalues and return sorted values with sort indices
   !> Helper subroutine for tetrahedron methods
   !---------------------------------------------------------------------------
   pure subroutine sort_eigenvalues(e_in, e_sorted, sort_idx)
      real(rp), dimension(4), intent(in) :: e_in
      real(rp), dimension(4), intent(out) :: e_sorted
      integer, dimension(4), intent(out) :: sort_idx
      
      ! Local variables
      integer :: i, j, min_idx
      real(rp) :: temp_e
      integer :: temp_idx
      
      ! Initialize
      e_sorted = e_in
      do i = 1, 4
         sort_idx(i) = i
      end do
      
      ! Simple selection sort
      do i = 1, 3
         min_idx = i
         do j = i + 1, 4
            if (e_sorted(j) < e_sorted(min_idx)) then
               min_idx = j
            end if
         end do
         
         if (min_idx /= i) then
            ! Swap energies
            temp_e = e_sorted(i)
            e_sorted(i) = e_sorted(min_idx)
            e_sorted(min_idx) = temp_e
            
            ! Swap indices
            temp_idx = sort_idx(i)
            sort_idx(i) = sort_idx(min_idx)
            sort_idx(min_idx) = temp_idx
         end if
      end do
   end subroutine sort_eigenvalues

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate Local Density Matrix (LDM) from k-space projected DOS
   !> 
   !> This subroutine calculates the local density matrix needed for LDA+U
   !> calculations by integrating the orbital and spin-resolved projected DOS
   !> up to the Fermi level. The LDM is stored in the symbolic_atom potentials.
   !>
   !> The projected DOS already contains the orbital character information, so
   !> we can directly integrate it to get the occupations. For the off-diagonal
   !> elements, we need to reconstruct them from the eigenvector projections.
   !>
   !> @param[in,out] this - Reciprocal object with projected DOS and eigenvalues
   !> @param[in,out] lattice_obj - Lattice object to store LDM in symbolic_atoms
   !>
   !> Created by Anders Bergman, 2025-10-21, for k-space LDA+U implementation
   !---------------------------------------------------------------------------
   subroutine calculate_ldm_from_projected_dos(this, lattice_obj)
      use lattice_mod
      class(reciprocal), intent(inout) :: this
      type(lattice), intent(inout) :: lattice_obj
      
      ! Local variables
      integer :: isite
      
      ! Orbital indexing: 
      ! iorb = 1 (s, 1 orbital), iorb = 2 (p, 3 orbitals), iorb = 3 (d, 5 orbitals)
      
      call g_logger%info('calculate_ldm_from_projected_dos: Computing LDM for LDA+U from k-space DOS', __FILE__, __LINE__)
      
      ! Initialize all LDM to zero
      do isite = 1, lattice_obj%nrec
         lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%ldm(:,:,:,:) = 0.0_rp
      end do
      
      ! Method 2: Full LDM (including off-diagonal) from eigenvector projections
      ! This properly captures orbital hybridization and is more accurate than diagonal-only
      call calculate_ldm_from_eigenvectors(this, lattice_obj)
      
      ! Save flattened copy for each site
      do isite = 1, lattice_obj%nrec
         call lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%flatten_ldm()
      end do
      
      ! Print LDM for debugging
      call g_logger%info('calculate_ldm_from_projected_dos: LDM calculation complete', __FILE__, __LINE__)
      
   end subroutine calculate_ldm_from_projected_dos

! DESCRIPTION:
   !> @brief
   !> Calculate full LDM including off-diagonal elements from eigenvectors
   !>
   !> This calculates: LDM(site, l, spin, m1, m2) = Σ_{k,n} f_{k,n} ψ*_{m1} ψ_{m2} w_k
   !> where f_{k,n} is the Fermi occupation (0 or 1 for T=0), ψ are eigenvector
   !> projections onto atomic orbitals, and w_k are k-point weights.
   !>
   !> @param[in] this - Reciprocal object with eigenvalues and eigenvectors
   !> @param[in,out] lattice_obj - Lattice object to store LDM
   !---------------------------------------------------------------------------
   subroutine calculate_ldm_from_eigenvectors(this, lattice_obj)
      use lattice_mod
      class(reciprocal), intent(inout) :: this
      type(lattice), intent(inout) :: lattice_obj
      
      ! Local variables
      integer :: ik, ib, isite, alpha, beta, iorb, norb_l
      integer :: ispin, basis_offset, site_offset, n_orb_site, n_orb_per_spin
      real(rp) :: kweight, fermi_occ
      complex(rp), allocatable :: DM(:,:)
      complex(rp) :: psi_a, psi_b
      integer :: nbands, nsites, i
      real(rp) :: trace_tot, trace_up, trace_dn, max_imag
      integer :: idx_a, idx_b

      call g_logger%info('calculate_ldm_from_eigenvectors: Computing site density matrices from eigenvectors', __FILE__, __LINE__)

      nbands = size(this%eigenvalues, 1)
      nsites = lattice_obj%nrec
      if (nsites <= 0) then
         call g_logger%fatal('calculate_ldm_from_eigenvectors: invalid number of sites', __FILE__, __LINE__)
      end if
      if (mod(size(this%eigenvectors, 1), nsites) /= 0) then
         call g_logger%fatal('calculate_ldm_from_eigenvectors: eigenvector basis size not divisible by number of sites', __FILE__, __LINE__)
      end if
      n_orb_site = size(this%eigenvectors, 1)/nsites
      if (mod(n_orb_site, this%n_spin_components) /= 0) then
         call g_logger%fatal('calculate_ldm_from_eigenvectors: per-site basis size incompatible with spin components', __FILE__, __LINE__)
      end if
      n_orb_per_spin = n_orb_site/this%n_spin_components

      ! Allocate temporary DM (n_orb_site x n_orb_site)
      allocate(DM(n_orb_site, n_orb_site))

      ! Zero per-site LDM storage before accumulation
      do isite = 1, nsites
         lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%ldm(:,:,:,:) = 0.0_rp
      end do

      ! Accumulate density matrix per site
      ! Loop order: k-points -> sites -> bands for better cache locality
      do ik = 1, this%nk_total
         kweight = this%k_weights(ik)
         
         do isite = 1, nsites
            ! Reset DM for this site at this k-point
            DM = (0.0_rp, 0.0_rp)
            site_offset = (isite - 1) * n_orb_site
            
            ! Accumulate contributions from all occupied bands at this k-point
            do ib = 1, nbands
               if (this%eigenvalues(ib, ik) > this%fermi_level) cycle
               fermi_occ = 1.0_rp
               
               ! Build density matrix: DM_{ab} += f * w_k * ψ*_a ψ_b
               do alpha = 1, n_orb_site
                  psi_a = this%eigenvectors(site_offset + alpha, ib, ik)
                  do beta = 1, n_orb_site
                     psi_b = this%eigenvectors(site_offset + beta, ib, ik)
                     DM(alpha, beta) = DM(alpha, beta) + conjg(psi_a) * psi_b * fermi_occ * kweight
                  end do
               end do
            end do ! bands
            
            ! Diagnostics: check DM properties for first k-point
            if (ik == 1) then
               trace_tot = 0.0_rp
               trace_up = 0.0_rp
               trace_dn = 0.0_rp
               max_imag = 0.0_rp
               
               do i = 1, n_orb_site
                  trace_tot = trace_tot + real(DM(i, i), rp)
                  ! Check for large imaginary parts (should be ~0 for physical DM)
                  max_imag = max(max_imag, abs(aimag(DM(i, i))))
                  
                  if (i <= n_orb_per_spin) then
                     trace_up = trace_up + real(DM(i, i), rp)
                  else
                     trace_dn = trace_dn + real(DM(i, i), rp)
                  end if
               end do
               
               call g_logger%info('calculate_ldm_from_eigenvectors: DM diagnostics site=' // trim(int2str(isite)) // &
                                 ', trace_tot=' // trim(real2str(trace_tot, '(F10.6)')) // &
                                 ', trace_up=' // trim(real2str(trace_up, '(F10.6)')) // &
                                 ', trace_dn=' // trim(real2str(trace_dn, '(F10.6)')) // &
                                 ', max_imag=' // trim(real2str(max_imag, '(E10.2)')), __FILE__, __LINE__)
               
               if (max_imag > 1.0e-6_rp) then
                  call g_logger%warning('calculate_ldm_from_eigenvectors: Large imaginary part in DM diagonal: ' // &
                                       trim(real2str(max_imag, '(E10.2)')), __FILE__, __LINE__)
               end if
            end if
            
            ! Project DM into ldm structure: iterate over angular momentum channels
            do iorb = 0, min(3, lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%lmax)
               norb_l = 2*iorb + 1
               basis_offset = iorb**2  ! 0, 1, 4 mapping
               
               do ispin = 1, this%n_spin_components
                  do alpha = 1, norb_l
                     do beta = 1, norb_l
                        ! Map to DM indices: spin_block_offset + basis_offset + m
                        idx_a = (ispin - 1) * n_orb_per_spin + basis_offset + alpha
                        idx_b = (ispin - 1) * n_orb_per_spin + basis_offset + beta
                        if (idx_a > n_orb_site .or. idx_b > n_orb_site) cycle
                        
                        ! Accumulate real part (imaginary should be ~0 for proper DM)
                        lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%ldm(iorb+1, ispin, alpha, beta) = &
                           lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)%potential%ldm(iorb+1, ispin, alpha, beta) + &
                           real(DM(idx_a, idx_b), rp)
                     end do
                  end do
               end do
            end do
            
         end do ! sites
      end do ! k-points

      deallocate(DM)
      call g_logger%info('calculate_ldm_from_eigenvectors: Eigenvector-based 18x18 DM -> ldm complete', __FILE__, __LINE__)
      
   end subroutine calculate_ldm_from_eigenvectors
end module reciprocal_mod
