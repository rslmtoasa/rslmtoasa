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
   use precision_mod, only: rp
   use math_mod
   use string_mod, only: int2str, real2str
   use logger_mod, only: g_logger
   use timer_mod, only: g_timer
   use symmetry_mod, only: symmetry
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
      complex(rp), dimension(:, :, :, :), allocatable :: hk_bulk
      !> Spin-orbit Hamiltonian in k-space
      complex(rp), dimension(:, :, :), allocatable :: hk_so
      !> Total k-space Hamiltonian (bulk + SO)
      complex(rp), dimension(:, :, :), allocatable :: hk_total
      !> Overlap matrix in k-space
      complex(rp), dimension(:, :, :), allocatable :: sk_overlap

      ! Band structure variables
      !> Maximum number of orbital channels per atom type
      integer :: max_orb_channels
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

   contains
      procedure :: generate_mp_mesh
      procedure :: generate_reciprocal_vectors
      procedure :: build_kspace_hamiltonian
      procedure :: build_kspace_hamiltonian_so
      procedure :: build_total_hamiltonian
      procedure :: calculate_structure_factors
      procedure :: fourier_transform_hamiltonian
      procedure :: set_basis_sizes
      procedure :: get_basis_type_from_size
      procedure :: diagonalize_hamiltonian
      procedure :: calculate_band_structure
      procedure :: calculate_density_of_states
      procedure :: calculate_dos_tetrahedron
      procedure :: calculate_dos_gaussian
      procedure :: setup_dos_energy_grid
      procedure :: setup_tetrahedra
      procedure :: tetrahedron_dos_contribution
      procedure :: get_kpoint_index
      procedure :: project_dos_orbitals
      procedure :: project_dos_orbitals_gaussian
      procedure :: project_dos_orbitals_tetrahedron
      procedure :: calculate_band_moments
      procedure :: find_fermi_level_from_dos
      procedure :: integrate_dos_up_to_energy
      procedure :: calculate_gaussian_weight_single
      procedure :: write_dos_to_file
      procedure :: restore_to_default
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
      call obj%generate_reciprocal_vectors()
      call obj%set_basis_sizes()
      call obj%symmetry_analysis%initialize(obj%lattice)
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
      if (allocated(this%eigenvalues)) deallocate (this%eigenvalues)
      if (allocated(this%eigenvalues_path)) deallocate (this%eigenvalues_path)
      if (allocated(this%eigenvectors)) deallocate (this%eigenvectors)
      if (allocated(this%eigenvectors_path)) deallocate (this%eigenvectors_path)
      ! DOS arrays
      if (allocated(this%dos_energy_grid)) deallocate (this%dos_energy_grid)
      if (allocated(this%total_dos)) deallocate (this%total_dos)
      if (allocated(this%projected_dos)) deallocate (this%projected_dos)
      if (allocated(this%band_moments)) deallocate (this%band_moments)
      if (allocated(this%tetrahedra)) deallocate (this%tetrahedra)
      if (allocated(this%tetrahedron_volumes)) deallocate (this%tetrahedron_volumes)
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
      this%k_offset = [0.0_rp, 0.0_rp, 0.0_rp]  ! No shift by default
      this%include_so = .false.
      this%max_orb_channels = 18  ! spd noncollinear default

   ! By default suppress internal verbose prints (can be enabled by user)
   this%suppress_internal_logs = .true.

      ! Initialize reciprocal lattice to zero
      this%reciprocal_vectors = 0.0_rp
      this%reciprocal_volume = 0.0_rp

      ! Default DOS settings
      this%n_energy_points = 1000
      this%dos_energy_range = [-10.0_rp, 10.0_rp]  ! Default energy range in eV
      this%dos_method = 'tetrahedron'  ! Default to tetrahedron method
      this%gaussian_sigma = 0.1_rp  ! Default Gaussian smearing in eV
      this%temperature = 300.0_rp  ! Default temperature in Kelvin
      this%fermi_level = 0.0_rp  ! Default Fermi level
      this%total_electrons = 0.0_rp  ! Default total electrons (will be set from input)
      this%auto_find_fermi = .false.  ! Default to using input Fermi level
      this%n_sites = 0
      this%n_orb_types = 4  ! s, p, d, f
      this%n_spin_components = 1  ! Default to non-spin-polarized
      this%n_tetrahedra = 0
   end subroutine restore_to_default

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

      this%nk_total = this%nk_mesh(1) * this%nk_mesh(2) * this%nk_mesh(3)

      ! Allocate k-point arrays
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.k_points', this%k_points, [3, this%nk_total])
      call g_safe_alloc%allocate('reciprocal.k_weights', this%k_weights, [this%nk_total])
#else
      allocate(this%k_points(3, this%nk_total))
      allocate(this%k_weights(this%nk_total))
#endif

      ! Generate Monkhorst-Pack mesh
      ik = 0
      do iz = 1, this%nk_mesh(3)
         do iy = 1, this%nk_mesh(2)
            do ix = 1, this%nk_mesh(1)
               ik = ik + 1
               
               ! Monkhorst-Pack formula: k_i = (2n_i - N_i - 1) / (2*N_i) + offset_i
               kx = (2.0_rp * ix - this%nk_mesh(1) - 1.0_rp) / (2.0_rp * this%nk_mesh(1)) + this%k_offset(1)
               ky = (2.0_rp * iy - this%nk_mesh(2) - 1.0_rp) / (2.0_rp * this%nk_mesh(2)) + this%k_offset(2)
               kz = (2.0_rp * iz - this%nk_mesh(3) - 1.0_rp) / (2.0_rp * this%nk_mesh(3)) + this%k_offset(3)

               this%k_points(1, ik) = kx
               this%k_points(2, ik) = ky
               this%k_points(3, ik) = kz

               ! Equal weights for now (should implement symmetry reduction)
               this%k_weights(ik) = 1.0_rp / real(this%nk_total, rp)
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
      allocate(this%basis_size(this%lattice%ntype))
#endif

      ! Determine basis size for each atom type
      ! This should ideally be read from input or deduced from potential info
      do ntype = 1, this%lattice%ntype
         ! For now, assume spd basis (9 orbitals) for all atom types
         ! TODO: Extend to read from symbolic_atoms potential configuration
         this%basis_size(ntype) = 9  ! spd basis: s(1) + p(3) + d(5) = 9 orbitals
         
         ! Alternative basis sizes:
         ! sp basis: s(1) + p(3) = 4 orbitals
         ! spdf basis: s(1) + p(3) + d(5) + f(7) = 16 orbitals
      end do

      ! Set maximum orbital channels (×2 for spin, stored as 18 for noncollinear spd)
      this%max_orb_channels = maxval(this%basis_size) * 2
      
      ! For current spd case: 9 orbitals × 2 = 18 for noncollinear
      if (maxval(this%basis_size) == 9) then
         this%max_orb_channels = 18
      else if (maxval(this%basis_size) == 4) then
         this%max_orb_channels = 8   ! sp basis
      else if (maxval(this%basis_size) == 16) then
         this%max_orb_channels = 32  ! spdf basis
      end if

      call g_logger%info('reciprocal%set_basis_sizes: Basis sizes set: max_orb_channels = ' // trim(int2str(this%max_orb_channels)), __FILE__, __LINE__)
   end subroutine set_basis_sizes

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate structure factors for Fourier transform
   !---------------------------------------------------------------------------
   subroutine calculate_structure_factors(this, k_vec, ntype, structure_factors)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      integer, intent(in) :: ntype
      complex(rp), dimension(:), intent(out) :: structure_factors
      ! Local variables
      integer :: ineigh, ia, nr
      real(rp) :: k_dot_r
      real(rp), dimension(3) :: r_vec
   logical :: debug_this_k

   ia = this%lattice%atlist(ntype)
   nr = this%lattice%nn(ia, 1)  ! Number of neighbors

   ! Only enable debug_this_k when user allows internal logs and for first atom type
   debug_this_k = .false.
   if (.not. this%suppress_internal_logs) debug_this_k = (ntype == 1)

   ! Calculate structure factors exp(i*k·R) for each neighbor
   do ineigh = 1, min(nr, size(structure_factors))
         if (ineigh == 1) then
            ! On-site term (R = 0)
            structure_factors(ineigh) = cmplx(1.0_rp, 0.0_rp, rp)
            if (debug_this_k) then
               call g_logger%info('fourier_transform_hamiltonian: On-site term (R=0)', __FILE__, __LINE__)
            end if
         else
            ! Off-site terms - get neighbor vector from lattice structure
            if (allocated(this%lattice%sbarvec_direct) .and. ineigh <= size(this%lattice%sbarvec_direct, 2)) then
               r_vec(1) = this%lattice%sbarvec_direct(1, ineigh)
               r_vec(2) = this%lattice%sbarvec_direct(2, ineigh)
               r_vec(3) = this%lattice%sbarvec_direct(3, ineigh)
               if (debug_this_k) then
                  call g_logger%info('fourier_transform_hamiltonian: Using sbarvec_direct neighbor vector', __FILE__, __LINE__)
               end if
            else if (allocated(this%lattice%sbarvec) .and. ineigh <= size(this%lattice%sbarvec, 2)) then
               r_vec(1) = this%lattice%sbarvec(1, ineigh)
               r_vec(2) = this%lattice%sbarvec(2, ineigh)
               r_vec(3) = this%lattice%sbarvec(3, ineigh)
               if (debug_this_k) then
                  call g_logger%info('fourier_transform_hamiltonian: Using sbarvec neighbor vector', __FILE__, __LINE__)
               end if
            else
               ! Fallback to zero vector if sbarvec not available
               r_vec = 0.0_rp
               if (debug_this_k) then
                  call g_logger%info('fourier_transform_hamiltonian: WARNING - Using zero fallback vector!', __FILE__, __LINE__)
               end if
            end if
            k_dot_r = 2.0_rp * 3.141592653589793_rp * dot_product(k_vec, r_vec)  ! 2π * k_frac · R_frac
            structure_factors(ineigh) = cmplx(cos(k_dot_r), sin(k_dot_r), rp)
         end if
      end do
   end subroutine calculate_structure_factors

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Fourier transform real-space Hamiltonian to k-space
   !---------------------------------------------------------------------------
   subroutine fourier_transform_hamiltonian(this, k_vec, ntype, hk_result)
      class(reciprocal), intent(in) :: this
      real(rp), dimension(3), intent(in) :: k_vec
      integer, intent(in) :: ntype
      complex(rp), dimension(:, :), intent(out) :: hk_result
      ! Local variables
      integer :: ineigh, ia, nr
      complex(rp), dimension(:), allocatable :: structure_factors
      complex(rp), dimension(18, 18) :: h_neighbor
      logical :: debug_this_k

      ia = this%lattice%atlist(ntype)
      nr = this%lattice%nn(ia, 1)  ! Number of neighbors

      ! Debug only for first k-point and first atom type to avoid spam
      debug_this_k = (ntype == 1)
      
      ! Only debug near gamma point
      if (sqrt(k_vec(1)**2 + k_vec(2)**2 + k_vec(3)**2) > 0.1_rp) debug_this_k = .false.

      if (.not. this%suppress_internal_logs) then
         if (debug_this_k) then
            call g_logger%info('fourier_transform_hamiltonian: Debug output for k-vector', __FILE__, __LINE__)
            if (allocated(this%lattice%sbarvec)) then
               call g_logger%info('fourier_transform_hamiltonian: sbarvec is allocated', __FILE__, __LINE__)
            else
               call g_logger%info('fourier_transform_hamiltonian: WARNING - sbarvec is NOT allocated!', __FILE__, __LINE__)
            end if
         end if

         ! Add simple check for ANY call during band structure calculation
         ! if (ntype == 1) then
         !    if (allocated(this%lattice%sbarvec)) then
         !       call g_logger%info('fourier_transform_hamiltonian: sbarvec allocated for band structure', __FILE__, __LINE__)
         !    else
         !       call g_logger%info('fourier_transform_hamiltonian: WARNING - sbarvec NOT allocated for band structure!', __FILE__, __LINE__)
         !    end if
         ! end if
      end if

      ! Allocate structure factors
      allocate(structure_factors(nr))

      ! Calculate structure factors
      call this%calculate_structure_factors(k_vec, ntype, structure_factors)
      
      ! Print first few structure factors for first k-point to debug
     ! if (.not. this%suppress_internal_logs) then
     !   if (ntype == 1 .and. size(structure_factors) >= 5) then
     !      write(*, '(A,2F10.6)') "First structure factor:", structure_factors(1)
     !      write(*, '(A,2F10.6)') "Second structure factor:", structure_factors(2)
     !      write(*, '(A,2F10.6)') "Third structure factor:", structure_factors(3)
     !      write(*, '(A,2F10.6)') "Fourth structure factor:", structure_factors(4)
     !      write(*, '(A,2F10.6)') "Fifth structure factor:", structure_factors(5)
     !   end if
     ! end if

      ! Initialize result
      hk_result = cmplx(0.0_rp, 0.0_rp, rp)

      ! Sum over neighbors: H(k) = Σ_R H(R) * exp(i*k·R)
      do ineigh = 1, nr
         h_neighbor = this%hamiltonian%ee(:, :, ineigh, ntype)
         hk_result = hk_result + h_neighbor * structure_factors(ineigh)
      end do

      ! if (debug_this_k) then
      !    call g_logger%info('fourier_transform_hamiltonian: Completed FT', __FILE__, __LINE__)
      ! end if

      deallocate(structure_factors)
      
   contains
      ! Helper function to calculate trace of real matrix
      function trace_real_matrix(mat) result(tr)
         complex(rp), dimension(:, :), intent(in) :: mat
         real(rp) :: tr
         integer :: i
         tr = 0.0_rp
         do i = 1, min(size(mat, 1), size(mat, 2))
            tr = tr + real(mat(i, i))
         end do
      end function trace_real_matrix
      
      ! Helper function to calculate trace of complex matrix
      function trace_complex_matrix(mat) result(tr)
         complex(rp), dimension(:, :), intent(in) :: mat
         complex(rp) :: tr
         integer :: i
         tr = cmplx(0.0_rp, 0.0_rp, rp)
         do i = 1, min(size(mat, 1), size(mat, 2))
            tr = tr + mat(i, i)
         end do
      end function trace_complex_matrix
   end subroutine fourier_transform_hamiltonian

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build k-space Hamiltonian for all k-points (bulk contribution)
   !---------------------------------------------------------------------------
   subroutine build_kspace_hamiltonian(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ik, ntype
      real(rp), dimension(3) :: k_cartesian
      character(len=200) :: debug_msg

      if (.not. allocated(this%k_points)) then
         call g_logger%error('reciprocal%build_kspace_hamiltonian: K-points not generated. Call generate_mp_mesh first.', __FILE__, __LINE__)
         return
      end if

      call g_logger%info('reciprocal%build_kspace_hamiltonian: Starting Fourier transform of real-space Hamiltonian', __FILE__, __LINE__)
      write(debug_msg, '(A,I0,A,I0,A)') 'build_kspace_hamiltonian: Building for ', this%nk_total, ' k-points and ', this%lattice%ntype, ' atom types'
      call g_logger%info(trim(debug_msg), __FILE__, __LINE__)

      ! Allocate k-space Hamiltonian
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.hk_bulk', this%hk_bulk, &
                                [this%max_orb_channels, this%max_orb_channels, this%nk_total, this%lattice%ntype])
#else
      allocate(this%hk_bulk(this%max_orb_channels, this%max_orb_channels, this%nk_total, this%lattice%ntype))
#endif

      ! Build Hamiltonian for each k-point and atom type
      do ntype = 1, this%lattice%ntype
         write(debug_msg, '(A,I0)') 'build_kspace_hamiltonian: Processing atom type ', ntype
         call g_logger%info(trim(debug_msg), __FILE__, __LINE__)
         
         do ik = 1, this%nk_total
            ! Convert k-point to Cartesian coordinates
            k_cartesian = matmul(this%reciprocal_vectors, this%k_points(:, ik))
            
            ! Debug first few k-points
            if (.not. this%suppress_internal_logs) then
               if (ik <= 3) then
                  write(debug_msg, '(A,I0,A,3F12.6,A,3F12.6,A)') 'build_kspace_hamiltonian: k-point ', ik, &
                        ' fractional=(', this%k_points(:, ik), ') cartesian=(', k_cartesian, ')'
                  call g_logger%info(trim(debug_msg), __FILE__, __LINE__)
               end if
            end if
            
            ! Fourier transform real-space Hamiltonian
            call this%fourier_transform_hamiltonian(k_cartesian, ntype, this%hk_bulk(:, :, ik, ntype))
         end do
      end do

      call g_logger%info('reciprocal%build_kspace_hamiltonian: K-space bulk Hamiltonian built', __FILE__, __LINE__)
   end subroutine build_kspace_hamiltonian

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build k-space spin-orbit Hamiltonian
   !---------------------------------------------------------------------------
   subroutine build_kspace_hamiltonian_so(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ik, ntype

      if (.not. this%include_so) then
         call g_logger%info('reciprocal%build_kspace_hamiltonian_so: Spin-orbit coupling disabled', __FILE__, __LINE__)
         return
      end if

      if (.not. allocated(this%k_points)) then
         call g_logger%error('reciprocal%build_kspace_hamiltonian_so: K-points not generated. Call generate_mp_mesh first.', __FILE__, __LINE__)
         return
      end if

      if (.not. allocated(this%hamiltonian%lsham)) then
         call g_logger%error('reciprocal%build_kspace_hamiltonian_so: Real-space SO Hamiltonian not built. Call hamiltonian%build_lsham first.', __FILE__, __LINE__)
         return
      end if

      ! Allocate k-space SO Hamiltonian
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.hk_so', this%hk_so, &
                                [this%max_orb_channels, this%max_orb_channels, this%lattice%ntype])
#else
      allocate(this%hk_so(this%max_orb_channels, this%max_orb_channels, this%lattice%ntype))
#endif

      ! SO coupling is local (on-site), so it's k-independent
      ! Copy from real-space lsham matrix
      do ntype = 1, this%lattice%ntype
         this%hk_so(:, :, ntype) = this%hamiltonian%lsham(:, :, ntype)
      end do

      call g_logger%info('reciprocal%build_kspace_hamiltonian_so: K-space SO Hamiltonian built', __FILE__, __LINE__)
   end subroutine build_kspace_hamiltonian_so

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build total k-space Hamiltonian (bulk + SO contributions)
   !---------------------------------------------------------------------------
   subroutine build_total_hamiltonian(this)
      class(reciprocal), intent(inout) :: this
      ! Local variables
      integer :: ik, ntype

      if (.not. allocated(this%hk_bulk)) then
         call g_logger%error('reciprocal%build_total_hamiltonian: Bulk k-space Hamiltonian not built. Call build_kspace_hamiltonian first.', __FILE__, __LINE__)
         return
      end if

      ! Allocate total Hamiltonian
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('reciprocal.hk_total', this%hk_total, &
                                [this%max_orb_channels, this%max_orb_channels, this%nk_total])
#else
      allocate(this%hk_total(this%max_orb_channels, this%max_orb_channels, this%nk_total))
#endif

      ! Combine bulk and SO contributions for each k-point
      do ik = 1, this%nk_total
         ! For now, use first atom type - should be generalized for multi-component systems
         ntype = 1
         this%hk_total(:, :, ik) = this%hk_bulk(:, :, ik, ntype)
         
         ! Add spin-orbit coupling if available
         if (this%include_so .and. allocated(this%hk_so)) then
            this%hk_total(:, :, ik) = this%hk_total(:, :, ik) + this%hk_so(:, :, ntype)
         end if
      end do

      call g_logger%info('reciprocal%build_total_hamiltonian: Total k-space Hamiltonian built', __FILE__, __LINE__)
   end subroutine build_total_hamiltonian

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
   !> Diagonalize Hamiltonian at k-points (wrapper for unified method)
   !---------------------------------------------------------------------------
   subroutine diagonalize_hamiltonian(this, ham, use_kpath)
      class(reciprocal), intent(inout) :: this
      class(hamiltonian), intent(in) :: ham
      logical, intent(in), optional :: use_kpath
      ! Local variables
      logical :: use_path
      integer :: nk, i, nmat, lwork, info
      complex(rp), dimension(:, :), allocatable :: h_k, h_k_copy
      real(rp), dimension(:), allocatable :: eigenvals
      complex(rp), dimension(:, :), allocatable :: eigenvecs
      complex(rp), dimension(:), allocatable :: work_complex
      real(rp), dimension(:), allocatable :: rwork
      character(len=100) :: info_msg

      use_path = .false.
      if (present(use_kpath)) use_path = use_kpath

      if (use_path) then
         if (.not. allocated(this%k_path)) then
            call g_logger%error('diagonalize_hamiltonian: k-path not generated', __FILE__, __LINE__)
            return
         end if
         nk = this%nk_path
         call g_logger%info('diagonalize_hamiltonian: Using k-path with ' // trim(int2str(nk)) // ' points', __FILE__, __LINE__)
      else
         if (.not. allocated(this%k_points)) then
            call g_logger%error('diagonalize_hamiltonian: k-mesh not generated', __FILE__, __LINE__)
            return
         end if
         nk = this%nk_total
         call g_logger%info('diagonalize_hamiltonian: Using k-mesh with ' // trim(int2str(nk)) // ' points', __FILE__, __LINE__)
      end if

      call g_logger%info('diagonalize_hamiltonian: Starting diagonalization of ' // trim(int2str(nk)) // ' k-points', __FILE__, __LINE__)

      ! Get matrix size from the actual Hamiltonian
      if (allocated(ham%ee)) then
         nmat = size(ham%ee, 1)  ! Use actual Hamiltonian dimension
      else
         nmat = this%max_orb_channels  ! Fallback to max orbital channels
      end if

      call g_logger%info('diagonalize_hamiltonian: Matrix size is ' // trim(int2str(nmat)) // 'x' // trim(int2str(nmat)), __FILE__, __LINE__)

      ! Allocate eigenvalue and eigenvector arrays
      if (allocated(this%eigenvalues)) deallocate(this%eigenvalues)
      if (allocated(this%eigenvectors)) deallocate(this%eigenvectors)
      
      allocate(this%eigenvalues(nmat, nk))
      allocate(this%eigenvectors(nmat, nmat, nk))

      ! Allocate work arrays
      allocate(h_k(nmat, nmat))
      allocate(h_k_copy(nmat, nmat))
      allocate(eigenvals(nmat))
      allocate(eigenvecs(nmat, nmat))
      
      ! Query optimal work array size
      lwork = -1
      allocate(work_complex(1))
      allocate(rwork(3*nmat-2))
      
      call zheev('V', 'U', nmat, h_k, nmat, eigenvals, work_complex, lwork, rwork, info)
      lwork = int(real(work_complex(1)))
      deallocate(work_complex)
      allocate(work_complex(lwork))

      ! Loop over k-points
      do i = 1, nk
         ! Build Hamiltonian at this k-point using unified Fourier transform
         if (use_path) then
            call this%fourier_transform_hamiltonian(this%k_path(:, i), 1, h_k)
            
            !!! ! Debug: Dump H(k) matrix for Γ and H points
            !!! if (all(abs(this%k_path(:, i) - [0.0_rp, 0.0_rp, 0.0_rp]) < 1.0e-6_rp)) then
            !!!    ! Γ point
            !!!    call dump_complex_matrix(h_k, 'H_k_Gamma.dat', this%k_path(:, i))
            !!!    call g_logger%info('diagonalize_hamiltonian: Dumped H(k) matrix for Γ point to H_k_Gamma.dat, k-vector: [' // &
            !!!                      trim(real2str(this%k_path(1,i))) // ', ' // trim(real2str(this%k_path(2,i))) // ', ' // &
            !!!                      trim(real2str(this%k_path(3,i))) // ']', __FILE__, __LINE__)
            !!! else if (all(abs(this%k_path(:, i) - [0.5_rp, -0.5_rp, 0.5_rp]) < 1.0e-6_rp)) then
            !!!    ! H point
            !!!    call dump_complex_matrix(h_k, 'H_k_H.dat', this%k_path(:, i))
            !!!    call g_logger%info('diagonalize_hamiltonian: Dumped H(k) matrix for H point to H_k_H.dat, k-vector: [' // &
            !!!                      trim(real2str(this%k_path(1,i))) // ', ' // trim(real2str(this%k_path(2,i))) // ', ' // &
            !!!                      trim(real2str(this%k_path(3,i))) // ']', __FILE__, __LINE__)
            !!! else if (all(abs(this%k_path(:, i) - [0.25_rp, -0.25_rp, 0.25_rp]) < 1.0e-6_rp)) then
            !!!    ! K/2 point (midpoint between Γ and H)
            !!!    call dump_complex_matrix(h_k, 'H_k_half.dat', this%k_path(:, i))
            !!!    call g_logger%info('diagonalize_hamiltonian: Dumped H(k) matrix for K/2 point to H_k_half.dat, k-vector: [' // &
            !!!                      trim(real2str(this%k_path(1,i))) // ', ' // trim(real2str(this%k_path(2,i))) // ', ' // &
            !!!                      trim(real2str(this%k_path(3,i))) // ']', __FILE__, __LINE__)
            !!! end if
         else
            call this%fourier_transform_hamiltonian(this%k_points(:, i), 1, h_k)
         end if

         ! Make a copy for diagonalization (LAPACK destroys input)
         h_k_copy = h_k

         ! Diagonalize using LAPACK ZHEEV
         call zheev('V', 'U', nmat, h_k_copy, nmat, eigenvals, work_complex, lwork, rwork, info)
         
         if (info /= 0) then
            write(info_msg, '(A,I0,A,I0)') 'diagonalize_hamiltonian: ZHEEV failed for k-point ', i, ' with info=', info
            call g_logger%error(trim(info_msg), __FILE__, __LINE__)
            cycle
         end if

         ! Store results
         this%eigenvalues(:, i) = eigenvals
         this%eigenvectors(:, :, i) = h_k_copy
      end do

      call g_logger%info('diagonalize_hamiltonian: Completed diagonalization', __FILE__, __LINE__)

      ! Clean up
      deallocate(h_k, h_k_copy, eigenvals, eigenvecs)
      deallocate(work_complex, rwork)
   end subroutine diagonalize_hamiltonian


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

      ! Diagonalize Hamiltonian along k-path using unified method
      call this%diagonalize_hamiltonian(ham, use_kpath=.true.)

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
   !> Generate symmetry-reduced k-point mesh using spglib
   !---------------------------------------------------------------------------
   subroutine generate_reduced_kpoint_mesh(this, mesh_dims, use_shift)
      class(reciprocal), intent(inout) :: this
      integer, intent(in) :: mesh_dims(3)
      logical, intent(in), optional :: use_shift
      integer :: shift(3)
      integer :: num_ir_kpoints
      logical :: do_shift

      do_shift = .false.
      if (present(use_shift)) do_shift = use_shift

      if (do_shift) then
         shift = [1, 1, 1]  ! Offset by half a mesh spacing
      else
         shift = [0, 0, 0]  ! No offset
      end if

      if (.not. this%symmetry_analysis%spglib%is_available()) then
         call g_logger%warning('generate_reduced_kpoint_mesh: spglib not available, using full mesh', __FILE__, __LINE__)
         call this%generate_mp_mesh()  ! Fall back to regular MP mesh
         return
      end if

      ! Get irreducible k-points using spglib
      num_ir_kpoints = this%symmetry_analysis%spglib%get_reduced_kpoint_mesh(mesh_dims, shift)

      call g_logger%info('generate_reduced_kpoint_mesh: Generated ' // trim(int2str(num_ir_kpoints)) // &
                        ' irreducible k-points from ' // trim(int2str(product(mesh_dims))) // ' total points', &
                        __FILE__, __LINE__)

      ! TODO: Implement actual k-point generation from spglib results
      ! For now, fall back to regular MP mesh
      call this%generate_mp_mesh()
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

      ! Diagonalize Hamiltonian on k-mesh if not already done
      if (.not. allocated(this%eigenvalues)) then
         call g_logger%info('calculate_density_of_states: Diagonalizing Hamiltonian on k-mesh', __FILE__, __LINE__)
         call this%diagonalize_hamiltonian(ham, use_kpath=.false.)
      end if

      ! Calculate DOS based on method
      select case (trim(this%dos_method))
      case ('tetrahedron')
         call g_logger%info('calculate_density_of_states: Using tetrahedron method', __FILE__, __LINE__)
         call this%calculate_dos_tetrahedron()
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
   !> Setup energy grid for DOS calculation
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

      ! Fill energy grid
      do i = 1, this%n_energy_points
         this%dos_energy_grid(i) = energy_min + real(i-1, rp) * delta_energy
      end do

      call g_logger%info('setup_dos_energy_grid: Created energy grid with ' // &
                        trim(int2str(this%n_energy_points)) // ' points from ' // &
                        trim(real2str(energy_min, '(F 8.5)')) // ' to ' // trim(real2str(energy_max, '(F 8.5)')) // ' eV', &
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
      integer :: i_energy, i_tet, i_corner, i_band
      real(rp) :: energy, dos_contrib
      real(rp), dimension(4) :: e_corners, sorted_e
      integer, dimension(4) :: sort_idx

      call g_logger%info('calculate_dos_tetrahedron: Calculating DOS using tetrahedron method', __FILE__, __LINE__)

      ! Setup tetrahedra if not already done
      if (.not. allocated(this%tetrahedra)) then
         call this%setup_tetrahedra()
      end if

      ! Allocate DOS arrays
      if (allocated(this%total_dos)) deallocate(this%total_dos)
      allocate(this%total_dos(this%n_energy_points))
      this%total_dos = 0.0_rp

      ! Loop over energy points
      do i_energy = 1, this%n_energy_points
         energy = this%dos_energy_grid(i_energy)

         ! Loop over tetrahedra
         do i_tet = 1, this%n_tetrahedra

            ! Loop over bands
            do i_band = 1, size(this%eigenvalues, 1)

               ! Get eigenvalues at tetrahedron corners
               do i_corner = 1, 4
                  e_corners(i_corner) = this%eigenvalues(i_band, this%tetrahedra(i_corner, i_tet))
               end do

               ! Sort eigenvalues
               call sort_real_array(e_corners, sorted_e, sort_idx)

               ! Calculate DOS contribution from this tetrahedron and band
               dos_contrib = this%tetrahedron_dos_contribution(energy, sorted_e)

               ! Add to total DOS (weight by tetrahedron volume and k-point weight)
               this%total_dos(i_energy) = this%total_dos(i_energy) + &
                                        dos_contrib * this%tetrahedron_volumes(i_tet)
            end do
         end do
      end do

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
   ! DESCRIPTION:
   !> @brief
   !> Calculate DOS using Gaussian smearing
   !---------------------------------------------------------------------------
   subroutine calculate_dos_gaussian(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i_energy, i_k, i_band
      real(rp) :: energy, weight, gaussian_factor
      real(rp) :: sigma_squared

      call g_logger%info('calculate_dos_gaussian: Calculating DOS with Gaussian smearing, sigma = ' // &
                        trim(real2str(this%gaussian_sigma)) // ' eV', __FILE__, __LINE__)

      ! Allocate DOS arrays
      if (allocated(this%total_dos)) deallocate(this%total_dos)
      allocate(this%total_dos(this%n_energy_points))
      this%total_dos = 0.0_rp

      sigma_squared = this%gaussian_sigma**2

      ! Loop over energy points
      do i_energy = 1, this%n_energy_points
         energy = this%dos_energy_grid(i_energy)

         ! Loop over k-points
         do i_k = 1, this%nk_total
            weight = this%k_weights(i_k)

            ! Loop over bands
            do i_band = 1, size(this%eigenvalues, 1)
               gaussian_factor = exp(-((energy - this%eigenvalues(i_band, i_k))**2) / (2.0_rp * sigma_squared))
               gaussian_factor = gaussian_factor / (this%gaussian_sigma * sqrt(2.0_rp * 3.141592653589793_rp))

               this%total_dos(i_energy) = this%total_dos(i_energy) + weight * gaussian_factor
            end do
         end do
      end do

      call g_logger%info('calculate_dos_gaussian: Gaussian DOS calculation completed', __FILE__, __LINE__)
   end subroutine calculate_dos_gaussian

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
   !> Project DOS onto orbitals, sites, and spin
   !---------------------------------------------------------------------------
   subroutine project_dos_orbitals(this)
      class(reciprocal), intent(inout) :: this

      call g_logger%info('project_dos_orbitals: Starting orbital projection calculation', __FILE__, __LINE__)

      ! Use tetrahedron or Gaussian method based on dos_method
      if (trim(this%dos_method) == 'tetrahedron') then
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
      integer :: ik, ib, ie, iorb, ispin, i
      integer :: n_orb_per_spin, orb_start
      real(rp) :: weight, orbital_char, energy
      complex(rp) :: psi_element

      call g_logger%info('project_dos_orbitals: Starting orbital projection calculation', __FILE__, __LINE__)

      ! Initialize dimensions
      this%n_sites = 1  ! For now, single site (can be extended)
      this%n_orb_types = 4  ! s, p, d, f
      this%n_spin_components = 2  ! spin up/down

      ! Allocate projected DOS array
      if (allocated(this%projected_dos)) deallocate(this%projected_dos)
      allocate(this%projected_dos(this%n_sites, this%n_orb_types, this%n_spin_components, this%n_energy_points))
      this%projected_dos = 0.0_rp

      ! Number of orbitals per spin (assuming spd basis: 9 orbitals per spin)
      n_orb_per_spin = this%max_orb_channels / 2

      ! For now, only implement Gaussian method for projections
      if (trim(this%dos_method) /= 'gaussian') then
         call g_logger%warning('project_dos_orbitals: Tetrahedron projections not yet implemented, using Gaussian', __FILE__, __LINE__)
      end if

      ! Loop over energy points
      do ie = 1, this%n_energy_points
         energy = this%dos_energy_grid(ie)

         ! Loop over k-points
         do ik = 1, this%nk_total
            ! Loop over bands
            do ib = 1, this%max_orb_channels
               ! Skip if eigenvalue is far from current energy
               if (abs(this%eigenvalues(ib, ik) - energy) > 5.0_rp * this%gaussian_sigma) cycle

               ! Calculate Gaussian weight
               weight = this%calculate_gaussian_weight_single(energy, this%eigenvalues(ib, ik))

               ! Skip if weight is negligible
               if (abs(weight) < 1.0e-10_rp) cycle

               ! Apply k-point weight
               weight = weight * this%k_weights(ik)

               ! Calculate orbital character for each orbital type and spin
               do ispin = 1, this%n_spin_components
                  ! Orbital range for this spin
                  orb_start = (ispin-1) * n_orb_per_spin + 1

                  ! s orbital (index 1 in each spin block)
                  iorb = 1
                  psi_element = this%eigenvectors(orb_start, ib, ik)
                  orbital_char = real(conjg(psi_element) * psi_element, rp)
                  this%projected_dos(1, iorb, ispin, ie) = this%projected_dos(1, iorb, ispin, ie) + &
                                                         orbital_char * weight

                  ! p orbitals (indices 2-4 in each spin block)
                  iorb = 2
                  orbital_char = 0.0_rp
                  do i = 2, 4
                     psi_element = this%eigenvectors(orb_start + i - 1, ib, ik)
                     orbital_char = orbital_char + real(conjg(psi_element) * psi_element, rp)
                  end do
                  this%projected_dos(1, iorb, ispin, ie) = this%projected_dos(1, iorb, ispin, ie) + &
                                                         orbital_char * weight

                  ! d orbitals (indices 5-9 in each spin block)
                  iorb = 3
                  orbital_char = 0.0_rp
                  do i = 5, 9
                     psi_element = this%eigenvectors(orb_start + i - 1, ib, ik)
                     orbital_char = orbital_char + real(conjg(psi_element) * psi_element, rp)
                  end do
                  this%projected_dos(1, iorb, ispin, ie) = this%projected_dos(1, iorb, ispin, ie) + &
                                                         orbital_char * weight

                  ! f orbitals (would be indices 10-16, but not present in spd basis)
                  iorb = 4
                  orbital_char = 0.0_rp  ! No f orbitals in current spd basis
                  this%projected_dos(1, iorb, ispin, ie) = this%projected_dos(1, iorb, ispin, ie) + &
                                                         orbital_char * weight
               end do
            end do
         end do
      end do

      call g_logger%info('project_dos_orbitals_gaussian: Gaussian orbital projection calculation completed', __FILE__, __LINE__)
   end subroutine project_dos_orbitals_gaussian

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Project DOS onto orbitals, sites, and spin using tetrahedron method
   !---------------------------------------------------------------------------
   subroutine project_dos_orbitals_tetrahedron(this)
      class(reciprocal), intent(inout) :: this

      ! Local variables
      integer :: i_energy, i_tet, i_corner, i_band, iorb, ispin, i
      integer :: n_orb_per_spin, orb_start, ik
      real(rp) :: energy, dos_contrib, orbital_char_avg, orbital_char
      real(rp), dimension(4) :: e_corners, sorted_e, orbital_chars
      integer, dimension(4) :: sort_idx
      complex(rp) :: psi_element

      call g_logger%info('project_dos_orbitals_tetrahedron: Starting tetrahedron orbital projection calculation', __FILE__, __LINE__)

      ! Initialize dimensions
      this%n_sites = 1  ! For now, single site (can be extended)
      this%n_orb_types = 4  ! s, p, d, f
      this%n_spin_components = 2  ! spin up/down

      ! Allocate projected DOS array
      if (allocated(this%projected_dos)) deallocate(this%projected_dos)
      allocate(this%projected_dos(this%n_sites, this%n_orb_types, this%n_spin_components, this%n_energy_points))
      this%projected_dos = 0.0_rp

      ! Number of orbitals per spin (assuming spd basis: 9 orbitals per spin)
      n_orb_per_spin = this%max_orb_channels / 2

      ! Setup tetrahedra if not already done
      if (.not. allocated(this%tetrahedra)) then
         call this%setup_tetrahedra()
      end if

      ! Loop over energy points
      do i_energy = 1, this%n_energy_points
         energy = this%dos_energy_grid(i_energy)

         ! Loop over tetrahedra
         do i_tet = 1, this%n_tetrahedra

            ! Loop over bands
            do i_band = 1, this%max_orb_channels

               ! Get eigenvalues at tetrahedron corners
               do i_corner = 1, 4
                  ik = this%tetrahedra(i_corner, i_tet)
                  e_corners(i_corner) = this%eigenvalues(i_band, ik)
               end do

               ! Sort eigenvalues
               call sort_real_array(e_corners, sorted_e, sort_idx)

               ! Calculate DOS contribution from this tetrahedron and band
               dos_contrib = this%tetrahedron_dos_contribution(energy, sorted_e)

               ! Skip if DOS contribution is negligible
               if (abs(dos_contrib) < 1.0e-12_rp) cycle

               ! Weight by tetrahedron volume
               dos_contrib = dos_contrib * this%tetrahedron_volumes(i_tet)

               ! Calculate orbital character for each orbital type and spin
               do ispin = 1, this%n_spin_components
                  ! Orbital range for this spin
                  orb_start = (ispin-1) * n_orb_per_spin + 1

                  ! s orbital (index 1 in each spin block)
                  iorb = 1
                  orbital_char_avg = 0.0_rp
                  do i_corner = 1, 4
                     ik = this%tetrahedra(i_corner, i_tet)
                     psi_element = this%eigenvectors(orb_start, i_band, ik)
                     orbital_chars(i_corner) = real(conjg(psi_element) * psi_element, rp)
                  end do
                  orbital_char_avg = sum(orbital_chars) / 4.0_rp  ! Average over tetrahedron corners
                  this%projected_dos(1, iorb, ispin, i_energy) = this%projected_dos(1, iorb, ispin, i_energy) + &
                                                               orbital_char_avg * dos_contrib

                  ! p orbitals (indices 2-4 in each spin block)
                  iorb = 2
                  orbital_char_avg = 0.0_rp
                  do i_corner = 1, 4
                     ik = this%tetrahedra(i_corner, i_tet)
                     orbital_char = 0.0_rp
                     do i = 2, 4
                        psi_element = this%eigenvectors(orb_start + i - 1, i_band, ik)
                        orbital_char = orbital_char + real(conjg(psi_element) * psi_element, rp)
                     end do
                     orbital_chars(i_corner) = orbital_char
                  end do
                  orbital_char_avg = sum(orbital_chars) / 4.0_rp
                  this%projected_dos(1, iorb, ispin, i_energy) = this%projected_dos(1, iorb, ispin, i_energy) + &
                                                               orbital_char_avg * dos_contrib

                  ! d orbitals (indices 5-9 in each spin block)
                  iorb = 3
                  orbital_char_avg = 0.0_rp
                  do i_corner = 1, 4
                     ik = this%tetrahedra(i_corner, i_tet)
                     orbital_char = 0.0_rp
                     do i = 5, 9
                        psi_element = this%eigenvectors(orb_start + i - 1, i_band, ik)
                        orbital_char = orbital_char + real(conjg(psi_element) * psi_element, rp)
                     end do
                     orbital_chars(i_corner) = orbital_char
                  end do
                  orbital_char_avg = sum(orbital_chars) / 4.0_rp
                  this%projected_dos(1, iorb, ispin, i_energy) = this%projected_dos(1, iorb, ispin, i_energy) + &
                                                               orbital_char_avg * dos_contrib

                  ! f orbitals (would be indices 10-16, but not present in spd basis)
                  iorb = 4
                  orbital_char_avg = 0.0_rp  ! No f orbitals in current spd basis
                  this%projected_dos(1, iorb, ispin, i_energy) = this%projected_dos(1, iorb, ispin, i_energy) + &
                                                               orbital_char_avg * dos_contrib
               end do
            end do
         end do
      end do

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

      ! Local variables
      integer :: isite, iorb, ispin, ie, n_energy
      real(rp) :: energy, dos_value, delta_energy, fermi_weight
      real(rp) :: m0, m1, m2, norm_factor
      real(rp), dimension(:), allocatable :: integrand
      real(rp) :: kT, fermi_arg

      call g_logger%info('calculate_band_moments: Starting band moments calculation', __FILE__, __LINE__)
      call g_logger%info('calculate_band_moments: Fermi level = ' // trim(real2str(this%fermi_level)) // ' Ry', __FILE__, __LINE__)
      call g_logger%info('calculate_band_moments: Temperature = ' // trim(real2str(this%temperature)) // ' K', __FILE__, __LINE__)

      ! Automatically find Fermi level from DOS if requested
      if (this%auto_find_fermi) then
         if (this%total_electrons > 0.0_rp) then
            this%fermi_level = this%find_fermi_level_from_dos(this%total_electrons)
            call g_logger%info('calculate_band_moments: Auto-found Fermi level = ' // trim(real2str(this%fermi_level, '(F 8.5)')) // ' Ry', __FILE__, __LINE__)
         else
            call g_logger%warning('calculate_band_moments: auto_find_fermi is true but total_electrons not set, using input Fermi level', __FILE__, __LINE__)
         end if
      end if

      ! Allocate band moments array if not already done
      if (allocated(this%band_moments)) deallocate(this%band_moments)
      allocate(this%band_moments(this%n_sites, this%n_orb_types, this%n_spin_components, 3))
      this%band_moments = 0.0_rp

      n_energy = this%n_energy_points

      ! Allocate temporary array for integration
      allocate(integrand(n_energy))

      ! Boltzmann constant in eV/K
      kT = this%temperature * 8.617333262145e-5_rp

      ! Calculate moments for each site, orbital type, and spin component
      do isite = 1, this%n_sites
         do iorb = 1, this%n_orb_types-1 ! Skip f orbitals if not present
            do ispin = 1, this%n_spin_components

               ! Initialize moment accumulators
               m0 = 0.0_rp
               m1 = 0.0_rp
               m2 = 0.0_rp

               ! Calculate zeroth moment (m0) = ∫ DOS(E) * f(E) dE
               ! where f(E) is the Fermi-Dirac distribution
               integrand = 0.0_rp
               do ie = 1, n_energy
                  energy = this%dos_energy_grid(ie)
                  fermi_arg = (energy - this%fermi_level) / kT
                  
                  ! Fermi-Dirac distribution: f(E) = 1 / (exp((E-Ef)/kT) + 1)
                  if (fermi_arg > 50.0_rp) then
                     fermi_weight = 0.0_rp  ! exp(-50) ≈ 0
                  else if (fermi_arg < -50.0_rp) then
                     fermi_weight = 1.0_rp  ! exp(50) ≈ very large
                  else
                     fermi_weight = 1.0_rp / (exp(fermi_arg) + 1.0_rp)
                  end if
                  
                  integrand(ie) = this%projected_dos(isite, iorb, ispin, ie) * fermi_weight
               end do
               m0 = trapezoidal_integral(this%dos_energy_grid, integrand)

               ! Calculate first moment (m1) = ∫ E * DOS(E) * f(E) dE / m0
               if (abs(m0) > 1.0e-12_rp) then
                  integrand = this%dos_energy_grid * this%projected_dos(isite, iorb, ispin, :) * &
                             [(1.0_rp / (exp((this%dos_energy_grid(ie) - this%fermi_level) / kT) + 1.0_rp), &
                               ie = 1, n_energy)]
                  m1 = trapezoidal_integral(this%dos_energy_grid, integrand) / m0
               else
                  m1 = 0.0_rp
               end if

               ! Calculate second moment (m2) = ∫ (E - m1)^2 * DOS(E) * f(E) dE / m0
               if (abs(m0) > 1.0e-12_rp) then
                  integrand = (this%dos_energy_grid - m1)**2 * this%projected_dos(isite, iorb, ispin, :) * &
                             [(1.0_rp / (exp((this%dos_energy_grid(ie) - this%fermi_level) / kT) + 1.0_rp), &
                               ie = 1, n_energy)]
                  m2 = trapezoidal_integral(this%dos_energy_grid, integrand) / m0
                  m2 = sqrt(max(m2, 0.0_rp))  ! Take square root to get width
               else
                  m2 = 0.0_rp
               end if

               ! Store moments
               this%band_moments(isite, iorb, ispin, 1) = m0  ! m0
               this%band_moments(isite, iorb, ispin, 2) = m1  ! m1
               this%band_moments(isite, iorb, ispin, 3) = m2  ! m2

               call g_logger%info('calculate_band_moments: Site ' // trim(int2str(isite)) // &
                                 ', Orbital ' // trim(int2str(iorb)) // ', Spin ' // trim(int2str(ispin)) // &
                                 ': m0= ' // trim(real2str(m0, '(F 8.5)')) // ', m1= ' // trim(real2str(m1, '(F 8.5)')) // &
                                 ', m2= ' // trim(real2str(m2, '(F 8.6)')), __FILE__, __LINE__)
            end do
         end do
      end do

      ! Clean up
      deallocate(integrand)

      call g_logger%info('calculate_band_moments: Band moments calculation completed', __FILE__, __LINE__)
   end subroutine calculate_band_moments

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Find Fermi level from calculated DOS by integrating to find electron count
   !---------------------------------------------------------------------------
   function find_fermi_level_from_dos(this, total_electrons) result(fermi_level)
      class(reciprocal), intent(in) :: this
      real(rp), intent(in) :: total_electrons
      real(rp) :: fermi_level

      ! Local variables
      integer :: ie, max_iter
      real(rp) :: integrated_dos, prev_integrated
      real(rp) :: energy, prev_energy, kT, fermi_weight
      real(rp) :: e_min, e_max, e_mid, electrons_at_e

      call g_logger%info('find_fermi_level_from_dos: Finding Fermi level for ' // &
                        trim(real2str(total_electrons, '(F 8.5)')) // ' electrons at T = ' // &
                        trim(real2str(this%temperature, '(F 8.5)')) // ' K', __FILE__, __LINE__)

      ! Check if DOS is calculated
      if (.not. allocated(this%total_dos)) then
         call g_logger%error('find_fermi_level_from_dos: Total DOS not calculated', __FILE__, __LINE__)
         fermi_level = 0.0_rp
         return
      end if

      ! Boltzmann constant in eV/K
      kT = this%temperature * 8.617333262145e-5_rp

      ! Use bisection method to find Fermi level
      e_min = this%dos_energy_grid(1)
      e_max = this%dos_energy_grid(this%n_energy_points)
      max_iter = 100

      do ie = 1, max_iter
         e_mid = (e_min + e_max) / 2.0_rp
         electrons_at_e = this%integrate_dos_up_to_energy(e_mid, kT)

         if (abs(electrons_at_e - total_electrons) < 1.0e-6_rp) then
            fermi_level = e_mid
            exit
         else if (electrons_at_e < total_electrons) then
            e_min = e_mid
         else
            e_max = e_mid
         end if
      end do

      fermi_level = e_mid

      ! Final check
      electrons_at_e = this%integrate_dos_up_to_energy(fermi_level, kT)
      call g_logger%info('find_fermi_level_from_dos: Found Fermi level at ' // &
                        trim(real2str(fermi_level, '(F 8.5)')) // ' Ry (integrated ' // &
                        trim(real2str(electrons_at_e, '(F 8.5)')) // ' electrons)', __FILE__, __LINE__)
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

      ! Local variables
      integer :: ie
      real(rp) :: e, fermi_weight, delta_e

      integral = 0.0_rp

      do ie = 1, this%n_energy_points - 1
         e = this%dos_energy_grid(ie)
         delta_e = this%dos_energy_grid(ie+1) - e

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

      call g_logger%info('write_dos_to_file: Writing DOS to ' // trim(filename), __FILE__, __LINE__)

      ! Write total DOS
      open(newunit=unit, file=trim(filename), status='replace', action='write')

      ! Write header
      write(unit, '(A)') '# Density of States'
      write(unit, '(A,A)') '# Method: ', trim(this%dos_method)
      if (trim(this%dos_method) == 'gaussian') then
         write(unit, '(A,F8.4)') '# Gaussian sigma: ', this%gaussian_sigma
      end if
      write(unit, '(A,I0)') '# Energy points: ', this%n_energy_points
      write(unit, '(A,2F8.3)') '# Energy range: ', this%dos_energy_range
      write(unit, '(A)') '# Energy (Ry)    Total DOS'

      ! Write DOS data
      do i_energy = 1, this%n_energy_points
         write(unit, '(2F12.6)') this%dos_energy_grid(i_energy), this%total_dos(i_energy)
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
            write(unit, '(A,F8.4)') '# Gaussian sigma: ', this%gaussian_sigma
         end if
         write(unit, '(A,I0)') '# Energy points: ', this%n_energy_points
         write(unit, '(A,2F8.3)') '# Energy range: ', this%dos_energy_range
         write(unit, '(A)') '# Columns: Energy, s_up, p_up, d_up, f_up, s_down, p_down, d_down, f_down'

         ! Write projected DOS data
         do i_energy = 1, this%n_energy_points
            write(unit, '(9F12.6)') this%dos_energy_grid(i_energy), &
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
            write(unit, '(A,F8.4)') '# Gaussian sigma: ', this%gaussian_sigma
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

end module reciprocal_mod