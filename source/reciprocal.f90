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
   use spglib_interface_mod, only: spglib_interface
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

      ! Eigenvalue/eigenvector storage
      !> Eigenvalues for k-mesh [nbands, nk_total]
      real(rp), dimension(:, :), allocatable :: eigenvalues
      !> Eigenvalues for k-path [nbands, nk_path]
      real(rp), dimension(:, :), allocatable :: eigenvalues_path
      !> Eigenvectors for k-mesh [max_orb_channels, nbands, nk_total]
      complex(rp), dimension(:, :, :), allocatable :: eigenvectors
      !> Eigenvectors for k-path [max_orb_channels, nbands, nk_path]
      complex(rp), dimension(:, :, :), allocatable :: eigenvectors_path

      ! Symmetry operations (spglib integration)
      !> spglib interface for crystallographic operations
      type(spglib_interface) :: spglib
      !> Space group number
      integer :: space_group_number
      !> Space group symbol
      character(len=20) :: space_group_symbol
      !> Crystal system (cubic, tetragonal, etc.)
      character(len=10) :: crystal_system

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
      procedure :: determine_crystal_structure
      procedure :: generate_kpath
      procedure :: get_high_symmetry_points
      procedure :: diagonalize_hamiltonian
      procedure :: diagonalize_at_kpoints
      procedure :: calculate_band_structure
      procedure :: restore_to_default
      ! Symmetry operations (spglib)
      procedure :: initialize_symmetry
      procedure :: get_symmetry_info
      procedure :: generate_reduced_kpoint_mesh
      procedure :: generate_symmetry_kpath
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
      call obj%initialize_symmetry()
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

      ! Initialize reciprocal lattice to zero
      this%reciprocal_vectors = 0.0_rp
      this%reciprocal_volume = 0.0_rp
   end subroutine restore_to_default

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
            if (allocated(this%lattice%sbarvec) .and. ineigh <= size(this%lattice%sbarvec, 2)) then
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
            k_dot_r = dot_product(k_vec, r_vec)
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

      if (debug_this_k) then
         call g_logger%info('fourier_transform_hamiltonian: Debug output for k-vector', __FILE__, __LINE__)
         if (allocated(this%lattice%sbarvec)) then
            call g_logger%info('fourier_transform_hamiltonian: sbarvec is allocated', __FILE__, __LINE__)
         else
            call g_logger%info('fourier_transform_hamiltonian: WARNING - sbarvec is NOT allocated!', __FILE__, __LINE__)
         end if
      end if

      ! Add simple check for ANY call during band structure calculation
      if (ntype == 1) then
         if (allocated(this%lattice%sbarvec)) then
            call g_logger%info('fourier_transform_hamiltonian: sbarvec allocated for band structure', __FILE__, __LINE__)
         else
            call g_logger%info('fourier_transform_hamiltonian: WARNING - sbarvec NOT allocated for band structure!', __FILE__, __LINE__)
         end if
      end if

      ! Allocate structure factors
      allocate(structure_factors(nr))

      ! Calculate structure factors
      call this%calculate_structure_factors(k_vec, ntype, structure_factors)
      
      ! Print first few structure factors for first k-point to debug
      if (ntype == 1 .and. size(structure_factors) >= 5) then
          write(*, '(A,2F10.6)') "First structure factor:", structure_factors(1)
          write(*, '(A,2F10.6)') "Second structure factor:", structure_factors(2)
          write(*, '(A,2F10.6)') "Third structure factor:", structure_factors(3)
          write(*, '(A,2F10.6)') "Fourth structure factor:", structure_factors(4)
          write(*, '(A,2F10.6)') "Fifth structure factor:", structure_factors(5)
      end if

      ! Initialize result
      hk_result = cmplx(0.0_rp, 0.0_rp, rp)

      ! Sum over neighbors: H(k) = Σ_R H(R) * exp(i*k·R)
      do ineigh = 1, nr
         h_neighbor = this%hamiltonian%ee(:, :, ineigh, ntype)
         hk_result = hk_result + h_neighbor * structure_factors(ineigh)
      end do

      if (debug_this_k) then
         call g_logger%info('fourier_transform_hamiltonian: Completed FT', __FILE__, __LINE__)
      end if

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
            if (ik <= 3) then
               write(debug_msg, '(A,I0,A,3F12.6,A,3F12.6,A)') 'build_kspace_hamiltonian: k-point ', ik, &
                     ' fractional=(', this%k_points(:, ik), ') cartesian=(', k_cartesian, ')'
               call g_logger%info(trim(debug_msg), __FILE__, __LINE__)
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
   !> Determine crystal structure from lattice vectors
   !---------------------------------------------------------------------------
   function determine_crystal_structure(this) result(crystal_type)
      class(reciprocal), intent(in) :: this
      character(len=10) :: crystal_type
      ! Local variables
      real(rp), dimension(3, 3) :: lattice_vectors
      real(rp) :: a, b, c, alpha, beta, gamma
      real(rp), parameter :: tol = 1.0e-6_rp
      
      ! Get the lattice vectors from lattice%a
      lattice_vectors = this%lattice%a
      
      ! Calculate lattice parameters
      a = norm2(lattice_vectors(:, 1))
      b = norm2(lattice_vectors(:, 2))
      c = norm2(lattice_vectors(:, 3))
      
      ! Calculate angles
      alpha = acos(dot_product(lattice_vectors(:, 2), lattice_vectors(:, 3)) / (b * c))
      beta = acos(dot_product(lattice_vectors(:, 1), lattice_vectors(:, 3)) / (a * c))
      gamma = acos(dot_product(lattice_vectors(:, 1), lattice_vectors(:, 2)) / (a * b))
      
      ! Convert to degrees
      alpha = alpha * 180.0_rp / acos(-1.0_rp)
      beta = beta * 180.0_rp / acos(-1.0_rp) 
      gamma = gamma * 180.0_rp / acos(-1.0_rp)
      
      call g_logger%info('determine_crystal_structure: Lattice parameters: a=' // &
         real2str(a) // ', b=' // real2str(b) // ', c=' // real2str(c), __FILE__, __LINE__)
      call g_logger%info('determine_crystal_structure: Angles: alpha=' // &
         real2str(alpha) // ', beta=' // real2str(beta) // ', gamma=' // real2str(gamma), __FILE__, __LINE__)
      
      ! Determine crystal system
      if (abs(a - b) < tol .and. abs(b - c) < tol .and. &
          abs(alpha - 90.0_rp) < tol .and. abs(beta - 90.0_rp) < tol .and. abs(gamma - 90.0_rp) < tol) then
         ! Cubic system - now determine if it's SC, BCC, or FCC
         ! This is a simplified approach - in practice, need to analyze the basis
         ! For typical LMTO calculations, BCC is common for Fe, Co, Cr
         crystal_type = 'bcc'
         call g_logger%info('determine_crystal_structure: Detected cubic system, assuming BCC', __FILE__, __LINE__)
      else if (abs(a - b) < tol .and. abs(alpha - 90.0_rp) < tol .and. &
               abs(beta - 90.0_rp) < tol .and. abs(gamma - 120.0_rp) < tol) then
         crystal_type = 'hexagonal'
         call g_logger%info('determine_crystal_structure: Detected hexagonal system', __FILE__, __LINE__)
      else if (abs(alpha - 90.0_rp) < tol .and. abs(beta - 90.0_rp) < tol .and. abs(gamma - 90.0_rp) < tol) then
         crystal_type = 'orthorhombic'
         call g_logger%info('determine_crystal_structure: Detected orthorhombic system', __FILE__, __LINE__)
      else
         crystal_type = 'bcc'  ! Default to BCC for LMTO
         call g_logger%warning('determine_crystal_structure: Unknown crystal system, defaulting to BCC', __FILE__, __LINE__)
      end if
      
   end function determine_crystal_structure

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Get high-symmetry points for different crystal systems
   !---------------------------------------------------------------------------
   subroutine get_high_symmetry_points(this, crystal_type, kpoints, labels)
      class(reciprocal), intent(in) :: this
      character(len=*), intent(in) :: crystal_type
      real(rp), dimension(:, :), allocatable, intent(out) :: kpoints
      character(len=10), dimension(:), allocatable, intent(out) :: labels
      ! Local variables
      integer :: npts

      ! For now, implement cubic (simple cubic, FCC, BCC) high-symmetry points
      ! This should be extended to use spglib for automatic detection
      select case (trim(crystal_type))
      case ('fcc', 'FCC')
         ! FCC high-symmetry points in units of 2π/a
         npts = 4
         allocate(kpoints(3, npts))
         allocate(labels(npts))
         
         ! Γ point
         kpoints(:, 1) = [0.0_rp, 0.0_rp, 0.0_rp]
         labels(1) = 'Γ'
         
         ! X point
         kpoints(:, 2) = [0.5_rp, 0.0_rp, 0.5_rp]
         labels(2) = 'X'
         
         ! L point  
         kpoints(:, 3) = [0.5_rp, 0.5_rp, 0.5_rp]
         labels(3) = 'L'
         
         ! W point
         kpoints(:, 4) = [0.5_rp, 0.25_rp, 0.75_rp]
         labels(4) = 'W'

      case ('bcc', 'BCC')
         ! BCC high-symmetry points
         npts = 4
         allocate(kpoints(3, npts))
         allocate(labels(npts))
         
         kpoints(:, 1) = [0.0_rp, 0.0_rp, 0.0_rp]    ! Γ
         labels(1) = 'Γ'
         kpoints(:, 2) = [0.0_rp, 0.5_rp, 0.5_rp]    ! H
         labels(2) = 'H'
         kpoints(:, 3) = [0.25_rp, 0.25_rp, 0.25_rp] ! P
         labels(3) = 'P'
         kpoints(:, 4) = [0.5_rp, 0.5_rp, 0.0_rp]    ! N
         labels(4) = 'N'

      case ('sc', 'simple_cubic', 'cubic')
         ! Simple cubic high-symmetry points
         npts = 4
         allocate(kpoints(3, npts))
         allocate(labels(npts))
         
         kpoints(:, 1) = [0.0_rp, 0.0_rp, 0.0_rp]    ! Γ
         labels(1) = 'Γ'
         kpoints(:, 2) = [0.5_rp, 0.0_rp, 0.0_rp]    ! X
         labels(2) = 'X'
         kpoints(:, 3) = [0.5_rp, 0.5_rp, 0.0_rp]    ! M
         labels(3) = 'M'
         kpoints(:, 4) = [0.5_rp, 0.5_rp, 0.5_rp]    ! R
         labels(4) = 'R'

      case default
         ! Default to simple cubic if unknown
         call g_logger%warning('get_high_symmetry_points: Unknown crystal type ' // trim(crystal_type) // ', using simple cubic', __FILE__, __LINE__)
         call get_high_symmetry_points(this, 'sc', kpoints, labels)
      end select

      call g_logger%info('get_high_symmetry_points: Generated ' // trim(int2str(npts)) // ' high-symmetry points for ' // trim(crystal_type), __FILE__, __LINE__)
   end subroutine get_high_symmetry_points

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Generate k-path through high-symmetry points
   !---------------------------------------------------------------------------
   subroutine generate_kpath(this, crystal_type, npts_per_segment, path_indices)
      class(reciprocal), intent(inout) :: this
      character(len=*), intent(in) :: crystal_type
      integer, intent(in) :: npts_per_segment
      integer, dimension(:), intent(in), optional :: path_indices
      ! Local variables
      real(rp), dimension(:, :), allocatable :: hs_kpoints
      character(len=10), dimension(:), allocatable :: hs_labels
      integer :: nhs, nseg, total_pts, i, j, k, idx
      integer, dimension(:), allocatable :: path
      real(rp) :: segment_length, total_distance
      real(rp), dimension(3) :: kvec_diff

      call g_logger%info('generate_kpath: Generating k-path for ' // trim(crystal_type), __FILE__, __LINE__)

      ! Get high-symmetry points
      call this%get_high_symmetry_points(crystal_type, hs_kpoints, hs_labels)
      nhs = size(hs_labels)

      ! Define default path if not provided
      if (present(path_indices)) then
         path = path_indices
      else
         ! Default path: connect all high-symmetry points in order
         allocate(path(nhs))
         do i = 1, nhs
            path(i) = i
         end do
      end if

      nseg = size(path) - 1
      total_pts = nseg * npts_per_segment + 1  ! +1 for the last point

      ! Allocate k-path arrays
      if (allocated(this%k_path)) deallocate(this%k_path)
      if (allocated(this%k_labels)) deallocate(this%k_labels)
      if (allocated(this%k_distances)) deallocate(this%k_distances)
      
      allocate(this%k_path(3, total_pts))
      allocate(this%k_labels(total_pts))
      allocate(this%k_distances(total_pts))

      ! Generate k-path
      idx = 1
      total_distance = 0.0_rp
      this%k_distances(1) = 0.0_rp

      do i = 1, nseg
         kvec_diff = hs_kpoints(:, path(i+1)) - hs_kpoints(:, path(i))
         segment_length = sqrt(sum(kvec_diff**2))
         
         do j = 1, npts_per_segment
            if (idx > total_pts) exit
            
            ! Linear interpolation
            this%k_path(:, idx) = hs_kpoints(:, path(i)) + &
                                  real(j-1, rp) / real(npts_per_segment-1, rp) * kvec_diff
            
            ! Set label for high-symmetry points
            if (j == 1) then
               this%k_labels(idx) = hs_labels(path(i))
            else
               this%k_labels(idx) = ''
            end if
            
            ! Calculate cumulative distance
            if (idx > 1) then
               total_distance = total_distance + sqrt(sum((this%k_path(:, idx) - this%k_path(:, idx-1))**2))
            end if
            this%k_distances(idx) = total_distance
            
            idx = idx + 1
         end do
      end do

      ! Add the final point
      if (idx <= total_pts) then
         this%k_path(:, idx) = hs_kpoints(:, path(nseg+1))
         this%k_labels(idx) = hs_labels(path(nseg+1))
         if (idx > 1) then
            total_distance = total_distance + sqrt(sum((this%k_path(:, idx) - this%k_path(:, idx-1))**2))
         end if
         this%k_distances(idx) = total_distance
      end if

      this%nk_path = total_pts
      
      call g_logger%info('generate_kpath: Generated k-path with ' // trim(int2str(total_pts)) // ' points', __FILE__, __LINE__)
      
      ! Clean up
      deallocate(hs_kpoints, hs_labels, path)
   end subroutine generate_kpath

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
      integer :: nk

      use_path = .false.
      if (present(use_kpath)) use_path = use_kpath

      if (use_path) then
         if (.not. allocated(this%k_path)) then
            call g_logger%error('diagonalize_hamiltonian: k-path not generated', __FILE__, __LINE__)
            return
         end if
         nk = this%nk_path
         call g_logger%info('diagonalize_hamiltonian: Using k-path with ' // trim(int2str(nk)) // ' points', __FILE__, __LINE__)
         call this%diagonalize_at_kpoints(ham, this%k_path, nk)
      else
         if (.not. allocated(this%k_points)) then
            call g_logger%error('diagonalize_hamiltonian: k-mesh not generated', __FILE__, __LINE__)
            return
         end if
         nk = this%nk_total
         call g_logger%info('diagonalize_hamiltonian: Using k-mesh with ' // trim(int2str(nk)) // ' points', __FILE__, __LINE__)
         call this%diagonalize_at_kpoints(ham, this%k_points, nk)
      end if
   end subroutine diagonalize_hamiltonian

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Unified diagonalization method for any set of k-points
   !---------------------------------------------------------------------------
   subroutine diagonalize_at_kpoints(this, ham, k_vectors, nk_points)
      class(reciprocal), intent(inout) :: this
      class(hamiltonian), intent(in) :: ham
      real(rp), dimension(:, :), intent(in) :: k_vectors  ! (3, nk_points)
      integer, intent(in) :: nk_points
      ! Local variables
      integer :: i, nmat, lwork, info
      complex(rp), dimension(:, :), allocatable :: h_k, h_k_copy
      real(rp), dimension(:), allocatable :: eigenvals
      complex(rp), dimension(:, :), allocatable :: eigenvecs
      complex(rp), dimension(:), allocatable :: work_complex
      real(rp), dimension(:), allocatable :: rwork
      character(len=100) :: info_msg

      call g_logger%info('diagonalize_at_kpoints: Starting diagonalization of ' // trim(int2str(nk_points)) // ' k-points', __FILE__, __LINE__)

      ! Get matrix size from the actual Hamiltonian
      if (allocated(ham%ee)) then
         nmat = size(ham%ee, 1)  ! Use actual Hamiltonian dimension
      else
         nmat = this%max_orb_channels  ! Fallback to max orbital channels
      end if

      call g_logger%info('diagonalize_at_kpoints: Matrix size is ' // trim(int2str(nmat)) // 'x' // trim(int2str(nmat)), __FILE__, __LINE__)

      ! Allocate eigenvalue and eigenvector arrays
      if (allocated(this%eigenvalues)) deallocate(this%eigenvalues)
      if (allocated(this%eigenvectors)) deallocate(this%eigenvectors)
      
      allocate(this%eigenvalues(nmat, nk_points))
      allocate(this%eigenvectors(nmat, nmat, nk_points))

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
      do i = 1, nk_points
         ! Build Hamiltonian at this k-point using unified Fourier transform
         call this%fourier_transform_hamiltonian(k_vectors(:, i), 1, h_k)

         ! Make a copy for diagonalization (LAPACK destroys input)
         h_k_copy = h_k

         ! Diagonalize using LAPACK ZHEEV
         call zheev('V', 'U', nmat, h_k_copy, nmat, eigenvals, work_complex, lwork, rwork, info)
         
         if (info /= 0) then
            write(info_msg, '(A,I0,A,I0)') 'diagonalize_at_kpoints: ZHEEV failed for k-point ', i, ' with info=', info
            call g_logger%error(trim(info_msg), __FILE__, __LINE__)
            cycle
         end if

         ! Store results
         this%eigenvalues(:, i) = eigenvals
         this%eigenvectors(:, :, i) = h_k_copy
      end do

      call g_logger%info('diagonalize_at_kpoints: Completed diagonalization', __FILE__, __LINE__)

      ! Clean up
      deallocate(h_k, h_k_copy, eigenvals, eigenvecs)
      deallocate(work_complex, rwork)
   end subroutine diagonalize_at_kpoints

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

      ! Generate k-path
      call this%generate_kpath(crystal_type, npts_per_segment)

      ! Diagonalize Hamiltonian along k-path using unified method
      call this%diagonalize_at_kpoints(ham, this%k_path, this%nk_path)

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
   !> Initialize spglib interface with crystal structure information
   !---------------------------------------------------------------------------
   subroutine initialize_symmetry(this)
      class(reciprocal), intent(inout) :: this
      real(rp), dimension(3,3) :: lattice_vectors
      real(rp), dimension(:,:), allocatable :: atomic_positions
      integer, dimension(:), allocatable :: atomic_types
      integer :: i

      call g_logger%info('initialize_symmetry: Setting up spglib interface', __FILE__, __LINE__)

      ! Check if spglib is available
      if (.not. this%spglib%is_available()) then
         call g_logger%warning('initialize_symmetry: spglib not available, using default symmetry', __FILE__, __LINE__)
         this%space_group_number = 0
         this%space_group_symbol = 'P1'
         this%crystal_system = 'unknown'
         return
      end if

      ! Get lattice vectors (convert from lattice object)
      do i = 1, 3
         lattice_vectors(i,:) = this%lattice%a(i,:)
      end do

      ! Get atomic positions and types from lattice
      allocate(atomic_positions(3, this%lattice%nrec))
      allocate(atomic_types(this%lattice%nrec))

      do i = 1, this%lattice%nrec
         atomic_positions(:,i) = this%lattice%cr(:,i)  ! Fractional coordinates
         atomic_types(i) = this%lattice%iz(i)          ! Atomic number
      end do

      ! Initialize spglib with crystal structure
      call this%spglib%initialize(lattice_vectors, atomic_positions, atomic_types)

      ! Get symmetry information
      call this%get_symmetry_info()

      ! Clean up temporary arrays
      deallocate(atomic_positions, atomic_types)
   end subroutine initialize_symmetry

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Get space group and crystal system information from spglib
   !---------------------------------------------------------------------------
   subroutine get_symmetry_info(this)
      class(reciprocal), intent(inout) :: this
      integer :: num_symmetries

      if (.not. this%spglib%is_available()) then
         return
      end if

      ! Get space group number and symbol (already done in spglib initialize)
      this%space_group_number = this%spglib%get_space_group_number()
      this%space_group_symbol = this%spglib%get_space_group_symbol()
      this%crystal_system = this%spglib%get_crystal_system_name()

      ! Get number of symmetry operations
      num_symmetries = this%spglib%get_symmetry_operations()

      call g_logger%info('get_symmetry_info: Crystal system: ' // trim(this%crystal_system), __FILE__, __LINE__)
      call g_logger%info('get_symmetry_info: Space group: ' // trim(this%space_group_symbol) // &
                        ' (#' // trim(int2str(this%space_group_number)) // ')', __FILE__, __LINE__)
      call g_logger%info('get_symmetry_info: Number of symmetry operations: ' // &
                        trim(int2str(num_symmetries)), __FILE__, __LINE__)
   end subroutine get_symmetry_info

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

      if (.not. this%spglib%is_available()) then
         call g_logger%warning('generate_reduced_kpoint_mesh: spglib not available, using full mesh', __FILE__, __LINE__)
         call this%generate_mp_mesh()  ! Fall back to regular MP mesh
         return
      end if

      ! Get irreducible k-points using spglib
      num_ir_kpoints = this%spglib%get_reduced_kpoint_mesh(mesh_dims, shift)

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
   !> Generate high-symmetry k-path using spglib crystal system information
   !---------------------------------------------------------------------------
   subroutine generate_symmetry_kpath(this, nk_path_in)
      class(reciprocal), intent(inout) :: this
      integer, intent(in), optional :: nk_path_in
      logical :: path_available

      if (.not. this%spglib%is_available()) then
         call g_logger%warning('generate_symmetry_kpath: spglib not available, using default k-path', __FILE__, __LINE__)
         call this%generate_kpath('cubic', 50)  ! Fall back to regular k-path
         return
      end if

      ! Check if band path generation is available for this crystal system
      path_available = this%spglib%get_band_path_points()

      if (path_available) then
         call g_logger%info('generate_symmetry_kpath: Using spglib-based k-path for ' // &
                           trim(this%crystal_system) // ' crystal system', __FILE__, __LINE__)
         
         ! TODO: Implement actual spglib-based k-path generation
         ! For now, fall back to regular k-path with enhanced system-specific logic
         call this%generate_kpath(this%crystal_system, 50)
      else
         call g_logger%warning('generate_symmetry_kpath: spglib k-path not implemented, using default', __FILE__, __LINE__)
         call this%generate_kpath('cubic', 50)
      end if
   end subroutine generate_symmetry_kpath

end module reciprocal_mod