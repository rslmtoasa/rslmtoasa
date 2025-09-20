!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief
!> Symmetry analysis module for crystallographic operations and k-path generation
!> 
!> This module provides symmetry analysis capabilities including:
!> - Crystal structure determination from lattice vectors
!> - High-symmetry point generation for various crystal systems
!> - K-path generation for band structure calculations
!> - Integration with spglib for space group operations
!>
!> @author
!> The RS-LMTO-ASA team
!---------------------------------------------------------------------------
module symmetry_mod
   use precision_mod
   use logger_mod, only: g_logger
   use string_mod, only: int2str, real2str
   use lattice_mod, only: lattice
#ifdef USE_SPGLIB
   use spglib_interface_mod, only: spglib_interface
#endif

   implicit none
   private

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Type for symmetry analysis and k-path generation
   !---------------------------------------------------------------------------
   type, public :: symmetry
      ! Lattice reference for crystallographic analysis
      type(lattice), pointer :: lattice => null()
      
      ! Crystal structure information
      character(len=32) :: crystal_system = ''
      character(len=32) :: space_group_symbol = ''
      integer :: space_group_number = 0
      integer :: num_symmetry_operations = 0
      
      ! K-path related arrays
      real(rp), dimension(:, :), allocatable :: k_path
      character(len=10), dimension(:), allocatable :: k_labels
      real(rp), dimension(:), allocatable :: k_distances
      integer :: nk_path = 0
      
      ! Control flags
      logical :: suppress_internal_logs = .true.
      
#ifdef USE_SPGLIB
      ! spglib interface for crystallographic operations
      type(spglib_interface) :: spglib
#endif

   contains
      ! Initialization and cleanup
      procedure :: initialize => initialize_symmetry
      procedure :: finalize => finalize_symmetry
      
      ! Crystal structure analysis
      procedure :: determine_crystal_structure
      procedure :: get_symmetry_info
      
      ! K-path generation
      procedure :: get_high_symmetry_points
      procedure :: generate_kpath
      procedure :: generate_symmetry_kpath
      procedure :: generate_canonical_kpath
      procedure :: get_canonical_kpath_for_spacegroup
      
      ! Utility procedures
      procedure :: set_lattice_reference
   end type symmetry

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Initialize symmetry analysis
   !---------------------------------------------------------------------------
   subroutine initialize_symmetry(this, lattice_ref)
      class(symmetry), intent(inout) :: this
      type(lattice), target, intent(in) :: lattice_ref
      
      ! Set lattice reference
      this%lattice => lattice_ref
      
#ifdef USE_SPGLIB
      ! Initialize spglib interface
      block
         real(rp), dimension(3,3) :: lattice_vectors
         real(rp), dimension(:,:), allocatable :: atomic_positions
         integer, dimension(:), allocatable :: atomic_types
         integer :: i
         
         ! Get lattice vectors (convert from units of alat to absolute Cartesian)
         do i = 1, 3
            lattice_vectors(i,:) = lattice_ref%a(i,:) * lattice_ref%alat
         end do
         
         ! Get atomic positions and types from lattice
         allocate(atomic_positions(3, lattice_ref%nrec))
         allocate(atomic_types(lattice_ref%nrec))
         
         do i = 1, lattice_ref%nrec
            atomic_positions(:,i) = lattice_ref%cr(:,i)  ! Fractional coordinates
            atomic_types(i) = lattice_ref%iz(i)          ! Atomic number
         end do
         
         call this%spglib%initialize(lattice_vectors, atomic_positions, atomic_types)
         
         ! Clean up temporary arrays
         deallocate(atomic_positions, atomic_types)
      end block
      
      if (this%spglib%is_available()) then
         ! Get space group information
         call this%get_symmetry_info()
         call g_logger%info('initialize_symmetry: Detected space group ' // &
                           trim(this%space_group_symbol) // ' ( #' // &
                           trim(int2str(this%space_group_number)) // ' ), crystal system: ' // &
                           trim(this%crystal_system), __FILE__, __LINE__)
      else
         call g_logger%warning('initialize_symmetry: spglib not available', __FILE__, __LINE__)
      end if
#else
      call g_logger%info('initialize_symmetry: Built without spglib support', __FILE__, __LINE__)
#endif
      
   end subroutine initialize_symmetry

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Finalize symmetry analysis
   !---------------------------------------------------------------------------
   subroutine finalize_symmetry(this)
      class(symmetry), intent(inout) :: this
      
      ! Deallocate arrays
      if (allocated(this%k_path)) deallocate(this%k_path)
      if (allocated(this%k_labels)) deallocate(this%k_labels)
      if (allocated(this%k_distances)) deallocate(this%k_distances)
      
      ! Reset values
      this%crystal_system = ''
      this%space_group_symbol = ''
      this%space_group_number = 0
      this%num_symmetry_operations = 0
      this%nk_path = 0
      
      ! Clear lattice reference
      this%lattice => null()
      
   end subroutine finalize_symmetry

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Set lattice reference
   !---------------------------------------------------------------------------
   subroutine set_lattice_reference(this, lattice_ref)
      class(symmetry), intent(inout) :: this
      type(lattice), target, intent(in) :: lattice_ref
      
      this%lattice => lattice_ref
      
   end subroutine set_lattice_reference

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Get symmetry information from spglib
   !---------------------------------------------------------------------------
   subroutine get_symmetry_info(this)
      class(symmetry), intent(inout) :: this
      
#ifdef USE_SPGLIB
      if (this%spglib%is_available()) then
         ! Get space group information
         this%space_group_number = this%spglib%get_space_group_number()
         this%space_group_symbol = this%spglib%get_space_group_symbol()
         this%crystal_system = this%spglib%get_crystal_system_name()
         this%num_symmetry_operations = this%spglib%get_symmetry_operations()
         
         if (.not. this%suppress_internal_logs) then
            call g_logger%info('get_symmetry_info: Crystal system: ' // trim(this%crystal_system), __FILE__, __LINE__)
            call g_logger%info('get_symmetry_info: Space group: ' // trim(this%space_group_symbol) // &
                              ' (#' // trim(int2str(this%space_group_number)) // ')', __FILE__, __LINE__)
            call g_logger%info('get_symmetry_info: Number of symmetry operations: ' // &
                              trim(int2str(this%num_symmetry_operations)), __FILE__, __LINE__)
         end if
      end if
#endif
      
   end subroutine get_symmetry_info

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Determine crystal structure from lattice vectors
   !---------------------------------------------------------------------------
   function determine_crystal_structure(this) result(crystal_type)
      class(symmetry), intent(in) :: this
      character(len=10) :: crystal_type
      ! Local variables
      real(rp), dimension(3, 3) :: lattice_vectors
      real(rp) :: a, b, c, alpha, beta, gamma
      real(rp) :: G11, G22, G33, G12, G13, G23
      real(rp) :: avg_diag, avg_off, maxd, maxo
      real(rp), parameter :: tol_rel = 1.0e-5_rp
      
      if (.not. associated(this%lattice)) then
         call g_logger%error('determine_crystal_structure: lattice not associated', __FILE__, __LINE__)
         crystal_type = 'unknown'
         return
      end if
      
      ! Get the lattice vectors from lattice%a (primitive vectors in units of alat)
      lattice_vectors = this%lattice%a

#ifdef USE_SPGLIB
      ! If spglib already ran during initialization, prefer its result
      if (this%spglib%is_available()) then
         ! space_group_symbol typically contains centering like 'Im-3m' or 'Fm-3m'
         if (index(this%space_group_symbol, 'I') /= 0) then
            crystal_type = 'bcc'
            return
         else if (index(this%space_group_symbol, 'F') /= 0) then
            crystal_type = 'fcc'
            return
         else
            crystal_type = 'cubic'
            return
         end if
      end if
#endif

      ! Calculate lattice parameter norms
      a = norm2(lattice_vectors(:, 1))
      b = norm2(lattice_vectors(:, 2))
      c = norm2(lattice_vectors(:, 3))

      ! Build metric elements G_ij = a_i · a_j
      G11 = dot_product(lattice_vectors(:, 1), lattice_vectors(:, 1))
      G22 = dot_product(lattice_vectors(:, 2), lattice_vectors(:, 2))
      G33 = dot_product(lattice_vectors(:, 3), lattice_vectors(:, 3))
      G12 = dot_product(lattice_vectors(:, 1), lattice_vectors(:, 2))
      G13 = dot_product(lattice_vectors(:, 1), lattice_vectors(:, 3))
      G23 = dot_product(lattice_vectors(:, 2), lattice_vectors(:, 3))

      ! Quick logs for diagnostics (suppressed unless internal logs enabled)
      if (.not. this%suppress_internal_logs) then
         call g_logger%info('determine_crystal_structure: G11=' // real2str(G11) // ', G22=' // real2str(G22) // &
                           ', G33=' // real2str(G33), __FILE__, __LINE__)
         call g_logger%info('determine_crystal_structure: G12=' // real2str(G12) // &
                           ', G13=' // real2str(G13) // ', G23=' // real2str(G23), __FILE__, __LINE__)
      end if

      ! Check for cubic by testing equality of diagonal elements and equality of off-diagonals
      avg_diag = (G11 + G22 + G33) / 3.0_rp
      avg_off = (G12 + G13 + G23) / 3.0_rp
      maxd = max(abs(G11 - avg_diag), abs(G22 - avg_diag), abs(G33 - avg_diag))
      maxo = max(abs(G12 - avg_off), abs(G13 - avg_off), abs(G23 - avg_off))

      if (abs(avg_diag) > 0.0_rp .and. maxd / abs(avg_diag) < tol_rel .and. &
          (abs(avg_off) > 0.0_rp .and. maxo / abs(avg_off) < tol_rel .or. abs(maxo) < tol_rel)) then
         ! It's a cubic lattice (primitive vectors may be non-orthogonal)
         ! Distinguish BCC vs FCC by sign of off-diagonal metric elements in common primitive choices
         if (avg_off < 0.0_rp) then
            crystal_type = 'bcc'
            call g_logger%info('determine_crystal_structure: Detected cubic primitive cell consistent with BCC', __FILE__, __LINE__)
         else
            crystal_type = 'fcc'
            call g_logger%info('determine_crystal_structure: Detected cubic primitive cell consistent with FCC', __FILE__, __LINE__)
         end if
      else
         ! Fall back to previous angle-based checks for hexagonal/orthorhombic
         alpha = acos(G23 / (b * c))
         beta  = acos(G13 / (a * c))
         gamma = acos(G12 / (a * b))
         alpha = alpha * 180.0_rp / acos(-1.0_rp)
         beta  = beta  * 180.0_rp / acos(-1.0_rp)
         gamma = gamma * 180.0_rp / acos(-1.0_rp)

         if (abs(a - b) < 1.0e-6_rp .and. abs(alpha - 90.0_rp) < 1.0e-6_rp .and. &
             abs(beta - 90.0_rp) < 1.0e-6_rp .and. abs(gamma - 120.0_rp) < 1.0e-6_rp) then
            crystal_type = 'hexagonal'
            call g_logger%info('determine_crystal_structure: Detected hexagonal system', __FILE__, __LINE__)
         else if (abs(alpha - 90.0_rp) < 1.0e-6_rp .and. abs(beta - 90.0_rp) < 1.0e-6_rp .and. &
                  abs(gamma - 90.0_rp) < 1.0e-6_rp) then
            crystal_type = 'orthorhombic'
            call g_logger%info('determine_crystal_structure: Detected orthorhombic system', __FILE__, __LINE__)
         else
            crystal_type = 'bcc'  ! Default to BCC for LMTO
            call g_logger%warning('determine_crystal_structure: Unknown crystal system, defaulting to BCC', __FILE__, __LINE__)
         end if
      end if
      
   end function determine_crystal_structure

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Get high-symmetry points for different crystal systems
   !---------------------------------------------------------------------------
   subroutine get_high_symmetry_points(this, crystal_type, kpoints, labels)
      class(symmetry), intent(in) :: this
      character(len=*), intent(in) :: crystal_type
      real(rp), dimension(:, :), allocatable, intent(out) :: kpoints
      character(len=10), dimension(:), allocatable, intent(out) :: labels
      ! Local variables
      integer :: npts

      ! For now, implement cubic (simple cubic, FCC, BCC) high-symmetry points
      ! This should be extended to use spglib for automatic detection
      select case (trim(crystal_type))
      case ('tetragonal', 'tetr', 'TET')
         ! Tetragonal high-symmetry points (conventional)
         npts = 5
         allocate(kpoints(3, npts))
         allocate(labels(npts))
         kpoints(:,1) = [0.0_rp, 0.0_rp, 0.0_rp]; labels(1) = 'Γ'
         kpoints(:,2) = [0.5_rp, 0.0_rp, 0.0_rp]; labels(2) = 'X'
         kpoints(:,3) = [0.0_rp, 0.5_rp, 0.0_rp]; labels(3) = 'M'
         kpoints(:,4) = [0.0_rp, 0.0_rp, 0.5_rp]; labels(4) = 'Z'
         kpoints(:,5) = [0.5_rp, 0.0_rp, 0.5_rp]; labels(5) = 'R'

      case ('hexagonal', 'trigonal', 'hex')
         ! Hexagonal/trigonal high-symmetry points (conventional)
         npts = 5
         allocate(kpoints(3, npts))
         allocate(labels(npts))
         kpoints(:,1) = [0.0_rp, 0.0_rp, 0.0_rp]; labels(1) = 'Γ'
         kpoints(:,2) = [1.0_rp/3.0_rp, 2.0_rp/3.0_rp, 0.0_rp]; labels(2) = 'K'
         kpoints(:,3) = [0.5_rp, 0.0_rp, 0.0_rp]; labels(3) = 'M'
         kpoints(:,4) = [0.0_rp, 0.0_rp, 0.5_rp]; labels(4) = 'A'
         kpoints(:,5) = [1.0_rp/3.0_rp, 2.0_rp/3.0_rp, 0.5_rp]; labels(5) = 'H'

      case ('orthorhombic', 'ortho')
         ! Orthorhombic standard points
         npts = 5
         allocate(kpoints(3, npts))
         allocate(labels(npts))
         kpoints(:,1) = [0.0_rp, 0.0_rp, 0.0_rp]; labels(1) = 'Γ'
         kpoints(:,2) = [0.5_rp, 0.0_rp, 0.0_rp]; labels(2) = 'X'
         kpoints(:,3) = [0.0_rp, 0.5_rp, 0.0_rp]; labels(3) = 'Y'
         kpoints(:,4) = [0.0_rp, 0.0_rp, 0.5_rp]; labels(4) = 'Z'
         kpoints(:,5) = [0.5_rp, 0.5_rp, 0.5_rp]; labels(5) = 'T'

      case ('monoclinic', 'mono')
         ! Monoclinic simple path (very basic)
         npts = 4
         allocate(kpoints(3, npts))
         allocate(labels(npts))
         kpoints(:,1) = [0.0_rp, 0.0_rp, 0.0_rp]; labels(1) = 'Γ'
         kpoints(:,2) = [0.5_rp, 0.0_rp, 0.0_rp]; labels(2) = 'A'
         kpoints(:,3) = [0.0_rp, 0.5_rp, 0.0_rp]; labels(3) = 'B'
         kpoints(:,4) = [0.0_rp, 0.0_rp, 0.5_rp]; labels(4) = 'C'

      case ('triclinic', 'tric')
         ! Triclinic: only Gamma is well-defined
         npts = 1
         allocate(kpoints(3, npts))
         allocate(labels(npts))
         kpoints(:,1) = [0.0_rp, 0.0_rp, 0.0_rp]; labels(1) = 'Γ'
         
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
         ! BCC high-symmetry points - canonical path Γ-H-N-Γ-P-H|P-N
         npts = 6
         allocate(kpoints(3, npts))
         allocate(labels(npts))
         
         kpoints(:, 1) = [0.0_rp, 0.0_rp, 0.0_rp]    ! Γ
         labels(1) = 'Γ'
         kpoints(:, 2) = [0.5_rp, -0.5_rp, 0.5_rp]   ! H
         labels(2) = 'H'
         kpoints(:, 3) = [0.5_rp, 0.5_rp, 0.0_rp]    ! N
         labels(3) = 'N'
         kpoints(:, 4) = [0.0_rp, 0.0_rp, 0.0_rp]    ! Γ (back)
         labels(4) = 'Γ'
         kpoints(:, 5) = [0.25_rp, 0.25_rp, 0.25_rp] ! P
         labels(5) = 'P'
         kpoints(:, 6) = [0.5_rp, -0.5_rp, 0.5_rp]   ! H (back)
         labels(6) = 'H'

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
      class(symmetry), intent(inout) :: this
      character(len=*), intent(in) :: crystal_type
      integer, intent(in) :: npts_per_segment
      integer, dimension(:), intent(in), optional :: path_indices
      ! Local variables
      real(rp), dimension(:, :), allocatable :: hs_kpoints
      character(len=10), dimension(:), allocatable :: hs_labels
      integer :: nhs, nseg, total_pts, i, j, k, idx
      character(len=256) :: hs_list
      integer :: il
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

      ! Generate the k-path
      idx = 1
      total_distance = 0.0_rp
      
      call g_logger%info('generate_kpath: Starting k-path generation with ' // &
                        trim(int2str(nseg)) // ' segments and ' // &
                        trim(int2str(this%nk_path)) // ' total points allocated', __FILE__, __LINE__)
      
      do i = 1, nseg
         ! Calculate segment length
         kvec_diff = hs_kpoints(:, path(i+1)) - hs_kpoints(:, path(i))
         segment_length = norm2(kvec_diff)
         
         call g_logger%info('generate_kpath: Segment ' // trim(int2str(i)) // ': ' // &
                           trim(hs_labels(path(i))) // ' -> ' // trim(hs_labels(path(i+1))) // &
                           ', length=' // real2str(segment_length), __FILE__, __LINE__)
         
         ! Generate points along this segment
         do j = 1, npts_per_segment
            this%k_path(:, idx) = hs_kpoints(:, path(i)) + &
                                  real(j-1, rp) / real(npts_per_segment, rp) * kvec_diff
            this%k_distances(idx) = total_distance + &
                                   real(j-1, rp) / real(npts_per_segment, rp) * segment_length
            this%k_labels(idx) = ''  ! Most points don't have labels
            
            ! Log first few k-points for verification
            if (idx <= 3 .or. mod(idx, 25) == 0) then
               call g_logger%info('generate_kpath: k-point ' // trim(int2str(idx)) // ': [' // &
                                 real2str(this%k_path(1,idx)) // ', ' // &
                                 real2str(this%k_path(2,idx)) // ', ' // &
                                 real2str(this%k_path(3,idx)) // ']', __FILE__, __LINE__)
            end if
            
            idx = idx + 1
         end do
         
         total_distance = total_distance + segment_length
      end do
      
      ! Add the final point
      this%k_path(:, total_pts) = hs_kpoints(:, path(nseg+1))
      this%k_distances(total_pts) = total_distance
      this%k_labels(total_pts) = ''
      
      ! Set labels for high-symmetry points
      idx = 1
      do i = 1, nseg + 1
         this%k_labels(idx) = hs_labels(path(i))
         if (i <= nseg) idx = idx + npts_per_segment
      end do
      
      this%nk_path = total_pts

      ! Log the symmetry point sequence
      hs_list = ''
      il = 0
      do i = 1, size(path)
         if (il > 0) then
            hs_list = hs_list(1:il) // ' -> '
            il = il + 4
         end if
         if (il + len_trim(hs_labels(path(i))) <= len(hs_list)) then
            hs_list = hs_list(1:il) // trim(hs_labels(path(i)))
            il = il + len_trim(hs_labels(path(i)))
         else
            exit  ! String too long, truncate
         end if
      end do
      
      if (il > 0) then
         call g_logger%info('generate_kpath: Symmetry points for path: ' // hs_list(1:il), __FILE__, __LINE__)
      else
         call g_logger%info('generate_kpath: No labelled symmetry points found', __FILE__, __LINE__)
      end if

      call g_logger%info('generate_kpath: Generated k-path with ' // trim(int2str(this%nk_path)) // ' points', __FILE__, __LINE__)
   end subroutine generate_kpath

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Generate symmetry-based k-path using spglib when available
   !---------------------------------------------------------------------------
   subroutine generate_symmetry_kpath(this, nk_path_in)
      class(symmetry), intent(inout) :: this
      integer, intent(in), optional :: nk_path_in

#ifdef USE_SPGLIB
      if (.not. this%spglib%is_available()) then
         call g_logger%info('generate_symmetry_kpath: spglib NOT available — using default k-path', __FILE__, __LINE__)
         call this%generate_kpath('cubic', 50)
         return
      end if

      ! Prefer spglib-detected crystal system but fall back when empty
      if (len_trim(this%crystal_system) == 0) then
         call g_logger%warning('generate_symmetry_kpath: spglib available but crystal_system empty, using cubic', __FILE__, __LINE__)
         call this%generate_kpath('cubic', 50)
         return
      end if

      call g_logger%info('generate_symmetry_kpath: Generating k-path for detected crystal system: ' // trim(this%crystal_system), __FILE__, __LINE__)
      call this%get_high_symmetry_points(this%crystal_system, this%k_path, this%k_labels)

      ! If get_high_symmetry_points returned a trivial set, fall back to default k-path construction
      if (.not. allocated(this%k_path) .or. size(this%k_path,2) == 0) then
         call g_logger%warning('generate_symmetry_kpath: Failed to get HS points, falling back to default k-path', __FILE__, __LINE__)
         call this%generate_kpath('cubic', 50)
      else
         ! Build full path by sampling each successive HS point segment
         ! For now sample 50 points per segment (configurable later)
         call this%generate_kpath(this%crystal_system, 50)
      end if
#else
      call g_logger%info('generate_symmetry_kpath: Built without spglib support, using default k-path', __FILE__, __LINE__)
      call this%generate_kpath('cubic', 50)
#endif
      
   end subroutine generate_symmetry_kpath

   !---------------------------------------------------------------------------
   !> @brief Generate canonical k-path using spglib dataset API and standard band paths
   !> This method uses the detailed space group information to select appropriate
   !> k-paths following crystallographic conventions
   !---------------------------------------------------------------------------
   subroutine generate_canonical_kpath(this, npts_per_segment)
      class(symmetry), intent(inout) :: this
      integer, intent(in), optional :: npts_per_segment
      
      integer :: npts, spg_number
      character(len=20) :: international_symbol
      character(len=100) :: hall_symbol
      logical :: dataset_success
      
      npts = 50
      if (present(npts_per_segment)) npts = npts_per_segment

      call g_logger%info('generate_canonical_kpath: Starting canonical k-path generation with ' // &
                        trim(int2str(npts)) // ' points per segment', __FILE__, __LINE__)

#ifdef USE_SPGLIB
      if (.not. this%spglib%is_available()) then
         call g_logger%warning('generate_canonical_kpath: spglib not available, falling back to basic k-path', __FILE__, __LINE__)
         call this%generate_kpath('cubic', npts)
         call g_logger%info('generate_canonical_kpath: FALLBACK used - basic k-path generated', __FILE__, __LINE__)
         return
      end if

      call g_logger%info('generate_canonical_kpath: spglib is available, attempting dataset retrieval', __FILE__, __LINE__)

      ! Get detailed dataset information
      dataset_success = this%spglib%get_dataset(space_group_number=spg_number, &
                                               international=international_symbol, &
                                               hall=hall_symbol)

      call g_logger%info('generate_canonical_kpath: Dataset retrieval success: ' // &
                        merge('TRUE ', 'FALSE', dataset_success), __FILE__, __LINE__)

      if (.not. dataset_success) then
         call g_logger%warning('generate_canonical_kpath: Failed to get spglib dataset, using default k-path', __FILE__, __LINE__)
         call this%generate_kpath(this%crystal_system, npts)
         call g_logger%info('generate_canonical_kpath: FALLBACK used - crystal system k-path generated', __FILE__, __LINE__)
         return
      end if

      call g_logger%info('generate_canonical_kpath: SUCCESS - Space group detected as #' // &
                        trim(int2str(spg_number)) // ' (' // trim(international_symbol) // ')', __FILE__, __LINE__)
      call g_logger%info('generate_canonical_kpath: Using crystal family-based k-path (not true spglib canonical)', __FILE__, __LINE__)

      ! Generate canonical k-path based on space group
      call this%get_canonical_kpath_for_spacegroup(spg_number, international_symbol, npts)
      call g_logger%info('generate_canonical_kpath: CRYSTAL FAMILY-BASED k-path generation completed', __FILE__, __LINE__)
#else
      call g_logger%info('generate_canonical_kpath: Built without spglib support, using default k-path', __FILE__, __LINE__)
      call this%generate_kpath('cubic', npts)
      call g_logger%info('generate_canonical_kpath: NO-SPGLIB fallback used', __FILE__, __LINE__)
#endif

   end subroutine generate_canonical_kpath

   !---------------------------------------------------------------------------
   !> @brief Get canonical k-path for specific space group number
   !> Implements standard high-symmetry paths for different crystal systems
   !---------------------------------------------------------------------------
   subroutine get_canonical_kpath_for_spacegroup(this, spg_number, international, npts_per_segment)
      class(symmetry), intent(inout) :: this
      integer, intent(in) :: spg_number
      character(len=*), intent(in) :: international
      integer, intent(in) :: npts_per_segment
      
      character(len=20) :: crystal_family
      
      ! Determine crystal family from space group number
      if (spg_number >= 1 .and. spg_number <= 2) then
         crystal_family = 'triclinic'
      else if (spg_number >= 3 .and. spg_number <= 15) then
         crystal_family = 'monoclinic'
      else if (spg_number >= 16 .and. spg_number <= 74) then
         crystal_family = 'orthorhombic'
      else if (spg_number >= 75 .and. spg_number <= 142) then
         crystal_family = 'tetragonal'
      else if (spg_number >= 143 .and. spg_number <= 167) then
         crystal_family = 'trigonal'
      else if (spg_number >= 168 .and. spg_number <= 194) then
         crystal_family = 'hexagonal'
      else if (spg_number >= 195 .and. spg_number <= 230) then
         crystal_family = 'cubic'
      else
         crystal_family = 'unknown'
      end if

      call g_logger%info('get_canonical_kpath_for_spacegroup: Crystal family: ' // trim(crystal_family) // &
                        ' for space group ' // trim(int2str(spg_number)), __FILE__, __LINE__)

      ! For now, use the existing crystal type detection
      ! A full implementation would need detailed k-path tables for each space group
      select case (trim(crystal_family))
      case ('cubic')
         ! For cubic systems, distinguish between primitive, FCC, and BCC based on centering
         if (index(international, 'P') == 1) then
            call this%generate_kpath('cubic', npts_per_segment)
         else if (index(international, 'F') == 1) then
            call this%generate_kpath('fcc', npts_per_segment)
         else if (index(international, 'I') == 1) then
            call this%generate_kpath('bcc', npts_per_segment)
         else
            call this%generate_kpath('cubic', npts_per_segment)
         end if
      case ('hexagonal', 'trigonal')
         call this%generate_kpath('hexagonal', npts_per_segment)
      case ('tetragonal')
         call this%generate_kpath('tetragonal', npts_per_segment)
      case ('orthorhombic')
         call this%generate_kpath('orthorhombic', npts_per_segment)
      case ('monoclinic')
         call this%generate_kpath('monoclinic', npts_per_segment)
      case ('triclinic')
         call this%generate_kpath('triclinic', npts_per_segment)
      case default
         call g_logger%warning('get_canonical_kpath_for_spacegroup: Unknown crystal family, using cubic', __FILE__, __LINE__)
         call this%generate_kpath('cubic', npts_per_segment)
      end select

   end subroutine get_canonical_kpath_for_spacegroup

end module symmetry_mod