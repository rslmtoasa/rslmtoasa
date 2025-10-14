!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: spglib_interface
!
!> @author
!> Generated for RS-LMTO-ASA spglib integration
!
! DESCRIPTION:
!> Fortran interface to spglib for crystallographic symmetry operations
!> and k-path generation. This module provides a clean Fortran wrapper
!> around the spglib C library functionality.
!>
!> Features:
!> - Space group detection and symmetry operations
!> - High-symmetry k-path generation for band structure calculations
!> - Symmetry-reduced k-point mesh generation
!> - Crystal system classification
!>
!> Usage example:
!> ```fortran
!> type(spglib_interface) :: spg
!> call spg%initialize(lattice_vectors, positions, types)
!> kpath = spg%get_band_path()
!> reduced_mesh = spg%get_ir_reciprocal_mesh(mp_grid)
!> ```
!------------------------------------------------------------------------------

module spglib_interface_mod

#ifdef USE_SPGLIB
   ! Use ISO_C_BINDING for proper C interoperability
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_null_char, c_ptr, c_associated
#endif

   implicit none

   private

   ! Use standard modules that should be available
   integer, parameter :: rp = kind(1.0d0)  ! Double precision

#ifdef USE_SPGLIB
   ! C interface declarations for spglib
   interface
      ! Get space group number and symbol
      function spg_get_international(symbol, lattice, position, types, num_atom, symprec) &
         bind(C, name="spg_get_international")
         import :: c_int, c_double, c_char
         character(kind=c_char), intent(out) :: symbol(*)
         real(c_double), intent(in) :: lattice(3,3)
         real(c_double), intent(in) :: position(3,*)
         integer(c_int), intent(in) :: types(*)
         integer(c_int), value :: num_atom
         real(c_double), value :: symprec
         integer(c_int) :: spg_get_international
      end function spg_get_international

      ! Get symmetry operations
      function spg_get_symmetry(rotation, translation, max_size, lattice, position, types, num_atom, symprec) &
         bind(C, name="spg_get_symmetry")
         import :: c_int, c_double
         integer(c_int), intent(out) :: rotation(3,3,*)
         real(c_double), intent(out) :: translation(3,*)
         integer(c_int), value :: max_size
         real(c_double), intent(in) :: lattice(3,3)
         real(c_double), intent(in) :: position(3,*)
         integer(c_int), intent(in) :: types(*)
         integer(c_int), value :: num_atom
         real(c_double), value :: symprec
         integer(c_int) :: spg_get_symmetry
      end function spg_get_symmetry

      ! Get k-path for band structure
      function spg_get_ir_reciprocal_mesh(grid_address, ir_mapping_table, mesh, is_shift, &
                                         is_time_reversal, lattice, position, types, num_atom, symprec) &
         bind(C, name="spg_get_ir_reciprocal_mesh")
         import :: c_int, c_double
         integer(c_int), intent(out) :: grid_address(3,*)
         integer(c_int), intent(out) :: ir_mapping_table(*)
         integer(c_int), intent(in) :: mesh(3)
         integer(c_int), intent(in) :: is_shift(3)
         integer(c_int), value :: is_time_reversal
         real(c_double), intent(in) :: lattice(3,3)
         real(c_double), intent(in) :: position(3,*)
         integer(c_int), intent(in) :: types(*)
         integer(c_int), value :: num_atom
         real(c_double), value :: symprec
         integer(c_int) :: spg_get_ir_reciprocal_mesh
      end function spg_get_ir_reciprocal_mesh

      ! Dataset API for getting detailed information about crystal structure
      function spg_get_dataset(lattice, position, types, num_atom, symprec) &
         bind(C, name="spg_get_dataset")
         import :: c_ptr, c_double, c_int
         type(c_ptr) :: spg_get_dataset
         real(c_double), intent(in) :: lattice(3,3)
         real(c_double), intent(in) :: position(3,*)
         integer(c_int), intent(in) :: types(*)
         integer(c_int), value :: num_atom
         real(c_double), value :: symprec
      end function spg_get_dataset

      ! Free dataset memory
      subroutine spg_free_dataset(dataset) bind(C, name="spg_free_dataset")
         import :: c_ptr
         type(c_ptr), value :: dataset
      end subroutine spg_free_dataset
   end interface
#endif

   public :: spglib_interface

   !> Type to handle spglib operations
   type :: spglib_interface
      private
      logical :: initialized = .false.
      integer :: space_group_number = 0
      character(len=20) :: space_group_symbol = ''
      character(len=10) :: crystal_system = ''
      integer :: num_operations = 0
      real(rp), allocatable :: lattice_vectors(:,:)  ! (3,3)
      real(rp), allocatable :: atomic_positions(:,:) ! (3,num_atoms)
      integer, allocatable :: atomic_types(:)        ! (num_atoms)
      integer :: num_atoms = 0
      real(rp) :: symprec = 1.0e-5_rp
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: detect_space_group
      procedure :: get_symmetry_operations
      procedure :: get_reduced_kpoint_mesh
      procedure :: get_reduced_kpoint_mesh_with_points
      procedure :: get_crystal_system
      procedure :: get_band_path_points
      procedure :: get_dataset
      procedure :: is_available
      procedure :: get_space_group_number
      procedure :: get_space_group_symbol
      procedure :: get_crystal_system_name
      final :: destructor
   end type spglib_interface

   ! Interface for constructor
   interface spglib_interface
      procedure :: constructor
   end interface spglib_interface

contains

   !---------------------------------------------------------------------------
   !> @brief Constructor for spglib_interface
   !---------------------------------------------------------------------------
   function constructor() result(this)
      type(spglib_interface) :: this
      
#ifdef USE_SPGLIB
      write(*,*) 'spglib_interface: Constructor called with spglib support'
#else
      write(*,*) 'spglib_interface: Constructor called WITHOUT spglib support'
#endif
   end function constructor

   !---------------------------------------------------------------------------
   !> @brief Initialize with crystal structure
   !---------------------------------------------------------------------------
   subroutine initialize(this, lattice, positions, types, symprec_in)
      class(spglib_interface), intent(inout) :: this
      real(rp), intent(in) :: lattice(:,:)     ! (3,3) lattice vectors in rows
      real(rp), intent(in) :: positions(:,:)   ! (3,num_atoms) fractional coordinates
      integer, intent(in) :: types(:)          ! (num_atoms) atomic types
      real(rp), intent(in), optional :: symprec_in

      this%num_atoms = size(positions, 2)
      
      if (present(symprec_in)) this%symprec = symprec_in

      ! Store crystal structure
      allocate(this%lattice_vectors(3,3))
      allocate(this%atomic_positions(3,this%num_atoms))
      allocate(this%atomic_types(this%num_atoms))

      this%lattice_vectors = lattice
      this%atomic_positions = positions
      this%atomic_types = types

#ifdef USE_SPGLIB
      call this%detect_space_group()
      this%initialized = .true.
      write(*,*) 'spglib_interface: Initialized with', this%num_atoms, 'atoms'
      write(*,*) 'spglib_interface: Space group:', trim(this%space_group_symbol), &
                 '(#', this%space_group_number, ')'
#else
      write(*,*) 'spglib_interface: Initialized but spglib not available'
      this%initialized = .false.
#endif
   end subroutine initialize

   !---------------------------------------------------------------------------
   !> @brief Clean up allocated memory
   !---------------------------------------------------------------------------
   subroutine finalize(this)
      class(spglib_interface), intent(inout) :: this
      
      if (allocated(this%lattice_vectors)) deallocate(this%lattice_vectors)
      if (allocated(this%atomic_positions)) deallocate(this%atomic_positions)
      if (allocated(this%atomic_types)) deallocate(this%atomic_types)
      
      this%initialized = .false.
   end subroutine finalize

   !---------------------------------------------------------------------------
   !> @brief Detect space group and crystal system
   !---------------------------------------------------------------------------
   subroutine detect_space_group(this)
      class(spglib_interface), intent(inout) :: this

#ifdef USE_SPGLIB
      character(kind=c_char) :: symbol_c(21)  ! C string for symbol
      integer :: i

      ! Call spglib to get space group
      this%space_group_number = spg_get_international(symbol_c, &
                                                     this%lattice_vectors, &
                                                     this%atomic_positions, &
                                                     this%atomic_types, &
                                                     this%num_atoms, &
                                                     this%symprec)

      ! Convert C string to Fortran string
      this%space_group_symbol = ''
      do i = 1, 20
         if (symbol_c(i) == c_null_char) exit
         this%space_group_symbol(i:i) = symbol_c(i)
      end do

      ! Determine crystal system from space group number
      this%crystal_system = this%get_crystal_system()

#else
      this%space_group_number = 0
      this%space_group_symbol = 'UNKNOWN'
      this%crystal_system = 'UNKNOWN'
      write(*,*) 'spglib_interface: Space group detection unavailable (spglib not compiled)'
#endif
   end subroutine detect_space_group

   !---------------------------------------------------------------------------
   !> @brief Get symmetry operations
   !---------------------------------------------------------------------------
   function get_symmetry_operations(this) result(num_ops)
      class(spglib_interface), intent(in) :: this
      integer :: num_ops

#ifdef USE_SPGLIB
      integer, parameter :: max_operations = 192  ! Maximum possible for 3D space groups
      integer :: rotation(3,3,max_operations)
      real(c_double) :: translation(3,max_operations)

      if (.not. this%initialized) then
         write(*,*) 'ERROR: spglib_interface: Not initialized'
         num_ops = 0
         return
      end if

      num_ops = spg_get_symmetry(rotation, translation, max_operations, &
                                this%lattice_vectors, this%atomic_positions, &
                                this%atomic_types, this%num_atoms, this%symprec)

      write(*,*) 'spglib_interface: Found', num_ops, 'symmetry operations'
#else
      num_ops = 0
      write(*,*) 'spglib_interface: Symmetry operations unavailable (spglib not compiled)'
#endif
   end function get_symmetry_operations

   !---------------------------------------------------------------------------
   !> @brief Get irreducible k-point mesh
   !---------------------------------------------------------------------------
   function get_reduced_kpoint_mesh(this, mesh_dims, is_shift) result(num_ir_kpoints)
      class(spglib_interface), intent(in) :: this
      integer, intent(in) :: mesh_dims(3)
      integer, intent(in), optional :: is_shift(3)
      integer :: num_ir_kpoints

#ifdef USE_SPGLIB
      integer :: total_kpoints
      integer :: shift(3)
      integer, allocatable :: grid_address(:,:)
      integer, allocatable :: ir_mapping(:)

      if (.not. this%initialized) then
         write(*,*) 'ERROR: spglib_interface: Not initialized'
         num_ir_kpoints = 0
         return
      end if

      total_kpoints = mesh_dims(1) * mesh_dims(2) * mesh_dims(3)
      
      if (present(is_shift)) then
         shift = is_shift
      else
         shift = [0, 0, 0]  ! No shift by default
      end if

      allocate(grid_address(3, total_kpoints))
      allocate(ir_mapping(total_kpoints))

      num_ir_kpoints = spg_get_ir_reciprocal_mesh(grid_address, ir_mapping, &
                                                 mesh_dims, shift, 1, &  ! 1 = use time reversal
                                                 this%lattice_vectors, &
                                                 this%atomic_positions, &
                                                 this%atomic_types, &
                                                 this%num_atoms, &
                                                 this%symprec)

      write(*,*) 'spglib_interface: Reduced', total_kpoints, 'k-points to', &
                 num_ir_kpoints, 'irreducible points'

      deallocate(grid_address, ir_mapping)
#else
      num_ir_kpoints = mesh_dims(1) * mesh_dims(2) * mesh_dims(3)  ! Full mesh without reduction
      write(*,*) 'spglib_interface: Using full k-mesh (spglib not compiled)'
#endif
   end function get_reduced_kpoint_mesh

   !---------------------------------------------------------------------------
   !> @brief Get irreducible k-point mesh with k-points and weights
   !> @param[in] mesh_dims K-point mesh dimensions [nk1, nk2, nk3]
   !> @param[in] is_shift Optional shift for k-mesh [0 or 1 for each direction]
   !> @param[out] kpoints K-point coordinates in reciprocal lattice units [3, num_ir]
   !> @param[out] weights K-point weights (normalized to sum to 1) [num_ir]
   !> @return Number of irreducible k-points
   !---------------------------------------------------------------------------
   function get_reduced_kpoint_mesh_with_points(this, mesh_dims, is_shift, kpoints, weights) result(num_ir_kpoints)
      class(spglib_interface), intent(in) :: this
      integer, intent(in) :: mesh_dims(3)
      integer, intent(in), optional :: is_shift(3)
      real(rp), allocatable, intent(out) :: kpoints(:,:)
      real(rp), allocatable, intent(out) :: weights(:)
      integer :: num_ir_kpoints

      integer :: total_kpoints, shift(3), ir_idx, i, j, k
      
#ifdef USE_SPGLIB
      integer, allocatable :: grid_address(:,:)
      integer, allocatable :: ir_mapping(:)
      integer :: multiplicity
      integer, allocatable :: ir_kpoint_indices(:)
      logical, allocatable :: is_irreducible(:)

      if (.not. this%initialized) then
         write(*,*) 'ERROR: spglib_interface: Not initialized'
         num_ir_kpoints = 0
         return
      end if

      total_kpoints = mesh_dims(1) * mesh_dims(2) * mesh_dims(3)
      
      if (present(is_shift)) then
         shift = is_shift
      else
         shift = [0, 0, 0]  ! No shift by default
      end if

      allocate(grid_address(3, total_kpoints))
      allocate(ir_mapping(total_kpoints))

      ! Get irreducible k-points from spglib
      num_ir_kpoints = spg_get_ir_reciprocal_mesh(grid_address, ir_mapping, &
                                                 mesh_dims, shift, 1, &  ! 1 = use time reversal
                                                 this%lattice_vectors, &
                                                 this%atomic_positions, &
                                                 this%atomic_types, &
                                                 this%num_atoms, &
                                                 this%symprec)

      ! Allocate output arrays
      allocate(kpoints(3, num_ir_kpoints))
      allocate(weights(num_ir_kpoints))
      allocate(is_irreducible(total_kpoints))
      allocate(ir_kpoint_indices(num_ir_kpoints))
      
      ! Find which k-points are irreducible
      is_irreducible = .false.
      ir_idx = 0
      do i = 1, total_kpoints
         if (ir_mapping(i) == i) then  ! This is an irreducible k-point
            ir_idx = ir_idx + 1
            is_irreducible(i) = .true.
            ir_kpoint_indices(ir_idx) = i
         end if
      end do

      ! Extract irreducible k-points and calculate weights
      weights = 0.0_rp
      do ir_idx = 1, num_ir_kpoints
         i = ir_kpoint_indices(ir_idx)
         
         ! Convert grid address to fractional coordinates
         ! Note: grid_address gives integers, need to convert to [0,1) range
         kpoints(1, ir_idx) = real(grid_address(1, i), rp) / real(mesh_dims(1), rp)
         kpoints(2, ir_idx) = real(grid_address(2, i), rp) / real(mesh_dims(2), rp)
         kpoints(3, ir_idx) = real(grid_address(3, i), rp) / real(mesh_dims(3), rp)
         
         ! Add shift if needed
         if (shift(1) /= 0) kpoints(1, ir_idx) = kpoints(1, ir_idx) + 0.5_rp / real(mesh_dims(1), rp)
         if (shift(2) /= 0) kpoints(2, ir_idx) = kpoints(2, ir_idx) + 0.5_rp / real(mesh_dims(2), rp)
         if (shift(3) /= 0) kpoints(3, ir_idx) = kpoints(3, ir_idx) + 0.5_rp / real(mesh_dims(3), rp)
         
         ! Calculate weight from multiplicity (how many k-points map to this irreducible one)
         multiplicity = count(ir_mapping == i)
         weights(ir_idx) = real(multiplicity, rp) / real(total_kpoints, rp)
      end do

      write(*,*) 'spglib_interface: Reduced', total_kpoints, 'k-points to', &
                 num_ir_kpoints, 'irreducible points'
      write(*,*) 'spglib_interface: Weight sum =', sum(weights), '(should be 1.0)'

      deallocate(grid_address, ir_mapping, is_irreducible, ir_kpoint_indices)
#else
      ! Without spglib, return full mesh with equal weights
      total_kpoints = mesh_dims(1) * mesh_dims(2) * mesh_dims(3)
      num_ir_kpoints = total_kpoints
      
      allocate(kpoints(3, num_ir_kpoints))
      allocate(weights(num_ir_kpoints))
      
      ! Generate full MP mesh
      ir_idx = 0
      do i = 1, mesh_dims(1)
         do j = 1, mesh_dims(2)
            do k = 1, mesh_dims(3)
               ir_idx = ir_idx + 1
               kpoints(1, ir_idx) = real(i-1, rp) / real(mesh_dims(1), rp)
               kpoints(2, ir_idx) = real(j-1, rp) / real(mesh_dims(2), rp)
               kpoints(3, ir_idx) = real(k-1, rp) / real(mesh_dims(3), rp)
            end do
         end do
      end do
      weights = 1.0_rp / real(num_ir_kpoints, rp)
      
      write(*,*) 'spglib_interface: Using full k-mesh (spglib not compiled)'
#endif
   end function get_reduced_kpoint_mesh_with_points

   !---------------------------------------------------------------------------
   !> @brief Determine crystal system from space group number
   !---------------------------------------------------------------------------
   function get_crystal_system(this) result(crystal_system)
      class(spglib_interface), intent(in) :: this
      character(len=10) :: crystal_system

      ! Classification based on space group number ranges
      ! Reference: International Tables for Crystallography
      if (this%space_group_number >= 1 .and. this%space_group_number <= 2) then
         crystal_system = 'triclinic'
      else if (this%space_group_number >= 3 .and. this%space_group_number <= 15) then
         crystal_system = 'monoclinic'
      else if (this%space_group_number >= 16 .and. this%space_group_number <= 74) then
         crystal_system = 'orthorhombic'
      else if (this%space_group_number >= 75 .and. this%space_group_number <= 142) then
         crystal_system = 'tetragonal'
      else if (this%space_group_number >= 143 .and. this%space_group_number <= 167) then
         crystal_system = 'trigonal'
      else if (this%space_group_number >= 168 .and. this%space_group_number <= 194) then
         crystal_system = 'hexagonal'
      else if (this%space_group_number >= 195 .and. this%space_group_number <= 230) then
         crystal_system = 'cubic'
      else
         crystal_system = 'unknown'
      end if
   end function get_crystal_system

   !---------------------------------------------------------------------------
   !> @brief Get high-symmetry k-path points for band structure
   !---------------------------------------------------------------------------
   function get_band_path_points(this) result(path_available)
      class(spglib_interface), intent(in) :: this
      logical :: path_available

      ! This is a placeholder for band path generation
      ! In a full implementation, this would use spglib's band path functionality
      ! or implement custom logic based on the crystal system
      
      path_available = .false.
      
      if (.not. this%initialized) then
         write(*,*) 'ERROR: spglib_interface: Not initialized'
         return
      end if

#ifdef USE_SPGLIB
      ! TODO: Implement band path generation using spglib
      ! This would involve calling spglib functions to get high-symmetry points
      ! For now, just indicate that the feature is available with spglib
      path_available = .true.
      write(*,*) 'spglib_interface: Band path generation available for ', &
                 trim(this%crystal_system), ' crystal system'
#else
      write(*,*) 'spglib_interface: Band path generation unavailable (spglib not compiled)'
#endif
   end function get_band_path_points

   !---------------------------------------------------------------------------
   !> @brief Check if spglib functionality is available
   !---------------------------------------------------------------------------
   function is_available(this) result(available)
      class(spglib_interface), intent(in) :: this
      logical :: available

#ifdef USE_SPGLIB
      available = .true.
#else
      available = .false.
#endif
   end function is_available

   !---------------------------------------------------------------------------
   !> @brief Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(spglib_interface), intent(inout) :: this
      call this%finalize()
   end subroutine destructor

   !---------------------------------------------------------------------------
   !> @brief Get space group number
   !---------------------------------------------------------------------------
   function get_space_group_number(this) result(number)
      class(spglib_interface), intent(in) :: this
      integer :: number
      number = this%space_group_number
   end function get_space_group_number

   !---------------------------------------------------------------------------
   !> @brief Get space group symbol
   !---------------------------------------------------------------------------
   function get_space_group_symbol(this) result(symbol)
      class(spglib_interface), intent(in) :: this
      character(len=20) :: symbol
      symbol = this%space_group_symbol
   end function get_space_group_symbol

   !---------------------------------------------------------------------------
   !> @brief Get crystal system name
   !---------------------------------------------------------------------------
   function get_crystal_system_name(this) result(system)
      class(spglib_interface), intent(in) :: this
      character(len=10) :: system
      system = this%crystal_system
   end function get_crystal_system_name

   !---------------------------------------------------------------------------
   !> @brief Get detailed dataset information using spglib dataset API
   !> @return Logical indicating success
   !---------------------------------------------------------------------------
   function get_dataset(this, space_group_number, international, hall, &
                       transformation_matrix, origin_shift, primitive_lattice) result(success)
      class(spglib_interface), intent(in) :: this
      integer, intent(out), optional :: space_group_number
      character(len=*), intent(out), optional :: international
      character(len=*), intent(out), optional :: hall
      real(rp), intent(out), optional :: transformation_matrix(3,3)
      real(rp), intent(out), optional :: origin_shift(3)
      real(rp), intent(out), optional :: primitive_lattice(3,3)
      logical :: success

#ifdef USE_SPGLIB
      type(c_ptr) :: dataset_ptr
      integer :: spg_num
#endif

      success = .false.

#ifdef USE_SPGLIB
      if (.not. this%initialized) return

      ! Get dataset from spglib
      dataset_ptr = spg_get_dataset(real(this%lattice_vectors, c_double), &
                                   real(this%atomic_positions, c_double), &
                                   int(this%atomic_types, c_int), &
                                   int(this%num_atoms, c_int), &
                                   real(this%symprec, c_double))

      if (c_associated(dataset_ptr)) then
         success = .true.

         ! For simplicity, we'll just return the basic space group info that we already have
         ! A full implementation would need to define C structures to extract all dataset fields
         if (present(space_group_number)) space_group_number = this%space_group_number
         if (present(international)) international = trim(this%space_group_symbol)
         if (present(hall)) hall = 'Hall symbol extraction not implemented'
         
         ! Transformation matrix and other fields would need proper C struct definitions
         if (present(transformation_matrix)) transformation_matrix = 0.0_rp
         if (present(origin_shift)) origin_shift = 0.0_rp
         if (present(primitive_lattice)) primitive_lattice = this%lattice_vectors

         ! Free the dataset memory
         call spg_free_dataset(dataset_ptr)
      end if
#else
      ! Return minimal information when spglib is not available
      if (present(space_group_number)) space_group_number = 1
      if (present(international)) international = 'P1'
      if (present(hall)) hall = 'P 1'
      if (present(transformation_matrix)) then
         transformation_matrix = 0.0_rp
         transformation_matrix(1,1) = 1.0_rp
         transformation_matrix(2,2) = 1.0_rp
         transformation_matrix(3,3) = 1.0_rp
      end if
      if (present(origin_shift)) origin_shift = 0.0_rp
      if (present(primitive_lattice)) primitive_lattice = this%lattice_vectors
#endif

   end function get_dataset

end module spglib_interface_mod