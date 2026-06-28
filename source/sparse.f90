!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Sparse
!
!> @brief
!> Sparse matrix algebra helper for TB-LMTO Hamiltonians.
!> Provides Block Sparse Row (BSR), Compressed Sparse Row (CSR), and
!> Coordinate (COO) formats compatible with oneMKL, cuSPARSE, SciPy, PyTorch.
!
!> Usage:
!>   hamiltonian_obj = hamiltonian(charge_obj)
!>   sparse_obj = sparse(hamiltonian_obj)
!>   recursion_obj = recursion(hamiltonian_obj, energy_obj, sparse_obj)
!
!------------------------------------------------------------------------------

module sparse_mod

   use precision_mod, only: rp
   use basis_mod, only: nb
   use iso_c_binding
   use hamiltonian_mod, only: hamiltonian
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Public interface for sparse type construction
   public :: sparse
   !> Public types for direct access when needed
   public :: bsr_matrix, csr_matrix, coo_matrix, bsr_matrix_c_handle

   !> Public utility routines
   public :: block_to_sparse_from_arrays, export_bsr_to_c_format
   public :: validate_bsr

   !---------------------------------------------------------------------------
   ! TYPE DEFINITIONS (must come before usage)
   !---------------------------------------------------------------------------

   !> BSR matrix (Block Sparse Row format)
   type :: bsr_matrix
      integer :: nblocks, nrows, ncols, blocksize
      complex(rp), allocatable :: values(:,:,:)
      integer, allocatable :: col_indices(:), row_ptr(:)
   contains
      procedure :: print_info_bsr => bsr_print_info
   end type bsr_matrix

   !> CSR matrix (Compressed Sparse Row format)
   type :: csr_matrix
      integer :: m, n, nnz
      complex(rp), allocatable :: values(:)
      integer, allocatable :: col_indices(:), row_ptr(:)
   end type csr_matrix

   !> COO matrix (Coordinate format)
   type :: coo_matrix
      integer :: m, n, nnz
      integer, allocatable :: row_indices(:), col_indices(:)
      complex(rp), allocatable :: values(:)
   end type coo_matrix

   !> C-compatible BSR handle for GPU/Python libraries
   type :: bsr_matrix_c_handle
      integer(c_int), allocatable :: row_ptr_c(:), col_idx_c(:)
      complex(c_double_complex), allocatable :: vals_c(:,:,:), vals_c_flat(:)
      integer(c_int) :: nrows, ncols, nblocks, blocksize
      character(len=1) :: block_order = 'F'
      type(c_ptr) :: onemkl_handle = c_null_ptr
      type(c_ptr) :: cusparse_handle = c_null_ptr
      logical :: is_initialized = .false.
   end type bsr_matrix_c_handle

   !---------------------------------------------------------------------------
   ! MAIN SPARSE TYPE
   !---------------------------------------------------------------------------

   type :: sparse
      type(hamiltonian), pointer :: hamiltonian => null()
      type(bsr_matrix) :: h_bsr
      type(bsr_matrix_c_handle) :: h_c_gpu
      type(bsr_matrix_c_handle) :: h_c_python
   contains
      procedure :: build_bsr
      procedure :: build_bsr_from_operator
      procedure :: export_gpu_onemkl
      procedure :: export_gpu_cusparse
      procedure :: export_python_scipy
      procedure :: export_python_torch
      procedure :: print_info
      procedure :: restore_to_default
   end type sparse

   !> Interface for sparse type construction
   interface sparse
      procedure :: constructor
   end interface sparse

contains

   !---------------------------------------------------------------------------
   ! CONSTRUCTOR
   !---------------------------------------------------------------------------
   function constructor(hamiltonian_obj) result(obj)
      type(sparse) :: obj
      type(hamiltonian), target, intent(in) :: hamiltonian_obj

      obj%hamiltonian => hamiltonian_obj
      call obj%restore_to_default()
   end function constructor

   !---------------------------------------------------------------------------
   ! DESTRUCTOR
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(sparse) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%h_bsr%values)) &
         call g_safe_alloc%deallocate('sparse.h_bsr.values', this%h_bsr%values)
      if (allocated(this%h_bsr%col_indices)) &
         call g_safe_alloc%deallocate('sparse.h_bsr.col_indices', this%h_bsr%col_indices)
      if (allocated(this%h_bsr%row_ptr)) &
         call g_safe_alloc%deallocate('sparse.h_bsr.row_ptr', this%h_bsr%row_ptr)
#else
      if (allocated(this%h_bsr%values)) deallocate(this%h_bsr%values)
      if (allocated(this%h_bsr%col_indices)) deallocate(this%h_bsr%col_indices)
      if (allocated(this%h_bsr%row_ptr)) deallocate(this%h_bsr%row_ptr)
#endif
   end subroutine destructor

   !---------------------------------------------------------------------------
   ! RESET TO DEFAULT
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      class(sparse), intent(inout) :: this

      if (allocated(this%h_bsr%values)) deallocate(this%h_bsr%values)
      if (allocated(this%h_bsr%col_indices)) deallocate(this%h_bsr%col_indices)
      if (allocated(this%h_bsr%row_ptr)) deallocate(this%h_bsr%row_ptr)

      this%h_bsr%nblocks = 0
      this%h_bsr%nrows = 0
      this%h_bsr%ncols = 0
      this%h_bsr%blocksize = nb

      this%h_c_gpu%is_initialized = .false.
      this%h_c_python%is_initialized = .false.
   end subroutine restore_to_default

   !---------------------------------------------------------------------------
   ! BUILD BSR
   !---------------------------------------------------------------------------
   subroutine build_bsr(this)
      class(sparse), intent(inout) :: this

      if (.not. associated(this%hamiltonian)) then
         print *, 'ERROR: hamiltonian not associated in sparse object'
         return
      end if

      call block_to_sparse_from_arrays( &
          this%h_bsr, this%hamiltonian%hall, this%hamiltonian%ee, &
          this%hamiltonian%lattice%nn, this%hamiltonian%lattice%iz, &
          this%hamiltonian%lattice%kk, this%hamiltonian%lattice%nmax, nb)
   end subroutine build_bsr

   subroutine build_bsr_from_operator(this, hall, ee)
      class(sparse), intent(inout) :: this
      complex(rp), intent(in) :: hall(:,:,:,:), ee(:,:,:,:)

      if (.not. associated(this%hamiltonian)) then
         print *, 'ERROR: hamiltonian not associated in sparse object'
         return
      end if

      call block_to_sparse_from_arrays( &
          this%h_bsr, hall, ee, &
          this%hamiltonian%lattice%nn, this%hamiltonian%lattice%iz, &
          this%hamiltonian%lattice%kk, this%hamiltonian%lattice%nmax, nb)

      this%h_c_gpu%is_initialized = .false.
      this%h_c_python%is_initialized = .false.
   end subroutine build_bsr_from_operator

   !---------------------------------------------------------------------------
   ! EXPORT ROUTINES
   !---------------------------------------------------------------------------
   subroutine export_gpu_onemkl(this, is_valid)
      class(sparse), intent(inout) :: this
      logical, intent(out) :: is_valid
      if (this%h_bsr%nblocks == 0) call this%build_bsr()
      call export_bsr_to_c_format(this%h_bsr, this%h_c_gpu, 'F', is_valid)
   end subroutine export_gpu_onemkl

   subroutine export_gpu_cusparse(this, is_valid)
      class(sparse), intent(inout) :: this
      logical, intent(out) :: is_valid
      if (this%h_bsr%nblocks == 0) call this%build_bsr()
      call export_bsr_to_c_format(this%h_bsr, this%h_c_gpu, 'F', is_valid)
   end subroutine export_gpu_cusparse

   subroutine export_python_scipy(this, is_valid)
      class(sparse), intent(inout) :: this
      logical, intent(out) :: is_valid
      if (this%h_bsr%nblocks == 0) call this%build_bsr()
      call export_bsr_to_c_format(this%h_bsr, this%h_c_python, 'C', is_valid)
   end subroutine export_python_scipy

   subroutine export_python_torch(this, is_valid)
      class(sparse), intent(inout) :: this
      logical, intent(out) :: is_valid
      if (this%h_bsr%nblocks == 0) call this%build_bsr()
      call export_bsr_to_c_format(this%h_bsr, this%h_c_python, 'C', is_valid)
   end subroutine export_python_torch

   !---------------------------------------------------------------------------
   ! PRINT INFO
   !---------------------------------------------------------------------------
   subroutine print_info(this)
      class(sparse), intent(in) :: this
      print '(A)', '=================================================='
      print '(A)', 'Sparse Matrix (BSR) Structure'
      print '(A)', '=================================================='
      if (this%h_bsr%nblocks == 0) then
         print '(A)', '(not built yet - call build_bsr())'
      else
         print '(A,I0,A,I0,A,I0)', 'BSR: ', this%h_bsr%nblocks, ' blocks, ', &
            this%h_bsr%nrows, ' block-rows, blocksize=', this%h_bsr%blocksize
         print '(A,I0)', 'Total elements: ', this%h_bsr%nblocks * this%h_bsr%blocksize * this%h_bsr%blocksize
      end if
      print '(A)', '=================================================='
   end subroutine print_info

   !---------------------------------------------------------------------------
   ! UTILITY SUBROUTINES
   !---------------------------------------------------------------------------

   subroutine allocate_bsr(bsr, nblocks, nrows, blocksize)
      type(bsr_matrix), intent(inout) :: bsr
      integer, intent(in) :: nblocks, nrows, blocksize
      bsr%nblocks = nblocks
      bsr%nrows = nrows
      bsr%ncols = nrows
      bsr%blocksize = blocksize
      allocate(bsr%values(blocksize, blocksize, nblocks))
      allocate(bsr%col_indices(nblocks))
      allocate(bsr%row_ptr(nrows + 1))
      bsr%values = (0.0_rp, 0.0_rp)
      bsr%col_indices = 0
      bsr%row_ptr = 0
   end subroutine allocate_bsr

   subroutine bsr_print_info(this)
      class(bsr_matrix), intent(in) :: this
      print '(A,I0,A,I0,A,I0)', 'BSR: ', this%nblocks, ' blocks, ', &
         this%nrows, ' block-rows, blocksize=', this%blocksize
      print '(A,I0)', 'Total elements: ', this%nblocks * this%blocksize * this%blocksize
   end subroutine bsr_print_info

   function validate_bsr(bsr) result(is_valid)
      type(bsr_matrix), intent(in) :: bsr
      logical :: is_valid
      integer :: i, k
      is_valid = .true.
      do i = 1, bsr%nrows
         if (bsr%row_ptr(i) > bsr%row_ptr(i + 1)) then
            print *, 'ERROR: row_ptr not strictly increasing'
            is_valid = .false.
            return
         end if
      end do
      if (bsr%row_ptr(1) /= 1 .or. bsr%row_ptr(bsr%nrows + 1) /= bsr%nblocks + 1) then
         print *, 'ERROR: invalid row_ptr bounds'
         is_valid = .false.
         return
      end if
      do k = 1, bsr%nblocks
         if (bsr%col_indices(k) < 1 .or. bsr%col_indices(k) > bsr%ncols) then
            print *, 'ERROR: invalid column index in block', k
            is_valid = .false.
            return
         end if
      end do
   end function validate_bsr

   subroutine block_to_sparse_from_arrays(bsr, hall, ee, nn, iz, kk, nmax, blocksize)
      type(bsr_matrix), intent(inout) :: bsr
      complex(rp), intent(in) :: hall(:,:,:,:), ee(:,:,:,:)
      integer, intent(in) :: nn(:,:), iz(:), kk, nmax, blocksize
      integer :: atom, neighbor_idx, neighbor, block_idx, nblocks, num_neighbors, block_col, ih
      logical :: is_impurity

      nblocks = 0
      do atom = 1, kk
         num_neighbors = nn(atom, 1)
         do neighbor_idx = 1, num_neighbors
            if (neighbor_idx == 1 .or. nn(atom, neighbor_idx) /= 0) nblocks = nblocks + 1
         end do
      end do

      call allocate_bsr(bsr, nblocks, kk, blocksize)
      block_idx = 0
      bsr%row_ptr(1) = 1
      is_impurity = (nmax > 0)

      do atom = 1, kk
         num_neighbors = nn(atom, 1)
         do neighbor_idx = 1, num_neighbors
            if (neighbor_idx == 1) then
               block_col = atom
            else
               neighbor = nn(atom, neighbor_idx)
               if (neighbor == 0) cycle
               block_col = neighbor
            end if
            block_idx = block_idx + 1
            if (is_impurity .and. atom <= nmax) then
               bsr%values(:, :, block_idx) = hall(:, :, neighbor_idx, atom)
            else
               ih = iz(atom)
               bsr%values(:, :, block_idx) = ee(:, :, neighbor_idx, ih)
            end if
            bsr%col_indices(block_idx) = block_col
         end do
         if (atom < kk) bsr%row_ptr(atom + 1) = block_idx + 1
      end do

      bsr%row_ptr(kk + 1) = nblocks + 1
   end subroutine block_to_sparse_from_arrays

   subroutine export_bsr_to_c_format(bsr_f, bsr_c, block_order, is_valid)
      type(bsr_matrix), intent(in) :: bsr_f
      type(bsr_matrix_c_handle), intent(out) :: bsr_c
      character(len=1), intent(in) :: block_order
      logical, intent(out) :: is_valid
      integer :: k, bi, bj, idx

      is_valid = validate_bsr(bsr_f)
      if (.not. is_valid) return

      bsr_c%nrows = int(bsr_f%nrows, c_int)
      bsr_c%ncols = int(bsr_f%ncols, c_int)
      bsr_c%nblocks = int(bsr_f%nblocks, c_int)
      bsr_c%blocksize = int(bsr_f%blocksize, c_int)
      bsr_c%block_order = block_order

      allocate(bsr_c%row_ptr_c(bsr_f%nrows + 1))
      allocate(bsr_c%col_idx_c(bsr_f%nblocks))

      do k = 1, bsr_f%nrows + 1
         bsr_c%row_ptr_c(k) = int(bsr_f%row_ptr(k) - 1, c_int)
      end do
      do k = 1, bsr_f%nblocks
         bsr_c%col_idx_c(k) = int(bsr_f%col_indices(k) - 1, c_int)
      end do

      if (block_order == 'F') then
         allocate(bsr_c%vals_c(bsr_f%blocksize, bsr_f%blocksize, bsr_f%nblocks))
         do k = 1, bsr_f%nblocks
            bsr_c%vals_c(:, :, k) = cmplx(bsr_f%values(:, :, k), kind=c_double_complex)
         end do
      else if (block_order == 'C') then
         allocate(bsr_c%vals_c_flat(bsr_f%nblocks * bsr_f%blocksize * bsr_f%blocksize))
         idx = 1
         do k = 1, bsr_f%nblocks
            do bi = 1, bsr_f%blocksize
               do bj = 1, bsr_f%blocksize
                  bsr_c%vals_c_flat(idx) = cmplx(bsr_f%values(bi, bj, k), kind=c_double_complex)
                  idx = idx + 1
               end do
            end do
         end do
      else
         print *, 'ERROR: block_order must be "F" or "C"'
         is_valid = .false.
         return
      end if

      bsr_c%is_initialized = .true.
   end subroutine export_bsr_to_c_format

end module sparse_mod
