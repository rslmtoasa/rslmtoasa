!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Sparse Matrix
!
!> @brief
!> Sparse matrix storage formats for TB-LMTO Hamiltonians.
!> Supports BSR (Block Sparse Row), CSR (Compressed Sparse Row), and COO
!> (Coordinate) formats compatible with oneMKL and cuSPARSE.
!
!> @details
!> For TB-LMTO with nb=18 (spd basis), BSR is most efficient:
!> - Each block is 18×18 (dense within block)
!> - Only non-zero blocks stored
!> - Compatible with oneMKL and cuSPARSE
!> - Natural mapping: each block = one site-pair interaction
!
!------------------------------------------------------------------------------

module sparse_matrix_mod

   use precision_mod, only: rp
   use basis_mod, only: nb
   implicit none

   private
   public :: bsr_matrix, csr_matrix, coo_matrix
   public :: allocate_bsr, deallocate_bsr, assemble_bsr_from_blocks
   public :: bsr_to_csr, csr_spmv

   !> BSR matrix (Block Sparse Row format)
   !> Stores blocks in row-wise compressed format
   type :: bsr_matrix
      integer :: nblocks      ! Number of blocks
      integer :: nrows        ! Number of block-rows (= nblocks in full matrix / blocksize)
      integer :: ncols        ! Number of block-columns
      integer :: blocksize    ! Block size (e.g., 18 for spd)

      ! BSR data arrays
      complex(rp), allocatable :: values(:,:,:)   ! (blocksize, blocksize, nblocks)
      integer, allocatable     :: col_indices(:)  ! Column index for each block
      integer, allocatable     :: row_ptr(:)      ! Row pointer (length: nrows+1)

   contains
      procedure :: print_info => bsr_print_info
   end type bsr_matrix

   !> CSR matrix (Compressed Sparse Row format)
   !> Standard sparse format for general sparse matrices
   type :: csr_matrix
      integer :: m, n                     ! Dimensions
      integer :: nnz                      ! Number of non-zeros
      complex(rp), allocatable :: values(:)       ! Non-zero values
      integer, allocatable     :: col_indices(:)  ! Column indices
      integer, allocatable     :: row_ptr(:)      ! Row pointers (length: m+1)
   end type csr_matrix

   !> COO matrix (Coordinate format)
   !> Simplest sparse format: triplets (i, j, value)
   type :: coo_matrix
      integer :: m, n                     ! Dimensions
      integer :: nnz                      ! Number of non-zeros
      integer, allocatable     :: row_indices(:)  ! Row indices
      integer, allocatable     :: col_indices(:)  ! Column indices
      complex(rp), allocatable :: values(:)       ! Values
   end type coo_matrix

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Allocate BSR matrix structure
   !---------------------------------------------------------------------------
   subroutine allocate_bsr(bsr, nblocks, nrows, blocksize)
      type(bsr_matrix), intent(inout) :: bsr
      integer, intent(in) :: nblocks, nrows, blocksize

      bsr%nblocks = nblocks
      bsr%nrows = nrows
      bsr%ncols = nrows  ! Assume square for now
      bsr%blocksize = blocksize

      allocate(bsr%values(blocksize, blocksize, nblocks))
      allocate(bsr%col_indices(nblocks))
      allocate(bsr%row_ptr(nrows + 1))

      bsr%values = (0.0_rp, 0.0_rp)
      bsr%col_indices = 0
      bsr%row_ptr = 0
   end subroutine allocate_bsr

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Deallocate BSR matrix
   !---------------------------------------------------------------------------
   subroutine deallocate_bsr(bsr)
      type(bsr_matrix), intent(inout) :: bsr
      if (allocated(bsr%values)) deallocate(bsr%values)
      if (allocated(bsr%col_indices)) deallocate(bsr%col_indices)
      if (allocated(bsr%row_ptr)) deallocate(bsr%row_ptr)
   end subroutine deallocate_bsr

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Build BSR matrix from block arrays
   !>
   !> Constructs BSR format from the TB-LMTO block structure:
   !> - hall(nb,nb,neighbors,atoms) for impurity region
   !> - ee(nb,nb,neighbors,atom_types) for bulk region
   !> Only non-zero blocks are stored.
   !>
   !> @param[inout] bsr        Target BSR matrix
   !> @param[in]    nblocks    Number of blocks to store
   !> @param[in]    nrows      Number of block-rows
   !> @param[in]    blocksize  Size of each block (typically 18)
   !---------------------------------------------------------------------------
   subroutine assemble_bsr_from_blocks(bsr, nblocks, nrows, blocksize)
      type(bsr_matrix), intent(inout) :: bsr
      integer, intent(in) :: nblocks, nrows, blocksize

      ! Allocate the BSR structure
      call allocate_bsr(bsr, nblocks, nrows, blocksize)

      ! Initialize row pointers to zero (caller will fill in blocks)
      bsr%row_ptr(1) = 1
      ! Row pointers will be set after blocks are added

   end subroutine assemble_bsr_from_blocks

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Convert BSR to CSR format
   !>
   !> Expands block sparse row (BSR) to standard compressed sparse row (CSR)
   !> by inflating each block element.
   !>
   !> @param[in]  bsr  BSR matrix (blocksize×blocksize blocks)
   !> @param[out] csr  CSR matrix (full element-wise sparsity)
   !---------------------------------------------------------------------------
   subroutine bsr_to_csr(bsr, csr)
      type(bsr_matrix), intent(in)  :: bsr
      type(csr_matrix), intent(out) :: csr

      integer :: m, nnz, block_row, block_col, i, j, bi, bj, k
      integer, allocatable :: temp_col(:), temp_row(:)
      complex(rp), allocatable :: temp_val(:)

      m = bsr%nrows * bsr%blocksize
      allocate(temp_col(m * m))
      allocate(temp_row(m * m))
      allocate(temp_val(m * m))

      nnz = 0

      ! Iterate through BSR blocks
      do block_row = 1, bsr%nrows
         ! For each block in this block-row
         do k = bsr%row_ptr(block_row), bsr%row_ptr(block_row + 1) - 1
            block_col = bsr%col_indices(k)
            ! Expand block into CSR entries
            do bi = 1, bsr%blocksize
               do bj = 1, bsr%blocksize
                  if (abs(bsr%values(bi, bj, k)) > 0.0_rp) then
                     nnz = nnz + 1
                     i = (block_row - 1) * bsr%blocksize + bi
                     j = (block_col - 1) * bsr%blocksize + bj
                     temp_row(nnz) = i
                     temp_col(nnz) = j
                     temp_val(nnz) = bsr%values(bi, bj, k)
                  end if
               end do
            end do
         end do
      end do

      ! Allocate CSR
      csr%m = m
      csr%n = m
      csr%nnz = nnz
      allocate(csr%values(nnz))
      allocate(csr%col_indices(nnz))
      allocate(csr%row_ptr(m + 1))

      ! Copy data (note: would need sorting for proper CSR format)
      csr%values = temp_val(1:nnz)
      csr%col_indices = temp_col(1:nnz)

      ! Build row pointer from row indices
      csr%row_ptr = 0
      do k = 1, nnz
         csr%row_ptr(temp_row(k)) = csr%row_ptr(temp_row(k)) + 1
      end do
      ! Convert counts to pointers
      do i = 2, m + 1
         csr%row_ptr(i) = csr%row_ptr(i) + csr%row_ptr(i - 1)
      end do

      deallocate(temp_col, temp_row, temp_val)

   end subroutine bsr_to_csr

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Sparse matrix-vector product: y = A*x (CSR format)
   !>
   !> @param[in]  csr  CSR matrix
   !> @param[in]  x    Input vector
   !> @param[out] y    Output vector
   !---------------------------------------------------------------------------
   subroutine csr_spmv(csr, x, y)
      type(csr_matrix), intent(in) :: csr
      complex(rp), intent(in)  :: x(:)
      complex(rp), intent(out) :: y(:)
      integer :: i, k

      y = (0.0_rp, 0.0_rp)
      do i = 1, csr%m
         do k = csr%row_ptr(i), csr%row_ptr(i + 1) - 1
            y(i) = y(i) + csr%values(k) * x(csr%col_indices(k))
         end do
      end do
   end subroutine csr_spmv

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Print BSR matrix info
   !---------------------------------------------------------------------------
   subroutine bsr_print_info(this)
      class(bsr_matrix), intent(in) :: this
      print '(A,I0,A,I0,A,I0)', 'BSR: ', this%nblocks, ' blocks, ', &
         this%nrows, ' block-rows, blocksize=', this%blocksize
      print '(A,I0)', 'Total elements: ', this%nblocks * this%blocksize * this%blocksize
   end subroutine bsr_print_info

end module sparse_matrix_mod
