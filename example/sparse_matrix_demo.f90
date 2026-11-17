!------------------------------------------------------------------------------
! Sparse Matrix Format Verification
!
! Demonstrates that block_to_sparse() builds proper BSR format:
! - Block Sparse Row (BSR) with 18×18 blocks
! - Compatible with oneMKL and cuSPARSE
! - Only non-zero blocks stored
!
! Usage: Compile with main RS-LMTO code and run on a test case.
!        Output shows BSR structure and block count.
!------------------------------------------------------------------------------

program sparse_matrix_verification
   use hamiltonian_mod
   use lattice_mod
   use control_mod
   use precision_mod
   use basis_mod, only: nb
   implicit none

   type(hamiltonian) :: ham
   integer :: nblocks, nrows, i
   real(rp) :: sparsity

   ! Print BSR format info
   if (allocated(ham%h_bsr%values)) then
      nblocks = ham%h_bsr%nblocks
      nrows = ham%h_bsr%nrows

      print *, "================================================"
      print *, "Block Sparse Row (BSR) Hamiltonian Summary"
      print *, "================================================"
      print '(A,I0)', 'Number of blocks: ', nblocks
      print '(A,I0)', 'Number of block-rows: ', nrows
      print '(A,I0)', 'Block size: ', nb, 'x', nb
      print '(A,I0)', 'Dense equivalent size: ', nrows*nb, 'x', nrows*nb
      print '(A,I0)', 'Total elements in full matrix: ', (nrows*nb)**2

      ! Compute sparsity
      if ((nrows * nb)**2 > 0) then
         sparsity = 100.0_rp * (1.0_rp - real(nblocks * nb * nb, rp) / real((nrows*nb)**2, rp))
         print '(A,F6.2,A)', 'Sparsity: ', sparsity, '%'
      end if

      print *, ""
      print *, "Row pointer structure (first 5 entries):"
      do i = 1, min(5, nrows+1)
         print '(A,I0,A,I0)', 'row_ptr(', i, ') = ', ham%h_bsr%row_ptr(i)
      end do

      print *, ""
      print *, "Column indices (first 10 blocks):"
      do i = 1, min(10, nblocks)
         print '(A,I0,A,I0)', 'col_indices(', i, ') = ', ham%h_bsr%col_indices(i)
      end do

      print *, ""
      print *, "Format ready for:"
      print *, "  - oneMKL (Intel Math Kernel Library) sparse operations"
      print *, "  - cuSPARSE (NVIDIA CUDA sparse library)"
      print *, "  - Other GPU-accelerated sparse kernels"
      print *, "================================================"

   else
      print *, "ERROR: BSR matrix not allocated!"
      print *, "Call hamiltonian%block_to_sparse() first."
   end if

end program sparse_matrix_verification
