module block_sparse_operator_mod

   use precision_mod, only: rp
   implicit none

   private

   type, public :: block_sparse_operator
      integer :: block_dim = 0
      integer :: n_sites = 0
      integer :: nnzb = 0
      character(len=32) :: kind = 'none'
      logical :: includes_overlap = .false.
      logical :: includes_soc = .false.
      logical :: includes_enim = .false.
      integer, allocatable :: row_ptr(:)
      integer, allocatable :: col_ind(:)
      integer, allocatable :: site_types(:)
      complex(rp), allocatable :: blocks(:, :, :)
   contains
      procedure :: clear
      procedure :: is_ready
      procedure :: write_text
   end type block_sparse_operator

contains

   subroutine clear(this)
      class(block_sparse_operator), intent(inout) :: this

      if (allocated(this%row_ptr)) deallocate (this%row_ptr)
      if (allocated(this%col_ind)) deallocate (this%col_ind)
      if (allocated(this%site_types)) deallocate (this%site_types)
      if (allocated(this%blocks)) deallocate (this%blocks)

      this%block_dim = 0
      this%n_sites = 0
      this%nnzb = 0
      this%kind = 'none'
      this%includes_overlap = .false.
      this%includes_soc = .false.
      this%includes_enim = .false.
   end subroutine clear

   logical function is_ready(this)
      class(block_sparse_operator), intent(in) :: this

      is_ready = allocated(this%row_ptr) .and. allocated(this%col_ind) .and. allocated(this%blocks)
   end function is_ready

   subroutine write_text(this, path)
      class(block_sparse_operator), intent(in) :: this
      character(len=*), intent(in) :: path
      integer :: unitno, ios, row, entry, i, j

      open (newunit=unitno, file=path, status='replace', action='write', iostat=ios)
      if (ios /= 0) return

      write (unitno, '(A)') '# block_sparse_operator'
      write (unitno, '(A,1X,A)') 'kind', trim(this%kind)
      write (unitno, '(A,1X,I0)') 'block_dim', this%block_dim
      write (unitno, '(A,1X,I0)') 'n_sites', this%n_sites
      write (unitno, '(A,1X,I0)') 'nnzb', this%nnzb
      write (unitno, '(A,1X,L1)') 'includes_overlap', this%includes_overlap
      write (unitno, '(A,1X,L1)') 'includes_soc', this%includes_soc
      write (unitno, '(A,1X,L1)') 'includes_enim', this%includes_enim

      if (allocated(this%site_types)) then
         write (unitno, '(A)') '# site_types(site, type)'
         do row = 1, size(this%site_types)
            write (unitno, '(I0,1X,I0)') row, this%site_types(row)
         end do
      end if

      write (unitno, '(A)') '# row_ptr'
      do row = 1, size(this%row_ptr)
         write (unitno, '(I0,1X,I0)') row, this%row_ptr(row)
      end do

      write (unitno, '(A)') '# entries: row col bi bj real imag'
      do row = 1, this%n_sites
         do entry = this%row_ptr(row), this%row_ptr(row + 1) - 1
            do i = 1, this%block_dim
               do j = 1, this%block_dim
                  write (unitno, '(4(I0,1X),2(ES24.16,1X))') row, this%col_ind(entry), i, j, &
                     real(this%blocks(i, j, entry)), aimag(this%blocks(i, j, entry))
               end do
            end do
         end do
      end do

      close (unitno)
   end subroutine write_text

end module block_sparse_operator_mod
