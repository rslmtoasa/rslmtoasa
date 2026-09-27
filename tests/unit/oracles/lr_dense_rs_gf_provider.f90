!------------------------------------------------------------------------------
! Dense coefficient-space GF provider for finite exact test oracles.
!------------------------------------------------------------------------------
module lr_dense_rs_gf_provider_mod
   use precision_mod, only: rp
   use linear_response_mod, only: lr_rs_gf_provider
   implicit none
   private

   real(rp), parameter :: translation_tolerance = 2.0e-13_rp

   type, extends(lr_rs_gf_provider), public :: lr_rs_dense_gf_provider
      complex(rp), allocatable :: hamiltonian(:,:), hgamma(:,:)
      integer :: nsite = 0, block_size = 0
      logical :: translation_independent = .false.
   contains
      procedure :: initialize => lr_rs_dense_initialize
      procedure :: get_pair => lr_rs_dense_get_pair
      procedure :: describe => lr_rs_dense_describe
   end type lr_rs_dense_gf_provider

contains

   subroutine lr_rs_dense_initialize(this, hamiltonian, hgamma, block_size, translation_independent)
      class(lr_rs_dense_gf_provider), intent(out) :: this
      complex(rp), intent(in) :: hamiltonian(:,:), hgamma(:,:)
      integer, intent(in) :: block_size
      logical, intent(in), optional :: translation_independent

      if (block_size < 1 .or. size(hamiltonian, 1) /= size(hamiltonian, 2) .or. &
          any(shape(hgamma) /= shape(hamiltonian)) .or. mod(size(hamiltonian, 1), block_size) /= 0) then
         error stop 'lr_rs_dense_gf_provider%initialize: inconsistent Hamiltonian/block dimensions'
      end if
      this%hamiltonian = hamiltonian
      this%hgamma = hgamma
      this%block_size = block_size
      this%nsite = size(hamiltonian, 1)/block_size
      this%nbasis_per_site = block_size
      this%translation_independent = .false.
      if (present(translation_independent)) this%translation_independent = translation_independent
      this%provider_kind = 'dense_exact_coefficient_inverse'
   end subroutine lr_rs_dense_initialize

   subroutine lr_rs_dense_get_pair(this, left_site, right_site, translation, z, gij, gji, hgamma_ij, hgamma_ji)
      class(lr_rs_dense_gf_provider), intent(inout) :: this
      integer, intent(in) :: left_site, right_site
      real(rp), intent(in) :: translation(3)
      complex(rp), intent(in) :: z
      complex(rp), intent(out) :: gij(:,:), gji(:,:), hgamma_ij(:,:), hgamma_ji(:,:)
      complex(rp), allocatable :: matrix(:,:), inverse(:,:)
      integer :: left_begin, right_begin, left_end, right_end

      if (.not. allocated(this%hamiltonian)) error stop 'dense coefficient provider: not initialized'
      if (left_site < 1 .or. left_site > this%nsite .or. right_site < 1 .or. right_site > this%nsite) then
         error stop 'dense coefficient provider: site index is outside the finite system'
      end if
      if (maxval(abs(translation)) > translation_tolerance .and. .not. this%translation_independent) then
         error stop 'dense coefficient provider: nonzero translation is unsupported by the finite oracle'
      end if
      if (any(shape(gij) /= [this%block_size, this%block_size]) .or. &
          any(shape(gji) /= [this%block_size, this%block_size]) .or. &
          any(shape(hgamma_ij) /= [this%block_size, this%block_size]) .or. &
          any(shape(hgamma_ji) /= [this%block_size, this%block_size])) then
         error stop 'dense coefficient provider: output block shape mismatch'
      end if

      allocate(matrix(size(this%hamiltonian, 1), size(this%hamiltonian, 2)), &
         inverse(size(this%hamiltonian, 1), size(this%hamiltonian, 2)))
      matrix = -this%hamiltonian
      do left_begin = 1, size(matrix, 1)
         matrix(left_begin, left_begin) = z - this%hamiltonian(left_begin, left_begin)
      end do
      call invert_dense_matrix(matrix, inverse)

      left_begin = (left_site - 1)*this%block_size + 1
      left_end = left_site*this%block_size
      right_begin = (right_site - 1)*this%block_size + 1
      right_end = right_site*this%block_size
      gij = inverse(left_begin:left_end, right_begin:right_end)
      gji = inverse(right_begin:right_end, left_begin:left_end)
      hgamma_ij = this%hgamma(left_begin:left_end, right_begin:right_end)
      hgamma_ji = this%hgamma(right_begin:right_end, left_begin:left_end)
      deallocate(matrix, inverse)
   end subroutine lr_rs_dense_get_pair

   function lr_rs_dense_describe(this) result(description)
      class(lr_rs_dense_gf_provider), intent(in) :: this
      character(len=256) :: description

      write (description, '(a,i0,a,i0,a,l1)') 'provider=dense_exact_coefficient_inverse nsite=', this%nsite, &
         ' block_size=', this%block_size, ' translation_independent=', this%translation_independent
   end function lr_rs_dense_describe

   subroutine invert_dense_matrix(matrix, inverse)
      complex(rp), intent(in) :: matrix(:,:)
      complex(rp), intent(out) :: inverse(:,:)
      complex(rp), allocatable :: work(:,:), row(:)
      complex(rp) :: pivot, factor
      integer :: n, i, j, pivot_row

      if (size(matrix, 1) /= size(matrix, 2) .or. any(shape(inverse) /= shape(matrix))) then
         error stop 'invert_dense_matrix: square shape required'
      end if
      n = size(matrix, 1)
      allocate(work(n,n), row(n))
      work = matrix
      inverse = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         inverse(i,i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      do i = 1, n
         pivot_row = i
         do j = i + 1, n
            if (abs(work(j,i)) > abs(work(pivot_row,i))) pivot_row = j
         end do
         if (abs(work(pivot_row,i)) <= 100.0_rp*epsilon(1.0_rp)) then
            error stop 'invert_dense_matrix: singular finite fixture Hamiltonian resolvent'
         end if
         if (pivot_row /= i) then
            row = work(i,:); work(i,:) = work(pivot_row,:); work(pivot_row,:) = row
            row = inverse(i,:); inverse(i,:) = inverse(pivot_row,:); inverse(pivot_row,:) = row
         end if
         pivot = work(i,i)
         work(i,:) = work(i,:)/pivot
         inverse(i,:) = inverse(i,:)/pivot
         do j = 1, n
            if (j == i) cycle
            factor = work(j,i)
            work(j,:) = work(j,:) - factor*work(i,:)
            inverse(j,:) = inverse(j,:) - factor*inverse(i,:)
         end do
      end do
      deallocate(work, row)
   end subroutine invert_dense_matrix

end module lr_dense_rs_gf_provider_mod
