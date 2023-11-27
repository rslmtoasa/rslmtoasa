!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Potential
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module manipulate arrays
!------------------------------------------------------------------------------

module array_mod
   use logger_mod, only: g_logger

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief pushes an element to array
   !
   !> @param[inout] arr vector push
   !> @param[in] new data to push
   !> @param[in, optional] index to push new element
   !---------------------------------------------------------------------------
   interface push
      procedure :: push_c
      procedure :: push_i
      procedure :: push_r4
      procedure :: push_r8
   end interface push

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief pops an element to array
   !
   !> @param[inout] arr vector pop
   !> @param[in, optional] index to pop new element
   !---------------------------------------------------------------------------
   interface pop
      procedure :: pop_c
      procedure :: pop_i
      procedure :: pop_r4
      procedure :: pop_r8
   end interface pop

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief swaps two elements in array
   !
   !> @param[inout] arr vector to swap indexes
   !> @param[in] i first index
   !> @param[in] j second index
   !---------------------------------------------------------------------------
   interface swap
      procedure :: swap_c
      procedure :: swap_i
      procedure :: swap_r4
      procedure :: swap_r8
   end interface swap

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief sorts array using quick-sort method
   !
   !> @param[inout] arr array to sort
   !---------------------------------------------------------------------------
   interface qsort
      procedure :: qsort_c
      procedure :: qsort_i
      procedure :: qsort_r4
      procedure :: qsort_r8
   end interface qsort

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief same as qsort subroutine, but a function instead
   !
   !> @param[in] arr array to sort
   !---------------------------------------------------------------------------
   interface fqsort
      procedure :: fqsort_c
      procedure :: fqsort_i
      procedure :: fqsort_r4
      procedure :: fqsort_r8
   end interface fqsort

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief returns indexes of sorted array using quick-sort method
   !
   !> @param[in] arr array to get ordered indexes
   !---------------------------------------------------------------------------
   interface fqsortloc
      procedure :: fqsortloc_c
      procedure :: fqsortloc_i
      procedure :: fqsortloc_r4
      procedure :: fqsortloc_r8
   end interface fqsortloc

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief returns uniques elements of array
   !
   !> @param[in] arr array to get uniques
   !---------------------------------------------------------------------------
   interface unique
      procedure :: unique_c
      procedure :: unique_i
      procedure :: unique_r4
      procedure :: unique_r8
   end interface unique

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief returns uniques elements of array
   !
   !> @param[in] arr array to get uniques
   !---------------------------------------------------------------------------
   interface funique
      procedure :: funique_c
      procedure :: funique_i
      procedure :: funique_r4
      procedure :: funique_r8
   end interface funique

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine push_c(arr, new, index)
      character(len=*), dimension(:), allocatable, intent(inout) :: arr
      character(len=:), dimension(:), allocatable :: tmp
      character(len=*), intent(in) :: new
      include 'include_codes/array/push.f90'
   end subroutine push_c

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine push_i(arr, new, index)
      integer, dimension(:), allocatable, intent(inout) :: arr
      integer, dimension(:), allocatable :: tmp
      integer, intent(in) :: new
      include 'include_codes/array/push.f90'
   end subroutine push_i

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine push_r4(arr, new, index)
      real(4), dimension(:), allocatable, intent(inout) :: arr
      real(4), dimension(:), allocatable :: tmp
      real(4), intent(in) :: new
      include 'include_codes/array/push.f90'
   end subroutine push_r4

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine push_r8(arr, new, index)
      real(8), dimension(:), allocatable, intent(inout) :: arr
      real(8), dimension(:), allocatable :: tmp
      real(8), intent(in) :: new
      include 'include_codes/array/push.f90'
   end subroutine push_r8

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine pop_c(arr, index)
      character(len=*), dimension(:), allocatable, intent(inout) :: arr
      character(len=:), dimension(:), allocatable :: tmp
      include 'include_codes/array/pop.f90'
   end subroutine pop_c

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine pop_i(arr, index)
      integer, dimension(:), allocatable, intent(inout) :: arr
      integer, dimension(:), allocatable :: tmp
      include 'include_codes/array/pop.f90'
   end subroutine pop_i

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine pop_r4(arr, index)
      real(4), dimension(:), allocatable, intent(inout) :: arr
      real(4), dimension(:), allocatable :: tmp
      include 'include_codes/array/pop.f90'
   end subroutine pop_r4

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine pop_r8(arr, index)
      real(8), dimension(:), allocatable, intent(inout) :: arr
      real(8), dimension(:), allocatable :: tmp
      include 'include_codes/array/pop.f90'
   end subroutine pop_r8

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine swap_c(arr, i, j)
      character(len=*), dimension(:), intent(inout) :: arr
      character(len=:), allocatable :: tmp
      include 'include_codes/array/swap.f90'
   end subroutine swap_c

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine swap_i(arr, i, j)
      integer, dimension(:), intent(inout) :: arr
      integer :: tmp
      include 'include_codes/array/swap.f90'
   end subroutine swap_i

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine swap_r4(arr, i, j)
      real(4), dimension(:), intent(inout) :: arr
      real(4) :: tmp
      include 'include_codes/array/swap.f90'
   end subroutine swap_r4

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine swap_r8(arr, i, j)
      real(8), dimension(:), intent(inout) :: arr
      real(8) :: tmp
      include 'include_codes/array/swap.f90'
   end subroutine swap_r8

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine qsort_c(arr)
      character(len=*), dimension(:), intent(inout) :: arr
      include 'include_codes/array/qsort.f90'
   end subroutine qsort_c

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine qsort_i(arr)
      integer, dimension(:), intent(inout) :: arr
      include 'include_codes/array/qsort.f90'
   end subroutine qsort_i

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine qsort_r4(arr)
      real(4), dimension(:), intent(inout) :: arr
      include 'include_codes/array/qsort.f90'
   end subroutine qsort_r4

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine qsort_r8(arr)
      real(8), dimension(:), intent(inout) :: arr
      include 'include_codes/array/qsort.f90'
   end subroutine qsort_r8

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   function fqsort_c(arr) result(arr_ord)
      character(len=*), dimension(:), intent(in) :: arr
      character(len=:), dimension(:), allocatable :: arr_ord
      call qsort(arr_ord)
   end function fqsort_c

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   function fqsort_i(arr) result(arr_ord)
      integer, dimension(:), intent(in) :: arr
      integer, dimension(:), allocatable :: arr_ord
      call qsort(arr_ord)
   end function fqsort_i

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   function fqsort_r4(arr) result(arr_ord)
      real(4), dimension(:), intent(in) :: arr
      real(4), dimension(:), allocatable :: arr_ord
      call qsort(arr_ord)
   end function fqsort_r4

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   function fqsort_r8(arr) result(arr_ord)
      real(8), dimension(:), intent(in) :: arr
      real(8), dimension(:), allocatable :: arr_ord
      call qsort(arr_ord)
   end function fqsort_r8

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   function fqsortloc_c(arr) result(arr_idx)
      character(len=*), allocatable, dimension(:), intent(in) :: arr
      include 'include_codes/array/qsortloc.f90'
   end function fqsortloc_c

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   function fqsortloc_i(arr) result(arr_idx)
      integer, dimension(:), intent(in) :: arr
      include 'include_codes/array/qsortloc.f90'
   end function fqsortloc_i

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   function fqsortloc_r4(arr) result(arr_idx)
      real(4), dimension(:), intent(in) :: arr
      include 'include_codes/array/qsortloc.f90'
   end function fqsortloc_r4

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   function fqsortloc_r8(arr) result(arr_idx)
      real(8), dimension(:), intent(in) :: arr
      include 'include_codes/array/qsortloc.f90'
   end function fqsortloc_r8

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine unique_c(arr)
      character(len=*), dimension(:), allocatable, intent(inout) :: arr
      character(len=:), dimension(:), allocatable :: tmp
      include 'include_codes/array/unique.f90'
   end subroutine unique_c

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine unique_i(arr)
      integer, dimension(:), allocatable, intent(inout) :: arr
      integer, dimension(:), allocatable :: tmp
      include 'include_codes/array/unique.f90'
   end subroutine unique_i

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine unique_r4(arr)
      real(4), dimension(:), allocatable, intent(inout) :: arr
      real(4), dimension(:), allocatable :: tmp
      include 'include_codes/array/unique.f90'
   end subroutine unique_r4

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   subroutine unique_r8(arr)
      real(8), dimension(:), allocatable, intent(inout) :: arr
      real(8), dimension(:), allocatable :: tmp
      include 'include_codes/array/unique.f90'
   end subroutine unique_r8

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   function funique_c(arr) result(tmp)
      character(len=*), dimension(:), allocatable, intent(in) :: arr
      character(len=:), dimension(:), allocatable :: tmp
      call g_logger%fatal('There is no overload for this function consider using the subroutine.', __FILE__, __LINE__)
   end function funique_c

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   function funique_i(arr) result(tmp)
      integer, dimension(:), intent(in) :: arr
      integer, dimension(:), allocatable :: tmp
      include 'include_codes/array/funique.f90'
   end function funique_i

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   function funique_r4(arr) result(tmp)
      real(4), dimension(:), intent(in) :: arr
      real(4), dimension(:), allocatable :: tmp
      include 'include_codes/array/funique.f90'
   end function funique_r4

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !---------------------------------------------------------------------------
   function funique_r8(arr) result(tmp)
      real(8), dimension(:), intent(in) :: arr
      real(8), dimension(:), allocatable :: tmp
      include 'include_codes/array/funique.f90'
   end function funique_r8
end module array_mod
