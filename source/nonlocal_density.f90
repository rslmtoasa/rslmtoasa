!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: nonlocal_density_mod
!
!> @brief Non-local density matrix from the inter-site Chebyshev moments of the SCF recursion.
!> @details Holds the (site, neighbour) pairs requested from the one-shot recursion
!>          when control%nonlocal_dm is set.
!------------------------------------------------------------------------------
module nonlocal_density_mod
   use lattice_mod, only: lattice
   implicit none
   private

   type, public :: nonlocal_density
      !> Requested (i, j) pairs of the local recursion sites, self pair first.
      integer, allocatable :: pairs(:, :)
   contains
      procedure :: neighbour_pairs
      procedure :: restore_to_default
      final :: destructor
   end type nonlocal_density

contains

   !> For each recursion site i = first..last: (irec(i), irec(i)), then (irec(i), j) for every
   !> neighbour j = nn(irec(i), 2:nn(irec(i),1)) in table order, zeros and repeated entries skipped.
   subroutine neighbour_pairs(this, lat, first, last)
      class(nonlocal_density), intent(inout) :: this
      type(lattice), intent(in) :: lat
      integer, intent(in) :: first, last
      integer :: i, m, k, n, site, j
      integer, allocatable :: tmp(:, :)

      call this%restore_to_default()
      allocate (tmp(max(0, last - first + 1)*size(lat%nn, 2), 2))
      n = 0
      do i = first, last
         site = lat%irec(i)
         n = n + 1
         tmp(n, :) = [site, site]
         k = n
         do m = 2, lat%nn(site, 1)
            j = lat%nn(site, m)
            if (j == 0) cycle
            if (any(tmp(k:n, 2) == j)) cycle
            n = n + 1
            tmp(n, :) = [site, j]
         end do
      end do
      this%pairs = tmp(1:n, :)
   end subroutine neighbour_pairs

   subroutine restore_to_default(this)
      class(nonlocal_density), intent(inout) :: this
      if (allocated(this%pairs)) deallocate (this%pairs)
   end subroutine restore_to_default

   subroutine destructor(this)
      type(nonlocal_density), intent(inout) :: this
      call this%restore_to_default()
   end subroutine destructor

end module nonlocal_density_mod
