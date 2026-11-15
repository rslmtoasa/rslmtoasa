!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: lmto_radial_augmentation
!
!> @brief Angular-momentum index of an orbital in the spherical ordering.
!------------------------------------------------------------------------------
module lmto_radial_augmentation_mod

   implicit none
   private

   public :: lmto_orbital_l

contains

   !> Orbital index in the live spherical ordering: s, p(-1:1), d(-2:2), ...
   pure integer function lmto_orbital_l(iorb) result(l)
      integer, intent(in) :: iorb
      integer :: trial

      if (iorb < 1) then
         l = -1
         return
      end if
      l = -1
      do trial = 0, 8
         if (trial*trial + 1 <= iorb .and. iorb <= (trial + 1)*(trial + 1)) then
            l = trial
            return
         end if
      end do
   end function lmto_orbital_l

end module lmto_radial_augmentation_mod
