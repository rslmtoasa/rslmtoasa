!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Canonical labels for charge-spin response components.
!>
!> Cartesian components are the primary response basis.  PLUS and MINUS are
!> circular views, with m^+ = m_x + i m_y and m^- = m_x - i m_y.
module response_components_mod
   implicit none

   private

   integer, parameter, public :: RESPONSE_CHARGE = 0
   integer, parameter, public :: RESPONSE_MX = 1
   integer, parameter, public :: RESPONSE_MY = 2
   integer, parameter, public :: RESPONSE_MZ = 3
   integer, parameter, public :: RESPONSE_PLUS = 4
   integer, parameter, public :: RESPONSE_MINUS = 5

   public :: response_component_name, is_response_component

contains

   pure function response_component_name(component) result(name)
      integer, intent(in) :: component
      character(len=7) :: name

      select case (component)
      case (RESPONSE_CHARGE)
         name = 'CHARGE'
      case (RESPONSE_MX)
         name = 'MX    '
      case (RESPONSE_MY)
         name = 'MY    '
      case (RESPONSE_MZ)
         name = 'MZ    '
      case (RESPONSE_PLUS)
         name = 'PLUS  '
      case (RESPONSE_MINUS)
         name = 'MINUS '
      case default
         name = 'INVALID'
      end select
   end function response_component_name

   pure logical function is_response_component(component)
      integer, intent(in) :: component
      is_response_component = component >= RESPONSE_CHARGE .and. component <= RESPONSE_MINUS
   end function is_response_component

end module response_components_mod
