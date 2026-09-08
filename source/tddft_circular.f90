!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Explicit circular-channel contract for collinear transverse TDDFT.
!>
!> A transverse response is not a single scalar object in a multi-sublattice
!> magnet.  The ordered pairs (+,-) and (-,+) have different Lehmann poles
!> and must retain their own source/measurement vertices and output records.
module tddft_circular_mod
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   implicit none

   private

   integer, parameter, public :: TDDFT_CIRCULAR_BOTH = 0
   integer, parameter, public :: TDDFT_CIRCULAR_PLUS_MINUS = 1
   integer, parameter, public :: TDDFT_CIRCULAR_MINUS_PLUS = 2

   public :: circular_channel_name
   public :: circular_channel_components
   public :: circular_channel_code
   public :: is_circular_channel_code
   public :: opposite_circular_channel
   public :: circular_channel_file_tag

contains

   pure function circular_channel_name(channel) result(name)
      integer, intent(in) :: channel
      character(len=16) :: name

      select case (channel)
      case (TDDFT_CIRCULAR_PLUS_MINUS)
         name = 'plus_minus'
      case (TDDFT_CIRCULAR_MINUS_PLUS)
         name = 'minus_plus'
      case default
         name = 'both'
      end select
   end function circular_channel_name

   !> Return the ordered measurement/source components for one channel.
   pure subroutine circular_channel_components(channel, left_component, right_component)
      integer, intent(in) :: channel
      integer, intent(out) :: left_component, right_component

      select case (channel)
      case (TDDFT_CIRCULAR_PLUS_MINUS)
         left_component = RESPONSE_PLUS
         right_component = RESPONSE_MINUS
      case (TDDFT_CIRCULAR_MINUS_PLUS)
         left_component = RESPONSE_MINUS
         right_component = RESPONSE_PLUS
      case default
         left_component = -1
         right_component = -1
      end select
   end subroutine circular_channel_components

   pure integer function circular_channel_code(name) result(channel)
      character(len=*), intent(in) :: name
      character(len=32) :: normalized

      normalized = lower_ascii(name)
      select case (trim(normalized))
      case ('plus_minus', 'plus-minus', '+-')
         channel = TDDFT_CIRCULAR_PLUS_MINUS
      case ('minus_plus', 'minus-plus', '-+')
         channel = TDDFT_CIRCULAR_MINUS_PLUS
      case ('both', 'all')
         channel = TDDFT_CIRCULAR_BOTH
      case default
         channel = -1
      end select
   end function circular_channel_code

   pure logical function is_circular_channel_code(channel) result(valid)
      integer, intent(in) :: channel
      valid = channel == TDDFT_CIRCULAR_BOTH .or. channel == TDDFT_CIRCULAR_PLUS_MINUS .or. &
         channel == TDDFT_CIRCULAR_MINUS_PLUS
   end function is_circular_channel_code

   !> Return the opposite ordered circular correlator.  BOTH is a request for
   !> the pair, rather than an ordered response object, and therefore maps to
   !> itself so callers can validate it before constructing vertices.
   pure integer function opposite_circular_channel(channel) result(opposite)
      integer, intent(in) :: channel

      select case (channel)
      case (TDDFT_CIRCULAR_PLUS_MINUS)
         opposite = TDDFT_CIRCULAR_MINUS_PLUS
      case (TDDFT_CIRCULAR_MINUS_PLUS)
         opposite = TDDFT_CIRCULAR_PLUS_MINUS
      case default
         opposite = TDDFT_CIRCULAR_BOTH
      end select
   end function opposite_circular_channel

   !> Return the filename tag for an ordered channel.
   !>
   !> The selected `plus_minus` response retains the historical empty tag for
   !> compatibility.  Every other ordered product is explicitly labelled from
   !> the channel code, so a selected `minus_plus` primary product cannot be
   !> mistaken for the historical `plus_minus` file.
   pure function circular_channel_file_tag(channel, historical_primary) result(tag)
      integer, intent(in) :: channel
      logical, intent(in), optional :: historical_primary
      character(len=16) :: tag
      logical :: keep_historical_name

      keep_historical_name = .false.
      if (present(historical_primary)) keep_historical_name = historical_primary
      tag = ''
      select case (channel)
      case (TDDFT_CIRCULAR_PLUS_MINUS)
         if (.not. keep_historical_name) tag = '_plus_minus'
      case (TDDFT_CIRCULAR_MINUS_PLUS)
         tag = '_minus_plus'
      case default
         tag = '_unspecified'
      end select
   end function circular_channel_file_tag

   pure function lower_ascii(input) result(output)
      character(len=*), intent(in) :: input
      character(len=len(input)) :: output
      integer :: i, code

      output = input
      do i = 1, len(input)
         code = iachar(output(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) output(i:i) = achar(code + 32)
      end do
   end function lower_ascii

end module tddft_circular_mod
