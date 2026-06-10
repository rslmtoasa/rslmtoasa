module strux_errors
!- Minimal fatal-error helpers kept separate from numerical kernels.
    implicit none
    private

    public :: rx
    public :: rxi

contains

subroutine rx(msg)
!- Abort execution with a plain string message.
!  This is the minimal replacement for the richer fatal-error helpers in the
!  original LMTO47/Questaal code base.
    implicit none
    character(len=*), intent(in) :: msg
    print *, 'ERROR: ', trim(msg)
    stop 1
end subroutine rx

subroutine rxi(msg, ival)
!- Abort execution with a string message followed by one integer value.
!  Used in places where the original code reported a failing size or index.
    implicit none
    character(len=*), intent(in) :: msg
    integer, intent(in) :: ival
    print *, 'ERROR: ', trim(msg), ival
    stop 1
end subroutine rxi

end module strux_errors
