!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Logger
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
!> Module to handle logging informations
!------------------------------------------------------------------------------

module logger_mod
   use iso_fortran_env, only: i8 => INT64, ERROR_UNIT, OUTPUT_UNIT
   use string_mod, only: sl, int2str
   implicit none

   type logger
      integer, private :: verbose
      character(len=245), private :: i_fmt, r4_fmt, r8_fmt

   contains
      procedure :: init => logger_init
      !Information which will only be useful for debugging, such as announcing when entering and exiting a procedure
      procedure :: debug => logger_debug
      !Information about the normal operation of the program
      procedure :: info => logger_info
      !Information produced when something happens which results in suboptimal completion of the program
      procedure :: warning => logger_warning
      !Information about something which will result in incorrect completion of the program
      procedure :: error => logger_error
      !Information that an event has occurred which will result in immediate termination of the program, without completion
      procedure :: fatal => logger_fatal
      procedure, private :: log_var => logger_log_var
      procedure, private :: message => logger_message
      procedure, private :: logger_log_integer, logger_log_real4, logger_log_real8, logger_log_character, logger_log_logical
      generic :: log => logger_log_integer, logger_log_real4, logger_log_real8, logger_log_character, logger_log_logical
   end type

   ! Only inside this module
   character(len=3), dimension(12), parameter, private :: months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

   ! Global
   type(logger) :: g_logger

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !---------------------------------------------------------------------------
   subroutine logger_init(this, verbose, i_fmt, r4_fmt, r8_fmt)
      class(logger), intent(inout) :: this
      integer, intent(in), optional :: verbose
      character(len=*), intent(in), optional :: i_fmt, r4_fmt, r8_fmt
      ! Variables
      

      ! default
      this%verbose = 0
      this%i_fmt = "(I5)"
      this%r4_fmt = "(F10.5)"
      this%r8_fmt = "(E15.10)"

      if (present(verbose)) this%verbose = verbose
      if (present(i_fmt)) this%i_fmt = i_fmt
      if (present(r4_fmt)) this%r4_fmt = r4_fmt
      if (present(r8_fmt)) this%r8_fmt = r8_fmt
   end subroutine logger_init

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !---------------------------------------------------------------------------
   subroutine logger_message(this, message, error, file, line)
      class(logger), intent(in) :: this
      character(len=*), intent(in) :: message
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      logical, intent(in), optional :: error

      integer, parameter :: file_width = 24
      integer :: unit
      character(len=file_width) :: spaces
      integer :: i_pad

      ! Additions to adjust the <info> printing
      i_pad = len(spaces) - len_trim(file) - len_trim(int2str(line))
      i_pad = max(1, i_pad)
      spaces = repeat(' ', i_pad)

      if (present(error) .and. error) then
         unit = ERROR_UNIT
      else
         unit = OUTPUT_UNIT
      end if
      write (unit, '("[", A24, "] [", A, ":", I0, "] ", A, A)') trim(current_time()), trim(file), line, spaces(1:i_pad), trim(message)
   end subroutine logger_message

   function current_time()
      character(len=24) :: current_time
      integer(i8), dimension(8) :: time_vals
      call date_and_time(values=time_vals)
      write (current_time, 1) months(time_vals(2)), time_vals(3), &
         time_vals(1), time_vals(5), time_vals(6), &
         time_vals(7), time_vals(8)
1     format(a3, 1X, i2, 1X, i4, 1X, i2.2, ":", i2.2, ":", i2.2, ".", i3.3)
   end function current_time

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief If compiled with flag -DDEBUG, prints a message to screen with level DEGUB
   !---------------------------------------------------------------------------
   subroutine logger_debug(this, message, file, line)
      class(logger), intent(in) :: this
      character(len=*), intent(in) :: message
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
#ifdef DEBUG
      call this%message(print("debug")//message, .false., file, line)
#endif
   end subroutine logger_debug

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Prints a message to screen with level INFO
   !---------------------------------------------------------------------------
   subroutine logger_info(this, message, file, line)
      class(logger), intent(in) :: this
      character(len=*), intent(in) :: message
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      character(24) :: r_file
      integer :: iidx

      ! Addition to remove full file paths
      iidx = index(file, '/', BACK=.true.)
      if (iidx > 0) then
         r_file = trim(file(iidx + 1:))
      else
         r_file = file
      end if

      call this%message(print("info")//message, .false., r_file, line)
   end subroutine logger_info

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Prints a message to screen with level WARNING
   !---------------------------------------------------------------------------
   subroutine logger_warning(this, message, file, line)
      class(logger), intent(in) :: this
      character(len=*), intent(in) :: message
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      call this%message(print("warning")//message, .false., file, line)
   end subroutine logger_warning

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Prints a message to screen with level ERROR
   !---------------------------------------------------------------------------
   subroutine logger_error(this, message, file, line)
      class(logger), intent(in) :: this
      character(len=*), intent(in) :: message
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      call this%message(print("error")//message, .true., file, line)
   end subroutine logger_error

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Prints a message to screen with level FATAL and stios the application
   !---------------------------------------------------------------------------
   subroutine logger_fatal(this, message, file, line)
      class(logger), intent(in) :: this
      character(len=*), intent(in) :: message
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      call this%message(print("fatal")//message, .true., file, line)
      error stop
   end subroutine logger_fatal

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief General subroutine to printo a variable to screen
   !---------------------------------------------------------------------------
   subroutine logger_log_var(this, varname, varval, file, line)
      ! parameters
      class(logger), intent(in) :: this
      character(len=*), intent(in) :: varname, varval
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      ! variables
      

      call this%message(print("log")//trim(varname)//" = "//trim(varval), .false., file, line)
   end subroutine logger_log_var

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Prints the value of a integer variable to screen
   !---------------------------------------------------------------------------
   subroutine logger_log_integer(this, varname, variable, file, line)
      ! parameters
      class(logger), intent(in) :: this
      character(len=*), intent(in) :: varname
      integer, intent(in) :: variable
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      ! variables
      character(len=25) :: varval

      write (varval, this%i_fmt) variable
      call this%log_var(varname, trim(varval), file, line)
   end subroutine logger_log_integer

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Prints the value of a real4 variable to screen
   !---------------------------------------------------------------------------
   subroutine logger_log_real4(this, varname, variable, file, line)
      ! parameters
      class(logger), intent(in) :: this
      character(len=*), intent(in) :: varname
      real(4), intent(in) :: variable
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      ! variables
      character(len=25) :: varval

      write (varval, this%r4_fmt) variable
      call this%log_var(varname, trim(varval), file, line)
   end subroutine logger_log_real4

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Prints the value of a real8 variable to screen
   !---------------------------------------------------------------------------
   subroutine logger_log_real8(this, varname, variable, file, line)
      ! parameters
      class(logger), intent(in) :: this
      character(len=*), intent(in) :: varname
      real(8), intent(in) :: variable
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      ! variables
      character(len=25) :: varval

      write (varval, this%r8_fmt) variable
      call this%log_var(varname, trim(varval), file, line)
   end subroutine logger_log_real8

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Prints the value of a character variable to screen
   !---------------------------------------------------------------------------
   subroutine logger_log_character(this, varname, variable, file, line)
      ! parameters
      class(logger), intent(in) :: this
      character(len=*), intent(in) :: varname
      character(len=*) :: variable
      character(len=*), intent(in) :: file
      integer, intent(in) :: line

      call this%log_var(varname, variable, file, line)
   end subroutine logger_log_character

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Prints the value of a logical variable to screen
   !---------------------------------------------------------------------------
   subroutine logger_log_logical(this, varname, variable, file, line)
      ! parameters
      class(logger), intent(in) :: this
      character(len=*), intent(in) :: varname
      logical :: variable
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      ! variables
      character(len=25) :: varval

      write (varval, *) variable
      call this%log_var(varname, trim(varval), file, line)
   end subroutine logger_log_logical

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief If compiled with flag -DCOLOR, add color to messagem
   !---------------------------------------------------------------------------
   function print(designator) result(cdesignator)
      use face, only: colorize
      implicit none
      character(len=*) :: designator
      character(len=:), allocatable :: cdesignator
      cdesignator = "<"//designator//"> "
#ifdef COLOR
      select case (cdesignator)
      case ('<debug>')
         cdesignator = colorize(cdesignator, color_fg='cyan', style='bold_on')
      case ('<log>')
         cdesignator = colorize(cdesignator, color_fg='blue', style='bold_on')
      case ('<info>')
         cdesignator = colorize(cdesignator, color_fg='green', style='bold_on')
      case ('<warning>')
         cdesignator = colorize(cdesignator, color_fg='yellow', style='bold_on')
      case ('<error>')
         cdesignator = colorize(cdesignator, color_fg='red', style='bold_on')
      case ('<fatal>')
         cdesignator = colorize(cdesignator, color_fg='red', style='bold_on')
      end select
#endif
   end function print
end module logger_mod
