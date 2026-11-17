! Runs the sub-cases of a consolidated unit test, each in its own process.
!
! Protocol of a consolidated executable:  exe [all|<case>] [extra arguments]
!   all (or no argument): the parent runs every case as `exe <case> <extra>`, prints
!       PASS/FAIL <case>, and stops with status 1 if any case failed.
!   <case>: the child runs that case only.
! A case that expects a fatal stop passes only if the child exits non-zero and its output
! does not contain 'UNEXPECTED'.
module subcase_runner_mod
   implicit none
   private
   public :: subcase_case_name, run_all_subcases

contains

   ! Name of the case requested on the command line; 'all' when none.
   function subcase_case_name() result(name)
      character(len=64) :: name
      call get_command_argument(1, name)
      if (len_trim(name) == 0) name = 'all'
   end function subcase_case_name

   ! workdirs(i) is a directory to enter before running case i ('' for the current one).
   subroutine run_all_subcases(cases, expect_fatal, workdirs)
      character(len=*), intent(in) :: cases(:), workdirs(:)
      logical, intent(in) :: expect_fatal(:)
      character(len=1024) :: exe, arg, extra, cmd, line
      character(len=:), allocatable :: logname, logpath
      integer :: i, k, es, cs, unit_log, ios, nfail
      logical :: ok, unexpected

      call get_command_argument(0, exe)
      extra = ''
      do k = 2, command_argument_count()
         call get_command_argument(k, arg)
         extra = trim(extra)//' "'//trim(arg)//'"'
      end do

      nfail = 0
      do i = 1, size(cases)
         logname = 'subcase_'//trim(cases(i))//'.log'
         logpath = logname
         cmd = '"'//trim(exe)//'" '//trim(cases(i))//trim(extra)//' > "'//logname//'" 2>&1'
         if (len_trim(workdirs(i)) > 0) then
            ! The child writes its log inside the directory it enters.
            cmd = 'cd "'//trim(workdirs(i))//'" && '//trim(cmd)
            logpath = trim(workdirs(i))//'/'//logname
         end if
         es = 0
         cs = 0
         call execute_command_line(trim(cmd), wait=.true., exitstat=es, cmdstat=cs)

         unexpected = .false.
         open (newunit=unit_log, file=logpath, status='old', action='read', iostat=ios)
         if (ios == 0) then
            do
               read (unit_log, '(a)', iostat=ios) line
               if (ios /= 0) exit
               if (index(line, 'UNEXPECTED') > 0) unexpected = .true.
            end do
            close (unit_log)
         end if

         if (expect_fatal(i)) then
            ok = (cs == 0 .and. es /= 0 .and. .not. unexpected)
         else
            ok = (cs == 0 .and. es == 0)
         end if

         if (ok) then
            write (*, '(a,a)') 'PASS ', trim(cases(i))
         else
            nfail = nfail + 1
            write (*, '(a,a,a,i0,a,i0,a)') 'FAIL ', trim(cases(i)), ' (exit ', es, ', cmdstat ', cs, ')'
            call show_log(logpath)
         end if
      end do

      write (*, '(a,i0,a,i0,a)') 'sub-cases: ', size(cases) - nfail, ' passed, ', nfail, ' failed'
      if (nfail > 0) then
         write (*, '(a)') 'RESULT: FAIL'
         error stop 1
      end if
      write (*, '(a)') 'RESULT: PASS'
   end subroutine run_all_subcases

   subroutine show_log(logpath)
      character(len=*), intent(in) :: logpath
      character(len=1024) :: line
      integer :: unit_log, ios
      open (newunit=unit_log, file=logpath, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read (unit_log, '(a)', iostat=ios) line
         if (ios /= 0) exit
         write (*, '(a,a)') '   | ', trim(line)
      end do
      close (unit_log)
   end subroutine show_log

end module subcase_runner_mod
