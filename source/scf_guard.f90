module scf_guard_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use precision_mod, only: rp
   use logger_mod, only: g_logger
   implicit none
   private
   public :: guard_q0_capacity, guard_nonfinite_scalar

contains

   subroutine guard_q0_capacity(atom, site, l, spin, iteration, q0)
      character(len=*), intent(in) :: atom
      integer, intent(in) :: site, l, spin, iteration
      real(rp), intent(in) :: q0
      real(rp) :: capacity
      character(len=24) :: value_text, limit_text
      character(len=8) :: shell

      capacity = real(2*l + 1, rp)
      if (.not. ieee_is_finite(q0)) then
         call guard_nonfinite_scalar(atom, site, l, spin, 'q0', iteration, q0)
      end if
      if (q0 <= capacity + 1.0e-6_rp) return

      write(value_text, '(es24.16e3)') q0
      write(limit_text, '(es24.16e3)') capacity + 1.0e-6_rp
      shell = shell_label(l)
      call g_logger%fatal('B1 occupation-capacity guard: atom='//trim(atom)//' site='//itoa(site)// &
         ' channel='//trim(shell)//' l='//itoa(l)//' spin='//itoa(spin)//' iteration='//iteration_label(iteration)// &
         ' q0='//trim(adjustl(value_text))//' capacity='//itoa(2*l + 1)// &
         ' allowed_max='//trim(adjustl(limit_text)), __FILE__, __LINE__)
   end subroutine guard_q0_capacity

   subroutine guard_nonfinite_scalar(atom, site, l, spin, field, iteration, value)
      character(len=*), intent(in) :: atom, field
      integer, intent(in) :: site, l, spin, iteration
      real(rp), intent(in) :: value
      character(len=24) :: value_text
      character(len=8) :: shell

      if (ieee_is_finite(value)) return
      write(value_text, '(es24.16e3)') value
      shell = shell_label(l)
      call g_logger%fatal('B1 non-finite guard: atom='//trim(atom)//' site='//itoa(site)// &
         ' channel='//trim(shell)//' l='//itoa(l)//' spin='//itoa(spin)//' iteration='//iteration_label(iteration)// &
         ' field='//trim(field)//' value='//trim(adjustl(value_text)), __FILE__, __LINE__)
   end subroutine guard_nonfinite_scalar

   function shell_label(l) result(label)
      integer, intent(in) :: l
      character(len=8) :: label
      select case (l)
      case (0)
         label = 's'
      case (1)
         label = 'p'
      case (2)
         label = 'd'
      case (3)
         label = 'f'
      case default
         label = 'l='//itoa(l)
      end select
   end function shell_label

   function itoa(value) result(text)
      integer, intent(in) :: value
      character(len=24) :: text
      write(text, '(i0)') value
      text = adjustl(text)
   end function itoa

   function iteration_label(iteration) result(text)
      integer, intent(in) :: iteration
      character(len=48) :: text
      if (iteration > 0) then
         write(text, '(i0)') iteration
         text = adjustl(text)
      else
         text = '0 (non-SCF/helper context)'
      end if
   end function iteration_label

end module scf_guard_mod
