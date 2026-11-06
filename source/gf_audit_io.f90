! Opt-in read-only diagnostics. No production physics lives in this module.
module gf_audit_io_mod
   use precision_mod, only: rp
   implicit none
   private
   public :: gf_audit_selected, gf_audit_stride, gf_audit_matrix, gf_audit_open, gf_audit_h_enabled
contains
   logical function gf_audit_selected(pair)
      integer, intent(in) :: pair
      character(32) :: value
      integer :: status, selected, ios
      gf_audit_selected = .false.
      call get_environment_variable('RSLMTO_GF_AUDIT', value, status=status)
      if (status /= 0) return
      if (trim(value) == 'all') then
         gf_audit_selected = pair >= 1 .and. pair <= 24
      else
         read(value, *, iostat=ios) selected
         if (ios == 0) gf_audit_selected = pair == selected .and. pair >= 1
      end if
   end function

   logical function gf_audit_h_enabled()
      character(8) :: value
      integer :: status
      call get_environment_variable('RSLMTO_GF_AUDIT_H', value, status=status)
      gf_audit_h_enabled = status == 0 .and. trim(value) == '1'
   end function

   integer function gf_audit_stride(ne)
      integer, intent(in) :: ne
      integer :: status, step, ios
      character(32) :: value
      gf_audit_stride = max(1, (ne + 62)/63)
      call get_environment_variable('RSLMTO_GF_AUDIT_STRIDE', value, status=status)
      if (status /= 0) return
      read(value, *, iostat=ios) step
      if (ios == 0) gf_audit_stride = max(gf_audit_stride, step)
   end function

   subroutine gf_audit_open(unit, kind, owner, pair, append)
      integer, intent(out) :: unit
      character(*), intent(in) :: kind
      integer, intent(in) :: owner, pair
      logical, optional, intent(in) :: append
      integer :: ios
      character(160) :: filename
      logical :: add
      add = .false.
      if (present(append)) add = append
      write(filename, '(a,a,a,i0,a,i0,a)') 'gf-audit-',kind,'-r',owner,'-p',pair,'.dat'
      if (add) then
         open(newunit=unit, file=trim(filename), status='old', position='append', action='write', iostat=ios)
      else
         open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
      end if
      if (ios /= 0) then
         unit = 0
         write(*, '(a,a)') 'GF audit could not open ', trim(filename)
      end if
   end subroutine

   subroutine gf_audit_matrix(unit, index, energy, label, matrix)
      integer, intent(in) :: unit, index
      real(rp), intent(in) :: energy
      character(*), intent(in) :: label
      complex(rp), intent(in) :: matrix(:, :)
      integer :: row, col
      write(unit, '(a,1x,i0,1x,es25.16e3,1x,a,1x,i0)') 'M',index,energy,label,size(matrix,1)
      do row = 1, size(matrix, 1)
         write(unit, '(*(es25.16e3,1x))') (real(matrix(row,col),rp), aimag(matrix(row,col)), col=1,size(matrix,2))
      end do
   end subroutine
end module gf_audit_io_mod
