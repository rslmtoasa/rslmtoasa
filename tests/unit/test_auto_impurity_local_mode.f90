program test_auto_impurity_local_mode
   use math_mod, only: auto_impurity_local_host_eligible
   implicit none
   logical :: eligible
   character(len=128) :: reason
   logical :: frozen_u_host(3)
   character(len=16) :: mode(3)
   character(len=512) :: deck_path, host_path, surface_deck_path

   call get_command_argument(1, deck_path)
   call get_command_argument(2, host_path)
   call get_command_argument(3, surface_deck_path)
   call check(len_trim(deck_path) > 0 .and. len_trim(host_path) > 0 .and. &
      len_trim(surface_deck_path) > 0, 'fixture paths supplied')
   call check(file_contains(trim(deck_path), "calctype = 'I'") .and. &
      file_contains(trim(deck_path), "hubbard_radial_moments = 'auto'") .and. &
      file_contains(trim(deck_path), 'hubbard_u_impurity(1,3)'), &
      'impurity override fixture selects auto mode')
   call check(file_contains(trim(host_path), "hubbard_radial_moments = 'local'") .and. &
      file_contains(trim(host_path), 'hubbard_u ='), 'frozen host fixture carries local U provenance')
   call check(file_contains(trim(surface_deck_path), "calctype = 'S'") .and. &
      file_contains(trim(surface_deck_path), "hubbard_radial_moments = 'auto'") .and. &
      file_contains(trim(surface_deck_path), 'hubbard_u(1,3)'), &
      'surface fixture carries an active d-shell U override in auto mode')
   frozen_u_host = [.true., .false., .false.]
   mode = ['local           ', 'spectral        ', 'spectral        ']
   call auto_impurity_local_host_eligible(.true., frozen_u_host, mode, eligible, reason)
   call check(eligible .and. index(reason, 'provenance') == 0, &
      'impurity override on a local U-carrying host resolves local without provenance error')

   ! Surface uses the same fixed-U frozen-host eligibility rule as impurity.
   call auto_impurity_local_host_eligible(.true., frozen_u_host, mode, eligible, reason)
   call check(eligible .and. index(reason, 'provenance') == 0, &
      'surface auto mode resolves local on a local U-carrying frozen host')

   frozen_u_host = [.true., .true., .false.]
   mode = ['local           ', 'local           ', 'spectral        ']
   call auto_impurity_local_host_eligible(.true., frozen_u_host, mode, eligible, reason)
   call check(eligible, 'every U-carrying frozen host is local')

   mode(2) = 'spectral'
   call auto_impurity_local_host_eligible(.true., frozen_u_host, mode, eligible, reason)
   call check(.not. eligible .and. index(reason, 'spectral') > 0, &
      'spectral frozen host keeps auto mode spectral with a reason')

   mode(2) = ''
   call auto_impurity_local_host_eligible(.true., frozen_u_host, mode, eligible, reason)
   call check(.not. eligible .and. index(reason, 'legacy provenance') > 0, &
      'legacy frozen host keeps auto mode spectral with a reason')

   mode(2) = 'local'
   call auto_impurity_local_host_eligible(.false., frozen_u_host, mode, eligible, reason)
   call check(.not. eligible .and. index(reason, 'scope') > 0, &
      'guarded or unsupported impurity scope remains ineligible')

   frozen_u_host = .false.
   call auto_impurity_local_host_eligible(.true., frozen_u_host, mode, eligible, reason)
   call check(.not. eligible .and. index(reason, 'no U-carrying') > 0, &
      'no U-carrying frozen host keeps the default spectral')

   write(*,'(a)') 'RESULT: PASS'
contains
   logical function file_contains(path, needle)
      character(len=*), intent(in) :: path, needle
      character(len=1024) :: line
      integer :: unit, ios
      file_contains = .false.
      open(newunit=unit, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, needle) > 0) file_contains = .true.
      end do
      close(unit)
   end function file_contains

   subroutine check(ok, label)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: label
      if (.not. ok) then
         write(*,'(a)') 'FAIL: '//label
         error stop 1
      end if
   end subroutine check
end program test_auto_impurity_local_mode
