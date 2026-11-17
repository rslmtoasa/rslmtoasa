module magnetic_representation_mod
   use precision_mod, only: rp
   implicit none
   private

   integer, parameter, public :: magnetic_representation_len = 24
   character(len=magnetic_representation_len), parameter, public :: periodic_nc = 'periodic_nc'
   character(len=magnetic_representation_len), parameter, public :: gbt_single_q = 'gbt_single_q'
   character(len=magnetic_representation_len), parameter, public :: explicit_texture = 'explicit_texture'

   public :: normalize_magnetic_representation
   public :: select_endpoint_moments
   public :: texture_bulk_reuse_is_valid

contains

   pure function normalize_magnetic_representation(value) result(normalized)
      character(len=*), intent(in) :: value
      character(len=magnetic_representation_len) :: normalized
      integer :: i, code

      normalized = adjustl(value)
      do i = 1, len_trim(normalized)
         code = iachar(normalized(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) normalized(i:i) = achar(code + 32)
      end do
   end function normalize_magnetic_representation

   pure subroutine select_endpoint_moments(representation, ia, ja, type_mom_i, type_mom_j, &
                                           site_moments, mom_i, mom_j, valid)
      character(len=*), intent(in) :: representation
      integer, intent(in) :: ia, ja
      real(rp), intent(in) :: type_mom_i(3), type_mom_j(3)
      real(rp), intent(in), optional :: site_moments(:, :)
      real(rp), intent(out) :: mom_i(3), mom_j(3)
      logical, intent(out) :: valid
      character(len=magnetic_representation_len) :: mode

      mode = normalize_magnetic_representation(representation)
      mom_i = type_mom_i
      mom_j = type_mom_j
      valid = .true.

      select case (trim(mode))
      case (periodic_nc)
         continue
      case (explicit_texture)
         if (.not. present(site_moments)) then
            valid = .false.
         else if (size(site_moments, 1) /= 3 .or. ia < 1 .or. ja < 1 .or. &
                  ia > size(site_moments, 2) .or. ja > size(site_moments, 2)) then
            valid = .false.
         else
            mom_i = site_moments(:, ia)
            mom_j = site_moments(:, ja)
         end if
      case default
         valid = .false.
      end select
   end subroutine select_endpoint_moments

   pure logical function texture_bulk_reuse_is_valid(site_types, representatives, site_moments, tolerance)
      integer, intent(in) :: site_types(:), representatives(:)
      real(rp), intent(in) :: site_moments(:, :), tolerance
      integer :: itype, ia, representative

      texture_bulk_reuse_is_valid = .false.
      if (size(site_moments, 1) /= 3 .or. size(site_moments, 2) < size(site_types)) return
      do itype = 1, size(representatives)
         representative = representatives(itype)
         if (representative < 1 .or. representative > size(site_moments, 2)) return
         do ia = 1, size(site_types)
            if (site_types(ia) /= itype) cycle
            if (maxval(abs(site_moments(:, ia) - site_moments(:, representative))) > tolerance) return
         end do
      end do
      texture_bulk_reuse_is_valid = .true.
   end function texture_bulk_reuse_is_valid

end module magnetic_representation_mod
