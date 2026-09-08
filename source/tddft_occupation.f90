!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Validated thermodynamic state shared by TDDFT response consumers.
module tddft_occupation_mod
   use precision_mod, only: rp
   implicit none

   private

   ! This sentinel is an implementation detail.  It must never reach response
   ! metadata because an unresolved occupation state is rejected first.
   real(rp), parameter, public :: tddft_unresolved_fermi_level = huge(1.0_rp)

   type, public :: tddft_response_occupation_state
      real(rp) :: fermi_level = tddft_unresolved_fermi_level
      real(rp) :: electronic_temperature = -1.0_rp
      real(rp) :: occupation_tolerance = -1.0_rp
      integer :: band_first = 0
      integer :: band_last = -1
      logical :: fermi_is_resolved = .false.
      character(len=48) :: fermi_source = 'unresolved'
      character(len=48) :: fermi_policy = 'unresolved'
   contains
      procedure :: set => set_response_occupation
      procedure :: validate => validate_response_occupation
      procedure :: same_as => same_response_occupation
   end type tddft_response_occupation_state

   public :: validate_response_occupation_fields

contains

   subroutine set_response_occupation(this, fermi_level, electronic_temperature, occupation_tolerance, &
      band_first, band_last, fermi_source, fermi_policy)
      class(tddft_response_occupation_state), intent(inout) :: this
      real(rp), intent(in) :: fermi_level, electronic_temperature, occupation_tolerance
      integer, intent(in) :: band_first, band_last
      character(len=*), intent(in) :: fermi_source, fermi_policy

      this%fermi_level = fermi_level
      this%electronic_temperature = electronic_temperature
      this%occupation_tolerance = occupation_tolerance
      this%band_first = band_first
      this%band_last = band_last
      this%fermi_is_resolved = .true.
      this%fermi_source = trim(fermi_source)
      this%fermi_policy = trim(fermi_policy)
      call this%validate('set_response_occupation')
   end subroutine set_response_occupation

   subroutine validate_response_occupation(this, context)
      class(tddft_response_occupation_state), intent(in) :: this
      character(len=*), intent(in) :: context

      if (.not. this%fermi_is_resolved .or. abs(this%fermi_level) >= tddft_unresolved_fermi_level/2.0_rp) then
         error stop trim(context)//': response occupation state is unresolved (Fermi level is missing)'
      end if
      call validate_response_occupation_fields(this%fermi_level, this%electronic_temperature, &
         this%occupation_tolerance, this%band_first, this%band_last, context)
      if (len_trim(this%fermi_source) == 0 .or. len_trim(this%fermi_policy) == 0) then
         error stop trim(context)//': response occupation state has no Fermi-level provenance'
      end if
   end subroutine validate_response_occupation

   logical function same_response_occupation(this, other) result(same)
      class(tddft_response_occupation_state), intent(in) :: this
      type(tddft_response_occupation_state), intent(in) :: other

      same = this%fermi_is_resolved .eqv. other%fermi_is_resolved .and. &
         this%fermi_level == other%fermi_level .and. &
         this%electronic_temperature == other%electronic_temperature .and. &
         this%occupation_tolerance == other%occupation_tolerance .and. &
         this%band_first == other%band_first .and. this%band_last == other%band_last .and. &
         trim(this%fermi_source) == trim(other%fermi_source) .and. &
         trim(this%fermi_policy) == trim(other%fermi_policy)
   end function same_response_occupation

   subroutine validate_response_occupation_fields(fermi_level, electronic_temperature, occupation_tolerance, &
      band_first, band_last, context)
      real(rp), intent(in) :: fermi_level, electronic_temperature, occupation_tolerance
      integer, intent(in) :: band_first, band_last
      character(len=*), intent(in) :: context

      if (abs(fermi_level) >= tddft_unresolved_fermi_level/2.0_rp) then
         error stop trim(context)//': response occupation state is unresolved (Fermi level is missing)'
      end if
      if (electronic_temperature < 0.0_rp) then
         error stop trim(context)//': response occupation state has an invalid electronic temperature'
      end if
      if (occupation_tolerance < 0.0_rp) then
         error stop trim(context)//': response occupation state has an invalid occupation tolerance'
      end if
      if (band_first < 1 .or. band_last < 0) then
         error stop trim(context)//': response occupation state has an invalid band-window policy'
      end if
   end subroutine validate_response_occupation_fields

end module tddft_occupation_mod
