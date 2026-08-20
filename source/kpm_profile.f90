!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: KPM profile
!
!> Small, transport-specific timing record for the KPM/Kubo-Bastin path.
!> The existing global timer remains the authoritative whole-program timer;
!> this record adds machine-readable intervals for the stages that matter when
!> comparing CPU and CUDA moment implementations.
!------------------------------------------------------------------------------

module kpm_profile_mod

   use iso_fortran_env, only: int64, output_unit
   use precision_mod, only: rp
   implicit none

   private

   integer, parameter, public :: kpm_stage_count = 10
   integer, parameter :: stage_operator = 1
   integer, parameter :: stage_trace_setup = 2
   integer, parameter :: stage_cheb_moments = 3
   integer, parameter :: stage_h2d = 4
   integer, parameter :: stage_d2h = 5
   integer, parameter :: stage_gamma = 6
   integer, parameter :: stage_mu_pack = 7
   integer, parameter :: stage_gamma_mu = 8
   integer, parameter :: stage_energy_integral = 9
   integer, parameter :: stage_transport_total = 10

   character(len=*), parameter :: stage_names(kpm_stage_count) = [ character(len=20) :: &
      'T_operator', 'T_trace_setup', 'T_cheb_moments', 'T_H2D', 'T_D2H', &
      'T_gamma', 'T_mu_pack', 'T_gamma_mu', 'T_energy_integral', 'T_transport_total' ]

   type, public :: kpm_profile
      private
      integer(int64) :: clock_rate = 0_int64
      integer(int64) :: stage_start(kpm_stage_count) = 0_int64
      real(rp) :: seconds(kpm_stage_count) = 0.0_rp
      logical :: active(kpm_stage_count) = .false.
      logical :: configured = .false.
      character(len=32) :: backend = 'unknown'
      character(len=48) :: precision = 'unknown'
      character(len=32) :: estimator = 'unknown'
      integer(int64) :: matrix_dimension = 0_int64
      integer(int64) :: nnz = 0_int64
      integer :: moments = 0
      integer :: lld = 0
      integer :: ntrace = 0
      integer(int64) :: bytes_h2d = 0_int64
      integer(int64) :: bytes_d2h = 0_int64
      integer(int64) :: bytes_gamma = 0_int64
      integer(int64) :: bytes_mu_pack = 0_int64
   contains
      procedure :: reset => kpm_profile_reset
      procedure :: configure => kpm_profile_configure
      procedure :: start => kpm_profile_start
      procedure :: stop => kpm_profile_stop
      procedure :: add_seconds => kpm_profile_add_seconds
      procedure :: add_bytes => kpm_profile_add_bytes
      procedure :: set_reconstruction_bytes => kpm_profile_set_reconstruction_bytes
      procedure :: emit => kpm_profile_emit
   end type kpm_profile

   type(kpm_profile), public, save :: g_kpm_profile

contains

   subroutine kpm_profile_reset(this)
      class(kpm_profile), intent(inout) :: this

      call system_clock(count_rate=this%clock_rate)
      this%stage_start = 0_int64
      this%seconds = 0.0_rp
      this%active = .false.
      this%configured = .false.
      this%backend = 'unknown'
      this%precision = 'unknown'
      this%estimator = 'unknown'
      this%matrix_dimension = 0_int64
      this%nnz = 0_int64
      this%moments = 0
      this%lld = 0
      this%ntrace = 0
      this%bytes_h2d = 0_int64
      this%bytes_d2h = 0_int64
      this%bytes_gamma = 0_int64
      this%bytes_mu_pack = 0_int64
   end subroutine kpm_profile_reset

   subroutine kpm_profile_configure(this, backend, precision, estimator, matrix_dimension, nnz, moments, lld, ntrace)
      class(kpm_profile), intent(inout) :: this
      character(len=*), intent(in) :: backend, precision, estimator
      integer(int64), intent(in) :: matrix_dimension, nnz
      integer, intent(in) :: moments, lld, ntrace

      if (this%clock_rate <= 0_int64) call system_clock(count_rate=this%clock_rate)
      this%configured = .true.
      this%backend = backend
      this%precision = precision
      this%estimator = estimator
      this%matrix_dimension = matrix_dimension
      this%nnz = nnz
      this%moments = moments
      this%lld = lld
      this%ntrace = ntrace
   end subroutine kpm_profile_configure

   integer function stage_index(name) result(index)
      character(len=*), intent(in) :: name
      integer :: i

      index = 0
      do i = 1, kpm_stage_count
         if (trim(stage_names(i)) == trim(name)) then
            index = i
            return
         end if
      end do
   end function stage_index

   subroutine kpm_profile_start(this, name)
      class(kpm_profile), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer :: index

      index = stage_index(name)
      if (index == 0) return
      if (this%clock_rate <= 0_int64) call system_clock(count_rate=this%clock_rate)
      call system_clock(this%stage_start(index))
      this%active(index) = .true.
   end subroutine kpm_profile_start

   subroutine kpm_profile_stop(this, name)
      class(kpm_profile), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer :: index
      integer(int64) :: now

      index = stage_index(name)
      if (index == 0 .or. .not. this%active(index)) return
      call system_clock(now)
      this%seconds(index) = this%seconds(index) + real(now - this%stage_start(index), rp) / &
         real(this%clock_rate, rp)
      this%active(index) = .false.
   end subroutine kpm_profile_stop

   subroutine kpm_profile_add_seconds(this, name, seconds)
      class(kpm_profile), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: seconds
      integer :: index

      index = stage_index(name)
      if (index > 0) this%seconds(index) = this%seconds(index) + max(0.0_rp, seconds)
   end subroutine kpm_profile_add_seconds

   subroutine kpm_profile_add_bytes(this, direction, bytes)
      class(kpm_profile), intent(inout) :: this
      character(len=*), intent(in) :: direction
      integer(int64), intent(in) :: bytes

      select case (trim(direction))
      case ('H2D')
         this%bytes_h2d = this%bytes_h2d + max(0_int64, bytes)
      case ('D2H')
         this%bytes_d2h = this%bytes_d2h + max(0_int64, bytes)
      end select
   end subroutine kpm_profile_add_bytes

   subroutine kpm_profile_set_reconstruction_bytes(this, gamma_bytes, mu_pack_bytes)
      class(kpm_profile), intent(inout) :: this
      integer(int64), intent(in) :: gamma_bytes, mu_pack_bytes

      this%bytes_gamma = max(0_int64, gamma_bytes)
      this%bytes_mu_pack = max(0_int64, mu_pack_bytes)
   end subroutine kpm_profile_set_reconstruction_bytes

   subroutine kpm_profile_emit(this)
      class(kpm_profile), intent(inout) :: this

      if (.not. this%configured) return
      do while (any(this%active))
         ! Close only malformed/incomplete intervals. Normal callers stop all
         ! stages explicitly; this keeps a profile useful on early exits.
         call close_active_stage(this)
      end do
      write (output_unit, *) 'KPM_PROFILE ', &
         'backend='//trim(this%backend)//' ', &
         'precision='//trim(this%precision)//' ', &
         'estimator='//trim(this%estimator)//' ', &
         'N= ', this%matrix_dimension, 'nnz= ', this%nnz, &
         'M= ', this%moments, 'lld= ', this%lld, 'Ntrace= ', this%ntrace, &
         'bytes_h2d= ', this%bytes_h2d, 'bytes_d2h= ', this%bytes_d2h, &
         'bytes_gamma= ', this%bytes_gamma, 'bytes_mu_pack= ', this%bytes_mu_pack, &
         'T_operator= ', this%seconds(stage_operator), &
         'T_trace_setup= ', this%seconds(stage_trace_setup), &
         'T_cheb_moments= ', this%seconds(stage_cheb_moments), &
         'T_H2D= ', this%seconds(stage_h2d), 'T_D2H= ', this%seconds(stage_d2h), &
         'T_gamma= ', this%seconds(stage_gamma), 'T_mu_pack= ', this%seconds(stage_mu_pack), &
         'T_gamma_mu= ', this%seconds(stage_gamma_mu), &
         'T_energy_integral= ', this%seconds(stage_energy_integral), &
         'T_transport_total= ', this%seconds(stage_transport_total)
   end subroutine kpm_profile_emit

   subroutine close_active_stage(this)
      class(kpm_profile), intent(inout) :: this
      integer :: i

      do i = 1, kpm_stage_count
         if (this%active(i)) then
            call kpm_profile_stop(this, stage_names(i))
            return
         end if
      end do
   end subroutine close_active_stage

end module kpm_profile_mod
