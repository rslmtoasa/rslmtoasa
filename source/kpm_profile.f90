!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: KPM profile
!
!> Small, transport-specific timing record for the KPM/Kubo-Bastin path.
!>
!> The profile deliberately has two fixed namespaces. ``P_*`` entries are
!> exclusive transport phases and ``D_*`` entries are implementation detail
!> timers nested inside one of those phases. The old ``T_*`` fields are kept
!> as deprecated compatibility aliases in the emitted record; they must not be
!> added to obtain a transport total.
!------------------------------------------------------------------------------

module kpm_profile_mod

   use iso_fortran_env, only: int64, output_unit
   use precision_mod, only: rp
   implicit none

   private

   integer, parameter, public :: kpm_phase_count = 8
   integer, parameter, public :: kpm_detail_count = 8

   integer, parameter :: phase_operator = 1
   integer, parameter :: phase_trace_setup = 2
   integer, parameter :: phase_moments_total = 3
   integer, parameter :: phase_gamma = 4
   integer, parameter :: phase_reconstruction_total = 5
   integer, parameter :: phase_energy_integration = 6
   integer, parameter :: phase_output_io = 7
   integer, parameter :: phase_other = 8

   integer, parameter :: detail_moment_h2d = 1
   integer, parameter :: detail_moment_gpu_kernel = 2
   integer, parameter :: detail_moment_d2h = 3
   integer, parameter :: detail_moment_conversion = 4
   integer, parameter :: detail_mu_pack = 5
   integer, parameter :: detail_reconstruction_blas = 6
   integer, parameter :: detail_gamma_basis = 7
   integer, parameter :: detail_gamma_fill = 8

   character(len=*), parameter :: phase_names(kpm_phase_count) = [ character(len=24) :: &
      'P_operator', 'P_trace_setup', 'P_moments_total', 'P_gamma', &
      'P_reconstruction_total', 'P_energy_integration', 'P_output_io', 'P_other' ]

   character(len=*), parameter :: detail_names(kpm_detail_count) = [ character(len=28) :: &
      'D_moment_H2D', 'D_moment_GPU_kernel', 'D_moment_D2H', 'D_conversion', &
      'D_mu_pack', 'D_reconstruction_BLAS', 'D_gamma_basis', 'D_gamma_fill' ]

   type, public :: kpm_profile
      private
      integer(int64) :: clock_rate = 0_int64
      integer(int64) :: phase_start(kpm_phase_count) = 0_int64
      integer(int64) :: detail_start(kpm_detail_count) = 0_int64
      integer(int64) :: transport_total_start = 0_int64
      real(rp) :: phase_seconds(kpm_phase_count) = 0.0_rp
      real(rp) :: detail_seconds(kpm_detail_count) = 0.0_rp
      real(rp) :: transport_total_seconds = 0.0_rp
      logical :: phase_active(kpm_phase_count) = .false.
      logical :: detail_active(kpm_detail_count) = .false.
      logical :: transport_total_active = .false.
      logical :: configured = .false.
      character(len=32) :: backend = 'unknown'
      character(len=32) :: moment_backend = 'unknown'
      character(len=16) :: moment_precision = 'unknown'
      character(len=32) :: reconstruction_backend = 'unknown'
      character(len=16) :: reconstruction_precision = 'unknown'
      character(len=48) :: precision = 'unknown'
      character(len=32) :: estimator = 'unknown'
      character(len=32) :: omp_threads = 'unset'
      character(len=32) :: blas_threads = 'unset'
      integer(int64) :: matrix_dimension = 0_int64
      integer(int64) :: nnz = 0_int64
      integer :: moments = 0
      integer :: lld = 0
      integer :: ntrace = 0
      integer(int64) :: bytes_h2d = 0_int64
      integer(int64) :: bytes_d2h = 0_int64
      integer(int64) :: bytes_gamma = 0_int64
      integer(int64) :: bytes_mu_pack = 0_int64
      real(rp) :: closure_error = 0.0_rp
      real(rp) :: child_error = 0.0_rp
      character(len=8) :: status = 'UNKNOWN'
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
      this%phase_start = 0_int64
      this%detail_start = 0_int64
      this%transport_total_start = 0_int64
      this%phase_seconds = 0.0_rp
      this%detail_seconds = 0.0_rp
      this%transport_total_seconds = 0.0_rp
      this%phase_active = .false.
      this%detail_active = .false.
      this%transport_total_active = .false.
      this%configured = .false.
      this%backend = 'unknown'
      this%moment_backend = 'unknown'
      this%moment_precision = 'unknown'
      this%reconstruction_backend = 'unknown'
      this%reconstruction_precision = 'unknown'
      this%precision = 'unknown'
      this%estimator = 'unknown'
      this%omp_threads = 'unset'
      this%blas_threads = 'unset'
      this%matrix_dimension = 0_int64
      this%nnz = 0_int64
      this%moments = 0
      this%lld = 0
      this%ntrace = 0
      this%bytes_h2d = 0_int64
      this%bytes_d2h = 0_int64
      this%bytes_gamma = 0_int64
      this%bytes_mu_pack = 0_int64
      this%closure_error = 0.0_rp
      this%child_error = 0.0_rp
      this%status = 'UNKNOWN'
   end subroutine kpm_profile_reset

   subroutine kpm_profile_configure(this, moment_backend, moment_precision, reconstruction_backend, &
                                    reconstruction_precision, estimator, matrix_dimension, nnz, moments, lld, ntrace)
      class(kpm_profile), intent(inout) :: this
      character(len=*), intent(in) :: moment_backend, moment_precision
      character(len=*), intent(in) :: reconstruction_backend, reconstruction_precision
      character(len=*), intent(in) :: estimator
      integer(int64), intent(in) :: matrix_dimension, nnz
      integer, intent(in) :: moments, lld, ntrace

      if (this%clock_rate <= 0_int64) call system_clock(count_rate=this%clock_rate)
      this%configured = .true.
      this%moment_backend = moment_backend
      this%moment_precision = moment_precision
      this%reconstruction_backend = reconstruction_backend
      this%reconstruction_precision = reconstruction_precision
      this%estimator = estimator
      if (trim(moment_backend) == 'cuda') then
         this%backend = 'cuda'
      else
         this%backend = 'cpu'
      end if
      this%precision = trim(moment_precision)//'_'//trim(reconstruction_precision)
      this%omp_threads = environment_value('OMP_NUM_THREADS')
      this%blas_threads = environment_value('BLAS_NUM_THREADS')
      if (trim(this%blas_threads) == 'unset') this%blas_threads = environment_value('MKL_NUM_THREADS')
      if (trim(this%blas_threads) == 'unset') this%blas_threads = environment_value('OPENBLAS_NUM_THREADS')
      this%matrix_dimension = matrix_dimension
      this%nnz = nnz
      this%moments = moments
      this%lld = lld
      this%ntrace = ntrace
   end subroutine kpm_profile_configure

   integer function phase_index(name) result(index)
      character(len=*), intent(in) :: name
      integer :: i

      index = 0
      do i = 1, kpm_phase_count
         if (trim(phase_names(i)) == trim(name)) then
            index = i
            return
         end if
      end do
      select case (trim(name))
      case ('T_operator')
         index = phase_operator
      case ('T_trace_setup')
         index = phase_trace_setup
      case ('T_cheb_moments')
         index = phase_moments_total
      case ('T_gamma')
         index = phase_gamma
      case ('T_energy_integral')
         index = phase_energy_integration
      end select
   end function phase_index

   integer function detail_index(name) result(index)
      character(len=*), intent(in) :: name
      integer :: i

      index = 0
      do i = 1, kpm_detail_count
         if (trim(detail_names(i)) == trim(name)) then
            index = i
            return
         end if
      end do
      select case (trim(name))
      case ('T_H2D')
         index = detail_moment_h2d
      case ('T_D2H')
         index = detail_moment_d2h
      case ('T_mu_pack')
         index = detail_mu_pack
      case ('T_gamma_mu')
         index = detail_reconstruction_blas
      end select
   end function detail_index

   subroutine kpm_profile_start(this, name)
      class(kpm_profile), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer :: index

      if (trim(name) == 'T_transport_total') then
         if (this%clock_rate <= 0_int64) call system_clock(count_rate=this%clock_rate)
         call system_clock(this%transport_total_start)
         this%transport_total_active = .true.
         return
      end if
      index = phase_index(name)
      if (index > 0 .and. index /= phase_other) then
         if (this%clock_rate <= 0_int64) call system_clock(count_rate=this%clock_rate)
         call system_clock(this%phase_start(index))
         this%phase_active(index) = .true.
         return
      end if
      index = detail_index(name)
      if (index > 0) then
         if (this%clock_rate <= 0_int64) call system_clock(count_rate=this%clock_rate)
         call system_clock(this%detail_start(index))
         this%detail_active(index) = .true.
      end if
   end subroutine kpm_profile_start

   subroutine kpm_profile_stop(this, name)
      class(kpm_profile), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer :: index
      integer(int64) :: now

      if (trim(name) == 'T_transport_total' .and. this%transport_total_active) then
         call system_clock(now)
         this%transport_total_seconds = this%transport_total_seconds + &
            real(now - this%transport_total_start, rp) / real(this%clock_rate, rp)
         this%transport_total_active = .false.
         return
      end if
      index = phase_index(name)
      if (index > 0 .and. index /= phase_other .and. this%phase_active(index)) then
         call system_clock(now)
         this%phase_seconds(index) = this%phase_seconds(index) + real(now - this%phase_start(index), rp) / &
            real(this%clock_rate, rp)
         this%phase_active(index) = .false.
         return
      end if
      index = detail_index(name)
      if (index > 0 .and. this%detail_active(index)) then
         call system_clock(now)
         this%detail_seconds(index) = this%detail_seconds(index) + real(now - this%detail_start(index), rp) / &
            real(this%clock_rate, rp)
         this%detail_active(index) = .false.
      end if
   end subroutine kpm_profile_stop

   subroutine kpm_profile_add_seconds(this, name, seconds)
      class(kpm_profile), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: seconds
      integer :: index

      index = detail_index(name)
      if (index > 0) then
         this%detail_seconds(index) = this%detail_seconds(index) + max(0.0_rp, seconds)
         return
      end if
      index = phase_index(name)
      if (index > 0 .and. index /= phase_other) then
         this%phase_seconds(index) = this%phase_seconds(index) + max(0.0_rp, seconds)
      end if
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
      real(rp) :: exclusive_sum, raw_other, total, tolerance
      real(rp) :: moment_children, reconstruction_children, gamma_children

      if (.not. this%configured) return
      do while (any(this%phase_active) .or. any(this%detail_active) .or. this%transport_total_active)
         call close_active_timer(this)
      end do

      total = this%transport_total_seconds
      raw_other = total - sum(this%phase_seconds(1:phase_output_io))
      this%phase_seconds(phase_other) = raw_other
      exclusive_sum = sum(this%phase_seconds)
      tolerance = max(0.05_rp * max(total, 1.0e-12_rp), 0.01_rp)
      this%closure_error = abs(exclusive_sum - total) / max(total, 1.0e-12_rp)

      moment_children = sum(this%detail_seconds(detail_moment_h2d:detail_moment_conversion))
      reconstruction_children = sum(this%detail_seconds(detail_mu_pack:detail_reconstruction_blas))
      gamma_children = sum(this%detail_seconds(detail_gamma_basis:detail_gamma_fill))
      this%child_error = max(0.0_rp, max(moment_children - this%phase_seconds(phase_moments_total), &
         max(reconstruction_children - this%phase_seconds(phase_reconstruction_total), &
             gamma_children - this%phase_seconds(phase_gamma))))
      if (this%closure_error <= 0.05_rp .and. this%phase_seconds(phase_other) >= -tolerance .and. &
          moment_children <= this%phase_seconds(phase_moments_total) + tolerance .and. &
          reconstruction_children <= this%phase_seconds(phase_reconstruction_total) + tolerance .and. &
          gamma_children <= this%phase_seconds(phase_gamma) + tolerance) then
         this%status = 'PASS'
      else
         this%status = 'FAIL'
      end if

      write (output_unit, *) 'KPM_PROFILE ', &
         'backend='//trim(this%backend)//' ', &
         'moment_backend='//trim(this%moment_backend)//' ', &
         'moment_precision='//trim(this%moment_precision)//' ', &
         'reconstruction_backend='//trim(this%reconstruction_backend)//' ', &
         'reconstruction_precision='//trim(this%reconstruction_precision)//' ', &
         'precision='//trim(this%precision)//' ', &
         'estimator='//trim(this%estimator)//' ', &
         'N= ', this%matrix_dimension, 'nnz= ', this%nnz, &
         'M= ', this%moments, 'lld= ', this%lld, 'Ntrace= ', this%ntrace, &
         'OMP_NUM_THREADS='//trim(this%omp_threads)//' ', &
         'BLAS_NUM_THREADS='//trim(this%blas_threads)//' ', &
         'omp_threads='//trim(this%omp_threads)//' ', &
         'blas_threads='//trim(this%blas_threads)//' ', &
         'clock_source=system_clock_wall ', &
         'bytes_h2d= ', this%bytes_h2d, 'bytes_d2h= ', this%bytes_d2h, &
         'bytes_gamma= ', this%bytes_gamma, 'bytes_mu_pack= ', this%bytes_mu_pack, &
         'P_operator= ', this%phase_seconds(phase_operator), &
         'P_trace_setup= ', this%phase_seconds(phase_trace_setup), &
         'P_moments_total= ', this%phase_seconds(phase_moments_total), &
         'P_gamma= ', this%phase_seconds(phase_gamma), &
         'P_reconstruction_total= ', this%phase_seconds(phase_reconstruction_total), &
         'P_energy_integration= ', this%phase_seconds(phase_energy_integration), &
         'P_output_io= ', this%phase_seconds(phase_output_io), &
         'P_other= ', this%phase_seconds(phase_other), &
         'D_moment_H2D= ', this%detail_seconds(detail_moment_h2d), &
         'D_moment_GPU_kernel= ', this%detail_seconds(detail_moment_gpu_kernel), &
         'D_moment_D2H= ', this%detail_seconds(detail_moment_d2h), &
         'D_conversion= ', this%detail_seconds(detail_moment_conversion), &
         'D_mu_pack= ', this%detail_seconds(detail_mu_pack), &
         'D_reconstruction_BLAS= ', this%detail_seconds(detail_reconstruction_blas), &
         'D_gamma_basis= ', this%detail_seconds(detail_gamma_basis), &
         'D_gamma_fill= ', this%detail_seconds(detail_gamma_fill), &
         'T_operator= ', this%phase_seconds(phase_operator), &
         'T_trace_setup= ', this%phase_seconds(phase_trace_setup), &
         'T_cheb_moments= ', this%phase_seconds(phase_moments_total), &
         'T_H2D= ', this%detail_seconds(detail_moment_h2d), &
         'T_D2H= ', this%detail_seconds(detail_moment_d2h), &
         'T_gamma= ', this%phase_seconds(phase_gamma), &
         'T_mu_pack= ', this%detail_seconds(detail_mu_pack), &
         'T_gamma_mu= ', this%detail_seconds(detail_reconstruction_blas), &
         'T_energy_integral= ', this%phase_seconds(phase_energy_integration), &
         'T_transport_total= ', total, &
         'profile_closure_error= ', this%closure_error, &
         'profile_child_error= ', this%child_error, &
         'PROFILE_STATUS='//trim(this%status)
   end subroutine kpm_profile_emit

   subroutine close_active_timer(this)
      class(kpm_profile), intent(inout) :: this
      integer :: i

      if (this%transport_total_active) then
         call kpm_profile_stop(this, 'T_transport_total')
         return
      end if
      do i = 1, kpm_phase_count
         if (this%phase_active(i)) then
            call kpm_profile_stop(this, phase_names(i))
            return
         end if
      end do
      do i = 1, kpm_detail_count
         if (this%detail_active(i)) then
            call kpm_profile_stop(this, detail_names(i))
            return
         end if
      end do
   end subroutine close_active_timer

   character(len=32) function environment_value(name) result(value)
      character(len=*), intent(in) :: name
      character(len=256) :: buffer
      integer :: length, status

      value = 'unset'
      buffer = ''
      call get_environment_variable(trim(name), buffer, length=length, status=status)
      if (status == 0 .and. length > 0) value = buffer(1:min(length, len(value)))
   end function environment_value

end module kpm_profile_mod
