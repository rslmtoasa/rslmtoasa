!------------------------------------------------------------------------------
! SCF-B0C benchmark profile collector.
!
! This module is deliberately a small observation layer.  It owns no SCF
! state and does not alter any numerical operation; callers mark the existing
! production stage boundaries and the collector emits machine-readable rows.
!------------------------------------------------------------------------------
module scf_benchmark_profile_mod

   use, intrinsic :: iso_fortran_env, only: int64
   use precision_mod, only: rp
   implicit none

   private

   integer, parameter, public :: scf_profile_stage_count = 19
   character(len=32), parameter, public :: scf_profile_stage_names(scf_profile_stage_count) = [ &
      character(len=32) :: 'P_hamiltonian_prepare', &
      'P_hk_assembly', &
      'P_eigensolver', &
      'P_eigenpair_transfer', &
      'P_occupations_fermi', &
      'P_density_build', &
      'P_charge_spin_accumulate', &
      'P_potential_update', &
      'P_mixing', &
      'P_scf_io', &
      'P_scf_misc', &
      'P_rs_hamiltonian_prepare', &
      'P_rs_solver_kernel', &
      'P_rs_green_function', &
      'P_rs_spectral_reconstruct', &
      'P_rs_energy_integration', &
      'P_rs_fermi', &
      'P_rs_density_build', &
      'P_rs_charge_spin_accumulate']

   type, public :: scf_benchmark_profile
      logical :: enabled = .false.
      integer(int64) :: clock_rate = 1_int64
      integer(int64) :: iteration_start = 0_int64
      integer(int64) :: stage_start = 0_int64
      integer :: active_stage = 0
      integer :: iteration = 0
      real(rp) :: stage_seconds(scf_profile_stage_count) = 0.0_rp
      real(rp) :: detail_h2d_seconds = 0.0_rp
      real(rp) :: detail_solver_seconds = 0.0_rp
      real(rp) :: detail_d2h_seconds = 0.0_rp
      real(rp) :: detail_total_seconds = 0.0_rp
      real(rp) :: detail_rs_h2d_seconds = 0.0_rp
      real(rp) :: detail_rs_kernel_seconds = 0.0_rp
      real(rp) :: detail_rs_d2h_seconds = 0.0_rp
      real(rp) :: detail_rs_sync_seconds = 0.0_rp
      real(rp) :: detail_rs_setup_seconds = 0.0_rp
      real(rp) :: previous_h2d_seconds = 0.0_rp
      real(rp) :: previous_solver_seconds = 0.0_rp
      real(rp) :: previous_d2h_seconds = 0.0_rp
      real(rp) :: previous_total_seconds = 0.0_rp
      real(rp) :: last_total_seconds = 0.0_rp
      real(rp) :: last_closure_error = 0.0_rp
      integer :: completed_iterations = 0
   contains
      procedure :: configure => scf_profile_configure
      procedure :: configure_from_environment => scf_profile_configure_from_environment
      procedure :: start_iteration => scf_profile_start_iteration
      procedure :: start_stage => scf_profile_start_stage
      procedure :: stop_stage => scf_profile_stop_stage
      procedure :: add_stage => scf_profile_add_stage
      procedure :: set_reciprocal_details => scf_profile_set_reciprocal_details
      procedure :: finish_iteration => scf_profile_finish_iteration
      procedure :: reset => scf_profile_reset
   end type scf_benchmark_profile

   type(scf_benchmark_profile), save, public :: g_scf_benchmark_profile

   public :: scf_profile_stage_index

contains

   subroutine scf_profile_configure(this, enabled)
      class(scf_benchmark_profile), intent(inout) :: this
      logical, intent(in) :: enabled

      call this%reset()
      this%enabled = enabled
      call system_clock(count_rate=this%clock_rate)
      if (this%clock_rate <= 0_int64) this%clock_rate = 1_int64
   end subroutine scf_profile_configure

   subroutine scf_profile_configure_from_environment(this)
      class(scf_benchmark_profile), intent(inout) :: this
      character(len=32) :: value
      integer :: status

      call get_environment_variable('RSLMTO_SCF_B0C_PROFILE', value, status=status)
      call this%configure(status == 0 .and. (trim(value) == '1' .or. trim(value) == 'true' .or. trim(value) == 'TRUE'))
   end subroutine scf_profile_configure_from_environment

   subroutine scf_profile_reset(this)
      class(scf_benchmark_profile), intent(inout) :: this

      this%iteration_start = 0_int64
      this%stage_start = 0_int64
      this%active_stage = 0
      this%iteration = 0
      this%stage_seconds = 0.0_rp
      this%detail_h2d_seconds = 0.0_rp
      this%detail_solver_seconds = 0.0_rp
      this%detail_d2h_seconds = 0.0_rp
      this%detail_total_seconds = 0.0_rp
      this%detail_rs_h2d_seconds = 0.0_rp
      this%detail_rs_kernel_seconds = 0.0_rp
      this%detail_rs_d2h_seconds = 0.0_rp
      this%detail_rs_sync_seconds = 0.0_rp
      this%detail_rs_setup_seconds = 0.0_rp
      this%previous_h2d_seconds = 0.0_rp
      this%previous_solver_seconds = 0.0_rp
      this%previous_d2h_seconds = 0.0_rp
      this%previous_total_seconds = 0.0_rp
      this%last_total_seconds = 0.0_rp
      this%last_closure_error = 0.0_rp
      this%completed_iterations = 0
   end subroutine scf_profile_reset

   integer function scf_profile_stage_index(name) result(index)
      character(len=*), intent(in) :: name
      integer :: i

      index = 0
      do i = 1, scf_profile_stage_count
         if (trim(name) == trim(scf_profile_stage_names(i))) then
            index = i
            return
         end if
      end do
   end function scf_profile_stage_index

   subroutine scf_profile_start_iteration(this, iteration)
      class(scf_benchmark_profile), intent(inout) :: this
      integer, intent(in) :: iteration

      if (.not. this%enabled) return
      if (this%active_stage /= 0) this%active_stage = 0
      this%iteration = iteration
      this%stage_seconds = 0.0_rp
      this%detail_h2d_seconds = 0.0_rp
      this%detail_solver_seconds = 0.0_rp
      this%detail_d2h_seconds = 0.0_rp
      this%detail_total_seconds = 0.0_rp
      this%detail_rs_h2d_seconds = 0.0_rp
      this%detail_rs_kernel_seconds = 0.0_rp
      this%detail_rs_d2h_seconds = 0.0_rp
      this%detail_rs_sync_seconds = 0.0_rp
      this%detail_rs_setup_seconds = 0.0_rp
      call system_clock(this%iteration_start)
   end subroutine scf_profile_start_iteration

   subroutine scf_profile_start_stage(this, name)
      class(scf_benchmark_profile), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer :: index

      if (.not. this%enabled) return
      index = scf_profile_stage_index(name)
      if (index == 0) return
      if (this%active_stage /= 0) then
         call this%stop_stage(scf_profile_stage_names(this%active_stage))
      end if
      this%active_stage = index
      call system_clock(this%stage_start)
   end subroutine scf_profile_start_stage

   subroutine scf_profile_stop_stage(this, name)
      class(scf_benchmark_profile), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer(int64) :: now
      integer :: index

      if (.not. this%enabled) return
      index = scf_profile_stage_index(name)
      if (index == 0 .or. this%active_stage == 0) return
      call system_clock(now)
      if (this%active_stage == index) then
         this%stage_seconds(index) = this%stage_seconds(index) + &
            real(max(0_int64, now - this%stage_start), rp) / real(this%clock_rate, rp)
         if (trim(name) == 'P_rs_solver_kernel') then
            ! The production RS entry point is the timer-equivalent kernel
            ! boundary for CPU and CUDA.  CUDA transfer/setup counters are
            ! optional nested diagnostics and are not fabricated here.
            this%detail_rs_kernel_seconds = this%detail_rs_kernel_seconds + &
               real(max(0_int64, now - this%stage_start), rp) / real(this%clock_rate, rp)
         end if
         this%active_stage = 0
      end if
   end subroutine scf_profile_stop_stage

   subroutine scf_profile_add_stage(this, name, seconds)
      class(scf_benchmark_profile), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: seconds
      integer :: index

      if (.not. this%enabled) return
      index = scf_profile_stage_index(name)
      if (index > 0) this%stage_seconds(index) = this%stage_seconds(index) + max(0.0_rp, seconds)
   end subroutine scf_profile_add_stage

   subroutine scf_profile_set_reciprocal_details(this, h2d, solver, d2h, total)
      class(scf_benchmark_profile), intent(inout) :: this
      real(rp), intent(in) :: h2d, solver, d2h, total
      real(rp) :: delta_h2d, delta_solver, delta_d2h, delta_total, transfer
      integer :: eigensolver_index, transfer_index

      if (.not. this%enabled) return
      delta_h2d = max(0.0_rp, h2d - this%previous_h2d_seconds)
      delta_solver = max(0.0_rp, solver - this%previous_solver_seconds)
      delta_d2h = max(0.0_rp, d2h - this%previous_d2h_seconds)
      delta_total = max(0.0_rp, total - this%previous_total_seconds)
      this%previous_h2d_seconds = h2d
      this%previous_solver_seconds = solver
      this%previous_d2h_seconds = d2h
      this%previous_total_seconds = total
      this%detail_h2d_seconds = delta_h2d
      this%detail_solver_seconds = delta_solver
      this%detail_d2h_seconds = delta_d2h
      this%detail_total_seconds = delta_total

      ! The outer diagonalization measurement is authoritative for exclusive
      ! closure.  Device transfer detail is nested beneath that interval, so
      ! split only the transfer portion and retain any launch/synchronization
      ! remainder in P_eigensolver.
      transfer = min(this%stage_seconds(3), delta_h2d + delta_d2h)
      eigensolver_index = scf_profile_stage_index('P_eigensolver')
      transfer_index = scf_profile_stage_index('P_eigenpair_transfer')
      this%stage_seconds(eigensolver_index) = max(0.0_rp, this%stage_seconds(eigensolver_index) - transfer)
      this%stage_seconds(transfer_index) = this%stage_seconds(transfer_index) + transfer
   end subroutine scf_profile_set_reciprocal_details

   subroutine scf_profile_finish_iteration(this, residual, total_energy, fermi_energy, magnetic_moment)
      class(scf_benchmark_profile), intent(inout) :: this
      real(rp), intent(in) :: residual, total_energy, fermi_energy, magnetic_moment
      integer(int64) :: now
      real(rp) :: total, exclusive_sum, closure_error
      character(len=64) :: fields(34)
      integer :: i, misc_index

      if (.not. this%enabled) return
      if (this%active_stage /= 0) call this%stop_stage(scf_profile_stage_names(this%active_stage))
      call system_clock(now)
      total = real(max(0_int64, now - this%iteration_start), rp) / real(this%clock_rate, rp)
      exclusive_sum = sum(this%stage_seconds)
      misc_index = scf_profile_stage_index('P_scf_misc')
      if (total > exclusive_sum) this%stage_seconds(misc_index) = this%stage_seconds(misc_index) + total - exclusive_sum
      exclusive_sum = sum(this%stage_seconds)
      closure_error = abs(exclusive_sum - total) / max(total, 1.0e-12_rp)
      this%last_total_seconds = total
      this%last_closure_error = closure_error
      this%completed_iterations = this%completed_iterations + 1

      do i = 1, scf_profile_stage_count
         fields(i) = profile_field(trim(scf_profile_stage_names(i)), this%stage_seconds(i))
      end do
      fields(20) = profile_field('P_scf_iteration_total', total)
      fields(21) = profile_field('profile_closure_error', closure_error)
      fields(22) = profile_field('residual', residual)
      fields(23) = profile_field('total_energy', total_energy)
      fields(24) = profile_field('fermi_energy', fermi_energy)
      fields(25) = profile_field('magnetic_moment', magnetic_moment)
      fields(26) = profile_field('T_H2D', this%detail_h2d_seconds)
      fields(27) = profile_field('T_solver', this%detail_solver_seconds)
      fields(28) = profile_field('T_D2H', this%detail_d2h_seconds)
      fields(29) = profile_field('T_total_steady', this%detail_total_seconds)
      fields(30) = profile_field('T_rs_H2D', this%detail_rs_h2d_seconds)
      fields(31) = profile_field('T_rs_kernel', this%detail_rs_kernel_seconds)
      fields(32) = profile_field('T_rs_D2H', this%detail_rs_d2h_seconds)
      fields(33) = profile_field('T_rs_sync', this%detail_rs_sync_seconds)
      fields(34) = profile_field('T_rs_setup', this%detail_rs_setup_seconds)
      write (*, '(A,1X,I0,1X,*(A,1X))') 'SCF_B0C_ITER', this%iteration, fields
   end subroutine scf_profile_finish_iteration

   function profile_field(name, value) result(token)
      character(len=*), intent(in) :: name
      real(rp), intent(in) :: value
      character(len=64) :: token
      character(len=32) :: number

      write(number, '(ES16.8)') value
      token = trim(name)//'='//trim(adjustl(number))
   end function profile_field

end module scf_benchmark_profile_mod
