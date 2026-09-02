!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Measurement and work-decomposition helpers for linear-response TDDFT.
!>
!> The planner is intentionally independent of MPI handles.  It describes the
!> ownership of response work and leaves collectives to the caller, which keeps
!> the response backends usable in serial builds and with the existing MPI
!> compatibility layer.  In particular, native real-space work is split over
!> R/energy blocks while every worker retains the complete q batch.  This is the
!> important distinction from q-parallel duplication of chi0(R,omega).
module tddft_performance_mod
   use, intrinsic :: iso_fortran_env, only: int64
   use precision_mod, only: rp
   implicit none
   private

   integer, parameter, public :: tddft_profile_stage_count = 7
   character(len=32), parameter, public :: tddft_profile_stage_names(tddft_profile_stage_count) = [ &
      character(len=32) :: 'H_G_construction', &
      'GF_contraction', &
      'energy_integration', &
      'R_q_FT', &
      'Dyson', &
      'spectral_analysis', &
      'MPI_reduction']

   type, public :: tddft_work_range
      integer :: first = 1
      integer :: last = 0
      integer :: count = 0
   contains
      procedure :: contains => tddft_work_range_contains
   end type tddft_work_range

   !> Hierarchical response-work plan.  `q` is the compute range and `owner_q`
   !> is the range allowed to write q-labelled products.  They differ for the
   !> native real-space route: all ranks compute the q batch after their local
   !> R/energy reduction, while only the q owners write files.
   type, public :: tddft_mpi_plan
      character(len=32) :: backend = ''
      character(len=24) :: integration = ''
      character(len=96) :: strategy = ''
      integer :: rank = 0
      integer :: size = 1
      integer :: nq = 0
      integer :: nw = 0
      integer :: nk = 0
      integer :: nr = 0
      integer :: ne = 0
      type(tddft_work_range) :: q
      type(tddft_work_range) :: owner_q
      type(tddft_work_range) :: omega
      type(tddft_work_range) :: k
      type(tddft_work_range) :: r
      type(tddft_work_range) :: energy
      logical :: preserves_realspace_reuse = .false.
      logical :: requires_collective_reduction = .false.
      logical :: q_fourier_is_batched = .false.
      real(rp) :: q_duplication_factor = 1.0_rp
   contains
      procedure :: emit => emit_tddft_mpi_plan
   end type tddft_mpi_plan

   !> Low-overhead wall-clock profile.  Disabled profiles have no timer calls
   !> in the production code; callers can add measured child timings when a
   !> backend already owns a more detailed timer.
   type, public :: tddft_performance_profile
      logical :: enabled = .false.
      integer(int64) :: clock_rate = 1_int64
      integer(int64) :: wall_start = 0_int64
      integer(int64) :: stage_start = 0_int64
      integer :: active_stage = 0
      integer :: sample_count = 0
      character(len=32) :: backend = ''
      character(len=24) :: integration = ''
      real(rp) :: stage_seconds(tddft_profile_stage_count) = 0.0_rp
      integer(int64) :: bytes_hg = 0_int64
      integer(int64) :: bytes_gf = 0_int64
      integer(int64) :: bytes_response = 0_int64
      integer(int64) :: bytes_peak = 0_int64
      integer(int64) :: response_points = 0_int64
      integer(int64) :: green_evaluations = 0_int64
      real(rp) :: scientific_checksum = 0.0_rp
   contains
      procedure :: configure => configure_tddft_profile
      procedure :: reset => reset_tddft_profile
      procedure :: start => start_tddft_profile
      procedure :: stop => stop_tddft_profile
      procedure :: add_seconds => add_tddft_profile_seconds
      procedure :: set_memory => set_tddft_profile_memory
      procedure :: set_workload => set_tddft_profile_workload
      procedure :: set_checksum => set_tddft_profile_checksum
      procedure :: emit => emit_tddft_profile
   end type tddft_performance_profile

   public :: make_tddft_mpi_plan
   public :: split_tddft_work
   public :: tddft_profile_stage_index
   public :: tddft_checksum_complex

contains

   !> Choose a decomposition from measured reuse opportunities.  q-parallel is
   !> the default for reciprocal backends.  If q does not fill the worker set,
   !> frequency parallelism is the next level and k parallelism is a fallback.
   !> Native R-GF deliberately keeps all q points local and partitions R first,
   !> then energy nodes; the caller reduces the partial response before the
   !> single batched q transform.
   function make_tddft_mpi_plan(backend, integration, nq, nw, nk, nr, rank_in, size_in, &
      allow_realspace_reuse, ne) result(plan)
      character(len=*), intent(in) :: backend, integration
      integer, intent(in) :: nq, nw, nk, nr, rank_in, size_in
      logical, intent(in) :: allow_realspace_reuse
      integer, intent(in), optional :: ne
      type(tddft_mpi_plan) :: plan
      integer :: effective_size

      plan%backend = lower_ascii(trim(backend))
      plan%integration = lower_ascii(trim(integration))
      plan%rank = max(0, rank_in)
      effective_size = max(1, size_in)
      plan%size = effective_size
      plan%nq = max(0, nq)
      plan%nw = max(0, nw)
      plan%nk = max(0, nk)
      plan%nr = max(0, nr)
      plan%ne = plan%nw
      if (present(ne)) plan%ne = max(0, ne)

      call split_tddft_work(plan%rank, plan%size, plan%nq, plan%owner_q)
      plan%q = plan%owner_q
      call split_tddft_work(plan%rank, plan%size, plan%nw, plan%omega)
      call split_tddft_work(plan%rank, plan%size, plan%nk, plan%k)
      call split_tddft_work(plan%rank, plan%size, plan%nr, plan%r)
      call split_tddft_work(plan%rank, plan%size, plan%ne, plan%energy)

      if (plan%backend == 'realspace_gf' .and. allow_realspace_reuse) then
         plan%q%first = 1
         plan%q%last = plan%nq
         plan%q%count = plan%nq
         plan%preserves_realspace_reuse = .true.
         plan%requires_collective_reduction = plan%size > 1
         plan%q_fourier_is_batched = .true.
         plan%q_duplication_factor = real(plan%size, rp)
         if (plan%nr >= plan%size) then
            plan%strategy = 'R-blocks -> collective chi0(R,w) -> batched q FT'
         else
            plan%strategy = 'energy-nodes -> collective chi0(R,w) -> batched q FT'
         end if
      else if (plan%nq >= plan%size) then
         plan%strategy = 'q-outer / omega-inner'
         plan%q_fourier_is_batched = .false.
      else if (plan%nw >= plan%size) then
         plan%strategy = 'omega-outer / q-inner'
         plan%q%first = 1
         plan%q%last = plan%nq
         plan%q%count = plan%nq
         plan%requires_collective_reduction = plan%size > 1
      else if (plan%nk >= plan%size) then
         plan%strategy = 'k-block fallback with collective response reduction'
         plan%q%first = 1
         plan%q%last = plan%nq
         plan%q%count = plan%nq
         plan%requires_collective_reduction = plan%size > 1
      else
         plan%strategy = 'q-outer / oversubscribed worker set'
      end if
   end function make_tddft_mpi_plan

   pure subroutine split_tddft_work(worker_rank, worker_size, n_items, range)
      integer, intent(in) :: worker_rank, worker_size, n_items
      type(tddft_work_range), intent(out) :: range
      integer :: base, remainder

      range%first = 1
      range%last = 0
      range%count = 0
      if (worker_size <= 0 .or. worker_rank < 0 .or. worker_rank >= worker_size .or. n_items <= 0) return
      base = n_items/worker_size
      remainder = mod(n_items, worker_size)
      if (worker_rank < remainder) then
         range%count = base + 1
         range%first = worker_rank*range%count + 1
      else
         range%count = base
         range%first = worker_rank*base + remainder + 1
      end if
      range%last = range%first + range%count - 1
      if (range%count == 0) then
         range%first = 1
         range%last = 0
      end if
   end subroutine split_tddft_work

   pure logical function tddft_work_range_contains(this, item)
      class(tddft_work_range), intent(in) :: this
      integer, intent(in) :: item
      tddft_work_range_contains = this%count > 0 .and. item >= this%first .and. item <= this%last
   end function tddft_work_range_contains

   subroutine emit_tddft_mpi_plan(this, unit, label)
      class(tddft_mpi_plan), intent(in) :: this
      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: label
      integer :: output_unit
      character(len=64) :: output_label

      output_unit = 6
      if (present(unit)) output_unit = unit
      output_label = 'response'
      if (present(label)) output_label = trim(label)
      write(output_unit, '(a,1x,a,1x,a,a,1x,a,a,1x,a,i0,1x,a,i0,1x,a,a,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,l1,1x,a,l1,1x,a,es12.4)') &
         'TDDFT_PERF_PLAN', trim(output_label), 'backend=', trim(this%backend), 'integration=', trim(this%integration), &
         'rank=', this%rank, 'size=', this%size, 'strategy=', trim(this%strategy), 'q_first=', this%q%first, &
         'q_last=', this%q%last, 'q_count=', this%q%count, 'owner_q_first=', this%owner_q%first, &
         'owner_q_last=', this%owner_q%last, 'r_first=', this%r%first, 'r_last=', this%r%last, &
         'energy_first=', this%energy%first, 'energy_last=', this%energy%last, 'reuse=', this%preserves_realspace_reuse, &
         'collective=', this%requires_collective_reduction, 'q_duplication=', this%q_duplication_factor
   end subroutine emit_tddft_mpi_plan

   integer function tddft_profile_stage_index(name) result(index)
      character(len=*), intent(in) :: name
      integer :: i

      index = 0
      do i = 1, tddft_profile_stage_count
         if (trim(lower_ascii(name)) == trim(lower_ascii(tddft_profile_stage_names(i)))) then
            index = i
            return
         end if
      end do
   end function tddft_profile_stage_index

   subroutine configure_tddft_profile(this, enabled, backend, integration)
      class(tddft_performance_profile), intent(inout) :: this
      logical, intent(in) :: enabled
      character(len=*), intent(in), optional :: backend, integration

      call this%reset()
      this%enabled = enabled
      if (present(backend)) this%backend = trim(backend)
      if (present(integration)) this%integration = trim(integration)
      call system_clock(count_rate=this%clock_rate)
      if (this%clock_rate <= 0_int64) this%clock_rate = 1_int64
   end subroutine configure_tddft_profile

   subroutine reset_tddft_profile(this)
      class(tddft_performance_profile), intent(inout) :: this

      this%wall_start = 0_int64
      this%stage_start = 0_int64
      this%active_stage = 0
      this%sample_count = 0
      this%stage_seconds = 0.0_rp
      this%bytes_hg = 0_int64
      this%bytes_gf = 0_int64
      this%bytes_response = 0_int64
      this%bytes_peak = 0_int64
      this%response_points = 0_int64
      this%green_evaluations = 0_int64
      this%scientific_checksum = 0.0_rp
   end subroutine reset_tddft_profile

   subroutine start_tddft_profile(this, stage)
      class(tddft_performance_profile), intent(inout) :: this
      character(len=*), intent(in) :: stage
      integer :: index

      if (.not. this%enabled) return
      index = tddft_profile_stage_index(stage)
      if (index == 0) return
      if (this%active_stage /= 0) call this%stop(tddft_profile_stage_names(this%active_stage))
      this%active_stage = index
      call system_clock(this%stage_start)
      if (this%wall_start == 0_int64) this%wall_start = this%stage_start
   end subroutine start_tddft_profile

   subroutine stop_tddft_profile(this, stage)
      class(tddft_performance_profile), intent(inout) :: this
      character(len=*), intent(in) :: stage
      integer :: index
      integer(int64) :: now

      if (.not. this%enabled) return
      index = tddft_profile_stage_index(stage)
      if (index == 0 .or. this%active_stage /= index) return
      call system_clock(now)
      this%stage_seconds(index) = this%stage_seconds(index) + real(max(0_int64, now-this%stage_start), rp)/real(this%clock_rate, rp)
      this%active_stage = 0
   end subroutine stop_tddft_profile

   subroutine add_tddft_profile_seconds(this, stage, seconds)
      class(tddft_performance_profile), intent(inout) :: this
      character(len=*), intent(in) :: stage
      real(rp), intent(in) :: seconds
      integer :: index

      if (.not. this%enabled) return
      index = tddft_profile_stage_index(stage)
      if (index > 0) this%stage_seconds(index) = this%stage_seconds(index) + max(0.0_rp, seconds)
   end subroutine add_tddft_profile_seconds

   subroutine set_tddft_profile_memory(this, bytes_hg, bytes_gf, bytes_response, bytes_peak)
      class(tddft_performance_profile), intent(inout) :: this
      integer(int64), intent(in) :: bytes_hg, bytes_gf, bytes_response, bytes_peak

      if (.not. this%enabled) return
      this%bytes_hg = max(0_int64, bytes_hg)
      this%bytes_gf = max(0_int64, bytes_gf)
      this%bytes_response = max(0_int64, bytes_response)
      this%bytes_peak = max(0_int64, bytes_peak)
   end subroutine set_tddft_profile_memory

   subroutine set_tddft_profile_workload(this, response_points, green_evaluations, sample_count)
      class(tddft_performance_profile), intent(inout) :: this
      integer(int64), intent(in) :: response_points, green_evaluations
      integer, intent(in), optional :: sample_count

      if (.not. this%enabled) return
      this%response_points = max(0_int64, response_points)
      this%green_evaluations = max(0_int64, green_evaluations)
      if (present(sample_count)) this%sample_count = max(0, sample_count)
   end subroutine set_tddft_profile_workload

   subroutine set_tddft_profile_checksum(this, checksum)
      class(tddft_performance_profile), intent(inout) :: this
      real(rp), intent(in) :: checksum

      if (this%enabled) this%scientific_checksum = checksum
   end subroutine set_tddft_profile_checksum

   subroutine emit_tddft_profile(this, unit, label)
      class(tddft_performance_profile), intent(inout) :: this
      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: label
      integer :: output_unit
      character(len=64) :: output_label
      real(rp) :: total

      if (.not. this%enabled) return
      if (this%active_stage /= 0) call this%stop(tddft_profile_stage_names(this%active_stage))
      output_unit = 6
      if (present(unit)) output_unit = unit
      output_label = 'response'
      if (present(label)) output_label = trim(label)
      total = sum(this%stage_seconds)
      write(output_unit, '(a,1x,a,1x,a,a,1x,a,a,1x,a,i0,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,es12.4,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,es16.8)') &
         'TDDFT_PERF_PROFILE', trim(output_label), 'backend=', trim(this%backend), 'integration=', trim(this%integration), &
         'samples=', this%sample_count, 'H_G_s=', this%stage_seconds(1), 'GF_s=', this%stage_seconds(2), &
         'energy_s=', this%stage_seconds(3), 'R_q_FT_s=', this%stage_seconds(4), 'Dyson_s=', this%stage_seconds(5), &
         'spectral_s=', this%stage_seconds(6), 'MPI_s=', this%stage_seconds(7), 'total_s=', total, &
         'bytes_H_G=', this%bytes_hg, 'bytes_GF=', this%bytes_gf, 'bytes_response=', this%bytes_response, &
         'bytes_peak=', this%bytes_peak, 'response_points=', this%response_points, 'green_evaluations=', this%green_evaluations, &
         'checksum=', this%scientific_checksum
   end subroutine emit_tddft_profile

   pure real(rp) function tddft_checksum_complex(values) result(checksum)
      complex(rp), intent(in) :: values(..)

      checksum = 0.0_rp
      select rank (values)
      rank (1)
         checksum = checksum_rank1(values)
      rank (2)
         checksum = checksum_rank1(reshape(values, [size(values)]))
      rank (3)
         checksum = checksum_rank1(reshape(values, [size(values)]))
      rank (4)
         checksum = checksum_rank1(reshape(values, [size(values)]))
      rank default
         checksum = 0.0_rp
      end select
   end function tddft_checksum_complex

   pure real(rp) function checksum_rank1(values) result(checksum)
      complex(rp), intent(in) :: values(:)
      integer :: i

      checksum = 0.0_rp
      do i = 1, size(values)
         checksum = checksum + real(i, rp)*real(values(i), rp) + real(i+1, rp)*aimag(values(i))
      end do
   end function checksum_rank1

   pure function lower_ascii(input) result(output)
      character(len=*), intent(in) :: input
      character(len=len(input)) :: output
      integer :: i, code

      output = input
      do i = 1, len(input)
         code = iachar(output(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) output(i:i) = achar(code+32)
      end do
   end function lower_ascii

end module tddft_performance_mod
