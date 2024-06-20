!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Timer
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module to mesure time between code
!------------------------------------------------------------------------------

module timer_mod
   use string_mod, only: sl, clean_str, split, endswith, fmt
   use precision_mod, only: rp, longint
   use logger_mod, only: g_logger
   use report_mod
   implicit none

   ! Private by default
   private

   type :: simple_timer
      character(len=sl) :: label
      integer(longint), private :: t_start, t_stop
      !real, private :: t_start, t_stop
   contains
      procedure :: start => simple_timer_start
      procedure :: stop => simple_timer_stop
      procedure :: has_stoped => simple_timer_has_stoped
      procedure :: delta => simpletimer_delta
      procedure :: report => simple_timer_report
   end type simple_timer

   type :: timer
      type(simple_timer), allocatable, dimension(:), private :: started_times
      type(simple_timer), allocatable, dimension(:), private :: finished_times
      type(simple_timer), private :: t_global
      logical :: auto_couple
      character(len=2*sl) :: coupled_label
      integer(kind=longint):: clock_rate
   contains
      procedure :: start => timer_start
      procedure :: stop => timer_stop
      procedure :: print_report => timer_report
      final :: timer_destructor
   end type timer

   interface timer
      procedure :: timer_constructor
   end interface timer

   ! Paramenters
   integer, dimension(8) :: values
   real(8), dimension(8), parameter :: seconds_of_values = (/3.154e+7, 2.628e+6, 1440.0, 0.0, 3600.0, 60.0, 1.0, 1e-3/)

   ! Global
   type(timer) :: g_timer

   ! Public
   public :: timer, g_timer

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> timer_constructor
   !
   !> @return type(timer)
   !---------------------------------------------------------------------------
   function timer_constructor(auto_couple) result(obj)
      type(timer) :: obj
      logical, optional :: auto_couple
      call obj%t_global%start('global')
      call system_clock(count_rate=obj%clock_rate)
      allocate (obj%started_times(0))
      allocate (obj%finished_times(0))
      !obj%auto_couple = .not.present(auto_couple) .or. auto_couple
      if (.not. present(auto_couple)) then
         obj%auto_couple = .true.
      else if (auto_couple) then
         obj%auto_couple = .true.
      else
         obj%auto_couple = .false.
      end if
      if (obj%auto_couple) obj%coupled_label = ''
   end function timer_constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine timer_destructor(this)
      type(timer) :: this
      if (allocated(this%started_times)) deallocate (this%started_times)
      if (allocated(this%finished_times)) deallocate (this%finished_times)
   end subroutine timer_destructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Measures the start time
   !---------------------------------------------------------------------------
   subroutine simple_timer_start(this, label)
      class(simple_timer) :: this
      character(len=*), intent(in) :: label
      this%label = label
      call system_clock(this%t_start)
      !call cpu_time(this%t_start)
      this%t_stop = this%t_start
   end subroutine simple_timer_start

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Measures the final time
   !---------------------------------------------------------------------------
   subroutine simple_timer_stop(this)
      class(simple_timer) :: this
      call system_clock(this%t_stop)
      !call cpu_time(this%t_stop)
   end subroutine simple_timer_stop

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Reports the time spent with label
   !---------------------------------------------------------------------------
   subroutine simple_timer_report(this)
      class(simple_timer) :: this
      call g_logger%log(this%label, this%delta(), __FILE__, __LINE__)
   end subroutine simple_timer_report

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Test if t_stop .gt. t_start
   !---------------------------------------------------------------------------
   function simple_timer_has_stoped(this) result(has_stoped)
      class(simple_timer) :: this
      logical :: has_stoped
      has_stoped = this%t_stop .gt. this%t_start
   end function simple_timer_has_stoped

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Returns the difference between t_stop and t_start
   !---------------------------------------------------------------------------
   function simpletimer_delta(this) result(delta)
      class(simple_timer) :: this
      real(rp) :: delta
      delta = this%t_stop - this%t_start
   end function simpletimer_delta

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Measures the start time
   !---------------------------------------------------------------------------
   subroutine timer_start(this, label)
      class(timer) :: this
      character(len=*), intent(in) :: label
      type(simple_timer) :: new_timer
      if (this%auto_couple) then
         if (len(trim(this%coupled_label)) == 0) then
            this%coupled_label = trim(label)
         else
            this%coupled_label = trim(this%coupled_label)//'.'//trim(label)
         end if
         call new_timer%start(this%coupled_label)
      else
         call new_timer%start(label)
      end if
      call push(this%started_times, new_timer)
   end subroutine timer_start

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Measures the final time
   !---------------------------------------------------------------------------
   subroutine timer_stop(this, label)
      class(timer) :: this
      character(len=*), intent(in) :: label
      logical :: has_found, use_auto_couple
      integer :: i
      has_found = .False.
      use_auto_couple = this%auto_couple .and. endswith(this%coupled_label, label)
      do i = 1, size(this%started_times)
      if ((trim(this%started_times(i)%label) == trim(label) .or. (use_auto_couple .and. endswith(this%started_times(i)%label, label))) .and. .not. this%started_times(i)%has_stoped()) then
         if (use_auto_couple) this%coupled_label = this%coupled_label(:len(trim(this%coupled_label)) - len(trim(label)) - 1)
         call this%started_times(i)%stop()
         call push(this%finished_times, this%started_times(i))
         call remove(this%started_times, i)
         has_found = .True.
         continue
      end if
      end do
      if (.not. has_found) then
         call g_logger%error('Finishing '//trim(label)//' without starting first', __FILE__, __LINE__)
      end if
   end subroutine timer_stop

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Reports the time spent with label
   !---------------------------------------------------------------------------
   subroutine timer_report(this)
      class(timer) :: this
      type(report) :: rep
      integer :: i
      character(len=:), allocatable :: label
      real(rp) :: delta
      rep = report('total')
      do i = 1, size(this%finished_times)
         if (.not. this%finished_times(i)%has_stoped()) then
            call this%finished_times(i)%stop()
            call g_logger%error('Label time "'//trim(this%finished_times(i)%label)//'" not closed before report!', __FILE__, __LINE__)
         end if
         label = trim(rep%label)//'.'//trim(this%finished_times(i)%label)
         delta = this%finished_times(i)%delta()/this%clock_rate
         call g_logger%debug(fmt('A', label)//' '//fmt('F10.2', delta), __FILE__, __LINE__)
         call rep%add_value(label, delta)
      end do
      call rep%print_report('TIMER REPORT')
   end subroutine timer_report

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Pushes new element to list (private)
   !> @param[inout] list list to be pushed
   !> @param[in] mew_element element to be added
   !---------------------------------------------------------------------------
   subroutine push(list, new_element)
      type(simple_timer), dimension(:), allocatable, intent(inout) :: list
      type(simple_timer), intent(in) :: new_element
      type(simple_timer), dimension(:), allocatable :: aux_list
      integer :: new_size
      new_size = size(list) + 1
      allocate (aux_list(new_size))
      aux_list(:new_size - 1) = list(:new_size - 1)
      aux_list(new_size) = new_element
      call move_alloc(aux_list, list)
   end subroutine push

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Removes element from list by index (private)
   !> @param[inout] list list to be pushed
   !> @param[in] index index of the element to be removed
   !---------------------------------------------------------------------------
   subroutine remove(list, index)
      type(simple_timer), dimension(:), allocatable, intent(inout) :: list
      integer, intent(in) :: index
      type(simple_timer), dimension(:), allocatable :: aux_list
      integer :: new_size
      new_size = size(list) - 1
      if (index < 1 .or. index > size(list)) call g_logger%fatal('wrong index value index='//fmt('I0', index)//' size(list)='//fmt('I0', size(list)), __FILE__, __LINE__)
      allocate (aux_list(new_size))
      aux_list(:index - 1) = list(:index - 1)
      aux_list(index:) = list(index + 1:)
      call move_alloc(aux_list, list)
   end subroutine remove

end module timer_mod

