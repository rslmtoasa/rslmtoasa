!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Frozen_magnon
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> @brief Sweep configuration for the `frozen_magnon` post-processing mode.
!> @details Holds the list of spin-spiral `q_ss` points to evaluate, either
!>          from a compact external text file or the legacy namelist q_ss_list,
!>          plus the q-coordinate convention ('cartesian': Cartesian units of
!>          2*pi/alat, matching &hamiltonian q_ss; 'direct': reciprocal-lattice
!>          coordinates), the execution mode ('scf': full SCF at every point;
!>          'mft': full SCF at the
!>          reference point only, single-iteration band-energy pass for the
!>          rest), and the output filename. The sweep itself is driven from
!>          calculation.f90::post_processing_frozen_magnon, which owns the
!>          control/lattice/hamiltonian/self object orchestration; this type
!>          is pure namelist-driven configuration, following the same
!>          restore_to_default/build_from_file lifecycle as control/hamiltonian.
!------------------------------------------------------------------------------

module frozen_magnon_mod

   use precision_mod, only: rp
   use string_mod, only: sl, int2str, lower
   use logger_mod, only: g_logger
   implicit none

   private

   type, public :: frozen_magnon
      !> Input file name
      character(len=sl) :: fname
      !> 'scf' or 'mft' (magnetic force theorem / single-shot)
      character(len=10) :: mode
      !> Number of q-points in the sweep (>= 2: reference + at least one q)
      integer :: n_q
      !> Sweep points as read, shape (3, n_q). Column 1 is the reference point.
      real(rp), dimension(:, :), allocatable :: q_ss_list
      !> Optional external q-point file. If set, it replaces n_q_points/q_ss_list.
      character(len=sl) :: q_file
      !> 'cartesian' (2*pi/alat) or 'direct' (reciprocal lattice coordinates).
      character(len=16) :: q_coordinates
      !> Output data file (one row per q: q1 q2 q3 etot eband mtot_1..mtot_nrec omega)
      character(len=sl) :: output_file
   contains
      procedure :: build_from_file
      procedure :: restore_to_default
   end type frozen_magnon

   interface frozen_magnon
      procedure :: constructor
   end interface frozen_magnon

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] fname Namelist file with a `&frozen_magnon` table
   !> @return type(frozen_magnon)
   !---------------------------------------------------------------------------
   function constructor(fname) result(obj)
      type(frozen_magnon) :: obj
      character(len=*), intent(in) :: fname

      call obj%restore_to_default()
      call obj%build_from_file(fname)
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      class(frozen_magnon), intent(inout) :: this

      this%mode = 'mft'
      this%n_q = 0
      this%q_file = ''
      this%q_coordinates = 'cartesian'
      this%output_file = 'frozen_magnon.dat'
   end subroutine restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read the `&frozen_magnon` namelist from `fname` and validate it.
   !
   !> @param[in] fname Namelist file
   !> @note This is a true input boundary and raises fatal diagnostics for
   !>       invalid/missing sweep configuration.
   !---------------------------------------------------------------------------
   subroutine build_from_file(this, fname)
      class(frozen_magnon), intent(inout) :: this
      character(len=*), intent(in) :: fname
      integer :: iostatus, funit

      include 'include_codes/namelists/frozen_magnon.f90'

      mode = this%mode
      output_file = this%output_file
      q_file = this%q_file
      q_coordinates = this%q_coordinates
      n_q_points = 0
      q_ss_list = 0.0_rp

      this%fname = fname

      open (newunit=funit, file=fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=frozen_magnon, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%fatal('[frozen_magnon.build_from_file]: error reading &frozen_magnon namelist, '// &
                              'iostatus = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      this%mode = lower(trim(mode))
      if (this%mode /= 'mft' .and. this%mode /= 'scf') then
         call g_logger%fatal("[frozen_magnon.build_from_file]: mode must be 'mft' or 'scf'", __FILE__, __LINE__)
      end if

      this%output_file = output_file
      this%q_file = trim(q_file)
      this%q_coordinates = lower(trim(q_coordinates))
      if (len_trim(this%q_coordinates) == 0) this%q_coordinates = 'cartesian'
      if (this%q_coordinates /= 'cartesian' .and. this%q_coordinates /= 'direct') then
         call g_logger%fatal("[frozen_magnon.build_from_file]: q_coordinates must be 'cartesian' or 'direct'", &
                              __FILE__, __LINE__)
      end if

      if (len_trim(this%q_file) > 0) then
         call read_q_points_file(this, trim(this%q_file))
         return
      end if

      if (n_q_points < 2) then
         call g_logger%fatal('[frozen_magnon.build_from_file]: n_q_points must be >= 2 '// &
                              '(a reference point plus at least one sweep point)', __FILE__, __LINE__)
      end if
      if (n_q_points > frozen_magnon_max_q) then
         call g_logger%fatal('[frozen_magnon.build_from_file]: n_q_points exceeds the compiled-in limit of '// &
                              int2str(frozen_magnon_max_q), __FILE__, __LINE__)
      end if

      this%n_q = n_q_points
      if (allocated(this%q_ss_list)) deallocate (this%q_ss_list)
      allocate (this%q_ss_list(3, this%n_q))
      this%q_ss_list(:, :) = q_ss_list(:, 1:this%n_q)
   end subroutine build_from_file

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read q-points from an external text file.
   !
   !> @param[in] q_file File with first line n_q, followed by n_q rows of qx qy qz.
   !---------------------------------------------------------------------------
   subroutine read_q_points_file(this, q_file)
      class(frozen_magnon), intent(inout) :: this
      character(len=*), intent(in) :: q_file
      integer :: funit, iostatus, iq

      open (newunit=funit, file=q_file, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('[frozen_magnon.read_q_points_file]: q_file '//trim(q_file)//' not found', &
                              __FILE__, __LINE__)
      end if

      read (funit, *, iostat=iostatus) this%n_q
      if (iostatus /= 0) then
         close (funit)
         call g_logger%fatal('[frozen_magnon.read_q_points_file]: could not read q-point count from '// &
                              trim(q_file)//', iostatus = '//int2str(iostatus), __FILE__, __LINE__)
      end if

      if (this%n_q < 2) then
         close (funit)
         call g_logger%fatal('[frozen_magnon.read_q_points_file]: q-point count must be >= 2 '// &
                              '(a reference point plus at least one sweep point)', __FILE__, __LINE__)
      end if

      if (allocated(this%q_ss_list)) deallocate (this%q_ss_list)
      allocate (this%q_ss_list(3, this%n_q))

      do iq = 1, this%n_q
         read (funit, *, iostat=iostatus) this%q_ss_list(:, iq)
         if (iostatus /= 0) then
            close (funit)
            call g_logger%fatal('[frozen_magnon.read_q_points_file]: could not read q-point '// &
                                 int2str(iq)//' from '//trim(q_file)//', iostatus = '//int2str(iostatus), &
                                 __FILE__, __LINE__)
         end if
      end do
      close (funit)
   end subroutine read_q_points_file

end module frozen_magnon_mod
