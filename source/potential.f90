
!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Potential
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
!> Module to hold potential's parameters
!------------------------------------------------------------------------------

module potential_mod
   use os_mod, only: exists
   use precision_mod, only: rp
   use globals_mod, only: GLOBAL_DATABASE_FOLDER, GLOBAL_CHAR_SIZE
   use string_mod, only: path_join, sl, fmt, real2str, int2str
   use logger_mod, only: g_logger
   use namelist_generator_mod, only: namelist_generator
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   ! public functions
   public :: array_of_potentials

   type, public :: potential
      !> Orbital index. Determines the size of the Hamiltonian
      integer :: lmax
      !> Potential parameters \f$ C \f$ and \f$ \sqrt{\Delta} \f$
      !> 1st index 1 = s-orbital, 2 = p-orbital, 3 = d-orbital
      !> 2nd index 1 = spin-up, 2 = spin-dw
      real(rp), dimension(:, :), allocatable :: center_band
      real(rp), dimension(:, :), allocatable :: width_band
      real(rp), dimension(:, :), allocatable :: gravity_center
      real(rp), dimension(:, :), allocatable :: shifted_band ! center_band - gravity_center
      real(rp), dimension(:, :), allocatable :: obar
      !> Potential parameters treated internally in the code
      !> cx -> center of the band
      !> wx -> width of the band
      !> cex -> center of the band - gravity center of the band
      !> obx -> o parameter for the 'hoh' calculation
      !> 1st index 1 = s-orbital, 2-4 = p-orbital, 5-9 = d-orbital
      !> 2nd index 1 = spin-up, 2 = spin-dw
      complex(rp), dimension(:, :), allocatable :: cx, wx, cex, obx

      !> cx0 -> cx-up + cx-down
      !> cx1 -> cx-up - cx-dwown
      !> wx0 -> wx-up + wx-down
      !> wx1 -> wx-up - wx-down
      complex(rp), dimension(:), allocatable :: cx0, cx1, wx0, wx1, cex0, cex1, obx0, obx1
      !> Potential parameters Pl's defined as \f$ P_l = 0.5 - \frac{1}{\pi}arctg(D_{l})\f$
      !> 1st index 1 = s-orbital, 2 = p-orbital, 3 = d-orbital
      !> 2nd index 1 = spin-up, 2 = spin-dw
      real(rp), dimension(:, :), allocatable :: pl
      !> Moments as defined in Eq. 48. of Phys. Rev. B 43, 9538 (1991).
      !> 1st index 1 = s-orbital, 2 = p-orbital, 3 = d-orbital
      !> 2nd index 1 = spin-up, 2 = spin-dw
      real(rp), dimension(:, :, :), allocatable :: ql
      !> Potential parameters on the orthogonal basis
      real(rp), dimension(:, :), allocatable :: c, enu, ppar, qpar, srdel, vl, pnu, qi, dele
      !> Normalized magnetic moments
      real(rp), dimension(:), allocatable :: mom
      !> Direction of zeroth order band moments (i.e. non-normalized mom)
      real(rp), dimension(:), allocatable :: mom0
      !> Direction of first order band moments (0th order = mom0, 1st order = mom1)
      real(rp), dimension(:), allocatable :: mom1
      !> Magnetic moments
      real(rp) :: mx, my, mz, mtot
      !> Orbital moments (non-normalized)
      real(rp), dimension(:), allocatable :: lmom
      !> Orbital moments (not used)
      real(rp) :: lx, ly, lz, ltot
      !> Band variables
      real(rp), dimension(:), allocatable :: cshi, dw_l
      !> Wignzer Seitz Radius
      real(rp) :: ws_r
      ! Energy variables
      real(rp) :: sumec, sumev, etot, utot, ekin, rhoeps
      ! Madelung potential
      real(rp) :: vmad
      ! Spin-orbit coupling coefficients and Racah parameters
      real(rp), dimension(2) :: xi_p, xi_d, rac
      ! Hyperfine fields
      real(rp), dimension(2) :: hyper_field(2)
   contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: print_state
      procedure :: print_state_full
      procedure :: print_state_formatted
      procedure :: sizeof_potential
      procedure :: sizeof_potential_lite
      procedure :: sizeof_potential_full
      procedure :: flatten_potential
      procedure :: flatten_potential_lite
      procedure :: flatten_potential_full
      procedure :: expand_potential
      procedure :: expand_potential_lite
      procedure :: expand_potential_full
      procedure :: print_hyperfine
      procedure :: copy_mom_to_scal
      final :: destructor
   end type potential

   interface potential
      procedure :: constructor
   end interface potential

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] potential Namelist file in database
   !> @param[in] database Directory to database files with 'potential' namelist
   !> @return type(control)
   !---------------------------------------------------------------------------
   function constructor(label, database, reload) result(obj)
      type(potential) :: obj
      character(len=*), intent(in) :: label
      character(len=*), intent(in), optional :: database
      logical, optional, intent(in) :: reload
      logical :: reload_, built
      character(len=:), allocatable :: path_to_file
      character(len=sl), dimension(2) :: lst_path_to_file

      !reload_ = merge(reload, .True., present(reload))
      if (present(reload)) then
         reload_ = reload
      else
         reload_ = .true.
      end if
      call obj%restore_to_default()
      if (present(database)) then
         lst_path_to_file(1) = database
      else
         lst_path_to_file(1) = './'
      end if
      built = .false.
      if (reload_) then
         lst_path_to_file(2) = trim(label)//'_out.nml'
         path_to_file = path_join(lst_path_to_file)
         if (exists(path_to_file)) then
            call g_logger%info('Reading potential from previous result '//trim(path_to_file), __FILE__, __LINE__)
            call obj%build_from_file(path_to_file)
            built = .true.
         end if
      end if
      if (.not. built) then
         lst_path_to_file(2) = trim(label)//'.nml'
         path_to_file = path_join(lst_path_to_file)
         if (exists(path_to_file)) then
            call obj%build_from_file(path_to_file)
         else
            lst_path_to_file(1) = GLOBAL_DATABASE_FOLDER
            path_to_file = path_join(lst_path_to_file)
            if (exists(path_to_file)) then
               call obj%build_from_file(path_to_file)
            else
               call g_logger%fatal('Element '//trim(label)//' not found in any database', __FILE__, __LINE__)
            end if
         end if
      end if
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(potential) :: this
   end subroutine destructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read parameters from input file
   !
   !> @param[in] fname Namelist file
   !---------------------------------------------------------------------------
   subroutine build_from_file(this, fname)
      class(potential), intent(inout) :: this
      character(len=*), intent(in) :: fname

      ! variables associated with the reading processes
      integer :: iostatus, funit

      ! Local variables


      include 'include_codes/namelists/potential.f90'

      ! Save previous values
      ws_r = this%ws_r
      sumec = this%sumec
      sumev = this%sumev
      etot = this%etot
      utot = this%utot
      ekin = this%ekin
      rhoeps = this%rhoeps
      lmax = this%lmax
      vmad = this%vmad
      xi_p = this%xi_p
      xi_d = this%xi_d
      rac = this%rac

      ! Normalize the magnetic moments
      this%mom(:) = this%mom(:) / norm2(this%mom(:))

      call move_alloc(this%ql, ql)
      call move_alloc(this%mom, mom)
      call move_alloc(this%lmom, lmom)
      call move_alloc(this%pl, pl)

      call move_alloc(this%center_band, center_band)
      call move_alloc(this%width_band, width_band)
      call move_alloc(this%gravity_center, gravity_center)
      call move_alloc(this%shifted_band, shifted_band)
      call move_alloc(this%obar, obar)

      call move_alloc(this%c, c)
      call move_alloc(this%enu, enu)
      call move_alloc(this%ppar, ppar)
      call move_alloc(this%qpar, qpar)
      call move_alloc(this%srdel, srdel)
      call move_alloc(this%vl, vl)

      open (newunit=funit, file=fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//fmt('A', fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=par, iostat=iostatus)
      if ( &
         iostatus /= 0 .and. &
         .not. IS_IOSTAT_END(iostatus) .and. &
         iostatus /= 5010 & ! LIBERROR_READ_VALUE (according to https://www.hep.manchester.ac.uk/u/samt/misc/gfortran_errors.html)
         ) then
         call g_logger%fatal('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%fatal('iostatus = '//fmt('I0', iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      ! Setting user values
      this%ws_r = ws_r
      this%sumec = sumec
      this%sumev = sumev
      this%etot = etot
      this%utot = utot
      this%ekin = ekin
      this%rhoeps = rhoeps
      this%lmax = lmax
      this%vmad = vmad
      this%xi_p = xi_p
      this%xi_d = xi_d
      this%rac = rac

      call move_alloc(center_band, this%center_band)
      call move_alloc(width_band, this%width_band)
      call move_alloc(gravity_center, this%gravity_center)
      call move_alloc(shifted_band, this%shifted_band)
      call move_alloc(obar, this%obar)

      mom(:) = mom(:) / norm2(mom(:))

      call move_alloc(mom, this%mom)
      call move_alloc(lmom, this%lmom)
      call move_alloc(pl, this%pl)
      call move_alloc(ql, this%ql)

      call move_alloc(c, this%c)
      call move_alloc(enu, this%enu)
      call move_alloc(ppar, this%ppar)
      call move_alloc(qpar, this%qpar)
      call move_alloc(srdel, this%srdel)
      call move_alloc(vl, this%vl)
   end subroutine build_from_file

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      implicit none
      class(potential), intent(out) :: this

      this%lmax = 2

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('potential.center_band', this%center_band, (/this%lmax + 1, 2/))
      call g_safe_alloc%allocate('potential.width_band', this%width_band, (/this%lmax + 1, 2/))
      allocate (this%pl(0:this%lmax, 2))
      call g_safe_alloc%report_allocate('potential.pl', this%pl)
      call g_safe_alloc%allocate('potential.gravity_center', this%gravity_center, (/this%lmax + 1, 2/))
      allocate (this%ql(3, 0:this%lmax, 2))
      call g_safe_alloc%report_allocate('potential.ql', this%ql)
      call g_safe_alloc%allocate('potential.cx', this%cx, (/(this%lmax + 1)**2, 2/))
      call g_safe_alloc%allocate('potential.wx', this%wx, (/(this%lmax + 1)**2, 2/))
      call g_safe_alloc%allocate('potential.cex', this%cex, (/(this%lmax + 1)**2, 2/))
      call g_safe_alloc%allocate('potential.obx', this%obx, (/(this%lmax + 1)**2, 2/))
      call g_safe_alloc%allocate('potential.cx0', this%cx0, (/(this%lmax + 1)**2/))
      call g_safe_alloc%allocate('potential.cx1', this%cx1, (/(this%lmax + 1)**2/))
      call g_safe_alloc%allocate('potential.wx0', this%wx0, (/(this%lmax + 1)**2/))
      call g_safe_alloc%allocate('potential.wx1', this%wx1, (/(this%lmax + 1)**2/))
      call g_safe_alloc%allocate('potential.cex0', this%cex0, (/(this%lmax + 1)**2/))
      call g_safe_alloc%allocate('potential.cex1', this%cex1, (/(this%lmax + 1)**2/))
      call g_safe_alloc%allocate('potential.obx0', this%obx0, (/(this%lmax + 1)**2/))
      call g_safe_alloc%allocate('potential.obx1', this%obx1, (/(this%lmax + 1)**2/))
      call g_safe_alloc%allocate('potential.mom', this%mom, (/3/))
      call g_safe_alloc%allocate('potential.lmom', this%lmom, (/3/))
      call g_safe_alloc%allocate('potential.mom0', this%mom0, (/3/))
      call g_safe_alloc%allocate('potential.mom1', this%mom1, (/3/))
      call g_safe_alloc%allocate('potential.cshi', this%cshi, (/(2 * (this%lmax + 1))**2/))
      call g_safe_alloc%allocate('potential.dw_l', this%dw_l, (/(2 * (this%lmax + 1))**2/))
      allocate (this%c(0:this%lmax, 2))
      call g_safe_alloc%report_allocate('potential.c', this%c)
      allocate (this%enu(0:this%lmax, 2))
      call g_safe_alloc%report_allocate('potential.enu', this%enu)
      allocate (this%ppar(0:this%lmax, 2))
      call g_safe_alloc%report_allocate('potential.ppar', this%ppar)
      allocate (this%qpar(0:this%lmax, 2))
      call g_safe_alloc%report_allocate('potential.qpar', this%qpar)
      allocate (this%srdel(0:this%lmax, 2))
      call g_safe_alloc%report_allocate('potential.srdel', this%srdel)
      allocate (this%vl(0:this%lmax, 2))
      call g_safe_alloc%report_allocate('potential.vl', this%vl)
      allocate (this%pnu(0:this%lmax, 2))
      call g_safe_alloc%report_allocate('potential.pnu', this%pnu)
      call g_safe_alloc%allocate('potential.shifted_band', this%shifted_band, (/this%lmax + 1, 2/))
      call g_safe_alloc%allocate('potential.obar', this%obar, (/this%lmax + 1, 2/))
#else
      allocate (this%center_band(this%lmax + 1, 2), this%width_band(this%lmax + 1, 2), this%pl(0:this%lmax, 2))
      allocate (this%gravity_center(this%lmax + 1, 2), this%ql(3, 0:this%lmax, 2), this%cx((this%lmax + 1)**2, 2), this%wx((this%lmax + 1)**2, 2))
      allocate (this%cex((this%lmax + 1)**2, 2), this%obx((this%lmax + 1)**2, 2), this%cx0((this%lmax + 1)**2))
      allocate (this%cx1((this%lmax + 1)**2), this%wx0((this%lmax + 1)**2), this%wx1((this%lmax + 1)**2))
      allocate (this%cex0((this%lmax + 1)**2), this%cex1((this%lmax + 1)**2), this%obx0((this%lmax + 1)**2), this%obx1((this%lmax + 1)**2))
      allocate (this%mom(3), this%lmom(3), this%cshi((2 * (this%lmax + 1))**2), this%dw_l((2 * (this%lmax + 1))**2))
      allocate (this%mom0(3), this%mom1(3))
      allocate (this%c(0:this%lmax, 2), this%enu(0:this%lmax, 2), this%ppar(0:this%lmax, 2), this%qpar(0:this%lmax, 2), &
                this%srdel(0:this%lmax, 2), this%vl(0:this%lmax, 2), this%pnu(0:this%lmax, 2), this%qi(0:this%lmax, 2), this%dele(0:this%lmax, 2))
      allocate (this%shifted_band(this%lmax + 1, 2), this%obar(this%lmax + 1, 2))
#endif

      this%xi_p(:) = 0.0D0
      this%xi_d(:) = 0.0D0
      this%rac(:) = 0.0D0
      this%ws_r = 0.0D0
      this%center_band(:, :) = 0.0D0
      this%width_band(:, :) = 0.0D0
      this%shifted_band(:, :) = 0.0D0
      this%obar(:, :) = 0.0D0
      this%sumec = 0.0D0
      this%sumev = 0.0D0
      this%etot = 0.0D0
      this%utot = 0.0D0
      this%ekin = 0.0D0
      this%rhoeps = 0.0D0
      this%vmad = 0.0D0
      this%pl(:, :) = 0.0D0
      this%ql(:, :, :) = 0.0D0
      this%cx(:, :) = 0.0D0
      this%wx(:, :) = 0.0D0
      this%cex(:, :) = 0.0D0
      this%obx(:, :) = 0.0D0
      this%cx0(:) = 0.0D0
      this%cx1(:) = 0.0D0
      this%wx0(:) = 0.0D0
      this%wx1(:) = 0.0D0
      this%mom(:) = [0.0D0, 0.0D0, 1.0D0]
      this%lmom(:) = [0.0D0, 0.0D0, 0.0D0]
      this%mom0(:) = [0.0D0, 0.0D0, 0.0D0]
      this%mom1(:) = [0.0D0, 0.0D0, 0.0D0]
      this%cshi(:) = 0.0D0
      this%dw_l(:) = 1.0D0
      this%c(:, :) = 0.0D0
      this%enu(:, :) = 0.0D0
      this%ppar(:, :) = 0.0D0
      this%qpar(:, :) = 0.0D0
      this%srdel(:, :) = 0.0D0
      this%vl(:, :) = 0.0D0
      this%pnu(:, :) = 0.0D0
      this%qi(:, :) = 0.0D0
      this%dele(:, :) = 0.0D0
   end subroutine restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
   !> @param[in] unit File unit used to write namelist
   !> @param[in] file File name used to write namelist
   !---------------------------------------------------------------------------
   subroutine print_state(this, unit, file)
      implicit none
      class(potential), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

      include 'include_codes/namelists/potential.f90'

      mom = this%mom
      lmom = this%lmom
      !mom0 = this%mom0
      !mom1 = this%mom1
      pl = this%pl
      ql = this%ql
      center_band = this%center_band
      width_band = this%width_band
      gravity_center = this%gravity_center
      sumec = this%sumec
      sumev = this%sumev
      etot = this%etot
      utot = this%utot
      ekin = this%ekin
      rhoeps = this%rhoeps
      vmad = this%vmad
      lmax = this%lmax
      ws_r = this%ws_r
      c = this%c
      enu = this%enu
      ppar = this%ppar
      qpar = this%qpar
      srdel = this%srdel
      vl = this%vl

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         write (unit, nml=par)
      else if (present(file)) then
         open (unit=newunit, file=file)
         write (newunit, nml=par)
         close (newunit)
      else
         write (*, nml=par)
      end if

   end subroutine print_state

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
   !> @param[in] unit File unit used to write namelist
   !> @param[in] file File name used to write namelist
   !---------------------------------------------------------------------------
   subroutine print_state_full(this, unit, file)
      implicit none
      class(potential), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

      include 'include_codes/namelists/potential.f90'

      center_band_s_up = this%center_band(1, 1)
      center_band_s_dw = this%center_band(1, 2)
      center_band_p_up = this%center_band(2, 1)
      center_band_p_dw = this%center_band(2, 2)
      center_band_d_up = this%center_band(3, 1)
      center_band_d_dw = this%center_band(3, 2)

      width_band_s_up = this%width_band(1, 1)
      width_band_s_dw = this%width_band(1, 2)
      width_band_p_up = this%width_band(2, 1)
      width_band_p_dw = this%width_band(2, 2)
      width_band_d_up = this%width_band(3, 1)
      width_band_d_dw = this%width_band(3, 2)

      mom = this%mom
      lmom = this%lmom
      !mom0 = this%mom0
      !mom1 = this%mom1

      cx = this%cx
      wx = this%wx
      cx0 = this%cx0
      wx0 = this%wx0
      cx1 = this%cx1
      wx1 = this%wx1

      pl = this%pl
      cshi = this%cshi
      dw_l = this%dw_l
      cex = this%cex
      obx = this%obx

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         write (unit, nml=par)
      else if (present(file)) then
         open (unit=newunit, file=file)
         write (newunit, nml=par)
         close (newunit)
      else
         write (*, nml=par)
      end if

   end subroutine print_state_full

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
   !> @param[in] unit File unit used to write namelist
   !> @param[in] file File name used to write namelist
   !---------------------------------------------------------------------------
   subroutine print_state_formatted(this, unit, file)
      implicit none
      class(potential), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file


      type(namelist_generator) :: nml

      nml = namelist_generator('par')

      call nml%add('center_band', this%center_band)
      call nml%add('width_band', this%width_band)
      call nml%add('gravity_center', this%gravity_center)
      call nml%add('shifted_band', this%shifted_band)
      call nml%add('obar', this%obar)
      call nml%add('sumec', this%sumec)
      call nml%add('sumev', this%sumev)
      call nml%add('etot', this%etot)
      call nml%add('utot', this%utot)
      call nml%add('ekin', this%ekin)
      call nml%add('rhoeps', this%rhoeps)
      call nml%add('c', this%c)
      call nml%add('enu', this%enu)
      call nml%add('ppar', this%ppar)
      call nml%add('qpar', this%qpar)
      call nml%add('srdel', this%srdel)
      call nml%add('vl', this%vl)
      call nml%add('pl', this%pl)
      call nml%add('mom', this%mom)
      call nml%add('lmom', this%lmom)
      call nml%add('ws_r', this%ws_r)
      call nml%add('ql', this%ql)
      call nml%add('lmax', this%lmax)
      call nml%add('vmad', this%vmad)
      call nml%add('xi_p', this%xi_p)
      call nml%add('xi_d', this%xi_d)
      call nml%add('rac', this%rac)

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         call nml%generate_namelist(unit=unit)
      else if (present(file)) then
         call nml%generate_namelist(file=file)
      else
         call nml%generate_namelist()
      end if
   end subroutine print_state_formatted

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Build an array of potentials
   !
   !> @param[in] potential List of labels in database to build potentials
   !> @param[in] database Directory to database files with 'potential' namelist
   !> @return type(potential), dimension(:), allocatable
   !---------------------------------------------------------------------------
   function array_of_potentials(potentials, database)
      type(potential), dimension(:), allocatable :: array_of_potentials
      character(len=*), dimension(:), intent(in) :: potentials
      character(len=*), intent(in), optional :: database
      integer :: i
      allocate (array_of_potentials(size(potentials)))
      if (present(database)) then
         do i = 1, size(potentials)
            array_of_potentials(i) = potential(potentials(i), database)
         end do
      else
         do i = 1, size(potentials)
            array_of_potentials(i) = potential(potentials(i))
         end do
      end if
   end function array_of_potentials

   subroutine flatten_potential_lite(this, flat_array)
      class(potential), intent(inout) :: this
      real(rp), dimension(:), intent(out) :: flat_array

      ! Determine the total size needed for the flat array
      integer :: total_size
      total_size = this%sizeof_potential_lite()
      !total_size = size(this%gravity_center) + size(this%ql) + size(this%mom) + 1 !mtot

      ! Flatten the data using PACK
      flat_array(1:total_size) = (/pack(this%gravity_center, .true.), &
                                   pack(this%ql, .true.), pack(this%mom, .true.), this%mtot/)
      return

   end subroutine flatten_potential_lite

   subroutine flatten_potential(this, flat_array)
      class(potential), intent(inout) :: this
      real(rp), dimension(:), intent(out) :: flat_array
      !real(rp), dimension(:), allocatable, intent(out) :: flat_array

      ! Determine the total size needed for the flat array
      integer :: total_size
      total_size = this%sizeof_potential()
      !total_size = size(this%center_band) + size(this%width_band) + size(this%gravity_center) + &
      !             size(this%shifted_band) + size(this%pl) + size(this%ql) + size(this%c) + &
      !             size(this%enu) + size(this%ppar) + size(this%qpar) + size(this%srdel) + &
      !             size(this%vl) + size(this%pnu) + size(this%qi) + size(this%dele) + &
      !             size(this%mom) + 1 !mtot

      ! Flatten the data using PACK
      flat_array(1:total_size) = (/pack(this%center_band, .true.), pack(this%width_band, .true.), &
                                   pack(this%gravity_center, .true.), pack(this%shifted_band, .true.), &
                                   pack(this%pl, .true.), pack(this%ql, .true.), pack(this%c, .true.), &
                                   pack(this%enu, .true.), pack(this%ppar, .true.), pack(this%qpar, .true.), &
                                   pack(this%srdel, .true.), pack(this%vl, .true.), pack(this%pnu, .true.), &
                                   pack(this%qi, .true.), pack(this%dele, .true.), pack(this%mom, .true.), &
                                   pack(this%lmom, .true.), this%mtot/)
      return

   end subroutine flatten_potential

   subroutine flatten_potential_full(this, flat_array)
      class(potential), intent(inout) :: this
      real(rp), dimension(:), intent(out) :: flat_array
      !real(rp), dimension(:), allocatable, intent(out) :: flat_array

      ! Determine the total size needed for the flat array
      integer :: total_size
      total_size = this%sizeof_potential_full()
      !total_size = size(this%center_band) + size(this%width_band) + size(this%gravity_center) + &
      !             size(this%shifted_band) + size(this%pl) + size(this%ql) + size(this%c) + &
      !             size(this%enu) + size(this%ppar) + size(this%qpar) + size(this%srdel) + &
      !             size(this%vl) + size(this%pnu) + size(this%qi) + size(this%dele) + &
      !             size(this%mom) + 1 !mtot

      ! Flatten the data using PACK
      flat_array(1:total_size) = (/pack(this%center_band, .true.), &
                                   pack(this%width_band, .true.), &
                                   pack(this%gravity_center, .true.), &
                                   pack(this%shifted_band, .true.), &
                                   pack(this%obar, .true.), &
                              !!! pack(real(this%cx), .true.), &
                              !!! pack(real(this%wx), .true.), &
                              !!! pack(aimag(this%wx), .true.), &
                              !!! pack(real(this%cex), .true.), &
                              !!! pack(aimag(this%cex), .true.), &
                              !!! pack(real(this%obx), .true.), &
                              !!! pack(aimag(this%obx), .true.), &
                                   pack(this%pl, .true.), &
                                   pack(this%ql, .true.), &
                                   pack(this%c, .true.), &
                                   pack(this%enu, .true.), &
                                   pack(this%ppar, .true.), &
                                   pack(this%qpar, .true.), &
                                   pack(this%srdel, .true.), &
                                   pack(this%vl, .true.), &
                                   pack(this%pnu, .true.), &
                                   pack(this%qi, .true.), &
                                   pack(this%dele, .true.), &
                                   pack(this%mom, .true.), &
                                   pack(this%lmom, .true.), &
                                   pack(this%xi_p, .true.), &
                                   pack(this%xi_d, .true.), &
                                   pack(this%rac, .true.), &
                                   this%sumec, &
                                   this%sumev, &
                                   this%etot, &
                                   this%utot, &
                                   this%ekin, &
                                   this%rhoeps, &
                                   this%vmad, &
                                   this%mtot/)
      return

   end subroutine flatten_potential_full

   function sizeof_potential_lite(this) result(total_size)
      class(potential), intent(inout) :: this
      integer :: total_size

      ! Determine the total size needed for the flat array
      total_size = size(this%gravity_center) + size(this%ql) &
                   + size(this%mom) + 1 !mtot

   end function sizeof_potential_lite

   function sizeof_potential(this) result(total_size)
      class(potential), intent(inout) :: this
      integer :: total_size

      ! Determine the total size needed for the flat array
      total_size = size(this%center_band) + size(this%width_band) + size(this%gravity_center) + &
                   size(this%shifted_band) + size(this%pl) + size(this%ql) + size(this%c) + &
                   size(this%enu) + size(this%ppar) + size(this%qpar) + size(this%srdel) + &
                   size(this%vl) + size(this%pnu) + size(this%qi) + size(this%dele) + &
                   size(this%mom) + size(this%lmom) + 1 !mtot

   end function sizeof_potential

   function sizeof_potential_full(this) result(total_size)
      class(potential), intent(inout) :: this
      integer :: total_size

      ! Determine the total size needed for the flat array
      total_size = size(this%center_band) &
                   + size(this%width_band) &
                   + size(this%gravity_center) &
                   + size(this%shifted_band) &
                   + size(this%obar) &
                !!! + size(this%cx)*2 &  !complex
                !!! + size(this%wx)*2 &  !complex
                !!! + size(this%cex)*2 &  !complex
                !!! + size(this%obx)*2 &  !complex
                   + size(this%pl) &
                   + size(this%ql) &
                   + size(this%c) &
                   + size(this%enu) &
                   + size(this%ppar) &
                   + size(this%qpar) &
                   + size(this%srdel) &
                   + size(this%vl) &
                   + size(this%pnu) &
                   + size(this%qi) &
                   + size(this%dele) &
                   + size(this%mom) &
                   + size(this%lmom) &
                   + size(this%xi_p) &
                   + size(this%xi_d) &
                   + size(this%rac) &
                   + 1 & ! sumec
                   + 1 & ! sumev
                   + 1 & ! etot
                   + 1 & ! utot
                   + 1 & ! ekin
                   + 1 & ! rhoeps
                   + 1 & ! vmad
                   + 1 !mtot

   end function sizeof_potential_full

   subroutine expand_potential_lite(this, flat_array)
      class(potential), intent(inout) :: this
      real(rp), dimension(*), intent(in) :: flat_array
      integer :: d_pos, d_size



      ! Determine the total size needed for the flat array
      integer :: total_size
      total_size = this%sizeof_potential_lite()
      !total_size = size(this%gravity_center) +  size(this%ql) &
      !            + size(this%mom) + 1 !mtot

      d_pos = 1
      ! gravity_center
      d_size = size(this%gravity_center)
      this%gravity_center = reshape(flat_array(d_pos:d_pos + d_size), shape(this%gravity_center))
      d_pos = d_pos + d_size
      ! ql
      d_size = size(this%ql)
      this%ql = reshape(flat_array(d_pos:d_pos + d_size), shape(this%ql))
      d_pos = d_pos + d_size
      ! mom
      d_size = size(this%mom)
      this%mom = reshape(flat_array(d_pos:d_pos + d_size), shape(this%mom))
      d_pos = d_pos + d_size
      ! mtot
      d_size = 1
      this%mtot = flat_array(d_pos)

      return

   end subroutine expand_potential_lite

   subroutine expand_potential(this, flat_array)
      class(potential), intent(inout) :: this
      real(rp), dimension(*), intent(in) :: flat_array
      integer :: d_pos, d_size



      ! Determine the total size needed for the flat array
      integer :: total_size

      total_size = this%sizeof_potential()
      !total_size = size(this%center_band) + size(this%width_band) + size(this%gravity_center) + &
      !             size(this%shifted_band) + size(this%pl) + size(this%ql) + size(this%c) + &
      !             size(this%enu) + size(this%ppar) + size(this%qpar) + size(this%srdel) + &
      !             size(this%vl) + size(this%pnu) + size(this%qi) + size(this%dele) + &
      !             size(this%mom) + 1 !mtot

      d_pos = 1
      ! center_band
      d_size = size(this%center_band)
      this%center_band = reshape(flat_array(d_pos:d_pos + d_size), shape(this%center_band))
      d_pos = d_pos + d_size
      ! width_band
      d_size = size(this%width_band)
      this%width_band = reshape(flat_array(d_pos:d_pos + d_size), shape(this%width_band))
      d_pos = d_pos + d_size
      ! gravity_center
      d_size = size(this%gravity_center)
      this%gravity_center = reshape(flat_array(d_pos:d_pos + d_size), shape(this%gravity_center))
      d_pos = d_pos + d_size
      ! shifted_band
      d_size = size(this%shifted_band)
      this%shifted_band = reshape(flat_array(d_pos:d_pos + d_size), shape(this%shifted_band))
      d_pos = d_pos + d_size
      ! pl
      d_size = size(this%pl)
      this%pl = reshape(flat_array(d_pos:d_pos + d_size), shape(this%pl))
      d_pos = d_pos + d_size
      ! ql
      d_size = size(this%ql)
      this%ql = reshape(flat_array(d_pos:d_pos + d_size), shape(this%ql))
      d_pos = d_pos + d_size
      ! c
      d_size = size(this%c)
      this%c = reshape(flat_array(d_pos:d_pos + d_size), shape(this%c))
      d_pos = d_pos + d_size
      ! enu
      d_size = size(this%enu)
      this%enu = reshape(flat_array(d_pos:d_pos + d_size), shape(this%enu))
      d_pos = d_pos + d_size
      ! ppar
      d_size = size(this%ppar)
      this%ppar = reshape(flat_array(d_pos:d_pos + d_size), shape(this%ppar))
      d_pos = d_pos + d_size
      ! qpar
      d_size = size(this%qpar)
      this%qpar = reshape(flat_array(d_pos:d_pos + d_size), shape(this%qpar))
      d_pos = d_pos + d_size
      ! srdel
      d_size = size(this%srdel)
      this%srdel = reshape(flat_array(d_pos:d_pos + d_size), shape(this%srdel))
      d_pos = d_pos + d_size
      ! vl
      d_size = size(this%vl)
      this%vl = reshape(flat_array(d_pos:d_pos + d_size), shape(this%vl))
      d_pos = d_pos + d_size
      ! pnu
      d_size = size(this%pnu)
      this%pnu = reshape(flat_array(d_pos:d_pos + d_size), shape(this%pnu))
      d_pos = d_pos + d_size
      ! qi
      d_size = size(this%qi)
      this%qi = reshape(flat_array(d_pos:d_pos + d_size), shape(this%qi))
      d_pos = d_pos + d_size
      ! dele
      d_size = size(this%dele)
      this%dele = reshape(flat_array(d_pos:d_pos + d_size), shape(this%dele))
      d_pos = d_pos + d_size
      ! mom
      d_size = size(this%mom)
      this%mom = reshape(flat_array(d_pos:d_pos + d_size), shape(this%mom))
      d_pos = d_pos + d_size
      ! lmom
      d_size = size(this%lmom)
      this%lmom = reshape(flat_array(d_pos:d_pos + d_size), shape(this%lmom))
      d_pos = d_pos + d_size
      ! mtot
      d_size = 1
      this%mtot = flat_array(d_pos)

      return

   end subroutine expand_potential

   subroutine expand_potential_full(this, flat_array)
      class(potential), intent(inout) :: this
      real(rp), dimension(*), intent(in) :: flat_array
      integer :: d_pos, d_size



      ! Determine the total size needed for the flat array
      integer :: total_size
      complex(rp)  :: jimag = (0.0_RP, 1.0_RP)

      total_size = this%sizeof_potential_full()
      !total_size = size(this%center_band) + size(this%width_band) + size(this%gravity_center) + &
      !             size(this%shifted_band) + size(this%pl) + size(this%ql) + size(this%c) + &
      !             size(this%enu) + size(this%ppar) + size(this%qpar) + size(this%srdel) + &
      !             size(this%vl) + size(this%pnu) + size(this%qi) + size(this%dele) + &
      !             size(this%mom) + 1 !mtot

      d_pos = 1
      ! center_band
      d_size = size(this%center_band)
      this%center_band = reshape(flat_array(d_pos:d_pos + d_size), shape(this%center_band))
      d_pos = d_pos + d_size
      ! width_band
      d_size = size(this%width_band)
      this%width_band = reshape(flat_array(d_pos:d_pos + d_size), shape(this%width_band))
      d_pos = d_pos + d_size
      ! gravity_center
      d_size = size(this%gravity_center)
      this%gravity_center = reshape(flat_array(d_pos:d_pos + d_size), shape(this%gravity_center))
      d_pos = d_pos + d_size
      ! shifted_band
      d_size = size(this%shifted_band)
      this%shifted_band = reshape(flat_array(d_pos:d_pos + d_size), shape(this%shifted_band))
      d_pos = d_pos + d_size
      ! obar
      d_size = size(this%obar)
      this%obar = reshape(flat_array(d_pos:d_pos + d_size), shape(this%obar))
      d_pos = d_pos + d_size
    !!! ! cx (complex)
    !!! d_size = size(this%cx)
    !!! this%cx = reshape(flat_array(d_pos:d_pos+d_size), shape(this%cx))
    !!! d_pos = d_pos + d_size
    !!! this%cx = this%cx + jimag * reshape(flat_array(d_pos:d_pos+d_size), shape(this%cx))
    !!! d_pos = d_pos + d_size
    !!! ! wx (complex)
    !!! d_size = size(this%wx)
    !!! this%wx = reshape(flat_array(d_pos:d_pos+d_size), shape(this%wx))
    !!! d_pos = d_pos + d_size
    !!! this%wx = this%wx + jimag * reshape(flat_array(d_pos:d_pos+d_size), shape(this%wx))
    !!! d_pos = d_pos + d_size
    !!! ! cex (complex)
    !!! d_size = size(this%cex)
    !!! this%cex = reshape(flat_array(d_pos:d_pos+d_size), shape(this%cex))
    !!! d_pos = d_pos + d_size
    !!! this%cex = this%cex + jimag * reshape(flat_array(d_pos:d_pos+d_size), shape(this%cex))
    !!! d_pos = d_pos + d_size
    !!! ! obx (complex)
    !!! d_size = size(this%obx)
    !!! this%obx = reshape(flat_array(d_pos:d_pos+d_size), shape(this%obx))
    !!! d_pos = d_pos + d_size
    !!! this%obx = this%obx + jimag * reshape(flat_array(d_pos:d_pos+d_size), shape(this%obx))
    !!! d_pos = d_pos + d_size
      ! pl
      d_size = size(this%pl)
      this%pl = reshape(flat_array(d_pos:d_pos + d_size), shape(this%pl))
      d_pos = d_pos + d_size
      ! ql
      d_size = size(this%ql)
      this%ql = reshape(flat_array(d_pos:d_pos + d_size), shape(this%ql))
      d_pos = d_pos + d_size
      ! c
      d_size = size(this%c)
      this%c = reshape(flat_array(d_pos:d_pos + d_size), shape(this%c))
      d_pos = d_pos + d_size
      ! enu
      d_size = size(this%enu)
      this%enu = reshape(flat_array(d_pos:d_pos + d_size), shape(this%enu))
      d_pos = d_pos + d_size
      ! ppar
      d_size = size(this%ppar)
      this%ppar = reshape(flat_array(d_pos:d_pos + d_size), shape(this%ppar))
      d_pos = d_pos + d_size
      ! qpar
      d_size = size(this%qpar)
      this%qpar = reshape(flat_array(d_pos:d_pos + d_size), shape(this%qpar))
      d_pos = d_pos + d_size
      ! srdel
      d_size = size(this%srdel)
      this%srdel = reshape(flat_array(d_pos:d_pos + d_size), shape(this%srdel))
      d_pos = d_pos + d_size
      ! vl
      d_size = size(this%vl)
      this%vl = reshape(flat_array(d_pos:d_pos + d_size), shape(this%vl))
      d_pos = d_pos + d_size
      ! pnu
      d_size = size(this%pnu)
      this%pnu = reshape(flat_array(d_pos:d_pos + d_size), shape(this%pnu))
      d_pos = d_pos + d_size
      ! qi
      d_size = size(this%qi)
      this%qi = reshape(flat_array(d_pos:d_pos + d_size), shape(this%qi))
      d_pos = d_pos + d_size
      ! dele
      d_size = size(this%dele)
      this%dele = reshape(flat_array(d_pos:d_pos + d_size), shape(this%dele))
      d_pos = d_pos + d_size
      ! mom
      d_size = size(this%mom)
      this%mom = reshape(flat_array(d_pos:d_pos + d_size), shape(this%mom))
      d_pos = d_pos + d_size
      ! lmom
      d_size = size(this%lmom)
      this%lmom = reshape(flat_array(d_pos:d_pos + d_size), shape(this%lmom))
      d_pos = d_pos + d_size
      ! xi_p
      d_size = size(this%xi_p)
      this%xi_p = reshape(flat_array(d_pos:d_pos + d_size), shape(this%xi_p))
      d_pos = d_pos + d_size
      ! xi_d
      d_size = size(this%xi_d)
      this%xi_d = reshape(flat_array(d_pos:d_pos + d_size), shape(this%xi_d))
      d_pos = d_pos + d_size
      ! xi_d
      d_size = size(this%rac)
      this%rac = reshape(flat_array(d_pos:d_pos + d_size), shape(this%rac))
      d_pos = d_pos + d_size
      ! sumec
      d_size = 1
      this%sumec = flat_array(d_pos)
      d_pos = d_pos + d_size
      ! sumev
      d_size = 1
      this%sumev = flat_array(d_pos)
      d_pos = d_pos + d_size
      ! etot
      d_size = 1
      this%etot = flat_array(d_pos)
      d_pos = d_pos + d_size
      ! utot
      d_size = 1
      this%utot = flat_array(d_pos)
      d_pos = d_pos + d_size
      ! ekin
      d_size = 1
      this%ekin = flat_array(d_pos)
      d_pos = d_pos + d_size
      ! rhoeps
      d_size = 1
      this%rhoeps = flat_array(d_pos)
      d_pos = d_pos + d_size
      ! vmad
      d_size = 1
      this%vmad = flat_array(d_pos)
      d_pos = d_pos + d_size
      ! mtot
      d_size = 1
      this%mtot = flat_array(d_pos)

      return

   end subroutine expand_potential_full

   subroutine print_hyperfine(this, ia)
      class(potential), intent(in) :: this
      integer, intent(in) :: ia

      real(rp) :: hyper_field_tot
      real(rp), dimension(2) :: hyper_field
      character*10 :: fmt

      hyper_field = this%hyper_field
      hyper_field_tot = sum(hyper_field)
      fmt = "(F 8.3)"

      call g_logger%info('Hyperfine field for atom '//int2str(ia)// &
                         ':  H_core='//real2str(hyper_field(1), fmt)//' T, H_val= '//real2str(hyper_field(2), fmt)//' T.' &
                         , __FILE__, __LINE__)

   end subroutine print_hyperfine

   subroutine copy_mom_to_scal(this)
      class(potential) :: this

      this%mx = this%mom(1)
      this%my = this%mom(2)
      this%mz = this%mom(3)

      this%lx = this%lmom(1)
      this%ly = this%lmom(2)
      this%lz = this%lmom(3)

   end subroutine copy_mom_to_scal

end module potential_mod
