!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Input boundary for the LR-TDDFT post-processing workflow.
!>
!> The configuration is deliberately separate from &reciprocal: reciprocal
!> owns the normal-state k integration, whereas this object owns response q,
!> frequency, projection, and output controls.  q points are always stored in
!> fractional reciprocal-lattice coordinates, the convention consumed by the
!> arbitrary-k eigenpair service.
module tddft_config_mod
   use precision_mod, only: rp
   use string_mod, only: lower, sl, int2str
   use logger_mod, only: g_logger
   implicit none

   private

   type, public :: tddft_config
      character(len=sl) :: fname
      logical :: enabled
      character(len=16) :: channel
      character(len=16) :: chi0_backend
      !> Self-enhancement construction.  The site-scalar route is retained as
      !> the backwards-compatible diagnostic baseline; pair_potential is the
      !> validated ham_only LMTO tangent route and compare writes both raw
      !> results from the same response inputs.
      character(len=24) :: xi_backend
      character(len=16) :: response_projection
      character(len=16) :: q_mode
      character(len=16) :: q_coordinates
      character(len=16) :: goldstone_mode
      ! TDDFT-02 policy: diagnose is the only implicit/default action;
      ! sum_rule and projected are explicit, auditable repairs.
      character(len=16) :: goldstone_policy
      logical :: goldstone_mode_migrated_from_sum_rule
      character(len=sl) :: output_prefix
      integer :: nomega
      real(rp) :: omega_min, omega_max, eta
      ! The response Fermi level is resolved internally on the complete
      ! response mesh from the ground-state target electron count.  Temperature
      ! may be explicitly overridden; all resolved values are emitted with
      ! every response product by the production driver.
      real(rp) :: electronic_temperature, fermi_level
      real(rp) :: ground_state_electronic_temperature, ground_state_fermi_level, ground_state_electron_count
      real(rp) :: response_electron_count
      logical :: electronic_temperature_overridden
      integer :: band_first, band_last
      real(rp) :: occupation_tolerance
      !> Real-axis GF bubble controls.  A zero green_eta means eta/2, and a
      !> reversed energy window requests a source-spectrum-derived window.
      real(rp) :: green_eta, green_energy_min, green_energy_max
      integer :: green_energy_points
      real(rp) :: realspace_rmax, realspace_tail_tolerance
      integer :: realspace_fourier_axes(3)
      character(len=16) :: realspace_representation
      character(len=sl) :: longitudinal_static_file
      real(rp) :: longitudinal_pair_tolerance, longitudinal_linearity_tolerance
      real(rp) :: longitudinal_static_agreement_tolerance, longitudinal_fit_omega_min, longitudinal_fit_omega_max
      logical :: output_chi0, output_xi, output_chi, output_modes, output_stoner
      real(rp), allocatable :: q_points(:, :)
   contains
      procedure :: restore_to_default
      procedure :: build_from_file
   end type tddft_config

   interface tddft_config
      procedure :: constructor
   end interface tddft_config

contains

   function constructor(fname) result(obj)
      type(tddft_config) :: obj
      character(len=*), intent(in) :: fname

      call obj%restore_to_default()
      call obj%build_from_file(fname)
   end function constructor

   subroutine restore_to_default(this)
      class(tddft_config), intent(inout) :: this

      this%enabled = .true.
      this%channel = 'transverse'
      this%chi0_backend = 'eigenpairs'
      this%xi_backend = 'legacy_site_scalar'
      this%response_projection = 'site'
      this%q_mode = 'list'
      this%q_coordinates = 'direct'
      this%goldstone_mode = 'diagnose'
      this%goldstone_policy = 'diagnose'
      this%goldstone_mode_migrated_from_sum_rule = .false.
      this%output_prefix = 'tddft'
      this%nomega = 201
      this%omega_min = 0.0_rp
      this%omega_max = 0.5_rp
      this%eta = 0.01_rp
      this%electronic_temperature = -1.0_rp
      this%fermi_level = huge(1.0_rp)
      this%ground_state_electronic_temperature = -1.0_rp
      this%ground_state_fermi_level = huge(1.0_rp)
      this%ground_state_electron_count = -1.0_rp
      this%response_electron_count = -1.0_rp
      this%electronic_temperature_overridden = .false.
      this%band_first = 1
      this%band_last = 0
      this%occupation_tolerance = 0.0_rp
      this%green_eta = 0.0_rp
      this%green_energy_min = huge(1.0_rp)
      this%green_energy_max = -huge(1.0_rp)
      this%green_energy_points = 2001
      this%realspace_rmax = huge(1.0_rp)
      this%realspace_tail_tolerance = 1.0e-3_rp
      this%realspace_fourier_axes = [1, 2, 3]
      this%realspace_representation = 'bulk'
      this%longitudinal_static_file = ''
      this%longitudinal_pair_tolerance = 1.0e-10_rp
      this%longitudinal_linearity_tolerance = 5.0e-2_rp
      this%longitudinal_static_agreement_tolerance = 5.0e-2_rp
      this%longitudinal_fit_omega_min = 0.0_rp
      this%longitudinal_fit_omega_max = huge(1.0_rp)
      this%output_chi0 = .true.
      this%output_xi = .false.
      this%output_chi = .false.
      this%output_modes = .false.
      this%output_stoner = .true.
      if (allocated(this%q_points)) deallocate(this%q_points)
      allocate(this%q_points(3, 1))
      this%q_points = 0.0_rp
   end subroutine restore_to_default

   subroutine build_from_file(this, fname)
      class(tddft_config), intent(inout) :: this
      character(len=*), intent(in) :: fname
      integer :: funit, iostatus, iq, ix, iy, iz, nmesh

      include 'include_codes/namelists/tddft.f90'

      enabled = this%enabled
      channel = this%channel
      chi0_backend = this%chi0_backend
      xi_backend = this%xi_backend
      response_projection = this%response_projection
      q_mode = this%q_mode
      q_coordinates = this%q_coordinates
      goldstone_mode = this%goldstone_mode
      goldstone_policy = this%goldstone_policy
      output_prefix = this%output_prefix
      q_file = ''
      n_q_points = 1
      q_mesh1 = 1; q_mesh2 = 1; q_mesh3 = 1
      q_start = 0.0_rp; q_end = 0.0_rp; q_list = 0.0_rp
      nomega = this%nomega
      omega_min = this%omega_min; omega_max = this%omega_max; eta = this%eta
      ! A negative temperature means inherit the reciprocal ground-state
      ! value.  The Fermi level is not an input: response occupations are
      ! always resolved from the target electron count on the response mesh.
      electronic_temperature = -1.0_rp
      band_first = this%band_first; band_last = this%band_last
      occupation_tolerance = this%occupation_tolerance
      green_eta = this%green_eta
      green_energy_min = this%green_energy_min; green_energy_max = this%green_energy_max
      green_energy_points = this%green_energy_points
      realspace_rmax = this%realspace_rmax
      realspace_tail_tolerance = this%realspace_tail_tolerance
      realspace_fourier_axes = this%realspace_fourier_axes
      realspace_representation = this%realspace_representation
      longitudinal_static_file = this%longitudinal_static_file
      longitudinal_pair_tolerance = this%longitudinal_pair_tolerance
      longitudinal_linearity_tolerance = this%longitudinal_linearity_tolerance
      longitudinal_static_agreement_tolerance = this%longitudinal_static_agreement_tolerance
      longitudinal_fit_omega_min = this%longitudinal_fit_omega_min
      longitudinal_fit_omega_max = this%longitudinal_fit_omega_max
      output_chi0 = this%output_chi0; output_xi = this%output_xi; output_chi = this%output_chi
      output_modes = this%output_modes; output_stoner = this%output_stoner

      open(newunit=funit, file=fname, action='read', status='old', iostat=iostatus)
      if (iostatus /= 0) call g_logger%fatal('[tddft_config.build_from_file]: input file not found', __FILE__, __LINE__)
      read(funit, nml=tddft, iostat=iostatus)
      close(funit)
      ! An absent &tddft group preserves the useful gamma/default response when
      ! the explicitly selected calculation route is used.
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%fatal('[tddft_config.build_from_file]: error reading &tddft, iostatus='//int2str(iostatus), &
            __FILE__, __LINE__)
      end if

      this%fname = fname
      this%enabled = enabled
      this%channel = lower(trim(channel))
      this%chi0_backend = lower(trim(chi0_backend))
      this%xi_backend = lower(trim(xi_backend))
      this%response_projection = lower(trim(response_projection))
      this%q_mode = lower(trim(q_mode))
      this%q_coordinates = lower(trim(q_coordinates))
      this%goldstone_mode = lower(trim(goldstone_mode))
      this%goldstone_policy = lower(trim(goldstone_policy))
      this%goldstone_mode_migrated_from_sum_rule = this%goldstone_mode == 'sum_rule'
      if (this%goldstone_mode_migrated_from_sum_rule) then
         this%goldstone_mode = 'correct'
         call g_logger%warning("[tddft_config]: goldstone_mode='sum_rule' is deprecated; using explicit pair-potential 'correct'.", &
            __FILE__, __LINE__)
      end if
      this%output_prefix = trim(output_prefix)
      this%nomega = nomega
      this%omega_min = omega_min; this%omega_max = omega_max; this%eta = eta
      this%electronic_temperature = electronic_temperature
      this%electronic_temperature_overridden = electronic_temperature >= 0.0_rp
      this%band_first = band_first; this%band_last = band_last
      this%occupation_tolerance = occupation_tolerance
      this%green_eta = green_eta
      this%green_energy_min = green_energy_min; this%green_energy_max = green_energy_max
      this%green_energy_points = green_energy_points
      this%realspace_rmax = realspace_rmax
      this%realspace_tail_tolerance = realspace_tail_tolerance
      this%realspace_fourier_axes = realspace_fourier_axes
      this%realspace_representation = lower(trim(realspace_representation))
      this%longitudinal_static_file = trim(longitudinal_static_file)
      this%longitudinal_pair_tolerance = longitudinal_pair_tolerance
      this%longitudinal_linearity_tolerance = longitudinal_linearity_tolerance
      this%longitudinal_static_agreement_tolerance = longitudinal_static_agreement_tolerance
      this%longitudinal_fit_omega_min = longitudinal_fit_omega_min
      this%longitudinal_fit_omega_max = longitudinal_fit_omega_max
      this%output_chi0 = output_chi0; this%output_xi = output_xi; this%output_chi = output_chi
      this%output_modes = output_modes; this%output_stoner = output_stoner

      call validate_scalar_settings(this)
      select case (this%q_mode)
      case ('list')
         if (len_trim(q_file) > 0) then
            call read_q_file(this, trim(q_file))
         else
            if (n_q_points < 1 .or. n_q_points > tddft_max_q) then
               call g_logger%fatal('[tddft_config.build_from_file]: list n_q_points must be 1..'//int2str(tddft_max_q), &
                  __FILE__, __LINE__)
            end if
            call replace_q_points(this, q_list(:, 1:n_q_points))
         end if
      case ('path')
         if (len_trim(q_file) > 0) then
            call read_q_file(this, trim(q_file))
         else
            if (n_q_points < 2 .or. n_q_points > tddft_max_q) then
               call g_logger%fatal('[tddft_config.build_from_file]: path n_q_points must be 2..'//int2str(tddft_max_q), &
                  __FILE__, __LINE__)
            end if
            if (allocated(this%q_points)) deallocate(this%q_points)
            allocate(this%q_points(3, n_q_points))
            do iq = 1, n_q_points
               this%q_points(:, iq) = q_start + real(iq-1, rp)*(q_end-q_start)/real(n_q_points-1, rp)
            end do
         end if
      case ('mesh')
         if (q_mesh1 < 1 .or. q_mesh2 < 1 .or. q_mesh3 < 1) then
            call g_logger%fatal('[tddft_config.build_from_file]: q_mesh1/2/3 must be positive', __FILE__, __LINE__)
         end if
         nmesh = q_mesh1*q_mesh2*q_mesh3
         if (allocated(this%q_points)) deallocate(this%q_points)
         allocate(this%q_points(3, nmesh))
         iq = 0
         do iz = 1, q_mesh3
            do iy = 1, q_mesh2
               do ix = 1, q_mesh1
                  iq = iq + 1
                  this%q_points(:, iq) = [real(2*ix-q_mesh1-1, rp)/(2.0_rp*real(q_mesh1, rp)), &
                     real(2*iy-q_mesh2-1, rp)/(2.0_rp*real(q_mesh2, rp)), &
                     real(2*iz-q_mesh3-1, rp)/(2.0_rp*real(q_mesh3, rp))]
               end do
            end do
         end do
      case default
         call g_logger%fatal("[tddft_config.build_from_file]: q_mode must be 'list', 'path', or 'mesh'", __FILE__, __LINE__)
      end select
   end subroutine build_from_file

   subroutine validate_scalar_settings(this)
      class(tddft_config), intent(in) :: this
      integer :: ix

      if (this%channel /= 'transverse' .and. this%channel /= 'longitudinal' .and. this%channel /= 'full') then
         call g_logger%fatal("[tddft_config]: channel must be 'transverse', 'longitudinal', or 'full'", __FILE__, __LINE__)
      end if
      if (this%chi0_backend /= 'eigenpairs' .and. this%chi0_backend /= 'green' .and. &
          this%chi0_backend /= 'kspace_lehmann' .and. this%chi0_backend /= 'kspace_gf' .and. &
          this%chi0_backend /= 'realspace_gf' .and. this%chi0_backend /= 'realspace' .and. &
          this%chi0_backend /= 'rs_gf' .and. this%chi0_backend /= 'real_space_gf') then
         call g_logger%fatal("[tddft_config]: chi0_backend must be 'eigenpairs', 'kspace_lehmann', or 'realspace_gf' "// &
            "(compatibility aliases: 'green', 'kspace_gf', 'realspace', 'rs_gf')", __FILE__, __LINE__)
      end if
      if (this%xi_backend /= 'legacy_site_scalar' .and. this%xi_backend /= 'pair_potential' .and. &
          this%xi_backend /= 'compare') then
         call g_logger%fatal("[tddft_config]: xi_backend must be 'legacy_site_scalar', 'pair_potential', or 'compare'", &
            __FILE__, __LINE__)
      end if
      if (this%response_projection /= 'site') then
         call g_logger%fatal("[tddft_config]: only response_projection='site' is currently implemented", __FILE__, __LINE__)
      end if
      if (this%q_coordinates /= 'direct') then
         call g_logger%fatal("[tddft_config]: q_coordinates must be 'direct' (fractional reciprocal coordinates)", &
            __FILE__, __LINE__)
      end if
      if (this%goldstone_mode /= 'off' .and. this%goldstone_mode /= 'diagnose' .and. this%goldstone_mode /= 'correct') then
         call g_logger%fatal("[tddft_config]: goldstone_mode must be 'off', 'diagnose', or 'correct' (deprecated 'sum_rule')", __FILE__, __LINE__)
      end if
      if (this%goldstone_policy /= 'diagnose' .and. this%goldstone_policy /= 'sum_rule' .and. &
          this%goldstone_policy /= 'projected') then
         call g_logger%fatal("[tddft_config]: goldstone_policy must be 'diagnose', 'sum_rule', or 'projected'", __FILE__, __LINE__)
      end if
      if (this%goldstone_mode == 'correct') then
         if (this%xi_backend /= 'pair_potential' .and. this%xi_backend /= 'compare') then
            call g_logger%fatal("[tddft_config]: goldstone_mode='correct' cannot be used with xi_backend='legacy_site_scalar'. "// &
               "Use goldstone_mode='diagnose' for the legacy raw spectrum, or select xi_backend='pair_potential' or 'compare' for a controlled correction.", &
               __FILE__, __LINE__)
         end if
         if (.not. (this%output_xi .or. this%output_chi)) then
            call g_logger%fatal("[tddft_config]: goldstone_mode='correct' requires output_xi=.true. or output_chi=.true. "// &
               'so the raw pair-potential Xi remains auditable.', __FILE__, __LINE__)
         end if
      end if
      if (this%nomega < 1 .or. this%omega_max < this%omega_min .or. this%eta <= 0.0_rp .or. &
          (this%electronic_temperature_overridden .and. this%electronic_temperature < 0.0_rp) .or. &
          this%band_first < 1 .or. this%band_last < 0 .or. &
          this%occupation_tolerance < 0.0_rp .or. len_trim(this%output_prefix) == 0) then
         call g_logger%fatal('[tddft_config]: invalid frequency, band, temperature, or output settings', __FILE__, __LINE__)
      end if
      if (this%green_eta < 0.0_rp .or. this%green_energy_points < 3 .or. &
          (this%green_energy_min < huge(1.0_rp)/2.0_rp .and. &
           this%green_energy_max > -huge(1.0_rp)/2.0_rp .and. this%green_energy_max <= this%green_energy_min)) then
         call g_logger%fatal('[tddft_config]: invalid Green-function energy integration settings', __FILE__, __LINE__)
      end if
      if (this%realspace_rmax <= 0.0_rp .or. this%realspace_tail_tolerance < 0.0_rp .or. &
          (this%realspace_representation /= 'bulk' .and. this%realspace_representation /= 'film' .and. &
           this%realspace_representation /= 'finite' .and. this%realspace_representation /= 'local') .or. &
          any(this%realspace_fourier_axes < 0) .or. any(this%realspace_fourier_axes > 3)) then
         call g_logger%fatal('[tddft_config]: invalid native real-space Green-function settings', __FILE__, __LINE__)
      end if
      if (this%realspace_representation == 'bulk' .and. any(this%realspace_fourier_axes == 0)) then
         call g_logger%fatal('[tddft_config]: bulk real-space representation requires all three Fourier axes', __FILE__, __LINE__)
      end if
      if (this%realspace_representation == 'film') then
         if (count(this%realspace_fourier_axes > 0) /= 2 .or. &
             any([(count(this%realspace_fourier_axes == ix), ix=1, 3)] > 1)) then
            call g_logger%fatal('[tddft_config]: film real-space representation requires two Fourier axes', __FILE__, __LINE__)
         end if
      end if
      if (this%channel == 'longitudinal') then
         if (len_trim(this%longitudinal_static_file) == 0 .or. this%longitudinal_pair_tolerance <= 0.0_rp .or. &
             this%longitudinal_linearity_tolerance < 0.0_rp .or. this%longitudinal_static_agreement_tolerance < 0.0_rp .or. &
             this%longitudinal_fit_omega_max < this%longitudinal_fit_omega_min) then
            call g_logger%fatal('[tddft_config]: longitudinal response requires a static +/- field file and valid fit tolerances', &
               __FILE__, __LINE__)
         end if
      end if
      if (this%output_modes .and. (this%nomega < 2 .or. this%omega_max <= this%omega_min)) then
         call g_logger%fatal('[tddft_config]: output_modes requires at least two strictly increasing frequencies', &
            __FILE__, __LINE__)
      end if
   end subroutine validate_scalar_settings

   subroutine read_q_file(this, filename)
      class(tddft_config), intent(inout) :: this
      character(len=*), intent(in) :: filename
      integer :: unit, ios, nq, iq
      real(rp), allocatable :: points(:, :)

      open(newunit=unit, file=filename, action='read', status='old', iostat=ios)
      if (ios /= 0) call g_logger%fatal('[tddft_config.read_q_file]: cannot open q_file '//trim(filename), __FILE__, __LINE__)
      read(unit, *, iostat=ios) nq
      if (ios /= 0 .or. nq < 1) then
         close(unit)
         call g_logger%fatal('[tddft_config.read_q_file]: first record must be a positive q-point count', __FILE__, __LINE__)
      end if
      allocate(points(3, nq))
      do iq = 1, nq
         read(unit, *, iostat=ios) points(:, iq)
         if (ios /= 0) then
            close(unit)
            call g_logger%fatal('[tddft_config.read_q_file]: unable to read q point '//int2str(iq), __FILE__, __LINE__)
         end if
      end do
      close(unit)
      call replace_q_points(this, points)
   end subroutine read_q_file

   subroutine replace_q_points(this, points)
      class(tddft_config), intent(inout) :: this
      real(rp), intent(in) :: points(:, :)

      if (size(points, 1) /= 3 .or. size(points, 2) < 1) error stop 'replace_q_points: invalid q-point array'
      if (allocated(this%q_points)) deallocate(this%q_points)
      allocate(this%q_points(3, size(points, 2)))
      this%q_points = points
   end subroutine replace_q_points

end module tddft_config_mod
