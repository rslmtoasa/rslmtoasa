!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Radial ground-state contract
!
!> This module owns the small, immutable-after-capture radial data contract
!> needed by a future response consumer.  It deliberately has no dependency
!> on the SCF solver or on XC implementation details.
!------------------------------------------------------------------------------

module radial_ground_state_mod
   use precision_mod, only: rp
   implicit none
   private

   real(rp), parameter, public :: RADIAL_PI = 3.1415926535897932384626433832795_rp

   type, public :: radial_xc_provenance
      character(len=32) :: backend_name = ''
      character(len=3) :: txch = ''
      character(len=512) :: functional_name = ''
      character(len=32) :: mapping_quality = ''
      character(len=16) :: internal_energy_units = 'Ry'
      integer :: txc = 0
      integer :: libxc_family = -1
      integer :: libxc_nspin = -1
      integer :: n_components = 0
      logical :: use_libxc = .false.
      logical :: spin_polarized = .true.
      logical :: radial_gga = .false.
      integer, allocatable :: component_ids(:)
      integer, allocatable :: component_families(:)
      integer, allocatable :: component_kinds(:)
   contains
      procedure :: clear => clear_xc_provenance
   end type radial_xc_provenance

   type, public :: radial_ground_state
      ! Lifecycle marker: valid is set by capture; accepted is set only after
      ! ATOMSC has computed the residual for the final, un-mixed inner step.
      logical :: valid = .false.
      logical :: accepted = .false.
      integer :: accepted_inner_iteration = 0
      real(rp) :: accepted_residual_control = huge(1.0_rp)

      ! Mesh metadata.  r is in bohr and A/B reproduce dr/di=A*(r+B).
      real(rp) :: a = 0.0_rp
      real(rp) :: b = 0.0_rp
      real(rp) :: rmax = 0.0_rp
      character(len=16) :: radial_units = 'bohr'
      character(len=32) :: density_units = 'bohr^-3'
      character(len=16) :: energy_units = 'Ry'
      character(len=128) :: channel_convention = &
         'up=RHO(:,1)=Vxc(:,1); down=RHO(:,2)=Vxc(:,2)'

      ! The exact arrays needed by a future response consumer.
      real(rp), allocatable :: r(:)
      real(rp), allocatable :: rho_weighted_up(:), rho_weighted_down(:)
      real(rp), allocatable :: n_up(:), n_down(:)
      real(rp), allocatable :: vxc_up(:), vxc_down(:)

      ! Convention-neutral and explicitly named derived fields.
      real(rp), allocatable :: delta_vxc(:)
      real(rp), allocatable :: vxc_scalar(:)
      real(rp), allocatable :: bxc_pauli(:)
      real(rp), allocatable :: total_ks_spin_splitting(:)
      real(rp) :: constraining_field_ry = 0.0_rp

      type(radial_xc_provenance) :: xc_provenance

      ! Native RS-LMTO magnetic moment is the positive up-minus-down spin
      ! number in units of mu_B.  These are diagnostics, not inputs to the
      ! XC field definition.
      real(rp) :: integrated_n_up = 0.0_rp
      real(rp) :: integrated_n_down = 0.0_rp
      real(rp) :: integrated_spin_number = 0.0_rp
      real(rp) :: integrated_moment_muB = 0.0_rp
      real(rp) :: reported_moment(3) = 0.0_rp
   contains
      procedure :: clear
      procedure :: capture
      procedure :: mark_accepted
      procedure :: set_reported_moment
      procedure :: log_mesh_integral
      procedure :: charge_integral
      procedure :: spin_number_integral
      procedure :: write_file
      final :: finalize_radial_ground_state
   end type radial_ground_state

   public :: radial_simpson_weight

contains

   pure function radial_simpson_weight(ir, nr) result(weight)
      integer, intent(in) :: ir, nr
      real(rp) :: weight

      weight = 2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp
      if (ir == 1 .or. ir == nr) weight = 1.0_rp/3.0_rp
   end function radial_simpson_weight

   subroutine clear_xc_provenance(this)
      class(radial_xc_provenance), intent(inout) :: this

      if (allocated(this%component_ids)) deallocate(this%component_ids)
      if (allocated(this%component_families)) deallocate(this%component_families)
      if (allocated(this%component_kinds)) deallocate(this%component_kinds)
      this%backend_name = ''
      this%txch = ''
      this%functional_name = ''
      this%mapping_quality = ''
      this%internal_energy_units = 'Ry'
      this%txc = 0
      this%libxc_family = -1
      this%libxc_nspin = -1
      this%n_components = 0
      this%use_libxc = .false.
      this%spin_polarized = .true.
      this%radial_gga = .false.
   end subroutine clear_xc_provenance

   subroutine clear(this)
      class(radial_ground_state), intent(inout) :: this

      if (allocated(this%r)) deallocate(this%r)
      if (allocated(this%rho_weighted_up)) deallocate(this%rho_weighted_up)
      if (allocated(this%rho_weighted_down)) deallocate(this%rho_weighted_down)
      if (allocated(this%n_up)) deallocate(this%n_up)
      if (allocated(this%n_down)) deallocate(this%n_down)
      if (allocated(this%vxc_up)) deallocate(this%vxc_up)
      if (allocated(this%vxc_down)) deallocate(this%vxc_down)
      if (allocated(this%delta_vxc)) deallocate(this%delta_vxc)
      if (allocated(this%vxc_scalar)) deallocate(this%vxc_scalar)
      if (allocated(this%bxc_pauli)) deallocate(this%bxc_pauli)
      if (allocated(this%total_ks_spin_splitting)) deallocate(this%total_ks_spin_splitting)
      call this%xc_provenance%clear()
      this%valid = .false.
      this%accepted = .false.
      this%accepted_inner_iteration = 0
      this%accepted_residual_control = huge(1.0_rp)
      this%a = 0.0_rp
      this%b = 0.0_rp
      this%rmax = 0.0_rp
      this%constraining_field_ry = 0.0_rp
      this%integrated_n_up = 0.0_rp
      this%integrated_n_down = 0.0_rp
      this%integrated_spin_number = 0.0_rp
      this%integrated_moment_muB = 0.0_rp
      this%reported_moment = 0.0_rp
   end subroutine clear

   subroutine capture(this, a, b, rofi, rho, density_origin, vxc_up, vxc_down, total_v, b_fsm)
      class(radial_ground_state), intent(inout) :: this
      real(rp), intent(in) :: a, b
      real(rp), intent(in) :: rofi(:)
      real(rp), intent(in) :: rho(:, :)
      real(rp), intent(in) :: density_origin(2)
      real(rp), intent(in) :: vxc_up(:), vxc_down(:)
      real(rp), intent(in) :: total_v(:, :)
      real(rp), intent(in) :: b_fsm

      integer :: ir, nr
      real(rp) :: four_pi_r2

      nr = size(rofi)
      if (nr < 3 .or. size(rho, 1) /= nr .or. size(rho, 2) /= 2 .or. &
          size(vxc_up) /= nr .or. size(vxc_down) /= nr .or. &
          size(total_v, 1) /= nr .or. size(total_v, 2) /= 2) then
         error stop 'radial_ground_state%capture: inconsistent radial dimensions'
      end if

      call this%clear()
      allocate(this%r(nr), this%rho_weighted_up(nr), this%rho_weighted_down(nr), &
               this%n_up(nr), this%n_down(nr), this%vxc_up(nr), this%vxc_down(nr), &
               this%delta_vxc(nr), this%vxc_scalar(nr), this%bxc_pauli(nr), &
               this%total_ks_spin_splitting(nr))

      this%a = a
      this%b = b
      this%r = rofi
      this%rmax = rofi(nr)
      this%rho_weighted_up = rho(:, 1)
      this%rho_weighted_down = rho(:, 2)
      this%vxc_up = vxc_up
      this%vxc_down = vxc_down
      this%delta_vxc = this%vxc_up - this%vxc_down
      this%vxc_scalar = 0.5_rp*(this%vxc_up + this%vxc_down)
      this%bxc_pauli = 0.5_rp*this%delta_vxc
      this%total_ks_spin_splitting = total_v(:, 1) - total_v(:, 2)
      this%constraining_field_ry = b_fsm
      ! Mark the arrays usable before evaluating the derived integrals.
      this%valid = .true.

      ! The first mesh point is r=0.  VXC0SP obtains the regular physical
      ! density there by extrapolation; every other point follows exactly
      ! n_sigma=RHO_sigma/(4*pi*r^2).
      this%n_up(1) = density_origin(1)
      this%n_down(1) = density_origin(2)
      do ir = 2, nr
         four_pi_r2 = 4.0_rp*RADIAL_PI*this%r(ir)**2
         this%n_up(ir) = this%rho_weighted_up(ir)/four_pi_r2
         this%n_down(ir) = this%rho_weighted_down(ir)/four_pi_r2
      end do

      this%integrated_n_up = this%charge_integral(1)
      this%integrated_n_down = this%charge_integral(2)
      this%integrated_spin_number = this%integrated_n_up - this%integrated_n_down
      this%integrated_moment_muB = this%integrated_spin_number
   end subroutine capture

   subroutine mark_accepted(this, inner_iteration, residual_control)
      class(radial_ground_state), intent(inout) :: this
      integer, intent(in) :: inner_iteration
      real(rp), intent(in) :: residual_control

      if (.not. this%valid) error stop 'radial_ground_state%mark_accepted: no captured state'
      this%accepted = .true.
      this%accepted_inner_iteration = inner_iteration
      this%accepted_residual_control = residual_control
   end subroutine mark_accepted

   subroutine set_reported_moment(this, moment)
      class(radial_ground_state), intent(inout) :: this
      real(rp), intent(in) :: moment(3)

      this%reported_moment = moment
   end subroutine set_reported_moment

   function log_mesh_integral(this, values) result(integral)
      class(radial_ground_state), intent(in) :: this
      real(rp), intent(in) :: values(:)
      real(rp) :: integral
      integer :: ir, nr

      if (.not. this%valid .or. size(values) /= size(this%r)) then
         error stop 'radial_ground_state%log_mesh_integral: invalid state or size'
      end if
      nr = size(this%r)
      integral = 0.0_rp
      do ir = 1, nr
         integral = integral + radial_simpson_weight(ir, nr)*this%a*(this%r(ir) + this%b)*values(ir)
      end do
   end function log_mesh_integral

   function charge_integral(this, channel) result(integral)
      class(radial_ground_state), intent(in) :: this
      integer, intent(in) :: channel
      real(rp) :: integral

      if (channel == 1) then
         integral = this%log_mesh_integral(this%rho_weighted_up)
      else if (channel == 2) then
         integral = this%log_mesh_integral(this%rho_weighted_down)
      else
         error stop 'radial_ground_state%charge_integral: channel must be 1 or 2'
      end if
   end function charge_integral

   function spin_number_integral(this) result(integral)
      class(radial_ground_state), intent(in) :: this
      real(rp) :: integral

      integral = this%charge_integral(1) - this%charge_integral(2)
   end function spin_number_integral

   subroutine write_file(this, filename, site_index)
      class(radial_ground_state), intent(in) :: this
      character(len=*), intent(in) :: filename
      integer, intent(in), optional :: site_index

      integer :: unit, ir, site
      character(len=16) :: logical_text

      site = 0
      if (present(site_index)) site = site_index
      open(newunit=unit, file=filename, status='replace', action='write')
      write(unit, '(a)') '# LR-01 converged radial ground-state snapshot'
      write(unit, '(a,i0)') '# site = ', site
      if (this%valid) then
         logical_text = 'true'
      else
         logical_text = 'false'
      end if
      write(unit, '(a,a)') '# valid = ', trim(logical_text)
      if (this%accepted) then
         logical_text = 'true'
      else
         logical_text = 'false'
      end if
      write(unit, '(a,a)') '# accepted = ', trim(logical_text)
      write(unit, '(a,i0)') '# accepted_inner_iteration = ', this%accepted_inner_iteration
      write(unit, '(a,es24.16)') '# accepted_residual_control = ', this%accepted_residual_control
      write(unit, '(a,2(es24.16,1x))') '# mesh_a mesh_b = ', this%a, this%b
      write(unit, '(a,es24.16)') '# rmax_bohr = ', this%rmax
      write(unit, '(a)') '# radial_units = bohr; density_units = bohr^-3; energy_units = Ry'
      write(unit, '(a,a)') '# channel_convention = ', trim(this%channel_convention)
      write(unit, '(a,es24.16)') '# constraining_field_ry = ', this%constraining_field_ry
      write(unit, '(a,a)') '# xc_backend = ', trim(this%xc_provenance%backend_name)
      write(unit, '(a,i0)') '# xc_txc = ', this%xc_provenance%txc
      write(unit, '(a,a)') '# xc_txch = ', trim(this%xc_provenance%txch)
      write(unit, '(a,a)') '# xc_functional = ', trim(this%xc_provenance%functional_name)
      write(unit, '(a,a)') '# xc_mapping_quality = ', trim(this%xc_provenance%mapping_quality)
      write(unit, '(a,l1)') '# xc_use_libxc = ', this%xc_provenance%use_libxc
      write(unit, '(a,l1)') '# xc_spin_polarized = ', this%xc_provenance%spin_polarized
      write(unit, '(a,l1)') '# xc_radial_gga = ', this%xc_provenance%radial_gga
      write(unit, '(a,i0)') '# xc_libxc_family = ', this%xc_provenance%libxc_family
      write(unit, '(a,i0)') '# xc_libxc_nspin = ', this%xc_provenance%libxc_nspin
      write(unit, '(a,i0)') '# xc_component_count = ', this%xc_provenance%n_components
      if (this%xc_provenance%n_components > 0) then
         write(unit, '(a,*(i0,1x))') '# xc_component_ids = ', this%xc_provenance%component_ids
         write(unit, '(a,*(i0,1x))') '# xc_component_families = ', this%xc_provenance%component_families
         write(unit, '(a,*(i0,1x))') '# xc_component_kinds = ', this%xc_provenance%component_kinds
      end if
      write(unit, '(a,3(es24.16,1x))') '# reported_moment_xyz_muB = ', this%reported_moment
      write(unit, '(a,es24.16)') '# integrated_n_up = ', this%integrated_n_up
      write(unit, '(a,es24.16)') '# integrated_n_down = ', this%integrated_n_down
      write(unit, '(a,es24.16)') '# integrated_spin_number = ', this%integrated_spin_number
      write(unit, '(a,es24.16)') '# integrated_moment_muB = ', this%integrated_moment_muB
      write(unit, '(a)') '# columns: i r rho_weighted_up rho_weighted_down n_up n_down vxc_up vxc_down '// &
                         'delta_vxc vxc_scalar bxc_pauli total_ks_spin_splitting'
      if (this%valid) then
         do ir = 1, size(this%r)
            write(unit, '(i8,1x,11(es24.16,1x))') ir, this%r(ir), this%rho_weighted_up(ir), &
               this%rho_weighted_down(ir), this%n_up(ir), this%n_down(ir), this%vxc_up(ir), &
               this%vxc_down(ir), this%delta_vxc(ir), this%vxc_scalar(ir), this%bxc_pauli(ir), &
               this%total_ks_spin_splitting(ir)
         end do
      end if
      close(unit)
   end subroutine write_file

   subroutine finalize_radial_ground_state(this)
      type(radial_ground_state) :: this

      call this%clear()
   end subroutine finalize_radial_ground_state

end module radial_ground_state_mod
