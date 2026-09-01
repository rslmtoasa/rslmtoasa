!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Ground-state XC data contract for future LR-TDDFT kernels.
!>
!> This module intentionally stores the radial spin-channel XC potentials that
!> the SCF functional returns, rather than a block extracted from the assembled
!> LMTO Hamiltonian.  A future site projection may choose its own radial weights
!> but must retain that provenance.
module xc_response_kernel_mod
   use precision_mod, only: rp
   use tddft_conventions_mod, only: tddft_circular_operator_factor, tddft_circular_source_factor
   use xc_mod, only: xc
   implicit none

   private

   type, public :: xc_response_sample
      ! Input densities use the SCF ordering: down first, up second.
      real(rp) :: rho_down = 0.0_rp
      real(rp) :: rho_up = 0.0_rp
      real(rp) :: vxc_down = 0.0_rp
      real(rp) :: vxc_up = 0.0_rp
      real(rp) :: exc = 0.0_rp
      ! Energy coefficients in V_xc = vxc_scalar sigma_0 + bxc_energy sigma_z.
      real(rp) :: vxc_scalar = 0.0_rp
      real(rp) :: bxc_energy = 0.0_rp
   end type xc_response_sample

   type, public :: xc_response_site
      ! Site-projected spin population n_up - n_down; it is not a Tesla field
      ! or a magnetic moment in SI units.
      real(rp) :: spin_population = 0.0_rp
      ! Signed collinear P_site sigma_z population.  It is separate from the
      ! legacy magnitude-oriented field above and is the Q_a normalization.
      real(rp) :: signed_spin_population = 0.0_rp
      ! Radial ASA spin population used to define the projected ALSDA
      ! stiffness.  It is retained separately from spin_population, whose
      ! value is the response-projector population (P_site sigma).
      real(rp) :: radial_spin_population = 0.0_rp
      ! Numerator of the radial ALSDA projection, retained until the
      ! P_site-sigma population has been supplied.
      real(rp) :: bxc_spin_moment = 0.0_rp
      ! Additional ordinary collinear diagnostics.  These are radial
      ! quadratures, not response-kernel parameters.
      real(rp) :: radial_spin_abs_population = 0.0_rp
      real(rp) :: radial_vxc_spin_difference = 0.0_rp
      real(rp) :: radial_vxc_spin_difference_abs = 0.0_rp
      real(rp) :: vxc_scalar = 0.0_rp
      real(rp) :: bxc_energy = 0.0_rp
      real(rp) :: dvxc_dn = 0.0_rp
      real(rp) :: dvxc_dm = 0.0_rp
      real(rp) :: dbxc_dn = 0.0_rp
      real(rp) :: dbxc_dm = 0.0_rp
      ! Circular (m+/- = mx +/- i my, sigma+/- = sigma_x +/- i sigma_y)
      ! transverse ALSDA scalar.  The unhalved circular operators require the
      ! explicit 1/2 in Bxc/(2 M).  It must not be used in Cartesian x/y.
      real(rp) :: k_perp_circular = 0.0_rp
      ! Unit vector giving the local ALSDA longitudinal axis in the global
      ! response frame.  It is explicitly retained instead of assuming z.
      real(rp) :: magnetization_direction(3) = [0.0_rp, 0.0_rp, 1.0_rp]
      logical :: has_dvxc_dn = .false.
      logical :: has_dvxc_dm = .false.
      logical :: has_dbxc_dn = .false.
      logical :: has_dbxc_dm = .false.
      logical :: has_k_perp_circular = .false.
      logical :: has_radial_projection = .false.
      logical :: has_magnetization_direction = .false.
   end type xc_response_site

   !> Accumulator for one radial VXC0SP evaluation.  `radial_weight` below is
   !> the SCF quadrature factor WGT*DRDI; its density arguments are the radial
   !> densities RHO used by the existing SCF integration.
   type, public :: xc_response_radial_projection
      real(rp) :: charge_population = 0.0_rp
      real(rp) :: spin_population = 0.0_rp
      real(rp) :: vxc_charge_moment = 0.0_rp
      real(rp) :: bxc_spin_moment = 0.0_rp
      real(rp) :: spin_abs_population = 0.0_rp
      real(rp) :: vxc_spin_difference = 0.0_rp
      real(rp) :: vxc_spin_difference_abs = 0.0_rp
   contains
      procedure :: clear => xc_radial_projection_clear
      procedure :: accumulate => xc_radial_projection_accumulate
   end type xc_response_radial_projection

   type, public :: xc_response_kernel_provider
      character(len=32) :: functional_label = 'unrecorded'
      ! Set only by a selected XC route after its full ALSDA derivative
      ! evaluator has passed its own derivative verification.  Merely filling
      ! synthetic slots is not production capability.
      logical :: full_alsda_derivatives_validated = .false.
      type(xc_response_site), allocatable :: site(:)
   contains
      procedure :: initialize => xc_kernel_initialize
      procedure :: clear => xc_kernel_clear
      procedure :: record_ground_state_site => xc_kernel_record_ground_state_site
      procedure :: record_radial_projection => xc_kernel_record_radial_projection
      procedure :: set_site_spin_population => xc_kernel_set_site_spin_population
      procedure :: set_site_signed_spin_population => xc_kernel_set_site_signed_spin_population
      procedure :: set_site_magnetization_direction => xc_kernel_set_site_magnetization_direction
      procedure :: set_site_derivatives => xc_kernel_set_site_derivatives
      procedure :: full_response_capability => xc_kernel_full_response_capability
   end type xc_response_kernel_provider

   public :: evaluate_ground_state_xc_sample
   public :: circular_transverse_kernel
   public :: cartesian_transverse_kernel

contains

   !> Evaluate the pointwise XC path called by self%VXC0SP.
   !> No constraining field is included here: B_fsm is added by VXC0SP after
   !> this functional call and is not an XC field.  A libXC GGA request is
   !> rejected by the pointwise wrapper; it must be evaluated by the complete
   !> radial helper instead.
   function evaluate_ground_state_xc_sample(functional, rho_down, rho_up, rho_total, &
      rho_gradient, rho_laplacian, radius) result(sample)
      class(xc), intent(in) :: functional
      real(rp), intent(in) :: rho_down, rho_up, rho_total, radius
      real(rp), intent(in) :: rho_gradient(2), rho_laplacian(2)
      type(xc_response_sample) :: sample
      real(rp) :: v_down, v_up

      call functional%XCPOT_hybrid(rho_down, rho_up, rho_total, rho_gradient, rho_laplacian, radius, &
         v_down, v_up, sample%exc)
      sample%rho_down = rho_down
      sample%rho_up = rho_up
      sample%vxc_down = v_down
      sample%vxc_up = v_up
      sample%vxc_scalar = 0.5_rp*(v_up + v_down)
      sample%bxc_energy = 0.5_rp*(v_up - v_down)
   end function evaluate_ground_state_xc_sample

   subroutine xc_radial_projection_clear(this)
      class(xc_response_radial_projection), intent(inout) :: this
      this%charge_population = 0.0_rp
      this%spin_population = 0.0_rp
      this%vxc_charge_moment = 0.0_rp
      this%bxc_spin_moment = 0.0_rp
      this%spin_abs_population = 0.0_rp
      this%vxc_spin_difference = 0.0_rp
      this%vxc_spin_difference_abs = 0.0_rp
   end subroutine xc_radial_projection_clear

   !> Accumulate an XC sample already evaluated by VXC0SP.  No second XC
   !> evaluation and no FSM/constraining field enter this record.
   subroutine xc_radial_projection_accumulate(this, radial_weight, rho_down, rho_up, vxc_down, vxc_up)
      class(xc_response_radial_projection), intent(inout) :: this
      real(rp), intent(in) :: radial_weight, rho_down, rho_up, vxc_down, vxc_up
      real(rp) :: charge_density, spin_density, vxc_scalar, bxc_energy

      charge_density = rho_up + rho_down
      spin_density = rho_up - rho_down
      vxc_scalar = 0.5_rp*(vxc_up + vxc_down)
      bxc_energy = 0.5_rp*(vxc_up - vxc_down)
      this%charge_population = this%charge_population + radial_weight*charge_density
      this%spin_population = this%spin_population + radial_weight*spin_density
      this%vxc_charge_moment = this%vxc_charge_moment + radial_weight*charge_density*vxc_scalar
      this%bxc_spin_moment = this%bxc_spin_moment + radial_weight*spin_density*bxc_energy
      this%spin_abs_population = this%spin_abs_population + radial_weight*abs(spin_density)
      this%vxc_spin_difference = this%vxc_spin_difference + radial_weight*(vxc_up - vxc_down)
      this%vxc_spin_difference_abs = this%vxc_spin_difference_abs + radial_weight*abs(vxc_up - vxc_down)
   end subroutine xc_radial_projection_accumulate

   subroutine xc_kernel_initialize(this, nsite, functional_label)
      class(xc_response_kernel_provider), intent(inout) :: this
      integer, intent(in) :: nsite
      character(len=*), intent(in), optional :: functional_label

      if (nsite < 0) error stop 'xc_response_kernel_provider%initialize: negative nsite'
      call this%clear()
      allocate(this%site(nsite))
      if (present(functional_label)) this%functional_label = functional_label
   end subroutine xc_kernel_initialize

   subroutine xc_kernel_clear(this)
      class(xc_response_kernel_provider), intent(inout) :: this
      if (allocated(this%site)) deallocate(this%site)
      this%functional_label = 'unrecorded'
      this%full_alsda_derivatives_validated = .false.
   end subroutine xc_kernel_clear

   subroutine xc_kernel_record_ground_state_site(this, isite, spin_population, sample)
      class(xc_response_kernel_provider), intent(inout) :: this
      integer, intent(in) :: isite
      real(rp), intent(in) :: spin_population
      type(xc_response_sample), intent(in) :: sample

      call require_site(this, isite, 'record_ground_state_site')
      this%site(isite)%spin_population = spin_population
      this%site(isite)%vxc_scalar = sample%vxc_scalar
      this%site(isite)%bxc_energy = sample%bxc_energy
   end subroutine xc_kernel_record_ground_state_site

   !> Store a radial ALSDA projection for one response site.
   !>
   !> The site response variable is the total P_site sigma population.  Taking
   !> its transverse fluctuation to preserve the SCF radial spin shape gives
   !>
   !>   K_perp(site) = integral B_xc(r)m(r) dr / (2 M_site^2),
   !>
   !> where M_site is the P_site sigma population supplied separately by
   !> set_site_spin_population.  Here B_xc is the *energy* coefficient
   !> (V_up-V_down)/2.  This is the
   !> direct rotational ALSDA derivative in the code’s m=n_up-n_down and
   !> H=H0+Hvec.sigma convention.  Since the response vertices measure
   !> m^+/-=sigma_x +/- i sigma_y while delta H=(delta B^+m^- +
   !> delta B^-m^+)/2, the circular kernel includes the explicit half.
   subroutine xc_kernel_record_radial_projection(this, isite, projection)
      class(xc_response_kernel_provider), intent(inout) :: this
      integer, intent(in) :: isite
      type(xc_response_radial_projection), intent(in) :: projection
      real(rp) :: mrad

      call require_site(this, isite, 'record_radial_projection')
      mrad = projection%spin_population
      this%site(isite)%radial_spin_population = mrad
      if (abs(projection%charge_population) > tiny(1.0_rp)) then
         this%site(isite)%vxc_scalar = projection%vxc_charge_moment/projection%charge_population
      else
         this%site(isite)%vxc_scalar = 0.0_rp
      end if
      this%site(isite)%bxc_spin_moment = projection%bxc_spin_moment
      this%site(isite)%radial_spin_abs_population = projection%spin_abs_population
      this%site(isite)%radial_vxc_spin_difference = projection%vxc_spin_difference
      this%site(isite)%radial_vxc_spin_difference_abs = projection%vxc_spin_difference_abs
      this%site(isite)%has_radial_projection = .true.
      call finalize_site_k_perp(this%site(isite))
   end subroutine xc_kernel_record_radial_projection

   !> Set the ground-state population in the same site projector used by the
   !> response vertices.  It is intentionally not substituted by a radial
   !> population, since those are distinct projections in the current code.
   subroutine xc_kernel_set_site_spin_population(this, isite, spin_population)
      class(xc_response_kernel_provider), intent(inout) :: this
      integer, intent(in) :: isite
      real(rp), intent(in) :: spin_population

      call require_site(this, isite, 'set_site_spin_population')
      this%site(isite)%spin_population = spin_population
      call finalize_site_k_perp(this%site(isite))
   end subroutine xc_kernel_set_site_spin_population

   subroutine xc_kernel_set_site_signed_spin_population(this, isite, spin_population)
      class(xc_response_kernel_provider), intent(inout) :: this
      integer, intent(in) :: isite
      real(rp), intent(in) :: spin_population
      call require_site(this, isite, 'set_site_signed_spin_population')
      this%site(isite)%signed_spin_population = spin_population
   end subroutine xc_kernel_set_site_signed_spin_population

   !> Record the direction of the ground-state magnetization for the site’s
   !> local ALSDA frame.  Its magnitude is deliberately not used here: the
   !> response projector population is stored separately as spin_population.
   subroutine xc_kernel_set_site_magnetization_direction(this, isite, direction)
      class(xc_response_kernel_provider), intent(inout) :: this
      integer, intent(in) :: isite
      real(rp), intent(in) :: direction(3)
      real(rp) :: norm_direction

      call require_site(this, isite, 'set_site_magnetization_direction')
      norm_direction = sqrt(sum(direction**2))
      if (norm_direction <= tiny(1.0_rp)) then
         error stop 'xc_response_kernel_provider: local magnetization direction is zero'
      end if
      this%site(isite)%magnetization_direction = direction/norm_direction
      this%site(isite)%has_magnetization_direction = .true.
   end subroutine xc_kernel_set_site_magnetization_direction

   subroutine finalize_site_k_perp(site)
      type(xc_response_site), intent(inout) :: site

      ! A radial m=0 record cannot define a rotation profile.  Likewise the
      ! stiffness cannot be normalized before the response projector’s M_site
      ! is known.  This deliberately leaves the provider invalid in either
      ! case rather than replacing M_site with an unrelated quantity.
      if (abs(site%radial_spin_population) <= tiny(1.0_rp) .or. &
         abs(site%spin_population) <= tiny(1.0_rp)) then
         site%bxc_energy = 0.0_rp
         site%k_perp_circular = 0.0_rp
         site%has_k_perp_circular = .false.
         return
      end if
      site%bxc_energy = site%bxc_spin_moment/site%spin_population
      site%k_perp_circular = tddft_circular_source_factor*site%bxc_energy/site%spin_population
      site%has_k_perp_circular = .true.
   end subroutine finalize_site_k_perp

   subroutine xc_kernel_set_site_derivatives(this, isite, dvxc_dn, dvxc_dm, dbxc_dn, dbxc_dm, k_perp_circular, &
      derivatives_validated)
      class(xc_response_kernel_provider), intent(inout) :: this
      integer, intent(in) :: isite
      real(rp), intent(in), optional :: dvxc_dn, dvxc_dm, dbxc_dn, dbxc_dm, k_perp_circular
      logical, intent(in), optional :: derivatives_validated

      call require_site(this, isite, 'set_site_derivatives')
      if (present(dvxc_dn)) then
         this%site(isite)%dvxc_dn = dvxc_dn
         this%site(isite)%has_dvxc_dn = .true.
      end if
      if (present(dvxc_dm)) then
         this%site(isite)%dvxc_dm = dvxc_dm
         this%site(isite)%has_dvxc_dm = .true.
      end if
      if (present(dbxc_dn)) then
         this%site(isite)%dbxc_dn = dbxc_dn
         this%site(isite)%has_dbxc_dn = .true.
      end if
      if (present(dbxc_dm)) then
         this%site(isite)%dbxc_dm = dbxc_dm
         this%site(isite)%has_dbxc_dm = .true.
      end if
      if (present(k_perp_circular)) then
         this%site(isite)%k_perp_circular = k_perp_circular
         this%site(isite)%has_k_perp_circular = .true.
      end if
      if (present(derivatives_validated)) this%full_alsda_derivatives_validated = derivatives_validated
   end subroutine xc_kernel_set_site_derivatives

   !> Return the circular ALSDA scalar for unhalved sigma+/- vertices.
   function circular_transverse_kernel(provider, isite) result(kernel)
      type(xc_response_kernel_provider), intent(in) :: provider
      integer, intent(in) :: isite
      real(rp) :: kernel

      call require_site(provider, isite, 'circular_transverse_kernel')
      if (.not. provider%site(isite)%has_k_perp_circular) then
         error stop 'circular_transverse_kernel: circular transverse kernel is absent'
      end if
      kernel = provider%site(isite)%k_perp_circular
   end function circular_transverse_kernel

   !> Return the Cartesian x/y derivative.  With the frozen unhalved circular
   !> convention it is exactly twice the circular scalar (the historical
   !> spelling is `2.0_rp*circular_transverse_kernel`).
   function cartesian_transverse_kernel(provider, isite) result(kernel)
      type(xc_response_kernel_provider), intent(in) :: provider
      integer, intent(in) :: isite
      real(rp) :: kernel

      kernel = tddft_circular_operator_factor*circular_transverse_kernel(provider, isite)
   end function cartesian_transverse_kernel

   subroutine xc_kernel_full_response_capability(this, supported, reason)
      class(xc_response_kernel_provider), intent(in) :: this
      logical, intent(out) :: supported
      character(len=*), intent(out) :: reason
      integer :: isite

      supported = .false.
      reason = 'full ALSDA derivatives have not been validated for the selected XC route'
      if (.not. allocated(this%site)) then
         reason = 'XC response provider is not initialized'
         return
      end if
      if (.not. this%full_alsda_derivatives_validated) return
      do isite = 1, size(this%site)
         if (.not. this%site(isite)%has_magnetization_direction) then
            write(reason, '(a,i0)') 'site lacks a ground-state magnetization direction: ', isite
            return
         end if
         if (.not. this%site(isite)%has_k_perp_circular .or. .not. this%site(isite)%has_dvxc_dn .or. &
             .not. this%site(isite)%has_dvxc_dm .or. .not. this%site(isite)%has_dbxc_dn .or. &
             .not. this%site(isite)%has_dbxc_dm) then
            write(reason, '(a,i0)') 'site lacks one or more full ALSDA derivatives: ', isite
            return
         end if
      end do
      supported = .true.
      reason = 'validated full ALSDA derivatives supplied by selected XC route'
   end subroutine xc_kernel_full_response_capability

   subroutine require_site(this, isite, caller)
      class(xc_response_kernel_provider), intent(in) :: this
      integer, intent(in) :: isite
      character(len=*), intent(in) :: caller

      if (.not. allocated(this%site)) error stop 'xc_response_kernel_provider: provider is not initialized'
      if (isite < 1 .or. isite > size(this%site)) error stop 'xc_response_kernel_provider: invalid site index'
   end subroutine require_site

end module xc_response_kernel_mod
