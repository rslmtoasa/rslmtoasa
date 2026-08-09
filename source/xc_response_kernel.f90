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
      real(rp) :: vxc_scalar = 0.0_rp
      real(rp) :: bxc_energy = 0.0_rp
      real(rp) :: dvxc_dn = 0.0_rp
      real(rp) :: dvxc_dm = 0.0_rp
      real(rp) :: dbxc_dn = 0.0_rp
      real(rp) :: dbxc_dm = 0.0_rp
      real(rp) :: k_perp = 0.0_rp
      logical :: has_dvxc_dn = .false.
      logical :: has_dvxc_dm = .false.
      logical :: has_dbxc_dn = .false.
      logical :: has_dbxc_dm = .false.
      logical :: has_k_perp = .false.
   end type xc_response_site

   type, public :: xc_response_kernel_provider
      character(len=32) :: functional_label = 'unrecorded'
      type(xc_response_site), allocatable :: site(:)
   contains
      procedure :: initialize => xc_kernel_initialize
      procedure :: clear => xc_kernel_clear
      procedure :: record_ground_state_site => xc_kernel_record_ground_state_site
      procedure :: set_site_derivatives => xc_kernel_set_site_derivatives
   end type xc_response_kernel_provider

   public :: evaluate_ground_state_xc_sample

contains

   !> Evaluate exactly the legacy XCPOT path called by self%VXC0SP.
   !> No constraining field is included here: B_fsm is added by VXC0SP after
   !> this functional call and is not an XC field.
   function evaluate_ground_state_xc_sample(functional, rho_down, rho_up, rho_total, &
      rho_gradient, rho_laplacian, radius) result(sample)
      class(xc), intent(in) :: functional
      real(rp), intent(in) :: rho_down, rho_up, rho_total, radius
      real(rp), intent(in) :: rho_gradient(2), rho_laplacian(2)
      type(xc_response_sample) :: sample
      real(rp) :: v_down, v_up

      call functional%XCPOT(rho_down, rho_up, rho_total, rho_gradient, rho_laplacian, radius, &
         v_down, v_up, sample%exc)
      sample%rho_down = rho_down
      sample%rho_up = rho_up
      sample%vxc_down = v_down
      sample%vxc_up = v_up
      sample%vxc_scalar = 0.5_rp*(v_up + v_down)
      sample%bxc_energy = 0.5_rp*(v_up - v_down)
   end function evaluate_ground_state_xc_sample

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

   subroutine xc_kernel_set_site_derivatives(this, isite, dvxc_dn, dvxc_dm, dbxc_dn, dbxc_dm, k_perp)
      class(xc_response_kernel_provider), intent(inout) :: this
      integer, intent(in) :: isite
      real(rp), intent(in), optional :: dvxc_dn, dvxc_dm, dbxc_dn, dbxc_dm, k_perp

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
      if (present(k_perp)) then
         this%site(isite)%k_perp = k_perp
         this%site(isite)%has_k_perp = .true.
      end if
   end subroutine xc_kernel_set_site_derivatives

   subroutine require_site(this, isite, caller)
      class(xc_response_kernel_provider), intent(in) :: this
      integer, intent(in) :: isite
      character(len=*), intent(in) :: caller

      if (.not. allocated(this%site)) error stop 'xc_response_kernel_provider: provider is not initialized'
      if (isite < 1 .or. isite > size(this%site)) error stop 'xc_response_kernel_provider: invalid site index'
   end subroutine require_site

end module xc_response_kernel_mod
