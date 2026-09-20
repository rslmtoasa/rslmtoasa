!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: lmto_radial_augmentation
!
!> @brief Small, baseline-only contract for physical LMTO augmentation.
!>
!> The legacy atomic solver owns the numerical construction of the radial
!> solutions.  This module only stores an audit snapshot of its
!> radial outputs and evaluates the linearized augmentation and its density
!> polynomial.  It deliberately does not build a Hamiltonian or solve a radial
!> equation.
!
!> Supported physical contract:
!>   orthogonal reciprocal ham_only, second-order/HOH, collinear, no SOC or
!>   additive operator, and sp/spd angular bases.
!------------------------------------------------------------------------------
module lmto_radial_augmentation_mod

   use precision_mod, only: rp
   implicit none
   private

   real(rp), parameter, public :: scalar_relativistic_c = 274.074_rp

   type, public :: lmto_radial_basis
      !> Number of points and highest represented angular momentum.
      integer :: npoint = 0
      integer :: lmax = -1
      integer :: nspin = 0
      !> Radial mesh and legacy logarithmic-mesh parameters.
      real(rp), allocatable :: rofi(:)
      real(rp) :: mesh_a = 0.0_rp
      real(rp) :: mesh_b = 0.0_rp
      real(rp) :: nuclear_z = 0.0_rp
      !> Immutable accepted-state provenance for the scalar-relativistic
      !> reduction.  The potential is channel dependent; TMC also carries
      !> the linearisation-energy dependence and is therefore l dependent.
      real(rp), allocatable :: potential(:, :), tmc(:, :, :)
      !> Large/small scalar-relativistic radial numerators and their energy
      !> derivatives.  Dimensions are (mesh,l+1,spin), with l=0 in slot 1.
      real(rp), allocatable :: phi_large(:, :, :), phi_small(:, :, :)
      real(rp), allocatable :: phidot_large(:, :, :), phidot_small(:, :, :)
      real(rp), allocatable :: phiddot_large(:, :, :), phiddot_small(:, :, :)
      !> GFAC is the legacy scalar-relativistic large-component metric factor.
      real(rp), allocatable :: gfac(:, :, :)
      !> Energies at which the stored radial derivatives were evaluated, and
      !> the working energies used by the reciprocal Hamiltonian.
      real(rp), allocatable :: enu_radial(:, :), enu_work(:, :)
      logical, allocatable :: channel_present(:, :)
   contains
      procedure :: initialize => lmto_radial_basis_initialize
      procedure :: capture_channel => lmto_radial_basis_capture_channel
      procedure :: supported => lmto_radial_basis_supported
      procedure :: require_supported => lmto_radial_basis_require_supported
      procedure :: reconstruct_state => lmto_reconstruct_state
      procedure :: state_density => lmto_state_density
      procedure :: density_from_moments => lmto_density_from_moments
      procedure :: normalization_integral => lmto_normalization_integral
      procedure :: overlap_integral => lmto_overlap_integral
   end type lmto_radial_basis

   public :: lmto_orbital_l

contains

   !> Initialise storage for one site/species radial basis.
   subroutine lmto_radial_basis_initialize(this, npoint, lmax, nspin)
      class(lmto_radial_basis), intent(inout) :: this
      integer, intent(in) :: npoint, lmax, nspin

      if (npoint < 3 .or. lmax < 0 .or. nspin < 1) then
         error stop 'lmto_radial_basis%initialize: invalid dimensions'
      end if

      if (allocated(this%rofi)) then
         deallocate(this%rofi, this%phi_large, this%phi_small, this%phidot_large, this%phidot_small, &
                    this%phiddot_large, this%phiddot_small, this%gfac, this%enu_radial, this%enu_work, &
                    this%channel_present, this%potential, this%tmc)
      end if

      this%npoint = npoint
      this%lmax = lmax
      this%nspin = nspin
      allocate(this%rofi(npoint))
      allocate(this%phi_large(npoint, lmax + 1, nspin), this%phi_small(npoint, lmax + 1, nspin))
      allocate(this%phidot_large(npoint, lmax + 1, nspin), this%phidot_small(npoint, lmax + 1, nspin))
      allocate(this%phiddot_large(npoint, lmax + 1, nspin), this%phiddot_small(npoint, lmax + 1, nspin))
      allocate(this%gfac(npoint, lmax + 1, nspin))
      allocate(this%enu_radial(lmax + 1, nspin), this%enu_work(lmax + 1, nspin))
      allocate(this%potential(npoint, nspin), this%tmc(npoint, lmax + 1, nspin))
      allocate(this%channel_present(lmax + 1, nspin))

      this%rofi = 0.0_rp
      this%phi_large = 0.0_rp
      this%phi_small = 0.0_rp
      this%phidot_large = 0.0_rp
      this%phidot_small = 0.0_rp
      this%phiddot_large = 0.0_rp
      this%phiddot_small = 0.0_rp
      this%gfac = 1.0_rp
      this%enu_radial = 0.0_rp
      this%enu_work = 0.0_rp
      this%potential = 0.0_rp
      this%tmc = scalar_relativistic_c
      this%channel_present = .false.
   end subroutine lmto_radial_basis_initialize

   !> Copy the outputs of the legacy RSEQSR/PHDFSR pair without changing them.
   !> `g` is the legacy (large,small) packed array of length 2*npoint.
   subroutine lmto_radial_basis_capture_channel(this, l, ispin, energy, rofi, potential, a, b, z, &
                                                 g, gp, gpp, enu_work)
      class(lmto_radial_basis), intent(inout) :: this
      integer, intent(in) :: l, ispin
      real(rp), intent(in) :: energy, rofi(:), potential(:), a, b, z
      real(rp), intent(in) :: g(:), gp(:), gpp(:)
      real(rp), intent(in), optional :: enu_work
      integer :: ir
      real(rp) :: r, tmc, fllp1

      if (l < 0 .or. l > this%lmax .or. ispin < 1 .or. ispin > this%nspin) then
         error stop 'lmto_radial_basis%capture_channel: channel out of range'
      end if
      if (size(rofi) /= this%npoint .or. size(potential) /= this%npoint .or. &
          size(g) < 2*this%npoint .or. size(gp) < 2*this%npoint .or. size(gpp) < 2*this%npoint) then
         error stop 'lmto_radial_basis%capture_channel: inconsistent radial extents'
      end if

      this%rofi = rofi
      this%mesh_a = a
      this%mesh_b = b
      this%nuclear_z = z
      this%potential(:, ispin) = potential
      fllp1 = real(l*(l + 1), rp)
      do ir = 1, this%npoint
         this%phi_large(ir, l + 1, ispin) = g(ir)
         this%phi_small(ir, l + 1, ispin) = g(ir + this%npoint)
         this%phidot_large(ir, l + 1, ispin) = gp(ir)
         this%phidot_small(ir, l + 1, ispin) = gp(ir + this%npoint)
         this%phiddot_large(ir, l + 1, ispin) = gpp(ir)
         this%phiddot_small(ir, l + 1, ispin) = gpp(ir + this%npoint)

         ! This is the same factor used by GINTSR and NEWRHO.  The origin is
         ! not integrated by the legacy Simpson loop, so its diagnostic value
         ! is set to one rather than evaluating the 0/0 expression.
         if (ir == 1 .or. abs(rofi(ir)) <= tiny(1.0_rp)) then
            this%gfac(ir, l + 1, ispin) = 1.0_rp
            this%tmc(ir, l + 1, ispin) = scalar_relativistic_c
         else
            r = rofi(ir)
            tmc = scalar_relativistic_c - (potential(ir) - 2.0_rp*z/r - energy) / scalar_relativistic_c
            this%tmc(ir, l + 1, ispin) = tmc
            this%gfac(ir, l + 1, ispin) = 1.0_rp + fllp1/(tmc*r)**2
         end if
      end do
      this%enu_radial(l + 1, ispin) = energy
      if (present(enu_work)) then
         this%enu_work(l + 1, ispin) = enu_work
      else
         this%enu_work(l + 1, ispin) = energy
      end if
      this%channel_present(l + 1, ispin) = .true.
   end subroutine lmto_radial_basis_capture_channel

   !> Return the production baseline capability decision.
   pure logical function lmto_radial_basis_supported(this, reciprocal_mode, hamiltonian_order, orthogonal, &
                                                     collinear, has_soc, has_extra_operator) result(ok)
      class(lmto_radial_basis), intent(in) :: this
      character(len=*), intent(in) :: reciprocal_mode, hamiltonian_order
      logical, intent(in) :: orthogonal, collinear, has_soc, has_extra_operator

      ok = trim(reciprocal_mode) == 'ham_only' .and. trim(hamiltonian_order) == 'second' .and. &
           orthogonal .and. collinear .and. (.not. has_soc) .and. (.not. has_extra_operator) .and. &
           (this%lmax == 1 .or. this%lmax == 2)
   end function lmto_radial_basis_supported

   !> Fail closed at the API boundary; callers must not silently use the
   !> baseline formula for SOC, overlap, noncollinear, additive, or spdf data.
   subroutine lmto_radial_basis_require_supported(this, reciprocal_mode, hamiltonian_order, orthogonal, &
                                                  collinear, has_soc, has_extra_operator)
      class(lmto_radial_basis), intent(in) :: this
      character(len=*), intent(in) :: reciprocal_mode, hamiltonian_order
      logical, intent(in) :: orthogonal, collinear, has_soc, has_extra_operator

      if (.not. this%supported(reciprocal_mode, hamiltonian_order, orthogonal, collinear, has_soc, &
                               has_extra_operator)) then
         error stop 'LMTO augmentation baseline unsupported for the supplied representation'
      end if
   end subroutine lmto_radial_basis_require_supported

   !> Orbital index in the live spherical ordering: s, p(-1:1), d(-2:2), ...
   pure integer function lmto_orbital_l(iorb) result(l)
      integer, intent(in) :: iorb
      integer :: trial

      if (iorb < 1) then
         l = -1
         return
      end if
      l = -1
      do trial = 0, 8
         if (trial*trial + 1 <= iorb .and. iorb <= (trial + 1)*(trial + 1)) then
            l = trial
            return
         end if
      end do
   end function lmto_orbital_l

   !> Reconstruct radial large/small coefficients for one state.
   !> The angular factor Y_lm is intentionally not multiplied in here.
   subroutine lmto_reconstruct_state(this, energy, coefficients, state_large, state_small)
      class(lmto_radial_basis), intent(in) :: this
      real(rp), intent(in) :: energy
      complex(rp), intent(in) :: coefficients(:, :)
      complex(rp), intent(out) :: state_large(:, :, :), state_small(:, :, :)
      integer :: iorb, ispin, l, ir
      real(rp) :: delta

      if (size(coefficients, 1) /= (this%lmax + 1)**2 .or. size(coefficients, 2) /= this%nspin) then
         error stop 'lmto_reconstruct_state: coefficient shape is not site-local orbital/spin'
      end if
      call lmto_require_complete_channels(this, 'lmto_reconstruct_state')
      if (size(state_large, 1) /= this%npoint .or. size(state_large, 2) /= size(coefficients, 1) .or. &
          size(state_large, 3) /= this%nspin .or. size(state_small, 1) /= this%npoint .or. &
          size(state_small, 2) /= size(coefficients, 1) .or. size(state_small, 3) /= this%nspin) then
         error stop 'lmto_reconstruct_state: output shape mismatch'
      end if

      state_large = (0.0_rp, 0.0_rp)
      state_small = (0.0_rp, 0.0_rp)
      do ispin = 1, this%nspin
         do iorb = 1, size(coefficients, 1)
            l = lmto_orbital_l(iorb)
            if (l < 0 .or. l > this%lmax .or. .not. this%channel_present(l + 1, ispin)) cycle
            delta = energy - this%enu_work(l + 1, ispin)
            do ir = 1, this%npoint
               state_large(ir, iorb, ispin) = coefficients(iorb, ispin) * &
                  (this%phi_large(ir, l + 1, ispin) + delta*this%phidot_large(ir, l + 1, ispin))
               state_small(ir, iorb, ispin) = coefficients(iorb, ispin) * &
                  (this%phi_small(ir, l + 1, ispin) + delta*this%phidot_small(ir, l + 1, ispin))
            end do
         end do
      end do
   end subroutine lmto_reconstruct_state

   !> State-wise radial density Taylor polynomial in the same scalar-relativistic
   !> metric as NEWRHO.  The result is (mesh,l+1,spin), summed over m for each
   !> l.  reconstruct_state is the separate first-order wavefunction oracle;
   !> this routine retains the legacy phi*phiddot term at second order.
   subroutine lmto_state_density(this, energy, coefficients, density)
      class(lmto_radial_basis), intent(in) :: this
      real(rp), intent(in) :: energy
      complex(rp), intent(in) :: coefficients(:, :)
      real(rp), intent(out) :: density(:, :, :)
      integer :: iorb, ispin, l, ir
      real(rp) :: delta, amplitude, cross0, cross1, cross2

      if (size(density, 1) /= this%npoint .or. size(density, 2) /= this%lmax + 1 .or. &
          size(density, 3) /= this%nspin) error stop 'lmto_state_density: output shape mismatch'
      if (size(coefficients, 1) /= (this%lmax + 1)**2 .or. size(coefficients, 2) /= this%nspin) then
         error stop 'lmto_state_density: coefficient shape is not site-local orbital/spin'
      end if
      call lmto_require_complete_channels(this, 'lmto_state_density')
      density = 0.0_rp
      do ispin = 1, this%nspin
         do iorb = 1, size(coefficients, 1)
            l = lmto_orbital_l(iorb)
            if (l < 0 .or. l > this%lmax) cycle
            if (.not. this%channel_present(l + 1, ispin)) cycle
            amplitude = real(coefficients(iorb, ispin)*conjg(coefficients(iorb, ispin)), rp)
            delta = energy - this%enu_work(l + 1, ispin)
            do ir = 1, this%npoint
               cross0 = this%gfac(ir, l + 1, ispin)*this%phi_large(ir, l + 1, ispin)**2 + &
                        this%phi_small(ir, l + 1, ispin)**2
               cross1 = this%gfac(ir, l + 1, ispin)*this%phi_large(ir, l + 1, ispin)* &
                        this%phidot_large(ir, l + 1, ispin) + this%phi_small(ir, l + 1, ispin)* &
                        this%phidot_small(ir, l + 1, ispin)
               cross2 = this%gfac(ir, l + 1, ispin)*(this%phidot_large(ir, l + 1, ispin)**2 + &
                        this%phi_large(ir, l + 1, ispin)*this%phiddot_large(ir, l + 1, ispin)) + &
                        this%phidot_small(ir, l + 1, ispin)**2 + this%phi_small(ir, l + 1, ispin)* &
                        this%phiddot_small(ir, l + 1, ispin)
               density(ir, l + 1, ispin) = density(ir, l + 1, ispin) + amplitude*(cross0 + &
                  2.0_rp*delta*cross1 + delta**2*cross2)
            end do
         end do
      end do
   end subroutine lmto_state_density

   !> Evaluate the legacy Q0/Q1/Q2 density polynomial from coefficient-space
   !> energy moments.  Q1/Q2 are moments of delta=epsilon-E_nu_work.
   subroutine lmto_density_from_moments(this, moments, density)
      class(lmto_radial_basis), intent(in) :: this
      real(rp), intent(in) :: moments(:, :, :) ! (moment, l+1, spin)
      real(rp), intent(out) :: density(:, :, :)
      integer :: ir, l, ispin
      real(rp) :: q0, q1, q2, cross0, cross1, cross2

      if (size(moments, 1) /= 3 .or. size(moments, 2) /= this%lmax + 1 .or. &
          size(moments, 3) /= this%nspin) error stop 'lmto_density_from_moments: moment shape mismatch'
      if (size(density, 1) /= this%npoint .or. size(density, 2) /= this%lmax + 1 .or. &
          size(density, 3) /= this%nspin) error stop 'lmto_density_from_moments: output shape mismatch'
      call lmto_require_complete_channels(this, 'lmto_density_from_moments')

      density = 0.0_rp
      do ispin = 1, this%nspin
         do l = 0, this%lmax
            q0 = moments(1, l + 1, ispin)
            q1 = moments(2, l + 1, ispin) - this%enu_work(l + 1, ispin)*q0
            q2 = moments(3, l + 1, ispin) - 2.0_rp*this%enu_work(l + 1, ispin)*moments(2, l + 1, ispin) + &
                 this%enu_work(l + 1, ispin)**2*q0
            do ir = 1, this%npoint
               cross0 = this%gfac(ir, l + 1, ispin)*this%phi_large(ir, l + 1, ispin)**2 + &
                        this%phi_small(ir, l + 1, ispin)**2
               cross1 = this%gfac(ir, l + 1, ispin)*this%phi_large(ir, l + 1, ispin)*this%phidot_large(ir, l + 1, ispin) + &
                        this%phi_small(ir, l + 1, ispin)*this%phidot_small(ir, l + 1, ispin)
               cross2 = this%gfac(ir, l + 1, ispin)*(this%phidot_large(ir, l + 1, ispin)**2 + &
                        this%phi_large(ir, l + 1, ispin)*this%phiddot_large(ir, l + 1, ispin)) + &
                        this%phidot_small(ir, l + 1, ispin)**2 + &
                        this%phi_small(ir, l + 1, ispin)*this%phiddot_small(ir, l + 1, ispin)
               density(ir, l + 1, ispin) = q0*cross0 + 2.0_rp*q1*cross1 + q2*cross2
            end do
         end do
      end do
   end subroutine lmto_density_from_moments

   subroutine lmto_require_complete_channels(this, caller)
      class(lmto_radial_basis), intent(in) :: this
      character(len=*), intent(in) :: caller

      if (.not. all(this%channel_present)) error stop trim(caller)//': radial channels are incomplete'
   end subroutine lmto_require_complete_channels

   !> Same Simpson metric as GINTSR, for one packed large/small radial pair.
   real(rp) function lmto_normalization_integral(this, l, ispin, which) result(value)
      class(lmto_radial_basis), intent(in) :: this
      integer, intent(in) :: l, ispin
      character(len=*), intent(in) :: which
      integer :: ir
      real(rp), allocatable :: first(:), second(:)

      allocate(first(this%npoint), second(this%npoint))
      select case (trim(which))
      case ('phi'); first = this%phi_large(:, l + 1, ispin); second = this%phi_small(:, l + 1, ispin)
      case ('phidot'); first = this%phidot_large(:, l + 1, ispin); second = this%phidot_small(:, l + 1, ispin)
      case ('phiddot'); first = this%phiddot_large(:, l + 1, ispin); second = this%phiddot_small(:, l + 1, ispin)
      case default; error stop 'lmto_normalization_integral: unknown radial member'
      end select
      value = 0.0_rp
      do ir = 2, this%npoint - 1, 2
         value = value + (this%rofi(ir) + this%mesh_b)*(this%gfac(ir, l + 1, ispin)*first(ir)**2 + second(ir)**2)
      end do
      value = 2.0_rp*value
      do ir = 3, this%npoint - 2, 2
         value = value + (this%rofi(ir) + this%mesh_b)*(this%gfac(ir, l + 1, ispin)*first(ir)**2 + second(ir)**2)
      end do
      value = 2.0_rp*value
      ir = this%npoint
      value = value + (this%rofi(ir) + this%mesh_b)*(this%gfac(ir, l + 1, ispin)*first(ir)**2 + second(ir)**2)
      value = value*this%mesh_a/3.0_rp
      deallocate(first, second)
   end function lmto_normalization_integral

   !> Same Simpson metric as GINTSR for phi against phidot/phiddot.
   real(rp) function lmto_overlap_integral(this, l, ispin, right) result(value)
      class(lmto_radial_basis), intent(in) :: this
      integer, intent(in) :: l, ispin
      character(len=*), intent(in) :: right
      integer :: ir
      value = 0.0_rp
      do ir = 2, this%npoint - 1, 2
         value = value + (this%rofi(ir) + this%mesh_b)*lmto_overlap_point(this, l, ispin, ir, right)
      end do
      value = 2.0_rp*value
      do ir = 3, this%npoint - 2, 2
         value = value + (this%rofi(ir) + this%mesh_b)*lmto_overlap_point(this, l, ispin, ir, right)
      end do
      value = 2.0_rp*value
      ir = this%npoint
      value = value + (this%rofi(ir) + this%mesh_b)*lmto_overlap_point(this, l, ispin, ir, right)
      value = value*this%mesh_a/3.0_rp
   end function lmto_overlap_integral

   pure real(rp) function lmto_overlap_point(this, l, ispin, ir, right) result(value)
      class(lmto_radial_basis), intent(in) :: this
      integer, intent(in) :: l, ispin, ir
      character(len=*), intent(in) :: right
      real(rp) :: rlarge, rsmall
      rlarge = 0.0_rp
      rsmall = 0.0_rp
      select case (trim(right))
      case ('phidot')
         rlarge = this%phidot_large(ir, l + 1, ispin)
         rsmall = this%phidot_small(ir, l + 1, ispin)
      case ('phiddot')
         rlarge = this%phiddot_large(ir, l + 1, ispin)
         rsmall = this%phiddot_small(ir, l + 1, ispin)
      case default
         error stop 'lmto_overlap_integral: unknown right member'
      end select
      value = this%gfac(ir, l + 1, ispin)*this%phi_large(ir, l + 1, ispin)*rlarge + &
              this%phi_small(ir, l + 1, ispin)*rsmall
   end function lmto_overlap_point

end module lmto_radial_augmentation_mod
