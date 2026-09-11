!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Numerical-coordinate primitives for the LR-02 response-space audit.
!>
!> This module deliberately stops short of a transition-density evaluator.
!> The live scalar-relativistic radial interface does not yet expose a
!> certified nonspherical lower-component spin-angular convention; the guard
!> below makes that limitation explicit to future response code.
!------------------------------------------------------------------------------
module response_basis_mapping_mod

   use precision_mod, only: rp
   use response_angular_basis_mod, only: response_lm_index, response_angular_pi
   implicit none
   private

   !> Sign of the live reciprocal/Fourier convention exp(+i 2*pi*k.R).
   integer, parameter, public :: response_fourier_phase_sign = 1

   type, public :: response_super_index
      integer :: site = 0
      integer :: response_l = 0
      integer :: response_m = 0
      integer :: radial_point = 0
      integer :: channel = 0
   end type response_super_index

   public :: response_superindex_size
   public :: response_flatten_superindex
   public :: response_unflatten_superindex
   public :: response_simpson_weight
   public :: response_log_mesh_jacobian
   public :: response_volume_measure
   public :: response_weighted_density
   public :: response_physical_density
   public :: response_log_mesh_integral
   public :: response_endpoint_phase
   public :: response_real_space_phase
   public :: response_site_gauge
   public :: response_apply_site_gauge
   public :: response_mapping_supported

contains

   pure integer function response_superindex_size(nsite, response_lmax, npoint, nchannel) result(size_out)
      integer, intent(in) :: nsite, response_lmax, npoint, nchannel
      size_out = nsite*(response_lmax + 1)**2*npoint*nchannel
   end function response_superindex_size

   !> Canonical order is site, response harmonic (L then M=-L..L), radial
   !> point, and explicit density/spin channel.
   subroutine response_flatten_superindex(item, nsite, response_lmax, npoint, nchannel, flat)
      type(response_super_index), intent(in) :: item
      integer, intent(in) :: nsite, response_lmax, npoint, nchannel
      integer, intent(out) :: flat
      integer :: harmonic

      call validate_layout(nsite, response_lmax, npoint, nchannel)
      if (item%site < 1 .or. item%site > nsite .or. item%response_l < 0 .or. &
          item%response_l > response_lmax .or. abs(item%response_m) > item%response_l .or. &
          item%radial_point < 1 .or. item%radial_point > npoint .or. &
          item%channel < 1 .or. item%channel > nchannel) then
         error stop 'response_flatten_superindex: index is outside the response layout'
      end if
      harmonic = response_lm_index(item%response_l, item%response_m) - 1
      flat = (((item%site - 1)*(response_lmax + 1)**2 + harmonic)*npoint + &
         item%radial_point - 1)*nchannel + item%channel
   end subroutine response_flatten_superindex

   subroutine response_unflatten_superindex(flat, nsite, response_lmax, npoint, nchannel, item)
      integer, intent(in) :: flat, nsite, response_lmax, npoint, nchannel
      type(response_super_index), intent(out) :: item
      integer :: zero_based, harmonic

      call validate_layout(nsite, response_lmax, npoint, nchannel)
      if (flat < 1 .or. flat > response_superindex_size(nsite, response_lmax, npoint, nchannel)) then
         error stop 'response_unflatten_superindex: flat index is outside the response layout'
      end if

      zero_based = flat - 1
      item%channel = mod(zero_based, nchannel) + 1
      zero_based = zero_based/nchannel
      item%radial_point = mod(zero_based, npoint) + 1
      zero_based = zero_based/npoint
      harmonic = mod(zero_based, (response_lmax + 1)**2) + 1
      item%site = zero_based/(response_lmax + 1)**2 + 1
      call response_lm_from_harmonic(harmonic, item%response_l, item%response_m)
   end subroutine response_unflatten_superindex

   pure real(rp) function response_simpson_weight(ir, npoint) result(weight)
      integer, intent(in) :: ir, npoint
      weight = 2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp
      if (ir == 1 .or. ir == npoint) weight = 1.0_rp/3.0_rp
   end function response_simpson_weight

   pure real(rp) function response_log_mesh_jacobian(a, b, radius) result(value)
      real(rp), intent(in) :: a, b, radius
      value = a*(radius + b)
   end function response_log_mesh_jacobian

   !> Radial part of r^2 dr on the production logarithmic mesh.  The angular
   !> integral is not included here.
   pure real(rp) function response_volume_measure(a, b, radius) result(value)
      real(rp), intent(in) :: a, b, radius
      value = radius**2*response_log_mesh_jacobian(a, b, radius)
   end function response_volume_measure

   pure real(rp) function response_weighted_density(radius, physical_density) result(value)
      real(rp), intent(in) :: radius, physical_density
      value = 4.0_rp*response_angular_pi*radius**2*physical_density
   end function response_weighted_density

   pure real(rp) function response_physical_density(radius, weighted_density) result(value)
      real(rp), intent(in) :: radius, weighted_density
      if (abs(radius) <= tiny(1.0_rp)) then
         value = 0.0_rp
      else
         value = weighted_density/(4.0_rp*response_angular_pi*radius**2)
      end if
   end function response_physical_density

   function response_log_mesh_integral(values, radius, a, b) result(integral)
      real(rp), intent(in) :: values(:), radius(:), a, b
      real(rp) :: integral
      integer :: ir

      if (size(values) /= size(radius) .or. size(radius) < 3 .or. mod(size(radius), 2) == 0) then
         error stop 'response_log_mesh_integral: Simpson mesh must have odd size'
      end if
      integral = 0.0_rp
      do ir = 1, size(radius)
         integral = integral + response_simpson_weight(ir, size(radius))* &
            response_log_mesh_jacobian(a, b, radius(ir))*values(ir)
      end do
   end function response_log_mesh_integral

   !> The production Fourier phase is exp(+i 2*pi*k.R), while the current
   !> cluster vector is an endpoint displacement R+tau_b-tau_a.  The site
   !> gauge under k -> k+G is therefore exp(-i 2*pi*G.tau_a).
   pure complex(rp) function response_endpoint_phase(site_tau, reciprocal_vector) result(phase)
      real(rp), intent(in) :: site_tau(3), reciprocal_vector(3)
      real(rp) :: angle
      angle = -2.0_rp*response_angular_pi*dot_product(reciprocal_vector, site_tau)
      phase = cmplx(cos(angle), sin(angle), rp)
   end function response_endpoint_phase

   !> Fourier phase for a real-space response pair.  `translation` is the
   !> direct-lattice pair translation and the endpoint displacement is
   !> R+tau_right-tau_left, matching reciprocal_fourier and LR-03.
   pure complex(rp) function response_real_space_phase(reciprocal_vector, translation, tau_left, tau_right) result(phase)
      real(rp), intent(in) :: reciprocal_vector(3), translation(3), tau_left(3), tau_right(3)
      real(rp) :: angle

      angle = real(response_fourier_phase_sign, rp)*2.0_rp*response_angular_pi*dot_product(reciprocal_vector, &
         translation + tau_right - tau_left)
      phase = cmplx(cos(angle), sin(angle), rp)
   end function response_real_space_phase

   subroutine response_site_gauge(site_tau, reciprocal_vector, gauge)
      real(rp), intent(in) :: site_tau(:, :), reciprocal_vector(3)
      complex(rp), intent(out) :: gauge(:)
      integer :: isite

      if (size(site_tau, 1) /= 3 .or. size(site_tau, 2) /= size(gauge)) then
         error stop 'response_site_gauge: site-coordinate shape mismatch'
      end if
      do isite = 1, size(gauge)
         gauge(isite) = response_endpoint_phase(site_tau(:, isite), reciprocal_vector)
      end do
   end subroutine response_site_gauge

   subroutine response_apply_site_gauge(coefficients, site_tau, reciprocal_vector, transformed)
      complex(rp), intent(in) :: coefficients(:, :)
      real(rp), intent(in) :: site_tau(:, :), reciprocal_vector(3)
      complex(rp), intent(out) :: transformed(:, :)
      complex(rp) :: gauge(size(coefficients, 2))
      integer :: isite

      if (any(shape(transformed) /= shape(coefficients)) .or. size(site_tau, 1) /= 3 .or. &
          size(site_tau, 2) /= size(coefficients, 2)) then
         error stop 'response_apply_site_gauge: coefficient/site-coordinate shape mismatch'
      end if
      call response_site_gauge(site_tau, reciprocal_vector, gauge)
      do isite = 1, size(gauge)
         transformed(:, isite) = gauge(isite)*coefficients(:, isite)
      end do
   end subroutine response_apply_site_gauge

   !> `full_scalar_relativistic_spin_angular` is intentionally an explicit
   !> input.  The current production radial contract supplies .false.; a
   !> future augmentation seam may set it only after proving the missing
   !> lower-component angular convention.
   pure logical function response_mapping_supported(orbital_lmax, scalar_relativistic, collinear, no_soc, &
                                                    no_extra_operator, full_scalar_relativistic_spin_angular) result(ok)
      integer, intent(in) :: orbital_lmax
      logical, intent(in) :: scalar_relativistic, collinear, no_soc, no_extra_operator
      logical, intent(in) :: full_scalar_relativistic_spin_angular

      ok = (orbital_lmax == 1 .or. orbital_lmax == 2) .and. scalar_relativistic .and. collinear .and. &
           no_soc .and. no_extra_operator .and. full_scalar_relativistic_spin_angular
   end function response_mapping_supported

   subroutine validate_layout(nsite, response_lmax, npoint, nchannel)
      integer, intent(in) :: nsite, response_lmax, npoint, nchannel
      if (nsite < 1 .or. response_lmax < 0 .or. npoint < 1 .or. nchannel < 1) then
         error stop 'response basis layout: dimensions must be positive'
      end if
   end subroutine validate_layout

   pure subroutine response_lm_from_harmonic(harmonic, l, m)
      integer, intent(in) :: harmonic
      integer, intent(out) :: l, m
      integer :: trial

      l = -1
      m = 0
      do trial = 0, 64
         if (trial*trial + 1 <= harmonic .and. harmonic <= (trial + 1)*(trial + 1)) then
            l = trial
            m = harmonic - trial*trial - trial - 1
            return
         end if
      end do
   end subroutine response_lm_from_harmonic

end module response_basis_mapping_mod
