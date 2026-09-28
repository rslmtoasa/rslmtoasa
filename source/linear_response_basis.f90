!------------------------------------------------------------------------------
! Linear-response basis, angular, mapping, endpoint, Pauli, product, GF, and
! projected-site services moved under linear_response_mod.
!------------------------------------------------------------------------------
submodule (linear_response_mod) linear_response_basis
   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use, intrinsic :: ieee_arithmetic
   use radial_ground_state_mod, only: radial_ground_state, radial_simpson_weight
   implicit none

   ! --- private state from response_basis_mapping_mod ---

   !> Sign of the live reciprocal/Fourier convention exp(+i 2*pi*k.R).

   ! --- private state from lr_pauli_transition_vertex_mod ---

   !> Capability tuple for the initial LR-05 production vertex.

   !> One reciprocal endpoint eigenstate in production site-blocked ordering.
   !>
   !> The coefficient vector is packed as
   !> `(site, spin-up orbitals, spin-down orbitals)`, with the orbital order
   !> `(s),(p,-1:1),(d,-2:2),...` supplied by `basis_mod`/LMTO.

   ! --- private state from lr_gf_endpoint_augmentation_mod ---


   !> Capability tuple accepted by this endpoint adapter.

   !> Callback for a production effective-Hamiltonian action.
   !>
   !> `seed` is a localized source-site block.  The implementation applies the
   !> same H_eff action used by the native solver and returns the destination
   !> site block.  For destination=source the caller subtracts E_nu_work to
   !> obtain h^gamma_aa; for an offsite block the returned H block is already
   !> h^gamma_ab.


   !> Simple dense provider used by small finite fixtures and by callers that
   !> already hold the completed H_eff matrix.  Production callers can provide
   !> a wrapper around their native two-sweep/action routine instead.

   !> Factored endpoint-augmented Green-function block.
   !>
   !> The four coefficient-space members are the explicit radial branches in
   !>   Phi G Phi^dagger,
   !>   Phidot (hG) Phi^dagger,
   !>   Phi (Gh) Phidot^dagger,
   !>   Phidot (hGh) Phidot^dagger.
   !>
   !> `phi_left/right` and `phidot_left/right` contain large-component radial
   !> numerators repeated in orbital/spin ordering.  Angular harmonics are not
   !> duplicated in storage; `evaluate_point` inserts the certified
   !> response_angular_basis convention and divides U_l by r at the endpoint.

   ! --- private state from lr_projected_site_spin_mod ---

   real(rp), parameter :: kB_Ry_per_K = 6.3336814e-6_rp
   real(rp), parameter :: occupation_kT_floor = 1.0e-10_rp
   real(rp), parameter :: mesh_tolerance = 2.0e-13_rp




   !> Immutable-after-initialization DRESP-01 contract metadata and operations.

   interface
      subroutine zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
         import :: rp
         character(len=1), intent(in) :: jobu, jobvt
         integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
         complex(rp), intent(inout) :: a(lda, *)
         real(rp), intent(out) :: s(*)
         complex(rp), intent(out) :: u(ldu, *), vt(ldvt, *), work(*)
         real(rp), intent(inout) :: rwork(*)
         integer, intent(out) :: info
      end subroutine zgesvd
   end interface

contains

   ! --- from response_angular_basis_mod ---

   !> The complete product of an orbital basis truncated at lmax has this
   !> response angular cutoff unless symmetry removes channels.
   module pure integer function response_lmax(orbital_lmax) result(value)
      integer, intent(in) :: orbital_lmax
      value = 2*orbital_lmax
   end function response_lmax

   module pure logical function response_product_space_complete(orbital_lmax, response_lmax_in) result(ok)
      integer, intent(in) :: orbital_lmax, response_lmax_in
      ok = orbital_lmax >= 0 .and. response_lmax_in >= response_lmax(orbital_lmax)
   end function response_product_space_complete

   !> One-based index in the live (l,m) order used by hcpx for complex sp/spd.
   module pure integer function response_lm_index(l, m) result(index)
      integer, intent(in) :: l, m
      if (l < 0 .or. abs(m) > l) then
         index = -1
      else
         index = l*l + l + m + 1
      end if
   end function response_lm_index

   module pure subroutine response_lm_from_index(index, l, m)
      integer, intent(in) :: index
      integer, intent(out) :: l, m
      integer :: trial

      l = -1
      m = 0
      if (index < 1) return
      do trial = 0, 64
         if (trial*trial + 1 <= index .and. index <= (trial + 1)*(trial + 1)) then
            l = trial
            m = index - trial*trial - trial - 1
            return
         end if
      end do
   end subroutine response_lm_from_index

   !> Complex spherical harmonic in the convention implied by hcpx.
   module pure function response_harmonic(l, m, theta, phi) result(value)
      integer, intent(in) :: l, m
      real(rp), intent(in) :: theta, phi
      complex(rp) :: value
      complex(rp) :: standard
      real(rp) :: x, p, norm
      integer :: ma

      value = cmplx(0.0_rp, 0.0_rp, rp)
      if (l < 0 .or. abs(m) > l) return

      x = cos(theta)
      ma = abs(m)
      p = associated_legendre_without_condon(l, ma, x)
      norm = sqrt(real(2*l + 1, rp)/(4.0_rp*response_angular_pi) * &
         factorial_real(l - ma)/factorial_real(l + ma))
      if (m >= 0) then
         standard = cmplx(norm*minus_one_power(ma)*p*cos(real(m, rp)*phi), &
            norm*minus_one_power(ma)*p*sin(real(m, rp)*phi), rp)
      else
         ! Condon--Shortley Y(l,-m)=(-1)^m conjugate(Y(l,m)).
         standard = minus_one_power(ma)*conjg(cmplx(norm*minus_one_power(ma)*p*cos(real(ma, rp)*phi), &
            norm*minus_one_power(ma)*p*sin(real(ma, rp)*phi), rp))
      end if

      ! hcpx columns are the conjugates of the usual Condon--Shortley
      ! harmonics; see the explicit p and d columns in source/math.f90.
      value = conjg(standard)
   end function response_harmonic

   !> Gaunt coefficient for the live hcpx harmonic convention:
   !>
   !>   integral Y(code)_LM^* Y(code)_lm^* Y(code)_l'm' dOmega.
   !>
   !> The expression is the Condon--Shortley 3j formula transformed by
   !> Y(code)=conjugate(Y(CS)).
   module pure real(rp) function response_gaunt(l1, m1, l2, m2, lout, mout) result(value)
      integer, intent(in) :: l1, m1, l2, m2, lout, mout
      real(rp) :: prefactor

      value = 0.0_rp
      if (l1 < 0 .or. l2 < 0 .or. lout < 0) return
      if (abs(m1) > l1 .or. abs(m2) > l2 .or. abs(mout) > lout) return
      if (mout /= m2 - m1) return
      if (lout < abs(l1 - l2) .or. lout > l1 + l2) return
      if (mod(lout + l1 + l2, 2) /= 0) return

      prefactor = sqrt(real((2*lout + 1)*(2*l1 + 1)*(2*l2 + 1), rp) / &
         (4.0_rp*response_angular_pi))
      value = minus_one_power(m2)*prefactor* &
         wigner_3j(lout, l1, l2, 0, 0, 0)*wigner_3j(lout, l1, l2, mout, m1, -m2)
   end function response_gaunt

   !> Fill coefficients for Y(code)_lm^* Y(code)_l'm' in response order.
   module subroutine response_product_coefficients(l, m, lp, mp, coefficients)
      integer, intent(in) :: l, m, lp, mp
      real(rp), intent(out) :: coefficients(:)
      integer :: lout, mout, expected_size

      expected_size = (l + lp + 1)**2
      if (l < 0 .or. lp < 0 .or. abs(m) > l .or. abs(mp) > lp .or. &
          size(coefficients) /= expected_size) then
         error stop 'response_product_coefficients: invalid product dimensions'
      end if
      coefficients = 0.0_rp
      do lout = 0, l + lp
         do mout = -lout, lout
            coefficients(response_lm_index(lout, mout)) = response_gaunt(l, m, lp, mp, lout, mout)
         end do
      end do
   end subroutine response_product_coefficients

   pure real(rp) function factorial_real(n) result(value)
      integer, intent(in) :: n
      integer :: i

      if (n < 0) then
         value = 0.0_rp
         return
      end if
      value = 1.0_rp
      do i = 2, n
         value = value*real(i, rp)
      end do
   end function factorial_real

   pure real(rp) function minus_one_power(n) result(value)
      integer, intent(in) :: n
      if (mod(abs(n), 2) == 0) then
         value = 1.0_rp
      else
         value = -1.0_rp
      end if
   end function minus_one_power

   pure real(rp) function associated_legendre_without_condon(l, m, x) result(value)
      integer, intent(in) :: l, m
      real(rp), intent(in) :: x
      real(rp) :: pmm, pmmp1, pll, somx2
      integer :: ell

      value = 0.0_rp
      if (m < 0 .or. l < m) return
      pmm = 1.0_rp
      if (m > 0) then
         somx2 = sqrt(max(0.0_rp, (1.0_rp - x)*(1.0_rp + x)))
         do ell = 1, m
            pmm = pmm*real(2*ell - 1, rp)*somx2
         end do
      end if
      if (l == m) then
         value = pmm
         return
      end if
      pmmp1 = x*real(2*m + 1, rp)*pmm
      if (l == m + 1) then
         value = pmmp1
         return
      end if
      do ell = m + 2, l
         pll = (real(2*ell - 1, rp)*x*pmmp1 - real(ell + m - 1, rp)*pmm)/real(ell - m, rp)
         pmm = pmmp1
         pmmp1 = pll
      end do
      value = pmmp1
   end function associated_legendre_without_condon

   pure real(rp) function wigner_3j(j1, j2, j3, m1, m2, m3) result(value)
      integer, intent(in) :: j1, j2, j3, m1, m2, m3
      integer :: z, zmin, zmax
      real(rp) :: delta, norm, denominator, sum_term

      value = 0.0_rp
      if (j1 < 0 .or. j2 < 0 .or. j3 < 0) return
      if (abs(m1) > j1 .or. abs(m2) > j2 .or. abs(m3) > j3) return
      if (m1 + m2 + m3 /= 0) return
      if (j3 < abs(j1 - j2) .or. j3 > j1 + j2) return

      delta = sqrt(factorial_real(j1 + j2 - j3)*factorial_real(j1 - j2 + j3)* &
         factorial_real(-j1 + j2 + j3)/factorial_real(j1 + j2 + j3 + 1))
      norm = sqrt(factorial_real(j1 + m1)*factorial_real(j1 - m1)* &
         factorial_real(j2 + m2)*factorial_real(j2 - m2)* &
         factorial_real(j3 + m3)*factorial_real(j3 - m3))
      zmin = max(0, max(-j3 + j2 - m1, -j3 + j1 + m2))
      zmax = min(j1 + j2 - j3, min(j1 - m1, j2 + m2))
      if (zmin > zmax) return

      sum_term = 0.0_rp
      do z = zmin, zmax
         denominator = factorial_real(z)*factorial_real(j1 + j2 - j3 - z)* &
            factorial_real(j1 - m1 - z)*factorial_real(j2 + m2 - z)* &
            factorial_real(j3 - j2 + m1 + z)*factorial_real(j3 - j1 - m2 + z)
         sum_term = sum_term + minus_one_power(z)/denominator
      end do
      value = minus_one_power(j1 - j2 - m3)*delta*norm*sum_term
   end function wigner_3j


   ! --- from response_basis_mapping_mod ---

   module pure integer function response_superindex_size(nsite, response_lmax, npoint, nchannel) result(size_out)
      integer, intent(in) :: nsite, response_lmax, npoint, nchannel
      size_out = nsite*(response_lmax + 1)**2*npoint*nchannel
   end function response_superindex_size

   !> Canonical order is site, response harmonic (L then M=-L..L), radial
   !> point, and explicit density/spin channel.
   module subroutine response_flatten_superindex(item, nsite, response_lmax, npoint, nchannel, flat)
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

   module subroutine response_unflatten_superindex(flat, nsite, response_lmax, npoint, nchannel, item)
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

   module pure real(rp) function response_simpson_weight(ir, npoint) result(weight)
      integer, intent(in) :: ir, npoint
      weight = 2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp
      if (ir == 1 .or. ir == npoint) weight = 1.0_rp/3.0_rp
   end function response_simpson_weight

   module pure real(rp) function response_log_mesh_jacobian(a, b, radius) result(value)
      real(rp), intent(in) :: a, b, radius
      value = a*(radius + b)
   end function response_log_mesh_jacobian

   !> Radial part of r^2 dr on the production logarithmic mesh.  The angular
   !> integral is not included here.
   module pure real(rp) function response_volume_measure(a, b, radius) result(value)
      real(rp), intent(in) :: a, b, radius
      value = radius**2*response_log_mesh_jacobian(a, b, radius)
   end function response_volume_measure

   module pure real(rp) function response_weighted_density(radius, physical_density) result(value)
      real(rp), intent(in) :: radius, physical_density
      value = 4.0_rp*response_angular_pi*radius**2*physical_density
   end function response_weighted_density

   module pure real(rp) function response_physical_density(radius, weighted_density) result(value)
      real(rp), intent(in) :: radius, weighted_density
      if (abs(radius) <= tiny(1.0_rp)) then
         value = 0.0_rp
      else
         value = weighted_density/(4.0_rp*response_angular_pi*radius**2)
      end if
   end function response_physical_density

   module function response_log_mesh_integral(values, radius, a, b) result(integral)
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
   module pure complex(rp) function response_endpoint_phase(site_tau, reciprocal_vector) result(phase)
      real(rp), intent(in) :: site_tau(3), reciprocal_vector(3)
      real(rp) :: angle
      angle = -2.0_rp*response_angular_pi*dot_product(reciprocal_vector, site_tau)
      phase = cmplx(cos(angle), sin(angle), rp)
   end function response_endpoint_phase

   !> Fourier phase for a real-space response pair.  `translation` is the
   !> direct-lattice pair translation and the endpoint displacement is
   !> R+tau_right-tau_left, matching reciprocal_fourier and LR-03.
   module pure complex(rp) function response_real_space_phase(reciprocal_vector, translation, tau_left, tau_right) result(phase)
      real(rp), intent(in) :: reciprocal_vector(3), translation(3), tau_left(3), tau_right(3)
      real(rp) :: angle

      angle = real(response_fourier_phase_sign, rp)*2.0_rp*response_angular_pi*dot_product(reciprocal_vector, &
         translation + tau_right - tau_left)
      phase = cmplx(cos(angle), sin(angle), rp)
   end function response_real_space_phase

   module subroutine response_site_gauge(site_tau, reciprocal_vector, gauge)
      real(rp), intent(in) :: site_tau(:, :), reciprocal_vector(3)
      complex(rp), intent(out) :: gauge(:)
      integer :: isite

      do isite = 1, size(gauge)
         gauge(isite) = response_endpoint_phase(site_tau(:, isite), reciprocal_vector)
      end do
   end subroutine response_site_gauge

   module subroutine response_apply_site_gauge(coefficients, site_tau, reciprocal_vector, transformed)
      complex(rp), intent(in) :: coefficients(:, :)
      real(rp), intent(in) :: site_tau(:, :), reciprocal_vector(3)
      complex(rp), intent(out) :: transformed(:, :)
      complex(rp) :: gauge(size(coefficients, 2))
      integer :: isite

      call response_site_gauge(site_tau, reciprocal_vector, gauge)
      do isite = 1, size(gauge)
         transformed(:, isite) = gauge(isite)*coefficients(:, isite)
      end do
   end subroutine response_apply_site_gauge

   !> `full_scalar_relativistic_spin_angular` is intentionally an explicit
   !> input.  The current production radial contract supplies .false.; a
   !> future augmentation seam may set it only after proving the missing
   !> lower-component angular convention.
   module pure logical function response_mapping_supported(orbital_lmax, scalar_relativistic, collinear, no_soc, &
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


   ! --- from lr_response_space_mod ---

   !> Build the physical radial metric for a pointwise density on the LR-01
   !> logarithmic mesh.  The angular integral is already represented by
   !> normalized response harmonics and is therefore not included here.
   module subroutine response_build_radial_metric(a, b, radius, weights)
      real(rp), intent(in) :: a, b, radius(:)
      real(rp), intent(out) :: weights(:)
      integer :: ir, npoint

      npoint = size(radius)
      if (size(weights) /= npoint .or. npoint < 3 .or. mod(npoint, 2) == 0) then
         error stop 'response_build_radial_metric: Simpson mesh must have odd size'
      end if
      if (a <= 0.0_rp .or. b <= 0.0_rp) then
         error stop 'response_build_radial_metric: logarithmic mesh parameters must be positive'
      end if
      if (radius(1) /= 0.0_rp) then
         error stop 'response_build_radial_metric: production response mesh must include r=0'
      end if
      do ir = 2, npoint
         if (radius(ir) <= radius(ir - 1)) then
            error stop 'response_build_radial_metric: radius must be strictly increasing'
         end if
      end do

      do ir = 1, npoint
         weights(ir) = response_simpson_weight(ir, npoint)*response_volume_measure(a, b, radius(ir))
      end do
   end subroutine response_build_radial_metric

   module subroutine response_space_initialize(this, nsite, response_lmax, radius, a, b, nchannel)
      class(response_space_layout), intent(out) :: this
      integer, intent(in) :: nsite, response_lmax, nchannel
      real(rp), intent(in) :: radius(:), a, b
      integer :: flat, radial_point

      if (nsite < 1 .or. response_lmax < 0 .or. nchannel < 1) then
         error stop 'response_space_layout%initialize: invalid response dimensions'
      end if
      if (size(radius) < 3 .or. mod(size(radius), 2) == 0) then
         error stop 'response_space_layout%initialize: invalid Simpson mesh size'
      end if

      this%nsite = nsite
      this%response_lmax = response_lmax
      this%npoint = size(radius)
      this%nchannel = nchannel
      this%ndim = response_superindex_size(nsite, response_lmax, size(radius), nchannel)
      this%active_dimension = response_superindex_size(nsite, response_lmax, size(radius) - 1, nchannel)
      this%a = a
      this%b = b
      allocate(this%radius(size(radius)), this%radial_weights(size(radius)), this%metric_weights(this%ndim))
      this%radius = radius
      call response_build_radial_metric(a, b, radius, this%radial_weights)

      ! The LR-02R order is site, harmonic, radial point, channel.  Repeating
      ! the radial metric in this order makes every contraction use the same
      ! super-index without introducing a second indexing implementation.
      do flat = 1, this%ndim
         radial_point = mod((flat - 1)/nchannel, this%npoint) + 1
         this%metric_weights(flat) = this%radial_weights(radial_point)
      end do
   end subroutine response_space_initialize

   module subroutine response_metric_weights(space, weights)
      type(response_space_layout), intent(in) :: space
      real(rp), intent(out) :: weights(:)

      call require_layout(space)
      weights = space%metric_weights
   end subroutine response_metric_weights

   module function response_vector_inner_product(space, left, right) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: left(:), right(:)
      complex(rp) :: value

      call require_layout(space)
      call require_vector(space, left, 'response_vector_inner_product: left vector shape mismatch')
      call require_vector(space, right, 'response_vector_inner_product: right vector shape mismatch')
      value = sum(conjg(left)*space%metric_weights*right)
   end function response_vector_inner_product

   module function response_vector_norm(space, vector) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: vector(:)
      real(rp) :: value
      complex(rp) :: inner

      inner = response_vector_inner_product(space, vector, vector)
      if (real(inner, rp) < -100.0_rp*epsilon(1.0_rp)*max(1.0_rp, abs(inner))) then
         error stop 'response_vector_norm: negative metric norm'
      end if
      value = sqrt(max(0.0_rp, real(inner, rp)))
   end function response_vector_norm

   module subroutine response_apply_operator(space, operator, field, result)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: operator(:, :), field(:)
      complex(rp), intent(out) :: result(:)

      call require_operator(space, operator, 'response_apply_operator: operator shape mismatch')
      call require_vector(space, field, 'response_apply_operator: field shape mismatch')
      call require_vector(space, result, 'response_apply_operator: result shape mismatch')
      result = matmul(operator, field)
   end subroutine response_apply_operator

   !> Composition follows directly from B=A W.  If C=A W and D=D_raw W,
   !> then the composite canonical matrix is C D; no additional local metric
   !> factor is allowed here.
   module subroutine response_compose_operators(space, first, second, composed)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: first(:, :), second(:, :)
      complex(rp), intent(out) :: composed(:, :)

      call require_operator(space, first, 'response_compose_operators: first shape mismatch')
      call require_operator(space, second, 'response_compose_operators: second shape mismatch')
      call require_operator(space, composed, 'response_compose_operators: result shape mismatch')
      composed = matmul(first, second)
   end subroutine response_compose_operators

   module subroutine response_identity_operator(space, identity)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(out) :: identity(:, :)
      integer :: flat

      call require_operator(space, identity, 'response_identity_operator: result shape mismatch')
      identity = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         identity(flat, flat) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end subroutine response_identity_operator

   !> Spherical local radial scalars are diagonal in site, L, M, radial point,
   !> and channel.  A rank-2 input is shared by all channels.
   module subroutine response_local_operator_site_radial(space, values, operator)
      type(response_space_layout), intent(in) :: space
      real(rp), intent(in) :: values(:, :)
      complex(rp), intent(out) :: operator(:, :)
      type(response_super_index) :: item
      integer :: flat

      call require_operator(space, operator, 'response_local_operator: result shape mismatch')
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, &
            space%npoint, space%nchannel, item)
         operator(flat, flat) = cmplx(values(item%site, item%radial_point), 0.0_rp, rp)
      end do
   end subroutine response_local_operator_site_radial

   module subroutine response_local_operator_site_radial_channel(space, values, operator)
      type(response_space_layout), intent(in) :: space
      real(rp), intent(in) :: values(:, :, :)
      complex(rp), intent(out) :: operator(:, :)
      type(response_super_index) :: item
      integer :: flat

      call require_operator(space, operator, 'response_local_operator: result shape mismatch')
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, &
            space%npoint, space%nchannel, item)
         operator(flat, flat) = cmplx(values(item%site, item%radial_point, item%channel), 0.0_rp, rp)
      end do
   end subroutine response_local_operator_site_radial_channel

   !> Complex-valued spherical local scalar.  This has exactly the same
   !> pointwise/canonical action as the real overload and is needed by static
   !> response routes evaluated with a finite retarded broadening.
   module subroutine response_local_operator_complex_site_radial(space, values, operator)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: values(:, :)
      complex(rp), intent(out) :: operator(:, :)
      type(response_super_index) :: item
      integer :: flat

      call require_operator(space, operator, 'response_local_operator: result shape mismatch')
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, &
            space%npoint, space%nchannel, item)
         operator(flat, flat) = values(item%site, item%radial_point)
      end do
   end subroutine response_local_operator_complex_site_radial

   module subroutine response_local_operator_complex_site_radial_channel(space, values, operator)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: values(:, :, :)
      complex(rp), intent(out) :: operator(:, :)
      type(response_super_index) :: item
      integer :: flat

      call require_operator(space, operator, 'response_local_operator: result shape mismatch')
      operator = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         call response_unflatten_superindex(flat, space%nsite, space%response_lmax, &
            space%npoint, space%nchannel, item)
         operator(flat, flat) = values(item%site, item%radial_point, item%channel)
      end do
   end subroutine response_local_operator_complex_site_radial_channel

   !> Metric adjoint is defined by <x,B y>_W=<B^dagger x,y>_W.
   !> The origin entries form a null subspace of the quadrature metric.  The
   !> canonical extension is the ordinary conjugate transpose within that
   !> null subspace, with no active/null mixing.  A finite quadrature operator
   !> must not map a null input at the origin into an active output.
   module subroutine response_operator_adjoint(space, operator, adjoint)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: operator(:, :)
      complex(rp), intent(out) :: adjoint(:, :)
      integer :: i, j

      call require_operator(space, operator, 'response_operator_adjoint: operator shape mismatch')
      call require_operator(space, adjoint, 'response_operator_adjoint: result shape mismatch')
      do j = 1, space%ndim
         if (space%metric_weights(j) == 0.0_rp) then
            do i = 1, space%ndim
               if (space%metric_weights(i) > 0.0_rp .and. operator(i, j) /= cmplx(0.0_rp, 0.0_rp, rp)) then
                  error stop 'response_operator_adjoint: null-origin input is not quadrature compatible'
               end if
            end do
         end if
      end do

      adjoint = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, space%ndim
         do j = 1, space%ndim
            if (space%metric_weights(i) > 0.0_rp .and. space%metric_weights(j) > 0.0_rp) then
               adjoint(i, j) = conjg(operator(j, i))*space%metric_weights(j)/space%metric_weights(i)
            elseif (space%metric_weights(i) == 0.0_rp .and. space%metric_weights(j) == 0.0_rp) then
               adjoint(i, j) = conjg(operator(j, i))
            end if
         end do
      end do
   end subroutine response_operator_adjoint

   !> Trace is the trace on the positive-measure response space.  Null origin
   !> point values are intentionally excluded.
   module function response_operator_trace(space, operator) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: operator(:, :)
      complex(rp) :: value
      integer :: flat

      call require_operator(space, operator, 'response_operator_trace: operator shape mismatch')
      value = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, space%ndim
         if (space%metric_weights(flat) > 0.0_rp) value = value + operator(flat, flat)
      end do
   end function response_operator_trace

   module function response_rigid_vector_norm(space, rigid_vector) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: rigid_vector(:)
      real(rp) :: value

      value = response_vector_norm(space, rigid_vector)
   end function response_rigid_vector_norm

   module function response_rigid_vector_overlap(space, vector, rigid_vector) result(value)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: vector(:), rigid_vector(:)
      complex(rp) :: value

      value = response_vector_inner_product(space, vector, rigid_vector)
   end function response_rigid_vector_overlap

   !> Convert a raw pointwise kernel A to canonical B=A W.  The origin column
   !> is exactly zero because its quadrature measure is exactly zero.
   module subroutine response_raw_to_canonical(space, raw, canonical)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: raw(:, :)
      complex(rp), intent(out) :: canonical(:, :)
      integer :: column

      call require_operator(space, raw, 'response_raw_to_canonical: raw shape mismatch')
      call require_operator(space, canonical, 'response_raw_to_canonical: result shape mismatch')
      do column = 1, space%ndim
         canonical(:, column) = raw(:, column)*space%metric_weights(column)
      end do
   end subroutine response_raw_to_canonical

   !> Convert canonical B back to the finite raw kernel on its active columns.
   !> A nonzero origin column has no finite pointwise-kernel representation and
   !> is rejected rather than hidden with a regularizer.
   module subroutine response_canonical_to_raw(space, canonical, raw)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: canonical(:, :)
      complex(rp), intent(out) :: raw(:, :)
      integer :: column

      call require_operator(space, canonical, 'response_canonical_to_raw: canonical shape mismatch')
      call require_operator(space, raw, 'response_canonical_to_raw: result shape mismatch')
      raw = cmplx(0.0_rp, 0.0_rp, rp)
      do column = 1, space%ndim
         if (space%metric_weights(column) > 0.0_rp) then
            raw(:, column) = canonical(:, column)/space%metric_weights(column)
         elseif (any(canonical(:, column) /= cmplx(0.0_rp, 0.0_rp, rp))) then
            error stop 'response_canonical_to_raw: nonzero null-origin column is not representable'
         end if
      end do
   end subroutine response_canonical_to_raw

   subroutine require_layout(space)
      type(response_space_layout), intent(in) :: space

      if (space%ndim < 1 .or. .not. allocated(space%metric_weights) .or. &
          size(space%metric_weights) /= space%ndim) then
         error stop 'LR response space: layout is not initialized'
      end if
   end subroutine require_layout

   subroutine require_vector(space, vector, message)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: vector(:)
      character(len=*), intent(in) :: message

      call require_layout(space)
   end subroutine require_vector

   subroutine require_operator(space, operator, message)
      type(response_space_layout), intent(in) :: space
      complex(rp), intent(in) :: operator(:, :)
      character(len=*), intent(in) :: message

      call require_layout(space)
   end subroutine require_operator


   ! --- from lr_lmto_endpoint_branches_mod ---

   module pure logical function lmto_product_branch_valid(branch) result(valid)
      integer, intent(in) :: branch

      valid = branch >= lmto_product_branch_00 .and. branch <= lmto_product_branch_02
   end function lmto_product_branch_valid

   module pure subroutine lmto_product_branch_powers(branch, p, q)
      integer, intent(in) :: branch
      integer, intent(out) :: p, q

      p = -1
      q = -1
      select case (branch)
      case (lmto_product_branch_00)
         p = 0; q = 0
      case (lmto_product_branch_10)
         p = 1; q = 0
      case (lmto_product_branch_01)
         p = 0; q = 1
      case (lmto_product_branch_11)
         p = 1; q = 1
      case (lmto_product_branch_20)
         p = 2; q = 0
      case (lmto_product_branch_02)
         p = 0; q = 2
      end select
   end subroutine lmto_product_branch_powers

   module pure character(len=2) function lmto_product_branch_label(branch) result(label)
      integer, intent(in) :: branch

      select case (branch)
      case (lmto_product_branch_00); label = '00'
      case (lmto_product_branch_10); label = '10'
      case (lmto_product_branch_01); label = '01'
      case (lmto_product_branch_11); label = '11'
      case (lmto_product_branch_20); label = '20'
      case (lmto_product_branch_02); label = '02'
      case default; label = '??'
      end select
   end function lmto_product_branch_label

   !> Stable real energy power for the endpoint and GF contracts.
   module pure real(rp) function lmto_product_energy_power(energy, power) result(value)
      real(rp), intent(in) :: energy
      integer, intent(in) :: power
      integer :: n

      if (power < 0 .or. power > lmto_product_max_gf_moment) then
         error stop 'lmto_product_energy_power: unsupported power'
      end if
      value = 1.0_rp
      do n = 1, power
         value = value*energy
      end do
   end function lmto_product_energy_power

   !> Apply one endpoint branch in coefficient space.
   !>
   !> For a matrix V whose band-state matrix element is the authoritative
   !> quantity
   !>
   !>   <n|delta H_pq|m> = E_n**p E_m**q <n|V_pq|m>,
   !>
   !> the direct action is H**p V H**q.  A stored-dual source is the
   !> transpose/adjoint of the measurement vertex; its endpoint labels remain
   !> attached to the original measurement endpoints, so its coefficient-space
   !> action is the adjoint-equivalent H**q V H**p.  Keeping both conventions
   !> here prevents callers from duplicating, and silently truncating, the
   !> Hamiltonian-power algebra.
   module subroutine lmto_product_apply_branch_action(hamiltonian, matrix, branch, action, stored_dual)
      complex(rp), intent(in) :: hamiltonian(:, :), matrix(:, :)
      integer, intent(in) :: branch
      complex(rp), intent(out) :: action(:, :)
      logical, intent(in), optional :: stored_dual
      integer :: p, q, left_power, right_power, power
      logical :: use_stored_dual

      if (.not. lmto_product_branch_valid(branch)) then
         error stop 'lmto_product_apply_branch_action: invalid branch'
      end if

      call lmto_product_branch_powers(branch, p, q)
      use_stored_dual = .false.
      if (present(stored_dual)) use_stored_dual = stored_dual
      if (use_stored_dual) then
         left_power = q
         right_power = p
      else
         left_power = p
         right_power = q
      end if

      action = matrix
      do power = 1, left_power
         action = matmul(hamiltonian, action)
      end do
      do power = 1, right_power
         action = matmul(action, hamiltonian)
      end do
   end subroutine lmto_product_apply_branch_action

   !> Evaluate one radial coefficient B_pq of the certified second-order
   !> absolute-energy product polynomial.  The expression is
   !>
   !>   R_L R_R |_2 = phi_L phi_R
   !>       + dL dotphi_L phi_R + phi_L dR dotphi_R
   !>       + dL dR dotphi_L dotphi_R
   !>       + 1/2 dL^2 phiddot_L phi_R
   !>       + 1/2 dR^2 phi_L phiddot_R,
   !>
   !> with dL=E_L-Enu_L and dR=E_R-Enu_R.  Expanding in absolute energies
   !> gives exactly the six branches 00/10/01/11/20/02.
   module recursive subroutine lmto_product_second_order_radial_branch(radial, ir, l, lp, spin_left, spin_right, branch, value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, lp, spin_left, spin_right, branch
      real(rp), intent(out) :: value
      real(rp) :: phi_l, phi_r, dot_l, dot_r, ddot_l, ddot_r
      real(rp) :: enu_l, enu_r, positive_value

      if (.not. lmto_product_branch_valid(branch)) error stop &
         'lmto_product_second_order_radial_branch: invalid branch'
      if (ir < 1 .or. ir > radial%npoint .or. l < 0 .or. l > radial%lmax .or. lp < 0 .or. lp > radial%lmax .or. &
          spin_left < 1 .or. spin_left > radial%nspin .or. spin_right < 1 .or. spin_right > radial%nspin) then
         error stop 'lmto_product_second_order_radial_branch: index outside radial basis'
      end if

      if (ir == 1) then
         if (l /= 0 .or. lp /= 0) then
            value = 0.0_rp
            return
         end if
         call lmto_product_second_order_radial_branch(radial, 2, l, lp, spin_left, spin_right, branch, positive_value)
         call lmto_product_second_order_radial_branch(radial, 3, l, lp, spin_left, spin_right, branch, value)
         ! The origin value is the same r^2 extrapolation used by the
         ! historical response-space product evaluator.
         value = (positive_value*radial%rofi(3)**2 - value*radial%rofi(2)**2)/ &
            (radial%rofi(3)**2 - radial%rofi(2)**2)
         return
      end if
      if (radial%rofi(ir) <= tiny(1.0_rp)) error stop &
         'lmto_product_second_order_radial_branch: invalid positive radius'

      phi_l = radial%phi_large(ir, l + 1, spin_left)
      phi_r = radial%phi_large(ir, lp + 1, spin_right)
      dot_l = radial%phidot_large(ir, l + 1, spin_left)
      dot_r = radial%phidot_large(ir, lp + 1, spin_right)
      ddot_l = radial%phiddot_large(ir, l + 1, spin_left)
      ddot_r = radial%phiddot_large(ir, lp + 1, spin_right)
      enu_l = radial%enu_work(l + 1, spin_left)
      enu_r = radial%enu_work(lp + 1, spin_right)

      select case (branch)
      case (lmto_product_branch_00)
         value = phi_l*phi_r - enu_l*dot_l*phi_r - enu_r*phi_l*dot_r + enu_l*enu_r*dot_l*dot_r + &
            0.5_rp*enu_l**2*ddot_l*phi_r + 0.5_rp*enu_r**2*phi_l*ddot_r
      case (lmto_product_branch_10)
         value = dot_l*phi_r - enu_r*dot_l*dot_r - enu_l*ddot_l*phi_r
      case (lmto_product_branch_01)
         value = phi_l*dot_r - enu_l*dot_l*dot_r - enu_r*phi_l*ddot_r
      case (lmto_product_branch_11)
         value = dot_l*dot_r
      case (lmto_product_branch_20)
         value = 0.5_rp*ddot_l*phi_r
      case (lmto_product_branch_02)
         value = 0.5_rp*phi_l*ddot_r
      end select
      value = value/radial%rofi(ir)**2
   end subroutine lmto_product_second_order_radial_branch


   ! --- from lr_pauli_transition_vertex_mod ---

   module subroutine pauli_endpoint_state_initialize(this, energy, coefficients)
      class(pauli_endpoint_state), intent(out) :: this
      real(rp), intent(in) :: energy
      complex(rp), intent(in) :: coefficients(:)

      this%energy = energy
      allocate(this%coefficients(size(coefficients)))
      this%coefficients = coefficients
   end subroutine pauli_endpoint_state_initialize

   module pure function pauli_charge_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      matrix(2, 2) = cmplx(1.0_rp, 0.0_rp, rp)
   end function pauli_charge_matrix

   module pure function pauli_sigma_x_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(1, 2) = cmplx(1.0_rp, 0.0_rp, rp)
      matrix(2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
   end function pauli_sigma_x_matrix

   module pure function pauli_sigma_y_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(1, 2) = cmplx(0.0_rp, -1.0_rp, rp)
      matrix(2, 1) = cmplx(0.0_rp, 1.0_rp, rp)
   end function pauli_sigma_y_matrix

   module pure function pauli_sigma_z_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(1, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      matrix(2, 2) = cmplx(-1.0_rp, 0.0_rp, rp)
   end function pauli_sigma_z_matrix

   !> The typed circular convention is encoded by matrices, not an integer
   !> selector: sigma_plus=(sigma_x+i sigma_y)/2 and sigma_minus=(sigma_x-i sigma_y)/2.
   module pure function pauli_sigma_plus_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(1, 2) = cmplx(1.0_rp, 0.0_rp, rp)
   end function pauli_sigma_plus_matrix

   module pure function pauli_sigma_minus_matrix() result(matrix)
      complex(rp) :: matrix(2, 2)

      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(2, 1) = cmplx(1.0_rp, 0.0_rp, rp)
   end function pauli_sigma_minus_matrix

   !> Evaluate one operator channel.  The response layout must have one
   !> channel; use the rank-3 overload for a vector containing several explicit
   !> Pauli operators.
   module subroutine evaluate_pauli_transition_vertex_one(space, radial_bases, left_state, right_state, operator_matrix, &
                                                   capabilities, transition_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(in) :: operator_matrix(:, :)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      complex(rp), intent(out) :: transition_vector(:)
      complex(rp), allocatable :: operator_channels(:, :, :)

      if (space%nchannel /= 1) then
         error stop 'evaluate_pauli_transition_vertex: one operator requires nchannel=1'
      end if

      allocate(operator_channels(2, 2, 1))
      operator_channels(:, :, 1) = operator_matrix
      call evaluate_pauli_transition_vertex_channels(space, radial_bases, left_state, right_state, operator_channels, &
                                                      capabilities, transition_vector)
   end subroutine evaluate_pauli_transition_vertex_one

   !> Evaluate one or more explicit 2x2 Pauli-space operators.  The third
   !> operator dimension is the LR-04 response channel, so the complete
   !> `(site,L,M,radial,mu)` vector is emitted in one call.
   module subroutine evaluate_pauli_transition_vertex_channels(space, radial_bases, left_state, right_state, operator_matrices, &
                                                        capabilities, transition_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(in) :: operator_matrices(:, :, :)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      complex(rp), intent(out) :: transition_vector(:)

      integer :: nsite, norb, npoint, nresponse, nchannel
      integer :: isite, ir, response_l, response_m, channel, iorb, jorb, ispin, jspin, branch, p, q
      integer :: orbital_l, orbital_lp, left_offset, right_offset, flat
      type(response_super_index) :: item
      complex(rp) :: value, coefficient_product
      real(rp) :: radial_product

      call validate_vertex_inputs(space, radial_bases, left_state, right_state, operator_matrices, capabilities, transition_vector)

      nsite = space%nsite
      norb = (radial_bases(1)%lmax + 1)**2
      npoint = space%npoint
      nresponse = space%response_lmax
      nchannel = space%nchannel
      transition_vector = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, nsite
         left_offset = (isite - 1)*2*norb
         right_offset = left_offset
         do response_l = 0, nresponse
            do response_m = -response_l, response_l
               do ir = 1, npoint
                  do channel = 1, nchannel
                     value = cmplx(0.0_rp, 0.0_rp, rp)
                     do iorb = 1, norb
                        orbital_l = lmto_orbital_l(iorb)
                       do jorb = 1, norb
                           orbital_lp = lmto_orbital_l(jorb)
                           do ispin = 1, 2
                              do jspin = 1, 2
                                 coefficient_product = conjg(left_state%coefficients(left_offset + &
                                    (ispin - 1)*norb + iorb))*right_state%coefficients(right_offset + &
                                    (jspin - 1)*norb + jorb)
                                 do branch = 1, lmto_product_nbranch
                                    call lmto_product_branch_powers(branch, p, q)
                                    call lmto_product_second_order_radial_branch(radial_bases(isite), ir, orbital_l, &
                                       orbital_lp, ispin, jspin, branch, radial_product)
                                    value = value + coefficient_product*operator_matrices(ispin, jspin, channel)* &
                                       lmto_product_energy_power(left_state%energy, p)* &
                                       lmto_product_energy_power(right_state%energy, q)*radial_product* &
                                       response_gaunt(orbital_l, orbital_m(iorb, orbital_l), orbital_lp, &
                                       orbital_m(jorb, orbital_lp), response_l, response_m)
                                 end do
                              end do
                           end do
                        end do
                     end do
                     item = response_super_index(isite, response_l, response_m, ir, channel)
                     call response_flatten_superindex(item, space%nsite, space%response_lmax, space%npoint, &
                                                      space%nchannel, flat)
                     transition_vector(flat) = value
                  end do
               end do
            end do
         end do
      end do

   end subroutine evaluate_pauli_transition_vertex_channels

   subroutine validate_vertex_inputs(space, radial_bases, left_state, right_state, operator_matrices, capabilities, transition_vector)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(in) :: operator_matrices(:, :, :)
      type(pauli_vertex_capabilities), intent(in) :: capabilities
      complex(rp), intent(out) :: transition_vector(:)

      integer :: isite, lmax, norb
      real(rp), parameter :: mesh_tolerance = 2.0e-13_rp

      ! This call is deliberately inside the public evaluator validation
      ! path.  Unsupported representations cannot be enabled by a caller-side
      ! convention or by accidentally supplying compatible-sized arrays.
      if (size(radial_bases) /= space%nsite) then
         error stop 'evaluate_pauli_transition_vertex: one radial basis is required per response site'
      end if
      do isite = 1, size(radial_bases)
         call radial_bases(isite)%require_supported(capabilities%reciprocal_mode, capabilities%hamiltonian_order, &
            capabilities%orthogonal, capabilities%collinear, capabilities%has_soc, capabilities%has_extra_operator)
      end do

      if (size(operator_matrices, 1) /= 2 .or. size(operator_matrices, 2) /= 2 .or. &
          size(operator_matrices, 3) /= space%nchannel) then
         error stop 'evaluate_pauli_transition_vertex: operator/channel shape mismatch'
      end if
      if (size(transition_vector) /= space%ndim) then
         error stop 'evaluate_pauli_transition_vertex: response vector shape mismatch'
      end if

      lmax = radial_bases(1)%lmax
      norb = (lmax + 1)**2
      if (space%response_lmax < 2*lmax) then
         error stop 'evaluate_pauli_transition_vertex: response angular product space is incomplete'
      end if
      if (.not. allocated(space%radius) .or. size(space%radius) /= radial_bases(1)%npoint .or. &
          space%npoint /= radial_bases(1)%npoint) then
         error stop 'evaluate_pauli_transition_vertex: response/radial mesh shape mismatch'
      end if
      do isite = 1, size(radial_bases)
         if (radial_bases(isite)%lmax /= lmax .or. radial_bases(isite)%nspin /= 2 .or. &
             radial_bases(isite)%npoint /= space%npoint .or. .not. allocated(radial_bases(isite)%rofi) .or. &
             .not. allocated(radial_bases(isite)%channel_present) .or. &
             .not. all(radial_bases(isite)%channel_present)) then
            error stop 'evaluate_pauli_transition_vertex: incomplete or inconsistent Pauli radial basis'
         end if
         if (abs(radial_bases(isite)%mesh_a - space%a) > mesh_tolerance*max(1.0_rp, abs(space%a)) .or. &
             abs(radial_bases(isite)%mesh_b - space%b) > mesh_tolerance*max(1.0_rp, abs(space%b)) .or. &
             maxval(abs(radial_bases(isite)%rofi - space%radius)) > mesh_tolerance* &
             max(1.0_rp, maxval(abs(space%radius)))) then
            error stop 'evaluate_pauli_transition_vertex: radial mesh provenance mismatch'
         end if
      end do

      if (.not. allocated(left_state%coefficients) .or. .not. allocated(right_state%coefficients) .or. &
          size(left_state%coefficients) /= 2*norb*space%nsite .or. &
          size(right_state%coefficients) /= 2*norb*space%nsite) then
         error stop 'evaluate_pauli_transition_vertex: endpoint coefficient shape mismatch'
      end if
   end subroutine validate_vertex_inputs

   pure integer function orbital_m(iorb, l) result(m)
      integer, intent(in) :: iorb, l
      m = iorb - l*l - l - 1
   end function orbital_m


   ! --- from lr_lmto_product_response_basis_mod ---

   module pure integer function lmto_product_candidate_count(lmax, response_l) result(count)
      integer, intent(in) :: lmax, response_l
      integer :: l, lp

      count = 0
      do l = 0, lmax
         do lp = 0, lmax
            if (allowed_pair(l, lp, response_l)) count = count + lmto_product_nbranch
         end do
      end do
   end function lmto_product_candidate_count

   module subroutine lmto_enumerate_product_candidates(lmax, response_l, candidates)
      integer, intent(in) :: lmax, response_l
      type(lmto_product_candidate), allocatable, intent(out) :: candidates(:)
      integer :: l, lp, p, q, branch, k

      allocate(candidates(lmto_product_candidate_count(lmax, response_l)))
      k = 0
      do l = 0, lmax
         do lp = 0, lmax
            if (.not. allowed_pair(l, lp, response_l)) cycle
            do branch = lmto_product_branch_00, lmto_product_branch_02
               call lmto_product_branch_powers(branch, p, q)
               k = k + 1
               candidates(k) = lmto_product_candidate(l, lp, p, q, branch)
            end do
         end do
      end do
   end subroutine lmto_enumerate_product_candidates

   pure logical function allowed_pair(l, lp, response_l) result(is_allowed)
      integer, intent(in) :: l, lp, response_l

      is_allowed = abs(l - lp) <= response_l .and. response_l <= l + lp .and. &
         mod(l + lp + response_l, 2) == 0
   end function allowed_pair

   module subroutine lmto_product_response_basis_initialize(this, space, radial_bases, circular_channel, strict_rank)
      class(lmto_product_response_basis), intent(out) :: this
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      integer, intent(in) :: circular_channel
      logical, intent(in), optional :: strict_rank
      logical :: strict
      integer :: isite, response_l

      strict = .true.
      if (present(strict_rank)) strict = strict_rank
      call validate_initialization_inputs(space, radial_bases, circular_channel)

      this%nsite = space%nsite
      this%orbital_lmax = radial_bases(1)%lmax
      this%response_lmax = space%response_lmax
      this%circular_channel = circular_channel
      this%npoint = space%npoint
      this%unpruned_dimension = 0
      this%product_dimension = 0
      allocate(this%blocks(this%nsite, 0:this%response_lmax))

      do isite = 1, this%nsite
         do response_l = 0, this%response_lmax
            this%unpruned_dimension = this%unpruned_dimension + &
               (2*response_l + 1)*lmto_product_candidate_count(this%orbital_lmax, response_l)
            call build_product_block(this%blocks(isite, response_l), space, radial_bases(isite), isite, &
               response_l, circular_channel, strict)
            this%product_dimension = this%product_dimension + &
               (2*response_l + 1)*this%blocks(isite, response_l)%rank
         end do
      end do
   end subroutine lmto_product_response_basis_initialize

   subroutine validate_initialization_inputs(space, radial_bases, circular_channel)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      integer, intent(in) :: circular_channel
      integer :: isite, lmax
      real(rp), parameter :: mesh_tolerance = 2.0e-13_rp

      if (space%nsite < 1 .or. space%response_lmax < 0 .or. space%nchannel /= 1 .or. &
          .not. allocated(space%radius) .or. .not. allocated(space%radial_weights)) then
         error stop 'lmto_product_response_basis: only one explicit response channel is supported'
      end if
      if (size(radial_bases) /= space%nsite) then
         error stop 'lmto_product_response_basis: one radial basis is required per site'
      end if
      if (circular_channel /= lmto_product_channel_plus .and. circular_channel /= lmto_product_channel_minus) then
         error stop 'lmto_product_response_basis: invalid circular channel'
      end if
      lmax = radial_bases(1)%lmax
      if (lmax < 0 .or. lmax > 2 .or. space%response_lmax > 2*lmax) then
         error stop 'lmto_product_response_basis: unsupported angular product space'
      end if
      if (space%npoint < 3 .or. size(space%radius) /= space%npoint .or. &
          size(space%radial_weights) /= space%npoint) then
         error stop 'lmto_product_response_basis: inconsistent radial layout'
      end if
      do isite = 1, space%nsite
         call radial_bases(isite)%require_supported('ham_only', 'second', .true., .true., .false., .false.)
         if (radial_bases(isite)%lmax /= lmax .or. radial_bases(isite)%nspin /= 2 .or. &
             radial_bases(isite)%npoint /= space%npoint .or. .not. allocated(radial_bases(isite)%rofi) .or. &
             .not. allocated(radial_bases(isite)%channel_present) .or. .not. all(radial_bases(isite)%channel_present)) then
            error stop 'lmto_product_response_basis: incomplete radial basis'
         end if
         if (abs(radial_bases(isite)%mesh_a - space%a) > mesh_tolerance*max(1.0_rp, abs(space%a)) .or. &
             abs(radial_bases(isite)%mesh_b - space%b) > mesh_tolerance*max(1.0_rp, abs(space%b)) .or. &
             maxval(abs(radial_bases(isite)%rofi - space%radius)) > mesh_tolerance* &
                max(1.0_rp, maxval(abs(space%radius)))) then
            error stop 'lmto_product_response_basis: radial mesh provenance mismatch'
         end if
      end do
      if (any(.not. ieee_is_finite(space%radial_weights)) .or. any(space%radial_weights < 0.0_rp)) then
         error stop 'lmto_product_response_basis: invalid LR-04 radial weights'
      end if
   end subroutine validate_initialization_inputs

   subroutine build_product_block(block, space, radial, site, response_l, circular_channel, strict)
      type(lmto_product_block), intent(out) :: block
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: site, response_l, circular_channel
      logical, intent(in) :: strict
      type(lmto_product_candidate), allocatable :: candidates(:)
      complex(rp), allocatable :: a(:, :), u(:, :), vt(:, :), work(:)
      complex(rp) :: work_query(1)
      real(rp), allocatable :: column_norms(:), singular_values(:), rwork(:)
      real(rp) :: tau1, tau10, tau100
      integer :: nr, ncandidate, nsv, lwork, info, ir, k, rank1, rank10, rank100

      nr = space%npoint
      call lmto_enumerate_product_candidates(radial%lmax, response_l, candidates)
      ncandidate = size(candidates)
      nsv = min(nr, ncandidate)
      allocate(a(nr, ncandidate), column_norms(ncandidate), singular_values(nsv), &
         u(nr, nsv), vt(nsv, ncandidate), rwork(max(1, 5*nsv)))
      a = cmplx(0.0_rp, 0.0_rp, rp)
      do k = 1, ncandidate
         do ir = 1, nr
            a(ir, k) = cmplx(radial_product(radial, ir, candidates(k)%l, candidates(k)%lp, &
               circular_channel, candidates(k)%branch), 0.0_rp, rp)
         end do
         column_norms(k) = sqrt(sum(space%radial_weights*real(a(:, k)*conjg(a(:, k)), rp)))
         if (.not. ieee_is_finite(column_norms(k)) .or. column_norms(k) <= tiny(1.0_rp)) then
            write (*, '(a,i0,a,i0,a,a,a,i0,a,i0,a,i0,a,i0,a,a,a)') &
               'LMTO product candidate: site=', site, ' response L=', response_l, ' M=all(-L:L) channel=', &
               merge('plus ', 'minus', circular_channel == lmto_product_channel_plus), ' candidate=', k, ' l=', &
               candidates(k)%l, ' l''=', candidates(k)%lp, ' branch=', candidates(k)%branch, ' (', &
               lmto_product_branch_label(candidates(k)%branch), ')'
            write (*, '(a,es16.8,a,es16.8)') '  candidate norm=', column_norms(k), &
               ' largest weighted scale=', maxval(abs(sqrt(space%radial_weights)*a(:, k)))
            error stop 'lmto_product_response_basis: zero or non-finite candidate norm'
         end if
         a(:, k) = sqrt(space%radial_weights)*a(:, k)/column_norms(k)
      end do

      ! Direct SVD of W^(1/2) Btilde.  The candidate Gram matrix is not formed.
      call zgesvd('S', 'S', nr, ncandidate, a, nr, singular_values, u, nr, vt, nsv, work_query, -1, rwork, info)
      if (info /= 0) then
         call report_product_svd_diagnostics(site, response_l, circular_channel, nr, ncandidate, column_norms, &
            singular_values, .false., info)
         error stop 'lmto_product_response_basis: zgesvd workspace query failed'
      end if
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      call zgesvd('S', 'S', nr, ncandidate, a, nr, singular_values, u, nr, vt, nsv, work, lwork, rwork, info)
      deallocate(work)
      if (info /= 0) then
         call report_product_svd_diagnostics(site, response_l, circular_channel, nr, ncandidate, column_norms, &
            singular_values, .false., info)
         error stop 'lmto_product_response_basis: zgesvd failed'
      end if
      if (.not. all(ieee_is_finite(singular_values))) then
         call report_product_svd_diagnostics(site, response_l, circular_channel, nr, ncandidate, column_norms, &
            singular_values, .false., info)
         error stop 'lmto_product_response_basis: non-finite singular value'
      end if

      tau1 = real(max(nr, ncandidate), rp)*epsilon(1.0_rp)*singular_values(1)
      tau10 = 10.0_rp*tau1
      tau100 = 100.0_rp*tau1
      rank1 = count(singular_values > tau1)
      rank10 = count(singular_values > tau10)
      rank100 = count(singular_values > tau100)
      if (rank1 < 1) error stop 'lmto_product_response_basis: no retained numerical product mode'
      if (strict .and. (rank1 /= rank10 .or. rank1 /= rank100)) then
         call report_product_svd_diagnostics(site, response_l, circular_channel, nr, ncandidate, column_norms, &
            singular_values, .true., info)
         write (*, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'LMTO product rank sensitivity site=', site, ' L=', response_l, &
            ' tau1/tau10/tau100=', rank1, '/', rank10, '/', rank100
         error stop 'lmto_product_response_basis: rank sensitivity requires review'
      end if

      block%site = site
      block%response_l = response_l
      block%channel = circular_channel
      block%npoint = nr
      block%ncandidate = ncandidate
      block%nsv = nsv
      block%rank = rank1
      block%rank_tau1 = rank1
      block%rank_tau10 = rank10
      block%rank_tau100 = rank100
      block%tau1 = tau1
      block%tau10 = tau10
      block%tau100 = tau100
      block%rank_stable = rank1 == rank10 .and. rank1 == rank100
      allocate(block%candidates(ncandidate), block%column_norms(ncandidate), block%singular_values(nsv), &
         block%weighted_modes(nr, rank1), block%right_modes(rank1, ncandidate), block%forward_transform(rank1, ncandidate))
      block%candidates = candidates
      block%column_norms = column_norms
      block%singular_values = singular_values
      block%weighted_modes = u(:, 1:rank1)
      block%right_modes = vt(1:rank1, :)
      do k = 1, ncandidate
         block%forward_transform(:, k) = singular_values(1:rank1)*vt(1:rank1, k)*column_norms(k)
      end do
      deallocate(candidates, a, u, vt, column_norms, singular_values, rwork)
   end subroutine build_product_block

   subroutine report_product_svd_diagnostics(site, response_l, circular_channel, nr, ncandidate, column_norms, &
                                             singular_values, spectrum_available, info)
      integer, intent(in) :: site, response_l, circular_channel, nr, ncandidate, info
      real(rp), intent(in) :: column_norms(:), singular_values(:)
      logical, intent(in) :: spectrum_available
      real(rp) :: minimum_singular, condition_estimate
      character(len=5) :: channel_label

      if (circular_channel == lmto_product_channel_plus) then
         channel_label = 'plus'
      else
         channel_label = 'minus'
      end if
      write (*, '(a,i0,a,i0,a,i0,a,i0,a,a,a,i0)') 'LMTO product SVD: matrix=', nr, 'x', ncandidate, &
         ' site=', site, ' L=', response_l, ' channel=', trim(channel_label), ' LAPACK info=', info
      write (*, '(a,es16.8)') '  max candidate norm=', maxval(column_norms)
      if (spectrum_available) then
         minimum_singular = minval(singular_values)
         if (minimum_singular > tiny(1.0_rp)) then
            condition_estimate = singular_values(1)/minimum_singular
         else
            condition_estimate = huge(1.0_rp)
         end if
         write (*, '(a,es16.8,a,es16.8)') '  smallest resolved singular value=', minimum_singular, &
            ' condition estimate=', condition_estimate
      else
         write (*, '(a)') '  smallest resolved singular value=unavailable condition estimate=unavailable'
      end if
   end subroutine report_product_svd_diagnostics

   module integer function lmto_product_flat_index(this, site, response_l, response_m, product_mode) result(flat)
      class(lmto_product_response_basis), intent(in) :: this
      integer, intent(in) :: site, response_l, response_m, product_mode
      integer :: isite, l, rank

      call validate_product_index_inputs(this, site, response_l, response_m, product_mode)
      flat = 0
      do isite = 1, site - 1
         do l = 0, this%response_lmax
            flat = flat + (2*l + 1)*this%blocks(isite, l)%rank
         end do
      end do
      do l = 0, response_l - 1
         flat = flat + (2*l + 1)*this%blocks(site, l)%rank
      end do
      rank = this%blocks(site, response_l)%rank
      flat = flat + (response_m + response_l)*rank + product_mode
   end function lmto_product_flat_index

   module subroutine lmto_product_unflatten_index(this, flat, site, response_l, response_m, product_mode)
      class(lmto_product_response_basis), intent(in) :: this
      integer, intent(in) :: flat
      integer, intent(out) :: site, response_l, response_m, product_mode
      integer :: remaining, block_size, isite, l

      if (flat < 1 .or. flat > this%product_dimension) then
         error stop 'lmto_product_response_basis: product flat index outside representation'
      end if
      remaining = flat
      do isite = 1, this%nsite
         do l = 0, this%response_lmax
            block_size = (2*l + 1)*this%blocks(isite, l)%rank
            if (remaining <= block_size) then
               site = isite
               response_l = l
               response_m = (remaining - 1)/this%blocks(isite, l)%rank - l
               product_mode = mod(remaining - 1, this%blocks(isite, l)%rank) + 1
               return
            end if
            remaining = remaining - block_size
         end do
      end do
      error stop 'lmto_product_response_basis: product unflatten failed'
   end subroutine lmto_product_unflatten_index

   subroutine validate_product_index_inputs(this, site, response_l, response_m, product_mode)
      class(lmto_product_response_basis), intent(in) :: this
      integer, intent(in) :: site, response_l, response_m, product_mode

      if (.not. allocated(this%blocks) .or. site < 1 .or. site > this%nsite .or. &
          response_l < 0 .or. response_l > this%response_lmax .or. abs(response_m) > response_l .or. &
          product_mode < 1 .or. product_mode > this%blocks(site, response_l)%rank) then
         error stop 'lmto_product_response_basis: product index outside representation'
      end if
   end subroutine validate_product_index_inputs

   module subroutine lmto_product_candidate_coefficients(this, site, response_l, response_m, left_state, right_state, coefficients, &
                                                  selected_l)
      class(lmto_product_response_basis), intent(in) :: this
      integer, intent(in) :: site, response_l, response_m
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: coefficients(:)
      logical, intent(in), optional :: selected_l(0:)
      integer :: norb, offset, iorb, jorb, spin_left, spin_right, k, ncandidate
      integer :: orbital_l, orbital_lp, orbital_m, orbital_mp
      real(rp) :: left_power, right_power
      integer :: branch

      if (site < 1 .or. site > this%nsite .or. response_l < 0 .or. response_l > this%response_lmax .or. &
          abs(response_m) > response_l) error stop 'lmto_product_response_basis: invalid candidate block index'
      ncandidate = this%blocks(site, response_l)%ncandidate
      norb = (this%orbital_lmax + 1)**2
      offset = (site - 1)*2*norb
      if (this%circular_channel == lmto_product_channel_plus) then
         spin_left = 1
         spin_right = 2
      else
         spin_left = 2
         spin_right = 1
      end if
      coefficients = cmplx(0.0_rp, 0.0_rp, rp)
      do k = 1, ncandidate
         branch = this%blocks(site, response_l)%candidates(k)%branch
         if (.not. lmto_product_branch_valid(branch)) then
            error stop 'lmto_product_response_basis: candidate has invalid endpoint branch'
         end if
         if (present(selected_l)) then
            if (this%blocks(site, response_l)%candidates(k)%l < lbound(selected_l, 1) .or. &
                this%blocks(site, response_l)%candidates(k)%l > ubound(selected_l, 1) .or. &
                this%blocks(site, response_l)%candidates(k)%lp < lbound(selected_l, 1) .or. &
                this%blocks(site, response_l)%candidates(k)%lp > ubound(selected_l, 1) .or. &
                .not. selected_l(this%blocks(site, response_l)%candidates(k)%l) .or. &
                .not. selected_l(this%blocks(site, response_l)%candidates(k)%lp)) cycle
         end if
         left_power = lmto_product_energy_power(left_state%energy, &
            this%blocks(site, response_l)%candidates(k)%p)
         right_power = lmto_product_energy_power(right_state%energy, &
            this%blocks(site, response_l)%candidates(k)%q)
         do iorb = 1, norb
            orbital_l = lmto_orbital_l(iorb)
            if (orbital_l /= this%blocks(site, response_l)%candidates(k)%l) cycle
            orbital_m = iorb - orbital_l*orbital_l - orbital_l - 1
            do jorb = 1, norb
               orbital_lp = lmto_orbital_l(jorb)
               if (orbital_lp /= this%blocks(site, response_l)%candidates(k)%lp) cycle
               orbital_mp = jorb - orbital_lp*orbital_lp - orbital_lp - 1
               coefficients(k) = coefficients(k) + left_power*right_power* &
                  conjg(left_state%coefficients(offset + (spin_left - 1)*norb + iorb))* &
                  right_state%coefficients(offset + (spin_right - 1)*norb + jorb)* &
                  response_gaunt(orbital_l, orbital_m, orbital_lp, orbital_mp, response_l, response_m)
            end do
         end do
      end do
   end subroutine lmto_product_candidate_coefficients

   module subroutine lmto_product_transition_coordinates(this, left_state, right_state, coordinates, selected_l, response_l_filter)
      class(lmto_product_response_basis), intent(in) :: this
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: coordinates(:)
      logical, intent(in), optional :: selected_l(0:)
      integer, intent(in), optional :: response_l_filter
      complex(rp), allocatable :: coefficients(:)
      integer :: site, response_l, response_m

      coordinates = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         do response_l = 0, this%response_lmax
            if (present(response_l_filter)) then
               if (response_l /= response_l_filter) cycle
            end if
            allocate(coefficients(this%blocks(site, response_l)%ncandidate))
            do response_m = -response_l, response_l
               if (present(selected_l)) then
                  call this%candidate_coefficients(site, response_l, response_m, left_state, right_state, coefficients, &
                     selected_l)
               else
                  call this%candidate_coefficients(site, response_l, response_m, left_state, right_state, coefficients)
               end if
               coordinates(this%flat_index(site, response_l, response_m, 1): &
                  this%flat_index(site, response_l, response_m, this%blocks(site, response_l)%rank)) = &
                  matmul(this%blocks(site, response_l)%forward_transform, coefficients)
            end do
            deallocate(coefficients)
         end do
      end do
   end subroutine lmto_product_transition_coordinates

   !> Construct the six energy-moment LMTO vertex components directly in the
   !> retained product representation.  The candidate map already contains
   !> the radial product information, so this routine adds only the orbital
   !> Gaunt factor and the selected circular spin block.
   !>
   !> Component ordering is the authoritative 00/10/01/11/20/02 branch map.
   !> The tensor is indexed as
   !> `(electronic-left, electronic-right, component, product-coordinate)` and
   !> contains no point-space response allocation.
   module subroutine lmto_product_component_vertex_tensor(this, vertices, selected_l)
      class(lmto_product_response_basis), intent(in) :: this
      complex(rp), allocatable, intent(out) :: vertices(:, :, :, :)
      logical, intent(in), optional :: selected_l(0:)

      integer :: norb, nbasis, flat, site, response_l, response_m, product_mode
      integer :: p, q, component, branch, k, iorb, jorb, spin_left, spin_right
      integer :: orbital_l, orbital_lp, orbital_m, orbital_mp, offset
      real(rp) :: gaunt

      if (.not. allocated(this%blocks) .or. this%nsite < 1 .or. this%product_dimension < 1) then
         error stop 'lmto_product_component_vertex_tensor: representation is uninitialized'
      end if
      norb = (this%orbital_lmax + 1)**2
      nbasis = 2*norb*this%nsite
      allocate(vertices(nbasis, nbasis, lmto_product_nbranch, this%product_dimension))
      vertices = cmplx(0.0_rp, 0.0_rp, rp)

      if (this%circular_channel == lmto_product_channel_plus) then
         spin_left = 1
         spin_right = 2
      else if (this%circular_channel == lmto_product_channel_minus) then
         spin_left = 2
         spin_right = 1
      else
         error stop 'lmto_product_component_vertex_tensor: invalid circular channel'
      end if

      do flat = 1, this%product_dimension
         call this%unflatten_index(flat, site, response_l, response_m, product_mode)
         offset = (site - 1)*2*norb
         do branch = lmto_product_branch_00, lmto_product_branch_02
            call lmto_product_branch_powers(branch, p, q)
            component = branch
               do k = 1, this%blocks(site, response_l)%ncandidate
                  if (this%blocks(site, response_l)%candidates(k)%branch /= branch) cycle
                  orbital_l = this%blocks(site, response_l)%candidates(k)%l
                  orbital_lp = this%blocks(site, response_l)%candidates(k)%lp
                  if (present(selected_l)) then
                     if (orbital_l < lbound(selected_l, 1) .or. orbital_l > ubound(selected_l, 1) .or. &
                         orbital_lp < lbound(selected_l, 1) .or. orbital_lp > ubound(selected_l, 1) .or. &
                         .not. selected_l(orbital_l) .or. .not. selected_l(orbital_lp)) cycle
                  end if
                  do iorb = 1, norb
                     if (lmto_orbital_l(iorb) /= orbital_l) cycle
                     orbital_m = iorb - orbital_l*orbital_l - orbital_l - 1
                     do jorb = 1, norb
                        if (lmto_orbital_l(jorb) /= orbital_lp) cycle
                        orbital_mp = jorb - orbital_lp*orbital_lp - orbital_lp - 1
                        gaunt = response_gaunt(orbital_l, orbital_m, orbital_lp, orbital_mp, &
                           response_l, response_m)
                        vertices(offset + (spin_left - 1)*norb + iorb, &
                           offset + (spin_right - 1)*norb + jorb, component, flat) = &
                           vertices(offset + (spin_left - 1)*norb + iorb, &
                           offset + (spin_right - 1)*norb + jorb, component, flat) + &
                           this%blocks(site, response_l)%forward_transform(product_mode, k)*cmplx(gaunt, 0.0_rp, rp)
                     end do
                  end do
               end do
         end do
      end do
   end subroutine lmto_product_component_vertex_tensor

   real(rp) function radial_product(radial, ir, l, lp, circular_channel, branch) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, lp, circular_channel, branch
      integer :: spin_left, spin_right

      if (circular_channel == lmto_product_channel_plus) then
         spin_left = 1
         spin_right = 2
      else
         spin_left = 2
         spin_right = 1
      end if
      call lmto_product_second_order_radial_branch(radial, ir, l, lp, spin_left, spin_right, branch, value)
   end function radial_product


   ! --- from lr_gf_endpoint_augmentation_mod ---

   module pure logical function lr_gf_endpoint_capabilities_supported(this) result(ok)
      class(lr_gf_endpoint_capabilities), intent(in) :: this

      ok = trim(this%reciprocal_mode) == 'ham_only' .and. trim(this%hamiltonian_order) == 'second' .and. &
           this%orthogonal .and. this%collinear .and. (.not. this%has_soc) .and. &
           (.not. this%has_extra_operator) .and. (.not. this%has_hubbard) .and. &
           (.not. this%has_additive_operator)
   end function lr_gf_endpoint_capabilities_supported

   module subroutine lr_gf_endpoint_require_supported(this, caller)
      class(lr_gf_endpoint_capabilities), intent(in) :: this
      character(len=*), intent(in), optional :: caller
      character(len=96) :: message

      if (.not. this%supported()) then
         if (present(caller)) then
            message = trim(caller)//': native-GF endpoint augmentation capability is unsupported'
         else
            message = 'native-GF endpoint augmentation capability is unsupported'
         end if
         error stop trim(message)
      end if
   end subroutine lr_gf_endpoint_require_supported

   module subroutine lr_gf_dense_provider_initialize(this, h_eff, block_size)
      class(lr_gf_dense_hamiltonian_provider), intent(out) :: this
      complex(rp), intent(in) :: h_eff(:, :)
      integer, intent(in) :: block_size

      if (block_size < 1 .or. size(h_eff, 1) /= size(h_eff, 2) .or. &
          mod(size(h_eff, 1), block_size) /= 0) then
         error stop 'lr_gf_dense_provider_initialize: inconsistent matrix/block dimensions'
      end if
      allocate(this%h_eff(size(h_eff, 1), size(h_eff, 2)))
      this%h_eff = h_eff
      this%block_size = block_size
   end subroutine lr_gf_dense_provider_initialize

   module subroutine lr_gf_dense_provider_apply(this, source_site, destination_site, seed, destination_block)
      class(lr_gf_dense_hamiltonian_provider), intent(inout) :: this
      integer, intent(in) :: source_site, destination_site
      complex(rp), intent(in) :: seed(:, :)
      complex(rp), intent(out) :: destination_block(:, :)
      integer :: nsite, source_first, destination_first

      if (.not. allocated(this%h_eff) .or. this%block_size < 1) then
         error stop 'lr_gf_dense_provider_apply: provider is not initialized'
      end if
      nsite = size(this%h_eff, 1)/this%block_size
      if (source_site < 1 .or. source_site > nsite .or. destination_site < 1 .or. &
          destination_site > nsite .or. size(seed, 1) /= this%block_size .or. &
          size(destination_block, 1) /= this%block_size .or. size(destination_block, 2) /= size(seed, 2)) then
         error stop 'lr_gf_dense_provider_apply: invalid site or seed shape'
      end if

      source_first = (source_site - 1)*this%block_size + 1
      destination_first = (destination_site - 1)*this%block_size + 1
      destination_block = matmul(this%h_eff(destination_first:destination_first + this%block_size - 1, &
                                             source_first:source_first + this%block_size - 1), seed)
   end subroutine lr_gf_dense_provider_apply

   !> Promote one native coefficient-space G_ab(z) block to a factored
   !> Pauli/no-SOC endpoint representation.
   !>
   !> hgamma_block is optional only because a provider may obtain the block
   !> from the production H_eff action.  Exactly one of hgamma_block and
   !> hamiltonian_provider must be supplied.  `reverse_coefficient_gf` and
   !> `reverse_hgamma_block` are accepted for pair callers that need G_ba; the
   !> forward result itself uses only the explicitly supplied forward blocks.
   module subroutine augment_lr_gf_endpoint(radial_left, radial_right, left_site, right_site, z, coefficient_gf, &
                                     augmented, capabilities, hgamma_block, hamiltonian_provider, &
                                     reverse_coefficient_gf, reverse_hgamma_block)
      type(lmto_radial_basis), intent(in) :: radial_left, radial_right
      integer, intent(in) :: left_site, right_site
      complex(rp), intent(in) :: z
      complex(rp), intent(in) :: coefficient_gf(:, :)
      type(lr_gf_augmented_block), intent(out) :: augmented
      type(lr_gf_endpoint_capabilities), intent(in), optional :: capabilities
      complex(rp), intent(in), optional :: hgamma_block(:, :)
      class(lr_gf_effective_hamiltonian_provider), intent(inout), optional :: hamiltonian_provider
      complex(rp), intent(in), optional :: reverse_coefficient_gf(:, :)
      complex(rp), intent(in), optional :: reverse_hgamma_block(:, :)

      type(lr_gf_endpoint_capabilities) :: caps
      complex(rp), allocatable :: hgamma(:, :)
      complex(rp), allocatable :: d_left(:, :), d_right(:, :), identity(:, :)
      complex(rp), allocatable :: seed(:, :), effective_block(:, :)
      integer :: nlocal, nspin, norb
      logical :: have_hgamma

      caps = lr_gf_endpoint_capabilities()
      if (present(capabilities)) caps = capabilities
      call caps%require_supported('augment_lr_gf_endpoint')

      call validate_radial_pair(radial_left, radial_right, coefficient_gf)
      if (left_site < 1 .or. right_site < 1) then
         error stop 'augment_lr_gf_endpoint: endpoint site indices must be positive'
      end if
      norb = (radial_left%lmax + 1)**2
      nspin = 2
      nlocal = nspin*norb
      have_hgamma = present(hgamma_block) .or. present(hamiltonian_provider)
      if (.not. have_hgamma) then
         error stop 'augment_lr_gf_endpoint: hgamma_ab requires an explicit block or effective-action provider'
      end if
      allocate(hgamma(nlocal, nlocal), d_left(nlocal, nlocal), d_right(nlocal, nlocal), &
               identity(nlocal, nlocal), seed(nlocal, nlocal), effective_block(nlocal, nlocal))
      if (present(hgamma_block)) then
         hgamma = hgamma_block
      else
         identity = cmplx(0.0_rp, 0.0_rp, rp)
         do nlocal = 1, size(identity, 1)
            identity(nlocal, nlocal) = cmplx(1.0_rp, 0.0_rp, rp)
         end do
         seed = identity
         call hamiltonian_provider%apply(right_site, left_site, seed, effective_block)
         hgamma = effective_block
         if (left_site == right_site) hgamma = hgamma - diagonal_enu(radial_left)
      end if

      d_left = diagonal_difference(radial_left, z)
      d_right = diagonal_difference(radial_right, z)
      identity = cmplx(0.0_rp, 0.0_rp, rp)
      do nlocal = 1, size(identity, 1)
         identity(nlocal, nlocal) = cmplx(1.0_rp, 0.0_rp, rp)
      end do

      augmented%left_site = left_site
      augmented%right_site = right_site
      augmented%norb = norb
      augmented%nspin = nspin
      augmented%z = z
      allocate(augmented%coefficient_gf(nlocal, nlocal), augmented%hgamma_ab(nlocal, nlocal), &
               augmented%h_g(nlocal, nlocal), augmented%g_h(nlocal, nlocal), augmented%h_g_h(nlocal, nlocal))
      augmented%coefficient_gf = coefficient_gf
      augmented%hgamma_ab = hgamma
      ! Keep the contact terms visible.  They are not optional onsite
      ! corrections and they vanish algebraically only for offsite blocks.
      augmented%h_g = matmul(d_left, coefficient_gf)
      augmented%g_h = matmul(coefficient_gf, d_right)
      if (left_site == right_site) then
         augmented%h_g = augmented%h_g - identity
         augmented%g_h = augmented%g_h - identity
      end if
      augmented%h_g_h = matmul(d_left, matmul(coefficient_gf, d_right))
      if (left_site == right_site) augmented%h_g_h = augmented%h_g_h - d_left
      augmented%h_g_h = augmented%h_g_h - hgamma

      call copy_radial_maps(radial_left, radial_right, augmented)
      deallocate(hgamma, d_left, d_right, identity, seed, effective_block)
   end subroutine augment_lr_gf_endpoint

   !> Convenience pair entry point.  It makes the reverse G_ba/hgamma_ba
   !> requirement explicit while retaining the same four-branch representation
   !> for each direction.
   module subroutine augment_lr_gf_endpoint_pair(radial_left, radial_right, left_site, right_site, z, coefficient_gf, &
                                          reverse_coefficient_gf, forward, reverse, capabilities, hgamma_block, &
                                          reverse_hgamma_block, hamiltonian_provider)
      type(lmto_radial_basis), intent(in) :: radial_left, radial_right
      integer, intent(in) :: left_site, right_site
      complex(rp), intent(in) :: z
      complex(rp), intent(in) :: coefficient_gf(:, :), reverse_coefficient_gf(:, :)
      type(lr_gf_augmented_block), intent(out) :: forward, reverse
      type(lr_gf_endpoint_capabilities), intent(in), optional :: capabilities
      complex(rp), intent(in), optional :: hgamma_block(:, :), reverse_hgamma_block(:, :)
      class(lr_gf_effective_hamiltonian_provider), intent(inout), optional :: hamiltonian_provider

      call augment_lr_gf_endpoint(radial_left, radial_right, left_site, right_site, z, coefficient_gf, forward, &
         capabilities, hgamma_block, hamiltonian_provider, reverse_coefficient_gf, reverse_hgamma_block)
      if (present(reverse_hgamma_block)) then
         call augment_lr_gf_endpoint(radial_right, radial_left, right_site, left_site, z, reverse_coefficient_gf, &
            reverse, capabilities, reverse_hgamma_block, reverse_coefficient_gf=coefficient_gf)
      else
         call augment_lr_gf_endpoint(radial_right, radial_left, right_site, left_site, z, reverse_coefficient_gf, &
            reverse, capabilities, hamiltonian_provider=hamiltonian_provider, reverse_coefficient_gf=coefficient_gf)
      end if
   end subroutine augment_lr_gf_endpoint_pair

   subroutine validate_radial_pair(radial_left, radial_right, coefficient_gf)
      type(lmto_radial_basis), intent(in) :: radial_left, radial_right
      complex(rp), intent(in) :: coefficient_gf(:, :)
      integer :: nlocal

      call radial_left%require_supported('ham_only', 'second', .true., .true., .false., .false.)
      call radial_right%require_supported('ham_only', 'second', .true., .true., .false., .false.)
      if (radial_left%nspin /= 2 .or. radial_right%nspin /= 2 .or. radial_left%lmax /= radial_right%lmax) then
         error stop 'augment_lr_gf_endpoint: endpoint radial bases are not a common Pauli sp/spd basis'
      end if
      if (.not. allocated(radial_left%channel_present) .or. .not. allocated(radial_right%channel_present)) then
         error stop 'augment_lr_gf_endpoint: endpoint radial channels are incomplete'
      end if
      if (.not. all(radial_left%channel_present) .or. .not. all(radial_right%channel_present)) then
         error stop 'augment_lr_gf_endpoint: endpoint radial channels are incomplete'
      end if
      nlocal = 2*(radial_left%lmax + 1)**2
      if (any(shape(coefficient_gf) /= [nlocal, nlocal])) then
         error stop 'augment_lr_gf_endpoint: coefficient GF block shape mismatch'
      end if
      if (.not. allocated(radial_left%rofi) .or. .not. allocated(radial_right%rofi)) then
         error stop 'augment_lr_gf_endpoint: endpoint radial meshes are not allocated'
      end if
      if (.not. allocated(radial_left%phi_large) .or. .not. allocated(radial_left%phidot_large) .or. &
          .not. allocated(radial_right%phi_large) .or. .not. allocated(radial_right%phidot_large) .or. &
          .not. allocated(radial_left%enu_work) .or. .not. allocated(radial_right%enu_work)) then
         error stop 'augment_lr_gf_endpoint: endpoint radial augmentation arrays are incomplete'
      end if
   end subroutine validate_radial_pair

   function diagonal_enu(radial) result(enu)
      type(lmto_radial_basis), intent(in) :: radial
      complex(rp) :: enu(2*(radial%lmax + 1)**2, 2*(radial%lmax + 1)**2)
      integer :: iorb, ispin, l, index

      enu = cmplx(0.0_rp, 0.0_rp, rp)
      do ispin = 1, 2
         do iorb = 1, (radial%lmax + 1)**2
            l = lmto_orbital_l(iorb)
            index = (ispin - 1)*(radial%lmax + 1)**2 + iorb
            enu(index, index) = cmplx(radial%enu_work(l + 1, ispin), 0.0_rp, rp)
         end do
      end do
   end function diagonal_enu

   function diagonal_difference(radial, z) result(difference)
      type(lmto_radial_basis), intent(in) :: radial
      complex(rp), intent(in) :: z
      complex(rp) :: difference(2*(radial%lmax + 1)**2, 2*(radial%lmax + 1)**2)
      integer :: i

      difference = -diagonal_enu(radial)
      do i = 1, size(difference, 1)
         difference(i, i) = z + difference(i, i)
      end do
   end function diagonal_difference

   subroutine copy_radial_maps(radial_left, radial_right, augmented)
      type(lmto_radial_basis), intent(in) :: radial_left, radial_right
      type(lr_gf_augmented_block), intent(inout) :: augmented
      integer :: iorb, ispin, l

      allocate(augmented%phi_left(radial_left%npoint, augmented%norb, 2), &
               augmented%phidot_left(radial_left%npoint, augmented%norb, 2), &
               augmented%phi_right(radial_right%npoint, augmented%norb, 2), &
               augmented%phidot_right(radial_right%npoint, augmented%norb, 2), &
               augmented%radius_left(radial_left%npoint), augmented%radius_right(radial_right%npoint))
      augmented%radius_left = radial_left%rofi
      augmented%radius_right = radial_right%rofi
      do ispin = 1, 2
         do iorb = 1, augmented%norb
            l = lmto_orbital_l(iorb)
            augmented%phi_left(:, iorb, ispin) = radial_left%phi_large(:, l + 1, ispin)
            augmented%phidot_left(:, iorb, ispin) = radial_left%phidot_large(:, l + 1, ispin)
            augmented%phi_right(:, iorb, ispin) = radial_right%phi_large(:, l + 1, ispin)
            augmented%phidot_right(:, iorb, ispin) = radial_right%phidot_large(:, l + 1, ispin)
         end do
      end do
   end subroutine copy_radial_maps

   !> Evaluate the full 2x2 Pauli spin GF at two radial/angular endpoints.
   module subroutine lr_gf_augmented_block_evaluate_point(this, left_radial_index, left_theta, left_phi, &
                                                   right_radial_index, right_theta, right_phi, point_gf)
      class(lr_gf_augmented_block), intent(in) :: this
      integer, intent(in) :: left_radial_index, right_radial_index
      real(rp), intent(in) :: left_theta, left_phi, right_theta, right_phi
      complex(rp), intent(out) :: point_gf(:, :)
      complex(rp) :: branches(2, 2, 4)

      call this%branch_at_point(left_radial_index, left_theta, left_phi, right_radial_index, right_theta, right_phi, branches)
      point_gf = sum(branches, dim=3)
   end subroutine lr_gf_augmented_block_evaluate_point

   !> Return the four radial/angular branch matrices separately.  This is the
   !> diagnostic form used by the dense and negative-control oracles.
   module subroutine lr_gf_augmented_block_branch_at_point(this, left_radial_index, left_theta, left_phi, &
                                                    right_radial_index, right_theta, right_phi, branches)
      class(lr_gf_augmented_block), intent(in) :: this
      integer, intent(in) :: left_radial_index, right_radial_index
      real(rp), intent(in) :: left_theta, left_phi, right_theta, right_phi
      complex(rp), intent(out) :: branches(:, :, :)
      complex(rp) :: phi_left(2, 2*this%norb), dot_left(2, 2*this%norb)
      complex(rp) :: phi_right(2, 2*this%norb), dot_right(2, 2*this%norb)
      integer :: iorb, ispin, l, index

      call validate_point_indices(this, left_radial_index, right_radial_index)
      phi_left = cmplx(0.0_rp, 0.0_rp, rp)
      dot_left = cmplx(0.0_rp, 0.0_rp, rp)
      phi_right = cmplx(0.0_rp, 0.0_rp, rp)
      dot_right = cmplx(0.0_rp, 0.0_rp, rp)
      do ispin = 1, 2
         do iorb = 1, this%norb
            l = lmto_orbital_l(iorb)
            index = (ispin - 1)*this%norb + iorb
            phi_left(ispin, index) = cmplx(radial_over_radius(this%phi_left(:, iorb, ispin), &
               this%radius_left, left_radial_index, l), 0.0_rp, rp)* &
               response_harmonic(l, orbital_m(iorb, l), left_theta, left_phi)
            dot_left(ispin, index) = cmplx(radial_over_radius(this%phidot_left(:, iorb, ispin), &
               this%radius_left, left_radial_index, l), 0.0_rp, rp)* &
               response_harmonic(l, orbital_m(iorb, l), left_theta, left_phi)
            phi_right(ispin, index) = cmplx(radial_over_radius(this%phi_right(:, iorb, ispin), &
               this%radius_right, right_radial_index, l), 0.0_rp, rp)* &
               response_harmonic(l, orbital_m(iorb, l), right_theta, right_phi)
            dot_right(ispin, index) = cmplx(radial_over_radius(this%phidot_right(:, iorb, ispin), &
               this%radius_right, right_radial_index, l), 0.0_rp, rp)* &
               response_harmonic(l, orbital_m(iorb, l), right_theta, right_phi)
         end do
      end do

      branches(:, :, 1) = matmul(phi_left, matmul(this%coefficient_gf, conjg(transpose(phi_right))))
      branches(:, :, 2) = matmul(dot_left, matmul(this%h_g, conjg(transpose(phi_right))))
      branches(:, :, 3) = matmul(phi_left, matmul(this%g_h, conjg(transpose(dot_right))))
      branches(:, :, 4) = matmul(dot_left, matmul(this%h_g_h, conjg(transpose(dot_right))))
   end subroutine lr_gf_augmented_block_branch_at_point

   subroutine validate_point_indices(this, left_radial_index, right_radial_index)
      class(lr_gf_augmented_block), intent(in) :: this
      integer, intent(in) :: left_radial_index, right_radial_index

      if (.not. allocated(this%phi_left) .or. .not. allocated(this%phi_right) .or. &
          left_radial_index < 1 .or. left_radial_index > size(this%radius_left) .or. &
          right_radial_index < 1 .or. right_radial_index > size(this%radius_right)) then
         error stop 'lr_gf_augmented_block: invalid or uninitialized radial endpoint'
      end if
   end subroutine validate_point_indices

   function radial_over_radius(values, radius, index, l) result(value)
      real(rp), intent(in) :: values(:), radius(:)
      integer, intent(in) :: index, l
      real(rp) :: value
      real(rp) :: first, second

      if (index /= 1) then
         if (radius(index) <= tiny(1.0_rp)) error stop 'lr_gf endpoint: nonpositive radial mesh point'
         value = values(index)/radius(index)
         return
      end if
      if (l /= 0) then
         value = 0.0_rp
      else
         if (size(radius) < 3 .or. radius(2) <= 0.0_rp .or. radius(3) <= radius(2)) then
            error stop 'lr_gf endpoint: cannot form the regular s-channel origin limit'
         end if
         first = values(2)/radius(2)
         second = values(3)/radius(3)
         value = (first*radius(3) - second*radius(2))/(radius(3) - radius(2))
      end if
   end function radial_over_radius



   ! --- from lr_projected_site_spin_mod ---

   module pure logical function projected_selector_supported(selection) result(ok)
      character(len=*), intent(in) :: selection

      select case (trim(selection))
      case (projected_selection_d, projected_selection_spd)
         ok = .true.
      case default
         ok = .false.
      end select
   end function projected_selector_supported

   module subroutine projected_site_spin_initialize(this, space, radial_bases, selection)
      class(projected_site_spin_contract), intent(out) :: this
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      character(len=*), intent(in) :: selection

      if (.not. projected_selector_supported(selection)) then
         error stop 'DRESP-01: selector must be exactly d or spd; spdf is unsupported'
      end if
      if (space%nsite < 1 .or. space%nchannel /= 1 .or. space%response_lmax < 4 .or. &
          space%npoint < 3 .or. .not. allocated(space%radius) .or. &
          .not. allocated(space%radial_weights)) then
         error stop 'DRESP-01: response space is not a complete spd one-channel layout'
      end if

      this%nsite = space%nsite
      this%orbital_lmax = radial_bases(1)%lmax
      this%response_lmax = space%response_lmax
      this%npoint = space%npoint
      this%mesh_a = space%a
      this%mesh_b = space%b
      this%selector = trim(selection)
      this%core_policy = projected_core_policy
      this%selected_l = .false.
      select case (trim(selection))
      case (projected_selection_d)
         this%selector_kind = projected_selector_d
         this%selected_l(2) = .true.
      case (projected_selection_spd)
         this%selector_kind = projected_selector_spd
         this%selected_l = .true.
      end select

      if (this%orbital_lmax /= 2 .or. this%response_lmax /= 2*this%orbital_lmax) then
         error stop 'DRESP-01: only the accepted lmax=2, response_lmax=4 spd contract is supported'
      end if
      allocate(this%radius(this%npoint), this%radial_weights(this%npoint))
      this%radius = space%radius
      this%radial_weights = space%radial_weights
      if (any(.not. ieee_is_finite(this%radial_weights)) .or. any(this%radial_weights < 0.0_rp)) then
         error stop 'DRESP-01: invalid response-space metric weights'
      end if
      this%angular_scalar_integral = 4.0_rp*response_angular_pi* &
         real(response_harmonic(0, 0, 0.0_rp, 0.0_rp), rp)
      if (.not. ieee_is_finite(this%angular_scalar_integral) .or. this%angular_scalar_integral <= 0.0_rp) then
         error stop 'DRESP-01: invalid Y00 full-sphere normalization'
      end if
      call validate_radial_bases(this, radial_bases)
   end subroutine projected_site_spin_initialize

   module subroutine projected_site_spin_functional(this, product, functionals)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(out) :: functionals(:, :)

      integer :: flat, site, response_l, response_m, product_mode

      call validate_product(this, product, 'DRESP-01 site integration')
      functionals = cmplx(0.0_rp, 0.0_rp, rp)
      do flat = 1, product%product_dimension
         call product%unflatten_index(flat, site, response_l, response_m, product_mode)
         if (response_l /= 0 .or. response_m /= 0) cycle
         ! The row integral is a^T z.  Store p=conjg(a), so p^H z is the
         ! Fortran DOT_PRODUCT(functional,z) contraction below.
         functionals(flat, site) = conjg(cmplx(this%angular_scalar_integral, 0.0_rp, rp)* &
            sum(sqrt(this%radial_weights)*product%blocks(site, response_l)%weighted_modes(:, product_mode)))
      end do
   end subroutine projected_site_spin_functional

   module subroutine projected_selected_coordinates(this, product, left_state, right_state, coordinates)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: coordinates(:)

      complex(rp), allocatable :: coefficients(:)
      integer :: site, response_l, response_m, k, first, last

      call validate_product(this, product, 'DRESP-01 selected coordinates')
      coordinates = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         do response_l = 0, product%response_lmax
            allocate(coefficients(product%blocks(site, response_l)%ncandidate))
            do response_m = -response_l, response_l
               call product%candidate_coefficients(site, response_l, response_m, left_state, right_state, coefficients)
               do k = 1, size(coefficients)
                  if (.not. this%selected_l(product%blocks(site, response_l)%candidates(k)%l) .or. &
                      .not. this%selected_l(product%blocks(site, response_l)%candidates(k)%lp)) then
                     coefficients(k) = cmplx(0.0_rp, 0.0_rp, rp)
                  end if
               end do
               first = product%flat_index(site, response_l, response_m, 1)
               last = product%flat_index(site, response_l, response_m, product%blocks(site, response_l)%rank)
               coordinates(first:last) = matmul(product%blocks(site, response_l)%forward_transform, coefficients)
            end do
            deallocate(coefficients)
         end do
      end do
   end subroutine projected_selected_coordinates

   module subroutine projected_transition_amplitudes(this, product, left_state, right_state, amplitudes)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      complex(rp), intent(out) :: amplitudes(:)

      complex(rp), allocatable :: coordinates(:), functionals(:, :)
      integer :: site

      call validate_product(this, product, 'DRESP-01 transition amplitudes')
      allocate(coordinates(product%product_dimension), functionals(product%product_dimension, this%nsite))
      call this%selected_transition_coordinates(product, left_state, right_state, coordinates)
      call this%site_integration_functional(product, functionals)
      amplitudes = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         amplitudes(site) = dot_product(functionals(:, site), coordinates)
      end do
   end subroutine projected_transition_amplitudes

   module subroutine projected_direct_operator_matrix(this, radial_bases, energy, operator_kind, matrix)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(in) :: energy
      integer, intent(in) :: operator_kind
      complex(rp), intent(out) :: matrix(:, :)

      integer :: norb, site, iorb, l, up, down, offset
      real(rp) :: overlap

      call validate_radial_bases(this, radial_bases)
      call validate_operator_kind(operator_kind, 'DRESP-01 direct operator')
      norb = (this%orbital_lmax + 1)**2
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         offset = (site - 1)*2*norb
         do iorb = 1, norb
            l = lmto_orbital_l(iorb)
            if (.not. this%selected_l(l)) cycle
            up = offset + iorb
            down = offset + norb + iorb
            select case (operator_kind)
            case (projected_operator_plus)
               overlap = endpoint_radial_overlap(radial_bases(site), l, 1, energy, l, 2, energy)
               matrix(up, down) = cmplx(overlap, 0.0_rp, rp)
            case (projected_operator_minus)
               overlap = endpoint_radial_overlap(radial_bases(site), l, 2, energy, l, 1, energy)
               matrix(down, up) = cmplx(overlap, 0.0_rp, rp)
            case (projected_operator_z)
               matrix(up, up) = cmplx(endpoint_radial_overlap(radial_bases(site), l, 1, energy, l, 1, energy), 0.0_rp, rp)
               matrix(down, down) = cmplx(-endpoint_radial_overlap(radial_bases(site), l, 2, energy, l, 2, energy), &
                                          0.0_rp, rp)
            end select
         end do
      end do
   end subroutine projected_direct_operator_matrix

   module subroutine projected_direct_transition_amplitudes(this, radial_bases, left_state, right_state, operator_kind, amplitudes)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      integer, intent(in) :: operator_kind
      complex(rp), intent(out) :: amplitudes(:)

      integer :: norb, site, iorb, l, offset
      complex(rp) :: left_coefficient, right_coefficient

      call validate_radial_bases(this, radial_bases)
      call validate_operator_kind(operator_kind, 'DRESP-01 direct transition')
      norb = (this%orbital_lmax + 1)**2
      call validate_endpoint_shapes(this%nsite, norb, left_state, right_state, 'DRESP-01 direct transition')
      amplitudes = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         offset = (site - 1)*2*norb
         do iorb = 1, norb
            l = lmto_orbital_l(iorb)
            if (.not. this%selected_l(l)) cycle
            select case (operator_kind)
            case (projected_operator_plus)
               left_coefficient = left_state%coefficients(offset + iorb)
               right_coefficient = right_state%coefficients(offset + norb + iorb)
               amplitudes(site) = amplitudes(site) + conjg(left_coefficient)*right_coefficient* &
                  endpoint_radial_overlap(radial_bases(site), l, 1, left_state%energy, l, 2, right_state%energy)
            case (projected_operator_minus)
               left_coefficient = left_state%coefficients(offset + norb + iorb)
               right_coefficient = right_state%coefficients(offset + iorb)
               amplitudes(site) = amplitudes(site) + conjg(left_coefficient)*right_coefficient* &
                  endpoint_radial_overlap(radial_bases(site), l, 2, left_state%energy, l, 1, right_state%energy)
            case (projected_operator_z)
               left_coefficient = left_state%coefficients(offset + iorb)
               right_coefficient = right_state%coefficients(offset + iorb)
               amplitudes(site) = amplitudes(site) + conjg(left_coefficient)*right_coefficient* &
                  endpoint_radial_overlap(radial_bases(site), l, 1, left_state%energy, l, 1, right_state%energy)
               left_coefficient = left_state%coefficients(offset + norb + iorb)
               right_coefficient = right_state%coefficients(offset + norb + iorb)
               amplitudes(site) = amplitudes(site) - conjg(left_coefficient)*right_coefficient* &
                  endpoint_radial_overlap(radial_bases(site), l, 2, left_state%energy, l, 2, right_state%energy)
            end select
         end do
      end do
   end subroutine projected_direct_transition_amplitudes

   !> Build the six second-order endpoint components of the selected, site-integrated
   !> DRESP-01 operator in coefficient space.  This is the GF-facing form of
   !> the same operator used by transition_amplitudes; it contains no energy
   !> denominator and no response accumulation.
   !>
   !> Component ordering is the authoritative `00,10,01,11,20,02` branch map.
   !> The returned tensor is indexed as
   !> `(coefficient-left, coefficient-right, component, site)`.
   module subroutine projected_site_component_vertex_tensor(this, product, vertices)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), allocatable, intent(out) :: vertices(:, :, :, :)

      complex(rp), allocatable :: product_vertices(:, :, :, :), functionals(:, :)
      integer :: nbasis, flat, site

      call validate_product(this, product, 'DRESP-01 site component vertices')
      call product%component_vertex_tensor(product_vertices, this%selected_l)
      allocate(functionals(product%product_dimension, this%nsite))
      call this%site_integration_functional(product, functionals)
      nbasis = size(product_vertices, 1)
      allocate(vertices(nbasis, nbasis, lmto_product_nbranch, this%nsite))
      vertices = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, this%nsite
         do flat = 1, product%product_dimension
            vertices(:, :, :, site) = vertices(:, :, :, site) + &
               conjg(functionals(flat, site))*product_vertices(:, :, :, flat)
         end do
      end do
      deallocate(product_vertices, functionals)
   end subroutine projected_site_component_vertex_tensor

   module subroutine projected_moment_from_density(this, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, &
                                            ground_states, moment)
      class(projected_site_spin_contract), intent(in) :: this
      real(rp), intent(in) :: eigenvalues(:, :), k_weights(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      type(radial_ground_state), intent(in) :: ground_states(:)
      real(rp), intent(out) :: moment(:)

      integer :: site, ik, ib, ir, spin, iorb, l, norb, nbands, nk, offset
      real(rp) :: wk, occupation, coefficient_weight, u_value
      real(rp), allocatable :: profile(:, :)

      call validate_moment_inputs(this, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, &
                                  ground_states, moment)
      norb = (this%orbital_lmax + 1)**2
      nbands = size(eigenvalues, 1)
      nk = size(eigenvalues, 2)
      moment = 0.0_rp
      do site = 1, this%nsite
         allocate(profile(this%npoint, 2))
         profile = 0.0_rp
         do ik = 1, nk
            wk = k_weights(ik)/sum(k_weights)
            do ib = 1, nbands
               occupation = fermi_dirac_occupation(eigenvalues(ib, ik), fermi_level, temperature)
               if (occupation == 0.0_rp) cycle
               offset = (site - 1)*2*norb
               do spin = 1, 2
                  do iorb = 1, norb
                     l = lmto_orbital_l(iorb)
                     if (.not. this%selected_l(l)) cycle
                     coefficient_weight = real(eigenvectors(offset + (spin - 1)*norb + iorb, ib, ik)* &
                        conjg(eigenvectors(offset + (spin - 1)*norb + iorb, ib, ik)), rp)
                     do ir = 1, this%npoint
                        u_value = ground_states(site)%pauli_large(ir, l + 1, spin) + &
                           (eigenvalues(ib, ik) - ground_states(site)%pauli_enu(l + 1, spin))* &
                           ground_states(site)%pauli_large_dot(ir, l + 1, spin)
                        profile(ir, spin) = profile(ir, spin) + wk*occupation*coefficient_weight*u_value**2
                     end do
                  end do
               end do
            end do
         end do
         moment(site) = ground_states(site)%log_mesh_integral(profile(:, 1) - profile(:, 2))
         deallocate(profile)
      end do
   end subroutine projected_moment_from_density

   module subroutine projected_moment_from_operator(this, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, &
                                             ground_states, moment)
      class(projected_site_spin_contract), intent(in) :: this
      real(rp), intent(in) :: eigenvalues(:, :), k_weights(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      type(radial_ground_state), intent(in) :: ground_states(:)
      real(rp), intent(out) :: moment(:)

      complex(rp), allocatable :: local_state(:), local_operator(:, :), action(:)
      integer :: norb, nbands, nk, site, ik, ib, iorb, l, spin, offset
      real(rp) :: wk, occupation, expectation

      call validate_moment_inputs(this, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, &
                                  ground_states, moment)
      norb = (this%orbital_lmax + 1)**2
      nbands = size(eigenvalues, 1)
      nk = size(eigenvalues, 2)
      moment = 0.0_rp
      allocate(local_state(2*norb), local_operator(2*norb, 2*norb), action(2*norb))
      do ik = 1, nk
         wk = k_weights(ik)/sum(k_weights)
         do ib = 1, nbands
            occupation = fermi_dirac_occupation(eigenvalues(ib, ik), fermi_level, temperature)
            if (occupation == 0.0_rp) cycle
            do site = 1, this%nsite
               offset = (site - 1)*2*norb
               local_state = eigenvectors(offset + 1:offset + 2*norb, ib, ik)
               local_operator = cmplx(0.0_rp, 0.0_rp, rp)
               do spin = 1, 2
                  do iorb = 1, norb
                     l = lmto_orbital_l(iorb)
                     if (.not. this%selected_l(l)) cycle
                     local_operator((spin - 1)*norb + iorb, (spin - 1)*norb + iorb) = &
                        cmplx(merge(1.0_rp, -1.0_rp, spin == 1)* &
                              ground_state_radial_norm(ground_states(site), l, spin, eigenvalues(ib, ik)), 0.0_rp, rp)
                  end do
               end do
               action = matmul(local_operator, local_state)
               expectation = real(dot_product(local_state, action), rp)
               moment(site) = moment(site) + wk*occupation*expectation
            end do
         end do
      end do
      deallocate(local_state, local_operator, action)
   end subroutine projected_moment_from_operator

   module subroutine projected_core_spin_number(this, ground_states, core_spin)
      class(projected_site_spin_contract), intent(in) :: this
      type(radial_ground_state), intent(in) :: ground_states(:)
      real(rp), intent(out) :: core_spin(:)
      integer :: site

      core_spin = 0.0_rp
      do site = 1, this%nsite
         if (.not. ground_states(site)%core_density_valid .or. &
             .not. allocated(ground_states(site)%core_pauli_weighted_up) .or. &
             .not. allocated(ground_states(site)%core_pauli_weighted_down)) cycle
         core_spin(site) = ground_states(site)%log_mesh_integral( &
            ground_states(site)%core_pauli_weighted_up - ground_states(site)%core_pauli_weighted_down)
      end do
   end subroutine projected_core_spin_number

   subroutine validate_product(this, product, caller)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_product_response_basis), intent(in) :: product
      character(len=*), intent(in) :: caller

      if (.not. allocated(product%blocks) .or. product%nsite /= this%nsite .or. &
          product%orbital_lmax /= this%orbital_lmax .or. product%response_lmax /= this%response_lmax .or. &
          product%npoint /= this%npoint .or. product%product_dimension < 1) then
         error stop trim(caller)//': product representation does not match the DRESP-01 contract'
      end if
      if (product%circular_channel /= lmto_product_channel_plus .and. &
          product%circular_channel /= lmto_product_channel_minus) then
         error stop trim(caller)//': invalid circular product channel'
      end if
   end subroutine validate_product

   subroutine validate_radial_bases(this, radial_bases)
      class(projected_site_spin_contract), intent(in) :: this
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      integer :: site

      if (size(radial_bases) /= this%nsite) then
         error stop 'DRESP-01 radial validation: site shape mismatch'
      end if
      do site = 1, this%nsite
         call radial_bases(site)%require_supported('ham_only', 'second', .true., .true., .false., .false.)
         if (radial_bases(site)%lmax /= this%orbital_lmax .or. radial_bases(site)%nspin /= 2 .or. &
             radial_bases(site)%npoint /= this%npoint .or. .not. allocated(radial_bases(site)%rofi) .or. &
             .not. allocated(radial_bases(site)%phi_large) .or. .not. allocated(radial_bases(site)%phidot_large) .or. &
             .not. allocated(radial_bases(site)%enu_work) .or. .not. allocated(radial_bases(site)%channel_present)) then
            error stop 'DRESP-01 radial validation: incomplete accepted radial basis'
         end if
         if (.not. all(radial_bases(site)%channel_present)) then
            error stop 'DRESP-01 radial validation: incomplete accepted radial basis'
         end if
         if (any(.not. ieee_is_finite(radial_bases(site)%rofi)) .or. &
             any(.not. ieee_is_finite(radial_bases(site)%phi_large)) .or. &
             any(.not. ieee_is_finite(radial_bases(site)%phidot_large)) .or. &
             any(.not. ieee_is_finite(radial_bases(site)%enu_work))) then
            error stop 'DRESP-01 radial validation: non-finite radial data is rejected'
         end if
         if (abs(radial_bases(site)%mesh_a - this%mesh_a) > mesh_tolerance*max(1.0_rp, abs(this%mesh_a)) .or. &
             abs(radial_bases(site)%mesh_b - this%mesh_b) > mesh_tolerance*max(1.0_rp, abs(this%mesh_b)) .or. &
             maxval(abs(radial_bases(site)%rofi - this%radius)) > mesh_tolerance* &
                max(1.0_rp, maxval(abs(this%radius)))) then
            error stop 'DRESP-01 radial validation: mesh provenance mismatch'
         end if
      end do
   end subroutine validate_radial_bases

   subroutine validate_operator_kind(operator_kind, caller)
      integer, intent(in) :: operator_kind
      character(len=*), intent(in) :: caller

      if (operator_kind < projected_operator_plus .or. operator_kind > projected_operator_z) then
         error stop trim(caller)//': operator kind must be plus, minus, or z'
      end if
   end subroutine validate_operator_kind

   subroutine validate_endpoint_shapes(nsite, norb, left_state, right_state, caller)
      integer, intent(in) :: nsite, norb
      type(pauli_endpoint_state), intent(in) :: left_state, right_state
      character(len=*), intent(in) :: caller

      if (.not. allocated(left_state%coefficients) .or. .not. allocated(right_state%coefficients) .or. &
          size(left_state%coefficients) /= 2*norb*nsite .or. size(right_state%coefficients) /= 2*norb*nsite) then
         error stop trim(caller)//': endpoint coefficient shape mismatch'
      end if
   end subroutine validate_endpoint_shapes

   subroutine validate_moment_inputs(this, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, &
                                     ground_states, moment)
      class(projected_site_spin_contract), intent(in) :: this
      real(rp), intent(in) :: eigenvalues(:, :), k_weights(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      type(radial_ground_state), intent(in) :: ground_states(:)
      real(rp), intent(out) :: moment(:)
      integer :: site
      integer :: expected_basis

      expected_basis = 2*((this%orbital_lmax + 1)**2)*this%nsite
      if (size(moment) /= this%nsite .or. size(ground_states) /= this%nsite .or. &
          size(eigenvalues, 1) < 1 .or. size(eigenvalues, 2) < 1 .or. size(k_weights) /= size(eigenvalues, 2) .or. &
          any(shape(eigenvectors) /= [expected_basis, size(eigenvalues, 1), size(eigenvalues, 2)]) .or. &
          sum(k_weights) <= tiny(1.0_rp) .or. any(k_weights < 0.0_rp)) then
         error stop 'DRESP-01 moment: accepted eigensystem or k-weight shape is invalid'
      end if
      do site = 1, this%nsite
         if (.not. ground_states(site)%valid .or. .not. ground_states(site)%accepted .or. &
             .not. ground_states(site)%pauli_basis_valid .or. ground_states(site)%pauli_lmax < this%orbital_lmax .or. &
             .not. allocated(ground_states(site)%r) .or. .not. allocated(ground_states(site)%pauli_large) .or. &
             .not. allocated(ground_states(site)%pauli_large_dot) .or. .not. allocated(ground_states(site)%pauli_enu)) then
            error stop 'DRESP-01 moment: accepted Pauli radial snapshot is required'
         end if
         if (size(ground_states(site)%r) /= this%npoint .or. &
             any(shape(ground_states(site)%pauli_large) /= [this%npoint, this%orbital_lmax + 1, 2]) .or. &
             any(shape(ground_states(site)%pauli_large_dot) /= [this%npoint, this%orbital_lmax + 1, 2]) .or. &
             any(shape(ground_states(site)%pauli_enu) /= [this%orbital_lmax + 1, 2])) then
            error stop 'DRESP-01 moment: accepted Pauli radial snapshot dimensions are invalid'
         end if
         if (abs(ground_states(site)%a - this%mesh_a) > mesh_tolerance*max(1.0_rp, abs(this%mesh_a)) .or. &
             abs(ground_states(site)%b - this%mesh_b) > mesh_tolerance*max(1.0_rp, abs(this%mesh_b)) .or. &
             maxval(abs(ground_states(site)%r - this%radius)) > mesh_tolerance*max(1.0_rp, maxval(abs(this%radius)))) then
            error stop 'DRESP-01 moment: radial snapshot mesh provenance mismatch'
         end if
      end do
   end subroutine validate_moment_inputs

   pure real(rp) function endpoint_radial_value(radial, ir, l, spin, energy) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin
      real(rp), intent(in) :: energy

      value = radial%phi_large(ir, l + 1, spin) + (energy - radial%enu_work(l + 1, spin))* &
         radial%phidot_large(ir, l + 1, spin)
   end function endpoint_radial_value

   real(rp) function endpoint_radial_overlap(radial, l_left, spin_left, energy_left, l_right, spin_right, energy_right) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: l_left, spin_left, l_right, spin_right
      real(rp), intent(in) :: energy_left, energy_right
      integer :: ir

      value = 0.0_rp
      do ir = 1, radial%npoint
         value = value + radial_simpson_weight(ir, radial%npoint)*radial%mesh_a* &
            (radial%rofi(ir) + radial%mesh_b)*endpoint_radial_value(radial, ir, l_left, spin_left, energy_left)* &
            endpoint_radial_value(radial, ir, l_right, spin_right, energy_right)
      end do
   end function endpoint_radial_overlap

   real(rp) function ground_state_radial_norm(ground_state, l, spin, energy) result(value)
      type(radial_ground_state), intent(in) :: ground_state
      integer, intent(in) :: l, spin
      real(rp), intent(in) :: energy
      integer :: ir
      real(rp) :: u_value

      value = 0.0_rp
      do ir = 1, size(ground_state%r)
         u_value = ground_state%pauli_large(ir, l + 1, spin) + (energy - ground_state%pauli_enu(l + 1, spin))* &
            ground_state%pauli_large_dot(ir, l + 1, spin)
         value = value + radial_simpson_weight(ir, size(ground_state%r))*ground_state%a* &
            (ground_state%r(ir) + ground_state%b)*u_value**2
      end do
   end function ground_state_radial_norm

   pure real(rp) function fermi_dirac_occupation(eigenvalue, fermi_level, temperature) result(occupation)
      real(rp), intent(in) :: eigenvalue, fermi_level, temperature
      real(rp) :: argument

      argument = (eigenvalue - fermi_level)/max(temperature*kB_Ry_per_K, occupation_kT_floor)
      if (argument >= 50.0_rp) then
         occupation = 0.0_rp
      else if (argument <= -50.0_rp) then
         occupation = 1.0_rp
      else
         occupation = 1.0_rp/(exp(argument) + 1.0_rp)
      end if
   end function fermi_dirac_occupation


end submodule linear_response_basis
