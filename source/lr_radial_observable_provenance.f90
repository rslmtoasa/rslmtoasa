!------------------------------------------------------------------------------
! DRESP-09W radial observable provenance services.
!
! The routines in this module are deliberately representation-explicit.  They
! provide the Pauli large-component polynomial, the complete scalar-
! relativistic polynomial, independent occupied-state oracles, and the Pauli
! L=0 response contraction.  No SCF, field mapping, kernel, or Dyson state is
! changed here.
!------------------------------------------------------------------------------
module lr_radial_observable_provenance_mod

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   use radial_ground_state_mod, only: radial_ground_state
   use reciprocal_mod, only: reciprocal
   use response_angular_basis_mod, only: response_gaunt
   use response_basis_mapping_mod, only: response_super_index, response_unflatten_superindex
   use lr_response_space_mod, only: response_space_layout
   implicit none
   private

   real(rp), parameter, public :: dresp09w_four_pi = 4.0_rp*acos(-1.0_rp)

   public :: channel_moments_from_matrix
   public :: pauli_density_from_moments
   public :: sr_density_from_moments
   public :: pauli_density_from_occupied_states
   public :: sr_density_from_occupied_states
   public :: pauli_l0_density_from_endpoint_branches
   public :: weighted_from_density
   public :: physical_from_weighted
   public :: volume_integral
   public :: volume_l2_norm
   public :: radial_relative_metrics

contains

   !> Extract equilibrium diagonal site/l/spin blocks from production M0/M1/M2.
   !> The result is real because the accepted state is collinear; the caller
   !> supplies the k-weight normalization explicitly through the matrix input.
   subroutine channel_moments_from_matrix(moment_matrices, nsite, moments)
      complex(rp), intent(in) :: moment_matrices(:, :, :)
      integer, intent(in) :: nsite
      real(rp), intent(out) :: moments(:, :, :, :) ! (moment,site,l+1,spin)
      integer :: n, norb, site, spin, l, iorb, first, last

      n = size(moment_matrices, 1)
      if (size(moment_matrices, 2) /= n .or. size(moment_matrices, 3) /= 3 .or. &
          nsite < 1 .or. mod(n, 2*nsite) /= 0) error stop 'DRESP-09W channel moments: shape mismatch'
      norb = n/(2*nsite)
      if (size(moments, 1) /= 3 .or. size(moments, 2) /= nsite .or. &
          size(moments, 3) < 1 .or. size(moments, 4) /= 2) then
         error stop 'DRESP-09W channel moments: output shape mismatch'
      end if
      moments = 0.0_rp
      do site = 1, nsite
         do spin = 1, 2
            first = (site - 1)*2*norb + (spin - 1)*norb + 1
            last = first + norb - 1
            do iorb = first, last
               l = lmto_orbital_l(iorb - first + 1)
               if (l < 0 .or. l + 1 > size(moments, 3)) cycle
               moments(:, site, l + 1, spin) = moments(:, site, l + 1, spin) + &
                  real([(moment_matrices(iorb, iorb, 1)), (moment_matrices(iorb, iorb, 2)), &
                        (moment_matrices(iorb, iorb, 3))], rp)
            end do
         end do
      end do
   end subroutine channel_moments_from_matrix

   !> Production Pauli moment polynomial.  This is the large-component-only
   !> observable, but it retains the legacy phi*phiddot endpoint term.
   subroutine pauli_density_from_moments(radial_bases, channel_moments, density, density_l)
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(in) :: channel_moments(:, :, :, :)
      real(rp), intent(out) :: density(:, :)
      real(rp), intent(out), optional :: density_l(:, :, :)
      integer :: site, ir, l, spin
      real(rp) :: q0, q1, q2, value, sign_spin

      call validate_channel_shapes(radial_bases, channel_moments, density)
      density = 0.0_rp
      if (present(density_l)) then
         if (any(shape(density_l) /= [size(density, 1), size(density, 2), size(channel_moments, 3)])) then
            error stop 'DRESP-09W Pauli moments: l-resolved output shape mismatch'
         end if
         density_l = 0.0_rp
      end if
      do site = 1, size(radial_bases)
         do l = 0, radial_bases(site)%lmax
            do spin = 1, 2
               q0 = channel_moments(1, site, l + 1, spin)
               q1 = channel_moments(2, site, l + 1, spin) - radial_bases(site)%enu_work(l + 1, spin)*q0
               q2 = channel_moments(3, site, l + 1, spin) - &
                    2.0_rp*radial_bases(site)%enu_work(l + 1, spin)*channel_moments(2, site, l + 1, spin) + &
                    radial_bases(site)%enu_work(l + 1, spin)**2*q0
               sign_spin = merge(1.0_rp, -1.0_rp, spin == 1)
               do ir = 1, radial_bases(site)%npoint
                  value = q0*radial_bases(site)%phi_large(ir, l + 1, spin)**2 + &
                     2.0_rp*q1*radial_bases(site)%phi_large(ir, l + 1, spin)* &
                     radial_bases(site)%phidot_large(ir, l + 1, spin) + &
                     q2*(radial_bases(site)%phidot_large(ir, l + 1, spin)**2 + &
                         radial_bases(site)%phi_large(ir, l + 1, spin)*radial_bases(site)%phiddot_large(ir, l + 1, spin))
                  density(site, ir) = density(site, ir) + sign_spin*value
                  if (present(density_l)) density_l(site, ir, l + 1) = density_l(site, ir, l + 1) + sign_spin*value
               end do
            end do
         end do
      end do
   end subroutine pauli_density_from_moments

   !> Full scalar-relativistic production polynomial, delegated to the live
   !> lmto_density_from_moments implementation so no second convention exists.
   subroutine sr_density_from_moments(radial_bases, channel_moments, density, density_l)
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(in) :: channel_moments(:, :, :, :)
      real(rp), intent(out) :: density(:, :)
      real(rp), intent(out), optional :: density_l(:, :, :)
      real(rp), allocatable :: local_moments(:, :, :), local_density(:, :, :)
      integer :: site, l, spin

      call validate_channel_shapes(radial_bases, channel_moments, density)
      allocate(local_moments(3, radial_bases(1)%lmax + 1, 2), local_density(radial_bases(1)%npoint, &
         radial_bases(1)%lmax + 1, 2))
      density = 0.0_rp
      if (present(density_l)) then
         if (any(shape(density_l) /= [size(density, 1), size(density, 2), size(channel_moments, 3)])) then
            error stop 'DRESP-09W SR moments: l-resolved output shape mismatch'
         end if
         density_l = 0.0_rp
      end if
      do site = 1, size(radial_bases)
         do l = 1, radial_bases(site)%lmax + 1
            do spin = 1, 2
               local_moments(:, l, spin) = channel_moments(:, site, l, spin)
            end do
         end do
         call radial_bases(site)%density_from_moments(local_moments, local_density)
         do spin = 1, 2
            density(site, :) = density(site, :) + merge(1.0_rp, -1.0_rp, spin == 1)*sum(local_density(:, :, spin), dim=2)
            if (present(density_l)) then
               do l = 1, radial_bases(site)%lmax + 1
                  density_l(site, :, l) = density_l(site, :, l) + &
                     merge(1.0_rp, -1.0_rp, spin == 1)*local_density(:, l, spin)
               end do
            end if
         end do
      end do
      deallocate(local_moments, local_density)
   end subroutine sr_density_from_moments

   !> Independent Pauli occupied-state oracle: explicit phi+DeltaE*phidot.
   subroutine pauli_density_from_occupied_states(reciprocal_obj, radial_bases, density, density_l)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(out) :: density(:, :)
      real(rp), intent(out), optional :: density_l(:, :, :)
      integer :: site, ir, l, spin, iorb, ib, ik, ik_global, norb, nsite
      real(rp) :: wk, occ, energy, delta, amplitude, value
      complex(rp) :: coefficient

      call validate_reciprocal_radial(reciprocal_obj, radial_bases, density)
      nsite = size(radial_bases); norb = (radial_bases(1)%lmax + 1)**2
      density = 0.0_rp
      if (present(density_l)) density_l = 0.0_rp
      do site = 1, nsite
         do ik = 1, size(reciprocal_obj%eigenvalues, 2)
            ik_global = ik
            if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
            wk = reciprocal_obj%k_weights(ik_global)
            do ib = 1, size(reciprocal_obj%eigenvalues, 1)
               energy = reciprocal_obj%eigenvalues(ib, ik)
               occ = fermi_occupation_local(energy, reciprocal_obj%fermi_level, reciprocal_obj%temperature)
               if (occ <= 1.0e-14_rp) cycle
               do spin = 1, 2
                  do iorb = 1, norb
                     l = lmto_orbital_l(iorb)
                     coefficient = reciprocal_obj%eigenvectors((site - 1)*2*norb + (spin - 1)*norb + iorb, ib, ik)
                     amplitude = real(coefficient*conjg(coefficient), rp)
                     delta = energy - radial_bases(site)%enu_work(l + 1, spin)
                     do ir = 1, radial_bases(site)%npoint
                        value = amplitude*(radial_bases(site)%phi_large(ir, l + 1, spin) + &
                           delta*radial_bases(site)%phidot_large(ir, l + 1, spin))**2
                        density(site, ir) = density(site, ir) + merge(1.0_rp, -1.0_rp, spin == 1)*wk*occ*value
                        if (present(density_l)) density_l(site, ir, l + 1) = density_l(site, ir, l + 1) + &
                           merge(1.0_rp, -1.0_rp, spin == 1)*wk*occ*value
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine pauli_density_from_occupied_states

   !> Independent full-SR occupied-state oracle using lmto_state_density.
   subroutine sr_density_from_occupied_states(reciprocal_obj, radial_bases, density, density_l)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(out) :: density(:, :)
      real(rp), intent(out), optional :: density_l(:, :, :)
      complex(rp), allocatable :: coefficients(:, :)
      real(rp), allocatable :: local_density(:, :, :)
      integer :: site, ik, ik_global, ib, spin, iorb, norb, nsite
      real(rp) :: wk, occ, energy

      call validate_reciprocal_radial(reciprocal_obj, radial_bases, density)
      nsite = size(radial_bases); norb = (radial_bases(1)%lmax + 1)**2
      allocate(coefficients(norb, 2), local_density(radial_bases(1)%npoint, radial_bases(1)%lmax + 1, 2))
      density = 0.0_rp
      if (present(density_l)) density_l = 0.0_rp
      do site = 1, nsite
         do ik = 1, size(reciprocal_obj%eigenvalues, 2)
            ik_global = ik
            if (allocated(reciprocal_obj%k_l2g_map)) ik_global = reciprocal_obj%k_l2g_map(ik)
            wk = reciprocal_obj%k_weights(ik_global)
            do ib = 1, size(reciprocal_obj%eigenvalues, 1)
               energy = reciprocal_obj%eigenvalues(ib, ik)
               occ = fermi_occupation_local(energy, reciprocal_obj%fermi_level, reciprocal_obj%temperature)
               if (occ <= 1.0e-14_rp) cycle
               do spin = 1, 2
                  do iorb = 1, norb
                     coefficients(iorb, spin) = reciprocal_obj%eigenvectors((site - 1)*2*norb + &
                        (spin - 1)*norb + iorb, ib, ik)
                  end do
               end do
               call radial_bases(site)%state_density(energy, coefficients, local_density)
               do spin = 1, 2
                  density(site, :) = density(site, :) + merge(1.0_rp, -1.0_rp, spin == 1)*wk*occ* &
                     sum(local_density(:, :, spin), dim=2)
                  if (present(density_l)) then
                     do iorb = 1, radial_bases(site)%lmax + 1
                        density_l(site, :, iorb) = density_l(site, :, iorb) + &
                           merge(1.0_rp, -1.0_rp, spin == 1)*wk*occ*local_density(:, iorb, spin)
                     end do
                  end if
               end do
            end do
         end do
      end do
      deallocate(coefficients, local_density)
   end subroutine sr_density_from_occupied_states

   !> Pauli large-component analogue of sr_l0_density_from_endpoint_branches.
   subroutine pauli_l0_density_from_endpoint_branches(space, radial_bases, endpoint_branches, circular_channel, density, &
                                                      l_first, l_last)
      type(response_space_layout), intent(in) :: space
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      complex(rp), intent(in) :: endpoint_branches(:, :, :)
      integer, intent(in) :: circular_channel
      complex(rp), intent(out) :: density(:, :)
      integer, intent(in), optional :: l_first, l_last
      integer :: nsite, norb, n, site, ir, iorb, l, m, p, q, branch, row, col, first_l, last_l
      integer :: spin_left, spin_right
      real(rp) :: gaunt, pair

      nsite = space%nsite; norb = (radial_bases(1)%lmax + 1)**2; n = 2*norb*nsite
      if (size(radial_bases) /= nsite .or. any(shape(endpoint_branches) /= [n, n, 4]) .or. &
          any(shape(density) /= [nsite, space%npoint])) error stop 'DRESP-09W Pauli response: shape mismatch'
      first_l = 0; last_l = radial_bases(1)%lmax
      if (present(l_first)) first_l = max(0, l_first)
      if (present(l_last)) last_l = min(radial_bases(1)%lmax, l_last)
      if (circular_channel == 1) then
         spin_left = 1; spin_right = 2
      else if (circular_channel == 2) then
         spin_left = 2; spin_right = 1
      else
         error stop 'DRESP-09W Pauli response: invalid circular channel'
      end if
      density = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, nsite
         do ir = 2, space%npoint
            do iorb = 1, norb
               l = lmto_orbital_l(iorb)
               if (l < first_l .or. l > last_l) cycle
               m = iorb - l*l - l - 1
               gaunt = response_gaunt(l, m, l, m, 0, 0)
               row = (site - 1)*2*norb + (spin_left - 1)*norb + iorb
               col = (site - 1)*2*norb + (spin_right - 1)*norb + iorb
               do p = 0, 1
                  do q = 0, 1
                     branch = 1 + p + 2*q
                     pair = pauli_endpoint_pair(radial_bases(site), ir, l, spin_left, p, spin_right, q)
                     density(site, ir) = density(site, ir) + gaunt*endpoint_branches(row, col, branch)*pair
                  end do
               end do
            end do
         end do
      end do
   end subroutine pauli_l0_density_from_endpoint_branches

   subroutine weighted_from_density(r, density, weighted)
      real(rp), intent(in) :: r(:), density(:, :)
      real(rp), intent(out) :: weighted(:, :)
      integer :: site, ir
      if (any(shape(weighted) /= shape(density)) .or. size(weighted, 2) /= size(r)) then
         error stop 'DRESP-09W weighted density: shape mismatch'
      end if
      weighted = 0.0_rp
      do site = 1, size(density, 1)
         do ir = 2, size(r)
            weighted(site, ir) = dresp09w_four_pi*r(ir)**2*density(site, ir)
         end do
      end do
   end subroutine weighted_from_density

   subroutine physical_from_weighted(r, weighted, density)
      real(rp), intent(in) :: r(:), weighted(:, :)
      real(rp), intent(out) :: density(:, :)
      integer :: site, ir
      if (any(shape(density) /= shape(weighted)) .or. size(density, 2) /= size(r)) then
         error stop 'DRESP-09W physical density: shape mismatch'
      end if
      density = 0.0_rp
      do site = 1, size(weighted, 1)
         do ir = 2, size(r)
            density(site, ir) = weighted(site, ir)/(dresp09w_four_pi*r(ir)**2)
         end do
      end do
   end subroutine physical_from_weighted

   real(rp) function volume_integral(state, weighted) result(value)
      type(radial_ground_state), intent(in) :: state
      real(rp), intent(in) :: weighted(:)
      integer :: ir
      value = 0.0_rp
      do ir = 2, size(state%r)
         value = value + radial_jacobian(state, ir)*weighted(ir)
      end do
   end function volume_integral

   real(rp) function volume_l2_norm(state, weighted) result(value)
      type(radial_ground_state), intent(in) :: state
      real(rp), intent(in) :: weighted(:)
      integer :: ir
      real(rp) :: r
      value = 0.0_rp
      do ir = 2, size(state%r)
         r = state%r(ir)
         value = value + radial_jacobian(state, ir)*weighted(ir)**2/(dresp09w_four_pi*r*r)
      end do
      value = sqrt(value)
   end function volume_l2_norm

   !> Relative volume L2, absolute volume L2, maximum physical absolute and
   !> maximum relative residual.  All diagnostics exclude the zero-measure
   !> origin point.
   subroutine radial_relative_metrics(values, target, states, relative, l2, max_absolute, max_relative, integral)
      real(rp), intent(in) :: values(:, :), target(:, :)
      type(radial_ground_state), intent(in) :: states(:)
      real(rp), intent(out) :: relative, l2, max_absolute, max_relative, integral
      integer :: site, ir
      real(rp) :: jac, r, num, den, difference, physical_difference, physical_target

      num = 0.0_rp; den = 0.0_rp; max_absolute = 0.0_rp; max_relative = 0.0_rp; integral = 0.0_rp
      do site = 1, size(states)
         do ir = 2, size(states(site)%r)
            r = states(site)%r(ir); jac = radial_jacobian(states(site), ir)
            difference = values(site, ir) - target(site, ir)
            num = num + jac*difference**2/(dresp09w_four_pi*r*r)
            den = den + jac*target(site, ir)**2/(dresp09w_four_pi*r*r)
            physical_difference = abs(difference)/(dresp09w_four_pi*r*r)
            physical_target = abs(target(site, ir))/(dresp09w_four_pi*r*r)
            max_absolute = max(max_absolute, physical_difference)
            max_relative = max(max_relative, physical_difference/max(physical_target, tiny(1.0_rp)))
            integral = integral + jac*values(site, ir)
         end do
      end do
      l2 = sqrt(num)
      relative = l2/max(sqrt(den), tiny(1.0_rp))
   end subroutine radial_relative_metrics

   pure real(rp) function radial_jacobian(state, ir) result(value)
      type(radial_ground_state), intent(in) :: state
      integer, intent(in) :: ir
      value = (2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp)*state%a*(state%r(ir) + state%b)
   end function radial_jacobian

   subroutine validate_channel_shapes(radial_bases, moments, density)
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(in) :: moments(:, :, :, :), density(:, :)
      if (size(radial_bases) < 1 .or. size(moments, 1) /= 3 .or. size(moments, 2) /= size(radial_bases) .or. &
          size(moments, 3) /= radial_bases(1)%lmax + 1 .or. size(moments, 4) /= 2 .or. &
          any(shape(density) /= [size(radial_bases), radial_bases(1)%npoint])) then
         error stop 'DRESP-09W channel density: shape mismatch'
      end if
   end subroutine validate_channel_shapes

   subroutine validate_reciprocal_radial(reciprocal_obj, radial_bases, density)
      type(reciprocal), intent(in) :: reciprocal_obj
      type(lmto_radial_basis), intent(in) :: radial_bases(:)
      real(rp), intent(in) :: density(:, :)
      integer :: nsite, nmat, norb
      nsite = size(radial_bases); norb = (radial_bases(1)%lmax + 1)**2
      if (nsite < 1 .or. any(shape(density) /= [nsite, radial_bases(1)%npoint]) .or. &
          .not. allocated(reciprocal_obj%eigenvalues) .or. .not. allocated(reciprocal_obj%eigenvectors) .or. &
          .not. allocated(reciprocal_obj%k_weights)) error stop 'DRESP-09W occupied oracle: incomplete state'
      nmat = size(reciprocal_obj%eigenvectors, 1)
      if (nmat /= 2*norb*nsite .or. size(reciprocal_obj%eigenvectors, 2) /= size(reciprocal_obj%eigenvalues, 1) .or. &
          size(reciprocal_obj%eigenvectors, 3) /= size(reciprocal_obj%eigenvalues, 2)) then
         error stop 'DRESP-09W occupied oracle: eigensystem/radial shape mismatch'
      end if
   end subroutine validate_reciprocal_radial

   pure real(rp) function pauli_endpoint_component(radial, ir, l, spin, power) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin, power
      if (power == 0) then
         value = radial%phi_large(ir, l + 1, spin) - radial%enu_work(l + 1, spin)* &
            radial%phidot_large(ir, l + 1, spin)
      else
         value = radial%phidot_large(ir, l + 1, spin)
      end if
   end function pauli_endpoint_component

   pure real(rp) function pauli_endpoint_pair(radial, ir, left_l, left_spin, left_power, right_spin, right_power) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, left_l, left_spin, left_power, right_spin, right_power
      real(rp) :: r
      if (ir == 1 .or. radial%rofi(ir) <= tiny(1.0_rp)) then
         value = 0.0_rp
         return
      end if
      r = radial%rofi(ir)
      value = pauli_endpoint_component(radial, ir, left_l, left_spin, left_power)* &
         pauli_endpoint_component(radial, ir, left_l, right_spin, right_power)/r**2
   end function pauli_endpoint_pair

   pure real(rp) function fermi_occupation_local(energy, fermi_level, temperature) result(value)
      real(rp), intent(in) :: energy, fermi_level, temperature
      real(rp) :: argument, kt
      kt = max(temperature*6.3336814e-6_rp, 1.0e-12_rp)
      argument = (energy - fermi_level)/kt
      if (argument >= 50.0_rp) then
         value = 0.0_rp
      else if (argument <= -50.0_rp) then
         value = 1.0_rp
      else
         value = 1.0_rp/(exp(argument) + 1.0_rp)
      end if
   end function fermi_occupation_local

end module lr_radial_observable_provenance_mod
