!------------------------------------------------------------------------------
! Independent direct-integration check for the literature-faithful Juelich-d
! projector, moment, spin-flip transition, angular factor, and eta bubble.
! Expected values below do not call any production projector routine.
!------------------------------------------------------------------------------
module test_juelich_d_reference_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, scalar_relativistic_c
   use linear_response_mod, only: lr_electronic_state, juelich_d_projector, juelich_d_state_projection, &
      evaluate_juelich_d_lehmann_chi0, extrapolate_juelich_static_chi0, lr_fermi_dirac_occupation, response_angular_pi
   implicit none
   private
   public :: run_test_juelich_d_reference

   integer, parameter :: nr = 9, lmax = 2, norb = 9, nbands = 18, nk = 1
   real(rp), parameter :: mesh_a = 0.18_rp, mesh_b = 0.12_rp, nuclear_z = 26.0_rp
   real(rp), parameter :: ef = 0.035_rp, temperature = 300.0_rp, tol = 3.0e-12_rp

contains

   subroutine run_test_juelich_d_reference()
      type(lmto_radial_basis) :: radial(1)
      type(juelich_d_projector) :: projector
      type(juelich_d_state_projection) :: projection
      type(lr_electronic_state) :: state
      real(rp) :: radius(nr), eigenvalues(nbands, nk), k_points(3, nk), k_weights(nk), occupations(nbands, nk)
      real(rp) :: expected_norm(2), expected_moment(1), moment(1), frequencies(2), eta, wrong_ef_moment, wrong_norm_moment
      complex(rp) :: eigenvectors(nbasis_value(), nbands, nk), chi(1, 1, 2), expected_chi(1, 1, 2)
      complex(rp) :: expected_projected(5, 2, nbands), endpoint_up, endpoint_down, transition
      real(rp) :: eta_fit(4), estimate, fit_error, eta_imaginary
      complex(rp) :: chi_fit(1, 1, 4), chi_static_fit(1, 1)
      real(rp) :: projector_energy_delta, frozen_transition, energy_dependent_transition, no_angular, double_spin
      integer :: ir, l, spin, ib, jb, im, iorb, ispin
      logical :: failed

      failed = .false.
      projector_energy_delta = 0.0_rp
      call basis_init(lmax)
      call build_mesh(radius)
      call build_radial(radial(1), radius)
      do ib = 1, nbands
         eigenvalues(ib, 1) = -0.24_rp + 0.026_rp*real(ib - 1, rp)
         occupations(ib, 1) = lr_fermi_dirac_occupation(eigenvalues(ib, 1), ef, temperature)
      end do
      eigenvectors = cmplx(0.0_rp, 0.0_rp, rp)
      do ib = 1, nbands
         eigenvectors(ib, ib, 1) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      k_points = 0.0_rp
      k_weights = 1.0_rp
      call state%initialize(eigenvalues, eigenvectors, k_points, k_weights, occupations, ef, temperature)

      call independent_projector_norms(radial(1), ef, expected_norm)
      call projector%initialize(radial, ef)
      if (maxval(abs(projector%normalization(1, :) - expected_norm)) > tol*maxval(expected_norm)) failed = .true.
      do ispin = 1, 2
         call independent_norm(radial(1), ef + 0.17_rp, ispin, wrong_ef_moment)
         projector_energy_delta = max(projector_energy_delta, abs(projector%normalization(1, ispin) - wrong_ef_moment))
      end do
      if (projector_energy_delta <= 1.0e-7_rp) failed = .true.

      call projector%project_state(radial, state, projection)
      call projector%moment_from_state(radial, state, moment)
      call independent_band_projection(radial(1), eigenvalues(:, 1), eigenvectors(:, :, 1), ef, expected_projected)
      expected_moment = 0.0_rp
      wrong_norm_moment = 0.0_rp
      do ib = 1, nbands
         do im = 1, 5
            expected_moment(1) = expected_moment(1) + occupations(ib, 1)*( &
               abs(expected_projected(im, 1, ib))**2 - abs(expected_projected(im, 2, ib))**2)
         end do
      end do
      if (abs(moment(1) - expected_moment(1)) > tol*max(1.0_rp, abs(expected_moment(1)))) failed = .true.
      do ib = 1, nbands
         do im = 1, 5
            wrong_norm_moment = wrong_norm_moment + occupations(ib, 1)*( &
               abs(expected_projected(im, 1, ib))**2*expected_norm(1) - &
               abs(expected_projected(im, 2, ib))**2*expected_norm(2))
         end do
      end do
      if (abs(moment(1) - wrong_norm_moment) <= 1.0e-5_rp) failed = .true.

      ! Compare one occupied-up/empty-down transition independently.  The
      ! expected d amplitudes are direct radial integrals against R_d(EF).
      ib = 6
      jb = norb + 6
      do im = 1, 5
         if (abs(projection%amplitudes(im, 1, 1, ib, 1) - expected_projected(im, 1, ib)) > tol .or. &
             abs(projection%amplitudes(im, 1, 2, jb, 1) - expected_projected(im, 2, jb)) > tol) failed = .true.
      end do
      endpoint_up = expected_projected(2, 1, ib)
      endpoint_down = expected_projected(2, 2, jb)
      frozen_transition = conjg(endpoint_up)*endpoint_down
      energy_dependent_transition = independent_energy_dependent_overlap(radial(1), 1, eigenvalues(ib, 1), &
         2, eigenvalues(jb, 1))
      if (abs(frozen_transition - energy_dependent_transition) <= 1.0e-7_rp) failed = .true.

      frequencies = [0.0_rp, 0.11_rp]
      eta = 0.027_rp
      call evaluate_juelich_d_lehmann_chi0(projector, radial, state, state, [0.0_rp, 0.0_rp, 0.0_rp], &
         frequencies, eta, chi)
      call independent_juelich_chi(radial(1), state, ef, frequencies, eta, expected_chi)
      if (maxval(abs(chi - expected_chi)) > tol*max(1.0_rp, maxval(abs(expected_chi)))) failed = .true.
      no_angular = maxval(abs(chi - expected_chi/(4.0_rp*response_angular_pi)))
      double_spin = maxval(abs(chi - 2.0_rp*expected_chi))
      if (no_angular <= 1.0e-6_rp .or. double_spin <= 1.0e-6_rp) failed = .true.

      eta_fit = [0.08_rp, 0.04_rp, 0.02_rp, 0.01_rp]
      chi_fit(1, 1, :) = cmplx(-2.3_rp + 0.7_rp*eta_fit**2, 1.4_rp*eta_fit, rp)
      call extrapolate_juelich_static_chi0(eta_fit, chi_fit, chi_static_fit, estimate, fit_error, eta_imaginary)
      if (abs(real(chi_static_fit(1, 1), rp) + 2.3_rp) > tol .or. abs(aimag(chi_static_fit(1, 1))) > tol .or. &
          estimate > tol .or. fit_error > tol .or. eta_imaginary <= 0.0_rp) failed = .true.

      write(*, '(a,2(es14.6,1x))') '  independent EF radial norms up/down = ', expected_norm
      write(*, '(a,es14.6,a,es14.6)') '  independent projected M_d = ', expected_moment(1), &
         ' production error = ', abs(moment(1) - expected_moment(1))
      write(*, '(a,2(es14.6,1x))') '  chi0(0,0; omega=0) Re/Im = ', real(chi(1, 1, 1), rp), aimag(chi(1, 1, 1))
      if (.not. all(ieee_is_finite(real(chi, rp))) .or. .not. all(ieee_is_finite(aimag(chi)))) failed = .true.
      if (failed) error stop 'UnitLrJuelichDReference: FAIL'
      write(*, '(a)') 'UnitLrJuelichDReference: PASS (independent EF norm, moment, transition, 4pi and unit spin factor)'
   end subroutine run_test_juelich_d_reference

   pure integer function nbasis_value() result(n)
      n = 2*norb
   end function nbasis_value

   subroutine build_mesh(r)
      real(rp), intent(out) :: r(:)
      integer :: ir
      r(1) = 0.0_rp
      do ir = 2, size(r)
         r(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine build_radial(radial, r)
      type(lmto_radial_basis), intent(out) :: radial
      real(rp), intent(in) :: r(:)
      real(rp) :: tmc, alpha, small_scale
      integer :: ir, l, spin

      call radial%initialize(size(r), lmax, 2)
      radial%rofi = r
      radial%mesh_a = mesh_a
      radial%mesh_b = mesh_b
      radial%nuclear_z = nuclear_z
      radial%channel_present = .true.
      do spin = 1, 2
         do ir = 1, size(r)
            radial%potential(ir, spin) = -0.21_rp + 0.08_rp*real(spin - 1, rp) + 0.013_rp*r(ir)
         end do
         do l = 0, lmax
            radial%enu_work(l + 1, spin) = -0.20_rp + 0.035_rp*real(l, rp) + &
               0.021_rp*real(spin - 1, rp)
            radial%enu_radial(l + 1, spin) = radial%enu_work(l + 1, spin) - 0.01_rp
            alpha = 0.38_rp + 0.07_rp*real(l, rp) + 0.04_rp*real(spin - 1, rp)
            small_scale = 0.07_rp + 0.015_rp*real(l + spin, rp)
            do ir = 1, size(r)
               radial%phi_large(ir, l + 1, spin) = (0.8_rp + 0.05_rp*spin)*r(ir)**l* &
                  (1.0_rp + 0.14_rp*r(ir))*exp(-alpha*r(ir))
               radial%phidot_large(ir, l + 1, spin) = (0.31_rp + 0.03_rp*l)*r(ir)**l* &
                  (1.0_rp + 0.22_rp*r(ir))*exp(-0.57_rp*alpha*r(ir))
               radial%phi_small(ir, l + 1, spin) = small_scale*r(ir)**(l + 1)*exp(-1.11_rp*alpha*r(ir))
               radial%phidot_small(ir, l + 1, spin) = 0.43_rp*small_scale*r(ir)**(l + 1)* &
                  (1.0_rp + 0.1_rp*r(ir))*exp(-0.77_rp*alpha*r(ir))
               if (ir == 1 .or. r(ir) <= tiny(1.0_rp)) then
                  radial%gfac(ir, l + 1, spin) = 1.0_rp
                  radial%tmc(ir, l + 1, spin) = scalar_relativistic_c
               else
                  tmc = scalar_relativistic_c - (radial%potential(ir, spin) - 2.0_rp*nuclear_z/r(ir) - &
                     radial%enu_radial(l + 1, spin))/scalar_relativistic_c
                  radial%tmc(ir, l + 1, spin) = tmc
                  radial%gfac(ir, l + 1, spin) = 1.0_rp + real(l*(l + 1), rp)/(tmc*r(ir))**2
               end if
            end do
         end do
      end do
   end subroutine build_radial

   subroutine independent_projector_norms(radial, energy, norms)
      type(lmto_radial_basis), intent(in) :: radial
      real(rp), intent(in) :: energy
      real(rp), intent(out) :: norms(2)
      integer :: spin
      do spin = 1, 2
         call independent_norm(radial, energy, spin, norms(spin))
      end do
   end subroutine independent_projector_norms

   subroutine independent_norm(radial, energy, spin, norm)
      type(lmto_radial_basis), intent(in) :: radial
      real(rp), intent(in) :: energy
      integer, intent(in) :: spin
      real(rp), intent(out) :: norm
      real(rp) :: measure, tmc, metric, large, small, r
      integer :: ir
      norm = 0.0_rp
      do ir = 1, radial%npoint
         r = radial%rofi(ir)
         measure = simpson_coefficient(ir, radial%npoint)*radial%mesh_a*(r + radial%mesh_b)
         if (ir == 1 .or. r <= tiny(1.0_rp)) then
            metric = 1.0_rp
         else
            tmc = scalar_relativistic_c - (radial%potential(ir, spin) - 2.0_rp*radial%nuclear_z/r - &
               energy)/scalar_relativistic_c
            metric = 1.0_rp + 6.0_rp/(tmc*r)**2
         end if
         large = radial%phi_large(ir, 3, spin) + (energy - radial%enu_work(3, spin))* &
            radial%phidot_large(ir, 3, spin)
         small = radial%phi_small(ir, 3, spin) + (energy - radial%enu_work(3, spin))* &
            radial%phidot_small(ir, 3, spin)
         norm = norm + measure*(metric*large**2 + small**2)
      end do
   end subroutine independent_norm

   subroutine independent_band_projection(radial, energies, vectors, energy, projected)
      type(lmto_radial_basis), intent(in) :: radial
      real(rp), intent(in) :: energies(:), energy
      complex(rp), intent(in) :: vectors(:, :)
      complex(rp), intent(out) :: projected(5, 2, size(energies))
      integer :: ib, spin, im, ir, iorb, basis_index
      real(rp) :: norm, overlap, tmc, metric, measure, r, frozen_large, frozen_small, state_large, state_small
      projected = cmplx(0.0_rp, 0.0_rp, rp)
      do spin = 1, 2
         call independent_norm(radial, energy, spin, norm)
         do ib = 1, size(energies)
            overlap = 0.0_rp
            do ir = 1, radial%npoint
               r = radial%rofi(ir)
               measure = simpson_coefficient(ir, radial%npoint)*radial%mesh_a*(r + radial%mesh_b)
               if (ir == 1 .or. r <= tiny(1.0_rp)) then
                  metric = 1.0_rp
               else
                  tmc = scalar_relativistic_c - (radial%potential(ir, spin) - 2.0_rp*radial%nuclear_z/r - &
                     energy)/scalar_relativistic_c
                  metric = 1.0_rp + 6.0_rp/(tmc*r)**2
               end if
               frozen_large = (radial%phi_large(ir, 3, spin) + (energy - radial%enu_work(3, spin))* &
                  radial%phidot_large(ir, 3, spin))/sqrt(norm)
               frozen_small = (radial%phi_small(ir, 3, spin) + (energy - radial%enu_work(3, spin))* &
                  radial%phidot_small(ir, 3, spin))/sqrt(norm)
               state_large = radial%phi_large(ir, 3, spin) + (energies(ib) - radial%enu_work(3, spin))* &
                  radial%phidot_large(ir, 3, spin)
               state_small = radial%phi_small(ir, 3, spin) + (energies(ib) - radial%enu_work(3, spin))* &
                  radial%phidot_small(ir, 3, spin)
               overlap = overlap + measure*(metric*frozen_large*state_large + frozen_small*state_small)
            end do
            do im = 1, 5
               iorb = 4 + im
               basis_index = (spin - 1)*norb + iorb
               projected(im, spin, ib) = vectors(basis_index, ib)*overlap
            end do
         end do
      end do
   end subroutine independent_band_projection

   subroutine independent_juelich_chi(radial, state, energy, frequencies, eta, chi)
      type(lmto_radial_basis), intent(in) :: radial
      type(lr_electronic_state), intent(in) :: state
      real(rp), intent(in) :: energy, frequencies(:), eta
      complex(rp), intent(out) :: chi(1, 1, size(frequencies))
      complex(rp) :: projected(5, 2, state%nbands), transition
      complex(rp) :: denominator
      real(rp) :: difference, factor
      integer :: ib, jb, iw, im
      call independent_band_projection(radial, state%eigenvalues(:, 1), state%eigenvectors(:, :, 1), energy, projected)
      chi = cmplx(0.0_rp, 0.0_rp, rp)
      do ib = 1, state%nbands
         do jb = 1, state%nbands
            difference = state%occupations(ib, 1) - state%occupations(jb, 1)
            if (difference == 0.0_rp) cycle
            transition = cmplx(0.0_rp, 0.0_rp, rp)
            do im = 1, 5
               transition = transition + conjg(projected(im, 1, ib))*projected(im, 2, jb)
            end do
            do iw = 1, size(frequencies)
               denominator = cmplx(frequencies(iw) + state%eigenvalues(ib, 1) - state%eigenvalues(jb, 1), eta, rp)
               factor = 4.0_rp*response_angular_pi*difference
               chi(1, 1, iw) = chi(1, 1, iw) + factor*transition*conjg(transition)/denominator
            end do
         end do
      end do
   end subroutine independent_juelich_chi

   real(rp) function independent_energy_dependent_overlap(radial, spin_left, energy_left, spin_right, energy_right) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: spin_left, spin_right
      real(rp), intent(in) :: energy_left, energy_right
      real(rp) :: r, measure, left, right
      integer :: ir
      value = 0.0_rp
      do ir = 1, radial%npoint
         r = radial%rofi(ir)
         measure = simpson_coefficient(ir, radial%npoint)*radial%mesh_a*(r + radial%mesh_b)
         left = radial%phi_large(ir, 3, spin_left) + (energy_left - radial%enu_work(3, spin_left))* &
            radial%phidot_large(ir, 3, spin_left)
         right = radial%phi_large(ir, 3, spin_right) + (energy_right - radial%enu_work(3, spin_right))* &
            radial%phidot_large(ir, 3, spin_right)
         value = value + measure*left*right
      end do
   end function independent_energy_dependent_overlap

   pure real(rp) function simpson_coefficient(ir, n)
      integer, intent(in) :: ir, n
      simpson_coefficient = 2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp
      if (ir == 1 .or. ir == n) simpson_coefficient = 1.0_rp/3.0_rp
   end function simpson_coefficient

end module test_juelich_d_reference_mod
