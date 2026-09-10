!------------------------------------------------------------------------------
! LR-02 angular/radial-coordinate and endpoint-gauge oracles.
!
! This test intentionally does not claim a physical scalar-relativistic
! transition-density oracle.  The production radial interface lacks the
! lower-component spin-angular contract required for that mapping; the guard
! at the end verifies that the capability is not silently promoted.
!------------------------------------------------------------------------------
program test_lr_response_basis
   use precision_mod, only: rp
   use response_angular_basis_mod, only: response_angular_pi, response_lmax, &
      response_product_space_complete, response_lm_index, response_lm_from_index, response_harmonic, &
      response_gaunt, response_product_coefficients
   use response_basis_mapping_mod, only: response_super_index, response_superindex_size, &
      response_flatten_superindex, response_unflatten_superindex, response_simpson_weight, &
      response_log_mesh_jacobian, response_volume_measure, response_weighted_density, &
      response_physical_density, response_log_mesh_integral, response_endpoint_phase, &
      response_apply_site_gauge, response_mapping_supported
   implicit none

   logical :: failed

   failed = .false.
   call angular_product_oracle(failed)
   call cutoff_oracle(failed)
   call charge_spin_l0_oracle(failed)
   call radial_measure_oracle(failed)
   call superindex_oracle(failed)
   call endpoint_gauge_oracle(failed)
   call capability_guard_oracle(failed)

   if (failed) then
      write (*, '(a)') 'UnitLrResponseBasis: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrResponseBasis: PASS (independent angular/radial/gauge oracles; physical mapping guarded)'

contains

   subroutine angular_product_oracle(test_failed)
      logical, intent(inout) :: test_failed
      call check_product(0, 0, 0, 0, test_failed) ! s x s
      call check_product(0, 0, 1, 0, test_failed) ! s x p
      call check_product(1, -1, 1, 1, test_failed) ! p x p off-diagonal
      call check_product(1, 0, 1, 0, test_failed) ! p x p diagonal
      call check_product(1, -1, 2, 1, test_failed) ! p x d off-diagonal
      call check_product(1, 0, 2, 0, test_failed) ! p x d diagonal
      call check_product(2, -2, 2, 2, test_failed) ! d x d off-diagonal
      call check_product(2, 1, 2, 1, test_failed) ! d x d diagonal
   end subroutine angular_product_oracle

   subroutine check_product(l, m, lp, mp, test_failed)
      integer, intent(in) :: l, m, lp, mp
      logical, intent(inout) :: test_failed
      integer, parameter :: nx = 32, nphi = 64
      real(rp) :: x(nx), wx(nx), theta, phi, weight
      real(rp) :: analytic((l + lp + 1)**2), numerical((l + lp + 1)**2)
      real(rp) :: max_coeff_error, max_point_error
      complex(rp) :: lhs, yi, yj
      integer :: ix, iphi, lout, mout, index

      call gauss_legendre(nx, x, wx)
      analytic = 0.0_rp
      numerical = 0.0_rp
      call fill_analytic_coefficients(l, m, lp, mp, analytic)
      max_point_error = 0.0_rp
      do ix = 1, nx
         theta = acos(x(ix))
         do iphi = 1, nphi
            phi = 2.0_rp*response_angular_pi*real(iphi - 1, rp)/real(nphi, rp)
            yi = independent_harmonic(l, m, theta, phi)
            yj = independent_harmonic(lp, mp, theta, phi)
            lhs = conjg(yi)*yj
            max_point_error = max(max_point_error, abs(lhs - &
               independent_harmonic_product_from_coefficients(analytic, l + lp, theta, phi)))
            weight = wx(ix)*(2.0_rp*response_angular_pi/real(nphi, rp))
            do lout = 0, l + lp
               do mout = -lout, lout
                  index = response_lm_index(lout, mout)
                  numerical(index) = numerical(index) + weight* &
                     real(conjg(independent_harmonic(lout, mout, theta, phi))*conjg(yi)*yj, rp)
               end do
            end do
         end do
      end do
      max_coeff_error = maxval(abs(analytic - numerical))
      write (*, '(a,4(i2,1x),a,2(es12.4,1x))') 'Gaunt pair ', l, m, lp, mp, &
         'coeff/point errors = ', max_coeff_error, max_point_error
      if (max_coeff_error > 2.0e-10_rp .or. max_point_error > 2.0e-10_rp) test_failed = .true.
   end subroutine check_product

   subroutine fill_analytic_coefficients(l, m, lp, mp, coefficients)
      integer, intent(in) :: l, m, lp, mp
      real(rp), intent(out) :: coefficients(:)
      call response_product_coefficients(l, m, lp, mp, coefficients)
   end subroutine fill_analytic_coefficients

   function independent_harmonic_product_from_coefficients(coefficients, lmax, theta, phi) result(value)
      real(rp), intent(in) :: coefficients(:), theta, phi
      integer, intent(in) :: lmax
      complex(rp) :: value
      integer :: L, M

      value = cmplx(0.0_rp, 0.0_rp, rp)
      do L = 0, lmax
         do M = -L, L
            value = value + coefficients(response_lm_index(L, M))*response_harmonic(L, M, theta, phi)
         end do
      end do
   end function independent_harmonic_product_from_coefficients

   subroutine cutoff_oracle(test_failed)
      logical, intent(inout) :: test_failed
      if (response_lmax(1) /= 2 .or. response_lmax(2) /= 4) test_failed = .true.
      if (.not. response_product_space_complete(2, 4) .or. response_product_space_complete(2, 3)) test_failed = .true.
      write (*, '(a,i0)') 'Derived spd response Lmax = ', response_lmax(2)
   end subroutine cutoff_oracle

   subroutine charge_spin_l0_oracle(test_failed)
      logical, intent(inout) :: test_failed
      real(rp) :: up_l0, down_l0, charge_l0, spin_z_l0, expected_l0

      up_l0 = response_gaunt(2, -2, 2, -2, 0, 0)
      down_l0 = response_gaunt(2, -2, 2, -2, 0, 0)
      charge_l0 = 0.75_rp*up_l0 + 0.25_rp*down_l0
      spin_z_l0 = 0.75_rp*up_l0 - 0.25_rp*down_l0
      expected_l0 = 1.0_rp/sqrt(4.0_rp*response_angular_pi)
      write (*, '(a,2(es12.4,1x))') 'Algebraic charge/spin-z L0 oracle = ', charge_l0, spin_z_l0
      if (abs(charge_l0 - expected_l0) > 2.0e-13_rp .or. &
          abs(spin_z_l0 - 0.5_rp*expected_l0) > 2.0e-13_rp) test_failed = .true.
   end subroutine charge_spin_l0_oracle

   subroutine radial_measure_oracle(test_failed)
      logical, intent(inout) :: test_failed
      integer, parameter :: nr = 101
      real(rp), parameter :: a = 0.03_rp, b = 0.10_rp
      real(rp) :: radius(nr), physical(nr), weighted(nr), recovered(nr)
      real(rp) :: direct, converted, jacobian, volume, max_recovery
      integer :: ir

      radius(1) = 0.0_rp
      do ir = 2, nr
         radius(ir) = b*(exp(a*real(ir - 1, rp)) - 1.0_rp)
      end do
      do ir = 1, nr
         physical(ir) = 0.4_rp*exp(-radius(ir))
         weighted(ir) = response_weighted_density(radius(ir), physical(ir))
         recovered(ir) = response_physical_density(radius(ir), weighted(ir))
      end do
      recovered(1) = physical(1)
      direct = response_log_mesh_integral(weighted, radius, a, b)
      converted = 0.0_rp
      do ir = 1, nr
         jacobian = response_log_mesh_jacobian(a, b, radius(ir))
         volume = response_volume_measure(a, b, radius(ir))
         converted = converted + response_simpson_weight(ir, nr)*jacobian*weighted(ir)
         if (ir > 1 .and. abs(volume - radius(ir)**2*jacobian) > 1.0e-15_rp) test_failed = .true.
      end do
      max_recovery = maxval(abs(recovered - physical))
      write (*, '(a,2(es12.4,1x),a,es12.4)') 'Radial weighted/log integral pair = ', direct, converted, &
         'density recovery error = ', max_recovery
      if (abs(direct - converted) > 2.0e-13_rp .or. max_recovery > 2.0e-13_rp) test_failed = .true.
   end subroutine radial_measure_oracle

   subroutine superindex_oracle(test_failed)
      logical, intent(inout) :: test_failed
      integer, parameter :: nsite = 2, lmax = 4, npoint = 7, nchannel = 4
      type(response_super_index) :: original, recovered
      integer :: site, L, M, ir, channel, flat

      if (response_superindex_size(nsite, lmax, npoint, nchannel) /= 2*25*7*4) test_failed = .true.
      do flat = 1, (lmax + 1)**2
         call response_lm_from_index(flat, L, M)
         if (response_lm_index(L, M) /= flat) test_failed = .true.
      end do
      do site = 1, nsite
         do L = 0, lmax
            do M = -L, L
               do ir = 1, npoint
                  do channel = 1, nchannel
                     original = response_super_index(site, L, M, ir, channel)
                     call response_flatten_superindex(original, nsite, lmax, npoint, nchannel, flat)
                     call response_unflatten_superindex(flat, nsite, lmax, npoint, nchannel, recovered)
                     if (original%site /= recovered%site .or. original%response_l /= recovered%response_l .or. &
                         original%response_m /= recovered%response_m .or. original%radial_point /= recovered%radial_point .or. &
                         original%channel /= recovered%channel) test_failed = .true.
                  end do
               end do
            end do
         end do
      end do
      write (*, '(a)') 'Response super-index round trip: checked'
   end subroutine superindex_oracle

   subroutine endpoint_gauge_oracle(test_failed)
      logical, intent(inout) :: test_failed
      real(rp) :: tau(3, 2), G(3), angle, k(3), q(3), kq_raw(3), kq_fold(3), Gcross(3)
      complex(rp) :: coefficients(2, 2), transformed(2, 2), expected(2, 2)
      complex(rp) :: H(2, 2), Hfold(2, 2), D(2, 2)

      tau(:, 1) = [0.0_rp, 0.0_rp, 0.0_rp]
      tau(:, 2) = [0.25_rp, 0.125_rp, 0.0_rp]
      G = [1.0_rp, 0.0_rp, 0.0_rp]
      k = [0.49_rp, 0.10_rp, 0.0_rp]
      q = [0.62_rp, 0.0_rp, 0.0_rp]
      kq_raw = k + q
      kq_fold = kq_raw - floor(kq_raw + 0.5_rp)
      Gcross = floor(kq_raw + 0.5_rp)
      if (maxval(abs(kq_raw - kq_fold - Gcross)) > 2.0e-14_rp .or. &
          maxval(abs(Gcross - G)) > 2.0e-14_rp) test_failed = .true.
      coefficients = reshape([cmplx(0.2_rp, 0.1_rp, rp), cmplx(-0.3_rp, 0.4_rp, rp), &
         cmplx(0.5_rp, -0.2_rp, rp), cmplx(0.7_rp, 0.3_rp, rp)], [2, 2])
      call response_apply_site_gauge(coefficients, tau, G, transformed)
      expected = coefficients
      expected(:, 2) = response_endpoint_phase(tau(:, 2), G)*expected(:, 2)
      if (maxval(abs(transformed - expected)) > 2.0e-14_rp) test_failed = .true.
      angle = -2.0_rp*response_angular_pi*0.25_rp
      if (abs(response_endpoint_phase(tau(:, 2), G) - cmplx(cos(angle), sin(angle), rp)) > 2.0e-14_rp) test_failed = .true.

      H = reshape([cmplx(1.0_rp, 0.0_rp, rp), cmplx(0.2_rp, -0.1_rp, rp), &
         cmplx(0.2_rp, 0.1_rp, rp), cmplx(2.0_rp, 0.0_rp, rp)], [2, 2])
      D = cmplx(0.0_rp, 0.0_rp, rp)
      D(1, 1) = response_endpoint_phase(tau(:, 1), G)
      D(2, 2) = response_endpoint_phase(tau(:, 2), G)
      Hfold = matmul(D, matmul(H, conjg(transpose(D))))
      if (abs(Hfold(1, 2) - H(1, 2)*conjg(D(2, 2))) > 2.0e-14_rp) test_failed = .true.
      write (*, '(a,2(es12.4,1x))') 'Two-site endpoint gauge phase/H relation = ', &
         abs(response_endpoint_phase(tau(:, 2), G)), abs(Hfold(1, 2) - H(1, 2)*conjg(D(2, 2)))
   end subroutine endpoint_gauge_oracle

   subroutine capability_guard_oracle(test_failed)
      logical, intent(inout) :: test_failed
      if (response_mapping_supported(2, .true., .true., .true., .true., .false.)) test_failed = .true.
      if (.not. response_mapping_supported(2, .true., .true., .true., .true., .true.)) test_failed = .true.
      if (response_mapping_supported(3, .true., .true., .true., .true., .true.)) test_failed = .true.
      write (*, '(a)') 'Unsupported scalar-relativistic angular completion: guarded'
   end subroutine capability_guard_oracle

   pure function independent_harmonic(l, m, theta, phi) result(value)
      integer, intent(in) :: l, m
      real(rp), intent(in) :: theta, phi
      complex(rp) :: value, standard
      real(rp) :: p, norm
      integer :: ma

      ma = abs(m)
      p = independent_legendre(l, ma, cos(theta))
      norm = sqrt(real(2*l + 1, rp)/(4.0_rp*response_angular_pi)* &
         independent_factorial(l - ma)/independent_factorial(l + ma))
      if (m >= 0) then
         standard = cmplx(norm*independent_sign(ma)*p*cos(real(m, rp)*phi), &
            norm*independent_sign(ma)*p*sin(real(m, rp)*phi), rp)
      else
         standard = independent_sign(ma)*conjg(cmplx(norm*independent_sign(ma)*p*cos(real(ma, rp)*phi), &
            norm*independent_sign(ma)*p*sin(real(ma, rp)*phi), rp))
      end if
      value = conjg(standard)
   end function independent_harmonic

   pure real(rp) function independent_factorial(n) result(value)
      integer, intent(in) :: n
      integer :: i
      value = 1.0_rp
      do i = 2, n
         value = value*real(i, rp)
      end do
   end function independent_factorial

   pure real(rp) function independent_sign(n) result(value)
      integer, intent(in) :: n
      value = merge(1.0_rp, -1.0_rp, mod(abs(n), 2) == 0)
   end function independent_sign

   pure real(rp) function independent_legendre(l, m, x) result(value)
      integer, intent(in) :: l, m
      real(rp), intent(in) :: x
      real(rp) :: pmm, pmmp1, pll, somx2
      integer :: ell
      pmm = 1.0_rp
      if (m > 0) then
         somx2 = sqrt(max(0.0_rp, (1.0_rp - x)*(1.0_rp + x)))
         do ell = 1, m
            pmm = pmm*real(2*ell - 1, rp)*somx2
         end do
      end if
      if (l == m) then
         value = pmm
      else
         pmmp1 = x*real(2*m + 1, rp)*pmm
         if (l == m + 1) then
            value = pmmp1
         else
            do ell = m + 2, l
               pll = (real(2*ell - 1, rp)*x*pmmp1 - real(ell + m - 1, rp)*pmm)/real(ell - m, rp)
               pmm = pmmp1
               pmmp1 = pll
            end do
            value = pmmp1
         end if
      end if
   end function independent_legendre

   subroutine gauss_legendre(n, nodes, weights)
      integer, intent(in) :: n
      real(rp), intent(out) :: nodes(n), weights(n)
      integer :: i, j, half
      real(rp) :: z, z1, p1, p2, p3, pp

      half = (n + 1)/2
      do i = 1, half
         z = cos(response_angular_pi*(real(i, rp) - 0.25_rp)/(real(n, rp) + 0.5_rp))
         do
            p1 = 1.0_rp
            p2 = 0.0_rp
            do j = 1, n
               p3 = p2
               p2 = p1
               p1 = ((2.0_rp*real(j, rp) - 1.0_rp)*z*p2 - real(j - 1, rp)*p3)/real(j, rp)
            end do
            pp = real(n, rp)*(z*p1 - p2)/(z*z - 1.0_rp)
            z1 = z
            z = z1 - p1/pp
            if (abs(z - z1) < 2.0e-15_rp) exit
         end do
         nodes(i) = -z
         nodes(n + 1 - i) = z
         weights(i) = 2.0_rp/((1.0_rp - z*z)*pp*pp)
         weights(n + 1 - i) = weights(i)
      end do
   end subroutine gauss_legendre

end program test_lr_response_basis
