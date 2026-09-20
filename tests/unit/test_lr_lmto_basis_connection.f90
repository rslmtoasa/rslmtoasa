! DRESP-09T radial moving-basis connection and finite-angle oracle.
!
! The fixture is intentionally independent of the native DRESP-08 tangent.
! It builds distinct up/down scalar-relativistic radial channels, retains all
! four energy-linearisation branches, and checks the raw cross-overlap
! derivative before testing the covariant S^{-1} Gamma connection.
program test_lr_lmto_basis_connection

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_response_space_mod, only: response_space_layout
   use lr_lmto_basis_connection_mod, only: lmto_build_basis_connection, lmto_build_finite_basis_overlap, &
      lmto_covariant_connection, lmto_basis_tangent, lmto_connection_relative_residual, &
      lmto_connection_metric_residual, lmto_build_second_order_basis_tangent
   implicit none

   type(lmto_radial_basis) :: radial
   type(response_space_layout) :: space
   complex(rp), allocatable :: gamma(:, :, :), overlap(:, :, :), plus(:, :, :), minus(:, :, :), fd(:, :, :)
   complex(rp), allocatable :: metric(:, :), raw(:, :), connection(:, :), h(:, :), o(:, :), enu(:, :), h2(:, :), &
      tangent(:, :), second_tangent(:, :), hzero(:, :), dzero(:, :), hp(:, :), hm(:, :)
   complex(rp), allocatable :: scratch1(:, :), scratch2(:, :), scratch3(:, :), scratch4(:, :), scratch5(:, :)
   real(rp), allocatable :: radius(:), potential(:)
   real(rp) :: theta, error_fd(4), error_ratio(4), metric_error, ordinary_error, finite_error(4)
   real(rp) :: hamiltonian_error, radial_value, product_error, nonmagnetic_error
   integer :: npoint, lmax, norb, n, ir, l, spin, p, q, branch, i, j
   real(rp) :: g(18), gp(18), gpp(18), r
   logical :: failed

   npoint = 9
   lmax = 2
   norb = (lmax + 1)**2
   n = 2*norb
   allocate(radius(npoint), potential(npoint))
   do ir = 1, npoint
      radius(ir) = real(ir - 1, rp)*0.12_rp
      potential(ir) = -1.5_rp + 0.04_rp*real(ir - 1, rp)
   end do
   call radial%initialize(npoint, lmax, 2)
   radial%rofi = radius
   radial%mesh_a = 0.08_rp
   radial%mesh_b = 0.7_rp
   radial%nuclear_z = 26.0_rp
   do spin = 1, 2
      do l = 0, lmax
         do ir = 1, npoint
            r = radius(ir)
            g = 0.0_rp
            gp = 0.0_rp
            gpp = 0.0_rp
            g(ir) = (0.82_rp + 0.07_rp*spin + 0.03_rp*l)*max(r, 0.02_rp)**(l + 1)
            gp(ir) = (0.14_rp + 0.02_rp*spin + 0.01_rp*l)*max(r, 0.02_rp)**(l + 1)
            gpp(ir) = (0.02_rp + 0.01_rp*spin)*max(r, 0.02_rp)**(l + 1)
            g(ir + npoint) = (0.11_rp + 0.02_rp*spin + 0.01_rp*l)*max(r, 0.02_rp)**(l + 1)
            gp(ir + npoint) = (0.03_rp + 0.004_rp*spin)*max(r, 0.02_rp)**(l + 1)
            gpp(ir + npoint) = (0.006_rp + 0.001_rp*spin)*max(r, 0.02_rp)**(l + 1)
         end do
         call radial%capture_channel(l, spin, merge(-0.30_rp, -0.17_rp, spin == 1), radius, &
            potential + 0.02_rp*real(spin - 1, rp), radial%mesh_a, radial%mesh_b, radial%nuclear_z, g, gp, gpp, &
            merge(-0.26_rp, -0.11_rp, spin == 1))
      end do
   end do
   call space%initialize(1, 0, radius, radial%mesh_a, radial%mesh_b, 1)

   allocate(gamma(n, n, 4), overlap(n, n, 4), plus(n, n, 4), minus(n, n, 4), fd(n, n, 4))
   call lmto_build_basis_connection(radial, space%radial_weights, gamma, overlap)
   do branch = 1, 4
      theta = 2.0e-3_rp
      call lmto_build_finite_basis_overlap(radial, space%radial_weights, theta, plus)
      call lmto_build_finite_basis_overlap(radial, space%radial_weights, -theta, minus)
      fd(:, :, branch) = (plus(:, :, branch) - minus(:, :, branch))/(2.0_rp*theta)
      error_fd(branch) = lmto_connection_relative_residual(fd(:, :, branch), gamma(:, :, branch))
      theta = theta/2.0_rp
      call lmto_build_finite_basis_overlap(radial, space%radial_weights, theta, plus)
      call lmto_build_finite_basis_overlap(radial, space%radial_weights, -theta, minus)
      fd(:, :, branch) = (plus(:, :, branch) - minus(:, :, branch))/(2.0_rp*theta)
      error_ratio(branch) = error_fd(branch)/max(lmto_connection_relative_residual(fd(:, :, branch), gamma(:, :, branch)), tiny(1.0_rp))
   end do

   ! The raw overlap derivative satisfies the metric identity after the
   ! covariant connection K=S^{-1} Gamma is formed.  This is the nonorthogonal
   ! check; no anti-Hermiticity is imposed on Gamma by projection.
   allocate(metric(n, n), raw(n, n), connection(n, n))
   metric = overlap(:, :, 1)
   raw = gamma(:, :, 1)
   call lmto_covariant_connection(metric, raw, connection)
   metric_error = lmto_connection_metric_residual(metric, connection, transpose(conjg(raw)) + raw)

   ! Equal radial channels must reduce K to the ordinary y-rotation generator.
   radial%phi_large(:, :, 2) = radial%phi_large(:, :, 1)
   radial%phi_small(:, :, 2) = radial%phi_small(:, :, 1)
   radial%phidot_large(:, :, 2) = radial%phidot_large(:, :, 1)
   radial%phidot_small(:, :, 2) = radial%phidot_small(:, :, 1)
   radial%tmc(:, :, 2) = radial%tmc(:, :, 1)
   radial%gfac(:, :, 2) = radial%gfac(:, :, 1)
   radial%enu_work(:, 2) = radial%enu_work(:, 1)
   call lmto_build_basis_connection(radial, space%radial_weights, gamma, overlap)
   metric = overlap(:, :, 1)
   raw = gamma(:, :, 1)
   call lmto_covariant_connection(metric, raw, connection)
   ordinary_error = 0.0_rp
   do i = 1, norb
      ordinary_error = max(ordinary_error, abs(connection(i, norb + i) + 0.5_rp))
      ordinary_error = max(ordinary_error, abs(connection(norb + i, i) - 0.5_rp))
   end do

   ! Independent finite-angle matrix oracle in the orthogonal equal-channel
   ! representation.  The matrix is deliberately noncommuting and complex.
   allocate(h(n, n), o(n, n), enu(n, n), h2(n, n), tangent(n, n), second_tangent(n, n), hzero(n, n), dzero(n, n), &
      hp(n, n), hm(n, n), scratch1(n, n), scratch2(n, n), scratch3(n, n), scratch4(n, n), scratch5(n, n))
   h = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, n
      h(i, i) = cmplx(0.11_rp + 0.017_rp*real(i, rp), 0.0_rp, rp)
      do j = i + 1, n
         h(i, j) = cmplx(0.008_rp*real(i + j, rp), 0.002_rp*real(j - i, rp), rp)
         h(j, i) = conjg(h(i, j))
      end do
   end do
   call lmto_basis_tangent(connection, h, tangent)
   o = cmplx(0.0_rp, 0.0_rp, rp)
   enu = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, n
      o(i,i) = cmplx(0.08_rp + 0.003_rp*real(i, rp), 0.0_rp, rp)
      enu(i,i) = cmplx(-0.17_rp + 0.011_rp*real(i, rp), 0.0_rp, rp)
   end do
   h2 = enu + h - matmul(matmul(h, o), h)
   call lmto_build_second_order_basis_tangent(connection, connection, connection, enu, h, o, second_tangent, &
      scratch1, scratch2, scratch3, scratch4, scratch5)
   product_error = lmto_connection_relative_residual(second_tangent, tangent_for_matrix(connection, h2))
   hzero = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, n/2
      hzero(i,i) = cmplx(0.4_rp + 0.01_rp*real(i, rp), 0.0_rp, rp)
      hzero(n/2+i,n/2+i) = hzero(i,i)
   end do
   call lmto_basis_tangent(connection, hzero, dzero)
   nonmagnetic_error = sqrt(sum(abs(dzero)**2))
   theta = 0.2_rp
   do i = 1, n
      do j = 1, n
         hp(i, j) = h(i, j)
         hm(i, j) = h(i, j)
      end do
   end do
   call rotate_by_connection(h, connection, theta, hp)
   call rotate_by_connection(h, connection, -theta, hm)
   hamiltonian_error = lmto_connection_relative_residual((hp - hm)/(2.0_rp*theta), tangent)
   radial_value = 0.0_rp
   do i = 1, 4
      radial_value = max(radial_value, error_fd(i))
   end do
   finite_error = 0.0_rp
   do i = 1, 4
      finite_error(i) = error_fd(i)
   end do

   failed = maxval(error_fd) > 2.0e-6_rp .or. minval(error_ratio) < 3.0_rp .or. metric_error > 2.0e-12_rp .or. &
      ordinary_error > 2.0e-12_rp .or. hamiltonian_error > 2.0e-2_rp .or. product_error > 2.0e-12_rp .or. &
      nonmagnetic_error > 2.0e-12_rp
   if (failed) then
      write(*, '(a)') 'UnitLrLmtoBasisConnection: FAIL'
      write(*, '(a,4es12.4)') '  branch_fd=', error_fd
      write(*, '(a,4es12.4)') '  branch_theta2_ratio=', error_ratio
      write(*, '(a,es12.4)') '  metric_covariant=', metric_error
      write(*, '(a,es12.4)') '  symmetric_generator=', ordinary_error
      write(*, '(a,es12.4)') '  finite_matrix=', hamiltonian_error
      write(*, '(a,es12.4)') '  product_rule=', product_error
      write(*, '(a,es12.4)') '  nonmagnetic_zero=', nonmagnetic_error
      error stop 1
   end if
   write(*, '(a,4es12.4)') '  branch_fd=', finite_error
   write(*, '(a,4es12.4)') '  branch_theta2_ratio=', error_ratio
   write(*, '(a,es12.4)') '  metric_covariant=', metric_error
   write(*, '(a,es12.4)') '  symmetric_generator=', ordinary_error
   write(*, '(a,es12.4)') '  finite_matrix=', hamiltonian_error
   write(*, '(a,es12.4)') '  product_rule=', product_error
   write(*, '(a,es12.4)') '  nonmagnetic_zero=', nonmagnetic_error
   write(*, '(a)') 'UnitLrLmtoBasisConnection: PASS (00/01/10/11, metric, symmetric, finite-angle)'

contains

   subroutine rotate_by_connection(reference, kappa, theta, rotated)
      complex(rp), intent(in) :: reference(:, :), kappa(:, :)
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: rotated(:, :)
      complex(rp), allocatable :: rotation(:, :)
      integer :: norb, i
      real(rp) :: c, s
      norb = size(reference, 1)/2
      c = cos(0.5_rp*theta)
      s = sin(0.5_rp*theta)
      allocate(rotation(2*norb, 2*norb))
      rotation = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, norb
         rotation(i, i) = cmplx(c, 0.0_rp, rp)
         rotation(norb + i, norb + i) = cmplx(c, 0.0_rp, rp)
         rotation(i, norb + i) = cmplx(-s, 0.0_rp, rp)
         rotation(norb + i, i) = cmplx(s, 0.0_rp, rp)
      end do
      rotated = matmul(transpose(conjg(rotation)), matmul(reference, rotation))
      deallocate(rotation)
   end subroutine rotate_by_connection

   function tangent_for_matrix(kappa, matrix) result(tangent_matrix)
      complex(rp), intent(in) :: kappa(:, :), matrix(:, :)
      complex(rp) :: tangent_matrix(size(matrix,1),size(matrix,2))
      tangent_matrix = matmul(transpose(conjg(kappa)), matrix) + matmul(matrix, kappa)
   end function tangent_for_matrix

end program test_lr_lmto_basis_connection
