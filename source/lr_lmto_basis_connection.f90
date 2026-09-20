!------------------------------------------------------------------------------
! DRESP-09T moving-basis / representation connection.
!
! This module owns the radial part of the connection independently of the
! native LMTO tangent.  In particular, it does not obtain a connection by
! subtracting a fixed-field result from a native result.
!
! The arrays returned by lmto_build_basis_connection are raw matrix elements
!
!     Gamma_raw(ab) = <chi_a | d chi_b / d theta>
!
! for the four linearisation endpoint branches.  If the retained basis has a
! non-unit Gram matrix S, the coefficient-space connection is K=S^{-1}Gamma_raw.
! This distinction is intentional: Gamma_raw is not in general
! anti-Hermitian, while K is anti-Hermitian only after the relevant
! representation has actually been orthogonalised.
!------------------------------------------------------------------------------
module lr_lmto_basis_connection_mod

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis, lmto_orbital_l
   implicit none
   private

   integer, parameter, public :: lmto_connection_branch_count = 4

   public :: lmto_build_basis_connection
   public :: lmto_build_finite_basis_overlap
   public :: lmto_build_effective_connection
   public :: lmto_build_effective_endpoint_matrix
   public :: lmto_covariant_connection
   public :: lmto_basis_tangent
   public :: lmto_build_second_order_basis_tangent
   public :: lmto_connection_metric_residual
   public :: lmto_connection_relative_residual

contains

   !> Build the raw connection and the branch Gram matrices.
   !>
   !> The lower scalar-relativistic contribution here is the overlap metric
   !> of the augmented spin-angular basis.  It contains the mixed-spin
   !> TMC-up/TMC-down angular term.  The -1/3 transverse factor belongs to the
   !> physical transverse field insertion and is deliberately not reused for
   !> a basis overlap.  With identical up/down radial channels this makes the
   !> coefficient connection the ordinary spin generator, as required.
   subroutine lmto_build_basis_connection(radial, radial_weights, gamma_raw, overlap)
      type(lmto_radial_basis), intent(in) :: radial
      real(rp), intent(in) :: radial_weights(:)
      complex(rp), intent(out) :: gamma_raw(:, :, :), overlap(:, :, :)
      integer :: norb, branch, p, q, iorb, l, up, down
      real(rp) :: same_up, same_down, mixed_ud, mixed_du

      norb = (radial%lmax + 1)**2
      if (size(radial_weights) /= radial%npoint .or. any(shape(gamma_raw) /= [2*norb, 2*norb, 4]) .or. &
          any(shape(overlap) /= [2*norb, 2*norb, 4])) then
         error stop 'DRESP-09T basis connection: inconsistent fixture dimensions'
      end if

      gamma_raw = cmplx(0.0_rp, 0.0_rp, rp)
      overlap = cmplx(0.0_rp, 0.0_rp, rp)
      do q = 0, 1
         do p = 0, 1
            branch = 1 + p + 2*q
            do iorb = 1, norb
               l = lmto_orbital_l(iorb)
               up = iorb
               down = norb + iorb
               same_up = radial_overlap(radial, radial_weights, l, 1, p, 1, q)
               same_down = radial_overlap(radial, radial_weights, l, 2, p, 2, q)
               mixed_ud = radial_overlap(radial, radial_weights, l, 1, p, 2, q)
               mixed_du = radial_overlap(radial, radial_weights, l, 2, p, 1, q)
               overlap(up, up, branch) = cmplx(same_up, 0.0_rp, rp)
               overlap(down, down, branch) = cmplx(same_down, 0.0_rp, rp)
               ! dU/dtheta for a positive y rotation is [[0,-1/2],[1/2,0]].
               gamma_raw(up, down, branch) = cmplx(-0.5_rp*mixed_ud, 0.0_rp, rp)
               gamma_raw(down, up, branch) = cmplx(0.5_rp*mixed_du, 0.0_rp, rp)
            end do
         end do
      end do
   end subroutine lmto_build_basis_connection

   !> Independent finite-angle cross-overlap oracle for all endpoint branches.
   !> Only the spin-angular part is rotated; no radial equation is rerun.
   subroutine lmto_build_finite_basis_overlap(radial, radial_weights, theta, overlap_theta)
      type(lmto_radial_basis), intent(in) :: radial
      real(rp), intent(in) :: radial_weights(:), theta
      complex(rp), intent(out) :: overlap_theta(:, :, :)
      integer :: norb, branch, p, q, iorb, l, up, down
      real(rp) :: c, s, same_up, same_down, mixed_ud, mixed_du

      norb = (radial%lmax + 1)**2
      if (size(radial_weights) /= radial%npoint .or. any(shape(overlap_theta) /= [2*norb, 2*norb, 4])) then
         error stop 'DRESP-09T finite basis overlap: inconsistent dimensions'
      end if
      c = cos(0.5_rp*theta)
      s = sin(0.5_rp*theta)
      overlap_theta = cmplx(0.0_rp, 0.0_rp, rp)
      do q = 0, 1
         do p = 0, 1
            branch = 1 + p + 2*q
            do iorb = 1, norb
               l = lmto_orbital_l(iorb)
               up = iorb
               down = norb + iorb
               same_up = radial_overlap(radial, radial_weights, l, 1, p, 1, q)
               same_down = radial_overlap(radial, radial_weights, l, 2, p, 2, q)
               mixed_ud = radial_overlap(radial, radial_weights, l, 1, p, 2, q)
               mixed_du = radial_overlap(radial, radial_weights, l, 2, p, 1, q)
               overlap_theta(up, up, branch) = cmplx(c*same_up, 0.0_rp, rp)
               overlap_theta(down, down, branch) = cmplx(c*same_down, 0.0_rp, rp)
               overlap_theta(up, down, branch) = cmplx(-s*mixed_ud, 0.0_rp, rp)
               overlap_theta(down, up, branch) = cmplx(s*mixed_du, 0.0_rp, rp)
            end do
         end do
      end do
   end subroutine lmto_build_finite_basis_overlap

   !> Contract endpoint branches into the coefficient-space connection.
   !> The order is the live linearisation order: H**q Gamma_pq H**p.
   subroutine lmto_build_effective_connection(gamma_raw, hamiltonian, gamma_effective)
      complex(rp), intent(in) :: gamma_raw(:, :, :), hamiltonian(:, :)
      complex(rp), intent(out) :: gamma_effective(:, :)
      complex(rp), allocatable :: term(:, :)
      integer :: n, p, q, branch

      n = size(hamiltonian, 1)
      if (size(hamiltonian, 2) /= n .or. any(shape(gamma_raw) /= [n, n, 4]) .or. any(shape(gamma_effective) /= [n, n])) then
         error stop 'DRESP-09T effective connection: inconsistent dimensions'
      end if
      allocate(term(n, n))
      gamma_effective = cmplx(0.0_rp, 0.0_rp, rp)
      do q = 0, 1
         do p = 0, 1
            branch = 1 + p + 2*q
            term = gamma_raw(:, :, branch)
            if (q == 1) term = matmul(hamiltonian, term)
            if (p == 1) term = matmul(term, hamiltonian)
            gamma_effective = gamma_effective + term
         end do
      end do
      deallocate(term)
   end subroutine lmto_build_effective_connection

   !> Contract any endpoint matrix with the same H**q X_pq H**p order used
   !> for the linearised augmented basis.  This is used for both S and Gamma.
   subroutine lmto_build_effective_endpoint_matrix(endpoint, hamiltonian, effective)
      complex(rp), intent(in) :: endpoint(:, :, :), hamiltonian(:, :)
      complex(rp), intent(out) :: effective(:, :)
      complex(rp), allocatable :: term(:, :)
      integer :: n, p, q, branch

      n = size(hamiltonian, 1)
      if (size(hamiltonian, 2) /= n .or. any(shape(endpoint) /= [n, n, 4]) .or. any(shape(effective) /= [n, n])) then
         error stop 'DRESP-09T endpoint contraction: inconsistent dimensions'
      end if
      allocate(term(n, n))
      effective = cmplx(0.0_rp, 0.0_rp, rp)
      do q = 0, 1
         do p = 0, 1
            branch = 1 + p + 2*q
            term = endpoint(:, :, branch)
            if (q == 1) term = matmul(hamiltonian, term)
            if (p == 1) term = matmul(term, hamiltonian)
            effective = effective + term
         end do
      end do
      deallocate(term)
   end subroutine lmto_build_effective_endpoint_matrix

   !> Convert raw overlap matrix elements into a coefficient connection.
   !> This solves S*K=Gamma_raw, column by column, and is valid for a
   !> genuinely nonorthogonal reference representation.
   subroutine lmto_covariant_connection(overlap, gamma_raw, connection)
      complex(rp), intent(in) :: overlap(:, :), gamma_raw(:, :)
      complex(rp), intent(out) :: connection(:, :)
      complex(rp), allocatable :: a(:, :), b(:, :)
      integer :: n

      n = size(overlap, 1)
      if (size(overlap, 2) /= n .or. any(shape(gamma_raw) /= [n, n]) .or. any(shape(connection) /= [n, n])) then
         error stop 'DRESP-09T covariant connection: inconsistent dimensions'
      end if
      allocate(a(n, n), b(n, n))
      a = overlap
      b = gamma_raw
      call solve_linear_system(a, b)
      connection = b
      deallocate(a, b)
   end subroutine lmto_covariant_connection

   !> Moving-basis contribution for a matrix in the same coefficient basis.
   !> For S=I this is the familiar Gamma^H H + H Gamma expression.
   subroutine lmto_basis_tangent(connection, matrix, tangent)
      complex(rp), intent(in) :: connection(:, :), matrix(:, :)
      complex(rp), intent(out) :: tangent(:, :)
      integer :: n

      n = size(matrix, 1)
      if (size(matrix, 2) /= n .or. any(shape(connection) /= [n, n]) .or. any(shape(tangent) /= [n, n])) then
         error stop 'DRESP-09T basis tangent: inconsistent dimensions'
      end if
      tangent = matmul(transpose(conjg(connection)), matrix) + matmul(matrix, connection)
   end subroutine lmto_basis_tangent

   !> Product-rule basis tangent for H2=E+h-h*o*h.
   subroutine lmto_build_second_order_basis_tangent(connection_e, connection_h, connection_o, enu, h, o, tangent, &
                                                    tangent_enu, tangent_h, tangent_left, tangent_middle, tangent_right)
      complex(rp), intent(in) :: connection_e(:, :), connection_h(:, :), connection_o(:, :), enu(:, :), h(:, :), o(:, :)
      complex(rp), intent(out) :: tangent(:, :), tangent_enu(:, :), tangent_h(:, :), tangent_left(:, :), &
         tangent_middle(:, :), tangent_right(:, :)
      complex(rp), allocatable :: de(:, :), dh(:, :), do_( :, :)
      integer :: n

      n = size(h, 1)
      if (size(h, 2) /= n .or. any(shape(enu) /= [n, n]) .or. any(shape(o) /= [n, n]) .or. &
          any(shape(connection_e) /= [n, n]) .or. any(shape(connection_h) /= [n, n]) .or. &
          any(shape(connection_o) /= [n, n]) .or. any(shape(tangent) /= [n, n])) then
         error stop 'DRESP-09T second-order basis tangent: inconsistent dimensions'
      end if
      allocate(de(n, n), dh(n, n), do_(n, n))
      call lmto_basis_tangent(connection_e, enu, de)
      call lmto_basis_tangent(connection_h, h, dh)
      call lmto_basis_tangent(connection_o, o, do_)
      tangent_enu = de
      tangent_h = dh
      tangent_left = matmul(matmul(dh, o), h)
      tangent_middle = matmul(matmul(h, do_), h)
      tangent_right = matmul(matmul(h, o), dh)
      tangent = de + dh - tangent_left - tangent_middle - tangent_right
      deallocate(de, dh, do_)
   end subroutine lmto_build_second_order_basis_tangent

   !> Check the metric-covariant identity dS=K^H S+S K.
   real(rp) function lmto_connection_metric_residual(overlap, connection, delta_overlap) result(value)
      complex(rp), intent(in) :: overlap(:, :), connection(:, :), delta_overlap(:, :)
      complex(rp), allocatable :: predicted(:, :)
      real(rp) :: denominator
      integer :: n

      n = size(overlap, 1)
      if (size(overlap, 2) /= n .or. any(shape(connection) /= [n, n]) .or. any(shape(delta_overlap) /= [n, n])) then
         error stop 'DRESP-09T metric residual: inconsistent dimensions'
      end if
      allocate(predicted(n, n))
      predicted = matmul(transpose(conjg(connection)), overlap) + matmul(overlap, connection)
      denominator = max(sqrt(sum(abs(delta_overlap)**2)), tiny(1.0_rp))
      value = sqrt(sum(abs(predicted - delta_overlap)**2))/denominator
      deallocate(predicted)
   end function lmto_connection_metric_residual

   real(rp) function lmto_connection_relative_residual(left, right) result(value)
      complex(rp), intent(in) :: left(:, :), right(:, :)
      real(rp) :: denominator
      if (any(shape(left) /= shape(right))) error stop 'DRESP-09T connection residual: shape mismatch'
      denominator = max(sqrt(sum(abs(right)**2)), tiny(1.0_rp))
      value = sqrt(sum(abs(left - right)**2))/denominator
   end function lmto_connection_relative_residual

   real(rp) function radial_overlap(radial, weights, l, spin_left, p, spin_right, q) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      real(rp), intent(in) :: weights(:)
      integer, intent(in) :: l, spin_left, p, spin_right, q
      integer :: ir
      real(rp) :: r, left_large, right_large, left_small, right_small, tleft, tright

      value = 0.0_rp
      do ir = 2, radial%npoint
         r = radial%rofi(ir)
         if (r <= tiny(1.0_rp)) cycle
         left_large = endpoint(radial%phi_large(:, l + 1, spin_left), radial%phidot_large(:, l + 1, spin_left), &
            radial%enu_work(l + 1, spin_left), ir, p)
         right_large = endpoint(radial%phi_large(:, l + 1, spin_right), radial%phidot_large(:, l + 1, spin_right), &
            radial%enu_work(l + 1, spin_right), ir, q)
         left_small = endpoint(radial%phi_small(:, l + 1, spin_left), radial%phidot_small(:, l + 1, spin_left), &
            radial%enu_work(l + 1, spin_left), ir, p)
         right_small = endpoint(radial%phi_small(:, l + 1, spin_right), radial%phidot_small(:, l + 1, spin_right), &
            radial%enu_work(l + 1, spin_right), ir, q)
         tleft = radial%tmc(ir, l + 1, spin_left)
         tright = radial%tmc(ir, l + 1, spin_right)
         value = value + weights(ir)/r**2*(left_large*right_large*(1.0_rp + real(l*(l + 1), rp)/ &
            (tleft*tright*r**2)) + left_small*right_small)
      end do
   end function radial_overlap

   pure real(rp) function endpoint(phi, phidot, enu, ir, power) result(value)
      real(rp), intent(in) :: phi(:), phidot(:), enu
      integer, intent(in) :: ir, power
      if (power == 0) then
         value = phi(ir) - enu*phidot(ir)
      else if (power == 1) then
         value = phidot(ir)
      else
         error stop 'DRESP-09T endpoint branch must be 0 or 1'
      end if
   end function endpoint

   ! Small dense complex Gaussian solve.  The production connection uses
   ! small site-local blocks; keeping this local avoids silently assuming a
   ! diagonal metric in the public covariant API.
   subroutine solve_linear_system(a, b)
      complex(rp), intent(inout) :: a(:, :), b(:, :)
      integer :: n, i, j, pivot
      complex(rp) :: factor, pivot_value
      real(rp) :: pivot_norm
      complex(rp), allocatable :: row(:)

      n = size(a, 1)
      allocate(row(n))
      do i = 1, n
         pivot = i
         pivot_norm = abs(a(i, i))
         do j = i + 1, n
            if (abs(a(j, i)) > pivot_norm) then
               pivot = j
               pivot_norm = abs(a(j, i))
            end if
         end do
         if (pivot_norm <= 100.0_rp*epsilon(1.0_rp)) error stop 'DRESP-09T covariant connection: singular overlap'
         if (pivot /= i) then
            row = a(i, :); a(i, :) = a(pivot, :); a(pivot, :) = row
            row = b(i, :); b(i, :) = b(pivot, :); b(pivot, :) = row
         end if
         pivot_value = a(i, i)
         a(i, :) = a(i, :)/pivot_value
         b(i, :) = b(i, :)/pivot_value
         do j = 1, n
            if (j == i) cycle
            factor = a(j, i)
            a(j, :) = a(j, :) - factor*a(i, :)
            b(j, :) = b(j, :) - factor*b(i, :)
         end do
      end do
      deallocate(row)
   end subroutine solve_linear_system

end module lr_lmto_basis_connection_mod
