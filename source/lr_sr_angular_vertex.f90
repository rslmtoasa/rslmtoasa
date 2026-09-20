!------------------------------------------------------------------------------
! DRESP-09Z arbitrary-(L,M) scalar-relativistic spin angular algebra.
!
! The response_gaunt convention is deliberately used as the only analytic
! Gaunt primitive.  In that convention a coefficient is
!
!   int Y(LM)^* Y(lm)^* Y(l'm') dOmega .
!
! This is the covariant coefficient used by the existing LR response space.
! The source rank (L,M) is kept distinct from the orbital-product rank (K,Q).
! The latter is produced when the source harmonic is recoupled with the
! Kölling-Harmon rank-2 tensor.
!
! No numerical quadrature is used by the production evaluators below.  The
! quadrature routine is intentionally independent and is exported only as an
! oracle for DRESP-09Z tests.
!------------------------------------------------------------------------------
module lr_sr_angular_vertex_mod

   use precision_mod, only: rp
   use response_angular_basis_mod, only: response_angular_pi, response_gaunt, response_harmonic
   implicit none
   private

   integer, parameter, public :: sr_vertex_component_x = 1
   integer, parameter, public :: sr_vertex_component_y = 2
   integer, parameter, public :: sr_vertex_component_z = 3
   integer, parameter, public :: sr_vertex_branch_00 = 1
   integer, parameter, public :: sr_vertex_branch_10 = 2
   integer, parameter, public :: sr_vertex_branch_01 = 3
   integer, parameter, public :: sr_vertex_branch_11 = 4
   integer, parameter, public :: sr_vertex_branch_20 = 5
   integer, parameter, public :: sr_vertex_branch_02 = 6
   real(rp), parameter, public :: sr_vertex_lower_rank0 = -1.0_rp/3.0_rp

   public :: sr_angular_upper_vertex
   public :: sr_angular_lower_vertex
   public :: sr_angular_lower_linear_vertex
   public :: sr_angular_rank2_integral
   public :: sr_angular_nn_rank2_coefficient
   public :: sr_angular_quadrature_upper
   public :: sr_angular_quadrature_lower
   public :: sr_angular_circular_vertex
   public :: sr_angular_circular_conjugate_residual
   public :: sr_angular_product_rank_amplitudes
   public :: sr_angular_l0_scalar_vertex
   public :: sr_angular_l0_projected_components
   public :: sr_angular_build_vertex
   public :: sr_vertex_branch_name

contains

   !> Ordinary live-harmonic upper component.
   pure complex(rp) function sr_angular_upper_vertex(l, m, lp, mp, source_l, source_m) result(value)
      integer, intent(in) :: l, m, lp, mp, source_l, source_m

      value = cmplx(response_gaunt(l, m, lp, mp, source_l, source_m), 0.0_rp, rp)
   end function sr_angular_upper_vertex

   !> Cartesian lower K-H operator for sigma_i.
   !>
   !>   (sigma.n) sigma_i (sigma.n) = -sigma_i/3
   !>       + 2 sum_jq c_ij,q Y_2q sigma_j .
   !>
   !> The four-harmonic rank-2 integral is recoupled in the order
   !> (source L) x (tensor 2) -> product K, and only then contracted with the
   !> orbital pair.  K is not silently clipped: response_gaunt supplies the
   !> exact zero when K exceeds the retained orbital-product space.
   subroutine sr_angular_lower_vertex(l, m, lp, mp, source_l, source_m, component, value)
      integer, intent(in) :: l, m, lp, mp, source_l, source_m, component
      complex(rp), intent(out) :: value(2, 2)
      complex(rp) :: sigma(2, 2), linear_weights(3)

      call cartesian_sigma(component, sigma)
      linear_weights = cmplx(0.0_rp, 0.0_rp, rp)
      linear_weights(component) = cmplx(1.0_rp, 0.0_rp, rp)
      call sr_angular_lower_linear_vertex(l, m, lp, mp, source_l, source_m, linear_weights, value)
   end subroutine sr_angular_lower_vertex

   !> Lower K-H operator for an arbitrary Cartesian linear combination.
   !> weights=(1/2,+i/2,0) and (1/2,-i/2,0) are sigma_+ and sigma_-.
   subroutine sr_angular_lower_linear_vertex(l, m, lp, mp, source_l, source_m, weights, value)
      integer, intent(in) :: l, m, lp, mp, source_l, source_m
      complex(rp), intent(in) :: weights(3)
      complex(rp), intent(out) :: value(2, 2)
      complex(rp) :: rank0(2, 2), rank2(2, 2)

      call sr_angular_lower_rank0_linear_vertex(l, m, lp, mp, source_l, source_m, weights, rank0)
      call sr_angular_lower_rank2_linear_vertex(l, m, lp, mp, source_l, source_m, weights, rank2)
      value = rank0 + rank2
   end subroutine sr_angular_lower_linear_vertex

   pure subroutine sr_angular_lower_rank0_linear_vertex(l, m, lp, mp, source_l, source_m, weights, value)
      integer, intent(in) :: l, m, lp, mp, source_l, source_m
      complex(rp), intent(in) :: weights(3)
      complex(rp), intent(out) :: value(2, 2)
      complex(rp) :: sigma(2, 2)
      integer :: i

      value = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, 3
         call cartesian_sigma(i, sigma)
         value = value + weights(i)*lower_rank0_matrix(&
            sr_angular_upper_vertex(l, m, lp, mp, source_l, source_m), sigma)
      end do
   end subroutine sr_angular_lower_rank0_linear_vertex

   pure subroutine sr_angular_lower_rank2_linear_vertex(l, m, lp, mp, source_l, source_m, weights, value)
      integer, intent(in) :: l, m, lp, mp, source_l, source_m
      complex(rp), intent(in) :: weights(3)
      complex(rp), intent(out) :: value(2, 2)
      complex(rp) :: sigma(2, 2, 3)
      integer :: i, j, q

      call cartesian_sigma(sr_vertex_component_x, sigma(:, :, 1))
      call cartesian_sigma(sr_vertex_component_y, sigma(:, :, 2))
      call cartesian_sigma(sr_vertex_component_z, sigma(:, :, 3))
      value = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, 3
         do j = 1, 3
            do q = -2, 2
               value = value + 2.0_rp*weights(i)*sr_angular_nn_rank2_coefficient(i, j, q)* &
                  sr_angular_rank2_integral(l, m, lp, mp, source_l, source_m, i, j, q)*sigma(:, :, j)
            end do
         end do
      end do
   end subroutine sr_angular_lower_rank2_linear_vertex

   !> The spatial four-harmonic integral multiplying Q^(2)_ij.
   pure complex(rp) function sr_angular_rank2_integral(l, m, lp, mp, source_l, source_m, i, j, q) result(value)
      integer, intent(in) :: l, m, lp, mp, source_l, source_m, i, j, q
      integer :: product_l, product_m
      complex(rp) :: source_tensor_coeff, orbital_integral

      value = cmplx(0.0_rp, 0.0_rp, rp)
      if (q < -2 .or. q > 2) return
      if (abs(m) > l .or. abs(mp) > lp .or. abs(source_m) > source_l) return
      if (abs(sr_angular_nn_rank2_coefficient(i, j, q)) <= tiny(1.0_rp)) return
      do product_l = 0, source_l + 2
         do product_m = -product_l, product_l
            source_tensor_coeff = cmplx(response_gaunt(source_l, source_m, 2, q, product_l, product_m), 0.0_rp, rp)
            if (abs(source_tensor_coeff) <= tiny(1.0_rp)) cycle
            ! Y_KQ = (-1)^Q conj(Y_K,-Q), which converts the remaining
            ! integral to the authoritative live-harmonic Gaunt primitive.
            orbital_integral = cmplx(real_power_minus_one(product_m), 0.0_rp, rp)* &
               cmplx(response_gaunt(product_l, -product_m, lp, mp, l, m), 0.0_rp, rp)
            value = value + source_tensor_coeff*orbital_integral
         end do
      end do
   end function sr_angular_rank2_integral

   !> Coefficient of the live Y_2q in the traceless part of n_i n_j.
   pure complex(rp) function sr_angular_nn_rank2_coefficient(i, j, q) result(value)
      integer, intent(in) :: i, j, q
      real(rp) :: a, b, c

      value = cmplx(0.0_rp, 0.0_rp, rp)
      if (q < -2 .or. q > 2 .or. i < 1 .or. i > 3 .or. j < 1 .or. j > 3) return
      a = sqrt(15.0_rp/(8.0_rp*response_angular_pi))
      b = sqrt(15.0_rp/(32.0_rp*response_angular_pi))
      c = 1.0_rp/(3.0_rp*sqrt(5.0_rp/(16.0_rp*response_angular_pi)))
      select case (i)
      case (1)
         select case (j)
         case (1); if (q == 0) value = cmplx(-0.5_rp*c, 0.0_rp, rp); if (abs(q) == 2) value = cmplx(1.0_rp/(4.0_rp*b), 0.0_rp, rp)
         case (2); if (q == -2) value = cmplx(0.0_rp, -1.0_rp/(4.0_rp*b), rp); if (q == 2) value = cmplx(0.0_rp, 1.0_rp/(4.0_rp*b), rp)
         case (3); if (q == -1) value = cmplx(1.0_rp/(2.0_rp*a), 0.0_rp, rp); if (q == 1) value = cmplx(-1.0_rp/(2.0_rp*a), 0.0_rp, rp)
         end select
      case (2)
         select case (j)
         case (1); if (q == -2) value = cmplx(0.0_rp, -1.0_rp/(4.0_rp*b), rp); if (q == 2) value = cmplx(0.0_rp, 1.0_rp/(4.0_rp*b), rp)
         case (2); if (q == 0) value = cmplx(-0.5_rp*c, 0.0_rp, rp); if (abs(q) == 2) value = cmplx(-1.0_rp/(4.0_rp*b), 0.0_rp, rp)
         case (3); if (q == -1) value = cmplx(0.0_rp, -1.0_rp/(2.0_rp*a), rp); if (q == 1) value = cmplx(0.0_rp, -1.0_rp/(2.0_rp*a), rp)
         end select
      case (3)
         select case (j)
         case (1); if (q == -1) value = cmplx(1.0_rp/(2.0_rp*a), 0.0_rp, rp); if (q == 1) value = cmplx(-1.0_rp/(2.0_rp*a), 0.0_rp, rp)
         case (2); if (q == -1) value = cmplx(0.0_rp, -1.0_rp/(2.0_rp*a), rp); if (q == 1) value = cmplx(0.0_rp, -1.0_rp/(2.0_rp*a), rp)
         case (3); if (q == 0) value = cmplx(c, 0.0_rp, rp)
         end select
      end select
   end function sr_angular_nn_rank2_coefficient

   !> Independent circular-channel construction from Cartesian operators.
   subroutine sr_angular_circular_vertex(l, m, lp, mp, source_l, source_m, circular_channel, value)
      integer, intent(in) :: l, m, lp, mp, source_l, source_m, circular_channel
      complex(rp), intent(out) :: value(2, 2)
      complex(rp) :: weights(3)

      weights = cmplx(0.0_rp, 0.0_rp, rp)
      if (circular_channel == 1) then
         weights(1) = cmplx(0.5_rp, 0.0_rp, rp)
         weights(2) = cmplx(0.0_rp, 0.5_rp, rp)
      else if (circular_channel == -1) then
         weights(1) = cmplx(0.5_rp, 0.0_rp, rp)
         weights(2) = cmplx(0.0_rp, -0.5_rp, rp)
      else
         error stop 'DRESP-09Z: circular channel must be +1 or -1'
      end if
      call sr_angular_lower_linear_vertex(l, m, lp, mp, source_l, source_m, weights, value)
   end subroutine sr_angular_circular_vertex

   !> Assemble the angular operator in the accepted site-local coefficient
   !> ordering (spin-up orbitals followed by spin-down orbitals).  `weights`
   !> describes sigma_i; use a unit Cartesian vector or the circular weights
   !> documented by sr_angular_circular_vertex.
   subroutine sr_angular_build_vertex(orbital_lmax, source_l, source_m, weights, upper_matrix, lower_matrix)
      integer, intent(in) :: orbital_lmax, source_l, source_m
      complex(rp), intent(in) :: weights(3)
      complex(rp), intent(out) :: upper_matrix(:, :), lower_matrix(:, :)
      integer :: norb, iorb, jorb, spin_left, spin_right, l, lp, m, mp, row, col
      complex(rp) :: sigma(2,2,3), lower(2,2)

      if (orbital_lmax < 0 .or. source_l < 0 .or. abs(source_m) > source_l) then
         error stop 'DRESP-09Z: invalid angular vertex rank'
      end if
      norb = (orbital_lmax+1)**2
      if (any(shape(upper_matrix) /= [2*norb,2*norb]) .or. any(shape(lower_matrix) /= [2*norb,2*norb])) then
         error stop 'DRESP-09Z: angular vertex matrix shape mismatch'
      end if
      call cartesian_sigma(1,sigma(:,:,1)); call cartesian_sigma(2,sigma(:,:,2)); call cartesian_sigma(3,sigma(:,:,3))
      upper_matrix = cmplx(0.0_rp,0.0_rp,rp); lower_matrix = upper_matrix
      do iorb = 1, norb
         l = orbital_l_from_index(iorb); m = iorb-l*l-l-1
         do jorb = 1, norb
            lp = orbital_l_from_index(jorb); mp = jorb-lp*lp-lp-1
            lower = cmplx(0.0_rp,0.0_rp,rp)
            call sr_angular_lower_linear_vertex(l,m,lp,mp,source_l,source_m,weights,lower)
            ! Build the spin matrix for the ordinary upper factor.
            do spin_left = 1, 2
               do spin_right = 1, 2
                  row = (spin_left-1)*norb+iorb; col = (spin_right-1)*norb+jorb
                  upper_matrix(row,col) = sr_angular_upper_vertex(l,m,lp,mp,source_l,source_m)* &
                     sum(weights*[(sigma(spin_left,spin_right,1)),(sigma(spin_left,spin_right,2)), &
                                  (sigma(spin_left,spin_right,3))])
                  lower_matrix(row,col) = lower(spin_left,spin_right)
               end do
            end do
         end do
      end do
   end subroutine sr_angular_build_vertex

   !> Maximum residual of O^-_{L,M}=(-1)^M (O^+_{L,-M})^dagger.
   pure real(rp) function sr_angular_circular_conjugate_residual(plus_lm, minus_lm, source_l, source_m) result(value)
      complex(rp), intent(in) :: plus_lm(2, 2), minus_lm(2, 2)
      integer, intent(in) :: source_l, source_m
      complex(rp) :: expected(2, 2)

      expected = cmplx(real_power_minus_one(source_m), 0.0_rp, rp)*transpose(conjg(plus_lm))
      value = maxval(abs(minus_lm - expected))
      if (source_l < 0) value = huge(1.0_rp)
   end function sr_angular_circular_conjugate_residual

   !> Rank-resolved lower contribution, useful for the K<=4 audit.
   subroutine sr_angular_product_rank_amplitudes(l, m, lp, mp, source_l, source_m, weights, amplitudes)
      integer, intent(in) :: l, m, lp, mp, source_l, source_m
      complex(rp), intent(in) :: weights(3)
      real(rp), intent(out) :: amplitudes(0:4)
      integer :: i, j, q, product_l, product_m
      complex(rp) :: coefficient, orbital_integral

      amplitudes = 0.0_rp
      do i = 1, 3
         do j = 1, 3
            do q = -2, 2
               do product_l = 0, min(4, source_l + 2)
                  do product_m = -product_l, product_l
                     coefficient = cmplx(response_gaunt(source_l, source_m, 2, q, product_l, product_m), 0.0_rp, rp)
                     orbital_integral = cmplx(real_power_minus_one(product_m), 0.0_rp, rp)* &
                        cmplx(response_gaunt(product_l, -product_m, lp, mp, l, m), 0.0_rp, rp)
                     amplitudes(product_l) = amplitudes(product_l) + abs(2.0_rp*weights(i)* &
                        sr_angular_nn_rank2_coefficient(i, j, q)*coefficient*orbital_integral)**2
                  end do
               end do
            end do
         end do
      end do
      amplitudes = sqrt(amplitudes)
   end subroutine sr_angular_product_rank_amplitudes

   !> The certified scalar L=0 reduction used by the DRESP-09X/Y gate.
   pure subroutine sr_angular_l0_scalar_vertex(l, m, lp, mp, sigma, value)
      integer, intent(in) :: l, m, lp, mp
      complex(rp), intent(in) :: sigma(2, 2)
      complex(rp), intent(out) :: value(2, 2)

      value = lower_rank0_matrix(sr_angular_upper_vertex(l, m, lp, mp, 0, 0), sigma)
   end subroutine sr_angular_l0_scalar_vertex

   !> Apply the orbital scalar projector to the complete L=0 operator in one
   !> l shell.  The rank-2 result is retained as a separate diagnostic and is
   !> expected to vanish only after the complete shell trace.
   subroutine sr_angular_l0_projected_components(l, component, upper, lower_small, lower_rank0, lower_rank2, total)
      integer, intent(in) :: l, component
      complex(rp), intent(out) :: upper(2,2), lower_small(2,2), lower_rank0(2,2), lower_rank2(2,2), total(2,2)
      complex(rp) :: sigma(2,2), rank0(2,2), rank2(2,2)
      integer :: m

      call cartesian_sigma(component, sigma)
      upper = cmplx(0.0_rp,0.0_rp,rp); lower_small = upper; lower_rank0 = upper; lower_rank2 = upper
      do m = -l, l
         upper = upper + sr_angular_upper_vertex(l,m,l,m,0,0)*sigma
         call sr_angular_lower_rank0_linear_vertex(l,m,l,m,0,0,unit_cartesian_weights(component),rank0)
         call sr_angular_lower_rank2_linear_vertex(l,m,l,m,0,0,unit_cartesian_weights(component),rank2)
         lower_small = lower_small + rank0
         lower_rank0 = lower_rank0 + real(l*(l+1),rp)*rank0
         lower_rank2 = lower_rank2 + rank2
      end do
      upper = upper/real(2*l+1,rp)
      lower_small = lower_small/real(2*l+1,rp)
      lower_rank0 = lower_rank0/real(2*l+1,rp)
      lower_rank2 = lower_rank2/real(2*l+1,rp)
      total = upper + lower_small + lower_rank0 + lower_rank2
   end subroutine sr_angular_l0_projected_components

   pure function lower_rank0_matrix(upper, sigma) result(value)
      complex(rp), intent(in) :: upper, sigma(2, 2)
      complex(rp) :: value(2, 2)
      value = sr_vertex_lower_rank0*upper*sigma
   end function lower_rank0_matrix

   pure function unit_cartesian_weights(component) result(weights)
      integer, intent(in) :: component
      complex(rp) :: weights(3)
      weights = cmplx(0.0_rp,0.0_rp,rp)
      if (component < 1 .or. component > 3) error stop 'DRESP-09Z: invalid Cartesian component'
      weights(component) = cmplx(1.0_rp,0.0_rp,rp)
   end function unit_cartesian_weights

   pure integer function real_power_minus_one(n) result(value)
      integer, intent(in) :: n
      value = merge(1, -1, mod(abs(n), 2) == 0)
   end function real_power_minus_one

   pure integer function orbital_l_from_index(iorb) result(l)
      integer, intent(in) :: iorb
      integer :: trial
      l = -1
      do trial = 0, 64
         if (trial*trial+1 <= iorb .and. iorb <= (trial+1)*(trial+1)) then
            l = trial
            return
         end if
      end do
   end function orbital_l_from_index

   pure subroutine cartesian_sigma(component, sigma)
      integer, intent(in) :: component
      complex(rp), intent(out) :: sigma(2, 2)
      sigma = cmplx(0.0_rp, 0.0_rp, rp)
      select case (component)
      case (1); sigma(1,2) = 1.0_rp; sigma(2,1) = 1.0_rp
      case (2); sigma(1,2) = cmplx(0.0_rp,-1.0_rp,rp); sigma(2,1) = cmplx(0.0_rp,1.0_rp,rp)
      case (3); sigma(1,1) = 1.0_rp; sigma(2,2) = -1.0_rp
      case default; error stop 'DRESP-09Z: Cartesian spin component must be 1, 2, or 3'
      end select
   end subroutine cartesian_sigma

   !> Independent Gauss-Legendre x trapezoid oracle for the upper vertex.
   function sr_angular_quadrature_upper(l, m, lp, mp, source_l, source_m) result(value)
      integer, intent(in) :: l, m, lp, mp, source_l, source_m
      complex(rp) :: value
      integer, parameter :: nx = 24, nphi = 96
      integer :: ix, ip
      real(rp) :: x(nx), w(nx), phi, pi2
      complex(rp) :: integrand

      call gauss_legendre_nodes(nx, x, w)
      pi2 = 2.0_rp*response_angular_pi
      value = cmplx(0.0_rp, 0.0_rp, rp)
      do ix = 1, nx
         do ip = 0, nphi - 1
            phi = pi2*real(ip, rp)/real(nphi, rp)
            integrand = conjg(response_harmonic(source_l, source_m, acos(x(ix)), phi))* &
               conjg(response_harmonic(l, m, acos(x(ix)), phi))* &
               response_harmonic(lp, mp, acos(x(ix)), phi)
            value = value + w(ix)*pi2/real(nphi, rp)*integrand
         end do
      end do
   end function sr_angular_quadrature_upper

   !> Independent numerical oracle for the full lower K-H Cartesian operator.
   subroutine sr_angular_quadrature_lower(l, m, lp, mp, source_l, source_m, component, value)
      integer, intent(in) :: l, m, lp, mp, source_l, source_m, component
      complex(rp), intent(out) :: value(2, 2)
      integer, parameter :: nx = 24, nphi = 96
      integer :: ix, ip, j
      real(rp) :: x(nx), w(nx), phi, pi2, nxv, nyv, nzv, n(3)
      complex(rp) :: bra, source, ket, sigma(2,2,3), operator(2,2), identity_part(2,2)

      call gauss_legendre_nodes(nx, x, w)
      call cartesian_sigma(1, sigma(:, :, 1)); call cartesian_sigma(2, sigma(:, :, 2)); call cartesian_sigma(3, sigma(:, :, 3))
      pi2 = 2.0_rp*response_angular_pi
      value = cmplx(0.0_rp, 0.0_rp, rp)
      do ix = 1, nx
         do ip = 0, nphi - 1
            phi = pi2*real(ip, rp)/real(nphi, rp)
            nxv = sqrt(max(0.0_rp,1.0_rp-x(ix)*x(ix)))*cos(phi)
            nyv = sqrt(max(0.0_rp,1.0_rp-x(ix)*x(ix)))*sin(phi)
            nzv = x(ix)
            n = [nxv, nyv, nzv]
            operator = cmplx(0.0_rp, 0.0_rp, rp)
            do j = 1, 3
               operator = operator + 2.0_rp*n(component)*n(j)*sigma(:, :, j)
            end do
            call cartesian_sigma(component, identity_part)
            operator = operator - identity_part
            source = conjg(response_harmonic(source_l, source_m, acos(x(ix)), phi))
            bra = conjg(response_harmonic(l, m, acos(x(ix)), phi))
            ket = response_harmonic(lp, mp, acos(x(ix)), phi)
            value = value + w(ix)*pi2/real(nphi, rp)*source*bra*ket*operator
         end do
      end do
   end subroutine sr_angular_quadrature_lower

   subroutine gauss_legendre_nodes(n, nodes, weights)
      integer, intent(in) :: n
      real(rp), intent(out) :: nodes(n), weights(n)
      integer :: i, j, half
      real(rp) :: z, z1, p1, p2, p3, pp

      half = (n + 1)/2
      do i = 1, half
         z = cos(response_angular_pi*(real(i, rp)-0.25_rp)/(real(n, rp)+0.5_rp))
         do
            p1 = 1.0_rp; p2 = 0.0_rp
            do j = 1, n
               p3 = p2; p2 = p1
               p1 = ((2.0_rp*real(j, rp)-1.0_rp)*z*p2 - real(j-1, rp)*p3)/real(j, rp)
            end do
            pp = real(n, rp)*(z*p1-p2)/(z*z-1.0_rp)
            z1 = z; z = z1-p1/pp
            if (abs(z-z1) < 4.0_rp*epsilon(1.0_rp)) exit
         end do
         nodes(i) = -z; nodes(n+1-i) = z
         weights(i) = 2.0_rp/((1.0_rp-z*z)*pp*pp)
         weights(n+1-i) = weights(i)
      end do
   end subroutine gauss_legendre_nodes

   pure function sr_vertex_branch_name(branch) result(name)
      integer, intent(in) :: branch
      character(len=2) :: name
      select case (branch)
      case (sr_vertex_branch_00); name = '00'
      case (sr_vertex_branch_10); name = '10'
      case (sr_vertex_branch_01); name = '01'
      case (sr_vertex_branch_11); name = '11'
      case (sr_vertex_branch_20); name = '20'
      case (sr_vertex_branch_02); name = '02'
      case default; name = '??'
      end select
   end function sr_vertex_branch_name

end module lr_sr_angular_vertex_mod
