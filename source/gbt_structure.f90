module gbt_structure_mod
   use precision_mod, only: rp
   implicit none
   private

   !> A single authoritative description of a GBT site frame.
   !>
   !> `U` maps a spinor from the cell-independent rotating/reference frame to
   !> the laboratory frame.  `R` is the corresponding SO(3) map for vectors.
   !> Both components use the same convention: a positive phase advances the
   !> transverse (x,y) components counter-clockwise about +z.
   type, public :: gbt_frame_t
      complex(rp) :: U(2, 2)
      real(rp) :: R(3, 3)
   end type gbt_frame_t

   public :: gbt_bond_phase
   public :: gbt_reference_frame
   public :: gbt_frame_from_phase
   public :: gbt_frame_for_cell
   public :: gbt_frame_link
   public :: gbt_endpoint_link
   public :: gbt_rotating_to_lab_vector
   public :: gbt_lab_to_rotating_vector
   public :: gbt_rotating_to_lab_spinor
   public :: gbt_lab_to_rotating_spinor
   public :: gbt_rotating_to_lab_density
   public :: gbt_lab_to_rotating_density
   public :: gbt_reference_moment_is_collinear
   public :: gbt_lift_orbital_block
   public :: gbt_contract_collinear

contains

   pure function gbt_bond_phase(q_cart_2pi_over_alat, d_cart, alat) result(alpha)
      real(rp), intent(in) :: q_cart_2pi_over_alat(3), d_cart(3), alat
      real(rp) :: alpha

      if (alat <= 0.0_rp) then
         alpha = huge(1.0_rp)
      else
         alpha = 2.0_rp*acos(-1.0_rp)*(q_cart_2pi_over_alat(1)*d_cart(1) + &
                  q_cart_2pi_over_alat(2)*d_cart(2) + q_cart_2pi_over_alat(3)*d_cart(3))/alat
      end if
   end function gbt_bond_phase

   pure subroutine gbt_reference_frame(theta, phi, frame)
      real(rp), intent(in) :: theta, phi
      complex(rp), intent(out) :: frame(2, 2)
      real(rp) :: c, s
      complex(rp) :: eiphi

      c = cos(0.5_rp*theta)
      s = sin(0.5_rp*theta)
      eiphi = cmplx(cos(phi), sin(phi), rp)
      frame(1, 1) = cmplx(c, 0.0_rp, rp)
      frame(1, 2) = -conjg(eiphi)*s
      frame(2, 1) = eiphi*s
      frame(2, 2) = cmplx(c, 0.0_rp, rp)
   end subroutine gbt_reference_frame

   !> Construct the frame at a translated cell with a supplied phase.
   !>
   !> The reference frame is
   !> `U0 = D(phi) Ry(theta) D(-phi)`.  Translation applies
   !> `D(phase)` on the left, so the associated vector rotation is
   !> `Rz(phase) R0`.  Keeping both representations here prevents the
   !> Hamiltonian, density, and output paths from inventing separate signs.
   pure subroutine gbt_frame_from_phase(theta, phi, phase, frame)
      real(rp), intent(in) :: theta, phi, phase
      type(gbt_frame_t), intent(out) :: frame
      complex(rp) :: reference_u(2, 2), translation_u(2, 2)
      real(rp) :: reference_r(3, 3), translation_r(3, 3)

      call gbt_reference_frame(theta, phi, reference_u)
      call gbt_spin_z_rotation(phase, translation_u)
      frame%U = matmul(translation_u, reference_u)

      call gbt_reference_rotation(theta, phi, reference_r)
      call gbt_z_rotation(phase, translation_r)
      frame%R = matmul(translation_r, reference_r)
   end subroutine gbt_frame_from_phase

   !> Construct the frame for a physical cell translation.
   !> `q_cart_2pi_over_alat` is in the repository’s Cartesian `2*pi/alat`
   !> convention and `r_cart` is a Cartesian length.
   pure subroutine gbt_frame_for_cell(q_cart_2pi_over_alat, r_cart, alat, theta, phi, frame)
      real(rp), intent(in) :: q_cart_2pi_over_alat(3), r_cart(3), alat, theta, phi
      type(gbt_frame_t), intent(out) :: frame
      real(rp) :: phase

      phase = gbt_bond_phase(q_cart_2pi_over_alat, r_cart, alat)
      call gbt_frame_from_phase(theta, phi, phase, frame)
   end subroutine gbt_frame_for_cell

   !> Form the directed link `U_a^dagger U_b` between two authoritative frames.
   pure subroutine gbt_frame_link(frame_a, frame_b, link)
      type(gbt_frame_t), intent(in) :: frame_a, frame_b
      complex(rp), intent(out) :: link(2, 2)

      link = matmul(transpose(conjg(frame_a%U)), frame_b%U)
   end subroutine gbt_frame_link

   pure subroutine gbt_endpoint_link(q_cart_2pi_over_alat, d_cart, alat, theta_a, phi_a, theta_b, phi_b, link)
      real(rp), intent(in) :: q_cart_2pi_over_alat(3), d_cart(3), alat
      real(rp), intent(in) :: theta_a, phi_a, theta_b, phi_b
      complex(rp), intent(out) :: link(2, 2)
      type(gbt_frame_t) :: frame_a, frame_b
      real(rp) :: alpha

      alpha = gbt_bond_phase(q_cart_2pi_over_alat, d_cart, alat)
      call gbt_frame_from_phase(theta_a, phi_a, 0.0_rp, frame_a)
      call gbt_frame_from_phase(theta_b, phi_b, alpha, frame_b)
      call gbt_frame_link(frame_a, frame_b, link)
   end subroutine gbt_endpoint_link

   !> Convert a vector from the rotating frame to the laboratory frame.
   pure subroutine gbt_rotating_to_lab_vector(frame, rotating, lab)
      type(gbt_frame_t), intent(in) :: frame
      real(rp), intent(in) :: rotating(3)
      real(rp), intent(out) :: lab(3)

      lab = matmul(frame%R, rotating)
   end subroutine gbt_rotating_to_lab_vector

   !> Convert a vector from the laboratory frame to the rotating frame.
   pure subroutine gbt_lab_to_rotating_vector(frame, lab, rotating)
      type(gbt_frame_t), intent(in) :: frame
      real(rp), intent(in) :: lab(3)
      real(rp), intent(out) :: rotating(3)

      rotating = matmul(transpose(frame%R), lab)
   end subroutine gbt_lab_to_rotating_vector

   !> Convert a two-component spinor from the rotating frame to the lab frame.
   pure subroutine gbt_rotating_to_lab_spinor(frame, rotating, lab)
      type(gbt_frame_t), intent(in) :: frame
      complex(rp), intent(in) :: rotating(2)
      complex(rp), intent(out) :: lab(2)

      lab = matmul(frame%U, rotating)
   end subroutine gbt_rotating_to_lab_spinor

   !> Convert a two-component spinor from the lab frame to the rotating frame.
   pure subroutine gbt_lab_to_rotating_spinor(frame, lab, rotating)
      type(gbt_frame_t), intent(in) :: frame
      complex(rp), intent(in) :: lab(2)
      complex(rp), intent(out) :: rotating(2)

      rotating = matmul(transpose(conjg(frame%U)), lab)
   end subroutine gbt_lab_to_rotating_spinor

   !> Convert a spin density matrix from the rotating frame to the lab frame.
   !> The density convention is `rho = psi psi^dagger`.
   pure subroutine gbt_rotating_to_lab_density(frame, rotating, lab)
      type(gbt_frame_t), intent(in) :: frame
      complex(rp), intent(in) :: rotating(2, 2)
      complex(rp), intent(out) :: lab(2, 2)

      lab = matmul(frame%U, matmul(rotating, transpose(conjg(frame%U))))
   end subroutine gbt_rotating_to_lab_density

   !> Convert a spin density matrix from the lab frame to the rotating frame.
   pure subroutine gbt_lab_to_rotating_density(frame, lab, rotating)
      type(gbt_frame_t), intent(in) :: frame
      complex(rp), intent(in) :: lab(2, 2)
      complex(rp), intent(out) :: rotating(2, 2)

      rotating = matmul(transpose(conjg(frame%U)), matmul(lab, frame%U))
   end subroutine gbt_lab_to_rotating_density

   !> Whether `moment` is a valid strict GBT rotating-frame reference axis.
   !> The GBT builder consumes only the sign of z for the collinear channels;
   !> transverse cone information belongs to `theta_ss` and the sublattice
   !> reference-angle arrays instead.
   pure logical function gbt_reference_moment_is_collinear(moment, tolerance)
      real(rp), intent(in) :: moment(3), tolerance

      gbt_reference_moment_is_collinear = sqrt(moment(1)**2 + moment(2)**2) <= tolerance
   end function gbt_reference_moment_is_collinear

   pure subroutine gbt_spin_z_rotation(phase, rotation)
      real(rp), intent(in) :: phase
      complex(rp), intent(out) :: rotation(2, 2)

      rotation = cmplx(0.0_rp, 0.0_rp, rp)
      rotation(1, 1) = cmplx(cos(-0.5_rp*phase), sin(-0.5_rp*phase), rp)
      rotation(2, 2) = cmplx(cos( 0.5_rp*phase), sin( 0.5_rp*phase), rp)
   end subroutine gbt_spin_z_rotation

   pure subroutine gbt_z_rotation(angle, rotation)
      real(rp), intent(in) :: angle
      real(rp), intent(out) :: rotation(3, 3)

      rotation = 0.0_rp
      rotation(1, 1) = cos(angle)
      rotation(1, 2) = -sin(angle)
      rotation(2, 1) = sin(angle)
      rotation(2, 2) = cos(angle)
      rotation(3, 3) = 1.0_rp
   end subroutine gbt_z_rotation

   pure subroutine gbt_y_rotation(angle, rotation)
      real(rp), intent(in) :: angle
      real(rp), intent(out) :: rotation(3, 3)

      rotation = 0.0_rp
      rotation(1, 1) = cos(angle)
      rotation(1, 3) = sin(angle)
      rotation(2, 2) = 1.0_rp
      rotation(3, 1) = -sin(angle)
      rotation(3, 3) = cos(angle)
   end subroutine gbt_y_rotation

   pure subroutine gbt_reference_rotation(theta, phi, rotation)
      real(rp), intent(in) :: theta, phi
      real(rp), intent(out) :: rotation(3, 3)
      real(rp) :: rz_phi(3, 3), ry_theta(3, 3), rz_minus_phi(3, 3)

      call gbt_z_rotation(phi, rz_phi)
      call gbt_y_rotation(theta, ry_theta)
      call gbt_z_rotation(-phi, rz_minus_phi)
      rotation = matmul(rz_phi, matmul(ry_theta, rz_minus_phi))
   end subroutine gbt_reference_rotation

   pure subroutine gbt_lift_orbital_block(s_orb, link, s_spinor)
      complex(rp), intent(in) :: s_orb(:, :), link(2, 2)
      complex(rp), intent(out) :: s_spinor(:, :)
      integer :: nrow, ncol, i, j, sa, sb

      nrow = size(s_orb, 1)
      ncol = size(s_orb, 2)
      s_spinor = cmplx(0.0_rp, 0.0_rp, rp)
      do sa = 1, 2
         do sb = 1, 2
            do i = 1, nrow
               do j = 1, ncol
                  s_spinor(i + (sa - 1)*nrow, j + (sb - 1)*ncol) = s_orb(i, j)*link(sa, sb)
               end do
            end do
         end do
      end do
   end subroutine gbt_lift_orbital_block

   pure subroutine gbt_contract_collinear(s_spinor, wx0_a, wx1_a, wx0_b, wx1_b, sign_a, sign_b, h_spinor)
      complex(rp), intent(in) :: s_spinor(:, :)
      complex(rp), intent(in) :: wx0_a(:), wx1_a(:), wx0_b(:), wx1_b(:)
      real(rp), intent(in) :: sign_a, sign_b
      complex(rp), intent(out) :: h_spinor(:, :)
      integer :: norb_a, norb_b, i, j, sa, sb
      complex(rp) :: wa, wb

      norb_a = size(wx0_a)
      norb_b = size(wx0_b)
      h_spinor = cmplx(0.0_rp, 0.0_rp, rp)
      do sa = 1, 2
         do sb = 1, 2
            do i = 1, norb_a
               wa = wx0_a(i) + merge(sign_a, -sign_a, sa == 1)*wx1_a(i)
               do j = 1, norb_b
                  wb = wx0_b(j) + merge(sign_b, -sign_b, sb == 1)*wx1_b(j)
                  h_spinor(i + (sa - 1)*norb_a, j + (sb - 1)*norb_b) = &
                     wa*s_spinor(i + (sa - 1)*norb_a, j + (sb - 1)*norb_b)*wb
               end do
            end do
         end do
      end do
   end subroutine gbt_contract_collinear

end module gbt_structure_mod
