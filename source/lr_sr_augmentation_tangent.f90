!------------------------------------------------------------------------------
! DRESP-09Y scalar-relativistic augmentation-frame tangent.
!
! The radial arrays are diagonal in the local collinear spin frame.  This
! module rotates those 2x2 matrices into the global coefficient frame and
! differentiates the physical operator before any density-side contraction.
! It deliberately contains no residual, target subtraction, SCF update, or
! Hamiltonian/potential tangent.
!------------------------------------------------------------------------------
module lr_sr_augmentation_tangent_mod

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   implicit none
   private

   real(rp), parameter, public :: sr_aug_lower_spin_factor = -1.0_rp/3.0_rp
   integer, parameter, public :: sr_aug_upper = 1
   integer, parameter, public :: sr_aug_lower_small = 2
   integer, parameter, public :: sr_aug_lower_angular = 3

   public :: sr_aug_generator_y
   public :: sr_aug_sigma_x
   public :: sr_aug_sigma_y
   public :: sr_aug_sigma_z
   public :: sr_aug_sigma_plus
   public :: sr_aug_sigma_minus
   public :: sr_aug_diagonal_matrix
   public :: sr_aug_rotate_matrix
   public :: sr_aug_matrix_tangent
   public :: sr_aug_observable_tangent
   public :: sr_aug_branch_observable
   public :: sr_aug_branch_observable_at_angle
   public :: sr_aug_branch_observable_tangent

contains

   subroutine sr_aug_generator_y(generator)
      complex(rp), intent(out) :: generator(2, 2)
      generator = cmplx(0.0_rp, 0.0_rp, rp)
      generator(1, 2) = cmplx(0.0_rp, -0.5_rp, rp)
      generator(2, 1) = cmplx(0.0_rp, 0.5_rp, rp)
   end subroutine sr_aug_generator_y

   subroutine sr_aug_sigma_x(sigma)
      complex(rp), intent(out) :: sigma(2, 2)
      sigma = cmplx(0.0_rp, 0.0_rp, rp)
      sigma(1, 2) = 1.0_rp
      sigma(2, 1) = 1.0_rp
   end subroutine sr_aug_sigma_x

   subroutine sr_aug_sigma_y(sigma)
      complex(rp), intent(out) :: sigma(2, 2)
      sigma = cmplx(0.0_rp, 0.0_rp, rp)
      sigma(1, 2) = cmplx(0.0_rp, -1.0_rp, rp)
      sigma(2, 1) = cmplx(0.0_rp, 1.0_rp, rp)
   end subroutine sr_aug_sigma_y

   subroutine sr_aug_sigma_z(sigma)
      complex(rp), intent(out) :: sigma(2, 2)
      sigma = cmplx(0.0_rp, 0.0_rp, rp)
      sigma(1, 1) = 1.0_rp
      sigma(2, 2) = -1.0_rp
   end subroutine sr_aug_sigma_z

   subroutine sr_aug_sigma_plus(sigma)
      complex(rp), intent(out) :: sigma(2, 2)
      sigma = cmplx(0.0_rp, 0.0_rp, rp)
      sigma(1, 2) = 1.0_rp
   end subroutine sr_aug_sigma_plus

   subroutine sr_aug_sigma_minus(sigma)
      complex(rp), intent(out) :: sigma(2, 2)
      sigma = cmplx(0.0_rp, 0.0_rp, rp)
      sigma(2, 1) = 1.0_rp
   end subroutine sr_aug_sigma_minus

   subroutine sr_aug_diagonal_matrix(up, down, matrix)
      real(rp), intent(in) :: up, down
      complex(rp), intent(out) :: matrix(2, 2)
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      matrix(1, 1) = cmplx(up, 0.0_rp, rp)
      matrix(2, 2) = cmplx(down, 0.0_rp, rp)
   end subroutine sr_aug_diagonal_matrix

   subroutine sr_aug_rotate_matrix(generator, theta, matrix, rotated)
      complex(rp), intent(in) :: generator(2, 2), matrix(2, 2)
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: rotated(2, 2)
      complex(rp) :: unitary(2, 2)

      call rotation_from_generator(generator, theta, unitary)
      rotated = matmul(unitary, matmul(matrix, conjg(transpose(unitary))))
   end subroutine sr_aug_rotate_matrix

   !> Analytic derivative of U A U^dagger at theta=0.
   subroutine sr_aug_matrix_tangent(generator, matrix, tangent)
      complex(rp), intent(in) :: generator(2, 2), matrix(2, 2)
      complex(rp), intent(out) :: tangent(2, 2)
      tangent = cmplx(0.0_rp, -1.0_rp, rp)* &
         (matmul(generator, matrix) - matmul(matrix, generator))
   end subroutine sr_aug_matrix_tangent

   !> O_i=A_L sigma_i A_R and its independent augmentation-frame tangent.
   subroutine sr_aug_observable_tangent(generator, left, right, sigma, observable, tangent)
      complex(rp), intent(in) :: generator(2, 2), left(2, 2), right(2, 2), sigma(2, 2)
      complex(rp), intent(out) :: observable(2, 2), tangent(2, 2)
      complex(rp) :: delta_left(2, 2), delta_right(2, 2)

      call sr_aug_matrix_tangent(generator, left, delta_left)
      call sr_aug_matrix_tangent(generator, right, delta_right)
      observable = matmul(left, matmul(sigma, right))
      tangent = matmul(delta_left, matmul(sigma, right)) + matmul(left, matmul(sigma, delta_right))
   end subroutine sr_aug_observable_tangent

   !> Build one absolute-energy radial branch of A_L sigma A_R.
   subroutine sr_aug_branch_observable(radial, ir, l, branch, sigma, observable_upper, observable_small, &
                                       observable_angular, observable_total)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, branch
      complex(rp), intent(in) :: sigma(2, 2)
      complex(rp), intent(out) :: observable_upper(2, 2), observable_small(2, 2), observable_angular(2, 2), &
         observable_total(2, 2)
      complex(rp) :: generator(2, 2)

      call sr_aug_generator_y(generator)
      call branch_observable_with_generator(radial, ir, l, branch, sigma, generator, .false., 0.0_rp, observable_upper, &
         observable_small, observable_angular, observable_total)
   end subroutine sr_aug_branch_observable

   subroutine sr_aug_branch_observable_at_angle(radial, ir, l, branch, theta, sigma, observable_upper, observable_small, &
                                                observable_angular, observable_total)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, branch
      real(rp), intent(in) :: theta
      complex(rp), intent(in) :: sigma(2, 2)
      complex(rp), intent(out) :: observable_upper(2, 2), observable_small(2, 2), observable_angular(2, 2), &
         observable_total(2, 2)
      complex(rp) :: generator(2, 2)

      call sr_aug_generator_y(generator)
      call branch_observable_with_generator(radial, ir, l, branch, sigma, generator, .false., theta, observable_upper, &
         observable_small, observable_angular, observable_total)
   end subroutine sr_aug_branch_observable_at_angle

   !> Analytic derivative of the absolute-energy radial branch observable.
   subroutine sr_aug_branch_observable_tangent(radial, ir, l, branch, sigma, tangent_upper, tangent_small, &
                                               tangent_angular, tangent_total)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, branch
      complex(rp), intent(in) :: sigma(2, 2)
      complex(rp), intent(out) :: tangent_upper(2, 2), tangent_small(2, 2), tangent_angular(2, 2), &
         tangent_total(2, 2)
      complex(rp) :: generator(2, 2)

      call sr_aug_generator_y(generator)
      call branch_observable_with_generator(radial, ir, l, branch, sigma, generator, .true., 0.0_rp, tangent_upper, &
         tangent_small, tangent_angular, tangent_total)
   end subroutine sr_aug_branch_observable_tangent

   subroutine branch_observable_with_generator(radial, ir, l, branch, sigma, generator, derivative, theta, upper, small, &
                                               angular, total)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, branch
      complex(rp), intent(in) :: sigma(2, 2), generator(2, 2)
      logical, intent(in) :: derivative
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: upper(2, 2), small(2, 2), angular(2, 2), total(2, 2)
      integer :: nterm, iterm
      integer :: left_power(6), right_power(6), left_energy_power(6), right_energy_power(6)
      real(rp) :: coefficient(6)
      complex(rp) :: left(2, 2), right(2, 2), observable(2, 2), tangent(2, 2)
      complex(rp) :: rotated_left(2, 2), rotated_right(2, 2)

      upper = cmplx(0.0_rp, 0.0_rp, rp)
      small = upper
      angular = upper
      if (ir == 1 .or. radial%rofi(ir) <= tiny(1.0_rp)) then
         total = upper
         return
      end if
      if (branch < 1 .or. branch > 6) error stop 'DRESP-09Y: invalid observable branch'

      call branch_terms(branch, nterm, left_power, right_power, left_energy_power, right_energy_power, coefficient)
      do iterm = 1, nterm
         call radial_spin_matrix(radial, ir, l, sr_aug_upper, left_power(iterm), left_energy_power(iterm), left)
         call radial_spin_matrix(radial, ir, l, sr_aug_upper, right_power(iterm), right_energy_power(iterm), right)
         if (derivative) then
            call sr_aug_observable_tangent(generator, left, right, sigma, observable, tangent)
            upper = upper + coefficient(iterm)*tangent/radial%rofi(ir)**2
         else
            if (theta /= 0.0_rp) then
               call sr_aug_rotate_matrix(generator, theta, left, rotated_left)
               call sr_aug_rotate_matrix(generator, theta, right, rotated_right)
               left = rotated_left; right = rotated_right
            end if
            observable = matmul(left, matmul(sigma, right))
            upper = upper + coefficient(iterm)*observable/radial%rofi(ir)**2
         end if

         call radial_spin_matrix(radial, ir, l, sr_aug_lower_small, left_power(iterm), left_energy_power(iterm), left)
         call radial_spin_matrix(radial, ir, l, sr_aug_lower_small, right_power(iterm), right_energy_power(iterm), right)
         if (derivative) then
            call sr_aug_observable_tangent(generator, left, right, sigma, observable, tangent)
            small = small + coefficient(iterm)*tangent/radial%rofi(ir)**2
         else
            if (theta /= 0.0_rp) then
               call sr_aug_rotate_matrix(generator, theta, left, rotated_left)
               call sr_aug_rotate_matrix(generator, theta, right, rotated_right)
               left = rotated_left; right = rotated_right
            end if
            observable = matmul(left, matmul(sigma, right))
            small = small + coefficient(iterm)*observable/radial%rofi(ir)**2
         end if

         call radial_spin_matrix(radial, ir, l, sr_aug_lower_angular, left_power(iterm), left_energy_power(iterm), left)
         call radial_spin_matrix(radial, ir, l, sr_aug_lower_angular, right_power(iterm), right_energy_power(iterm), right)
         if (derivative) then
            call sr_aug_observable_tangent(generator, left, right, sigma, observable, tangent)
            angular = angular + coefficient(iterm)*tangent/radial%rofi(ir)**2
         else
            if (theta /= 0.0_rp) then
               call sr_aug_rotate_matrix(generator, theta, left, rotated_left)
               call sr_aug_rotate_matrix(generator, theta, right, rotated_right)
               left = rotated_left; right = rotated_right
            end if
            observable = matmul(left, matmul(sigma, right))
            angular = angular + coefficient(iterm)*observable/radial%rofi(ir)**2
         end if
      end do
      angular = real(l*(l + 1), rp)*angular
      total = upper + sr_aug_lower_spin_factor*(small + angular)

   end subroutine branch_observable_with_generator

   subroutine radial_spin_matrix(radial, ir, l, component, raw_power, energy_power, matrix)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, component, raw_power, energy_power
      complex(rp), intent(out) :: matrix(2, 2)
      real(rp) :: up, down, radius

      if (raw_power < 0 .or. raw_power > 2 .or. energy_power < 0 .or. energy_power > 2) then
         error stop 'DRESP-09Y: invalid radial branch power'
      end if
      radius = radial%rofi(ir)
      up = radial_value(radial, ir, l, 1, component, raw_power)
      down = radial_value(radial, ir, l, 2, component, raw_power)
      if (component == sr_aug_lower_angular) then
         up = up/(radial%tmc(ir, l + 1, 1)*radius)
         down = down/(radial%tmc(ir, l + 1, 2)*radius)
      end if
      up = up*radial%enu_work(l + 1, 1)**energy_power
      down = down*radial%enu_work(l + 1, 2)**energy_power
      call sr_aug_diagonal_matrix(up, down, matrix)
   end subroutine radial_spin_matrix

   pure real(rp) function radial_value(radial, ir, l, spin, component, power) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, spin, component, power
      select case (component)
      case (sr_aug_upper, sr_aug_lower_angular)
         select case (power)
         case (0); value = radial%phi_large(ir, l + 1, spin)
         case (1); value = radial%phidot_large(ir, l + 1, spin)
         case (2); value = radial%phiddot_large(ir, l + 1, spin)
         end select
      case (sr_aug_lower_small)
         select case (power)
         case (0); value = radial%phi_small(ir, l + 1, spin)
         case (1); value = radial%phidot_small(ir, l + 1, spin)
         case (2); value = radial%phiddot_small(ir, l + 1, spin)
         end select
      case default
         error stop 'DRESP-09Y: invalid radial component'
      end select
   end function radial_value

   subroutine branch_terms(branch, nterm, left_power, right_power, left_energy_power, right_energy_power, coefficient)
      integer, intent(in) :: branch
      integer, intent(out) :: nterm
      integer, intent(out) :: left_power(6), right_power(6), left_energy_power(6), right_energy_power(6)
      real(rp), intent(out) :: coefficient(6)

      nterm = 0
      left_power = 0; right_power = 0; left_energy_power = 0; right_energy_power = 0; coefficient = 0.0_rp
      select case (branch)
      case (1) ! 00
         nterm = 6
         left_power(1:6) = [0, 1, 0, 1, 2, 0]
         right_power(1:6) = [0, 0, 1, 1, 0, 2]
         left_energy_power(1:6) = [0, 1, 0, 1, 2, 0]
         right_energy_power(1:6) = [0, 0, 1, 1, 0, 2]
         coefficient(1:6) = [1.0_rp, -1.0_rp, -1.0_rp, 1.0_rp, 0.5_rp, 0.5_rp]
      case (2) ! 10
         nterm = 3
         left_power(1:3) = [1, 1, 2]
         right_power(1:3) = [0, 1, 0]
         left_energy_power(1:3) = [0, 0, 1]
         right_energy_power(1:3) = [0, 1, 0]
         coefficient(1:3) = [1.0_rp, -1.0_rp, -1.0_rp]
      case (3) ! 01
         nterm = 3
         left_power(1:3) = [0, 1, 0]
         right_power(1:3) = [1, 1, 2]
         left_energy_power(1:3) = [0, 1, 0]
         right_energy_power(1:3) = [0, 0, 1]
         coefficient(1:3) = [1.0_rp, -1.0_rp, -1.0_rp]
      case (4) ! 11
         nterm = 1
         left_power(1) = 1; right_power(1) = 1
         coefficient(1) = 1.0_rp
      case (5) ! 20
         nterm = 1
         left_power(1) = 2; right_power(1) = 0
         coefficient(1) = 0.5_rp
      case (6) ! 02
         nterm = 1
         left_power(1) = 0; right_power(1) = 2
         coefficient(1) = 0.5_rp
      case default
         error stop 'DRESP-09Y: invalid branch term request'
      end select
   end subroutine branch_terms

   subroutine rotation_from_generator(generator, theta, unitary)
      complex(rp), intent(in) :: generator(2, 2)
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: unitary(2, 2)
      integer :: i

      unitary = -cmplx(0.0_rp, 2.0_rp*sin(theta/2.0_rp), rp)*generator
      do i = 1, 2
         unitary(i, i) = unitary(i, i) + cmplx(cos(theta/2.0_rp), 0.0_rp, rp)
      end do
   end subroutine rotation_from_generator

end module lr_sr_augmentation_tangent_mod
