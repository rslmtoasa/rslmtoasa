!------------------------------------------------------------------------------
! DRESP-09Z arbitrary-(L,M) augmentation-frame observable tangent.
!
! This module is the spatial companion to lr_sr_augmentation_tangent_mod.
! The laboratory harmonic is held fixed.  Only the spin-dependent radial
! augmentation matrices and their coefficient-frame representation are
! rotated.  The six endpoint branches are explicit and ordered 00/10/01/
! 11/20/02.
!------------------------------------------------------------------------------
module lr_sr_spatial_augmentation_tangent_mod

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_sr_augmentation_tangent_mod, only: sr_aug_generator_y, sr_aug_rotate_matrix, sr_aug_matrix_tangent, &
      sr_aug_lower_spin_factor
   use lr_sr_angular_vertex_mod, only: sr_angular_upper_vertex, sr_angular_rank2_integral, &
      sr_angular_nn_rank2_coefficient
   implicit none
   private

   public :: sr_spatial_branch_observable
   public :: sr_spatial_branch_observable_at_angle
   public :: sr_spatial_branch_observable_tangent
   public :: sr_spatial_six_branch_point

contains

   !> Static arbitrary spatial vertex, split into the four requested pieces.
   subroutine sr_spatial_branch_observable(radial, ir, l, m, lp, mp, branch, source_l, source_m, component, upper, &
                                           lower_small, lower_rank0, lower_rank2, total)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, m, lp, mp, branch, source_l, source_m, component
      complex(rp), intent(out) :: upper(2,2), lower_small(2,2), lower_rank0(2,2), lower_rank2(2,2), total(2,2)
      complex(rp) :: generator(2,2)

      call sr_aug_generator_y(generator)
      call spatial_branch_with_generator(radial, ir, l, m, lp, mp, branch, source_l, source_m, component, &
         generator, .false., 0.0_rp, upper, lower_small, lower_rank0, lower_rank2, total)
   end subroutine sr_spatial_branch_observable

   subroutine sr_spatial_branch_observable_at_angle(radial, ir, l, m, lp, mp, branch, source_l, source_m, component, &
                                                    theta, upper, lower_small, lower_rank0, lower_rank2, total)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, m, lp, mp, branch, source_l, source_m, component
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: upper(2,2), lower_small(2,2), lower_rank0(2,2), lower_rank2(2,2), total(2,2)
      complex(rp) :: generator(2,2)

      call sr_aug_generator_y(generator)
      call spatial_branch_with_generator(radial, ir, l, m, lp, mp, branch, source_l, source_m, component, &
         generator, .false., theta, upper, lower_small, lower_rank0, lower_rank2, total)
   end subroutine sr_spatial_branch_observable_at_angle

   !> Analytic delta-O at theta=0.  No spatial harmonic is differentiated.
   subroutine sr_spatial_branch_observable_tangent(radial, ir, l, m, lp, mp, branch, source_l, source_m, component, &
                                                   tangent_upper, tangent_small, tangent_rank0, tangent_rank2, tangent_total)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, m, lp, mp, branch, source_l, source_m, component
      complex(rp), intent(out) :: tangent_upper(2,2), tangent_small(2,2), tangent_rank0(2,2), tangent_rank2(2,2), &
         tangent_total(2,2)
      complex(rp) :: generator(2,2)

      call sr_aug_generator_y(generator)
      call spatial_branch_with_generator(radial, ir, l, m, lp, mp, branch, source_l, source_m, component, &
         generator, .true., 0.0_rp, tangent_upper, tangent_small, tangent_rank0, tangent_rank2, tangent_total)
   end subroutine sr_spatial_branch_observable_tangent

   !> Pointwise scalar branch for one spin matrix element.  This helper is
   !> useful to callers assembling a six-branch raw coefficient vertex.
   real(rp) function sr_spatial_six_branch_point(radial, ir, l, lp, spin_left, spin_right, branch, lower_factor) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, lp, spin_left, spin_right, branch
      real(rp), intent(in) :: lower_factor
      integer :: nterm, iterm
      integer :: left_power(6), right_power(6), left_energy_power(6), right_energy_power(6)
      real(rp) :: coefficient(6), left, right, left_small, right_small, left_angular, right_angular, radius

      value = 0.0_rp
      if (ir == 1 .or. radial%rofi(ir) <= tiny(1.0_rp)) return
      call branch_terms(branch, nterm, left_power, right_power, left_energy_power, right_energy_power, coefficient)
      radius = radial%rofi(ir)
      do iterm = 1, nterm
         left = radial_value(radial, ir, l, spin_left, 1, left_power(iterm))* &
            radial%enu_work(l+1,spin_left)**left_energy_power(iterm)
         right = radial_value(radial, ir, lp, spin_right, 1, right_power(iterm))* &
            radial%enu_work(lp+1,spin_right)**right_energy_power(iterm)
         left_small = radial_value(radial, ir, l, spin_left, 2, left_power(iterm))* &
            radial%enu_work(l+1,spin_left)**left_energy_power(iterm)
         right_small = radial_value(radial, ir, lp, spin_right, 2, right_power(iterm))* &
            radial%enu_work(lp+1,spin_right)**right_energy_power(iterm)
         left_angular = left/(radial%tmc(ir,l+1,spin_left)*radius)
         right_angular = right/(radial%tmc(ir,lp+1,spin_right)*radius)
         value = value + coefficient(iterm)*(left*right + lower_factor*(left_small*right_small + &
            sqrt(real(l*(l+1)*lp*(lp+1),rp))*left_angular*right_angular))/radius**2
      end do
   end function sr_spatial_six_branch_point

   subroutine spatial_branch_with_generator(radial, ir, l, m, lp, mp, branch, source_l, source_m, component, &
                                            generator, derivative, theta, upper, lower_small, lower_rank0, lower_rank2, total)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir, l, m, lp, mp, branch, source_l, source_m, component
      complex(rp), intent(in) :: generator(2,2)
      logical, intent(in) :: derivative
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: upper(2,2), lower_small(2,2), lower_rank0(2,2), lower_rank2(2,2), total(2,2)
      integer :: nterm, iterm, p, q, j, tensor_q
      integer :: left_power(6), right_power(6), left_energy_power(6), right_energy_power(6)
      real(rp) :: coefficient(6), source_upper, angular_factor
      complex(rp) :: angular_coefficient
      complex(rp) :: left(2,2), right(2,2), dleft(2,2), dright(2,2), observable(2,2), tangent(2,2)
      complex(rp) :: left_small_m(2,2), right_small_m(2,2), left_angular_m(2,2), right_angular_m(2,2)
      complex(rp) :: sigma(2,2), sigma_j(2,2), rotated(2,2)

      upper = cmplx(0.0_rp,0.0_rp,rp); lower_small = upper; lower_rank0 = upper; lower_rank2 = upper
      total = upper
      if (ir == 1 .or. radial%rofi(ir) <= tiny(1.0_rp)) return
      if (branch < 1 .or. branch > 6) error stop 'DRESP-09Z: invalid six-branch request'
      call branch_terms(branch, nterm, left_power, right_power, left_energy_power, right_energy_power, coefficient)
      call cartesian_sigma_local(component, sigma)
      source_upper = real(sr_angular_upper_vertex(l,m,lp,mp,source_l,source_m),rp)
      angular_factor = sqrt(real(l*(l+1)*lp*(lp+1),rp))
      do iterm = 1, nterm
         p = left_power(iterm); q = right_power(iterm)
         angular_factor = sqrt(real(l*(l+1)*lp*(lp+1),rp))
         call radial_matrix(radial,ir,l,1,p,left_energy_power(iterm),left)
         call radial_matrix(radial,ir,lp,1,q,right_energy_power(iterm),right)
         call radial_matrix(radial,ir,l,2,p,left_energy_power(iterm),left_small_m)
         call radial_matrix(radial,ir,lp,2,q,right_energy_power(iterm),right_small_m)
         call radial_matrix(radial,ir,l,3,p,left_energy_power(iterm),left_angular_m)
         call radial_matrix(radial,ir,lp,3,q,right_energy_power(iterm),right_angular_m)
         if (derivative) then
            call product_tangent(generator,left,right,sigma,observable,tangent)
            upper = upper + coefficient(iterm)*source_upper*tangent/radial%rofi(ir)**2
            call product_tangent(generator,left_small_m,right_small_m,sigma,observable,tangent)
            lower_small = lower_small + coefficient(iterm)*sr_aug_lower_spin_factor*source_upper*tangent/radial%rofi(ir)**2
            call product_tangent(generator,left_angular_m,right_angular_m,sigma,observable,tangent)
            lower_rank0 = lower_rank0 + coefficient(iterm)*sr_aug_lower_spin_factor*source_upper*angular_factor*tangent/&
               radial%rofi(ir)**2
            do j = 1, 3
               call cartesian_sigma_local(j,sigma_j)
               do tensor_q = -2, 2
                  angular_coefficient = 2.0_rp*sr_angular_rank2_integral(l,m,lp,mp,source_l,source_m,component,j,tensor_q)* &
                     sr_angular_nn_rank2_coefficient(component,j,tensor_q)
                  if (abs(angular_coefficient) <= tiny(1.0_rp)) cycle
                  call product_tangent(generator,left_angular_m,right_angular_m,sigma_j,observable,tangent)
                  lower_rank2 = lower_rank2 + coefficient(iterm)*angular_coefficient*tangent/radial%rofi(ir)**2
               end do
            end do
         else
            if (theta /= 0.0_rp) then
               call sr_aug_rotate_matrix(generator,theta,left,rotated); left=rotated
               call sr_aug_rotate_matrix(generator,theta,right,rotated); right=rotated
               call sr_aug_rotate_matrix(generator,theta,left_small_m,rotated); left_small_m=rotated
               call sr_aug_rotate_matrix(generator,theta,right_small_m,rotated); right_small_m=rotated
               call sr_aug_rotate_matrix(generator,theta,left_angular_m,rotated); left_angular_m=rotated
               call sr_aug_rotate_matrix(generator,theta,right_angular_m,rotated); right_angular_m=rotated
            end if
            upper = upper + coefficient(iterm)*source_upper*matmul(left,matmul(sigma,right))/radial%rofi(ir)**2
            lower_small = lower_small + coefficient(iterm)*sr_aug_lower_spin_factor*source_upper*&
               matmul(left_small_m,matmul(sigma,right_small_m))/radial%rofi(ir)**2
            lower_rank0 = lower_rank0 + coefficient(iterm)*sr_aug_lower_spin_factor*source_upper*angular_factor*&
               matmul(left_angular_m,matmul(sigma,right_angular_m))/radial%rofi(ir)**2
            do j = 1, 3
               call cartesian_sigma_local(j,sigma_j)
               do tensor_q = -2, 2
                  angular_coefficient = 2.0_rp*sr_angular_rank2_integral(l,m,lp,mp,source_l,source_m,component,j,tensor_q)* &
                     sr_angular_nn_rank2_coefficient(component,j,tensor_q)
                  if (abs(angular_coefficient) <= tiny(1.0_rp)) cycle
                  lower_rank2 = lower_rank2 + coefficient(iterm)*angular_coefficient*&
                     matmul(left_angular_m,matmul(sigma_j,right_angular_m))/radial%rofi(ir)**2
               end do
            end do
         end if
      end do
      total = upper + lower_small + lower_rank0 + lower_rank2
   end subroutine spatial_branch_with_generator

   subroutine product_tangent(generator,left,right,sigma,observable,tangent)
      complex(rp), intent(in) :: generator(2,2), left(2,2), right(2,2), sigma(2,2)
      complex(rp), intent(out) :: observable(2,2), tangent(2,2)
      complex(rp) :: dl(2,2), dr(2,2)
      call sr_aug_matrix_tangent(generator,left,dl); call sr_aug_matrix_tangent(generator,right,dr)
      observable = matmul(left,matmul(sigma,right))
      tangent = matmul(dl,matmul(sigma,right)) + matmul(left,matmul(sigma,dr))
   end subroutine product_tangent

   subroutine radial_matrix(radial,ir,l,component,power,energy_power,matrix)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir,l,component,power,energy_power
      complex(rp), intent(out) :: matrix(2,2)
      real(rp) :: up,down,radius
      radius = radial%rofi(ir)
      up = radial_value(radial,ir,l,1,component,power)*radial%enu_work(l+1,1)**energy_power
      down = radial_value(radial,ir,l,2,component,power)*radial%enu_work(l+1,2)**energy_power
      if (component == 3) then
         up = up/(radial%tmc(ir,l+1,1)*radius); down = down/(radial%tmc(ir,l+1,2)*radius)
      end if
      matrix = cmplx(0.0_rp,0.0_rp,rp); matrix(1,1)=up; matrix(2,2)=down
   end subroutine radial_matrix

   pure real(rp) function radial_value(radial,ir,l,spin,component,power) result(value)
      type(lmto_radial_basis), intent(in) :: radial
      integer, intent(in) :: ir,l,spin,component,power
      select case(component)
      case(1,3)
         select case(power); case(0); value=radial%phi_large(ir,l+1,spin); case(1); value=radial%phidot_large(ir,l+1,spin); case(2); value=radial%phiddot_large(ir,l+1,spin); end select
      case(2)
         select case(power); case(0); value=radial%phi_small(ir,l+1,spin); case(1); value=radial%phidot_small(ir,l+1,spin); case(2); value=radial%phiddot_small(ir,l+1,spin); end select
      case default; error stop 'DRESP-09Z: invalid radial component'
      end select
   end function radial_value

   subroutine branch_terms(branch,nterm,left_power,right_power,left_energy_power,right_energy_power,coefficient)
      integer,intent(in)::branch
      integer,intent(out)::nterm,left_power(6),right_power(6),left_energy_power(6),right_energy_power(6)
      real(rp),intent(out)::coefficient(6)
      nterm=0; left_power=0; right_power=0; left_energy_power=0; right_energy_power=0; coefficient=0.0_rp
      select case(branch)
      case(1); nterm=6; left_power(1:6)=[0,1,0,1,2,0]; right_power(1:6)=[0,0,1,1,0,2]; left_energy_power(1:6)=left_power(1:6); right_energy_power(1:6)=right_power(1:6); coefficient(1:6)=[1.0_rp,-1.0_rp,-1.0_rp,1.0_rp,0.5_rp,0.5_rp]
      case(2); nterm=3; left_power(1:3)=[1,1,2]; right_power(1:3)=[0,1,0]; left_energy_power(1:3)=[0,0,1]; right_energy_power(1:3)=[0,1,0]; coefficient(1:3)=[1.0_rp,-1.0_rp,-1.0_rp]
      case(3); nterm=3; left_power(1:3)=[0,1,0]; right_power(1:3)=[1,1,2]; left_energy_power(1:3)=[0,1,0]; right_energy_power(1:3)=[0,0,1]; coefficient(1:3)=[1.0_rp,-1.0_rp,-1.0_rp]
      case(4); nterm=1; left_power(1)=1; right_power(1)=1; coefficient(1)=1.0_rp
      case(5); nterm=1; left_power(1)=2; right_power(1)=0; coefficient(1)=0.5_rp
      case(6); nterm=1; left_power(1)=0; right_power(1)=2; coefficient(1)=0.5_rp
      case default; error stop 'DRESP-09Z: invalid branch term request'
      end select
   end subroutine branch_terms

   pure subroutine cartesian_sigma_local(component,sigma)
      integer,intent(in)::component
      complex(rp),intent(out)::sigma(2,2)
      sigma=cmplx(0.0_rp,0.0_rp,rp)
      select case(component)
      case(1); sigma(1,2)=1.0_rp; sigma(2,1)=1.0_rp
      case(2); sigma(1,2)=cmplx(0.0_rp,-1.0_rp,rp); sigma(2,1)=cmplx(0.0_rp,1.0_rp,rp)
      case(3); sigma(1,1)=1.0_rp; sigma(2,2)=-1.0_rp
      case default; error stop 'DRESP-09Z: invalid spin component'
      end select
   end subroutine cartesian_sigma_local

end module lr_sr_spatial_augmentation_tangent_mod
