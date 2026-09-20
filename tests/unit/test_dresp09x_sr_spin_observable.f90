! DRESP-09X focused observable/algebra oracle.
program test_dresp09x_sr_spin_observable

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_dresp09s_scalar_relativistic_mod, only: sr_second_order_branch_point
   use lr_radial_observable_provenance_mod, only: pauli_density_from_moments, sr_density_from_moments, &
      sr_spin_density_from_moments, pauli_density_from_second_order_endpoints, &
      sr_density_from_second_order_endpoints, sr_spin_density_from_second_order_endpoints
   use lr_lmto_density_moment_tangent_mod, only: endpoint_tangent_branches_second_order, relative_matrix_residual
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator
   implicit none

   type(lmto_radial_basis) :: radial(1)
   complex(rp) :: sx(2,2), sy(2,2), sz(2,2), nmat(2,2), product(2,2), average(2,2)
   complex(rp) :: h(4,4), rho(4,4), dh(4,4), drho(4,4), generator(4,4), tangent(4,4,6)
   complex(rp) :: endpoint(4,4,6), endpoint_plus(4,4,6), endpoint_minus(4,4,6), rotation(4,4), temp(4,4)
   complex(rp) :: endpoint2(2,2,6)
   complex(rp) :: moments(2,2,3)
   real(rp) :: channel_moments(3,1,1,2), p3(1,7), sr2(1,7), srspin2(1,7)
   real(rp) :: p3_six(1,7), sr2_six(1,7), srspin2_six(1,7)
   real(rp) :: d_en, d_em, left, right, expected, polynomial, b20, b02, ratio1, ratio2
   real(rp) :: tol, error, theta, mx_plus, mx_minus, mz, up, down
   complex(rp) :: rho_phys(2,2), rho_rot(2,2), pauli_x(2,2)
   integer :: i, ir, branch, npoint

   tol = 3.0e-12_rp
   npoint = 7
   call radial(1)%initialize(npoint, 0, 2)
   radial(1)%rofi = [0.0_rp, 0.08_rp, 0.16_rp, 0.24_rp, 0.32_rp, 0.40_rp, 0.48_rp]
   radial(1)%enu_work(:,1) = 0.11_rp
   radial(1)%enu_work(:,2) = -0.07_rp
   radial(1)%enu_radial = radial(1)%enu_work
   radial(1)%channel_present = .true.
   do ir = 1, npoint
      radial(1)%tmc(ir,:,1) = 274.0_rp
      radial(1)%tmc(ir,:,2) = 275.0_rp
      radial(1)%phi_large(ir,1,1) = 0.82_rp + 0.01_rp*ir
      radial(1)%phi_large(ir,1,2) = 0.57_rp + 0.015_rp*ir
      radial(1)%phidot_large(ir,1,1) = 0.13_rp - 0.002_rp*ir
      radial(1)%phidot_large(ir,1,2) = 0.09_rp + 0.001_rp*ir
      radial(1)%phiddot_large(ir,1,1) = 0.031_rp + 0.001_rp*ir
      radial(1)%phiddot_large(ir,1,2) = -0.027_rp + 0.002_rp*ir
      radial(1)%phi_small(ir,1,:) = [0.17_rp + 0.003_rp*ir, 0.09_rp + 0.002_rp*ir]
      radial(1)%phidot_small(ir,1,:) = [0.021_rp, -0.013_rp]
      radial(1)%phiddot_small(ir,1,:) = [0.008_rp, 0.006_rp]
   end do

   ! Independent Pauli-matrix/angular-average oracle for all three axes.
   sx = cmplx(0.0_rp,0.0_rp,rp); sx(1,2)=1.0_rp; sx(2,1)=1.0_rp
   sy = cmplx(0.0_rp,0.0_rp,rp); sy(1,2)=cmplx(0.0_rp,-1.0_rp,rp); sy(2,1)=cmplx(0.0_rp,1.0_rp,rp)
   sz = cmplx(0.0_rp,0.0_rp,rp); sz(1,1)=1.0_rp; sz(2,2)=-1.0_rp
   do i = 1, 3
      average = cmplx(0.0_rp,0.0_rp,rp)
      do branch = 1, 6
         select case (branch)
         case (1); nmat = sx
         case (2); nmat = -sx
         case (3); nmat = sy
         case (4); nmat = -sy
         case (5); nmat = sz
         case (6); nmat = -sz
         end select
         if (i == 1) then
            product = matmul(nmat,matmul(sx,nmat))
         else if (i == 2) then
            product = matmul(nmat,matmul(sy,nmat))
         else
            product = matmul(nmat,matmul(sz,nmat))
         end if
         average = average + product/6.0_rp
      end do
      if (i == 1) error = maxval(abs(average + sx/3.0_rp))
      if (i == 2) error = maxval(abs(average + sy/3.0_rp))
      if (i == 3) error = maxval(abs(average + sz/3.0_rp))
      if (error > tol) error stop 'DRESP-09X Pauli angular average failed'
   end do

   ! Explicit Taylor-product and pure 20/02 fixtures, with lower part
   ! removed so the check isolates the shifted endpoint polynomial.
   ir = 3; d_en = 0.37_rp; d_em = -0.22_rp
   left = radial(1)%phi_large(ir,1,1) + (d_en-radial(1)%enu_work(1,1))*radial(1)%phidot_large(ir,1,1) + &
      0.5_rp*(d_en-radial(1)%enu_work(1,1))**2*radial(1)%phiddot_large(ir,1,1)
   right = radial(1)%phi_large(ir,1,1) + (d_em-radial(1)%enu_work(1,1))*radial(1)%phidot_large(ir,1,1) + &
      0.5_rp*(d_em-radial(1)%enu_work(1,1))**2*radial(1)%phiddot_large(ir,1,1)
   expected = (radial(1)%phi_large(ir,1,1)**2 + (d_en-radial(1)%enu_work(1,1))* &
      radial(1)%phidot_large(ir,1,1)*radial(1)%phi_large(ir,1,1) + &
      (d_em-radial(1)%enu_work(1,1))*radial(1)%phi_large(ir,1,1)*radial(1)%phidot_large(ir,1,1) + &
      (d_en-radial(1)%enu_work(1,1))*(d_em-radial(1)%enu_work(1,1))*radial(1)%phidot_large(ir,1,1)**2 + &
      0.5_rp*(d_en-radial(1)%enu_work(1,1))**2*radial(1)%phiddot_large(ir,1,1)*radial(1)%phi_large(ir,1,1) + &
      0.5_rp*(d_em-radial(1)%enu_work(1,1))**2*radial(1)%phi_large(ir,1,1)*radial(1)%phiddot_large(ir,1,1)) / radial(1)%rofi(ir)**2
   polynomial = 0.0_rp
   polynomial = polynomial + sr_second_order_branch_point(radial(1),ir,0,1,1,1,0.0_rp)
   polynomial = polynomial + d_en*sr_second_order_branch_point(radial(1),ir,0,1,1,2,0.0_rp)
   polynomial = polynomial + d_em*sr_second_order_branch_point(radial(1),ir,0,1,1,3,0.0_rp)
   polynomial = polynomial + d_en*d_em*sr_second_order_branch_point(radial(1),ir,0,1,1,4,0.0_rp)
   polynomial = polynomial + d_en*d_en*sr_second_order_branch_point(radial(1),ir,0,1,1,5,0.0_rp)
   polynomial = polynomial + d_em*d_em*sr_second_order_branch_point(radial(1),ir,0,1,1,6,0.0_rp)
   if (abs(polynomial-expected) > tol) then
      error stop 'DRESP-09X six-branch Taylor product failed'
   end if
   b20 = sr_second_order_branch_point(radial(1),ir,0,1,1,5,0.0_rp)
   b02 = sr_second_order_branch_point(radial(1),ir,0,1,1,6,0.0_rp)
   if (abs(b20 - 0.5_rp*radial(1)%phiddot_large(ir,1,1)*radial(1)%phi_large(ir,1,1)/radial(1)%rofi(ir)**2) > tol .or. &
       abs(b02 - 0.5_rp*radial(1)%phi_large(ir,1,1)*radial(1)%phiddot_large(ir,1,1)/radial(1)%rofi(ir)**2) > tol) then
      error stop 'DRESP-09X pure 20/02 branch failed'
   end if

   ! Equilibrium six-branch targets must reproduce the independent production
   ! moment polynomials for Pauli, SR probability, and physical SR spin.
   channel_moments = 0.0_rp
   channel_moments(:,1,1,1) = [2.0_rp, 0.31_rp, 0.071_rp]
   channel_moments(:,1,1,2) = [1.0_rp, -0.12_rp, 0.019_rp]
   moments = cmplx(0.0_rp,0.0_rp,rp)
   moments(1,1,1) = cmplx(channel_moments(1,1,1,1),0.0_rp,rp)
   moments(2,2,1) = cmplx(channel_moments(1,1,1,2),0.0_rp,rp)
   moments(1,1,2) = cmplx(channel_moments(2,1,1,1),0.0_rp,rp)
   moments(2,2,2) = cmplx(channel_moments(2,1,1,2),0.0_rp,rp)
   moments(1,1,3) = cmplx(channel_moments(3,1,1,1),0.0_rp,rp)
   moments(2,2,3) = cmplx(channel_moments(3,1,1,2),0.0_rp,rp)
   endpoint2 = cmplx(0.0_rp,0.0_rp,rp)
   endpoint2(:,:,1)=moments(:,:,1); endpoint2(:,:,2)=moments(:,:,2); endpoint2(:,:,3)=moments(:,:,2)
   endpoint2(:,:,4)=moments(:,:,3); endpoint2(:,:,5)=moments(:,:,3); endpoint2(:,:,6)=moments(:,:,3)
   call pauli_density_from_moments(radial,channel_moments,p3)
   call sr_density_from_moments(radial,channel_moments,sr2)
   call sr_spin_density_from_moments(radial,channel_moments,srspin2)
   call pauli_density_from_second_order_endpoints(radial,endpoint2,p3_six)
   call sr_density_from_second_order_endpoints(radial,endpoint2,sr2_six)
   call sr_spin_density_from_second_order_endpoints(radial,endpoint2,srspin2_six)
   if (maxval(abs(p3(:,2:)-p3_six(:,2:))) > 1.0e-11_rp .or. maxval(abs(sr2(:,2:)-sr2_six(:,2:))) > 1.0e-11_rp .or. &
       maxval(abs(srspin2(:,2:)-srspin2_six(:,2:))) > 1.0e-11_rp) then
      write(*,'(a,3es16.8)') 'six target errors = ', maxval(abs(p3(:,2:)-p3_six(:,2:))), &
         maxval(abs(sr2(:,2:)-sr2_six(:,2:))), maxval(abs(srspin2(:,2:)-srspin2_six(:,2:)))
      error stop 'DRESP-09X equilibrium six-branch reconstruction failed'
   end if

   ! Six endpoint product-rule tangent and central finite-angle oracle.
   h = cmplx(0.0_rp,0.0_rp,rp); rho = h
   h(1,1)=0.10_rp; h(2,2)=0.21_rp; h(3,3)=-0.13_rp; h(4,4)=0.29_rp
   rho(1,1)=0.80_rp; rho(2,2)=0.55_rp; rho(3,3)=0.30_rp; rho(4,4)=0.12_rp
   call exact_ks_spin_rotation_generator(2,1,generator)
   call commutator(generator,h,dh); call commutator(generator,rho,drho)
   call endpoint_tangent_branches_second_order(h,rho,dh,drho,tangent)
   endpoint(:,:,1)=rho; endpoint(:,:,2)=matmul(h,rho); endpoint(:,:,3)=matmul(rho,h)
   endpoint(:,:,4)=matmul(h,matmul(rho,h)); endpoint(:,:,5)=matmul(h,matmul(h,rho)); endpoint(:,:,6)=matmul(rho,matmul(h,h))
   do branch = 1, 6
      call commutator(generator,endpoint(:,:,branch),temp)
      if (relative_matrix_residual(tangent(:,:,branch),temp) > tol) error stop 'DRESP-09X endpoint commutator tangent failed'
   end do
   theta=1.0e-4_rp
   call spin_rotation(generator,theta,rotation)
   do branch=1,6
      endpoint_plus(:,:,branch)=matmul(rotation,matmul(endpoint(:,:,branch),conjg(transpose(rotation))))
   end do
   call spin_rotation(generator,-theta,rotation)
   do branch=1,6
      endpoint_minus(:,:,branch)=matmul(rotation,matmul(endpoint(:,:,branch),conjg(transpose(rotation))))
   end do
   if (maxval([(relative_matrix_residual((endpoint_plus(:,:,branch)-endpoint_minus(:,:,branch))/(2.0_rp*theta),tangent(:,:,branch)),branch=1,6)]) > 1.0e-7_rp) &
      error stop 'DRESP-09X finite-angle endpoint oracle failed'

   ! Physical-space finite-angle spin oracle: d/dtheta <sigma_x>=<sigma_z>.
   pauli_x=sx; up=1.7_rp; down=-0.4_rp; mz=up-down; rho_phys=cmplx(0.0_rp,0.0_rp,rp); rho_phys(1,1)=up; rho_phys(2,2)=down
   theta=1.0e-6_rp; call spin_rotation_2(theta,rotation(1:2,1:2)); rho_rot=matmul(rotation(1:2,1:2),matmul(rho_phys,conjg(transpose(rotation(1:2,1:2))))); mx_plus=real(sum(pauli_x*transpose(rho_rot)),rp)
   call spin_rotation_2(-theta,rotation(1:2,1:2)); rho_rot=matmul(rotation(1:2,1:2),matmul(rho_phys,conjg(transpose(rotation(1:2,1:2))))); mx_minus=real(sum(pauli_x*transpose(rho_rot)),rp)
   if (abs((mx_plus-mx_minus)/(2.0_rp*theta)-mz) > 1.0e-8_rp) error stop 'DRESP-09X physical spin rotation oracle failed'

   write(*,'(a)') 'UnitDresp09XSrSpinObservable: PASS (KH x/y/z, six branches, targets, tangents, physical rotation)'

contains

   subroutine commutator(generator,matrix,delta)
      complex(rp), intent(in) :: generator(:,:), matrix(:,:)
      complex(rp), intent(out) :: delta(:,:)
      delta=cmplx(0.0_rp,-1.0_rp,rp)*(matmul(generator,matrix)-matmul(matrix,generator))
   end subroutine commutator

   subroutine spin_rotation(generator,theta,rotation)
      complex(rp), intent(in) :: generator(:,:)
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: rotation(:,:)
      integer :: i
      rotation=-cmplx(0.0_rp,2.0_rp*sin(theta/2.0_rp),rp)*generator
      do i=1,size(rotation,1); rotation(i,i)=rotation(i,i)+cmplx(cos(theta/2.0_rp),0.0_rp,rp); end do
   end subroutine spin_rotation

   subroutine spin_rotation_2(theta,rotation)
      real(rp), intent(in) :: theta
      complex(rp), intent(out) :: rotation(:,:)
      rotation=cmplx(0.0_rp,0.0_rp,rp); rotation(1,1)=cmplx(cos(theta/2.0_rp),0.0_rp,rp); rotation(2,2)=rotation(1,1)
      rotation(1,2)=cmplx(-sin(theta/2.0_rp),0.0_rp,rp); rotation(2,1)=cmplx(sin(theta/2.0_rp),0.0_rp,rp)
   end subroutine spin_rotation_2

end program test_dresp09x_sr_spin_observable
