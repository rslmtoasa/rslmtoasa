! DRESP-09Y focused augmentation-frame tangent and branch oracle.
program test_dresp09y_augmentation_tangent

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_sr_augmentation_tangent_mod, only: sr_aug_generator_y, sr_aug_sigma_x, sr_aug_sigma_y, sr_aug_sigma_z, &
      sr_aug_sigma_plus, sr_aug_sigma_minus, sr_aug_diagonal_matrix, sr_aug_rotate_matrix, sr_aug_matrix_tangent, &
      sr_aug_observable_tangent, sr_aug_branch_observable, sr_aug_branch_observable_at_angle, &
      sr_aug_branch_observable_tangent, sr_aug_lower_angular
   implicit none

   type(lmto_radial_basis) :: radial(1)
   complex(rp) :: generator(2,2), sx(2,2), sy(2,2), sz(2,2), sp(2,2), sm(2,2)
   complex(rp) :: a(2,2), b(2,2), da(2,2), db(2,2), oa(2,2), tangent(2,2), op(2,2), om(2,2)
   complex(rp) :: rotated_a(2,2), rotated_b(2,2), rotated_a_minus(2,2), rotated_b_minus(2,2), oa_plus(2,2), oa_minus(2,2)
   complex(rp) :: upper(2,2), small(2,2), angular(2,2), total(2,2)
   complex(rp) :: tupper(2,2), tsmall(2,2), tangular(2,2), ttotal(2,2)
   complex(rp) :: d(2,2), dd(2,2), dplus(2,2), dminus(2,2)
   real(rp) :: theta, error, error_half, ratio, max_error, max_ratio
   real(rp) :: up, down, r
   integer :: i, ir, l, branch, npoint

   call sr_aug_sigma_x(sx); call sr_aug_sigma_y(sy); call sr_aug_sigma_z(sz)
   call sr_aug_sigma_plus(sp); call sr_aug_sigma_minus(sm)

   ! A non-axis generator: G=(0.31 sigma_x+0.73 sigma_y+0.62 sigma_z)/2.
   generator = 0.31_rp*sx + 0.73_rp*sy + 0.62_rp*sz
   generator = generator/(2.0_rp*sqrt(0.31_rp**2 + 0.73_rp**2 + 0.62_rp**2))
   call sr_aug_diagonal_matrix(2.1_rp, -0.7_rp, a)
   call sr_aug_diagonal_matrix(-0.4_rp, 1.3_rp, b)
   theta = 2.0e-5_rp
   call sr_aug_matrix_tangent(generator, a, da)
   call sr_aug_matrix_tangent(generator, b, db)
   call sr_aug_rotate_matrix(generator, theta, a, rotated_a)
   call sr_aug_rotate_matrix(generator, theta, b, rotated_b)
   call sr_aug_rotate_matrix(generator, -theta, a, rotated_a_minus)
   call sr_aug_rotate_matrix(generator, -theta, b, rotated_b_minus)
   if (maxval(abs((rotated_a - rotated_a_minus)/(2.0_rp*theta) - da)) > 2.0e-9_rp .or. &
       maxval(abs((rotated_b - rotated_b_minus)/(2.0_rp*theta) - db)) > 2.0e-9_rp) then
      error stop 'DRESP-09Y arbitrary augmentation deltaA failed'
   end if

   ! Independent 2x2 O=A_L sigma A_R oracle for all live Cartesian/circular channels.
   max_error = 0.0_rp
   max_ratio = 0.0_rp
   do i = 1, 5
      select case (i)
      case (1); call sr_aug_observable_tangent(generator, a, b, sx, oa, tangent)
      case (2); call sr_aug_observable_tangent(generator, a, b, sy, oa, tangent)
      case (3); call sr_aug_observable_tangent(generator, a, b, sz, oa, tangent)
      case (4); call sr_aug_observable_tangent(generator, a, b, sp, oa, tangent)
      case (5); call sr_aug_observable_tangent(generator, a, b, sm, oa, tangent)
      end select
      theta = 1.0e-2_rp
      call sr_aug_rotate_matrix(generator, theta, a, rotated_a)
      call sr_aug_rotate_matrix(generator, theta, b, rotated_b)
      if (i == 1) oa_plus = matmul(rotated_a, matmul(sx, rotated_b))
      if (i == 2) oa_plus = matmul(rotated_a, matmul(sy, rotated_b))
      if (i == 3) oa_plus = matmul(rotated_a, matmul(sz, rotated_b))
      if (i == 4) oa_plus = matmul(rotated_a, matmul(sp, rotated_b))
      if (i == 5) oa_plus = matmul(rotated_a, matmul(sm, rotated_b))
      call sr_aug_rotate_matrix(generator, -theta, a, rotated_a)
      call sr_aug_rotate_matrix(generator, -theta, b, rotated_b)
      if (i == 1) oa_minus = matmul(rotated_a, matmul(sx, rotated_b))
      if (i == 2) oa_minus = matmul(rotated_a, matmul(sy, rotated_b))
      if (i == 3) oa_minus = matmul(rotated_a, matmul(sz, rotated_b))
      if (i == 4) oa_minus = matmul(rotated_a, matmul(sp, rotated_b))
      if (i == 5) oa_minus = matmul(rotated_a, matmul(sm, rotated_b))
      error = maxval(abs((oa_plus - oa_minus)/(2.0_rp*theta) - tangent))
      theta = theta/2.0_rp
      call sr_aug_rotate_matrix(generator, theta, a, rotated_a)
      call sr_aug_rotate_matrix(generator, theta, b, rotated_b)
      if (i == 1) op = matmul(rotated_a, matmul(sx, rotated_b))
      if (i == 2) op = matmul(rotated_a, matmul(sy, rotated_b))
      if (i == 3) op = matmul(rotated_a, matmul(sz, rotated_b))
      if (i == 4) op = matmul(rotated_a, matmul(sp, rotated_b))
      if (i == 5) op = matmul(rotated_a, matmul(sm, rotated_b))
      call sr_aug_rotate_matrix(generator, -theta, a, rotated_a)
      call sr_aug_rotate_matrix(generator, -theta, b, rotated_b)
      if (i == 1) om = matmul(rotated_a, matmul(sx, rotated_b))
      if (i == 2) om = matmul(rotated_a, matmul(sy, rotated_b))
      if (i == 3) om = matmul(rotated_a, matmul(sz, rotated_b))
      if (i == 4) om = matmul(rotated_a, matmul(sp, rotated_b))
      if (i == 5) om = matmul(rotated_a, matmul(sm, rotated_b))
      error_half = maxval(abs((op - om)/(2.0_rp*theta) - tangent))
      ratio = error/max(error_half, tiny(1.0_rp))
      max_error = max(max_error, error)
      max_ratio = max(max_ratio, ratio)
      if (ratio < 3.0_rp .or. ratio > 5.0_rp) error stop 'DRESP-09Y 2x2 theta^2 convergence failed'
   end do
   if (max_error > 2.0e-4_rp) error stop 'DRESP-09Y 2x2 analytic observable tangent failed'

   ! Spin-independent anti-circular control.
   call sr_aug_diagonal_matrix(1.7_rp, 1.7_rp, a)
   call sr_aug_diagonal_matrix(-0.2_rp, -0.2_rp, b)
   call sr_aug_matrix_tangent(generator, a, da)
   call sr_aug_matrix_tangent(generator, b, db)
   call sr_aug_observable_tangent(generator, a, b, sx, oa, tangent)
   if (maxval(abs(da)) > 1.0e-14_rp .or. maxval(abs(db)) > 1.0e-14_rp .or. maxval(abs(tangent)) > 1.0e-14_rp) then
      error stop 'DRESP-09Y spin-independent deltaO control failed'
   end if

   ! Radial branch fixture, including both TMC angular amplitudes.
   npoint = 7
   call radial(1)%initialize(npoint, 1, 2)
   radial(1)%rofi = [0.0_rp, 0.08_rp, 0.16_rp, 0.24_rp, 0.32_rp, 0.40_rp, 0.48_rp]
   radial(1)%enu_work(:,1) = [0.11_rp, -0.03_rp]
   radial(1)%enu_work(:,2) = [-0.07_rp, 0.05_rp]
   do ir = 1, npoint
      r = radial(1)%rofi(ir)
      radial(1)%tmc(ir,:,1) = 274.0_rp + 0.2_rp*r
      radial(1)%tmc(ir,:,2) = 275.0_rp + 0.3_rp*r
      do l = 0, 1
         radial(1)%phi_large(ir,l+1,:) = [0.82_rp + 0.02_rp*l + 0.01_rp*ir, 0.57_rp + 0.01_rp*l + 0.015_rp*ir]
         radial(1)%phidot_large(ir,l+1,:) = [0.13_rp - 0.002_rp*ir, 0.09_rp + 0.001_rp*ir]
         radial(1)%phiddot_large(ir,l+1,:) = [0.031_rp + 0.001_rp*ir, -0.027_rp + 0.002_rp*ir]
         radial(1)%phi_small(ir,l+1,:) = [0.17_rp + 0.003_rp*ir, 0.09_rp + 0.002_rp*ir]
         radial(1)%phidot_small(ir,l+1,:) = [0.021_rp, -0.013_rp]
         radial(1)%phiddot_small(ir,l+1,:) = [0.008_rp, 0.006_rp]
      end do
   end do

   call sr_aug_sigma_x(sx); call sr_aug_sigma_plus(sp); call sr_aug_sigma_minus(sm)
   max_error = 0.0_rp
   do l = 0, 1
      do ir = 2, npoint
         do branch = 1, 6
            call sr_aug_branch_observable_tangent(radial(1), ir, l, branch, sx, tupper, tsmall, tangular, ttotal)
            theta = 2.0e-5_rp
            call sr_aug_branch_observable_at_angle(radial(1), ir, l, branch, theta, sx, upper, small, angular, total)
            call sr_aug_branch_observable_at_angle(radial(1), ir, l, branch, -theta, sx, op, om, dplus, dminus)
            error = maxval(abs((upper-op)/(2.0_rp*theta)-tupper))
            error = max(error, maxval(abs((small-om)/(2.0_rp*theta)-tsmall)))
            error = max(error, maxval(abs((angular-dplus)/(2.0_rp*theta)-tangular)))
            error = max(error, maxval(abs((total-dminus)/(2.0_rp*theta)-ttotal)))
            max_error = max(max_error, error)
            call sr_aug_branch_observable_tangent(radial(1), ir, l, branch, sp, tupper, tsmall, tangular, ttotal)
            call sr_aug_branch_observable_tangent(radial(1), ir, l, branch_swap(branch), sm, upper, small, angular, total)
            if (maxval(abs(upper - conjg(transpose(tupper)))) > 2.0e-12_rp) then
               error stop 'DRESP-09Y circular Hermitian relation failed'
            end if
            if (l == 0 .and. maxval(abs(tangular)) > 1.0e-14_rp) then
               error stop 'DRESP-09Y s-channel angular tangent is nonzero'
            end if
         end do
      end do
   end do
   if (max_error > 3.0e-8_rp) error stop 'DRESP-09Y radial branch tangent finite difference failed'

   ! Direct finite-dimensional trace oracle: d Tr[D O] = Tr[dD O]+Tr[D dO].
   call sr_aug_diagonal_matrix(1.8_rp, -0.4_rp, a)
   call sr_aug_diagonal_matrix(0.6_rp, 1.2_rp, b)
   d = reshape([cmplx(0.7_rp,0.0_rp,rp), cmplx(0.2_rp,-0.1_rp,rp), &
      cmplx(0.2_rp,0.1_rp,rp), cmplx(0.3_rp,0.0_rp,rp)], [2,2])
   call sr_aug_matrix_tangent(generator, d, dd)
   call sr_aug_observable_tangent(generator, a, b, sx, oa, tangent)
   theta = 1.0e-5_rp
   call sr_aug_rotate_matrix(generator, theta, d, dplus)
   call sr_aug_rotate_matrix(generator, -theta, d, dminus)
   call sr_aug_rotate_matrix(generator, theta, a, rotated_a)
   call sr_aug_rotate_matrix(generator, theta, b, rotated_b)
   op = matmul(rotated_a, matmul(sx, rotated_b))
   call sr_aug_rotate_matrix(generator, -theta, a, rotated_a)
   call sr_aug_rotate_matrix(generator, -theta, b, rotated_b)
   om = matmul(rotated_a, matmul(sx, rotated_b))
   error = abs(sum((dplus*transpose(op)-dminus*transpose(om))/(2.0_rp*theta)) - &
      (sum(dd*transpose(oa)) + sum(d*transpose(tangent))))
   if (error > 2.0e-9_rp) then
      write(*,'(a,es16.8)') 'DRESP-09Y trace error = ', error
      write(*,'(a,es16.8)') 'DRESP-09Y deltaD FD error = ', maxval(abs((dplus-dminus)/(2.0_rp*theta)-dd))
      write(*,'(a,es16.8)') 'DRESP-09Y deltaO FD error = ', maxval(abs((op-om)/(2.0_rp*theta)-tangent))
      write(*,'(a,2(es16.8,1x))') 'DRESP-09Y trace FD/RHS = ', &
         real(sum((dplus*transpose(op)-dminus*transpose(om))/(2.0_rp*theta)),rp), &
         real(sum(dd*transpose(oa)) + sum(d*transpose(tangent)),rp)
      error stop 'DRESP-09Y finite-dimensional trace oracle failed'
   end if

   write(*,'(a,es16.8)') 'DRESP-09Y 2x2 analytic-vs-FD max = ', max_error
   write(*,'(a,es16.8)') 'DRESP-09Y theta^2 ratio max = ', max_ratio
   write(*,'(a)') 'UnitDresp09YAugmentationTangent: PASS (deltaA, deltaO, branches, controls, trace)'

contains

   integer function branch_swap(branch) result(swapped)
      integer, intent(in) :: branch
      select case (branch)
      case (2); swapped = 3
      case (3); swapped = 2
      case (5); swapped = 6
      case (6); swapped = 5
      case default; swapped = branch
      end select
   end function branch_swap

end program test_dresp09y_augmentation_tangent
