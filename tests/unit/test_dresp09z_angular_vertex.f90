program test_dresp09z_angular_vertex

   use precision_mod, only: rp
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_sr_angular_vertex_mod, only: sr_angular_upper_vertex, sr_angular_lower_vertex, &
      sr_angular_quadrature_upper, sr_angular_quadrature_lower, sr_angular_circular_vertex, &
      sr_angular_circular_conjugate_residual, sr_angular_product_rank_amplitudes, sr_angular_l0_scalar_vertex
   use lr_sr_augmentation_tangent_mod, only: sr_aug_sigma_x, sr_aug_sigma_y, sr_aug_sigma_z
   use lr_sr_spatial_augmentation_tangent_mod, only: sr_spatial_branch_observable, &
      sr_spatial_branch_observable_at_angle, sr_spatial_branch_observable_tangent
   use lr_dresp09z_bridge_mod, only: dresp09z_historical_candidate_count, &
      dresp09z_shadow_six_branch_candidate_count, dresp09z_source_product_ranks
   implicit none

   integer :: l, m, lp, mp, source_l, source_m, component
   real(rp) :: upper_error, lower_error, hermitian_error, fixture_error, amplitudes(0:4)
   real(rp) :: forbidden_error, l0_reduction_error
   complex(rp) :: analytic, oracle, lower(2,2), lower_oracle(2,2), plus(2,2), minus(2,2), lower_scalar(2,2), sigma(2,2)
   complex(rp) :: weights(3)
   type(lmto_radial_basis) :: radial
   complex(rp) :: u0(2,2), s0(2,2), r0(2,2), k0(2,2), t0(2,2)
   complex(rp) :: up(2,2), sp(2,2), rp0(2,2), kp(2,2), tplus(2,2), atot(2,2)
   complex(rp) :: um(2,2), sm(2,2), rm(2,2), km(2,2), tm(2,2)
   real(rp) :: tangent_error, theta, radius
   real(rp) :: branch_error(6), augmentation_l_error(0:4)
   logical :: ranks(0:4)
   integer :: ir, branch, fixture, fixture_l(0:4), fixture_m(0:4), fixture_mp(0:4), fixture_source_m(0:4)

   upper_error = 0.0_rp; lower_error = 0.0_rp; hermitian_error = 0.0_rp; fixture_error = 0.0_rp
   forbidden_error = 0.0_rp; l0_reduction_error = 0.0_rp
   do l = 0, 2
      do m = -l, l
         do lp = 0, 2
            do mp = -lp, lp
               do source_l = 0, 4
                  do source_m = -source_l, source_l
                     analytic = sr_angular_upper_vertex(l,m,lp,mp,source_l,source_m)
                     oracle = sr_angular_quadrature_upper(l,m,lp,mp,source_l,source_m)
                     upper_error = max(upper_error,abs(analytic-oracle))
                     if (.not. allowed_upper_pair(l,m,lp,mp,source_l,source_m)) then
                        forbidden_error = max(forbidden_error,abs(analytic))
                     end if
                     do component = 1, 3
                        call sr_angular_lower_vertex(l,m,lp,mp,source_l,source_m,component,lower)
                        call sr_angular_quadrature_lower(l,m,lp,mp,source_l,source_m,component,lower_oracle)
                        lower_error = max(lower_error,maxval(abs(lower-lower_oracle)))
                        if (source_l == 0 .and. source_m == 0) then
                           select case (component)
                           case (1); call sr_aug_sigma_x(sigma)
                           case (2); call sr_aug_sigma_y(sigma)
                           case (3); call sr_aug_sigma_z(sigma)
                           end select
                           call sr_angular_l0_scalar_vertex(l,m,lp,mp,sigma,lower_scalar)
                           l0_reduction_error = max(l0_reduction_error,maxval(abs(lower-lower_scalar)))
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do

   ! The circular convention is checked with orbital bra/ket exchange and
   ! the live harmonic M -> -M covariance.
   do l = 0, 2
      do m = -l, l
         do lp = 0, 2
            do mp = -lp, lp
               do source_l = 0, 4
                  do source_m = -source_l, source_l
                     call sr_angular_circular_vertex(l,m,lp,mp,source_l,source_m,1,plus)
                     call sr_angular_circular_vertex(lp,mp,l,m,source_l,-source_m,-1,minus)
                     hermitian_error = max(hermitian_error, sr_angular_circular_conjugate_residual(plus, minus, source_l, source_m))
                  end do
               end do
            end do
         end do
      end do
   end do

   weights = cmplx(0.0_rp,0.0_rp,rp); weights(1)=0.5_rp; weights(2)=cmplx(0.0_rp,0.5_rp,rp)
   call sr_angular_product_rank_amplitudes(2,-2,2,2,4,4,weights,amplitudes)
   fixture_error = maxval(amplitudes)
   if (dresp09z_historical_candidate_count(2,4) /= 232 .or. &
       dresp09z_shadow_six_branch_candidate_count(2,4) /= 348) then
      error stop 'DRESP-09Z compact candidate inventory failed'
   end if
   call dresp09z_source_product_ranks(4,4,ranks)
   if (.not. ranks(2) .or. .not. ranks(4) .or. ranks(0) .or. ranks(1) .or. ranks(3)) then
      error stop 'DRESP-09Z source/product rank map failed'
   end if

   ! Six-branch arbitrary-L augmentation tangent and finite-angle oracle.
   call radial%initialize(7,2,2)
   radial%rofi = [0.0_rp,0.08_rp,0.16_rp,0.24_rp,0.32_rp,0.40_rp,0.48_rp]
   radial%enu_work(:,1) = [0.11_rp,-0.03_rp,0.07_rp]
   radial%enu_work(:,2) = [-0.07_rp,0.05_rp,-0.02_rp]
   do ir = 1, radial%npoint
      radius = radial%rofi(ir)
      radial%tmc(ir,:,:) = 274.0_rp + 0.2_rp*radius
      do l = 0, 2
         radial%phi_large(ir,l+1,:) = [0.82_rp+0.02_rp*l+0.01_rp*ir,0.57_rp+0.01_rp*l+0.015_rp*ir]
         radial%phidot_large(ir,l+1,:) = [0.13_rp-0.002_rp*ir,0.09_rp+0.001_rp*ir]
         radial%phiddot_large(ir,l+1,:) = [0.031_rp+0.001_rp*ir,-0.027_rp+0.002_rp*ir]
         radial%phi_small(ir,l+1,:) = [0.17_rp+0.003_rp*ir,0.09_rp+0.002_rp*ir]
         radial%phidot_small(ir,l+1,:) = [0.021_rp,-0.013_rp]
         radial%phiddot_small(ir,l+1,:) = [0.008_rp,0.006_rp]
      end do
   end do
   tangent_error = 0.0_rp
   branch_error = 0.0_rp; augmentation_l_error = 0.0_rp
   theta = 2.0e-5_rp
   fixture_l = [0,1,2,3,4]
   fixture_m = [0,0,-1,2,-2]
   fixture_mp = [0,0,0,0,2]
   fixture_source_m = [0,0,1,-2,4]
   do ir = 2, radial%npoint
      do fixture = 0, 4
         do branch = 1, 6
            call sr_spatial_branch_observable_tangent(radial,ir,2,fixture_m(fixture),2,fixture_mp(fixture),branch, &
               fixture_l(fixture),fixture_source_m(fixture),1,t0,s0,r0,k0,atot)
            call sr_spatial_branch_observable_at_angle(radial,ir,2,fixture_m(fixture),2,fixture_mp(fixture),branch, &
               fixture_l(fixture),fixture_source_m(fixture),1,theta,up,sp,rp0,kp,tplus)
            call sr_spatial_branch_observable_at_angle(radial,ir,2,fixture_m(fixture),2,fixture_mp(fixture),branch, &
               fixture_l(fixture),fixture_source_m(fixture),1,-theta,um,sm,rm,km,tm)
            branch_error(branch) = max(branch_error(branch),maxval(abs((up-um)/(2.0_rp*theta)-t0)))
            branch_error(branch) = max(branch_error(branch),maxval(abs((sp-sm)/(2.0_rp*theta)-s0)))
            branch_error(branch) = max(branch_error(branch),maxval(abs((rp0-rm)/(2.0_rp*theta)-r0)))
            branch_error(branch) = max(branch_error(branch),maxval(abs((kp-km)/(2.0_rp*theta)-k0)))
            branch_error(branch) = max(branch_error(branch),maxval(abs((tplus-tm)/(2.0_rp*theta)-atot)))
            augmentation_l_error(fixture_l(fixture)) = max(augmentation_l_error(fixture_l(fixture)),branch_error(branch))
         end do
      end do
   end do
   tangent_error = maxval(branch_error)

   write (*,'(a,es16.8)') 'DRESP-09Z upper quadrature max error = ', upper_error
   write (*,'(a,es16.8)') 'DRESP-09Z lower KH quadrature max error = ', lower_error
   write (*,'(a,es16.8)') 'DRESP-09Z forbidden upper-channel amplitude = ', forbidden_error
   write (*,'(a,es16.8)') 'DRESP-09Z raw L0-vs-certified reduction residual = ', l0_reduction_error
   write (*,'(a,es16.8)') 'DRESP-09Z circular covariance max error = ', hermitian_error
   write (*,'(a,es16.8)') 'DRESP-09Z L4 rank-mixing fixture norm = ', fixture_error
   write (*,'(a,es16.8)') 'DRESP-09Z arbitrary-L augmentation FD max error = ', tangent_error
   write (*,'(a,6(es12.4,1x))') 'DRESP-09Z branch FD errors 00/10/01/11/20/02 = ', branch_error
   write (*,'(a,5(es12.4,1x))') 'DRESP-09Z augmentation FD errors L0..L4 = ', augmentation_l_error
   if (upper_error > 2.0e-11_rp) error stop 'DRESP-09Z upper live-harmonic oracle failed'
   if (lower_error > 2.0e-11_rp) error stop 'DRESP-09Z lower KH quadrature oracle failed'
   if (hermitian_error > 2.0e-11_rp) error stop 'DRESP-09Z circular covariance failed'
   if (fixture_error <= 1.0e-12_rp) error stop 'DRESP-09Z rank-2 fixture did not mix'
   if (tangent_error > 3.0e-7_rp) error stop 'DRESP-09Z arbitrary-L augmentation tangent failed'
   write (*,'(a)') 'UnitDresp09ZAngularVertex: PASS'

contains

   pure logical function allowed_upper_pair(l,m,lp,mp,source_l,source_m) result(value)
      integer, intent(in) :: l,m,lp,mp,source_l,source_m
      value = abs(m) <= l .and. abs(mp) <= lp .and. abs(source_m) <= source_l .and. &
         source_m == mp-m .and. abs(l-lp) <= source_l .and. source_l <= l+lp .and. &
         mod(l+lp+source_l,2) == 0
   end function allowed_upper_pair

end program test_dresp09z_angular_vertex
