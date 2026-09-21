! DRESP-10R Cartesian/circular transverse convention oracle.
program test_lr_full_spatial_circular_convention
   use precision_mod, only: rp
   use lr_full_spatial_alsda_mod, only: lr_full_spatial_cartesian_to_circular, &
      lr_full_spatial_circular_to_cartesian
   implicit none

   real(rp), parameter :: tolerance = 2.0e-13_rp
   complex(rp) :: sx(2,2), sy(2,2), plus(2,2), minus(2,2), x_back(2,2), y_back(2,2)
   complex(rp) :: generator(2,2), sz(2,2), native_x(2,2), native_x_back(2,2), native_y_back(2,2)
   real(rp) :: error

   sx = cmplx(0.0_rp,0.0_rp,rp); sx(1,2)=1.0_rp; sx(2,1)=1.0_rp
   sy = cmplx(0.0_rp,0.0_rp,rp); sy(1,2)=cmplx(0.0_rp,-1.0_rp,rp); sy(2,1)=cmplx(0.0_rp,1.0_rp,rp)
   sz = cmplx(0.0_rp,0.0_rp,rp); sz(1,1)=1.0_rp; sz(2,2)=-1.0_rp
   call lr_full_spatial_cartesian_to_circular(sx, sy, plus, minus)
   call lr_full_spatial_circular_to_cartesian(plus, minus, x_back, y_back)
   error = maxval(abs(x_back-sx)) + maxval(abs(y_back-sy))

   generator = 0.5_rp*sy
   native_x = cmplx(0.0_rp,-1.0_rp,rp)*(matmul(generator,sz)-matmul(sz,generator))
   call lr_full_spatial_cartesian_to_circular(native_x, sy, plus, minus)
   call lr_full_spatial_circular_to_cartesian(plus, minus, native_x_back, native_y_back)
   error = max(error, maxval(abs(native_x_back-native_x)), maxval(abs(native_y_back-sy)))
   write (*,'(a,es12.4)') '  cartesian_circular_max_error=', error
   if (error > tolerance) error stop 'UnitLrFullSpatialCircularConvention: FAIL'
   write (*,'(a)') 'UnitLrFullSpatialCircularConvention: PASS (native SU(2) x tangent; exact +/- reconstruction)'
end program test_lr_full_spatial_circular_convention
