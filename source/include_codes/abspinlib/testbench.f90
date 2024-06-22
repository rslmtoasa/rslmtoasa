program test_abSpinlib
   !
   use Constants
   use abSpinLib
   !
   implicit none
   !
   integer :: na
   integer, parameter :: ldimx = 4
   integer, parameter :: ldimy = 4
   integer :: alloc_flag
   integer :: i, j, iiter, niter, x_plus, x_min, y_plus, y_min, x, y
   real(selected_real_kind(15)) :: delta
   real(selected_real_kind(15)) :: temp = 0.0d0
   real(selected_real_kind(15)) :: dt = 1d-3
   !
   real(selected_real_kind(15)), dimension(:, :), allocatable :: moments
   real(selected_real_kind(15)), dimension(:, :), allocatable :: fields
   real(selected_real_kind(15)), dimension(:, :), allocatable :: exchange
   !  real(selected_real_kind(15)), dimension(:,:), allocatable :: exchange
   integer, dimension(:, :), allocatable :: nnlist
   integer, dimension(:, :), allocatable :: atomgrid
   !
   na = ldimx*ldimy
   !
   call allocate_asd(na, 1)
   !
   call read_asd_input()
   temp = asd_temp
   dt = asd_dt
   !
   allocate (atomgrid(1:ldimx, 1:ldimy))
   allocate (moments(3, na))
   allocate (fields(3, na))
   allocate (nnlist(4, na))
   nnlist = 0
   allocate (exchange(4, na))
   !exchange=1.2d-3/2.0d0
   ! J1=1.2mRy, factor 2 for double application of exchange field then conversion to Tesla
   exchange = 2.0d0*2.0d0*1.2d0*mry/mub
   i = 0
   do x = 1, ldimx
      do y = 1, ldimy
         i = i + 1
         atomgrid(x, y) = i
      end do
   end do
   moments(:, :) = 0.0d0
   moments(3, :) = 1.0d0
   moments(1, 8) = 1.0d0
   moments(3, 8) = 0.0d0
   i = 0
   do x = 1, ldimx
      do y = 1, ldimy
         i = i + 1
         ! Fill moments
         !call random_number(delta)
         !moments(1,i)=0.1*delta
         !call random_number(delta)
         !moments(2,i)=0.1*delta
         !moments(3,i)=1.0
         print '(a,3i4,a,4f10.4)', 'Atom: ', i, x, y, ' Moment    : ', moments(:, i)
         moments(:, i) = moments(:, i)/norm2(moments(:, i))
         ! Find neighbours
         x_plus = x + 1
         if (x_plus > ldimx) x_plus = x_plus - ldimx
         x_min = x - 1
         if (x_min < 1) x_min = x_min + ldimx
         y_plus = y + 1
         if (y_plus > ldimy) y_plus = y_plus - ldimy
         y_min = y - 1
         if (y_min < 1) y_min = y_min + ldimy
         !print ´(a,3i4,a,4i4)´,´Atom: ´,i,x,y,´ Neighbours: ´,x_plus,x_min,y_plus,y_min
         nnlist(1, i) = atomgrid(x_plus, y)
         nnlist(2, i) = atomgrid(x_min, y)
         nnlist(3, i) = atomgrid(x, y_plus)
         nnlist(4, i) = atomgrid(x, y_min)
         print '(a,3i4,a,4i4)', 'Atom: ', i, x, y, ' Neighbours: ', nnlist(1:4, i)
      end do
   end do
   !
   niter = 1000
   do iiter = 1, niter
      ! Calculate fields for predictor step
      fields = 0.0d0
      do i = 1, na
         do j = 1, 4
            fields(:, i) = fields(:, i) - moments(:, nnlist(j, i))*exchange(j, i)
         end do
      end do
      ! Then do predictor step
      call asd_pred(moments, fields, temp, dt, na)
      ! Then calculate fields for corrector step
      fields = 0.0d0
      do i = 1, na
         do j = 1, 4
            fields(:, i) = fields(:, i) - moments(:, nnlist(j, i))*exchange(j, i)
         end do
      end do
      ! Then do corrector step
      call asd_corr(moments, fields, temp, dt, na)
      print *, 'Step', iiter, " finished", sum(moments(3, :))/na
      !do i=1,na
      !   print ´(2x,i4,3f12.6)´, i,moments(:,i)
      !end do
   end do
   !
   call allocate_asd(na, -1)
end program test_abSpinlib
