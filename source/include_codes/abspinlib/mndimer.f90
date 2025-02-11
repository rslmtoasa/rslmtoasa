program test_abSpinlib
   !
   use Constants
   use abSpinLib
   !
   implicit none
   !
   integer :: na
   integer, parameter :: ldimx = 2
   integer, parameter :: ldimy = 1
   integer :: alloc_flag
   integer :: i, j, iiter, niter, x_plus, x_min, y_plus, y_min, x, y
   real(selected_real_kind(15)) :: delta
   real(selected_real_kind(15)) :: temp = 0.0d0
   real(selected_real_kind(15)) :: dt = 1d-3
   !
   real(selected_real_kind(15)), dimension(:, :), allocatable :: moments
   real(selected_real_kind(15)), dimension(:), allocatable :: mnorms
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
   allocate (mnorms(na))
   allocate (fields(3, na))
   allocate (nnlist(1, na))
   nnlist = 0
   allocate (exchange(1, na))
   !exchange=1.2d-3/2.0d0
   ! J1=1.2mRy, factor 2 for double application of exchange field then conversion to Tesla
   exchange = -3.4d0*mry/mub

   !positions on cartesian grid
   atomgrid = 0
   !atom 1
   atomgrid(1, 1) = 1
   !atom 1
   atomgrid(2, 1) = 1
   !
   ! moments
   moments(1, 1) = 0.00d0
   moments(2, 1) = 0.0d0
   moments(3, 1) = 5.0d0
   mnorms(1) = sqrt(sum(moments(:, 1)**2))
   moments(1, 2) = 5.00d0
   moments(2, 2) = 0.0d0
   moments(3, 2) = 0.0d0
   mnorms(2) = sqrt(sum(moments(:, 2)**2))
   !
   exchange = exchange/mnorms(1)*mnorms(2)
   !
   ! neighbour list
   nnlist(1, 1) = atomgrid(2, 1)
   nnlist(1, 2) = atomgrid(1, 1)
   i = 0
   !
   niter = 5000
   do iiter = 1, niter
      ! Calculate fields for predictor step
      fields = 0.0d0

      !field on atom 1
      fields(:, 1) = fields(:, 1) - moments(:, 2)*exchange(1, 1)/2.0d0/mnorms(2)
      !field on atom 1
      fields(:, 2) = fields(:, 2) - moments(:, 1)*exchange(1, 2)/2.0d0/mnorms(1)

      ! Then do predictor step
      call asd_pred(moments, fields, temp, dt, na)

      ! Then calculate new  fields for corrector step
      fields = 0.0d0

      !field on atom 1
      fields(:, 1) = fields(:, 1) - moments(:, 2)*exchange(1, 1)/2.0d0/mnorms(2)
      !field on atom 1
      fields(:, 2) = fields(:, 2) - moments(:, 1)*exchange(1, 2)/2.0d0/mnorms(1)

      ! Then do corrector step
      call asd_corr(moments, fields, temp, dt, na)
      print *, 'Step', iiter, " finished", sum(moments(3, :))/na
      print *, 'Energy:', 2.0d0*exchange(1, 1)*sum(moments(1:3, 1)*moments(1:3, 2))/(mry/mub), exchange(1, 1)/(mry/mub)
      !do i=1,na
      !   print ´(2x,i4,3f12.6)´, i,moments(:,i)
      !end do
   end do
   !
   call allocate_asd(na, -1)
end program test_abSpinlib
