program test_ldau_fixed_u_residual_mpi
   use precision_mod, only: rp
   use mix_mod, only: fixed_u_residual_metrics
   use mpi_mod, only: rank, numprocs
   use mpi
   implicit none

   complex(rp) :: pin(1,2,2,2), pout(1,2,2,2), vin(1,2,2,2), vout(1,2,2,2)
   real(rp) :: rp_abs, rp_rel, rv_abs, rv_rel, ahp, ahv
   logical :: physical
   integer :: worst_site, ierr, myrank, nprocs, contributing_rank

   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
   rank = myrank
   numprocs = nprocs
   contributing_rank = min(1, nprocs - 1)

   pin = (0.0_rp, 0.0_rp); pout = pin; vin = pin; vout = pin
   pin(1,1,1,1) = 0.3_rp; pin(1,1,2,2) = 0.2_rp
   pin(1,2,1,1) = 0.4_rp; pin(1,2,2,2) = 0.1_rp
   pout = pin
   vin(1,1,1,1) = -0.2_rp; vin(1,1,2,2) = 0.2_rp
   vin(1,2,1,1) = -0.1_rp; vin(1,2,2,2) = 0.1_rp
   vout = vin

   ! Only rank 1 contributes an imaginary off-diagonal density change and
   ! a potential change. Every rank must receive the global maxima.
   if (myrank == contributing_rank) then
      pout(1,1,1,2) = cmplx(0.0_rp, 0.05_rp, rp)
      pout(1,1,2,1) = conjg(pout(1,1,1,2))
      vout(1,2,1,1) = vout(1,2,1,1) + 0.125_rp
   end if
   call fixed_u_residual_metrics(pin,pout,vin,vout,rp_abs,rp_rel,rv_abs,rv_rel, &
                                 physical,ahp,ahv,worst_site)
   call check(abs(rp_abs-sqrt(2.0_rp)*0.05_rp) < 1.0e-12_rp, &
              'MPI global max includes rank-local imaginary P changes')
   call check(abs(rv_abs-0.125_rp) < 1.0e-12_rp, &
              'MPI global max includes rank-local V changes in Ry')
   call check(physical .and. worst_site == 1, 'valid global residual reports physical site')

   ! A nonphysical matrix on only one rank must fail the global physicality
   ! reduction, even when all other ranks have valid occupations.
   pout = pin
   if (myrank == contributing_rank) pout(1,1,1,1) = -0.1_rp
   call fixed_u_residual_metrics(pin,pout,vin,vin,rp_abs,rp_rel,rv_abs,rv_rel, &
                                 physical,ahp,ahv,worst_site)
   call check(.not. physical, 'MPI global min rejects nonphysical occupation on one rank')

   ! A one-rank run exercises the same API and confirms it does not require
   ! remote contributions to produce finite metrics.
   pout = pin
   call fixed_u_residual_metrics(pin,pout,vin,vin,rp_abs,rp_rel,rv_abs,rv_rel, &
                                 physical,ahp,ahv,worst_site)
   call check(physical .and. rp_abs == 0.0_rp .and. rv_abs == 0.0_rp, &
              'serial-equivalent residual evaluation is zero at a fixed point')

   call MPI_FINALIZE(ierr)

contains

   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      integer :: local_fail, global_fail
      local_fail = 0
      if (.not. condition) then
         local_fail = 1
         write (*, '(a,i0,2a)') 'rank ', myrank, ': FAIL: ', message
      end if
      call MPI_ALLREDUCE(local_fail, global_fail, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      if (global_fail /= 0) then
         call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      end if
   end subroutine check

end program test_ldau_fixed_u_residual_mpi
