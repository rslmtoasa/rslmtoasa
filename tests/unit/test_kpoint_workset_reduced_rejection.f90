program test_kpoint_workset_reduced_rejection
   use precision_mod, only: rp
   use mpi_mod, only: parallel_context
   use kpoint_workset_mod, only: kpoint_workset, make_kpoint_workset
   use logger_mod, only: g_logger
   implicit none

   type(parallel_context) :: context
   type(kpoint_workset) :: reduced, shifted
   real(rp) :: points(3, 2), weights(2)

   call g_logger%init()
   context%rank = 0
   context%size = 1
   points(:, 1) = [0.0_rp, 0.0_rp, 0.0_rp]
   points(:, 2) = [0.25_rp, 0.0_rp, 0.0_rp]
   weights = [0.5_rp, 0.5_rp]
   reduced = make_kpoint_workset(points, weights, context, .false., .false.)
   shifted = reduced%shifted([0.2_rp, 0.0_rp, 0.0_rp])

   ! Reaching this line means the finite-q reduced-workset rejection was
   ! removed.  The CTest entry is marked WILL_FAIL and therefore accepts only
   ! the fatal contract violation above.
   write(*, '(a)') 'UNEXPECTED: reduced finite-q workset was accepted'
end program test_kpoint_workset_reduced_rejection
