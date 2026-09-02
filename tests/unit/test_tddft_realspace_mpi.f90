!------------------------------------------------------------------------------
! TDDFT-R2-01 -- deterministic native real-space GF MPI reduction regression.
!
! The root rank builds the serial reference from all real-space pairs.  Every
! rank then builds only its contiguous local pair block, Fourier transforms its
! local contribution at q=0 and at a generic finite q, and collectively reduces
! the complex chi0(q,omega).  Thus the test exercises the same ownership and
! reduction contract as the production native-RGF route without constructing a
! reciprocal-space Green function.
!------------------------------------------------------------------------------
program test_tddft_realspace_mpi
   use precision_mod, only: rp
   use response_components_mod, only: RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel
   use tddft_chi0_mod, only: tddft_chi0_request, tddft_chi0_batch_result
   use tddft_chi0_realspace_mod, only: tddft_realspace_chi0_options, tddft_native_realspace_gf_provider, &
      reduce_realspace_chi0_batch
   use tddft_performance_mod, only: tddft_work_range, split_tddft_work
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   integer, parameter :: nmat = 2, nsite = 2, nresponse = 2, ne = 257, nr = 8, nq = 2, nw = 2
   real(rp), parameter :: absolute_tolerance = 2.0e-11_rp
   real(rp), parameter :: relative_tolerance = 2.0e-12_rp
   integer :: mpi_rank, mpi_size, mpi_ierr, ip, ie, nlocal
   type(tddft_native_realspace_gf_provider) :: reference_provider, local_provider
   type(tddft_realspace_chi0_options) :: options
   type(tddft_chi0_request) :: request
   type(tddft_chi0_batch_result) :: reference_batch, local_batch
   type(tddft_work_range) :: local_range
   type(response_channel) :: left_channels(nresponse), right_channels(nresponse)
   real(rp) :: energy(ne), r_vectors_all(3, nr), q_points(3, nq), omega(nw)
   complex(rp) :: reference_chi(nresponse, nresponse, nw, nq)
   integer :: pair_sites_all(nr, 2)
   integer, allocatable :: local_pair_sites(:, :)
   complex(rp) :: g_ab_all(nmat, nmat, ne, nr), g_ba_all(nmat, nmat, ne, nr)
   complex(rp), allocatable :: local_g_ab(:, :, :, :), local_g_ba(:, :, :, :)
   real(rp), allocatable :: local_r_vectors(:, :)
   real(rp) :: difference, reference_scale, max_reference_scale, local_partial_difference
   real(rp) :: max_partial_difference, max_absolute_difference, max_relative_difference
   logical :: failed

   mpi_rank = 0
   mpi_size = 1
#ifdef USE_MPI
   call MPI_INIT(mpi_ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpi_ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, mpi_ierr)
#endif

   call build_fixture(energy, r_vectors_all, pair_sites_all, g_ab_all, g_ba_all, q_points, omega, left_channels, &
      right_channels, options)
   allocate(request%q_points(3, nq), request%omega(nw))
   request%q_points = q_points
   request%omega = omega

   reference_chi = cmplx(0.0_rp, 0.0_rp, rp)
   if (mpi_rank == 0) then
      call reference_provider%initialize(energy, g_ab_all, g_ba_all, r_vectors_all, pair_sites_all, [1, 1], &
         left_channels, right_channels, options)
      call reference_provider%evaluate_realspace(request, reference_batch)
      do ip = 1, nq
         reference_chi(:, :, :, ip) = reference_batch%q_response(ip)%chi
      end do
   end if
#ifdef USE_MPI
   call MPI_BCAST(reference_chi, size(reference_chi), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, mpi_ierr)
#endif

   call split_tddft_work(mpi_rank, mpi_size, nr, local_range)
   nlocal = local_range%count
   allocate(local_g_ab(nmat, nmat, ne, nlocal), local_g_ba(nmat, nmat, ne, nlocal), local_r_vectors(3, nlocal), &
      local_pair_sites(nlocal, 2))
   if (nlocal > 0) then
      local_g_ab = g_ab_all(:, :, :, local_range%first:local_range%last)
      local_g_ba = g_ba_all(:, :, :, local_range%first:local_range%last)
      local_r_vectors = r_vectors_all(:, local_range%first:local_range%last)
      local_pair_sites = pair_sites_all(local_range%first:local_range%last, :)
   end if
   call local_provider%initialize(energy, local_g_ab, local_g_ba, local_r_vectors, local_pair_sites, [1, 1], &
      left_channels, right_channels, options)
   call local_provider%evaluate_realspace(request, local_batch)
   local_partial_difference = 0.0_rp
   do ip = 1, nq
      local_partial_difference = max(local_partial_difference, &
         maxval(abs(local_batch%q_response(ip)%chi-reference_chi(:, :, :, ip))))
   end do
#ifdef USE_MPI
   call MPI_ALLREDUCE(local_partial_difference, max_partial_difference, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_ierr)
#else
   max_partial_difference = local_partial_difference
#endif
   call reduce_realspace_chi0_batch(local_batch)

   max_absolute_difference = 0.0_rp
   max_relative_difference = 0.0_rp
   max_reference_scale = maxval(abs(reference_chi))
   do ip = 1, nq
      difference = maxval(abs(local_batch%q_response(ip)%chi-reference_chi(:, :, :, ip)))
      reference_scale = maxval(abs(reference_chi(:, :, :, ip)))
      max_absolute_difference = max(max_absolute_difference, difference)
      max_relative_difference = max(max_relative_difference, difference/max(1.0e-30_rp, reference_scale))
   end do
#ifdef USE_MPI
   call MPI_ALLREDUCE(MPI_IN_PLACE, max_absolute_difference, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, max_relative_difference, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_ierr)
#endif

   failed = max_absolute_difference > absolute_tolerance .or. max_relative_difference > relative_tolerance
   if (mpi_size > 1 .and. max_partial_difference <= 1.0e-8_rp) failed = .true.
   if (local_batch%q_response(1)%metadata%real_space_points /= nr) failed = .true.
   if (local_batch%q_response(1)%metadata%real_space_omitted_points /= 0) failed = .true.
   if (mpi_rank == 0) then
      write(*, '(a,i0,a,es24.16,a,es24.16,a,es24.16,a,es24.16)') 'TDDFT_R2_01 ranks=', mpi_size, ' max_abs=', &
         max_absolute_difference, ' max_rel=', max_relative_difference, ' chi_scale=', max_reference_scale, &
         ' max_partial_abs=', max_partial_difference
      if (max_absolute_difference > absolute_tolerance) write(*, '(a,es16.8)') &
         'TDDFT_R2_01 absolute tolerance exceeded: ', absolute_tolerance
      if (max_relative_difference > relative_tolerance) write(*, '(a,es16.8)') &
         'TDDFT_R2_01 relative tolerance exceeded: ', relative_tolerance
   end if
#ifdef USE_MPI
   call MPI_ALLREDUCE(MPI_IN_PLACE, failed, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, mpi_ierr)
   call MPI_FINALIZE(mpi_ierr)
#endif
   if (failed) error stop 'FAIL: native real-space GF MPI reduction changed complex chi0'

contains

   subroutine build_fixture(energy, r_vectors, pair_sites, g_ab, g_ba, q_points, omega, left, right, options)
      real(rp), intent(out) :: energy(:), r_vectors(:, :), q_points(:, :), omega(:)
      integer, intent(out) :: pair_sites(:, :)
      complex(rp), intent(out) :: g_ab(:, :, :, :), g_ba(:, :, :, :)
      type(response_channel), intent(out) :: left(:), right(:)
      type(tddft_realspace_chi0_options), intent(out) :: options
      integer :: i, j, ie, ip
      real(rp) :: pole, width
      complex(rp) :: denominator, numerator

      do ie = 1, size(energy)
         energy(ie) = -1.0_rp + 2.0_rp*real(ie-1, rp)/real(size(energy)-1, rp)
      end do
      r_vectors = 0.0_rp
      r_vectors(:, 2) = [1.0_rp, 0.0_rp, 0.0_rp]
      r_vectors(:, 3) = [-1.0_rp, 0.0_rp, 0.0_rp]
      r_vectors(:, 4) = [0.0_rp, 1.0_rp, 0.0_rp]
      r_vectors(:, 5) = [0.0_rp, -1.0_rp, 0.0_rp]
      r_vectors(:, 6) = [1.0_rp, 1.0_rp, 0.0_rp]
      r_vectors(:, 7) = [-1.0_rp, 1.0_rp, 0.0_rp]
      r_vectors(:, 8) = [2.0_rp, 0.0_rp, 0.0_rp]
      pair_sites = reshape([1, 1, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2, 2], [nr, 2])
      do ip = 1, size(pair_sites, 1)
         do ie = 1, size(energy)
            do i = 1, size(g_ab, 1)
               do j = 1, size(g_ab, 2)
                  pole = -0.35_rp + 0.025_rp*real(ip, rp) + 0.011_rp*real(i+j, rp)
                  width = 0.045_rp + 0.002_rp*real(mod(i+j+ip, 3), rp)
                  denominator = cmplx(energy(ie)-pole, width, rp)
                  numerator = cmplx(0.08_rp*real(i+ip, rp)+0.03_rp*real(j, rp), &
                     -0.01_rp*real(i-j+ip, rp), rp)
                  g_ab(i, j, ie, ip) = numerator/denominator
                  g_ba(i, j, ie, ip) = cmplx(0.91_rp, 0.07_rp, rp)*g_ab(i, j, ie, ip) + &
                     cmplx(0.002_rp*real(i-j, rp), -0.001_rp*real(ip, rp), rp)/denominator
               end do
            end do
         end do
      end do

      left(1) = response_channel(1, RESPONSE_PLUS)
      left(2) = response_channel(2, RESPONSE_PLUS)
      right(1) = response_channel(1, RESPONSE_MINUS)
      right(2) = response_channel(2, RESPONSE_MINUS)
      q_points = 0.0_rp
      q_points(:, 2) = [0.173_rp, -0.219_rp, 0.071_rp]
      omega = [0.07_rp, 0.23_rp]
      options = tddft_realspace_chi0_options()
      options%eta = 0.04_rp
      options%green_eta = 0.02_rp
      options%fermi_level = 0.0_rp
      options%electronic_temperature = 0.0_rp
      options%rmax = huge(1.0_rp)
      options%tail_tolerance = 1.0e-12_rp
   end subroutine build_fixture

end program test_tddft_realspace_mpi
