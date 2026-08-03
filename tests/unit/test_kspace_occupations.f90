!------------------------------------------------------------------------------
! RS-LMTO-ASA -- WP2 canonical k-space occupation/EBAND tests
!------------------------------------------------------------------------------
program test_kspace_occupations
   use precision_mod, only: rp
   use reciprocal_mod, only: reciprocal
   use logger_mod, only: g_logger
   use mpi_mod, only: rank, numprocs, ierr, get_mpi_range
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   integer, parameter :: nbands = 2, nk_global = 4
   real(rp), parameter :: target_electrons = 1.15_rp
   real(rp), parameter :: temperature = 500.0_rp
   real(rp), parameter :: tol = 2.0e-10_rp
   real(rp), parameter :: kB_Ry_per_K = 6.3336814e-6_rp
   real(rp), parameter :: eigenvalues_full(nbands, nk_global) = reshape([ &
      -1.00_rp, 0.15_rp,  -0.60_rp, 0.45_rp, &
      -0.25_rp, 0.80_rp,   0.05_rp, 1.10_rp], [nbands, nk_global])
   real(rp), parameter :: weights_full(nk_global) = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]

   type(reciprocal) :: obj
   integer :: k_start, k_end, nk_local, ik, ie, ne
   integer, allocatable :: l2g(:), g2l(:)
   real(rp) :: ef, ef_ref, nelect, eband, eband_ref, raw_weight_sum
   real(rp) :: nelect_scaled, eband_scaled, eband_grid, grid_error
   real(rp) :: nelect_again, eband_again
   real(rp) :: e_min, e_max, de, energy, sigma
   logical :: failed

#ifdef USE_MPI
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
#else
   rank = 0
   numprocs = 1
   ierr = 0
#endif
   call g_logger%init()
   failed = .false.

   call get_mpi_range(rank, nk_global, k_start, k_end, nk_local, l2g, g2l, 'WP2 occupations')
   allocate(obj%eigenvalues(nbands, nk_local))
   allocate(obj%k_weights(nk_global))
   allocate(obj%k_l2g_map(nk_local))
   allocate(obj%k_g2l_map(nk_global))
   obj%k_weights = weights_full
   obj%k_l2g_map = l2g
   obj%k_g2l_map = g2l
   do ik = 1, nk_local
      obj%eigenvalues(:, ik) = eigenvalues_full(:, l2g(ik))
   end do
   obj%nk_total = nk_global
   obj%nk_local = nk_local
   obj%k_start = k_start
   obj%k_end = k_end
   obj%k_mesh_distributed_active = numprocs > 1
   obj%temperature = temperature
   obj%total_electrons = target_electrons
   obj%auto_find_fermi = .true.

   ! Independent full-array reference (does not call the production evaluator).
   call independent_fermi_root(eigenvalues_full, weights_full, target_electrons, temperature, ef_ref)
   call independent_accumulate(eigenvalues_full, weights_full, ef_ref, temperature, nelect, eband_ref)

   ef = obj%find_fermi_level_from_eigenvalues(target_electrons)
   eband = obj%calculate_canonical_band_energy(find_fermi=.true., electron_count=nelect)
   raw_weight_sum = obj%canonical_weight_sum
   call check_close('electron conservation', nelect, target_electrons, tol)
   call check_close('EF against independent root', ef, ef_ref, tol)
   call check_close('canonical EF field', obj%fermi_level, ef_ref, tol)
   call check_close('canonical EBAND', eband, eband_ref, tol)
   call check_close('raw k-weight sum', raw_weight_sum, sum(weights_full), tol)

   ! Explicit normalization contract: raw multiplicities may be rescaled.
   obj%k_weights = 7.0_rp*weights_full
   call obj%evaluate_eigenvalue_occupations(ef, nelect_scaled, eband_scaled, raw_weight_sum)
   call check_close('weight-scale invariant N', nelect_scaled, nelect, tol)
   call check_close('weight-scale invariant EBAND', eband_scaled, eband, tol)
   call check_close('scaled raw k-weight sum', raw_weight_sum, 7.0_rp*sum(weights_full), 1.0e-9_rp)
   obj%k_weights = weights_full

   ! DOS output controls must be spectators for the canonical evaluator.
   obj%n_energy_points = 17
   obj%dos_energy_range = [-0.4_rp, 0.2_rp]
   call obj%evaluate_eigenvalue_occupations(ef, nelect_again, eband_again)
   call check_close('DOS-window independent N', nelect_again, nelect, tol)
   call check_close('DOS-window independent EBAND', eband_again, eband, tol)
   obj%n_energy_points = 200003
   obj%dos_energy_range = [-1.5_rp, 1.5_rp]
   call obj%evaluate_eigenvalue_occupations(ef, nelect_again, eband_again)
   call check_close('DOS-grid independent EBAND', eband_again, eband, tol)

   ! Converged total-DOS diagnostic: narrow Gaussian representation of the same
   ! eigenvalues must approach the canonical occupied-eigenvalue sum.
   ne = obj%n_energy_points
   allocate(obj%dos_energy_grid(ne), obj%total_dos(ne))
   e_min = obj%dos_energy_range(1)
   e_max = obj%dos_energy_range(2)
   de = (e_max - e_min)/real(ne - 1, rp)
   sigma = 1.0e-4_rp
   do ie = 1, ne
      energy = e_min + de*real(ie - 1, rp)
      obj%dos_energy_grid(ie) = energy
      obj%total_dos(ie) = 0.0_rp
      do ik = 1, nk_global
         obj%total_dos(ie) = obj%total_dos(ie) + weights_full(ik)*sum(exp( &
            -0.5_rp*((energy - eigenvalues_full(:, ik))/sigma)**2))/(sigma*sqrt(2.0_rp*acos(-1.0_rp)))
      end do
      obj%total_dos(ie) = obj%total_dos(ie)/sum(weights_full)
   end do
   eband_grid = obj%calculate_band_energy_from_total_dos()
   grid_error = abs(eband_grid - eband)
   if (rank == 0) write (*, '(a,es12.4)') 'Converged total-DOS EBAND error = ', grid_error
   if (grid_error > 2.0e-5_rp) failed = .true.

   ! q=0 global spin rotation is isospectral.  Diagonalize a spin-split 2x2
   ! block before/after a dense SU(2)-equivalent rotation and feed both spectra
   ! through the canonical evaluator.
   call test_global_rotation_invariance(obj, failed)

   ! Probe-boundary invalidation must discard all spectrum-derived state while
   ! preserving the reusable k mesh/weights.
   if (.not. allocated(obj%total_dos)) allocate(obj%total_dos(3))
   allocate(obj%band_moments(1, 1, 2, 3))
   obj%canonical_energy_valid = .true.
   call obj%invalidate_spectral_cache()
   if (allocated(obj%eigenvalues) .or. allocated(obj%total_dos) .or. allocated(obj%band_moments)) failed = .true.
   if (obj%canonical_energy_valid) failed = .true.
   if (.not. allocated(obj%k_weights)) failed = .true.

#ifdef USE_MPI
   call reduce_failure(failed)
#endif
   if (rank == 0) then
      if (failed) then
         write (*, '(a)') 'RESULT: FAIL'
      else
         write (*, '(a)') 'RESULT: PASS'
      end if
   end if
#ifdef USE_MPI
   call MPI_FINALIZE(ierr)
#endif
   if (failed) error stop 1

contains

   pure real(rp) function independent_occupation(e, mu, temp) result(value)
      real(rp), intent(in) :: e, mu, temp
      real(rp) :: arg, kt_local
      kt_local = max(temp*kB_Ry_per_K, 1.0e-10_rp)
      arg = (e - mu)/kt_local
      if (arg >= 50.0_rp) then
         value = 0.0_rp
      else if (arg <= -50.0_rp) then
         value = 1.0_rp
      else
         value = 1.0_rp/(exp(arg) + 1.0_rp)
      end if
   end function independent_occupation

   subroutine independent_accumulate(evals, weights, mu, temp, count, energy_sum)
      real(rp), intent(in) :: evals(:, :), weights(:), mu, temp
      real(rp), intent(out) :: count, energy_sum
      integer :: i, j
      real(rp) :: f
      count = 0.0_rp
      energy_sum = 0.0_rp
      do j = 1, size(evals, 2)
         do i = 1, size(evals, 1)
            f = independent_occupation(evals(i, j), mu, temp)
            count = count + weights(j)*f
            energy_sum = energy_sum + weights(j)*f*evals(i, j)
         end do
      end do
      count = count/sum(weights)
      energy_sum = energy_sum/sum(weights)
   end subroutine independent_accumulate

   subroutine independent_fermi_root(evals, weights, target, temp, mu)
      real(rp), intent(in) :: evals(:, :), weights(:), target, temp
      real(rp), intent(out) :: mu
      real(rp) :: lo, hi, count, energy_sum
      integer :: iteration
      lo = minval(evals) - 1.0_rp
      hi = maxval(evals) + 1.0_rp
      do iteration = 1, 256
         mu = 0.5_rp*(lo + hi)
         call independent_accumulate(evals, weights, mu, temp, count, energy_sum)
         if (abs(count - target) < 1.0e-13_rp) exit
         if (count < target) then
            lo = mu
         else
            hi = mu
         end if
      end do
   end subroutine independent_fermi_root

   subroutine check_close(label, actual, expected, tolerance)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual, expected, tolerance
      real(rp) :: error_value
      error_value = abs(actual - expected)
      if (rank == 0) write (*, '(a,a,a,es12.4)') 'Check ', trim(label), ' error = ', error_value
      if (error_value > tolerance) failed = .true.
   end subroutine check_close

   subroutine test_global_rotation_invariance(recip, test_failed)
      type(reciprocal), intent(inout) :: recip
      logical, intent(inout) :: test_failed
      complex(rp) :: h0(2, 2), hrot(2, 2), u(2, 2), work(8)
      real(rp) :: eval0(2), evalrot(2), rwork(4), theta, e0, erot, n0, nrot
      integer :: info

      theta = 0.73_rp
      h0 = cmplx(0.0_rp, 0.0_rp, rp)
      h0(1, 1) = cmplx(-0.7_rp, 0.0_rp, rp)
      h0(2, 2) = cmplx(0.3_rp, 0.0_rp, rp)
      u = cmplx(0.0_rp, 0.0_rp, rp)
      u(1, 1) = cmplx(cos(0.5_rp*theta), 0.0_rp, rp)
      u(1, 2) = cmplx(-sin(0.5_rp*theta), 0.0_rp, rp)
      u(2, 1) = cmplx(sin(0.5_rp*theta), 0.0_rp, rp)
      u(2, 2) = cmplx(cos(0.5_rp*theta), 0.0_rp, rp)
      hrot = matmul(u, matmul(h0, transpose(conjg(u))))
      call zheev('N', 'U', 2, h0, 2, eval0, work, size(work), rwork, info)
      if (info /= 0) test_failed = .true.
      call zheev('N', 'U', 2, hrot, 2, evalrot, work, size(work), rwork, info)
      if (info /= 0) test_failed = .true.

      if (allocated(recip%eigenvalues)) deallocate(recip%eigenvalues)
      if (allocated(recip%k_weights)) deallocate(recip%k_weights)
      if (allocated(recip%k_l2g_map)) deallocate(recip%k_l2g_map)
      allocate(recip%eigenvalues(2, 1), recip%k_weights(1), recip%k_l2g_map(1))
      recip%k_weights = 9.0_rp
      recip%k_l2g_map = 1
      recip%k_mesh_distributed_active = .false.
      recip%eigenvalues(:, 1) = eval0
      call recip%evaluate_eigenvalue_occupations(0.0_rp, n0, e0)
      recip%eigenvalues(:, 1) = evalrot
      call recip%evaluate_eigenvalue_occupations(0.0_rp, nrot, erot)
      call check_close('q=0 rotation electron count', nrot, n0, 2.0e-13_rp)
      call check_close('q=0 rotation EBAND', erot, e0, 2.0e-13_rp)
   end subroutine test_global_rotation_invariance

#ifdef USE_MPI
   subroutine reduce_failure(test_failed)
      logical, intent(inout) :: test_failed
      integer :: local_failure, global_failure
      local_failure = merge(1, 0, test_failed)
      call MPI_ALLREDUCE(local_failure, global_failure, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      test_failed = global_failure /= 0
   end subroutine reduce_failure
#endif

end program test_kspace_occupations
