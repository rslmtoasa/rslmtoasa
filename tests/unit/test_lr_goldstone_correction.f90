!------------------------------------------------------------------------------
! GCR-01 BES Goldstone eigenvalue-correction oracles.
!
! The fixtures construct canonical LR-04 chiKS and Kxc matrices directly.
! They exercise all three denominator classes and keep the optional route
! separate from both the susceptibility and the direct KXC evaluator.
!------------------------------------------------------------------------------
program test_lr_goldstone_correction
   use precision_mod, only: rp
   use lr_response_space_mod, only: response_space_layout
   use lr_goldstone_correction_mod, only: lr_goldstone_correction_request, lr_goldstone_correction_result, &
      lr_gcr_units, lr_gcr_kxc_provenance, lr_gcr_class_hermitian, lr_gcr_class_metric_hermitian, &
      lr_gcr_class_nonhermitian, lr_gcr_solver_hermitian, lr_gcr_solver_metric, lr_gcr_solver_nonhermitian, &
      evaluate_lr_goldstone_correction
   implicit none

   integer, parameter :: npoint = 5
   real(rp), parameter :: mesh_a = 0.04_rp, mesh_b = 0.09_rp
   real(rp), parameter :: tolerance = 2.0e-9_rp
   type(response_space_layout), target :: space
   type(lr_goldstone_correction_request) :: request
   type(lr_goldstone_correction_result) :: exact_result, hermitian_result, perturbed_result, wrong_result, &
      ambiguous_result, nonhermitian_result, provenance_result, disabled_result
   real(rp), allocatable :: radius(:)
   integer, allocatable :: active(:)
   complex(rp), allocatable :: chi(:, :), kernel(:, :), rigid(:)
   logical :: failed

   failed = .false.
   allocate(radius(npoint))
   call build_mesh(radius)
   call space%initialize(1, 0, radius, mesh_a, mesh_b, 1)
   call collect_active(space, active)

   ! Optional means a default request is an observable no-op.
   call evaluate_lr_goldstone_correction(request, disabled_result)
   if (disabled_result%corrected .or. disabled_result%selected .or. &
       index(disabled_result%status, 'DISABLED') == 0) then
      failed = .true.
      write (*, '(a)') 'Disabled correction was not a clean no-op'
   end if
   write (*, '(a)') 'Optional-disabled route: PASS'

   ! 1. Exact ordinary-Hermitian Goldstone system.
   call build_fixture(space, active, 'hermitian', 0.0_rp, chi, kernel, rigid)
   call prepare_request(request, space, chi, kernel, rigid, .true.)
   call evaluate_lr_goldstone_correction(request, hermitian_result)
   call check_pass(hermitian_result, 'exact Hermitian Goldstone correction', failed)
   if (trim(hermitian_result%representation_class) /= lr_gcr_class_hermitian .or. &
       trim(hermitian_result%eigensolver) /= lr_gcr_solver_hermitian) then
      failed = .true.
      write (*, '(a)') 'Hermitian denominator did not use zheev'
   end if

   ! 2. Exact metric-Hermitian Goldstone system.  The first mode is exact and
   ! the second/third modes are mixed, so ordinary Hermiticity is not assumed.
   call build_fixture(space, active, 'metric', 0.0_rp, chi, kernel, rigid)
   call prepare_request(request, space, chi, kernel, rigid, .true.)
   call evaluate_lr_goldstone_correction(request, exact_result)
   call check_pass(exact_result, 'exact metric-Hermitian Goldstone correction', failed)
   if (trim(exact_result%representation_class) /= lr_gcr_class_metric_hermitian .or. &
       trim(exact_result%eigensolver) /= lr_gcr_solver_metric) then
      failed = .true.
      write (*, '(a)') 'Metric-Hermitian denominator did not use the metric solver'
   end if
   call check_real(exact_result%corrected_rigid_residual, 0.0_rp, tolerance, &
      'exact metric-aware rigid mode', failed)
   call check_real(exact_result%untouched_eigenvalue_change, 0.0_rp, tolerance, &
      'exact untouched eigenvalues', failed)
   write (*, '(a,es12.4,a,es12.4)') 'Exact metric Goldstone |lambda|/overlap = ', &
      exact_result%goldstone_eigenvalue_abs, ' ', exact_result%rigid_overlap

   ! 3. Controlled perturbation of the single Goldstone eigenvalue.
   call build_fixture(space, active, 'metric', 2.0e-3_rp, chi, kernel, rigid)
   call prepare_request(request, space, chi, kernel, rigid, .true., 2.0e-3_rp)
   call evaluate_lr_goldstone_correction(request, perturbed_result)
   call check_pass(perturbed_result, 'controlled single-mode perturbation', failed)
   call check_real(perturbed_result%goldstone_eigenvalue_abs, 2.0e-3_rp, tolerance, &
      'controlled raw Goldstone eigenvalue', failed)
   call check_real(perturbed_result%untouched_eigenvalue_change, 0.0_rp, tolerance, &
      'controlled untouched eigenvalues', failed)
   call check_real(perturbed_result%corrected_rigid_residual, 0.0_rp, tolerance, &
      'controlled metric-aware rigid mode', failed)
   call check_kernel_reconstruction(perturbed_result, tolerance, 'controlled BES kernel reconstruction', failed)

   ! 4. A small eigenvalue with the wrong eigenvector character is refused.
   call build_fixture(space, active, 'wrong', 0.0_rp, chi, kernel, rigid)
   call prepare_request(request, space, chi, kernel, rigid, .true.)
   call evaluate_lr_goldstone_correction(request, wrong_result)
   if (.not. wrong_result%blocked .or. wrong_result%corrected .or. &
       index(wrong_result%status, 'weak') == 0) then
      failed = .true.
      write (*, '(a)') 'Wrong-small-eigenvalue fixture was not refused'
   end if
   write (*, '(a,a)') 'Wrong-small-eigenvalue fixture: ', trim(wrong_result%status)

   ! 5. Two nearby small eigenvalues are ambiguous even when one rigid overlap
   ! would otherwise be strong.
   call build_fixture(space, active, 'ambiguous', 0.0_rp, chi, kernel, rigid)
   call prepare_request(request, space, chi, kernel, rigid, .true.)
   call evaluate_lr_goldstone_correction(request, ambiguous_result)
   if (.not. ambiguous_result%blocked .or. ambiguous_result%corrected .or. &
       index(ambiguous_result%status, 'isolated') == 0) then
      failed = .true.
      write (*, '(a)') 'Two-small-eigenvalue fixture was not refused'
   end if
   write (*, '(a,a)') 'Two-small-eigenvalue fixture: ', trim(ambiguous_result%status)

   ! 6. Complex non-Hermitian D uses zgeev and an explicit right-eigenvector
   ! inverse in the BES reconstruction.
   call build_fixture(space, active, 'nonhermitian', 2.0e-3_rp, chi, kernel, rigid)
   call prepare_request(request, space, chi, kernel, rigid, .true., 2.0e-3_rp)
   call evaluate_lr_goldstone_correction(request, nonhermitian_result)
   call check_pass(nonhermitian_result, 'complex non-Hermitian correction', failed)
   if (trim(nonhermitian_result%representation_class) /= lr_gcr_class_nonhermitian .or. &
       trim(nonhermitian_result%eigensolver) /= lr_gcr_solver_nonhermitian) then
      failed = .true.
      write (*, '(a)') 'Non-Hermitian denominator did not use zgeev/right-inverse'
   end if
   call check_kernel_reconstruction(nonhermitian_result, tolerance, 'non-Hermitian BES kernel reconstruction', failed)

   ! 7. A Kxc provenance mismatch is a blocker, not a silently corrected route.
   call build_fixture(space, active, 'metric', 2.0e-3_rp, chi, kernel, rigid)
   call prepare_request(request, space, chi, kernel, rigid, .true., 2.0e-3_rp)
   request%kernel_provenance = 'unrelated kernel'
   call evaluate_lr_goldstone_correction(request, provenance_result)
   if (.not. provenance_result%blocked .or. provenance_result%corrected .or. &
       index(provenance_result%status, 'provenance') == 0) then
      failed = .true.
      write (*, '(a)') 'Kxc provenance mismatch was not refused'
   end if
   write (*, '(a,a)') 'Kxc provenance mismatch: ', trim(provenance_result%status)

   if (failed) then
      write (*, '(a)') 'UnitLrGoldstoneCorrection: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrGoldstoneCorrection: PASS (optional, metric, non-Hermitian, gates, BES reconstruction)'

contains

   subroutine build_mesh(mesh)
      real(rp), intent(out) :: mesh(:)
      integer :: i

      mesh(1) = 0.0_rp
      do i = 2, size(mesh)
         mesh(i) = mesh_b*(exp(mesh_a*real(i - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine collect_active(layout, indices)
      type(response_space_layout), intent(in) :: layout
      integer, allocatable, intent(out) :: indices(:)
      integer :: i, n

      allocate(indices(count(layout%metric_weights > 0.0_rp)))
      n = 0
      do i = 1, layout%ndim
         if (layout%metric_weights(i) > 0.0_rp) then
            n = n + 1
            indices(n) = i
         end if
      end do
   end subroutine collect_active

   subroutine build_fixture(layout, active_indices, fixture, goldstone_value, susceptibility, kernel, rigid)
      type(response_space_layout), intent(in) :: layout
      integer, intent(in) :: active_indices(:)
      character(len=*), intent(in) :: fixture
      real(rp), intent(in) :: goldstone_value
      complex(rp), allocatable, intent(out) :: susceptibility(:, :), kernel(:, :), rigid(:)
      complex(rp), allocatable :: denominator(:, :), metric_form(:, :)
      real(rp), allocatable :: sqrtw(:)
      integer :: n, i, j

      n = layout%ndim
      allocate(susceptibility(n, n), kernel(n, n), rigid(n), denominator(n, n))
      susceptibility = cmplx(0.0_rp, 0.0_rp, rp)
      denominator = cmplx(0.0_rp, 0.0_rp, rp)
      rigid = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         susceptibility(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
         denominator(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
      rigid(active_indices(1)) = cmplx(1.0_rp, 0.0_rp, rp)

      select case (trim(fixture))
      case ('hermitian')
         denominator(active_indices, active_indices) = cmplx(0.0_rp, 0.0_rp, rp)
         denominator(active_indices(1), active_indices(1)) = cmplx(goldstone_value, 0.0_rp, rp)
         denominator(active_indices(2), active_indices(2)) = cmplx(0.8_rp, 0.0_rp, rp)
         denominator(active_indices(3), active_indices(3)) = cmplx(1.2_rp, 0.0_rp, rp)
         denominator(active_indices(4), active_indices(4)) = cmplx(1.6_rp, 0.0_rp, rp)
      case ('metric')
         allocate(sqrtw(size(active_indices)), metric_form(size(active_indices), size(active_indices)))
         do i = 1, size(active_indices)
            sqrtw(i) = sqrt(layout%metric_weights(active_indices(i)))
         end do
         metric_form = cmplx(0.0_rp, 0.0_rp, rp)
         metric_form(1, 1) = cmplx(goldstone_value, 0.0_rp, rp)
         metric_form(2, 2) = cmplx(0.8_rp, 0.0_rp, rp)
         metric_form(3, 3) = cmplx(1.2_rp, 0.0_rp, rp)
         metric_form(4, 4) = cmplx(1.6_rp, 0.0_rp, rp)
         metric_form(2, 3) = cmplx(0.15_rp, 0.0_rp, rp)
         metric_form(3, 2) = metric_form(2, 3)
         do i = 1, size(active_indices)
            do j = 1, size(active_indices)
               denominator(active_indices(i), active_indices(j)) = metric_form(i, j)*sqrtw(j)/sqrtw(i)
            end do
         end do
      case ('wrong')
         denominator(active_indices, active_indices) = cmplx(0.0_rp, 0.0_rp, rp)
         denominator(active_indices(1), active_indices(1)) = cmplx(5.0e-4_rp, 0.0_rp, rp)
         denominator(active_indices(2), active_indices(2)) = cmplx(0.8_rp, 0.0_rp, rp)
         denominator(active_indices(3), active_indices(3)) = cmplx(1.2_rp, 0.0_rp, rp)
         denominator(active_indices(4), active_indices(4)) = cmplx(1.6_rp, 0.0_rp, rp)
         rigid(active_indices(1)) = cmplx(0.0_rp, 0.0_rp, rp)
         rigid(active_indices(2)) = cmplx(1.0_rp, 0.0_rp, rp)
      case ('ambiguous')
         denominator(active_indices, active_indices) = cmplx(0.0_rp, 0.0_rp, rp)
         denominator(active_indices(1), active_indices(1)) = cmplx(1.0e-3_rp, 0.0_rp, rp)
         denominator(active_indices(2), active_indices(2)) = cmplx(2.0e-3_rp, 0.0_rp, rp)
         denominator(active_indices(3), active_indices(3)) = cmplx(1.0_rp, 0.0_rp, rp)
         denominator(active_indices(4), active_indices(4)) = cmplx(1.6_rp, 0.0_rp, rp)
      case ('nonhermitian')
         denominator(active_indices, active_indices) = cmplx(0.0_rp, 0.0_rp, rp)
         denominator(active_indices(1), active_indices(1)) = cmplx(goldstone_value, 0.0_rp, rp)
         denominator(active_indices(2), active_indices(2)) = cmplx(0.8_rp, 0.0_rp, rp)
         denominator(active_indices(3), active_indices(3)) = cmplx(1.2_rp, 0.0_rp, rp)
         denominator(active_indices(4), active_indices(4)) = cmplx(1.6_rp, 0.0_rp, rp)
         denominator(active_indices(1), active_indices(2)) = cmplx(0.2_rp, 0.05_rp, rp)
         denominator(active_indices(2), active_indices(3)) = cmplx(0.1_rp, -0.03_rp, rp)
      case default
         error stop 'build_fixture: unknown fixture'
      end select
      kernel = susceptibility - denominator
   end subroutine build_fixture

   subroutine prepare_request(request, layout, susceptibility, kernel, rigid, enable, goldstone_value)
      type(lr_goldstone_correction_request), intent(inout) :: request
      type(response_space_layout), target, intent(in) :: layout
      complex(rp), intent(in) :: susceptibility(:, :), kernel(:, :), rigid(:)
      logical, intent(in) :: enable
      real(rp), intent(in), optional :: goldstone_value
      real(rp) :: sample_goldstone

      sample_goldstone = 0.0_rp
      if (present(goldstone_value)) sample_goldstone = goldstone_value

      if (allocated(request%static_susceptibility)) deallocate(request%static_susceptibility)
      if (allocated(request%xc_kernel)) deallocate(request%xc_kernel)
      if (allocated(request%rigid_rotation_vector)) deallocate(request%rigid_rotation_vector)
      if (allocated(request%consistency_samples)) deallocate(request%consistency_samples)
      request%response_space => layout
      allocate(request%static_susceptibility(size(susceptibility, 1), size(susceptibility, 2)), &
         request%xc_kernel(size(kernel, 1), size(kernel, 2)), request%rigid_rotation_vector(size(rigid)))
      request%static_susceptibility = susceptibility
      request%xc_kernel = kernel
      request%rigid_rotation_vector = rigid
      request%enabled = enable
      request%lr06_converged = .true.
      request%static_q_zero = .true.
      request%static_frequency = 0.0_rp
      request%eta = 0.02_rp
      request%k_mesh_points = 8
      request%kernel_units = lr_gcr_units
      request%kernel_provenance = lr_gcr_kxc_provenance
      request%kernel_provenance_verified = .true.
      request%soc_present = .false.
      request%external_field_present = .false.
      request%constraining_field_present = .false.
      allocate(request%consistency_samples(2))
      request%consistency_samples(1)%eta = request%eta
      request%consistency_samples(1)%k_mesh_points = request%k_mesh_points
      request%consistency_samples(2)%eta = 0.03_rp
      request%consistency_samples(2)%k_mesh_points = 10
      request%consistency_samples(:)%goldstone_eigenvalue = cmplx(sample_goldstone, 0.0_rp, rp)
      request%consistency_samples(:)%competing_eigenvalue_abs = 0.8_rp
      request%consistency_samples(:)%normalized_rigid_overlap = 1.0_rp
   end subroutine prepare_request

   subroutine check_pass(result, label, test_failed)
      type(lr_goldstone_correction_result), intent(in) :: result
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed

      if (result%blocked .or. .not. result%corrected) then
         test_failed = .true.
         write (*, '(a,a)') trim(label)//': ', trim(result%status)
      else
         write (*, '(a,a)') trim(label)//': ', trim(result%status)
      end if
   end subroutine check_pass

   subroutine check_real(actual, expected, tolerance_value, label, test_failed)
      real(rp), intent(in) :: actual, expected, tolerance_value
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed

      write (*, '(a,es12.4)') trim(label)//' error = ', abs(actual - expected)
      if (abs(actual - expected) > tolerance_value) test_failed = .true.
   end subroutine check_real

   subroutine check_kernel_reconstruction(result, tolerance_value, label, test_failed)
      type(lr_goldstone_correction_result), intent(in) :: result
      real(rp), intent(in) :: tolerance_value
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed

      write (*, '(a,es12.4)') trim(label)//' error = ', result%kernel_reconstruction_residual
      if (result%kernel_reconstruction_residual > tolerance_value) test_failed = .true.
   end subroutine check_kernel_reconstruction

end program test_lr_goldstone_correction
