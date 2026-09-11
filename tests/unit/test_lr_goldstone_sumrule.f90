!------------------------------------------------------------------------------
! GSR-01 independent Lounis/Costa/Muniz/Mills sum-rule oracles.
!
! The fixtures construct canonical LR-04 chiKS matrices directly.  They do
! not call the direct ALSDA implementation to obtain the sum-rule solution.
!------------------------------------------------------------------------------
program test_lr_goldstone_sumrule
   use precision_mod, only: rp
   use response_angular_basis_mod, only: response_angular_pi
   use response_basis_mapping_mod, only: response_super_index, response_flatten_superindex
   use lr_response_space_mod, only: response_space_layout, response_apply_operator
   use lr_goldstone_sumrule_mod, only: lr_goldstone_sumrule_result, lr_sumrule_linear_solve_result, &
      lr_gsr_magnetization_pauli, evaluate_lr_goldstone_sumrule, solve_lr_goldstone_equation
   implicit none

   integer, parameter :: nr = 5, nsite = 2
   real(rp), parameter :: mesh_a = 0.04_rp, mesh_b = 0.09_rp
   real(rp), parameter :: tolerance = 2.0e-9_rp
   real(rp), parameter :: four_pi = 4.0_rp*response_angular_pi
   type(response_space_layout) :: space
   type(lr_goldstone_sumrule_result) :: result, one_site_result, multisite_result, &
      manufactured_result, metric_result
   type(lr_sumrule_linear_solve_result) :: analytic_solve, singular_solve
   complex(rp), allocatable :: gamma(:, :), rhs(:), expected(:), singular_gamma(:, :), singular_rhs(:)
   complex(rp), allocatable :: chi(:, :), target(:), generated(:), expected_residual(:)
   real(rp), allocatable :: radius(:), magnetization(:, :), direct_u(:, :)
   complex(rp), allocatable :: expected_u(:, :), manufactured_gamma(:, :)
   integer, allocatable :: active(:)
   integer :: n, i, j, isite, ir, flat, active_index
   real(rp) :: metric_expected, direct_difference
   complex(rp) :: off_diagonal
   logical :: failed

   failed = .false.
   allocate(radius(nr))
   call build_mesh(radius)

   ! 1. Independent analytic finite Gamma/U/m equation.
   allocate(gamma(3, 3), rhs(3), expected(3))
   gamma = reshape([cmplx(1.4_rp, 0.2_rp, rp), cmplx(-0.3_rp, 0.1_rp, rp), cmplx(0.2_rp, -0.1_rp, rp), &
                    cmplx(0.4_rp, -0.2_rp, rp), cmplx(1.7_rp, 0.0_rp, rp), cmplx(0.1_rp, 0.3_rp, rp), &
                    cmplx(-0.2_rp, 0.1_rp, rp), cmplx(0.5_rp, -0.2_rp, rp), cmplx(1.1_rp, 0.4_rp, rp)], [3, 3])
   expected = [cmplx(0.4_rp, -0.2_rp, rp), cmplx(-0.7_rp, 0.3_rp, rp), cmplx(0.9_rp, 0.1_rp, rp)]
   rhs = matmul(gamma, expected)
   call solve_lr_goldstone_equation(gamma, rhs, analytic_solve)
   call check_vector(analytic_solve%solution, expected, tolerance, 'analytic Gamma/U/m solve', failed)
   if (analytic_solve%rank /= 3 .or. analytic_solve%blocked) then
      failed = .true.
      write (*, '(a)') 'Analytic full-rank solve did not pass its rank gate'
   end if
   write (*, '(a,es12.4,a,es12.4,a,i0)') 'Analytic solve condition/residual/rank = ', &
      analytic_solve%condition_number, ' ', analytic_solve%residual_norm, ' ', analytic_solve%rank

   ! Explicit rank deficiency is a diagnostic blocker, never a smoothed answer.
   allocate(singular_gamma(2, 2), singular_rhs(2))
   singular_gamma = reshape([cmplx(1.0_rp, 0.0_rp, rp), cmplx(2.0_rp, 0.0_rp, rp), &
                              cmplx(2.0_rp, 0.0_rp, rp), cmplx(4.0_rp, 0.0_rp, rp)], [2, 2])
   singular_rhs = [cmplx(1.0_rp, 0.0_rp, rp), cmplx(2.0_rp, 0.0_rp, rp)]
   call solve_lr_goldstone_equation(singular_gamma, singular_rhs, singular_solve)
   if (.not. singular_solve%rank_deficient .or. .not. singular_solve%blocked) then
      failed = .true.
      write (*, '(a)') 'Rank-deficient Gamma was not declared BLOCKED'
   end if
   write (*, '(a,i0,a,a)') 'Rank-deficient solve rank = ', singular_solve%rank, ', status = ', trim(singular_solve%status)

   ! 2. One-site radial fixture with a known local LCMM U.
   call space%initialize(1, 0, radius, mesh_a, mesh_b, 1)
   n = space%ndim
   allocate(magnetization(1, nr), expected_u(1, nr), chi(n, n))
   call fill_magnetization(magnetization)
   expected_u = cmplx(0.0_rp, 0.0_rp, rp)
   expected_u(1, 2:nr) = [cmplx(0.31_rp, 0.02_rp, rp), cmplx(0.36_rp, -0.01_rp, rp), &
                           cmplx(0.43_rp, 0.03_rp, rp), cmplx(0.49_rp, -0.02_rp, rp)]
   chi = cmplx(0.0_rp, 0.0_rp, rp)
   do ir = 2, nr
      chi(ir, ir) = 1.0_rp/(four_pi*expected_u(1, ir))
   end do
   call evaluate_lr_goldstone_sumrule(space, chi, magnetization, one_site_result, lr_gsr_magnetization_pauli)
   call check_matrix(one_site_result%u_lcmm, expected_u, tolerance, 'one-site radial U', failed)
   call check_real(one_site_result%metric_relative_residual, 0.0_rp, tolerance, &
      'one-site metric Goldstone residual', failed)
   if (one_site_result%blocked .or. one_site_result%rank /= nr - 1) then
      failed = .true.
      write (*, '(a)') 'One-site radial sum-rule fixture did not pass'
   end if
   call check_canonical_interaction(space, one_site_result, tolerance, 'one-site canonical K_eff action', failed)
   write (*, '(a,es12.4,a,es12.4)') 'One-site metric residual/condition = ', &
      one_site_result%metric_residual_norm, ' ', one_site_result%condition_number

   ! 3. Multisite fixture with unequal accepted moments and local U values.
   call space%initialize(nsite, 0, radius, mesh_a, mesh_b, 1)
   n = space%ndim
   deallocate(magnetization, expected_u, chi)
   allocate(magnetization(nsite, nr), expected_u(nsite, nr), chi(n, n))
   call fill_multisite_magnetization(magnetization)
   expected_u = cmplx(0.0_rp, 0.0_rp, rp)
   do isite = 1, nsite
      do ir = 2, nr
         expected_u(isite, ir) = cmplx(0.18_rp + 0.07_rp*real(isite, rp) + 0.01_rp*real(ir, rp), &
                                      0.005_rp*real(isite - ir, rp), rp)
      end do
   end do
   chi = cmplx(0.0_rp, 0.0_rp, rp)
   do isite = 1, nsite
      do ir = 2, nr
         flat = (isite - 1)*nr + ir
         chi(flat, flat) = 1.0_rp/(four_pi*expected_u(isite, ir))
      end do
   end do
   call evaluate_lr_goldstone_sumrule(space, chi, magnetization, multisite_result, lr_gsr_magnetization_pauli)
   call check_matrix(multisite_result%u_lcmm, expected_u, tolerance, 'multisite unequal-moment U', failed)
   if (multisite_result%blocked .or. multisite_result%rank /= nsite*(nr - 1)) then
      failed = .true.
      write (*, '(a)') 'Multisite unequal-moment fixture did not pass'
   end if
   write (*, '(a,es12.4)') 'Multisite metric Goldstone residual = ', multisite_result%metric_residual_norm

   ! 4. Exact manufactured nonlocal rigid-rotation problem.  This checks that
   ! the route solves Gamma built from chiKS rather than assuming a local chi.
   allocate(manufactured_gamma(nsite*(nr - 1), nsite*(nr - 1)), active(nsite*(nr - 1)))
   active = 0
   do isite = 1, nsite
      do ir = 2, nr
         active_index = (isite - 1)*(nr - 1) + ir - 1
         active(active_index) = (isite - 1)*nr + ir
      end do
   end do
   manufactured_gamma = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, size(active)
      do j = 1, size(active)
         if (i /= j) then
            manufactured_gamma(i, j) = cmplx(0.01_rp*real(i + 2*j, rp), &
                                               -0.003_rp*real(i - j, rp), rp)
         end if
      end do
      off_diagonal = cmplx(0.0_rp, 0.0_rp, rp)
      do j = 1, size(active)
         if (j /= i) off_diagonal = off_diagonal + manufactured_gamma(i, j)*expected_u( &
            (active(j) - 1)/nr + 1, mod(active(j) - 1, nr) + 1)
      end do
      manufactured_gamma(i, i) = (sqrt(four_pi)*magnetization((active(i) - 1)/nr + 1, &
         mod(active(i) - 1, nr) + 1) - off_diagonal)/ &
         expected_u((active(i) - 1)/nr + 1, mod(active(i) - 1, nr) + 1)
   end do
   chi = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, size(active)
      do j = 1, size(active)
         chi(active(i), active(j)) = manufactured_gamma(i, j)/ &
            (four_pi*sqrt(four_pi)*magnetization((active(j) - 1)/nr + 1, mod(active(j) - 1, nr) + 1))
      end do
   end do
   call evaluate_lr_goldstone_sumrule(space, chi, magnetization, manufactured_result, lr_gsr_magnetization_pauli)
   call check_matrix(manufactured_result%u_lcmm, expected_u, tolerance, 'manufactured nonlocal rigid U', failed)
   if (manufactured_result%blocked .or. manufactured_result%metric_residual_norm > tolerance) then
      failed = .true.
      write (*, '(a)') 'Manufactured rigid-rotation identity did not pass'
   end if
   write (*, '(a,es12.4,a,es12.4)') 'Manufactured identity residual/overlap = ', &
      manufactured_result%metric_residual_norm, ' ', abs(manufactured_result%rigid_overlap)

   ! 5. Deliberately non-spherical output tests the LR-04 metric-aware residual.
   call space%initialize(1, 1, radius, mesh_a, mesh_b, 1)
   n = space%ndim
   deallocate(magnetization, expected_u, chi, active, manufactured_gamma)
   allocate(magnetization(1, nr), expected_u(1, nr), chi(n, n), target(n), generated(n), expected_residual(n))
   call fill_magnetization(magnetization)
   expected_u = cmplx(0.0_rp, 0.0_rp, rp)
   expected_u(1, 2:nr) = [cmplx(0.27_rp, 0.0_rp, rp), cmplx(0.29_rp, 0.0_rp, rp), &
                           cmplx(0.33_rp, 0.0_rp, rp), cmplx(0.38_rp, 0.0_rp, rp)]
   chi = cmplx(0.0_rp, 0.0_rp, rp)
   do ir = 2, nr
      chi(ir, ir) = 1.0_rp/(four_pi*expected_u(1, ir))
      ! The first L=1 row is coupled to the spherical source.  The equation
      ! solve remains exact in L=0, but the full response identity is not.
      chi(nr + ir, ir) = 0.21_rp + 0.01_rp*real(ir, rp)
   end do
   call evaluate_lr_goldstone_sumrule(space, chi, magnetization, metric_result, lr_gsr_magnetization_pauli)
   expected_residual = metric_result%residual_vector
   metric_expected = sqrt(sum(space%metric_weights*abs(expected_residual)**2))
   call check_real(metric_result%metric_residual_norm, metric_expected, tolerance, &
      'metric-aware full-response residual', failed)
   if (.not. metric_result%blocked .or. metric_result%metric_residual_norm <= tolerance) then
      failed = .true.
      write (*, '(a)') 'Non-spherical manufactured residual was not reported as BLOCKED'
   end if
   write (*, '(a,es12.4,a,a)') 'Metric-aware residual/identity status = ', metric_result%metric_residual_norm, ' ', &
      trim(metric_result%status)

   ! 6. Direct ALSDA is a comparison-only diagnostic.  It is not used to
   ! construct any Gamma, and the independent values intentionally differ.
   allocate(direct_u(1, nr))
   do ir = 1, nr
      direct_u(1, ir) = 0.0_rp
      if (ir > 1) direct_u(1, ir) = real(-2.0_rp*0.06_rp*exp(-0.2_rp*radius(ir))/ &
         (four_pi*magnetization(1, ir)), rp)
   end do
   direct_difference = maxval(abs(real(metric_result%u_lcmm, rp) - direct_u))
   write (*, '(a,es12.4)') 'Direct ALSDA comparison diagnostic (not used in solve) = ', direct_difference

   if (failed) then
      write (*, '(a)') 'UnitLrGoldstoneSumrule: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitLrGoldstoneSumrule: PASS (analytic, radial, multisite, rigid, metric, independence)'

contains

   subroutine build_mesh(mesh)
      real(rp), intent(out) :: mesh(:)
      integer :: ir

      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine fill_magnetization(values)
      real(rp), intent(out) :: values(:, :)
      integer :: ir

      do ir = 1, size(values, 2)
         values(1, ir) = 0.8_rp*exp(-0.4_rp*radius(ir))
      end do
      values(1, 1) = 0.0_rp
   end subroutine fill_magnetization

   subroutine fill_multisite_magnetization(values)
      real(rp), intent(out) :: values(:, :)
      integer :: isite, ir

      do isite = 1, size(values, 1)
         do ir = 1, size(values, 2)
            values(isite, ir) = (0.55_rp + 0.37_rp*real(isite, rp))*exp(-0.35_rp*radius(ir))
         end do
         values(isite, 1) = 0.0_rp
      end do
   end subroutine fill_multisite_magnetization

   subroutine check_vector(actual, expected_value, tolerance_value, label, test_failed)
      complex(rp), intent(in) :: actual(:), expected_value(:)
      real(rp), intent(in) :: tolerance_value
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed
      real(rp) :: error_value

      error_value = maxval(abs(actual - expected_value))
      write (*, '(a,es12.4)') trim(label)//' max error = ', error_value
      if (error_value > tolerance_value) test_failed = .true.
   end subroutine check_vector

   subroutine check_matrix(actual, expected_value, tolerance_value, label, test_failed)
      complex(rp), intent(in) :: actual(:, :), expected_value(:, :)
      real(rp), intent(in) :: tolerance_value
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed
      real(rp) :: error_value

      error_value = maxval(abs(actual - expected_value))
      write (*, '(a,es12.4)') trim(label)//' max error = ', error_value
      if (error_value > tolerance_value) test_failed = .true.
   end subroutine check_matrix

   subroutine check_real(actual, expected_value, tolerance_value, label, test_failed)
      real(rp), intent(in) :: actual, expected_value, tolerance_value
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed

      write (*, '(a,es12.4)') trim(label)//' error = ', abs(actual - expected_value)
      if (abs(actual - expected_value) > tolerance_value) test_failed = .true.
   end subroutine check_real

   subroutine check_canonical_interaction(layout, gsr, tolerance_value, label, test_failed)
      type(response_space_layout), intent(in) :: layout
      type(lr_goldstone_sumrule_result), intent(in) :: gsr
      real(rp), intent(in) :: tolerance_value
      character(len=*), intent(in) :: label
      logical, intent(inout) :: test_failed
      complex(rp), allocatable :: source(:), action(:)
      integer :: site, radial_point, flat_index
      type(response_super_index) :: item

      allocate(source(layout%ndim), action(layout%ndim))
      source = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, layout%nsite
         do radial_point = 2, layout%npoint
            item = response_super_index(site, 0, 0, radial_point, 1)
            call response_flatten_superindex(item, layout%nsite, layout%response_lmax, layout%npoint, &
                                             layout%nchannel, flat_index)
            source(flat_index) = gsr%magnetization_response(flat_index)
         end do
      end do
      call response_apply_operator(layout, gsr%canonical_interaction, source, action)
      call check_vector(action, gsr%generated_field, tolerance_value, label, test_failed)
   end subroutine check_canonical_interaction

end program test_lr_goldstone_sumrule
