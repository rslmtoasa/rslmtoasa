! DRESP-07R radial SVD and exact local spherical projector regression.
program test_dresp07_radial_oracle
   use precision_mod, only: rp
   use lr_dresp07_radial_oracle_mod, only: dresp07_project_local_spherical_operator, dresp07_solve_real_least_squares
   implicit none

   integer, parameter :: norb_site = 9, nsite = 1, lmax = 2
   real(rp), parameter :: tolerance = 2.0e-11_rp
   real(rp) :: design(3, 4), v_true(4), rhs(3), solution(4), singular_values(3)
   real(rp) :: condition, residual, relative_residual, solution_norm
   real(rp) :: bad_design(4, 2), bad_rhs(4), bad_solution(2), bad_singular_values(2)
   real(rp) :: bad_condition, bad_residual, bad_relative_residual, bad_solution_norm
   real(rp) :: p_diagonal(3)
   real(rp) :: local_residual, local_relative, offdiag_norm, anisotropy_norm, cross_l_norm, intersite_norm
   real(rp) :: candidate_residual, candidate_two_residual
   complex(rp) :: target(norb_site, norb_site), best(norb_site, norb_site), candidate(norb_site, norb_site)
   complex(rp) :: candidate_two(norb_site, norb_site), coefficients(1, lmax + 1)
   integer :: rank, bad_rank, i

   design(1, :) = [1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
   design(2, :) = [0.0_rp, 1.0_rp, 1.0_rp, 0.0_rp]
   design(3, :) = [1.0_rp, 0.0_rp, 1.0_rp, 1.0_rp]
   v_true = [0.3_rp, -0.2_rp, 0.4_rp, 0.1_rp]
   rhs = matmul(design, v_true)
   call dresp07_solve_real_least_squares(design, rhs, solution, rank, singular_values, condition, residual, &
      relative_residual, solution_norm)
   if (rank /= 3 .or. relative_residual > tolerance .or. residual > tolerance .or. &
       solution_norm > sqrt(sum(v_true**2)) + tolerance .or. condition <= 1.0_rp) then
      error stop 'DRESP-07R unit: nonorthogonal representable SVD solve failed'
   end if

   bad_design(1, :) = [1.0_rp, 0.0_rp]
   bad_design(2, :) = [0.0_rp, 1.0_rp]
   bad_design(3, :) = [1.0_rp, 1.0_rp]
   bad_design(4, :) = [2.0_rp, 1.0_rp]
   bad_rhs = [1.0_rp, 2.0_rp, 4.0_rp, 8.0_rp]
   call dresp07_solve_real_least_squares(bad_design, bad_rhs, bad_solution, bad_rank, bad_singular_values, &
      bad_condition, bad_residual, bad_relative_residual, bad_solution_norm)
   if (bad_rank /= 2 .or. bad_relative_residual < 1.0e-3_rp .or. bad_residual < 1.0e-3_rp) then
      error stop 'DRESP-07R unit: nonrepresentable negative SVD fixture was accepted'
   end if

   target = cmplx(0.0_rp, 0.0_rp, rp)
   p_diagonal = [-0.2_rp, 0.0_rp, -0.4_rp]
   do i = 1, norb_site
      select case (i)
      case (1)
         target(i, i) = cmplx(0.4_rp, 0.0_rp, rp)
      case (2:4)
         target(i, i) = cmplx(p_diagonal(i - 1), 0.0_rp, rp)
      case default
         target(i, i) = cmplx(0.7_rp, 0.0_rp, rp)
      end select
   end do
   target(2, 3) = cmplx(0.05_rp, 0.02_rp, rp)
   target(3, 2) = conjg(target(2, 3))
   target(1, 2) = cmplx(0.03_rp, 0.0_rp, rp)
   target(2, 1) = conjg(target(1, 2))

   call dresp07_project_local_spherical_operator(target, norb_site, nsite, lmax, best, coefficients, local_residual, &
      local_relative, offdiag_norm, anisotropy_norm, cross_l_norm, intersite_norm)
   if (abs(coefficients(1, 1) - cmplx(0.4_rp, 0.0_rp, rp)) > tolerance .or. &
       abs(coefficients(1, 2) - cmplx(-0.2_rp, 0.0_rp, rp)) > tolerance .or. &
       abs(coefficients(1, 3) - cmplx(0.7_rp, 0.0_rp, rp)) > tolerance) then
      error stop 'DRESP-07R unit: local spherical projector coefficients are not exact shell averages'
   end if
   if (abs(best(2, 2) - best(3, 3)) > tolerance .or. abs(best(2, 2) - best(4, 4)) > tolerance) then
      error stop 'DRESP-07R unit: local spherical projector is not m-degenerate'
   end if
   if (offdiag_norm <= 0.0_rp .or. anisotropy_norm <= 0.0_rp .or. cross_l_norm <= 0.0_rp .or. intersite_norm > tolerance) then
      error stop 'DRESP-07R unit: projector-space defect decomposition failed'
   end if

   candidate = cmplx(0.0_rp, 0.0_rp, rp)
   candidate_two = cmplx(0.0_rp, 0.0_rp, rp)
   candidate(1, 1) = cmplx(0.1_rp, 0.0_rp, rp)
   candidate(2:4, 2:4) = cmplx(-0.3_rp, 0.0_rp, rp)
   candidate(5:9, 5:9) = cmplx(0.5_rp, 0.0_rp, rp)
   candidate_two(1, 1) = cmplx(0.8_rp, 0.0_rp, rp)
   candidate_two(2:4, 2:4) = cmplx(0.4_rp, 0.0_rp, rp)
   candidate_two(5:9, 5:9) = cmplx(0.9_rp, 0.0_rp, rp)
   candidate_residual = sqrt(sum(abs(target - candidate)**2))
   candidate_two_residual = sqrt(sum(abs(target - candidate_two)**2))
   if (local_residual > candidate_residual + tolerance .or. local_residual > candidate_two_residual + tolerance) then
      error stop 'DRESP-07R unit: best local spherical operator violates variational inequality'
   end if

   write (*, '(a,es12.4)') '  nonorthogonal_svd_relative=', relative_residual
   write (*, '(a,es12.4)') '  projector_relative_residual=', local_relative
   write (*, '(a,es12.4)') '  negative_fixture_relative_residual=', bad_relative_residual
   write (*, '(a)') 'UnitDresp07RadialOracle: PASS (SVD, projector, variational, and negative fixtures)'
end program test_dresp07_radial_oracle
