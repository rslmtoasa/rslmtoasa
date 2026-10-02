!------------------------------------------------------------------------------
! Independent oracles for the native RS-LMTO -> Mills-1U model projection.
!------------------------------------------------------------------------------
module test_mills_1u_mod
   use precision_mod, only: rp
   use symbolic_atom_mod, only: symbolic_atom
   use linear_response_mod, only: lr_electronic_state, mills_1u_bare_request, mills_1u_bare_result, &
      mills_1u_model_reduction_result, projected_dyson_result, mills_1u_splitting_from_centers, &
      mills_1u_coefficient_moment, mills_1u_transition_vertex, evaluate_mills_1u_bare, &
      evaluate_mills_1u_model_reduction, evaluate_mills_1u_dyson
   implicit none
   private
   public :: run_test_mills_1u

contains

   subroutine run_test_mills_1u()
      logical :: failed
      failed = .false.
      call native_splitting_oracle(failed)
      call coefficient_moment_oracle(failed)
      call transition_vertex_oracle(failed)
      call analytic_dyson_oracle(failed)
      call model_reduction_oracle(failed)
      if (failed) then
         write(*, '(a)') 'UnitLrMills1U: FAIL'
         error stop 1
      end if
      write(*, '(a)') 'UnitLrMills1U: PASS (native centers, coefficient moment/vertex, analytic Dyson, model residual)'
   end subroutine run_test_mills_1u

   subroutine native_splitting_oracle(failed)
      logical, intent(inout) :: failed
      type(symbolic_atom) :: atom
      real(rp) :: centers(3, 2, 1), delta(1), scalarity(1), expected, full_d_average
      complex(rp) :: cx(9, 2, 1), full_spin_difference(18, 18)
      integer :: m

      atom%potential%lmax = 2
      allocate(atom%potential%center_band(3, 2), atom%potential%width_band(3, 2), &
         atom%potential%shifted_band(3, 2), atom%potential%obar(3, 2), &
         atom%potential%cx(9, 2), atom%potential%wx(9, 2), atom%potential%cex(9, 2), atom%potential%obx(9, 2), &
         atom%potential%cx0(9), atom%potential%cx1(9), atom%potential%wx0(9), atom%potential%wx1(9), &
         atom%potential%cex0(9), atom%potential%cex1(9), atom%potential%obx0(9), atom%potential%obx1(9))
      atom%potential%center_band(:, 1) = [0.82_rp, 0.31_rp, 0.04_rp]
      atom%potential%center_band(:, 2) = [1.16_rp, 0.57_rp, 0.18_rp]
      atom%potential%width_band = 0.2_rp
      atom%potential%shifted_band = 0.0_rp
      atom%potential%obar = 0.0_rp
      call atom%build_pot()
      centers(:, :, 1) = atom%potential%center_band
      cx(:, :, 1) = atom%potential%cx
      expected = atom%potential%center_band(3, 2) - atom%potential%center_band(3, 1)
      call mills_1u_splitting_from_centers(centers, cx, delta, scalarity)
      call check('native center sign C_down-C_up', abs(delta(1) - expected) < 1.0e-14_rp, failed)
      call check('reversed sign negative control', abs(delta(1) + expected) > 0.1_rp, failed)
      call check('factor one-half negative control', abs(delta(1) - 0.5_rp*expected) > 0.05_rp, failed)
      call check('s center negative control', abs(delta(1) - (centers(1,2,1)-centers(1,1,1))) > 0.05_rp, failed)
      call check('p center negative control', abs(delta(1) - (centers(2,2,1)-centers(2,1,1))) > 0.05_rp, failed)
      do m = 5, 9
         call check('build_pot spherical d identity', &
            abs((atom%potential%cx(m,2)-atom%potential%cx(m,1))-cmplx(expected,0.0_rp,rp)) < 1.0e-14_rp, failed)
      end do
      call check('native d scalarity diagnostic', scalarity(1) < 1.0e-14_rp, failed)

      ! A complete H_down-H_up can contain unrelated terms; it is not an input
      ! to the direct center mapping and therefore cannot change native Delta.
      full_spin_difference = cmplx(0.0_rp, 0.0_rp, rp)
      do m = 5, 9
         full_spin_difference(m,m) = cmplx(expected, 0.0_rp, rp)
      end do
      full_spin_difference(5,5) = full_spin_difference(5,5) + cmplx(0.31_rp, 0.0_rp, rp)
      full_spin_difference(1,1) = cmplx(-0.22_rp, 0.0_rp, rp)
      full_spin_difference(5, 14) = cmplx(0.31_rp, 0.0_rp, rp)
      full_d_average = sum([(real(full_spin_difference(m,m),rp), m=5,9)])/5.0_rp
      call check('full H difference cannot redefine the native center splitting', &
         abs(delta(1)-expected) < 1.0e-14_rp .and. abs(delta(1)-full_d_average) > 0.05_rp .and. &
         abs(full_spin_difference(5,14)) > 0.0_rp, failed)
   end subroutine native_splitting_oracle

   subroutine coefficient_moment_oracle(failed)
      logical, intent(inout) :: failed
      type(lr_electronic_state) :: state
      real(rp) :: eigenvalues(2,1), kpoints(3,1), weights(1), occupations(2,1), moment(1), direct, radial_weighted
      complex(rp) :: vectors(18,2,1)
      integer :: band, orbital

      eigenvalues(:,1) = [-0.4_rp, 0.2_rp]
      kpoints = 0.0_rp
      weights = [1.0_rp]
      occupations(:,1) = [1.0_rp, 1.0_rp]
      vectors = cmplx(0.0_rp, 0.0_rp, rp)
      vectors(5,1,1) = cmplx(sqrt(0.64_rp), 0.0_rp, rp)
      vectors(1,1,1) = cmplx(sqrt(0.36_rp), 0.0_rp, rp)
      vectors(15,2,1) = cmplx(sqrt(0.25_rp), 0.0_rp, rp)
      vectors(12,2,1) = cmplx(sqrt(0.75_rp), 0.0_rp, rp)
      call state%initialize(eigenvalues, vectors, kpoints, weights, occupations, 0.0_rp, 300.0_rp, &
         reciprocal_mode='ham_only', hamiltonian_order='second', orthogonal=.true., collinear=.true., &
         has_soc=.false., has_extra_operator=.false.)
      direct = 0.0_rp
      do band = 1, 2
         do orbital = 5, 9
            direct = direct + weights(1)*occupations(band,1)* &
               (abs(vectors(orbital,band,1))**2-abs(vectors(orbital+9,band,1))**2)
         end do
      end do
      call mills_1u_coefficient_moment(state, 1, 2, moment)
      call check('direct coefficient-space d moment', abs(moment(1)-direct) < 1.0e-14_rp, failed)
      call check('moment fixture value', abs(moment(1)-0.39_rp) < 1.0e-14_rp, failed)
      radial_weighted = 0.64_rp*0.55_rp**2 - 0.25_rp*1.30_rp**2
      call check('non-unit radial norm rejects radial moment weighting', &
         abs(moment(1)-radial_weighted) > 0.1_rp, failed)
   end subroutine coefficient_moment_oracle

   subroutine transition_vertex_oracle(failed)
      logical, intent(inout) :: failed
      complex(rp) :: left(18), right(18), expected, vertex, radial_weighted
      integer :: m

      left = cmplx(0.0_rp, 0.0_rp, rp)
      right = cmplx(0.0_rp, 0.0_rp, rp)
      left(1) = cmplx(0.7_rp, 0.0_rp, rp)
      right(10) = cmplx(0.8_rp, 0.0_rp, rp)
      left(5) = cmplx(0.3_rp, 0.2_rp, rp)
      left(6) = cmplx(-0.1_rp, 0.4_rp, rp)
      left(9) = cmplx(0.2_rp, -0.3_rp, rp)
      right(14) = cmplx(0.5_rp, -0.1_rp, rp)
      right(15) = cmplx(0.2_rp, 0.6_rp, rp)
      right(18) = cmplx(-0.4_rp, 0.3_rp, rp)
      expected = cmplx(0.0_rp, 0.0_rp, rp)
      do m = 5, 9
         expected = expected + conjg(left(m))*right(m+9)
      end do
      vertex = mills_1u_transition_vertex(left, right, 1, 2, 1, 'chi_plus')
      call check('explicit local d transition vertex', abs(vertex-expected) < 1.0e-14_rp, failed)
      radial_weighted = conjg(left(5))*right(14)*0.55_rp + conjg(left(6))*right(15)*1.30_rp + &
         conjg(left(9))*right(18)*0.75_rp
      call check('transition vertex has no radial overlap', abs(vertex-radial_weighted) > 0.05_rp, failed)
      call check('s coefficient excluded from d vertex', &
         abs(vertex-(expected+conjg(left(1))*right(10))) > 0.1_rp, failed)
   end subroutine transition_vertex_oracle

   subroutine analytic_dyson_oracle(failed)
      logical, intent(inout) :: failed
      type(lr_electronic_state), target :: left_state, right_state
      type(mills_1u_bare_request) :: request, reverse_request
      type(mills_1u_bare_result) :: bare, reverse_bare
      type(projected_dyson_result) :: response, reverse_response
      real(rp) :: eigenvalues(2,1), kpoints(3,1), weights(1), occupations(2,1), delta, physical_u
      real(rp), parameter :: omega = 0.16_rp, eta = 0.015_rp
      complex(rp) :: vectors(18,2,1), transition, chi0_direct, denominator_direct, response_direct
      complex(rp) :: denominator_wrong_sign, denominator_wrong_factor

      delta = 0.7_rp
      physical_u = delta
      eigenvalues(:,1) = [-0.5_rp*delta, 0.5_rp*delta]
      kpoints = 0.0_rp
      weights = [1.0_rp]
      occupations(:,1) = [1.0_rp, 0.0_rp]
      vectors = cmplx(0.0_rp, 0.0_rp, rp)
      vectors(5,1,1) = cmplx(1.0_rp, 0.0_rp, rp)
      vectors(14,2,1) = cmplx(1.0_rp, 0.0_rp, rp)
      call left_state%initialize(eigenvalues, vectors, kpoints, weights, occupations, 0.0_rp, 0.0_rp, &
         reciprocal_mode='ham_only', hamiltonian_order='second', orthogonal=.true., collinear=.true., &
         has_soc=.false., has_extra_operator=.false.)
      call right_state%initialize(eigenvalues, vectors, kpoints, weights, occupations, 0.0_rp, 0.0_rp, &
         reciprocal_mode='ham_only', hamiltonian_order='second', orthogonal=.true., collinear=.true., &
         has_soc=.false., has_extra_operator=.false.)
      request%q = 0.0_rp
      request%frequencies = [omega]
      request%eta = eta
      request%channel = 'chi_plus'
      request%nsite = 1
      request%orbital_lmax = 2
      request%electronic_state => left_state
      request%q_endpoint_state => right_state
      call evaluate_mills_1u_bare(request, bare)
      transition = cmplx(1.0_rp, 0.0_rp, rp)
      chi0_direct = 2.0_rp*(occupations(1,1)-occupations(2,1))*transition*conjg(transition) / &
         cmplx(omega+eigenvalues(1,1)-eigenvalues(2,1), eta, rp)
      call check('analytic one-site Mills chi0', abs(bare%susceptibility(1,1,1)-chi0_direct) < 1.0e-13_rp, failed)
      call evaluate_mills_1u_dyson(bare, [physical_u], response)
      denominator_direct = 1.0_rp + physical_u / &
         cmplx(omega+eigenvalues(1,1)-eigenvalues(2,1), eta, rp)
      response_direct = chi0_direct/denominator_direct
      call check('physical RPA denominator 1+U chi_Kubo', &
         abs(response%denominator(1,1,1)-denominator_direct) < 1.0e-13_rp, failed)
      call check('analytic Mills Dyson response', abs(response%enhanced_chi(1,1,1)-response_direct) < 1.0e-13_rp, failed)
      call check('reported interaction remains physical U', abs(response%interaction_U(1)-physical_u) < 1.0e-14_rp, failed)
      denominator_wrong_sign = 1.0_rp - chi0_direct*(physical_u/2.0_rp)
      denominator_wrong_factor = 1.0_rp + chi0_direct*physical_u
      call check('wrong interaction sign negative control', &
         abs(response%denominator(1,1,1)-denominator_wrong_sign) > 0.1_rp, failed)
      call check('missing circular factor negative control', &
         abs(response%denominator(1,1,1)-denominator_wrong_factor) > 0.1_rp, failed)

      reverse_request = request
      reverse_request%channel = 'chi_minus'
      reverse_request%frequencies = [-omega]
      call evaluate_mills_1u_bare(reverse_request, reverse_bare)
      call evaluate_mills_1u_dyson(reverse_bare, [physical_u], reverse_response)
      call check('q/-q channel covariance at Gamma', &
         abs(reverse_response%enhanced_chi(1,1,1)-conjg(response%enhanced_chi(1,1,1))) < 1.0e-13_rp, failed)
   end subroutine analytic_dyson_oracle

   subroutine model_reduction_oracle(failed)
      logical, intent(inout) :: failed
      type(mills_1u_model_reduction_result) :: reduction
      complex(rp) :: difference(18,18)
      real(rp) :: delta(2), moment(2), physical_u(2)
      integer :: site, orbital, index

      delta = [0.30_rp, 0.28_rp]
      moment = [0.50_rp, 0.40_rp]
      physical_u = delta/moment
      difference = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, 2
         do orbital = 5, 9
            index = (site-1)*9 + orbital
            difference(index,index) = cmplx(delta(site), 0.0_rp, rp)
         end do
      end do
      call evaluate_mills_1u_model_reduction(difference, delta, 2, 9, reduction)
      call check('exact one-U model has zero reduction residual', reduction%residual_norm < 1.0e-14_rp, failed)
      call check('exact one-U model retains local d norm', abs(reduction%mills_local_d_norm- &
         sqrt(5.0_rp*sum(delta**2))) < 1.0e-14_rp, failed)

      difference(5, 14) = cmplx(0.17_rp, 0.0_rp, rp)
      difference(14, 5) = cmplx(0.17_rp, 0.0_rp, rp)
      call evaluate_mills_1u_model_reduction(difference, delta, 2, 9, reduction)
      call check('spin-dependent hopping increases residual', reduction%residual_norm > 0.2_rp, failed)
      call check('hopping appears in nonlocal component', reduction%site_offdiagonal_norm > 0.2_rp, failed)
      call check('hopping does not refit native U', maxval(abs(delta/moment-physical_u)) < 1.0e-14_rp, failed)

      difference(1,1) = cmplx(0.09_rp, 0.0_rp, rp)
      call evaluate_mills_1u_model_reduction(difference, delta, 2, 9, reduction)
      call check('s/p splitting appears in non-d component', reduction%non_d_sector_norm > 0.08_rp, failed)
      call check('s/p splitting does not refit native U', maxval(abs(delta/moment-physical_u)) < 1.0e-14_rp, failed)
   end subroutine model_reduction_oracle

   subroutine check(label, condition, failed)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      logical, intent(inout) :: failed
      if (.not. condition) then
         write(*, '(a)') '  FAIL: '//trim(label)
         failed = .true.
      end if
   end subroutine check

end module test_mills_1u_mod
