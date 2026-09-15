!------------------------------------------------------------------------------
! TDVK-07 compact Dyson representation and loss oracle.
!
! The point-side solve below is deliberately independent of tddft_dyson_mod.
! It uses a small canonical point-space matrix, an explicit metric projection,
! and a hand-written pivoted Gaussian solve.  The compact side uses the live
! Dyson service on the same physical operator expressed in an orthonormal
! transformed basis.  The fixture has two sites, four angular sectors per
! site, four active radial sectors, complex off-diagonal chiKS, and a local
! point interaction which becomes non-diagonal after compact projection.
!------------------------------------------------------------------------------
program test_tddft_compact_dyson_oracle
   use precision_mod, only: rp
   use math_mod, only: pi
   use lr_response_space_mod, only: response_space_layout
   use tddft_dyson_mod, only: tddft_dyson_request, tddft_dyson_result, evaluate_tddft_dyson, &
      tddft_loss_matrix, lr_dyson_route_direct_alsda
   implicit none

   integer, parameter :: nsite = 2, response_lmax = 1, npoint = 5
   integer, parameter :: nactive = nsite*(response_lmax + 1)**2*(npoint - 1)
   real(rp), parameter :: mesh_a = 0.04_rp, mesh_b = 0.08_rp
   real(rp), parameter :: tolerance = 3.0e-9_rp
   type(response_space_layout), target :: space
   type(tddft_dyson_request) :: request
   type(tddft_dyson_result) :: compact_result
   real(rp) :: radius(npoint), d_f, relative_d_f, d_inf
   real(rp) :: loss_d_f, loss_relative_d_f, loss_d_inf, compact_residual
   complex(rp), allocatable :: b_point(:, :), k_point(:, :), b_interacting(:, :), point_loss(:, :)
   complex(rp), allocatable :: compact_bare(:, :), compact_kernel(:, :), projected(:, :), projected_loss(:, :)
   complex(rp), allocatable :: wrong(:, :), denominator(:, :), identity(:, :), rhs(:, :)
   complex(rp), allocatable :: u(:, :), transformed_bare(:, :), transformed_kernel(:, :)
   integer :: ir, i, j, a, b, ia, ib
   logical :: failed

   failed = .false.
   radius(1) = 0.0_rp
   do ir = 2, npoint
      radius(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
   end do
   call space%initialize(nsite, response_lmax, radius, mesh_a, mesh_b, 1)
   if (space%active_dimension /= nactive) error stop 'compact oracle: active dimension changed unexpectedly'

   allocate(b_point(space%ndim, space%ndim), k_point(space%ndim, space%ndim), &
      b_interacting(space%ndim, space%ndim), point_loss(space%ndim, space%ndim), &
      u(nactive, nactive), compact_bare(nactive, nactive), compact_kernel(nactive, nactive), &
      projected(nactive, nactive), projected_loss(nactive, nactive), wrong(nactive, nactive), &
      identity(nactive, nactive), transformed_bare(nactive, nactive))
   u = cmplx(0.0_rp, 0.0_rp, rp)
   do a = 1, nactive
      do b = 1, nactive
         u(a, b) = cmplx(cos(2.0_rp*pi*real((a - 1)*(b - 1), rp)/real(nactive, rp)), &
            sin(2.0_rp*pi*real((a - 1)*(b - 1), rp)/real(nactive, rp)), rp)/sqrt(real(nactive, rp))
      end do
   end do

   ! The compact bare matrix is complex and non-Hermitian, with finite
   ! retarded imaginary parts and all sectors coupled.
   compact_bare = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, nactive
      do j = 1, nactive
         compact_bare(i, j) = cmplx(0.0015_rp*real(1 + mod(i + 2*j, 17), rp), &
            0.0004_rp*real(i - j, rp), rp)
      end do
      compact_bare(i, i) = compact_bare(i, i) + cmplx(0.035_rp + 0.0002_rp*real(i, rp), &
         -0.002_rp - 0.00003_rp*real(i, rp), rp)
   end do

   ! A local point operator.  Its compact projection is intentionally dense
   ! and non-diagonal because U mixes site/angular/radial sectors.
   k_point = cmplx(0.0_rp, 0.0_rp, rp)
   a = 0
   do i = 1, space%ndim
      if (space%metric_weights(i) <= 0.0_rp) cycle
      a = a + 1
      k_point(i, i) = cmplx(0.045_rp + 0.002_rp*sin(real(3*a, rp)), 0.0_rp, rp)
   end do

   ! Build the canonical point matrix from B_t=U*B_c*U^H.  B_t is the
   ! orthonormal transformed form of the canonical point operator B=A*W.
   transformed_bare = matmul(u, matmul(compact_bare, conjg(transpose(u))))
   b_point = cmplx(0.0_rp, 0.0_rp, rp)
   a = 0
   do i = 1, space%ndim
      if (space%metric_weights(i) <= 0.0_rp) cycle
      a = a + 1
      b = 0
      do j = 1, space%ndim
         if (space%metric_weights(j) <= 0.0_rp) cycle
         b = b + 1
         b_point(i, j) = transformed_bare(a, b)*sqrt(space%metric_weights(j)/space%metric_weights(i))
      end do
   end do

   ! Independent point-space legacy convention: (I-B*K) B_int=B.
   ! The origin is the first point in every sector.  Solve the complete
   ! point-space matrix, retaining its null-measure rows and columns.
   allocate(denominator(space%ndim, space%ndim), rhs(space%ndim, space%ndim))
   call set_identity(denominator)
   do i = 1, space%ndim
      do j = 1, space%ndim
         denominator(i, j) = denominator(i, j) - b_point(i, j)*k_point(j, j)
      end do
   end do
   rhs = b_point
   call independent_solve(denominator, rhs)
   b_interacting = rhs

   ! Explicit point -> compact canonical projection: U^H W^(1/2) B
   ! W^(-1/2) U.  This is not a call to the production projection helper.
   projected = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, nactive
      ia = active_index(space, i)
      do j = 1, nactive
         ib = active_index(space, j)
         do a = 1, nactive
            do b = 1, nactive
               projected(i, j) = projected(i, j) + conjg(u(a, i))*sqrt(space%metric_weights(active_index(space, a)))* &
                  b_interacting(active_index(space, a), active_index(space, b))/sqrt(space%metric_weights(active_index(space, b)))*u(b, j)
            end do
         end do
      end do
   end do

   compact_kernel = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, nactive
      do j = 1, nactive
         do a = 1, nactive
            compact_kernel(i, j) = compact_kernel(i, j) + conjg(u(a, i))*k_point(active_index(space, a), active_index(space, a))*u(a, j)
         end do
      end do
   end do

   request%response_space => space
   request%compact_orthonormal = .true.
   request%frequencies = [0.17_rp]
   request%eta = 0.004_rp
   request%channel = 'chi_plus'
   request%ks_susceptibility = reshape(compact_bare, [nactive, nactive, 1])
   request%canonical_interaction = compact_kernel
   request%interaction_route = lr_dyson_route_direct_alsda
   request%interaction_provenance = 'KXC-01 direct ALSDA TDVK-07 independent point oracle'
   request%electronic_state_provenance = 'TDVK-07 synthetic accepted compact fixture'
   call evaluate_tddft_dyson(request, compact_result)

   d_f = sqrt(sum(abs(projected - compact_result%enhanced_susceptibility(:, :, 1))**2))
   relative_d_f = d_f/max(sqrt(sum(abs(projected)**2)), &
      sqrt(sum(abs(compact_result%enhanced_susceptibility(:, :, 1))**2)), tiny(1.0_rp))
   d_inf = maxval(abs(projected - compact_result%enhanced_susceptibility(:, :, 1)))

   ! Independent point loss uses the W-adjoint; its compact projection must
   ! equal the ordinary Hermitian loss in the orthonormal product basis.
   point_loss = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, space%ndim
      do j = 1, space%ndim
         if (space%metric_weights(i) > 0.0_rp .and. space%metric_weights(j) > 0.0_rp) then
            point_loss(i, j) = -(b_interacting(i, j) - conjg(b_interacting(j, i))* &
               space%metric_weights(j)/space%metric_weights(i))/cmplx(0.0_rp, 2.0_rp*pi, rp)
         end if
      end do
   end do
   projected_loss = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, nactive
      ia = active_index(space, i)
      do j = 1, nactive
         ib = active_index(space, j)
         do a = 1, nactive
            do b = 1, nactive
               projected_loss(i, j) = projected_loss(i, j) + conjg(u(a, i))*sqrt(space%metric_weights(active_index(space, a)))* &
                  point_loss(active_index(space, a), active_index(space, b))/sqrt(space%metric_weights(active_index(space, b)))*u(b, j)
            end do
         end do
      end do
   end do
   loss_d_f = sqrt(sum(abs(projected_loss - compact_result%loss_matrix(:, :, 1))**2))
   loss_relative_d_f = loss_d_f/max(sqrt(sum(abs(projected_loss)**2)), &
      sqrt(sum(abs(compact_result%loss_matrix(:, :, 1))**2)), tiny(1.0_rp))
   loss_d_inf = maxval(abs(projected_loss - compact_result%loss_matrix(:, :, 1)))
   compact_residual = compact_result%dyson_residual_frobenius(1)

   call check('fixture has complex off-diagonal chiKS', maxval(abs(aimag(compact_bare))) > 1.0e-5_rp .and. &
      maxval(abs(compact_bare - transpose(compact_bare))) > 1.0e-4_rp, failed)
   call check('local point interaction projects to non-diagonal compact operator', &
      maxval(abs(compact_kernel - diagonal_part(compact_kernel))) > 1.0e-5_rp, failed)
   call check('point-space oracle solve succeeds', compact_result%solve_succeeded(1), failed)
   call check('complete projected Dyson matrix agrees', d_f < tolerance .and. relative_d_f < tolerance .and. d_inf < tolerance, failed)
   call check('loss oracle agrees', loss_d_f < tolerance .and. loss_relative_d_f < tolerance .and. loss_d_inf < tolerance, failed)
   call check('Dyson algebraic residual is controlled', compact_residual < tolerance .and. &
      compact_result%dyson_residual_relative(1) < tolerance .and. compact_result%dyson_residual_infinity(1) < tolerance, failed)

   call solve_compact_variant(compact_bare, compact_kernel, .false., .false., wrong)
   call check('wrong compact matrix ordering is detectable', maxval(abs(wrong - compact_result%enhanced_susceptibility(:, :, 1))) > 1.0e-6_rp, failed)
   call solve_compact_variant(compact_bare, compact_kernel, .true., .false., wrong)
   call check('wrong Dyson sign is detectable', maxval(abs(wrong - compact_result%enhanced_susceptibility(:, :, 1))) > 1.0e-6_rp, failed)
   call solve_compact_variant(compact_bare, compact_kernel, .false., .true., wrong)
   call check('wrong interaction projection/conjugation is detectable', maxval(abs(wrong - compact_result%enhanced_susceptibility(:, :, 1))) > 1.0e-6_rp, failed)

   write (*, '(a,i0)') 'TDVK-07 compact Dyson point-space oracle dimension = ', nactive
   write (*, '(a,es24.16)') '  point projected chi Frobenius = ', sqrt(sum(abs(projected)**2))
   write (*, '(a,es24.16)') '  compact chi Frobenius = ', sqrt(sum(abs(compact_result%enhanced_susceptibility(:, :, 1))**2))
   write (*, '(a,3(es24.16,1x))') '  oracle dF relative_dF dInf = ', d_f, relative_d_f, d_inf
   write (*, '(a,3(es24.16,1x))') '  loss dF relative_dF dInf = ', loss_d_f, loss_relative_d_f, loss_d_inf
   write (*, '(a,3(es24.16,1x))') '  Dyson residual Frobenius relative dInf = ', compact_residual, &
      compact_result%dyson_residual_relative(1), compact_result%dyson_residual_infinity(1)
   write (*, '(a,3(es24.16,1x))') '  denominator min_sv max_sv condition = ', compact_result%denominator_min_singular_value(1), &
      compact_result%denominator_max_singular_value(1), compact_result%denominator_condition_number(1)
   if (failed) then
      write (*, '(a)') 'UnitTddftCompactDysonOracle: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'UnitTddftCompactDysonOracle: PASS (independent point Dyson, compact projection, residual, loss, complex ordering guards)'

contains

   integer function active_index(layout, active) result(index_value)
      type(response_space_layout), intent(in) :: layout
      integer, intent(in) :: active
      integer :: flat, count_active
      count_active = 0
      index_value = 0
      do flat = 1, layout%ndim
         if (layout%metric_weights(flat) <= 0.0_rp) cycle
         count_active = count_active + 1
         if (count_active == active) then
            index_value = flat
            return
         end if
      end do
      error stop 'compact oracle: active index is outside the response space'
   end function active_index

   subroutine set_identity(matrix)
      complex(rp), intent(out) :: matrix(:, :)
      integer :: index
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do index = 1, size(matrix, 1)
         matrix(index, index) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end subroutine set_identity

   function diagonal_part(matrix) result(diagonal)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp) :: diagonal(size(matrix, 1), size(matrix, 2))
      integer :: index
      diagonal = cmplx(0.0_rp, 0.0_rp, rp)
      do index = 1, size(matrix, 1)
         diagonal(index, index) = matrix(index, index)
      end do
   end function diagonal_part

   subroutine independent_solve(matrix, rhs_matrix)
      complex(rp), intent(in) :: matrix(:, :)
      complex(rp), intent(inout) :: rhs_matrix(:, :)
      complex(rp), allocatable :: work(:, :), row(:), rhs_row(:)
      complex(rp) :: factor
      integer :: n, nrhs_local, row_index, column_index, rhs_index, pivot_index

      n = size(matrix, 1)
      nrhs_local = size(rhs_matrix, 2)
      allocate(work(n, n), row(n), rhs_row(nrhs_local))
      work = matrix
      do column_index = 1, n - 1
         pivot_index = column_index
         do row_index = column_index + 1, n
            if (abs(work(row_index, column_index)) > abs(work(pivot_index, column_index))) pivot_index = row_index
         end do
         if (abs(work(pivot_index, column_index)) <= tiny(1.0_rp)) error stop 'compact oracle: singular independent point solve'
         if (pivot_index /= column_index) then
            row = work(column_index, :); work(column_index, :) = work(pivot_index, :); work(pivot_index, :) = row
            rhs_row = rhs_matrix(column_index, :); rhs_matrix(column_index, :) = rhs_matrix(pivot_index, :); &
               rhs_matrix(pivot_index, :) = rhs_row
         end if
         do row_index = column_index + 1, n
            factor = work(row_index, column_index)/work(column_index, column_index)
            work(row_index, column_index:n) = work(row_index, column_index:n) - factor*work(column_index, column_index:n)
            rhs_matrix(row_index, :) = rhs_matrix(row_index, :) - factor*rhs_matrix(column_index, :)
         end do
      end do
      if (abs(work(n, n)) <= tiny(1.0_rp)) error stop 'compact oracle: singular independent back solve'
      do column_index = n, 1, -1
         rhs_matrix(column_index, :) = rhs_matrix(column_index, :)/work(column_index, column_index)
         do row_index = 1, column_index - 1
            rhs_matrix(row_index, :) = rhs_matrix(row_index, :) - work(row_index, column_index)*rhs_matrix(column_index, :)
         end do
      end do
      deallocate(work, row, rhs_row)
   end subroutine independent_solve

   subroutine solve_compact_variant(bare, kernel, wrong_sign, wrong_projection, solution)
      complex(rp), intent(in) :: bare(:, :), kernel(:, :)
      logical, intent(in) :: wrong_sign, wrong_projection
      complex(rp), intent(out) :: solution(:, :)
      complex(rp), allocatable :: d(:, :), rhs_variant(:, :), effective_kernel(:, :)

      allocate(d(size(bare, 1), size(bare, 2)), rhs_variant(size(bare, 1), size(bare, 2)), &
         effective_kernel(size(kernel, 1), size(kernel, 2)))
      effective_kernel = kernel
      ! A transposed (conjugation-dropped) projection is a deliberately wrong
      ! compact interaction for this complex unitary fixture.
      if (wrong_projection) effective_kernel = transpose(kernel)
      call set_identity(d)
      if (wrong_sign) then
         d = d + matmul(bare, effective_kernel)
      else if (wrong_projection) then
         d = d - matmul(bare, effective_kernel)
      else
         ! Deliberately reverse the operator order for this guard.
         d = d - matmul(effective_kernel, bare)
      end if
      rhs_variant = bare
      call independent_solve(d, rhs_variant)
      solution = rhs_variant
      deallocate(d, rhs_variant, effective_kernel)
   end subroutine solve_compact_variant

   subroutine check(label, condition, test_failed)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      logical, intent(inout) :: test_failed
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL', trim(label)
         test_failed = .true.
      end if
   end subroutine check

end program test_tddft_compact_dyson_oracle
