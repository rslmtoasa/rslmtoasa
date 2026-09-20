! DRESP-09V finite-dimensional production moment tangent fixture.
!
! The fixture is deliberately non-diagonal in orbital space and has unequal
! spin blocks.  It checks the complete endpoint product rule, the independent
! divided-difference Frechet route, rigid commutators, finite-angle theta^2
! convergence, the near-degenerate limit, production s/p block accumulation,
! and the spin-degenerate zero-response control.
program test_lr_lmto_density_moment_tangent

   use precision_mod, only: rp
   use spin_density_mod, only: spin_density
   use lr_exact_ks_ward_mod, only: exact_ks_spin_rotation_generator
   use lr_lmto_density_moment_tangent_mod, only: density_from_eigensystem, moment_matrices_from_eigensystem, &
      moment_tangent_product_rule, moment_tangent_frechet, endpoint_tangent_branches, endpoint_matrices_from_moments, &
      commutator_tangent, fill_production_spin_density_from_moments, relative_matrix_residual
   implicit none

   integer, parameter :: norb = 3, n = 2*norb, nl = 4
   real(rp), parameter :: fermi = 0.03_rp, temperature = 300.0_rp, tol = 3.0e-10_rp
   complex(rp) :: h(n,n), h0(n,n), eigenvectors(n,n), generator(n,n), dh(n,n), rho(n,n), drho(n,n)
   complex(rp) :: moments(n,n,3), dm_product(n,n,3), dm_frechet(n,n,3), dbranches(n,n,4)
   complex(rp) :: moments_plus(n,n,3), moments_minus(n,n,3), endpoints_plus(n,n,4), endpoints_minus(n,n,4)
   complex(rp) :: rotation(n,n), hrot(n,n), rhot(n,n), endpoints(n,n,4), fd(n,n,4)
   complex(rp) :: near_h(n,n), near_dh(n,n), near_v(n,n), near_plus(n,n), near_minus(n,n)
   complex(rp) :: near_tangent(n,n,3), near_fd(n,n,3)
   real(rp) :: evals(n), near_evals(n), near_plus_evals(n), near_minus_evals(n)
   real(rp) :: product_error(3), commutator_error(3), endpoint_error(4), finite_error(3), finite_ratio(2)
   real(rp) :: near_error(3), mapping_error, nonmagnetic_error, theta, tolerance_near
   integer :: i, j, p, branch, info
   logical :: failed
   type(spin_density) :: mapped
   external :: zheev

   call build_fixture(h)
   eigenvectors = h
   call diagonalize(eigenvectors, evals, info)
   if (info /= 0) error stop 'DRESP-09V unit: fixture eigensystem failed'
   call density_from_eigensystem(evals, eigenvectors, fermi, temperature, rho)
   call exact_ks_spin_rotation_generator(norb, 1, generator)
   call commutator_tangent(generator, h, dh)
   call commutator_tangent(generator, rho, drho)
   call moment_matrices_from_eigensystem(evals, eigenvectors, fermi, temperature, moments)
   call moment_tangent_product_rule(h, rho, dh, drho, dm_product)
   call moment_tangent_frechet(evals, eigenvectors, dh, fermi, temperature, dm_frechet)
   call endpoint_tangent_branches(h, rho, dh, drho, dbranches)

   do p = 1, 3
      product_error(p) = relative_matrix_residual(dm_product(:,:,p), dm_frechet(:,:,p))
      call commutator_tangent(generator, moments(:,:,p), rotation)
      commutator_error(p) = relative_matrix_residual(dm_product(:,:,p), rotation)
   end do
   call direct_endpoints(h, rho, endpoints)
   endpoint_error = 0.0_rp
   endpoint_error(1) = relative_matrix_residual(dbranches(:,:,1), dm_product(:,:,1))
   endpoint_error(2) = relative_matrix_residual(dbranches(:,:,2), dm_product(:,:,2))
   endpoint_error(3) = relative_matrix_residual(dbranches(:,:,3), dm_product(:,:,2))
   endpoint_error(4) = relative_matrix_residual(dbranches(:,:,4), dm_product(:,:,3))

   finite_error = 0.0_rp
   do i = 1, 3
      theta = 1.0e-2_rp/2.0_rp**real(i-1,rp)
      call build_rotation(theta, generator, rotation)
      hrot = matmul(rotation, matmul(h, transpose(conjg(rotation))))
      rhot = matmul(rotation, matmul(rho, transpose(conjg(rotation))))
      call direct_endpoints(hrot, rhot, endpoints_plus)
      hrot = matmul(transpose(conjg(rotation)), matmul(h, rotation))
      rhot = matmul(transpose(conjg(rotation)), matmul(rho, rotation))
      call direct_endpoints(hrot, rhot, endpoints_minus)
      fd = (endpoints_plus-endpoints_minus)/(2.0_rp*theta)
      finite_error(i) = maxval([relative_matrix_residual(fd(:,:,1), dbranches(:,:,1)), &
         relative_matrix_residual(fd(:,:,2), dbranches(:,:,2)), relative_matrix_residual(fd(:,:,3), dbranches(:,:,3)), &
         relative_matrix_residual(fd(:,:,4), dbranches(:,:,4))])
   end do
   finite_ratio = [finite_error(1)/max(finite_error(2),tiny(1.0_rp)), &
      finite_error(2)/max(finite_error(3),tiny(1.0_rp))]

   ! Near-degenerate divided differences are compared with a direct finite H
   ! derivative.  The two first levels differ by 1e-13 Ry.
   near_h = cmplx(0.0_rp, 0.0_rp, rp)
   near_h(1,1) = 0.1200000000000_rp
   near_h(2,2) = 0.1200000000001_rp
   near_h(3,3) = -0.11_rp
   near_h(4,4) = 0.21_rp
   near_h(5,5) = -0.07_rp
   near_h(6,6) = 0.31_rp
   near_h(1,2) = cmplx(0.013_rp, 0.004_rp, rp); near_h(2,1) = conjg(near_h(1,2))
   near_h(3,4) = cmplx(-0.009_rp, 0.003_rp, rp); near_h(4,3) = conjg(near_h(3,4))
   near_dh = cmplx(0.0_rp, 0.0_rp, rp)
   near_dh(1,2) = cmplx(0.017_rp, -0.006_rp, rp); near_dh(2,1) = conjg(near_dh(1,2))
   near_dh(2,3) = cmplx(-0.008_rp, 0.005_rp, rp); near_dh(3,2) = conjg(near_dh(2,3))
   near_v = near_h
   call diagonalize(near_v, near_evals, info)
   if (info /= 0) error stop 'DRESP-09V unit: near-degenerate diagonalization failed'
   call moment_tangent_frechet(near_evals, near_v, near_dh, fermi, temperature, near_tangent)
   near_plus = near_h + 1.0e-6_rp*near_dh
   near_minus = near_h - 1.0e-6_rp*near_dh
   call finite_moments(near_plus, near_plus_evals, moments_plus)
   call finite_moments(near_minus, near_minus_evals, moments_minus)
   near_fd = (moments_plus-moments_minus)/(2.0e-6_rp)
   do p = 1, 3
      near_error(p) = relative_matrix_residual(near_tangent(:,:,p), near_fd(:,:,p))
   end do
   tolerance_near = maxval(near_error)

   mapped = spin_density(1, nl)
   call fill_production_spin_density_from_moments(moments, mapped)
   mapping_error = maxval(abs(mapped%rho(1,1,:,:,:) - 0.0_rp))
   mapping_error = mapping_error/max(maxval(abs(mapped%rho)), tiny(1.0_rp))
   ! The block object must contain the direct s and p diagonal orbital sums.
   if (abs(mapped%rho(1,1,1,1,1) - moments(1,1,1)) > tol .or. &
       abs(mapped%rho(1,2,1,1,1) - moments(1,4,1)) > tol .or. &
       abs(mapped%rho(1,1,2,1,1) - moments(2,2,1) - moments(3,3,1)) > tol) then
      mapping_error = 1.0_rp
   else
      mapping_error = 0.0_rp
   end if

   h0 = cmplx(0.0_rp, 0.0_rp, rp)
   h0(1:norb,1:norb) = h(1:norb,1:norb)
   h0(norb+1:n,norb+1:n) = h(1:norb,1:norb)
   call commutator_tangent(generator, h0, dh)
   call diagonalize(h0, evals, info)
   call density_from_eigensystem(evals, h0, fermi, temperature, rho)
   call commutator_tangent(generator, rho, drho)
   nonmagnetic_error = max(sqrt(sum(abs(dh)**2)), sqrt(sum(abs(drho)**2)))

   failed = maxval(product_error) > tol .or. maxval(commutator_error) > tol .or. maxval(endpoint_error) > tol .or. &
      .not. (finite_error(1) > finite_error(2) .and. finite_error(2) > finite_error(3)) .or. &
      minval(finite_ratio) < 3.0_rp .or. tolerance_near > 2.0e-6_rp .or. mapping_error > tol .or. &
      nonmagnetic_error > tol
   if (failed) then
      write (*, '(a)') 'UnitLrLmtoDensityMomentTangent: FAIL'
      write (*, '(a,3es12.4)') '  product_vs_frechet=', product_error
      write (*, '(a,3es12.4)') '  commutator=', commutator_error
      write (*, '(a,4es12.4)') '  endpoint=', endpoint_error
      write (*, '(a,3es12.4)') '  finite=', finite_error
      write (*, '(a,2es12.4)') '  finite_theta2_ratio=', finite_ratio
      write (*, '(a,3es12.4)') '  near_degenerate=', near_error
      write (*, '(a,es12.4)') '  mapping=', mapping_error
      write (*, '(a,es12.4)') '  nonmagnetic_zero=', nonmagnetic_error
      error stop 1
   end if
   write (*, '(a,3es12.4)') '  product_vs_frechet=', product_error
   write (*, '(a,3es12.4)') '  commutator=', commutator_error
   write (*, '(a,4es12.4)') '  endpoint=', endpoint_error
   write (*, '(a,3es12.4)') '  finite=', finite_error
   write (*, '(a,3es12.4)') '  near_degenerate=', near_error
   write (*, '(a)') 'UnitLrLmtoDensityMomentTangent: PASS (M0/M1/M2, endpoints, production blocks, nonmagnetic control)'

contains

   subroutine build_fixture(matrix)
      complex(rp), intent(out) :: matrix(:,:)
      complex(rp) :: up(norb,norb), down(norb,norb)
      up = cmplx(0.0_rp,0.0_rp,rp); down = up
      do i = 1, norb
         up(i,i) = 0.07_rp + 0.031_rp*real(i,rp)
         down(i,i) = -0.12_rp + 0.047_rp*real(i,rp)
      end do
      do i = 1, norb
         do j = i+1, norb
            up(i,j) = cmplx(0.012_rp*real(i+j,rp), 0.003_rp*real(j-i,rp), rp)
            down(i,j) = cmplx(-0.009_rp*real(i+j,rp), 0.004_rp*real(j-i,rp), rp)
            up(j,i) = conjg(up(i,j)); down(j,i) = conjg(down(i,j))
         end do
      end do
      matrix = cmplx(0.0_rp,0.0_rp,rp)
      matrix(1:norb,1:norb) = up
      matrix(norb+1:n,norb+1:n) = down
   end subroutine build_fixture

   subroutine build_rotation(angle, gen, rot)
      real(rp), intent(in) :: angle
      complex(rp), intent(in) :: gen(:,:)
      complex(rp), intent(out) :: rot(:,:)
      integer :: ii
      rot = -cmplx(0.0_rp, 2.0_rp*sin(angle/2.0_rp), rp)*gen
      do ii = 1, size(rot,1)
         rot(ii,ii) = rot(ii,ii) + cmplx(cos(angle/2.0_rp),0.0_rp,rp)
      end do
   end subroutine build_rotation

   subroutine direct_endpoints(matrix, density_matrix, out)
      complex(rp), intent(in) :: matrix(:,:), density_matrix(:,:)
      complex(rp), intent(out) :: out(:,:,:)
      complex(rp) :: h2(size(matrix,1),size(matrix,2))
      out(:,:,1) = density_matrix
      out(:,:,2) = matmul(matrix,density_matrix)
      out(:,:,3) = matmul(density_matrix,matrix)
      h2 = matmul(matrix,matrix)
      out(:,:,4) = matmul(matrix,matmul(density_matrix,matrix))
   end subroutine direct_endpoints

   subroutine finite_moments(matrix, values, out)
      complex(rp), intent(in) :: matrix(:,:)
      real(rp), intent(out) :: values(:)
      complex(rp), intent(out) :: out(:,:,:)
      complex(rp) :: vec(size(matrix,1),size(matrix,2))
      integer :: inf
      vec = matrix
      call diagonalize(vec, values, inf)
      if (inf /= 0) error stop 'DRESP-09V unit: finite near-degenerate diagonalization failed'
      call moment_matrices_from_eigensystem(values, vec, fermi, temperature, out)
   end subroutine finite_moments

   subroutine diagonalize(matrix, values, inf)
      complex(rp), intent(inout) :: matrix(:,:)
      real(rp), intent(out) :: values(:)
      integer, intent(out) :: inf
      complex(rp) :: query(1)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      integer :: lwork
      allocate(rwork(max(1,3*size(values)-2)))
      call zheev('V','U',size(values),matrix,size(values),values,query,-1,rwork,inf)
      lwork = max(1,nint(real(query(1),rp)))
      allocate(work(lwork))
      call zheev('V','U',size(values),matrix,size(values),values,work,lwork,rwork,inf)
      deallocate(work,rwork)
   end subroutine diagonalize

end program test_lr_lmto_density_moment_tangent
