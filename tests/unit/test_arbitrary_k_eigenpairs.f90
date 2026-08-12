!------------------------------------------------------------------------------
! TDDFT-01 -- arbitrary-k reciprocal eigenpair service
!------------------------------------------------------------------------------
program test_arbitrary_k_eigenpairs
   use precision_mod, only: rp
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   implicit none

   type(reciprocal) :: recip
   type(hamiltonian), target :: ham
   type(lattice), target :: lat
   real(rp) :: k0(3), k_off(3)
   real(rp), allocatable :: evals(:, :), evals_qg(:, :), folded(:, :), evals_second(:, :), direct_evals(:)
   complex(rp), allocatable :: evecs(:, :, :), evecs_qg(:, :, :), evecs_second(:, :, :)
   complex(rp), allocatable :: h_direct(:, :), direct_vecs(:, :)
   logical :: failed

   call g_logger%init()
   failed = .false.
   call setup_one_site_model(recip, ham, lat)
   allocate(h_direct(nb, nb), direct_evals(nb), direct_vecs(nb, nb))

   k0 = [0.125_rp, 0.0_rp, 0.0_rp]
   k_off = [0.371_rp, -0.217_rp, 0.083_rp]

   ! q=0 must be the identity, including a meaningful gauge comparison for
   ! this non-degenerate diagonal model.
   call recip%calculate_eigenpairs_at_kpoints(reshape(k0, [3, 1]), evals, evecs, folded)
   call recip%calculate_eigenpairs_at_kpoints(reshape(k0 + 0.0_rp, [3, 1]), evals_qg, evecs_qg)
   call check_close('q=0 spectrum', maxval(abs(evals - evals_qg)), 1.0e-12_rp, failed)
   call check_gauge('q=0 eigenvectors', evecs(:, :, 1), evecs_qg(:, :, 1), failed)

   ! The service is caller-owned: it must not change the standard mesh or an
   ! already assembled mesh Hamiltonian as a side effect.
   if (recip%nk_total /= 1 .or. abs(recip%k_points(1, 1) - 0.25_rp) > 0.0_rp .or. &
       maxval(abs(recip%hk_bulk(:, :, 1) - cmplx(7.0_rp, 0.0_rp, rp))) > 0.0_rp) then
      write (*, '(a)') 'FAIL arbitrary-k service changed standard reciprocal mesh state'
      failed = .true.
   end if

   ! A q=0 response call is also the normal mesh eigensystem at the same k.
   recip%k_points(:, 1) = k0
   call recip%build_kspace_hamiltonian()
   call recip%diagonalize_hamiltonian()
   call check_close('q=0 normal-mesh spectrum', maxval(abs(evals(:, 1) - recip%eigenvalues(:, 1))), 1.0e-12_rp, failed)
   call check_gauge('q=0 normal-mesh eigenvectors', evecs(:, :, 1), recip%eigenvectors(:, :, 1), failed)

   ! k+q and k+q+G are the same physical point.  The second column also
   ! exercises exact duplicate folded-point caching in one service call.
   deallocate(evals, evecs, folded, evals_qg, evecs_qg)
   call recip%calculate_eigenpairs_at_kpoints(reshape([k0, k0 + [1.0_rp, 0.0_rp, 0.0_rp]], [3, 2]), &
                                                evals, evecs, folded)
   call check_close('q+G spectrum', maxval(abs(evals(:, 1) - evals(:, 2))), 1.0e-12_rp, failed)
   call check_close('q+G folded point', maxval(abs(folded(:, 1) - folded(:, 2))), 0.0_rp, failed)
   call check_gauge('q+G eigenvectors', evecs(:, :, 1), evecs(:, :, 2), failed)

   ! An off-mesh service point agrees with the direct normal Fourier H(k).
   deallocate(evals, evecs, folded)
   call recip%calculate_eigenpairs_at_kpoints(reshape(k_off, [3, 1]), evals, evecs, folded)
   call recip%fourier_transform_hamiltonian(folded(:, 1), h_direct)
   call check_hermiticity('first-order off-mesh Hamiltonian', h_direct, failed)
   call check_eigenpairs('first-order off-mesh service eigenpairs', h_direct, evals(:, 1), evecs(:, :, 1), failed)
   call hermitian_eig(h_direct, direct_evals, direct_vecs)
   call check_close('off-mesh direct Hamiltonian spectrum', maxval(abs(evals(:, 1) - direct_evals)), 1.0e-12_rp, failed)
   call check_gauge('off-mesh direct Hamiltonian eigenvectors', evecs(:, :, 1), direct_vecs, failed)

   ! Switch the same completed normal Hamiltonian model to the second-order
   ! route with an onsite spin-orbit term.  The arbitrary-k route must inherit
   ! exactly the same assembly as the direct normal reciprocal transform.
   call enable_second_order_terms(recip, ham)
   deallocate(evals, evecs, folded)
   call recip%calculate_eigenpairs_at_kpoints(reshape(k_off, [3, 1]), evals_second, evecs_second, folded)
   call recip%fourier_transform_hamiltonian_second_order(folded(:, 1), h_direct)
   call check_hermiticity('second-order/SOC off-mesh Hamiltonian', h_direct, failed)
   call check_eigenpairs('second-order/SOC off-mesh service eigenpairs', h_direct, evals_second(:, 1), &
                         evecs_second(:, :, 1), failed)
   call hermitian_eig(h_direct, direct_evals, direct_vecs)
   call check_close('second-order/SOC direct spectrum', maxval(abs(evals_second(:, 1) - direct_evals)), 1.0e-12_rp, failed)
   call check_gauge('second-order/SOC eigenvectors', evecs_second(:, :, 1), direct_vecs, failed)

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine setup_one_site_model(obj, ham_obj, lat_obj)
      type(reciprocal), intent(out) :: obj
      type(hamiltonian), target, intent(out) :: ham_obj
      type(lattice), target, intent(out) :: lat_obj
      integer :: i

      lat_obj%nrec = 1
      lat_obj%ntype = 1
      lat_obj%nn_max = 3
      lat_obj%kk = 1
      allocate(lat_obj%ib(1), lat_obj%atlist(1), lat_obj%iz(1), lat_obj%nn(1, 3))
      lat_obj%ib = 1
      lat_obj%atlist = 1
      lat_obj%iz = 1
      lat_obj%nn(1, :) = [3, 1, 1]

      ham_obj%lattice => lat_obj
      ham_obj%hoh = .false.
      ham_obj%ccor_2c = .false.
      allocate(ham_obj%ee(nb, nb, 3, 1))
      ham_obj%ee = (0.0_rp, 0.0_rp)
      do i = 1, nb
         ham_obj%ee(i, i, 1, 1) = cmplx(0.05_rp*real(i, rp), 0.0_rp, rp)
         ham_obj%ee(i, i, 2, 1) = cmplx(0.01_rp*real(i, rp), 0.0_rp, rp)
         ham_obj%ee(i, i, 3, 1) = cmplx(0.01_rp*real(i, rp), 0.0_rp, rp)
      end do

      obj%hamiltonian => ham_obj
      obj%lattice => lat_obj
      obj%reciprocal_mode = 'ham_only'
      obj%kspace_ham_order = 'first'
      obj%max_orbs = nb
      obj%dos_method = 'tetrahedron'
      obj%cached_operator_generation = 0
      ham_obj%operator_generation = 0
      ham_obj%magnetic_representation = 'periodic_nc'
      ! Pre-populate the geometry cache.  The service must not need a mesh to
      ! use it, and this avoids invoking lattice cluster construction here.
      allocate(obj%ham_vec_type(3, 3, 1), obj%ham_vec_type_direct(3, 3, 1))
      obj%ham_vec_type = 0.0_rp
      obj%ham_vec_type_direct = 0.0_rp
      obj%ham_vec_type_direct(1, 2, 1) = 1.0_rp
      obj%ham_vec_type_direct(1, 3, 1) = -1.0_rp

      obj%nk_total = 1
      allocate(obj%k_points(3, 1), obj%hk_bulk(nb, nb, 1))
      obj%k_points(:, 1) = [0.25_rp, 0.0_rp, 0.0_rp]
      obj%hk_bulk = cmplx(7.0_rp, 0.0_rp, rp)
   end subroutine setup_one_site_model

   subroutine enable_second_order_terms(obj, ham_obj)
      type(reciprocal), intent(inout) :: obj
      type(hamiltonian), intent(inout) :: ham_obj
      integer :: i

      ham_obj%hoh = .true.
      obj%kspace_ham_order = 'second'
      allocate(ham_obj%eeo(nb, nb, 3, 1), ham_obj%enim(nb, nb, 1), ham_obj%lsham(nb, nb, 1))
      ham_obj%eeo = (0.0_rp, 0.0_rp)
      ham_obj%enim = (0.0_rp, 0.0_rp)
      ham_obj%lsham = (0.0_rp, 0.0_rp)
      do i = 1, nb
         ham_obj%lsham(i, i, 1) = cmplx(0.001_rp*real(i, rp), 0.0_rp, rp)
      end do
   end subroutine enable_second_order_terms

   subroutine hermitian_eig(h, w, vecs)
      complex(rp), intent(in) :: h(:, :)
      real(rp), intent(out) :: w(:)
      complex(rp), intent(out) :: vecs(:, :)
      complex(rp) :: work_query(1)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      integer :: n, info, lwork

      n = size(h, 1)
      vecs = h
      allocate(rwork(max(1, 3*n - 2)))
      call zheev('V', 'U', n, vecs, n, w, work_query, -1, rwork, info)
      lwork = max(1, int(real(work_query(1), rp)))
      allocate(work(lwork))
      call zheev('V', 'U', n, vecs, n, w, work, lwork, rwork, info)
      if (info /= 0) error stop 'test_arbitrary_k_eigenpairs: zheev failed'
      deallocate(work, rwork)
   end subroutine hermitian_eig

   subroutine check_close(label, value, tolerance, test_failed)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: value, tolerance
      logical, intent(inout) :: test_failed
      write (*, '(a,a,es12.4)') 'Test ', trim(label)//' max_err = ', value
      if (value > tolerance) test_failed = .true.
   end subroutine check_close

   subroutine check_gauge(label, left, right, test_failed)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: left(:, :), right(:, :)
      logical, intent(inout) :: test_failed
      integer :: iband
      real(rp) :: max_error

      max_error = 0.0_rp
      do iband = 1, size(left, 2)
         max_error = max(max_error, abs(abs(sum(conjg(left(:, iband))*right(:, iband))) - 1.0_rp))
      end do
      call check_close(label, max_error, 1.0e-12_rp, test_failed)
   end subroutine check_gauge

   subroutine check_hermiticity(label, h, test_failed)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: h(:, :)
      logical, intent(inout) :: test_failed

      call check_close(label, maxval(abs(h - transpose(conjg(h)))), 1.0e-12_rp, test_failed)
   end subroutine check_hermiticity

   subroutine check_eigenpairs(label, h, eigenvalues, eigenvectors, test_failed)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: h(:, :), eigenvectors(:, :)
      real(rp), intent(in) :: eigenvalues(:)
      logical, intent(inout) :: test_failed
      complex(rp), allocatable :: identity(:, :), residual(:, :), overlap(:, :)
      integer :: n

      n = size(eigenvalues)
      allocate(identity(n, n), residual(n, n), overlap(n, n))
      identity = diag_identity(n)
      residual = matmul(h, eigenvectors) - eigenvectors*spread(eigenvalues, 1, n)
      overlap = matmul(conjg(transpose(eigenvectors)), eigenvectors) - identity
      call check_close(trim(label)//' residual', maxval(abs(residual)), 1.0e-12_rp, test_failed)
      call check_close(trim(label)//' orthonormality', maxval(abs(overlap)), 1.0e-12_rp, test_failed)
      deallocate(identity, residual, overlap)
   end subroutine check_eigenpairs

   function diag_identity(n) result(identity)
      integer, intent(in) :: n
      complex(rp) :: identity(n, n)
      integer :: i

      identity = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         identity(i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      end do
   end function diag_identity

end program test_arbitrary_k_eigenpairs
