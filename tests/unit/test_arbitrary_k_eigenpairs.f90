!------------------------------------------------------------------------------
! TDDFT-01 -- arbitrary-k reciprocal eigenpair service
!------------------------------------------------------------------------------
program test_arbitrary_k_eigenpairs
   use precision_mod, only: rp
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal, reciprocal_execution_capabilities, reciprocal_execution_request, &
                             reciprocal_execution_result, reciprocal_residency_host, lapack_reciprocal_backend
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   implicit none

   type(reciprocal) :: recip
   type(reciprocal) :: recip_two_site
   type(reciprocal) :: recip_empty
   type(hamiltonian), target :: ham
   type(hamiltonian), target :: ham_two_site
   type(hamiltonian), target :: ham_empty
   type(lattice), target :: lat
   type(lattice), target :: lat_two_site
   type(lattice), target :: lat_empty
   real(rp) :: k0(3), k_off(3)
   real(rp), allocatable :: evals(:, :), evals_qg(:, :), folded(:, :), evals_second(:, :), direct_evals(:)
   real(rp), allocatable :: evals_tile1(:, :), evals_tile2(:, :)
   complex(rp), allocatable :: evecs(:, :, :), evecs_qg(:, :, :), evecs_second(:, :, :)
   complex(rp), allocatable :: evecs_tile1(:, :, :), evecs_tile2(:, :, :)
   complex(rp), allocatable :: h_direct(:, :), direct_vecs(:, :)
   logical :: failed
   integer :: workspace_capacity
   integer :: allocation_count, query_count
   integer :: execute_before, combined_before, assemble_before, input_solve_before
   integer :: execute_after, combined_after, assemble_after, input_solve_after
   type(reciprocal_execution_capabilities) :: backend_caps
   type(reciprocal_execution_request) :: backend_request
   type(reciprocal_execution_result) :: backend_result

   call g_logger%init()
   failed = .false.
   call setup_one_site_model(recip, ham, lat)
   allocate(h_direct(nb, nb), direct_evals(nb), direct_vecs(nb, nb))

   k0 = [0.125_rp, 0.0_rp, 0.0_rp]
   k_off = [0.371_rp, -0.217_rp, 0.083_rp]

   ! RF-05 contract: the default factory is CPU/LAPACK, advertises the full
   ! current solve surface, and can return an assembled diagnostic matrix
   ! while keeping an eigenvalues-only tile free of eigenvector ownership.
   call recip%make_execution_backend('lapack')
   call recip%execution_backend%capabilities(backend_caps)
   if (.not. backend_caps%standard_hermitian .or. .not. backend_caps%generalized_hermitian .or. &
       .not. backend_caps%eigenvalues_only .or. .not. backend_caps%eigenvectors .or. &
       .not. backend_caps%first_order_assembly .or. .not. backend_caps%second_order_assembly .or. &
       .not. backend_caps%overlap .or. backend_caps%residency /= reciprocal_residency_host) then
      write (*, '(a)') 'FAIL RF-05 LAPACK backend capability contract'
      failed = .true.
   end if
   allocate(backend_request%k_points(3,1))
   backend_request%k_points(:,1) = k0
   backend_request%solve_eigensystem = .true.
   backend_request%request_eigenvectors = .false.
   backend_request%request_assembled_hamiltonian = .true.
   backend_request%operator_generation = ham%operator_generation
   call recip%execution_backend%execute_batch(backend_request, backend_result)
   if (.not. backend_result%eigenvalues_valid .or. backend_result%eigenvectors_valid .or. &
       .not. backend_result%assembled_hamiltonian_valid .or. allocated(backend_result%eigenvectors) .or. &
       backend_result%local_point_count /= 1 .or. backend_result%operator_generation /= ham%operator_generation) then
      write (*, '(a)') 'FAIL RF-05 eigenvalues-only result ownership contract'
      failed = .true.
   end if
   call check_hermiticity('RF-05 assembled backend H(k)', backend_result%hamiltonian(:,:,1), failed)
   call recip%execution_backend%release()
   deallocate(recip%execution_backend)
   call recip%make_execution_backend('cpu')

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
   recip%reciprocal_tile_size = 1
   call recip%execution_backend%execution_metrics(execute_before, combined_before, assemble_before, input_solve_before)
   call recip%build_kspace_hamiltonian()
   call recip%execution_backend%execution_metrics(execute_after, combined_after, assemble_after, input_solve_after)
#ifdef RSLMTO_DISABLE_FUSED_RECIPROCAL
   if (execute_after /= execute_before + 1 .or. combined_after /= combined_before .or. &
       assemble_after /= assemble_before + 1 .or. input_solve_after /= input_solve_before) then
      write (*, '(a)') 'FAIL compatibility normal mesh did not use one assembly request'
      failed = .true.
   end if
#else
   if (execute_after /= execute_before + 1 .or. combined_after /= combined_before + 1 .or. &
       assemble_after /= assemble_before .or. input_solve_after /= input_solve_before) then
      write (*, '(a)') 'FAIL GC-04 normal mesh did not use one combined backend request'
      failed = .true.
   end if
#endif
   call recip%diagonalize_hamiltonian()
   call recip%execution_backend%execution_metrics(execute_after, combined_after, assemble_after, input_solve_after)
#ifdef RSLMTO_DISABLE_FUSED_RECIPROCAL
   if (execute_after /= execute_before + 2 .or. combined_after /= combined_before .or. &
       assemble_after /= assemble_before + 1 .or. input_solve_after /= input_solve_before + 1) then
      write (*, '(a)') 'FAIL compatibility diagonalize adapter did not submit one input-H solve'
      failed = .true.
   end if
#else
   if (execute_after /= execute_before + 1 .or. combined_after /= combined_before + 1 .or. &
       assemble_after /= assemble_before .or. input_solve_after /= input_solve_before) then
      write (*, '(a)') 'FAIL GC-04 diagonalize adapter submitted a second backend request'
      failed = .true.
   end if
#endif
   call recip%fourier_transform_hamiltonian(k0, h_direct)
   call check_close('GC-04 normal-mesh H(k) compatibility cache', maxval(abs(recip%hk_bulk(:,:,1) - h_direct)), 1.0e-12_rp, failed)
   call check_eigenpairs('GC-04 normal-mesh combined eigenpairs', recip%hk_bulk(:,:,1), recip%eigenvalues(:,1), &
                         recip%eigenvectors(:,:,1), failed)
   call check_close('q=0 normal-mesh spectrum', maxval(abs(evals(:, 1) - recip%eigenvalues(:, 1))), 1.0e-12_rp, failed)
   call check_gauge('q=0 normal-mesh eigenvectors', evecs(:, :, 1), recip%eigenvectors(:, :, 1), failed)

   ! A generation change must invalidate both compatibility caches and rebuild
   ! them through one new combined tile, not an input-H solve.
   ham%operator_generation = ham%operator_generation + 1
   call recip%diagonalize_hamiltonian()
   call recip%execution_backend%execution_metrics(execute_after, combined_after, assemble_after, input_solve_after)
#ifdef RSLMTO_DISABLE_FUSED_RECIPROCAL
   if (execute_after /= execute_before + 4 .or. combined_after /= combined_before .or. &
       assemble_after /= assemble_before + 2 .or. input_solve_after /= input_solve_before + 2 .or. &
       recip%cached_operator_generation /= ham%operator_generation) then
      write (*, '(a)') 'FAIL compatibility stale normal mesh did not rebuild and re-solve'
      failed = .true.
   end if
#else
   if (execute_after /= execute_before + 2 .or. combined_after /= combined_before + 2 .or. &
       assemble_after /= assemble_before .or. input_solve_after /= input_solve_before .or. &
       recip%cached_operator_generation /= ham%operator_generation) then
      write (*, '(a)') 'FAIL GC-04 stale normal-mesh generation did not rebuild as a combined tile'
      failed = .true.
   end if
#endif
   call check_close('GC-04 generation-refresh spectrum', maxval(abs(evals(:,1) - recip%eigenvalues(:,1))), 1.0e-12_rp, failed)

   ! A k-path follows the same resident tile contract while retaining the
   ! ordinary public build/diagonalize entry points.
   allocate(recip%k_path(3,2))
   recip%k_path(:,1) = k0
   recip%k_path(:,2) = k_off
   recip%nk_path = 2
   call recip%execution_backend%execution_metrics(execute_before, combined_before, assemble_before, input_solve_before)
   call recip%build_kspace_hamiltonian()
   call recip%diagonalize_hamiltonian()
   call recip%execution_backend%execution_metrics(execute_after, combined_after, assemble_after, input_solve_after)
#ifdef RSLMTO_DISABLE_FUSED_RECIPROCAL
   if (execute_after /= execute_before + 3 .or. combined_after /= combined_before .or. &
       assemble_after /= assemble_before + 2 .or. input_solve_after /= input_solve_before + 1) then
      write (*, '(a)') 'FAIL compatibility k-path did not use assembly plus input-H solve'
      failed = .true.
   end if
#else
   if (execute_after /= execute_before + 2 .or. combined_after /= combined_before + 2 .or. &
       assemble_after /= assemble_before .or. input_solve_after /= input_solve_before) then
      write (*, '(a)') 'FAIL GC-04 k-path did not use one combined request per tile'
      failed = .true.
   end if
#endif
   call check_eigenpairs('GC-04 k-path combined eigenpairs', recip%hk_bulk(:,:,2), recip%eigenvalues(:,2), &
                         recip%eigenvectors(:,:,2), failed)
   deallocate(recip%k_path)
   recip%nk_path = 0

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

   ! RF-04: one-point, partial, and multi-tile batches must agree while exact
   ! duplicates still scatter back into the caller's original order.
   recip%reciprocal_tile_size = 1
   call recip%calculate_eigenpairs_at_kpoints(reshape([k0, k_off, k0 + [1.0_rp, 0.0_rp, 0.0_rp], &
                                                         [0.0_rp, 0.125_rp, 0.0_rp]], [3, 4]), &
                                              evals_tile1, evecs_tile1)
   workspace_capacity = recip%workspace%tile_capacity
   recip%reciprocal_tile_size = 3
   call recip%calculate_eigenpairs_at_kpoints(reshape([k0, k_off, k0 + [1.0_rp, 0.0_rp, 0.0_rp], &
                                                         [0.0_rp, 0.125_rp, 0.0_rp]], [3, 4]), &
                                              evals_tile2, evecs_tile2)
   call check_close('RF-04 tiled eigenvalue equivalence', maxval(abs(evals_tile1-evals_tile2)), 1.0e-12_rp, failed)
   call check_gauge('RF-04 tiled eigenvector equivalence', evecs_tile1(:,:,2), evecs_tile2(:,:,2), failed)
   if (recip%workspace%tile_capacity < workspace_capacity) then
      write (*, '(a)') 'FAIL RF-04 workspace capacity shrank unexpectedly'
      failed = .true.
   end if
   ham%operator_generation = ham%operator_generation + 1
   call recip%calculate_eigenpairs_at_kpoints(reshape(k0, [3, 1]), evals_qg, evecs_qg)
   if (recip%workspace%cached_operator_generation /= ham%operator_generation) then
      write (*, '(a)') 'FAIL RF-04 workspace did not refresh operator-generation fingerprint'
      failed = .true.
   end if
   deallocate(evals_tile1, evecs_tile1, evals_tile2, evecs_tile2, evals_qg, evecs_qg)

   ! RF-04 generalized/multi-site coverage.  A two-site model with a positive
   ! overlap proxy verifies tile equivalence in the zhegv path and proves that
   ! a prepared repeat performs neither another workspace allocation nor a
   ! LAPACK workspace query.
   call setup_two_site_generalized_model(recip_two_site, ham_two_site, lat_two_site)
   recip_two_site%reciprocal_tile_size = 3
   call recip_two_site%calculate_eigenpairs_at_kpoints(reshape([k0, k_off, [0.0_rp, 0.125_rp, 0.0_rp], &
                                                                   [0.25_rp, -0.125_rp, 0.0_rp]], [3, 4]), &
                                                       evals_tile1, evecs_tile1)
   if (allocated(backend_request%k_points)) deallocate(backend_request%k_points)
   backend_request%assemble_hamiltonian = .true.
   backend_request%assemble_overlap = .true.
   backend_request%solve_eigensystem = .true.
   backend_request%generalized = .true.
   backend_request%request_eigenvectors = .true.
   backend_request%request_assembled_hamiltonian = .true.
   backend_request%request_assembled_overlap = .true.
   backend_request%operator_generation = ham_two_site%operator_generation
   allocate(backend_request%k_points(3,1)); backend_request%k_points(:,1) = k0
   call recip_two_site%execution_backend%execute_batch(backend_request, backend_result)
   call check_generalized_eigenpairs('GC-02 generalized residual and metric orthogonality', &
                                     backend_result%hamiltonian(:,:,1), backend_result%overlap(:,:,1), &
                                     backend_result%eigenvalues(:,1), backend_result%eigenvectors(:,:,1), failed)
   allocation_count = recip_two_site%workspace%storage_allocations
   query_count = recip_two_site%workspace%lapack_workspace_queries
   recip_two_site%reciprocal_tile_size = 2
   call recip_two_site%calculate_eigenpairs_at_kpoints(reshape([k0, k_off, [0.0_rp, 0.125_rp, 0.0_rp], &
                                                                   [0.25_rp, -0.125_rp, 0.0_rp]], [3, 4]), &
                                                       evals_tile2, evecs_tile2)
   call check_close('RF-04 multi-site generalized tiled eigenvalue equivalence', &
                    maxval(abs(evals_tile1-evals_tile2)), 1.0e-12_rp, failed)
   call check_gauge_normalized('RF-04 multi-site generalized tiled eigenvector equivalence', &
                               evecs_tile1(:,:,3), evecs_tile2(:,:,3), failed)
   if (recip_two_site%workspace%storage_allocations /= allocation_count .or. &
       recip_two_site%workspace%lapack_workspace_queries /= query_count) then
      write (*, '(a)') 'FAIL RF-04 prepared generalized tile repeated an allocation or LAPACK query'
      failed = .true.
   end if
   select type (backend => recip_two_site%execution_backend)
   type is (lapack_reciprocal_backend)
      if (.not. allocated(backend%workspace%overlap_cholesky) .or. &
          any(shape(backend%workspace%overlap_cholesky) /= [2*nb, 2*nb])) then
         write (*, '(a)') 'FAIL GC-02 generalized overlap validation scratch is not prepared'
         failed = .true.
      end if
   class default
      write (*, '(a)') 'FAIL GC-02 generalized test did not retain the LAPACK backend'
      failed = .true.
   end select

   ! GC-04 generalized normal-mesh tiles must assemble H/S and solve together.
   recip_two_site%nk_total = 2
   allocate(recip_two_site%k_points(3,2))
   recip_two_site%k_points(:,1) = k0
   recip_two_site%k_points(:,2) = k_off
   recip_two_site%reciprocal_tile_size = 1
   call recip_two_site%execution_backend%execution_metrics(execute_before, combined_before, assemble_before, input_solve_before)
   call recip_two_site%build_kspace_hamiltonian()
   call recip_two_site%diagonalize_hamiltonian()
   call recip_two_site%execution_backend%execution_metrics(execute_after, combined_after, assemble_after, input_solve_after)
#ifdef RSLMTO_DISABLE_FUSED_RECIPROCAL
   if (execute_after /= execute_before + 3 .or. combined_after /= combined_before .or. &
       assemble_after /= assemble_before + 2 .or. input_solve_after /= input_solve_before + 1) then
      write (*, '(a)') 'FAIL compatibility generalized normal mesh did not assemble H/S and solve'
      failed = .true.
   end if
#else
   if (execute_after /= execute_before + 2 .or. combined_after /= combined_before + 2 .or. &
       assemble_after /= assemble_before .or. input_solve_after /= input_solve_before) then
      write (*, '(a)') 'FAIL GC-04 generalized normal mesh did not use one combined request per tile'
      failed = .true.
   end if
#endif
   call check_generalized_eigenpairs('GC-04 generalized normal-mesh eigenpairs', recip_two_site%hk_bulk(:,:,1), &
                                     recip_two_site%sk_overlap(:,:,1), recip_two_site%eigenvalues(:,1), &
                                     recip_two_site%eigenvectors(:,:,1), failed)
   call check_close('GC-04 generalized normal-mesh spectrum', &
                    maxval(abs(recip_two_site%eigenvalues(:,1) - evals_tile1(:,1))), 1.0e-12_rp, failed)
   deallocate(evals_tile1, evecs_tile1, evals_tile2, evecs_tile2)

   ! Empty local ownership is a defined no-op: caches have zero local extent
   ! and no backend request is created merely to represent an empty tile.
   call setup_one_site_model(recip_empty, ham_empty, lat_empty)
   deallocate(recip_empty%k_points, recip_empty%hk_bulk)
   allocate(recip_empty%k_points(3,0))
   recip_empty%nk_total = 0
   recip_empty%kanpur_diagnostics = .false.
   call recip_empty%build_kspace_hamiltonian()
   call recip_empty%diagonalize_hamiltonian()
   if (size(recip_empty%hk_bulk,3) /= 0 .or. size(recip_empty%eigenvalues,2) /= 0 .or. &
       size(recip_empty%eigenvectors,3) /= 0 .or. allocated(recip_empty%execution_backend)) then
      write (*, '(a)') 'FAIL GC-04 empty local mesh did not remain a no-op cache state'
      failed = .true.
   end if

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

      call obj%restore_to_default()
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

   subroutine setup_two_site_generalized_model(obj, ham_obj, lat_obj)
      type(reciprocal), intent(out) :: obj
      type(hamiltonian), target, intent(out) :: ham_obj
      type(lattice), target, intent(out) :: lat_obj
      integer :: i

      call obj%restore_to_default()
      lat_obj%nrec = 2; lat_obj%ntype = 2; lat_obj%nn_max = 2; lat_obj%kk = 2
      allocate(lat_obj%ib(2), lat_obj%atlist(2), lat_obj%iz(2), lat_obj%nn(2,2))
      lat_obj%ib = [1,2]; lat_obj%atlist = [1,2]; lat_obj%iz = [1,2]
      lat_obj%nn(1,:) = [2,2]; lat_obj%nn(2,:) = [2,1]
      ham_obj%lattice => lat_obj; ham_obj%hoh = .false.; ham_obj%ccor_2c = .false.
      ham_obj%operator_generation = 0; ham_obj%magnetic_representation = 'periodic_nc'
      allocate(ham_obj%ee(nb,nb,2,2), ham_obj%eeo(nb,nb,2,2))
      ham_obj%ee = cmplx(0.0_rp,0.0_rp,rp); ham_obj%eeo = cmplx(0.0_rp,0.0_rp,rp)
      do i=1,nb
         ham_obj%ee(i,i,1,1) = cmplx(0.04_rp*real(i,rp),0.0_rp,rp)
         ham_obj%ee(i,i,1,2) = cmplx(0.06_rp*real(i,rp),0.0_rp,rp)
         ham_obj%ee(i,i,2,1) = cmplx(0.003_rp,0.0_rp,rp)
         ham_obj%ee(i,i,2,2) = cmplx(0.003_rp,0.0_rp,rp)
         ham_obj%eeo(i,i,1,1) = cmplx(0.02_rp,0.0_rp,rp)
         ham_obj%eeo(i,i,1,2) = cmplx(0.02_rp,0.0_rp,rp)
         ham_obj%eeo(i,i,2,1) = cmplx(0.001_rp,0.0_rp,rp)
         ham_obj%eeo(i,i,2,2) = cmplx(0.001_rp,0.0_rp,rp)
      end do
      obj%hamiltonian => ham_obj; obj%lattice => lat_obj
      obj%reciprocal_mode = 'generalized_overlap_proxy'; obj%kspace_ham_order = 'first'
      obj%max_orbs = nb; obj%reciprocal_tile_size = 3; obj%cached_operator_generation = 0
      allocate(obj%ham_vec_type(3,2,2), obj%ham_vec_type_direct(3,2,2))
      obj%ham_vec_type = 0.0_rp; obj%ham_vec_type_direct = 0.0_rp
      obj%ham_vec_type_direct(1,2,1) = 0.25_rp
      obj%ham_vec_type_direct(1,2,2) = -0.25_rp
   end subroutine setup_two_site_generalized_model

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

   subroutine check_gauge_normalized(label, left, right, test_failed)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: left(:, :), right(:, :)
      logical, intent(inout) :: test_failed
      integer :: iband
      real(rp) :: max_error, denominator

      max_error = 0.0_rp
      do iband = 1, size(left, 2)
         denominator = sqrt(sum(abs(left(:,iband))**2)*sum(abs(right(:,iband))**2))
         max_error = max(max_error, abs(abs(sum(conjg(left(:,iband))*right(:,iband)))/denominator - 1.0_rp))
      end do
      call check_close(label, max_error, 1.0e-12_rp, test_failed)
   end subroutine check_gauge_normalized

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

   subroutine check_generalized_eigenpairs(label, h, s, eigenvalues, eigenvectors, test_failed)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: h(:, :), s(:, :), eigenvectors(:, :)
      real(rp), intent(in) :: eigenvalues(:)
      logical, intent(inout) :: test_failed
      complex(rp), allocatable :: identity(:, :), residual(:, :), metric_overlap(:, :)
      integer :: n

      n = size(eigenvalues)
      allocate(identity(n,n), residual(n,n), metric_overlap(n,n))
      identity = diag_identity(n)
      residual = matmul(h, eigenvectors) - matmul(s, eigenvectors)*spread(eigenvalues, 1, n)
      metric_overlap = matmul(conjg(transpose(eigenvectors)), matmul(s, eigenvectors)) - identity
      call check_close(trim(label)//' residual', maxval(abs(residual)), 1.0e-12_rp, test_failed)
      call check_close(trim(label)//' metric orthogonality', maxval(abs(metric_overlap)), 1.0e-12_rp, test_failed)
      deallocate(identity, residual, metric_overlap)
   end subroutine check_generalized_eigenpairs

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
