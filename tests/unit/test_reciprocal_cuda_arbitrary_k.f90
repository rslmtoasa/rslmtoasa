!------------------------------------------------------------------------------
! ACC-04 -- production arbitrary-k CPU/CUDA reciprocal integration
!------------------------------------------------------------------------------
program test_reciprocal_cuda_arbitrary_k
   use, intrinsic :: iso_c_binding, only: c_int
   use precision_mod, only: rp
   use basis_mod, only: nb
   use reciprocal_mod, only: reciprocal
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use logger_mod, only: g_logger
   implicit none

   interface
      function cuda_device_count(count) bind(C, name='rslmto_reciprocal_cuda_device_count') result(status)
         import c_int
         integer(c_int), intent(out) :: count
         integer(c_int) :: status
      end function cuda_device_count
   end interface

   type(reciprocal) :: recip
   type(hamiltonian), target :: ham
   type(lattice), target :: lat
   real(rp), parameter :: points(3,4) = reshape([ &
      0.125_rp, 0.0_rp, 0.0_rp, &
      0.371_rp, -0.217_rp, 0.083_rp, &
      1.125_rp, 0.0_rp, 0.0_rp, &
      -0.125_rp, 0.0_rp, 0.0_rp], [3,4])
   real(rp), allocatable :: cpu_values(:, :), cuda_values(:, :)
   complex(rp), allocatable :: cpu_vectors(:, :, :), cuda_vectors(:, :, :)
   complex(rp), allocatable :: h(:, :), cpu_projector(:, :), cuda_projector(:, :)
   integer(c_int) :: device_count, status
   integer :: ik
   real(rp) :: eigenvalue_error, projector_error, residual_error, orthogonality_error

   call g_logger%init()
   status = cuda_device_count(device_count)
   if (status /= 0_c_int .or. device_count <= 0_c_int) then
      write (*, '(a)') 'SKIP: no CUDA device is available'
      stop 77
   end if

   call setup_one_site_model(recip, ham, lat)
   recip%reciprocal_tile_size = 2

   call recip%make_execution_backend('lapack')
   call recip%calculate_eigenpairs_at_kpoints(points, cpu_values, cpu_vectors)

   call recip%make_execution_backend('cuda')
   call recip%calculate_eigenpairs_at_kpoints(points, cuda_values, cuda_vectors)

   eigenvalue_error = maxval(abs(cpu_values - cuda_values))
   if (eigenvalue_error > 5.0e-11_rp) then
      write (*, '(a,es14.6)') 'FAIL: CPU/CUDA arbitrary-k eigenvalue error = ', eigenvalue_error
      error stop 1
   end if

   allocate(h(nb,nb), cpu_projector(nb,nb), cuda_projector(nb,nb))
   residual_error = 0.0_rp
   orthogonality_error = 0.0_rp
   projector_error = 0.0_rp
   do ik = 1, size(points,2)
      call recip%build_hamiltonian_at_kpoint(points(:,ik), h)
      residual_error = max(residual_error, maxval(abs(matmul(h, cuda_vectors(:,:,ik)) - &
         cuda_vectors(:,:,ik)*spread(cuda_values(:,ik), 1, nb))))
      orthogonality_error = max(orthogonality_error, maxval(abs(&
         matmul(conjg(transpose(cuda_vectors(:,:,ik))), cuda_vectors(:,:,ik)) - identity(nb))))
      cpu_projector = matmul(cpu_vectors(:,:,ik), conjg(transpose(cpu_vectors(:,:,ik))))
      cuda_projector = matmul(cuda_vectors(:,:,ik), conjg(transpose(cuda_vectors(:,:,ik))))
      projector_error = max(projector_error, maxval(abs(cpu_projector - cuda_projector)))
   end do

   if (residual_error > 5.0e-11_rp .or. orthogonality_error > 5.0e-11_rp .or. projector_error > 5.0e-10_rp) then
      write (*, '(a,3(es14.6,1x))') 'FAIL: CUDA residual/orthogonality/projector errors = ', &
         residual_error, orthogonality_error, projector_error
      error stop 1
   end if
   if (maxval(abs(cuda_values(:,1) - cuda_values(:,3))) > 5.0e-11_rp .or. &
       maxval(abs(cuda_values(:,1) - cuda_values(:,4))) > 5.0e-11_rp) then
      write (*, '(a)') 'FAIL: CUDA folded/duplicate arbitrary-k behavior changed'
      error stop 1
   end if

   write (*, '(a,3(es14.6,1x))') 'PASS: CUDA arbitrary-k eigenpairs; eigenvalue/projector/residual errors = ', &
      eigenvalue_error, projector_error, residual_error

contains

   subroutine setup_one_site_model(obj, ham_obj, lat_obj)
      type(reciprocal), intent(out) :: obj
      type(hamiltonian), target, intent(out) :: ham_obj
      type(lattice), target, intent(out) :: lat_obj
      integer :: i

      call obj%restore_to_default()
      lat_obj%nrec = 1; lat_obj%ntype = 1; lat_obj%nn_max = 3; lat_obj%kk = 1
      allocate(lat_obj%ib(1), lat_obj%atlist(1), lat_obj%iz(1), lat_obj%nn(1,3))
      lat_obj%ib = 1; lat_obj%atlist = 1; lat_obj%iz = 1; lat_obj%nn(1,:) = [3,1,1]
      ham_obj%lattice => lat_obj; ham_obj%hoh = .false.; ham_obj%ccor_2c = .false.
      allocate(ham_obj%ee(nb,nb,3,1)); ham_obj%ee = (0.0_rp,0.0_rp)
      do i = 1, nb
         ham_obj%ee(i,i,1,1) = cmplx(0.05_rp*real(i,rp),0.0_rp,rp)
         ham_obj%ee(i,i,2,1) = cmplx(0.01_rp*real(i,rp),0.0_rp,rp)
         ham_obj%ee(i,i,3,1) = cmplx(0.01_rp*real(i,rp),0.0_rp,rp)
      end do
      obj%hamiltonian => ham_obj; obj%lattice => lat_obj
      obj%reciprocal_mode = 'ham_only'; obj%kspace_ham_order = 'first'; obj%max_orbs = nb
      obj%dos_method = 'tetrahedron'; obj%cached_operator_generation = 0
      ham_obj%operator_generation = 0; ham_obj%magnetic_representation = 'periodic_nc'
      allocate(obj%ham_vec_type(3,3,1), obj%ham_vec_type_direct(3,3,1))
      obj%ham_vec_type = 0.0_rp; obj%ham_vec_type_direct = 0.0_rp
      obj%ham_vec_type_direct(1,2,1) = 1.0_rp; obj%ham_vec_type_direct(1,3,1) = -1.0_rp
   end subroutine setup_one_site_model

   function identity(n) result(value)
      integer, intent(in) :: n
      complex(rp) :: value(n,n)
      integer :: i
      value = cmplx(0.0_rp,0.0_rp,rp)
      do i = 1, n
         value(i,i) = cmplx(1.0_rp,0.0_rp,rp)
      end do
   end function identity

end program test_reciprocal_cuda_arbitrary_k
