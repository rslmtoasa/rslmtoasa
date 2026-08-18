! ACC-06 benchmark driver.
!
! This is deliberately a benchmark-only executable.  It reuses the production
! reciprocal assembler and typed backend request, but is not a correctness
! test and is never registered as a default CTest test.
program reciprocal_crossover
   use, intrinsic :: iso_fortran_env, only: int64
   use precision_mod, only: rp
   use basis_mod, only: basis_init, nb, norb
   use math_mod, only: init_math_operators
   use reciprocal_mod, only: reciprocal, reciprocal_execution_request, reciprocal_execution_result
   use hamiltonian_mod, only: hamiltonian
   use lattice_mod, only: lattice
   use charge_mod, only: charge
   use control_mod, only: control
   use logger_mod, only: g_logger
   implicit none

   character(len=32) :: backend, strategy, fixture
   integer :: sites, nk, tile_size, lmax, eigenvectors
   type(reciprocal) :: recip
   type(hamiltonian), target :: ham
   type(lattice), target :: lat
   type(charge), target :: chg
   type(control), target :: ctl
   real(rp), allocatable :: points(:, :)
   complex(rp), allocatable :: hamiltonians(:, :, :)
   real(rp) :: assembly_seconds, solve_seconds, total_seconds
   integer :: status

   backend = 'lapack'
   strategy = 'backend'
   sites = 1
   nk = 8
   tile_size = 16
   lmax = 2
   eigenvectors = 1
   call parse_arguments(backend, strategy, sites, nk, tile_size, lmax, eigenvectors)
   if (sites < 1 .or. nk < 1 .or. tile_size < 1 .or. lmax < 1 .or. lmax > 2) then
      error stop 'ACC06: invalid benchmark dimensions'
   end if
   if (trim(backend) == 'cuda' .and. trim(strategy) /= 'backend') then
      error stop 'ACC06: CUDA supports only strategy=backend'
   end if

   call basis_init(lmax)
   call init_math_operators()
   call g_logger%init()
   call setup_model(recip, ham, lat, chg, ctl, sites)
   call make_points(points, nk)
   allocate(hamiltonians(nb*sites, nb*sites, nk))

   call wall_clock(assembly_seconds)
   call assemble_all(recip, points, hamiltonians)
   call wall_clock(total_seconds)
   assembly_seconds = total_seconds - assembly_seconds

   call wall_clock(solve_seconds)
   if (trim(strategy) == 'parallel') then
      call solve_parallel(hamiltonians, nk, tile_size, eigenvectors == 1, status)
   else
      call solve_backend(recip, hamiltonians, nk, tile_size, eigenvectors == 1, trim(backend), status)
   end if
   call wall_clock(total_seconds)
   solve_seconds = total_seconds - solve_seconds
   if (status /= 0) error stop 'ACC06: eigensolver benchmark failed'

   total_seconds = assembly_seconds + solve_seconds
   fixture = fixture_name(sites, lmax)
   write (*, '(a,1x,a,1x,a,1x,a,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0)') &
      'ACC06_DIMENSIONS', 'fixture='//trim(fixture), 'backend='//trim(backend), &
      'strategy='//trim(strategy), 'sites=', sites, 'matrix_dimension=', nb*sites, &
      'nk=', nk, 'tile_size=', tile_size, 'eigenvectors=', eigenvectors, 'lmax=', lmax
   write (*, '(a,1x,a,1x,a,es16.8,1x,a,es16.8,1x,a,es16.8)') &
      'ACC06_TIMING', 'fixture='//trim(fixture), 'assembly_s=', assembly_seconds, &
      'solve_s=', solve_seconds, 'total_s=', total_seconds
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine parse_arguments(backend, strategy, sites, nk, tile_size, lmax, eigenvectors)
      character(len=*), intent(inout) :: backend, strategy
      integer, intent(inout) :: sites, nk, tile_size, lmax, eigenvectors
      integer :: i, narg
      character(len=128) :: key, value

      narg = command_argument_count()
      i = 1
      do while (i <= narg)
         call get_command_argument(i, key)
         if (i == narg) error stop 'ACC06: option requires a value'
         call get_command_argument(i+1, value)
         select case (trim(key))
         case ('--backend'); read(value, *) backend
         case ('--strategy'); read(value, *) strategy
         case ('--sites'); read(value, *) sites
         case ('--nk'); read(value, *) nk
         case ('--tile-size'); read(value, *) tile_size
         case ('--lmax'); read(value, *) lmax
         case ('--eigenvectors'); read(value, *) eigenvectors
         case default; error stop 'ACC06: unknown option '//trim(key)
         end select
         i = i + 2
      end do
      if (trim(backend) /= 'lapack' .and. trim(backend) /= 'cuda') then
         error stop 'ACC06: backend must be lapack or cuda'
      end if
      if (trim(strategy) /= 'backend' .and. trim(strategy) /= 'parallel') then
         error stop 'ACC06: strategy must be backend or parallel'
      end if
      if (eigenvectors /= 0 .and. eigenvectors /= 1) error stop 'ACC06: eigenvectors must be 0 or 1'
   end subroutine parse_arguments

   subroutine setup_model(recip, ham, lat, chg, ctl, nsite)
      type(reciprocal), intent(out) :: recip
      type(hamiltonian), target, intent(out) :: ham
      type(lattice), target, intent(out) :: lat
      type(charge), target, intent(out) :: chg
      type(control), target, intent(out) :: ctl
      integer, intent(in) :: nsite
      integer :: isite, i, next_site, previous_site

      call recip%restore_to_default()
      call ctl%restore_to_default()
      lat%nrec = nsite; lat%ntype = nsite; lat%nn_max = 3; lat%kk = nsite; lat%nmax = nsite
      allocate(lat%ib(nsite), lat%atlist(nsite), lat%iz(nsite), lat%num(nsite), lat%nn(nsite,3), &
         lat%sbar(norb,norb,3,nsite), lat%symbolic_atoms(nsite))
      lat%ib = [(isite, isite=1,nsite)]
      lat%atlist = [(isite, isite=1,nsite)]
      lat%iz = [(isite, isite=1,nsite)]
      lat%num = [(isite, isite=1,nsite)]
      lat%sbar = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, nsite
         next_site = mod(isite, nsite) + 1
         previous_site = mod(isite - 2 + nsite, nsite) + 1
         if (nsite == 1) then
            lat%nn(isite,:) = [3, 1, 1]
         else
            ! Pair each forward bond with its reverse bond so the assembled
            ! Fourier matrix is Hermitian at every k point.
            lat%nn(isite,:) = [3, next_site, previous_site]
         end if
         call lat%symbolic_atoms(isite)%restore_to_default()
         lat%symbolic_atoms(isite)%potential%lmax = lmax
         lat%symbolic_atoms(isite)%potential%mom = [0.0_rp, 0.0_rp, 1.0_rp]
         lat%symbolic_atoms(isite)%potential%wx0 = cmplx(1.0_rp, 0.0_rp, rp)
         lat%symbolic_atoms(isite)%potential%wx1 = cmplx(0.31_rp, 0.0_rp, rp)
         lat%symbolic_atoms(isite)%potential%cx1 = cmplx(0.23_rp, 0.0_rp, rp)
         do i = 1, norb
            lat%sbar(i,i,1,isite) = cmplx(0.12_rp, 0.0_rp, rp)
            lat%sbar(i,i,2,isite) = cmplx(0.37_rp, 0.0_rp, rp)
            lat%sbar(i,i,3,isite) = cmplx(0.37_rp, 0.0_rp, rp)
         end do
      end do

      chg%lattice => lat
      chg%symbolic_atom => lat%symbolic_atoms
      ham%charge => chg; ham%lattice => lat; ham%control => ctl
      ham%hoh = .false.; ham%ccor_2c = .false.; ham%hubbard_u_general_check = .false.
      ham%hubbard_v_check = .false.; ham%hubbard_u_impurity_check = .false.; ham%local_axis = .false.
      ham%magnetic_representation = 'periodic_nc'; ham%operator_generation = 1
      allocate(ham%ee(nb,nb,3,nsite)); ham%ee = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, nsite
         do i = 1, nb
            ham%ee(i,i,1,isite) = cmplx(0.05_rp*real(i,rp) + 0.001_rp*real(isite-1,rp), 0.0_rp, rp)
            ham%ee(i,i,2,isite) = cmplx(0.01_rp*real(i,rp), 0.0_rp, rp)
            ham%ee(i,i,3,isite) = cmplx(0.01_rp*real(i,rp), 0.0_rp, rp)
         end do
      end do

      recip%hamiltonian => ham; recip%lattice => lat; recip%control => ctl
      recip%reciprocal_mode = 'ham_only'; recip%kspace_ham_order = 'first'; recip%max_orbs = nb
      recip%dos_method = 'tetrahedron'; recip%cached_operator_generation = 0
      allocate(recip%ham_vec_type(3,3,nsite), recip%ham_vec_type_direct(3,3,nsite))
      recip%ham_vec_type = 0.0_rp; recip%ham_vec_type_direct = 0.0_rp
      do isite = 1, nsite
         recip%ham_vec_type(1,2,isite) = 0.25_rp
         recip%ham_vec_type_direct(1,2,isite) = recip%ham_vec_type(1,2,isite)
         recip%ham_vec_type(1,3,isite) = -recip%ham_vec_type(1,2,isite)
         recip%ham_vec_type_direct(1,3,isite) = recip%ham_vec_type(1,3,isite)
      end do
   end subroutine setup_model

   subroutine make_points(points, nk)
      real(rp), allocatable, intent(out) :: points(:, :)
      integer, intent(in) :: nk
      integer :: ik
      allocate(points(3,nk))
      do ik = 1, nk
         points(:,ik) = [mod(real(ik-1,rp), 7.0_rp)/7.0_rp - 0.5_rp, &
            mod(real(3*(ik-1),rp), 5.0_rp)/5.0_rp - 0.5_rp, &
            mod(real(5*(ik-1),rp), 3.0_rp)/3.0_rp - 0.5_rp]
      end do
   end subroutine make_points

   subroutine assemble_all(recip, points, hamiltonians)
      type(reciprocal), intent(inout) :: recip
      real(rp), intent(in) :: points(:, :)
      complex(rp), intent(out) :: hamiltonians(:, :, :)
      integer :: ik
      do ik = 1, size(points,2)
         call recip%build_hamiltonian_at_kpoint(points(:,ik), hamiltonians(:,:,ik))
      end do
   end subroutine assemble_all

   subroutine solve_backend(recip, hamiltonians, nk, tile_size, want_vectors, backend, status)
      type(reciprocal), intent(inout) :: recip
      complex(rp), intent(in) :: hamiltonians(:, :, :)
      integer, intent(in) :: nk, tile_size
      logical, intent(in) :: want_vectors
      character(len=*), intent(in) :: backend
      integer, intent(out) :: status
      type(reciprocal_execution_request) :: request
      type(reciprocal_execution_result) :: result
      integer :: first, last, length

      status = 0
      call recip%make_execution_backend(backend)
      do first = 1, nk, tile_size
         last = min(nk, first + tile_size - 1); length = last - first + 1
         request%assemble_hamiltonian = .false.; request%assemble_overlap = .false.
         request%solve_eigensystem = .true.; request%generalized = .false.
         request%request_eigenvectors = want_vectors; request%operator_generation = 1
         if (allocated(request%input_hamiltonian)) deallocate(request%input_hamiltonian)
         allocate(request%input_hamiltonian(size(hamiltonians,1), size(hamiltonians,2), length), &
            source=hamiltonians(:,:,first:last))
         call recip%execution_backend%execute_batch(request, result)
         call recip%execution_backend%synchronize()
         if (.not. result%eigenvalues_valid) status = 1
         if (result%local_point_count /= length) status = 1
         if (want_vectors .and. .not. result%eigenvectors_valid) status = 1
         if ((.not. want_vectors) .and. result%eigenvectors_valid) status = 1
         if (allocated(result%eigenvalues)) deallocate(result%eigenvalues)
         if (allocated(result%eigenvectors)) deallocate(result%eigenvectors)
      end do
   end subroutine solve_backend

   subroutine solve_parallel(hamiltonians, nk, tile_size, want_vectors, status)
      complex(rp), intent(in) :: hamiltonians(:, :, :)
      integer, intent(in) :: nk, tile_size
      logical, intent(in) :: want_vectors
      integer, intent(out) :: status
      integer :: first, last, ik, n, info
      complex(rp), allocatable :: matrix(:, :), work(:)
      real(rp), allocatable :: eval(:), rwork(:)
      external :: zheev

      ! This is intentionally benchmark-only.  The production backend owns a
      ! reusable workspace and is not made OpenMP-shared by this experiment.
      status = 0; n = size(hamiltonians,1)
      do first = 1, nk, tile_size
         last = min(nk, first + tile_size - 1)
!$omp parallel do default(none) shared(hamiltonians,first,last,n,want_vectors,status) &
!$omp& private(ik,matrix,work,eval,rwork,info)
         do ik = first, last
            allocate(matrix(n,n), eval(n), work(max(1,2*n)), rwork(max(1,3*n-2)))
            matrix = hamiltonians(:,:,ik)
            call zheev(merge('V','N',want_vectors), 'U', n, matrix, n, eval, work, size(work), rwork, info)
            if (info /= 0) then
!$omp critical(acc06_status)
               status = 1
!$omp end critical(acc06_status)
            end if
            deallocate(matrix, eval, work, rwork)
         end do
!$omp end parallel do
      end do
   end subroutine solve_parallel

   function fixture_name(nsite, lmax) result(name)
      integer, intent(in) :: nsite, lmax
      character(len=32) :: name
      if (lmax == 1 .and. nsite == 1) then
         name = 'Si_sp'
      else if (lmax == 2 .and. nsite == 1) then
         name = 'bccFe_spd'
      else if (lmax == 2 .and. nsite == 2) then
         name = 'two_site_spd'
      else
         write(name, '(a,i0,a,i0)') 'multisite_', nsite, '_lmax_', lmax
      end if
   end function fixture_name

   subroutine wall_clock(seconds)
      real(rp), intent(out) :: seconds
      integer(int64) :: count, rate
      call system_clock(count, rate)
      seconds = real(count, rp) / real(rate, rp)
   end subroutine wall_clock

end program reciprocal_crossover
