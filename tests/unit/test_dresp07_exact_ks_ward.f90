! DRESP-07 exact finite-H Ward-closure regression.
!
! This test is independent of the production driver.  It locks the accepted
! site-major spin ordering, the commutator/divided-difference identity, the
! circular factor-of-two convention, the exact B_H extraction, and the full
! spd product contraction.
program test_dresp07_exact_ks_ward
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use logger_mod, only: g_logger
   use math_mod, only: init_math_operators
   use self_mod, only: legacy_radial_fixture
   use lmto_radial_augmentation_mod, only: lmto_radial_basis
   use lr_response_space_mod, only: response_space_layout
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis, lmto_product_channel_plus
   use lr_exact_ks_ward_mod, only: exact_ks_ward_metrics, exact_ks_global_rotation_oracle, &
      exact_ks_product_rotation_sweep, exact_ks_extract_collinear_field, exact_ks_build_transverse_field, &
      exact_ks_fermi_occupation, exact_ks_fermi_derivative
   implicit none

   integer, parameter :: nr = 51, lmax = 2, norb = (lmax + 1)**2, n = 2*norb
   real(rp), parameter :: mesh_a = 0.03_rp, mesh_b = 0.10_rp, nuclear_z = 1.0_rp
   real(rp), parameter :: fermi = 0.0_rp, temperature = 300.0_rp, eta = 0.01_rp, tolerance = 2.0e-10_rp
   real(rp) :: radius(nr), evals(n), weights(1), max_eigenpair, max_density, max_density_relative
   complex(rp) :: hamiltonian(n, n), eigenvectors(n, n), generator(n, n), delta_h(n, n), density(n, n)
   complex(rp) :: delta_rho_rot(n, n), delta_rho_spec(n, n), b_orbital(norb, norb), b_spin(n, n)
   complex(rp) :: delta_h_field(n, n), sigma_x(n, n), sigma_plus(n, n), operator_eigen(n, n), h_eigen(n, n)
   complex(rp) :: circular_response, mx_response
   complex(rp), allocatable :: response_rot(:), response_spec(:), response_retarded(:)
   type(lmto_radial_basis) :: radial
   type(response_space_layout) :: response_space
   type(lmto_product_response_basis) :: product
   type(exact_ks_ward_metrics) :: metrics
   real(rp) :: onsite_norm, nonlocal_norm, spin_offdiag_norm, kernel, product_error_relative
   integer :: i, j, info
   logical :: failed
   external :: zheev

   call basis_init(lmax)
   call init_math_operators()
   call g_logger%init()
   call build_mesh(radius)
   call setup_radial_basis(radial, radius)
   call response_space%initialize(1, 2*lmax, radius, mesh_a, mesh_b, 1)
   call product%initialize(response_space, [radial], lmto_product_channel_plus, .false.)
   if (product%product_dimension < 1) error stop 'DRESP-07 unit: empty spd product representation'
   call build_fixture_hamiltonian(hamiltonian)
   eigenvectors = hamiltonian
   call diagonalize(eigenvectors, evals, info)
   if (info /= 0) error stop 'DRESP-07 unit: eigensystem failed'

   call exact_ks_global_rotation_oracle(hamiltonian, evals, eigenvectors, fermi, temperature, generator, delta_h, density, &
      delta_rho_rot, delta_rho_spec, metrics, norb, 1)
   failed = metrics%max_eigenpair_residual > tolerance .or. metrics%hermiticity_residual > tolerance .or. &
      metrics%rotation_field_relation_residual > tolerance .or. metrics%density_rotation_spectral_relative > tolerance
   call exact_ks_extract_collinear_field(hamiltonian, norb, 1, b_orbital, b_spin, onsite_norm, nonlocal_norm, &
      spin_offdiag_norm)
   if (maxval(abs(b_orbital - 0.5_rp*(hamiltonian(1:norb, 1:norb) - hamiltonian(norb+1:n, norb+1:n)))) > tolerance) failed = .true.
   call exact_ks_build_transverse_field(b_orbital, norb, 1, delta_h_field)
   if (maxval(abs(delta_h_field - delta_h)) > tolerance .or. spin_offdiag_norm > tolerance) failed = .true.

   sigma_x = cmplx(0.0_rp, 0.0_rp, rp)
   sigma_plus = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, norb
      sigma_x(i, norb + i) = cmplx(1.0_rp, 0.0_rp, rp)
      sigma_x(norb + i, i) = cmplx(1.0_rp, 0.0_rp, rp)
      sigma_plus(i, norb + i) = cmplx(1.0_rp, 0.0_rp, rp)
   end do
   mx_response = sum(sigma_x*transpose(delta_rho_spec))
   operator_eigen = matmul(conjg(transpose(eigenvectors)), matmul(sigma_plus, eigenvectors))
   h_eigen = matmul(conjg(transpose(eigenvectors)), matmul(delta_h, eigenvectors))
   circular_response = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, n
      do j = 1, n
         if (abs(evals(i) - evals(j)) <= 1.0e-11_rp*max(1.0_rp, abs(evals(i)), abs(evals(j)))) then
            kernel = exact_ks_fermi_derivative(evals(i), fermi, temperature)
         else
            kernel = (exact_ks_fermi_occupation(evals(i), fermi, temperature) - &
               exact_ks_fermi_occupation(evals(j), fermi, temperature))/(evals(i) - evals(j))
         end if
         circular_response = circular_response + operator_eigen(i, j)*kernel*h_eigen(j, i)
      end do
   end do
   if (abs(real(mx_response, rp) - 2.0_rp*real(circular_response, rp)) > tolerance .or. &
       abs(aimag(mx_response) - 2.0_rp*aimag(circular_response)) > tolerance) failed = .true.

   allocate(response_rot(product%product_dimension), response_spec(product%product_dimension), &
      response_retarded(product%product_dimension))
   weights = 1.0_rp
   call exact_ks_product_rotation_sweep(product, reshape(hamiltonian, [n, n, 1]), reshape(evals, [n, 1]), &
      reshape(eigenvectors, [n, n, 1]), weights, fermi, temperature, eta, response_rot, response_spec, response_retarded, &
      max_eigenpair, max_density, max_density_relative)
   product_error_relative = sqrt(sum(abs(response_rot - response_spec)**2))/max(sqrt(sum(abs(response_spec)**2)), tiny(1.0_rp))
   if (max_eigenpair > tolerance .or. max_density_relative > tolerance .or. product_error_relative > tolerance) failed = .true.

   if (failed) then
      write (*, '(a)') 'UnitDresp07ExactKsWard: FAIL'
      error stop 1
   end if
   write (*, '(a,es12.4)') '  exact_density_relative=', metrics%density_rotation_spectral_relative
   write (*, '(a,es12.4)') '  circular_factor_two_error=', abs(mx_response - 2.0_rp*circular_response)
   write (*, '(a,es12.4)') '  product_rotation_spectral_relative=', product_error_relative
   write (*, '(a)') 'UnitDresp07ExactKsWard: PASS (H contract, divided difference, circular factor, B_H, product closure)'

contains

   subroutine build_mesh(mesh)
      real(rp), intent(out) :: mesh(:)
      integer :: ir
      mesh(1) = 0.0_rp
      do ir = 2, size(mesh)
         mesh(ir) = mesh_b*(exp(mesh_a*real(ir - 1, rp)) - 1.0_rp)
      end do
   end subroutine build_mesh

   subroutine setup_radial_basis(basis, mesh)
      type(lmto_radial_basis), intent(out) :: basis
      real(rp), intent(in) :: mesh(:)
      real(rp) :: potential(size(mesh)), energy
      real(rp), allocatable :: g(:, :), gp(:, :), gpp(:, :), gpack(:), gdotpack(:), gddotpack(:)
      integer :: ispin, l
      potential = 0.0_rp
      call basis%initialize(size(mesh), lmax, 2)
      do ispin = 1, 2
         do l = 0, lmax
            call legacy_radial_fixture(nuclear_z, l, mesh_a, mesh_b, mesh, potential, energy, g, gp, gpp, 2.0_rp)
            gpack = reshape(g, [2*size(mesh)])
            gdotpack = reshape(gp, [2*size(mesh)])
            gddotpack = reshape(gpp, [2*size(mesh)])
            call basis%capture_channel(l, ispin, energy, mesh, potential, mesh_a, mesh_b, nuclear_z, gpack, gdotpack, &
               gddotpack, energy)
         end do
      end do
   end subroutine setup_radial_basis

   subroutine build_fixture_hamiltonian(hamiltonian)
      complex(rp), intent(out) :: hamiltonian(:, :)
      complex(rp) :: up(norb, norb), down(norb, norb)
      integer :: i, j
      up = cmplx(0.0_rp, 0.0_rp, rp)
      down = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, norb
         up(i, i) = cmplx(-0.48_rp + 0.13_rp*real(i, rp), 0.0_rp, rp)
         down(i, i) = cmplx(-0.31_rp + 0.12_rp*real(i, rp), 0.0_rp, rp)
         do j = i + 1, norb
            up(i, j) = cmplx(0.006_rp*real(i + j, rp), 0.002_rp*real(j - i, rp), rp)
            down(i, j) = cmplx(-0.004_rp*real(i + j, rp), 0.001_rp*real(j - i, rp), rp)
            up(j, i) = conjg(up(i, j))
            down(j, i) = conjg(down(i, j))
         end do
      end do
      hamiltonian = cmplx(0.0_rp, 0.0_rp, rp)
      hamiltonian(1:norb, 1:norb) = up
      hamiltonian(norb+1:n, norb+1:n) = down
   end subroutine build_fixture_hamiltonian

   subroutine diagonalize(matrix, eigenvalues, info)
      complex(rp), intent(inout) :: matrix(:, :)
      real(rp), intent(out) :: eigenvalues(:)
      integer, intent(out) :: info
      complex(rp) :: work_query(1)
      complex(rp), allocatable :: work(:)
      real(rp), allocatable :: rwork(:)
      integer :: lwork
      allocate(rwork(max(1, 3*size(eigenvalues) - 2)))
      call zheev('V', 'U', size(eigenvalues), matrix, size(eigenvalues), eigenvalues, work_query, -1, rwork, info)
      lwork = max(1, nint(real(work_query(1), rp)))
      allocate(work(lwork))
      call zheev('V', 'U', size(eigenvalues), matrix, size(eigenvalues), eigenvalues, work, lwork, rwork, info)
      deallocate(work, rwork)
   end subroutine diagonalize

end program test_dresp07_exact_ks_ward
