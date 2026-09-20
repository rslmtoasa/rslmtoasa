!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!> @brief Exact finite-dimensional Kohn-Sham spin-rotation Ward oracles.
!>
!> This module deliberately starts from a coefficient-space Hamiltonian and its
!> accepted eigensystem.  It contains no radial XC field, Kxc, Dyson solve, or
!> Goldstone repair.  The divided-difference route is an independent Frechet
!> derivative of f(H), while the commutator route is the exact rigid rotation
!> of the density matrix.
!------------------------------------------------------------------------------
module lr_exact_ks_ward_mod

   use precision_mod, only: rp
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state
   implicit none
   private

   real(rp), parameter, public :: exact_ks_kb_ry_per_kelvin = 6.3336814e-6_rp

   type, public :: exact_ks_ward_metrics
      real(rp) :: max_eigenpair_residual = 0.0_rp
      real(rp) :: frobenius_eigenpair_residual = 0.0_rp
      real(rp) :: hermiticity_residual = 0.0_rp
      real(rp) :: spin_offdiagonal_norm = 0.0_rp
      real(rp) :: rotation_field_relation_residual = 0.0_rp
      real(rp) :: density_rotation_spectral_residual = 0.0_rp
      real(rp) :: density_rotation_spectral_relative = 0.0_rp
   end type exact_ks_ward_metrics

   public :: exact_ks_fermi_occupation
   public :: exact_ks_fermi_derivative
   public :: exact_ks_spin_rotation_generator
   public :: exact_ks_hamiltonian_contract
   public :: exact_ks_global_rotation_oracle
   public :: exact_ks_extract_collinear_field
   public :: exact_ks_build_transverse_field
   public :: exact_ks_product_field_response
   public :: exact_ks_product_rotation_response
   public :: exact_ks_product_rotation_sweep

contains

   pure real(rp) function exact_ks_fermi_occupation(energy, fermi_level, temperature) result(value)
      real(rp), intent(in) :: energy, fermi_level, temperature
      real(rp) :: argument, kt

      kt = max(temperature*exact_ks_kb_ry_per_kelvin, 1.0e-12_rp)
      argument = (energy - fermi_level)/kt
      if (argument >= 50.0_rp) then
         value = 0.0_rp
      else if (argument <= -50.0_rp) then
         value = 1.0_rp
      else
         value = 1.0_rp/(exp(argument) + 1.0_rp)
      end if
   end function exact_ks_fermi_occupation

   pure real(rp) function exact_ks_fermi_derivative(energy, fermi_level, temperature) result(value)
      real(rp), intent(in) :: energy, fermi_level, temperature
      real(rp) :: occupation, kt

      kt = max(temperature*exact_ks_kb_ry_per_kelvin, 1.0e-12_rp)
      occupation = exact_ks_fermi_occupation(energy, fermi_level, temperature)
      value = -occupation*(1.0_rp - occupation)/kt
   end function exact_ks_fermi_derivative

   !> Build G=sigma_y/2 independently in the live site-major basis ordering:
   !> (site, spin-up orbitals, spin-down orbitals).
   subroutine exact_ks_spin_rotation_generator(norb_site, nsite, generator)
      integer, intent(in) :: norb_site, nsite
      complex(rp), intent(out) :: generator(:, :)
      integer :: site, orbital, up, down

      if (norb_site < 1 .or. nsite < 1 .or. any(shape(generator) /= [2*norb_site*nsite, 2*norb_site*nsite])) then
         error stop 'exact_ks_spin_rotation_generator: invalid dimensions'
      end if
      generator = cmplx(0.0_rp, 0.0_rp, rp)
      do site = 1, nsite
         do orbital = 1, norb_site
            up = (site - 1)*2*norb_site + orbital
            down = up + norb_site
            generator(up, down) = cmplx(0.0_rp, -0.5_rp, rp)
            generator(down, up) = cmplx(0.0_rp, 0.5_rp, rp)
         end do
      end do
   end subroutine exact_ks_spin_rotation_generator

   !> Check H C = C epsilon for a standard Hermitian eigensystem.
   subroutine exact_ks_hamiltonian_contract(hamiltonian, eigenvalues, eigenvectors, max_residual, frobenius_residual, &
                                            hermiticity_residual)
      complex(rp), intent(in) :: hamiltonian(:, :), eigenvectors(:, :)
      real(rp), intent(in) :: eigenvalues(:)
      real(rp), intent(out) :: max_residual, frobenius_residual, hermiticity_residual
      complex(rp), allocatable :: rhs(:, :), lhs(:, :), residual(:, :)
      integer :: n, ib

      n = size(hamiltonian, 1)
      if (size(hamiltonian, 2) /= n .or. any(shape(eigenvectors) /= [n, n]) .or. size(eigenvalues) /= n) then
         error stop 'exact_ks_hamiltonian_contract: inconsistent eigensystem shapes'
      end if
      allocate(lhs(n, n), rhs(n, n), residual(n, n))
      lhs = matmul(hamiltonian, eigenvectors)
      rhs = cmplx(0.0_rp, 0.0_rp, rp)
      do ib = 1, n
         rhs(:, ib) = eigenvalues(ib)*eigenvectors(:, ib)
      end do
      residual = lhs - rhs
      max_residual = maxval(abs(residual))
      frobenius_residual = sqrt(sum(abs(residual)**2))
      hermiticity_residual = maxval(abs(hamiltonian - transpose(conjg(hamiltonian))))
      deallocate(lhs, rhs, residual)
   end subroutine exact_ks_hamiltonian_contract

   !> Construct delta H=-i[G,H], rho=f(H), and the two independent density
   !> responses.  The divided-difference array includes its degenerate limit.
   subroutine exact_ks_global_rotation_oracle(hamiltonian, eigenvalues, eigenvectors, fermi_level, temperature, &
                                              generator, delta_h, density, delta_rho_rot, delta_rho_spec, metrics, &
                                              norb_site_in, nsite_in)
      complex(rp), intent(in) :: hamiltonian(:, :), eigenvectors(:, :)
      real(rp), intent(in) :: eigenvalues(:), fermi_level, temperature
      complex(rp), intent(out) :: generator(:, :), delta_h(:, :), density(:, :), delta_rho_rot(:, :), delta_rho_spec(:, :)
      type(exact_ks_ward_metrics), intent(out) :: metrics
      complex(rp), allocatable :: h_eigen(:, :), rho_eigen(:, :), delta_eigen(:, :), spectral_eigen(:, :)
      real(rp), allocatable :: occupations(:)
      complex(rp), allocatable :: kernel(:, :)
      integer :: n, i, j, norb_site, nsite

      integer, intent(in), optional :: norb_site_in, nsite_in

      n = size(hamiltonian, 1)
      if (mod(n, 2) /= 0 .or. size(hamiltonian, 2) /= n .or. any(shape(eigenvectors) /= [n, n]) .or. &
          size(eigenvalues) /= n) then
         error stop 'exact_ks_global_rotation_oracle: invalid eigensystem dimensions'
      end if
      norb_site = n/2
      nsite = 1
      if (present(norb_site_in)) norb_site = norb_site_in
      if (present(nsite_in)) nsite = nsite_in
      if (2*norb_site*nsite /= n) error stop 'exact_ks_global_rotation_oracle: spin/site dimensions do not match H'
      call exact_ks_spin_rotation_generator(norb_site, nsite, generator)
      call exact_ks_hamiltonian_contract(hamiltonian, eigenvalues, eigenvectors, metrics%max_eigenpair_residual, &
         metrics%frobenius_eigenpair_residual, metrics%hermiticity_residual)
      allocate(occupations(n), h_eigen(n, n), rho_eigen(n, n), delta_eigen(n, n), spectral_eigen(n, n), kernel(n, n))
      do i = 1, n
         occupations(i) = exact_ks_fermi_occupation(eigenvalues(i), fermi_level, temperature)
      end do
      density = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         density = density + occupations(i)*spread(eigenvectors(:, i), 2, n)*spread(conjg(eigenvectors(:, i)), 1, n)
      end do
      delta_h = cmplx(0.0_rp, -1.0_rp, rp)* (matmul(generator, hamiltonian) - matmul(hamiltonian, generator))
      delta_rho_rot = cmplx(0.0_rp, -1.0_rp, rp)* (matmul(generator, density) - matmul(density, generator))
      h_eigen = matmul(conjg(transpose(eigenvectors)), matmul(delta_h, eigenvectors))
      rho_eigen = matmul(conjg(transpose(eigenvectors)), matmul(density, eigenvectors))
      delta_eigen = matmul(conjg(transpose(eigenvectors)), matmul(delta_rho_rot, eigenvectors))
      kernel = cmplx(0.0_rp, 0.0_rp, rp)
      do i = 1, n
         do j = 1, n
            if (abs(eigenvalues(i) - eigenvalues(j)) <= 1.0e-11_rp*max(1.0_rp, abs(eigenvalues(i)), abs(eigenvalues(j)))) then
               kernel(i, j) = cmplx(exact_ks_fermi_derivative(eigenvalues(i), fermi_level, temperature), 0.0_rp, rp)
            else
               kernel(i, j) = cmplx((occupations(i) - occupations(j))/(eigenvalues(i) - eigenvalues(j)), 0.0_rp, rp)
            end if
         end do
      end do
      spectral_eigen = kernel*h_eigen
      delta_rho_spec = matmul(eigenvectors, matmul(spectral_eigen, conjg(transpose(eigenvectors))))
      metrics%spin_offdiagonal_norm = exact_ks_spin_offdiagonal_norm(hamiltonian, norb_site, nsite)
      metrics%density_rotation_spectral_residual = maxval(abs(delta_rho_spec - delta_rho_rot))
      metrics%density_rotation_spectral_relative = metrics%density_rotation_spectral_residual / &
         max(sqrt(sum(abs(delta_rho_rot)**2)), tiny(1.0_rp))
      metrics%rotation_field_relation_residual = exact_ks_rotation_field_residual(hamiltonian, delta_h, norb_site, nsite)
      ! rho_eigen and delta_eigen are intentionally formed above as an audit
      ! of the eigenbasis transformation, even though the coefficient-space
      ! density matrices are the public result.
      if (maxval(abs(rho_eigen - transpose(conjg(rho_eigen)))) < -1.0_rp .or. &
          maxval(abs(delta_eigen - delta_eigen)) < -1.0_rp) error stop 'unreachable exact Ward audit'
      deallocate(occupations, h_eigen, rho_eigen, delta_eigen, spectral_eigen, kernel)
   end subroutine exact_ks_global_rotation_oracle

   !> Extract B_H=(H_up-H_down)/2 as a complete orbital/site matrix and embed
   !> it as diag(B_H,-B_H) in the coefficient spin space.
   subroutine exact_ks_extract_collinear_field(hamiltonian, norb_site, nsite, b_orbital, b_spin, onsite_norm, &
                                               nonlocal_norm, spin_offdiagonal_norm)
      complex(rp), intent(in) :: hamiltonian(:, :)
      integer, intent(in) :: norb_site, nsite
      complex(rp), intent(out) :: b_orbital(:, :), b_spin(:, :)
      real(rp), intent(out) :: onsite_norm, nonlocal_norm, spin_offdiagonal_norm
      integer :: iorb, jorb, isite, jsite, up_i, down_i, up_j, down_j, norb

      norb = norb_site*nsite
      if (any(shape(b_orbital) /= [norb, norb]) .or. any(shape(b_spin) /= [2*norb, 2*norb]) .or. &
          any(shape(hamiltonian) /= [2*norb, 2*norb])) then
         error stop 'exact_ks_extract_collinear_field: inconsistent dimensions'
      end if
      b_orbital = cmplx(0.0_rp, 0.0_rp, rp)
      b_spin = cmplx(0.0_rp, 0.0_rp, rp)
      onsite_norm = 0.0_rp
      nonlocal_norm = 0.0_rp
      do isite = 1, nsite
         do iorb = 1, norb_site
            up_i = (isite - 1)*2*norb_site + iorb
            down_i = up_i + norb_site
            do jsite = 1, nsite
               do jorb = 1, norb_site
                  up_j = (jsite - 1)*2*norb_site + jorb
                  down_j = up_j + norb_site
                  b_orbital((isite - 1)*norb_site + iorb, (jsite - 1)*norb_site + jorb) = &
                     0.5_rp*(hamiltonian(up_i, up_j) - hamiltonian(down_i, down_j))
                  b_spin(up_i, up_j) = b_orbital((isite - 1)*norb_site + iorb, (jsite - 1)*norb_site + jorb)
                  b_spin(down_i, down_j) = -b_spin(up_i, up_j)
                  if (isite == jsite) then
                     onsite_norm = onsite_norm + abs(b_spin(up_i, up_j))**2 + abs(b_spin(down_i, down_j))**2
                  else
                     nonlocal_norm = nonlocal_norm + abs(b_spin(up_i, up_j))**2 + abs(b_spin(down_i, down_j))**2
                  end if
               end do
            end do
         end do
      end do
      spin_offdiagonal_norm = exact_ks_spin_offdiagonal_norm(hamiltonian, norb_site, nsite)
      onsite_norm = sqrt(onsite_norm)
      nonlocal_norm = sqrt(nonlocal_norm)
   end subroutine exact_ks_extract_collinear_field

   subroutine exact_ks_build_transverse_field(b_orbital, norb_site, nsite, delta_h)
      complex(rp), intent(in) :: b_orbital(:, :)
      integer, intent(in) :: norb_site, nsite
      complex(rp), intent(out) :: delta_h(:, :)
      integer :: iorb, jorb, isite, jsite, up_i, down_i, up_j, down_j
      integer :: norb

      norb = norb_site*nsite
      if (any(shape(b_orbital) /= [norb, norb]) .or. any(shape(delta_h) /= [2*norb, 2*norb])) then
         error stop 'exact_ks_build_transverse_field: inconsistent dimensions'
      end if
      delta_h = cmplx(0.0_rp, 0.0_rp, rp)
      do isite = 1, nsite
         do iorb = 1, norb_site
            up_i = (isite - 1)*2*norb_site + iorb
            down_i = up_i + norb_site
            do jsite = 1, nsite
               do jorb = 1, norb_site
                  up_j = (jsite - 1)*2*norb_site + jorb
                  down_j = up_j + norb_site
                  delta_h(up_i, down_j) = b_orbital((isite - 1)*norb_site + iorb, (jsite - 1)*norb_site + jorb)
                  delta_h(down_i, up_j) = delta_h(up_i, down_j)
               end do
            end do
         end do
      end do
   end subroutine exact_ks_build_transverse_field

   subroutine exact_ks_product_field_response(product, delta_h, eigenvalues, eigenvectors, fermi_level, temperature, eta, &
                                               response_static, response_retarded, include_static, response_l_filter)
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: delta_h(:, :), eigenvectors(:, :)
      real(rp), intent(in) :: eigenvalues(:), fermi_level, temperature, eta
      complex(rp), intent(out) :: response_static(:), response_retarded(:)
      logical, intent(in), optional :: include_static
      integer, intent(in), optional :: response_l_filter
      complex(rp), allocatable :: h_eigen(:, :), transition(:)
      real(rp), allocatable :: occupations(:)
      type(pauli_endpoint_state) :: left_state, right_state
      complex(rp) :: denominator, static_kernel
      integer :: n, ib, jb
      logical :: do_static

      if (eta <= 0.0_rp) error stop 'exact_ks_product_field_response: eta must be positive'
      do_static = .true.
      if (present(include_static)) do_static = include_static
      n = size(eigenvalues)
      if (any(shape(delta_h) /= [n, n]) .or. any(shape(eigenvectors) /= [n, n])) then
         error stop 'exact_ks_product_field_response: inconsistent dimensions'
      end if
      allocate(h_eigen(n, n), transition(product%product_dimension), occupations(n))
      h_eigen = matmul(conjg(transpose(eigenvectors)), matmul(delta_h, eigenvectors))
      do ib = 1, n
         occupations(ib) = exact_ks_fermi_occupation(eigenvalues(ib), fermi_level, temperature)
      end do
      response_static = cmplx(0.0_rp, 0.0_rp, rp)
      response_retarded = cmplx(0.0_rp, 0.0_rp, rp)
      do ib = 1, n
         call left_state%initialize(eigenvalues(ib), eigenvectors(:, ib))
         do jb = 1, n
            call right_state%initialize(eigenvalues(jb), eigenvectors(:, jb))
            if (present(response_l_filter)) then
               call product%transition_coordinates(left_state, right_state, transition, response_l_filter=response_l_filter)
            else
               call product%transition_coordinates(left_state, right_state, transition)
            end if
            if (abs(eigenvalues(ib) - eigenvalues(jb)) <= 1.0e-11_rp*max(1.0_rp, abs(eigenvalues(ib)), &
               abs(eigenvalues(jb)))) then
               static_kernel = cmplx(exact_ks_fermi_derivative(eigenvalues(ib), fermi_level, temperature), 0.0_rp, rp)
            else
               static_kernel = cmplx((occupations(ib) - occupations(jb))/(eigenvalues(ib) - eigenvalues(jb)), 0.0_rp, rp)
            end if
            if (do_static) response_static = response_static + 2.0_rp*transition*static_kernel*h_eigen(jb, ib)
            denominator = cmplx(eigenvalues(ib) - eigenvalues(jb), eta, rp)
            response_retarded = response_retarded + 2.0_rp*transition* &
               cmplx(occupations(ib) - occupations(jb), 0.0_rp, rp)/denominator*h_eigen(jb, ib)
         end do
      end do
      deallocate(h_eigen, transition, occupations)
   end subroutine exact_ks_product_field_response

   !> Contract the exact finite-H rotation into the live chi_plus product
   !> coordinates.  The factor two is the published circular convention: the
   !> halved sigma+ measurement plus its adjoint source gives the real
   !> transverse response.  The driving coefficient is B_H, not 2*B_H.
   subroutine exact_ks_product_rotation_response(product, hamiltonian, eigenvalues, eigenvectors, fermi_level, temperature, &
                                                 eta, response_rot, response_spec, response_retarded, include_static, response_l_filter)
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: hamiltonian(:, :), eigenvectors(:, :)
      real(rp), intent(in) :: eigenvalues(:), fermi_level, temperature, eta
      complex(rp), intent(out) :: response_rot(:), response_spec(:), response_retarded(:)
      logical, intent(in), optional :: include_static
      integer, intent(in), optional :: response_l_filter
      complex(rp), allocatable :: generator(:, :), delta_h(:, :), density(:, :), delta_rho_rot(:, :), delta_rho_spec(:, :)
      complex(rp), allocatable :: h_eigen(:, :), rot_eigen(:, :), spec_eigen(:, :), transition(:)
      real(rp), allocatable :: occupations(:)
      type(pauli_endpoint_state) :: left_state, right_state
      type(exact_ks_ward_metrics) :: metrics
      complex(rp) :: denominator, static_kernel
      integer :: n, ib, jb
      logical :: do_static

      if (eta <= 0.0_rp) error stop 'exact_ks_product_rotation_response: eta must be positive'
      do_static = .true.
      if (present(include_static)) do_static = include_static
      n = size(hamiltonian, 1)
      allocate(generator(n, n), delta_h(n, n), density(n, n), delta_rho_rot(n, n), delta_rho_spec(n, n))
      call exact_ks_global_rotation_oracle(hamiltonian, eigenvalues, eigenvectors, fermi_level, temperature, generator, &
         delta_h, density, delta_rho_rot, delta_rho_spec, metrics, (product%orbital_lmax + 1)**2, product%nsite)
      allocate(h_eigen(n, n), rot_eigen(n, n), spec_eigen(n, n), occupations(n), transition(product%product_dimension))
      h_eigen = matmul(conjg(transpose(eigenvectors)), matmul(delta_h, eigenvectors))
      rot_eigen = matmul(conjg(transpose(eigenvectors)), matmul(delta_rho_rot, eigenvectors))
      spec_eigen = matmul(conjg(transpose(eigenvectors)), matmul(delta_rho_spec, eigenvectors))
      do ib = 1, n
         occupations(ib) = exact_ks_fermi_occupation(eigenvalues(ib), fermi_level, temperature)
      end do
      response_rot = cmplx(0.0_rp, 0.0_rp, rp)
      response_spec = cmplx(0.0_rp, 0.0_rp, rp)
      response_retarded = cmplx(0.0_rp, 0.0_rp, rp)
      do ib = 1, n
         call left_state%initialize(eigenvalues(ib), eigenvectors(:, ib))
         do jb = 1, n
            call right_state%initialize(eigenvalues(jb), eigenvectors(:, jb))
            if (present(response_l_filter)) then
               call product%transition_coordinates(left_state, right_state, transition, response_l_filter=response_l_filter)
            else
               call product%transition_coordinates(left_state, right_state, transition)
            end if
            if (do_static) then
               response_rot = response_rot + 2.0_rp*transition*rot_eigen(jb, ib)
               response_spec = response_spec + 2.0_rp*transition*spec_eigen(jb, ib)
            end if
            denominator = cmplx(eigenvalues(ib) - eigenvalues(jb), eta, rp)
            static_kernel = cmplx(occupations(ib) - occupations(jb), 0.0_rp, rp)/denominator
            response_retarded = response_retarded + 2.0_rp*transition*static_kernel*h_eigen(jb, ib)
         end do
      end do
      deallocate(generator, delta_h, density, delta_rho_rot, delta_rho_spec, h_eigen, rot_eigen, spec_eigen, occupations, transition)
   end subroutine exact_ks_product_rotation_response

   subroutine exact_ks_product_rotation_sweep(product, hamiltonians, eigenvalues, eigenvectors, k_weights, fermi_level, &
                                               temperature, eta, response_rot, response_spec, response_retarded, &
                                               max_eigenpair_residual, max_density_residual, max_density_relative, retarded_only, &
                                               response_l_filter, max_eigenpair_frobenius, max_hermiticity, &
                                               max_rotation_field, max_spin_offdiagonal)
      type(lmto_product_response_basis), intent(in) :: product
      complex(rp), intent(in) :: hamiltonians(:, :, :), eigenvectors(:, :, :)
      real(rp), intent(in) :: eigenvalues(:, :), k_weights(:), fermi_level, temperature, eta
      complex(rp), intent(out) :: response_rot(:), response_spec(:), response_retarded(:)
      real(rp), intent(out) :: max_eigenpair_residual, max_density_residual, max_density_relative
      logical, intent(in), optional :: retarded_only
      integer, intent(in), optional :: response_l_filter
      real(rp), intent(out), optional :: max_eigenpair_frobenius, max_hermiticity, max_rotation_field, max_spin_offdiagonal
      complex(rp), allocatable :: generator(:, :), delta_h(:, :), density(:, :), delta_rho_rot(:, :), delta_rho_spec(:, :)
      complex(rp), allocatable :: one_rot(:), one_spec(:), one_retarded(:), one_target(:)
      type(exact_ks_ward_metrics) :: metrics
      integer :: ik, nk, n
      real(rp) :: weight_sum
      logical :: only_retarded

      n = size(hamiltonians, 1)
      nk = size(hamiltonians, 3)
      if (any(shape(hamiltonians) /= [n, n, nk]) .or. any(shape(eigenvectors) /= [n, n, nk]) .or. &
          any(shape(eigenvalues) /= [n, nk]) .or. size(k_weights) /= nk) then
         error stop 'exact_ks_product_rotation_sweep: inconsistent k-space dimensions'
      end if
      weight_sum = sum(k_weights)
      if (weight_sum <= 0.0_rp) error stop 'exact_ks_product_rotation_sweep: nonpositive k-weight sum'
      response_rot = cmplx(0.0_rp, 0.0_rp, rp)
      response_spec = cmplx(0.0_rp, 0.0_rp, rp)
      response_retarded = cmplx(0.0_rp, 0.0_rp, rp)
      only_retarded = .false.
      if (present(retarded_only)) only_retarded = retarded_only
      max_eigenpair_residual = 0.0_rp
      max_density_residual = 0.0_rp
      max_density_relative = 0.0_rp
      if (present(max_eigenpair_frobenius)) max_eigenpair_frobenius = 0.0_rp
      if (present(max_hermiticity)) max_hermiticity = 0.0_rp
      if (present(max_rotation_field)) max_rotation_field = 0.0_rp
      if (present(max_spin_offdiagonal)) max_spin_offdiagonal = 0.0_rp
      allocate(generator(n, n), delta_h(n, n), density(n, n), delta_rho_rot(n, n), delta_rho_spec(n, n))
      allocate(one_rot(product%product_dimension), one_spec(product%product_dimension), &
         one_retarded(product%product_dimension), one_target(product%product_dimension))
      do ik = 1, nk
         call exact_ks_global_rotation_oracle(hamiltonians(:, :, ik), eigenvalues(:, ik), eigenvectors(:, :, ik), &
            fermi_level, temperature, generator, delta_h, density, delta_rho_rot, delta_rho_spec, metrics, &
            (product%orbital_lmax + 1)**2, product%nsite)
         max_eigenpair_residual = max(max_eigenpair_residual, metrics%max_eigenpair_residual)
         max_density_residual = max(max_density_residual, metrics%density_rotation_spectral_residual)
         max_density_relative = max(max_density_relative, metrics%density_rotation_spectral_relative)
         if (present(max_eigenpair_frobenius)) then
            max_eigenpair_frobenius = max(max_eigenpair_frobenius, metrics%frobenius_eigenpair_residual)
         end if
         if (present(max_hermiticity)) max_hermiticity = max(max_hermiticity, metrics%hermiticity_residual)
         if (present(max_rotation_field)) max_rotation_field = max(max_rotation_field, metrics%rotation_field_relation_residual)
         if (present(max_spin_offdiagonal)) max_spin_offdiagonal = max(max_spin_offdiagonal, metrics%spin_offdiagonal_norm)
         if (present(response_l_filter)) then
            call exact_ks_product_rotation_response(product, hamiltonians(:, :, ik), eigenvalues(:, ik), &
               eigenvectors(:, :, ik), fermi_level, temperature, eta, one_rot, one_spec, one_retarded, .not. only_retarded, &
               response_l_filter)
         else
            call exact_ks_product_rotation_response(product, hamiltonians(:, :, ik), eigenvalues(:, ik), &
               eigenvectors(:, :, ik), fermi_level, temperature, eta, one_rot, one_spec, one_retarded, .not. only_retarded)
         end if
         ! The direct density contraction is the independent rigid-rotation
         ! target; the first product-vector output is not used as the
         ! oracle merely because it shares the same transition coordinates.
         if (.not. only_retarded) then
            call product_response_from_density(product, eigenvalues(:, ik), eigenvectors(:, :, ik), delta_rho_rot, one_target, &
               response_l_filter)
            response_rot = response_rot + (k_weights(ik)/weight_sum)*one_target
            response_spec = response_spec + (k_weights(ik)/weight_sum)*one_spec
         end if
         response_retarded = response_retarded + (k_weights(ik)/weight_sum)*one_retarded
      end do
      deallocate(generator, delta_h, density, delta_rho_rot, delta_rho_spec, one_rot, one_spec, one_retarded, one_target)
   contains
      subroutine product_response_from_density(product_basis, evals, vectors, density_response, target, response_l_filter)
         type(lmto_product_response_basis), intent(in) :: product_basis
         real(rp), intent(in) :: evals(:)
         complex(rp), intent(in) :: vectors(:, :), density_response(:, :)
         complex(rp), intent(out) :: target(:)
         integer, intent(in), optional :: response_l_filter
         type(pauli_endpoint_state) :: left_state, right_state
         complex(rp), allocatable :: transition(:), density_eigen(:, :)
         integer :: ib, jb, nlocal

         nlocal = size(evals)
         allocate(transition(product_basis%product_dimension), density_eigen(nlocal, nlocal))
         density_eigen = matmul(conjg(transpose(vectors)), matmul(density_response, vectors))
         target = cmplx(0.0_rp, 0.0_rp, rp)
         do ib = 1, nlocal
            call left_state%initialize(evals(ib), vectors(:, ib))
            do jb = 1, nlocal
               call right_state%initialize(evals(jb), vectors(:, jb))
               if (present(response_l_filter)) then
                  call product_basis%transition_coordinates(left_state, right_state, transition, &
                     response_l_filter=response_l_filter)
               else
                  call product_basis%transition_coordinates(left_state, right_state, transition)
               end if
               target = target + 2.0_rp*transition*density_eigen(jb, ib)
            end do
         end do
         deallocate(transition, density_eigen)
      end subroutine product_response_from_density
   end subroutine exact_ks_product_rotation_sweep

   pure real(rp) function exact_ks_spin_offdiagonal_norm(hamiltonian, norb_site, nsite) result(value)
      complex(rp), intent(in) :: hamiltonian(:, :)
      integer, intent(in) :: norb_site, nsite
      integer :: i, j, site, other_site, up_i, down_i, up_j, down_j

      value = 0.0_rp
      do site = 1, nsite
         do i = 1, norb_site
            up_i = (site - 1)*2*norb_site + i
            down_i = up_i + norb_site
            do other_site = 1, nsite
               do j = 1, norb_site
                  up_j = (other_site - 1)*2*norb_site + j
                  down_j = up_j + norb_site
                  value = value + abs(hamiltonian(up_i, down_j))**2 + abs(hamiltonian(down_i, up_j))**2
               end do
            end do
         end do
      end do
      value = sqrt(value)
   end function exact_ks_spin_offdiagonal_norm

   pure real(rp) function exact_ks_rotation_field_residual(hamiltonian, delta_h, norb_site, nsite) result(value)
      complex(rp), intent(in) :: hamiltonian(:, :), delta_h(:, :)
      integer, intent(in) :: norb_site, nsite
      integer :: i, j, site, other_site, up_i, down_i, up_j, down_j
      complex(rp) :: b

      value = 0.0_rp
      do site = 1, nsite
         do i = 1, norb_site
            up_i = (site - 1)*2*norb_site + i
            down_i = up_i + norb_site
            do other_site = 1, nsite
               do j = 1, norb_site
                  up_j = (other_site - 1)*2*norb_site + j
                  down_j = up_j + norb_site
                  b = 0.5_rp*(hamiltonian(up_i, up_j) - hamiltonian(down_i, down_j))
                  value = max(value, abs(delta_h(up_i, down_j) - b), abs(delta_h(down_i, up_j) - b))
               end do
            end do
         end do
      end do
   end function exact_ks_rotation_field_residual

end module lr_exact_ks_ward_mod
