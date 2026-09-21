!------------------------------------------------------------------------------
! DRESP-10 exact static product response.
!
! This module is intentionally independent of the assembled compact chi0
! services.  It evaluates the Frechet derivative of f(H) in the accepted
! finite-temperature occupation convention and contracts the resulting density
! matrix with the live six-branch product observable.
!------------------------------------------------------------------------------
module lr_static_frechet_product_response_mod

   use precision_mod, only: rp
   use lr_exact_ks_ward_mod, only: exact_ks_fermi_occupation, exact_ks_fermi_derivative
   use lr_lmto_product_response_basis_mod, only: lmto_product_response_basis
   use lr_pauli_transition_vertex_mod, only: pauli_endpoint_state
   implicit none
   private

   real(rp), parameter, public :: lr_static_default_degeneracy_tolerance = 1.0e-11_rp

   public :: lr_static_divided_difference
   public :: lr_static_frechet_density
   public :: lr_static_product_action
   public :: lr_static_product_matrix

contains

   !> Finite-temperature Frechet/divided-difference coefficient.
   pure real(rp) function lr_static_divided_difference(left_energy, right_energy, fermi_level, temperature, &
                                                       degeneracy_tolerance) result(value)
      real(rp), intent(in) :: left_energy, right_energy, fermi_level, temperature
      real(rp), intent(in), optional :: degeneracy_tolerance
      real(rp) :: tolerance

      tolerance = lr_static_default_degeneracy_tolerance
      if (present(degeneracy_tolerance)) tolerance = degeneracy_tolerance
      if (abs(left_energy - right_energy) <= tolerance*max(1.0_rp, abs(left_energy), abs(right_energy))) then
         value = exact_ks_fermi_derivative(left_energy, fermi_level, temperature)
      else
         value = (exact_ks_fermi_occupation(left_energy, fermi_level, temperature) - &
            exact_ks_fermi_occupation(right_energy, fermi_level, temperature))/(left_energy - right_energy)
      end if
   end function lr_static_divided_difference

   !> Direct coefficient-space Frechet action delta rho = L_f(H)[delta H].
   subroutine lr_static_frechet_density(eigenvalues, eigenvectors, fermi_level, temperature, delta_h, delta_rho, &
                                         degeneracy_tolerance)
      real(rp), intent(in) :: eigenvalues(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :), delta_h(:, :)
      complex(rp), intent(out) :: delta_rho(:, :)
      real(rp), intent(in), optional :: degeneracy_tolerance
      complex(rp), allocatable :: delta_eigen(:, :), delta_rho_eigen(:, :)
      integer :: n, i, j

      n = size(eigenvalues)
      if (size(eigenvectors, 1) /= n .or. size(eigenvectors, 2) /= n .or. any(shape(delta_h) /= [n, n]) .or. &
          any(shape(delta_rho) /= [n, n])) then
         error stop 'lr_static_frechet_density: inconsistent eigensystem or perturbation shape'
      end if
      allocate(delta_eigen(n, n), delta_rho_eigen(n, n))
      delta_eigen = matmul(conjg(transpose(eigenvectors)), matmul(delta_h, eigenvectors))
      do i = 1, n
         do j = 1, n
            delta_rho_eigen(i, j) = cmplx(lr_static_divided_difference(eigenvalues(i), eigenvalues(j), &
               fermi_level, temperature, degeneracy_tolerance), 0.0_rp, rp)*delta_eigen(i, j)
         end do
      end do
      delta_rho = matmul(eigenvectors, matmul(delta_rho_eigen, conjg(transpose(eigenvectors))))
      deallocate(delta_eigen, delta_rho_eigen)
   end subroutine lr_static_frechet_density

   !> Independent static action oracle.  No compact susceptibility matrix is
   !> formed or used here.
   subroutine lr_static_product_action(product, eigenvalues, eigenvectors, fermi_level, temperature, delta_h, response, &
                                       selected_l, degeneracy_tolerance)
      type(lmto_product_response_basis), intent(in) :: product
      real(rp), intent(in) :: eigenvalues(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :), delta_h(:, :)
      complex(rp), intent(out) :: response(:)
      logical, intent(in), optional :: selected_l(0:)
      real(rp), intent(in), optional :: degeneracy_tolerance
      complex(rp), allocatable :: delta_rho(:, :), delta_rho_eigen(:, :), transition(:)
      type(pauli_endpoint_state) :: left_state, right_state
      integer :: n, ib, jb

      n = size(eigenvalues)
      if (size(response) /= product%product_dimension) error stop 'lr_static_product_action: response shape mismatch'
      if (size(eigenvectors, 1) /= n .or. size(eigenvectors, 2) /= n .or. any(shape(delta_h) /= [n, n])) then
         error stop 'lr_static_product_action: eigensystem or perturbation shape mismatch'
      end if
      allocate(delta_rho(n, n), delta_rho_eigen(n, n), transition(product%product_dimension))
      call lr_static_frechet_density(eigenvalues, eigenvectors, fermi_level, temperature, delta_h, delta_rho, &
         degeneracy_tolerance)
      delta_rho_eigen = matmul(conjg(transpose(eigenvectors)), matmul(delta_rho, eigenvectors))
      response = cmplx(0.0_rp, 0.0_rp, rp)
      do ib = 1, n
         call left_state%initialize(eigenvalues(ib), eigenvectors(:, ib))
         do jb = 1, n
            call right_state%initialize(eigenvalues(jb), eigenvectors(:, jb))
            if (present(selected_l)) then
               call product%transition_coordinates(left_state, right_state, transition, selected_l)
            else
               call product%transition_coordinates(left_state, right_state, transition)
            end if
            response = response + 2.0_rp*transition*delta_rho_eigen(jb, ib)
         end do
      end do
      deallocate(delta_rho, delta_rho_eigen, transition)
   end subroutine lr_static_product_action

   !> Assemble the exact static product matrix from transition amplitudes and
   !> divided differences.  This is deliberately a separate implementation
   !> from both the retarded Lehmann and reciprocal-GF accumulators.
   subroutine lr_static_product_matrix(product, eigenvalues, eigenvectors, k_weights, fermi_level, temperature, matrix, &
                                       degeneracy_tolerance)
      type(lmto_product_response_basis), intent(in) :: product
      real(rp), intent(in) :: eigenvalues(:, :), k_weights(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      complex(rp), intent(out) :: matrix(:, :)
      real(rp), intent(in), optional :: degeneracy_tolerance
      complex(rp), allocatable :: transition(:)
      type(pauli_endpoint_state) :: left_state, right_state
      real(rp) :: weight_sum, dd
      integer :: n, nk, ik, ib, jb, i, j

      nk = size(eigenvalues, 2)
      n = size(eigenvalues, 1)
      if (size(k_weights) /= nk .or. size(eigenvectors, 1) /= n .or. size(eigenvectors, 2) /= n .or. &
          size(eigenvectors, 3) /= nk .or. any(shape(matrix) /= [product%product_dimension, product%product_dimension])) then
         error stop 'lr_static_product_matrix: inconsistent k-space or product shapes'
      end if
      weight_sum = sum(k_weights)
      if (weight_sum <= 0.0_rp) error stop 'lr_static_product_matrix: nonpositive k-weight sum'
      allocate(transition(product%product_dimension))
      matrix = cmplx(0.0_rp, 0.0_rp, rp)
      do ik = 1, nk
         do ib = 1, n
            call left_state%initialize(eigenvalues(ib, ik), eigenvectors(:, ib, ik))
            do jb = 1, n
               call right_state%initialize(eigenvalues(jb, ik), eigenvectors(:, jb, ik))
               call product%transition_coordinates(left_state, right_state, transition)
               dd = lr_static_divided_difference(eigenvalues(ib, ik), eigenvalues(jb, ik), fermi_level, temperature, &
                  degeneracy_tolerance)
               do j = 1, product%product_dimension
                  do i = 1, product%product_dimension
                     matrix(i, j) = matrix(i, j) + (2.0_rp*k_weights(ik)/weight_sum)*dd*transition(i)*conjg(transition(j))
                  end do
               end do
            end do
         end do
      end do
      deallocate(transition)
   end subroutine lr_static_product_matrix

end module lr_static_frechet_product_response_mod
