!------------------------------------------------------------------------------
! DRESP-09V production energy-moment density tangent.
!
! The production reciprocal density is not defined by an unspecified
! coefficient transformation.  It is the three-moment object
!
!     M0 = rho, M1 = H rho, M2 = H**2 rho,
!
! accumulated into spin_density after the k-point/eigenstate sum.  This
! module contains the small, representation-explicit algebra used by the
! DRESP-09V bridge.  In particular, endpoint_tangent_branches differentiates
! the complete D_pq = H**p rho H**q object; the older DRESP-09S fixed-H route
! remains a separate oracle.
!------------------------------------------------------------------------------
module lr_lmto_density_moment_tangent_mod

   use precision_mod, only: rp
   use spin_density_mod, only: spin_density, sd_orders
   implicit none
   private

   real(rp), parameter, public :: density_moment_kb_ry_per_kelvin = 6.3336814e-6_rp

   public :: density_from_eigensystem
   public :: moment_matrices_from_eigensystem
   public :: moment_tangent_product_rule
   public :: moment_tangent_frechet
   public :: endpoint_tangent_branches
   public :: endpoint_tangent_branches_second_order
   public :: endpoint_fixed_h_branches
   public :: endpoint_fixed_h_branches_second_order
   public :: endpoint_matrices_from_moments
   public :: endpoint_matrices_from_moments_second_order
   public :: commutator_tangent
   public :: fill_production_spin_density_from_moments
   public :: relative_matrix_residual

contains

   !> Build rho=f(H) from the accepted standard Hermitian eigensystem.
   subroutine density_from_eigensystem(eigenvalues, eigenvectors, fermi_level, temperature, density)
      real(rp), intent(in) :: eigenvalues(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :)
      complex(rp), intent(out) :: density(:, :)
      integer :: n, ib, j
      real(rp) :: occupation

      n = size(eigenvalues)
      if (any(shape(eigenvectors) /= [n, n]) .or. any(shape(density) /= [n, n])) then
         error stop 'DRESP-09V density: eigensystem shape mismatch'
      end if
      density = cmplx(0.0_rp, 0.0_rp, rp)
      do ib = 1, n
         occupation = fermi_occupation(eigenvalues(ib), fermi_level, temperature)
         if (occupation <= 1.0e-14_rp) cycle
         do j = 1, n
            density(:, j) = density(:, j) + occupation*eigenvectors(:, ib)*conjg(eigenvectors(j, ib))
         end do
      end do
   end subroutine density_from_eigensystem

   !> Fill M0/M1/M2 from the production eigenvector definition for one k point.
   subroutine moment_matrices_from_eigensystem(eigenvalues, eigenvectors, fermi_level, temperature, moments)
      real(rp), intent(in) :: eigenvalues(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :)
      complex(rp), intent(out) :: moments(:, :, :)
      integer :: n, ib, j, iorder
      real(rp) :: occupation, epow

      n = size(eigenvalues)
      if (any(shape(eigenvectors) /= [n, n]) .or. any(shape(moments) /= [n, n, sd_orders])) then
         error stop 'DRESP-09V moments: eigensystem shape mismatch'
      end if
      moments = cmplx(0.0_rp, 0.0_rp, rp)
      do ib = 1, n
         occupation = fermi_occupation(eigenvalues(ib), fermi_level, temperature)
         if (occupation <= 1.0e-14_rp) cycle
         epow = 1.0_rp
         do iorder = 1, sd_orders
            do j = 1, n
               moments(:, j, iorder) = moments(:, j, iorder) + occupation*epow*eigenvectors(:, ib)* &
                  conjg(eigenvectors(j, ib))
            end do
            epow = epow*eigenvalues(ib)
         end do
      end do
   end subroutine moment_matrices_from_eigensystem

   !> Complete product-rule tangent of M0=rho, M1=H rho, M2=H^2 rho.
   subroutine moment_tangent_product_rule(hamiltonian, density, delta_h, delta_density, delta_moments)
      complex(rp), intent(in) :: hamiltonian(:, :), density(:, :), delta_h(:, :), delta_density(:, :)
      complex(rp), intent(out) :: delta_moments(:, :, :)
      integer :: n
      complex(rp), allocatable :: h2(:, :)

      n = size(hamiltonian, 1)
      if (any(shape(hamiltonian) /= [n, n]) .or. any(shape(density) /= [n, n]) .or. &
          any(shape(delta_h) /= [n, n]) .or. any(shape(delta_density) /= [n, n]) .or. &
          any(shape(delta_moments) /= [n, n, sd_orders])) then
         error stop 'DRESP-09V product tangent: shape mismatch'
      end if
      allocate(h2(n, n))
      h2 = matmul(hamiltonian, hamiltonian)
      delta_moments(:, :, 1) = delta_density
      delta_moments(:, :, 2) = matmul(delta_h, density) + matmul(hamiltonian, delta_density)
      delta_moments(:, :, 3) = matmul(delta_h, matmul(hamiltonian, density)) + &
         matmul(hamiltonian, matmul(delta_h, density)) + matmul(h2, delta_density)
      deallocate(h2)
   end subroutine moment_tangent_product_rule

   !> Independent Frechet tangent of g_p(H)=H^p f(H), p=0,1,2.
   !> The divided difference uses the analytic g_p derivative in the
   !> equal/near-degenerate limit.
   subroutine moment_tangent_frechet(eigenvalues, eigenvectors, delta_h, fermi_level, temperature, delta_moments)
      real(rp), intent(in) :: eigenvalues(:), fermi_level, temperature
      complex(rp), intent(in) :: eigenvectors(:, :), delta_h(:, :)
      complex(rp), intent(out) :: delta_moments(:, :, :)
      complex(rp), allocatable :: delta_eigen(:, :), tangent_eigen(:, :), work(:, :)
      real(rp) :: denominator, scale, kernel
      integer :: n, i, j, p

      n = size(eigenvalues)
      if (any(shape(eigenvectors) /= [n, n]) .or. any(shape(delta_h) /= [n, n]) .or. &
          any(shape(delta_moments) /= [n, n, sd_orders])) then
         error stop 'DRESP-09V Frechet tangent: shape mismatch'
      end if
      allocate(delta_eigen(n, n), tangent_eigen(n, n), work(n, n))
      delta_eigen = matmul(conjg(transpose(eigenvectors)), matmul(delta_h, eigenvectors))
      do p = 0, sd_orders - 1
         do i = 1, n
            do j = 1, n
               denominator = eigenvalues(i) - eigenvalues(j)
               scale = max(1.0_rp, abs(eigenvalues(i)), abs(eigenvalues(j)))
               if (abs(denominator) <= 1.0e-11_rp*scale) then
                  kernel = gp_derivative(p, eigenvalues(i), fermi_level, temperature)
               else
                  kernel = (gp_value(p, eigenvalues(i), fermi_level, temperature) - &
                            gp_value(p, eigenvalues(j), fermi_level, temperature))/denominator
               end if
               tangent_eigen(i, j) = cmplx(kernel, 0.0_rp, rp)*delta_eigen(i, j)
            end do
         end do
         work = matmul(eigenvectors, matmul(tangent_eigen, conjg(transpose(eigenvectors))))
         delta_moments(:, :, p + 1) = work
      end do
      deallocate(delta_eigen, tangent_eigen, work)
   end subroutine moment_tangent_frechet

   !> Complete tangent of D_pq=H^p rho H^q in branch order 00,10,01,11.
   subroutine endpoint_tangent_branches(hamiltonian, density, delta_h, delta_density, delta_branches)
      complex(rp), intent(in) :: hamiltonian(:, :), density(:, :), delta_h(:, :), delta_density(:, :)
      complex(rp), intent(out) :: delta_branches(:, :, :)
      integer :: n

      n = size(hamiltonian, 1)
      if (any(shape(hamiltonian) /= [n, n]) .or. any(shape(density) /= [n, n]) .or. &
          any(shape(delta_h) /= [n, n]) .or. any(shape(delta_density) /= [n, n]) .or. &
          any(shape(delta_branches) /= [n, n, 4])) then
         error stop 'DRESP-09V endpoint tangent: shape mismatch'
      end if
      delta_branches(:, :, 1) = delta_density
      delta_branches(:, :, 2) = matmul(delta_h, density) + matmul(hamiltonian, delta_density)
      delta_branches(:, :, 3) = matmul(delta_density, hamiltonian) + matmul(density, delta_h)
      delta_branches(:, :, 4) = matmul(delta_h, matmul(density, hamiltonian)) + &
         matmul(hamiltonian, matmul(delta_density, hamiltonian)) + &
         matmul(hamiltonian, matmul(density, delta_h))
   end subroutine endpoint_tangent_branches

   !> Complete second-order endpoint tangent in branch order
   !> 00,10,01,11,20,02.  The final two branches are deliberately written
   !> as product-rule expressions; they are not obtained by multiplying a
   !> generic endpoint component of power two.  This preserves the Taylor
   !> truncation of the radial observable.
   subroutine endpoint_tangent_branches_second_order(hamiltonian, density, delta_h, delta_density, delta_branches)
      complex(rp), intent(in) :: hamiltonian(:, :), density(:, :), delta_h(:, :), delta_density(:, :)
      complex(rp), intent(out) :: delta_branches(:, :, :)
      integer :: n
      complex(rp), allocatable :: h2(:, :)

      n = size(hamiltonian, 1)
      if (any(shape(hamiltonian) /= [n, n]) .or. any(shape(density) /= [n, n]) .or. &
          any(shape(delta_h) /= [n, n]) .or. any(shape(delta_density) /= [n, n]) .or. &
          any(shape(delta_branches) /= [n, n, 6])) then
         error stop 'DRESP-09X endpoint tangent: shape mismatch'
      end if
      allocate(h2(n, n))
      h2 = matmul(hamiltonian, hamiltonian)
      delta_branches(:, :, 1) = delta_density
      delta_branches(:, :, 2) = matmul(delta_h, density) + matmul(hamiltonian, delta_density)
      delta_branches(:, :, 3) = matmul(delta_density, hamiltonian) + matmul(density, delta_h)
      delta_branches(:, :, 4) = matmul(delta_h, matmul(density, hamiltonian)) + &
         matmul(hamiltonian, matmul(delta_density, hamiltonian)) + &
         matmul(hamiltonian, matmul(density, delta_h))
      delta_branches(:, :, 5) = matmul(delta_h, matmul(hamiltonian, density)) + &
         matmul(hamiltonian, matmul(delta_h, density)) + matmul(h2, delta_density)
      delta_branches(:, :, 6) = matmul(delta_density, h2) + &
         matmul(density, matmul(delta_h, hamiltonian)) + matmul(density, matmul(hamiltonian, delta_h))
      deallocate(h2)
   end subroutine endpoint_tangent_branches_second_order

   !> DRESP-09S fixed-H density derivative, retained as the comparison oracle.
   subroutine endpoint_fixed_h_branches(hamiltonian, delta_density, delta_branches)
      complex(rp), intent(in) :: hamiltonian(:, :), delta_density(:, :)
      complex(rp), intent(out) :: delta_branches(:, :, :)
      integer :: n

      n = size(hamiltonian, 1)
      if (any(shape(hamiltonian) /= [n, n]) .or. any(shape(delta_density) /= [n, n]) .or. &
          any(shape(delta_branches) /= [n, n, 4])) then
         error stop 'DRESP-09V fixed-H tangent: shape mismatch'
      end if
      delta_branches(:, :, 1) = delta_density
      delta_branches(:, :, 2) = matmul(hamiltonian, delta_density)
      delta_branches(:, :, 3) = matmul(delta_density, hamiltonian)
      delta_branches(:, :, 4) = matmul(hamiltonian, matmul(delta_density, hamiltonian))
   end subroutine endpoint_fixed_h_branches

   !> Frozen-endpoint density tangent for the complete second-order branch
   !> contract.  This is deliberately separate from the complete tangent:
   !> only delta_density is varied, while every Hamiltonian endpoint is held
   !> fixed.  Branch order is 00,10,01,11,20,02.
   subroutine endpoint_fixed_h_branches_second_order(hamiltonian, delta_density, delta_branches)
      complex(rp), intent(in) :: hamiltonian(:, :), delta_density(:, :)
      complex(rp), intent(out) :: delta_branches(:, :, :)
      integer :: n
      complex(rp), allocatable :: h2(:, :)

      n = size(hamiltonian, 1)
      if (any(shape(hamiltonian) /= [n, n]) .or. any(shape(delta_density) /= [n, n]) .or. &
          any(shape(delta_branches) /= [n, n, 6])) then
         error stop 'DRESP-12 frozen-H endpoint tangent: shape mismatch'
      end if
      allocate(h2(n, n))
      h2 = matmul(hamiltonian, hamiltonian)
      delta_branches(:, :, 1) = delta_density
      delta_branches(:, :, 2) = matmul(hamiltonian, delta_density)
      delta_branches(:, :, 3) = matmul(delta_density, hamiltonian)
      delta_branches(:, :, 4) = matmul(hamiltonian, matmul(delta_density, hamiltonian))
      delta_branches(:, :, 5) = matmul(h2, delta_density)
      delta_branches(:, :, 6) = matmul(delta_density, h2)
      deallocate(h2)
   end subroutine endpoint_fixed_h_branches_second_order

   !> Equilibrium endpoint objects from M0/M1/M2.
   subroutine endpoint_matrices_from_moments(moments, endpoints)
      complex(rp), intent(in) :: moments(:, :, :)
      complex(rp), intent(out) :: endpoints(:, :, :)
      integer :: n

      n = size(moments, 1)
      if (any(shape(moments) /= [n, n, sd_orders]) .or. any(shape(endpoints) /= [n, n, 4])) then
         error stop 'DRESP-09V endpoint moments: shape mismatch'
      end if
      endpoints(:, :, 1) = moments(:, :, 1)
      endpoints(:, :, 2) = moments(:, :, 2)
      endpoints(:, :, 3) = moments(:, :, 2)
      endpoints(:, :, 4) = moments(:, :, 3)
   end subroutine endpoint_matrices_from_moments

   !> Equilibrium six-branch endpoint representation from certified M0/M1/M2.
   !> In the accepted equilibrium state M2=H**2 rho=rho H**2 up to the
   !> independently audited [H,rho] residual, so both pure second-order
   !> endpoints are populated from the same production M2 moment.
   subroutine endpoint_matrices_from_moments_second_order(moments, endpoints)
      complex(rp), intent(in) :: moments(:, :, :)
      complex(rp), intent(out) :: endpoints(:, :, :)
      integer :: n

      n = size(moments, 1)
      if (any(shape(moments) /= [n, n, sd_orders]) .or. any(shape(endpoints) /= [n, n, 6])) then
         error stop 'DRESP-09X endpoint moments: shape mismatch'
      end if
      endpoints(:, :, 1) = moments(:, :, 1)
      endpoints(:, :, 2) = moments(:, :, 2)
      endpoints(:, :, 3) = moments(:, :, 2)
      endpoints(:, :, 4) = moments(:, :, 3)
      endpoints(:, :, 5) = moments(:, :, 3)
      endpoints(:, :, 6) = moments(:, :, 3)
   end subroutine endpoint_matrices_from_moments_second_order

   pure subroutine commutator_tangent(generator, matrix, delta_matrix)
      complex(rp), intent(in) :: generator(:, :), matrix(:, :)
      complex(rp), intent(out) :: delta_matrix(:, :)
      integer :: n

      n = size(matrix, 1)
      if (any(shape(generator) /= [n, n]) .or. any(shape(delta_matrix) /= [n, n])) then
         error stop 'DRESP-09V commutator tangent: shape mismatch'
      end if
      delta_matrix = cmplx(0.0_rp, -1.0_rp, rp)*(matmul(generator, matrix) - matmul(matrix, generator))
   end subroutine commutator_tangent

   !> Accumulate full coefficient-space moment matrices using the exact
   !> site/l grouping used by reciprocal_spin_density.f90.
   subroutine fill_production_spin_density_from_moments(moment_matrices, density)
      complex(rp), intent(in) :: moment_matrices(:, :, :)
      type(spin_density), intent(inout) :: density
      integer :: nsite, norb_site, n, site, il, io, iup, idn, iorder
      integer, parameter :: lstart(4) = [1, 2, 5, 10]
      integer, parameter :: lend(4) = [1, 4, 9, 16]
      complex(rp) :: block(2, 2)

      nsite = density%nsite
      if (density%nl < 4 .or. size(moment_matrices, 3) /= sd_orders) then
         error stop 'DRESP-09V spin_density mapping: incomplete channel dimensions'
      end if
      n = size(moment_matrices, 1)
      if (any(shape(moment_matrices) /= [n, n, sd_orders]) .or. mod(n, 2*nsite) /= 0) then
         error stop 'DRESP-09V spin_density mapping: matrix shape mismatch'
      end if
      norb_site = n/(2*nsite)
      call density%zero_density()
      do site = 1, nsite
         do il = 1, min(4, density%nl)
            if (lstart(il) > norb_site) cycle
            do iorder = 1, sd_orders
               block = cmplx(0.0_rp, 0.0_rp, rp)
               do io = lstart(il), min(lend(il), norb_site)
                  iup = (site - 1)*2*norb_site + io
                  idn = iup + norb_site
                  block(1, 1) = block(1, 1) + moment_matrices(iup, iup, iorder)
                  block(1, 2) = block(1, 2) + moment_matrices(iup, idn, iorder)
                  block(2, 1) = block(2, 1) + moment_matrices(idn, iup, iorder)
                  block(2, 2) = block(2, 2) + moment_matrices(idn, idn, iorder)
               end do
               call density%accumulate_block(site, il, iorder, block)
            end do
         end do
      end do
   end subroutine fill_production_spin_density_from_moments

   pure real(rp) function relative_matrix_residual(left, right) result(value)
      complex(rp), intent(in) :: left(:, :), right(:, :)
      value = sqrt(sum(abs(left - right)**2))/max(sqrt(sum(abs(right)**2)), tiny(1.0_rp))
   end function relative_matrix_residual

   pure real(rp) function fermi_occupation(energy, fermi_level, temperature) result(value)
      real(rp), intent(in) :: energy, fermi_level, temperature
      real(rp) :: argument, kt

      kt = max(temperature*density_moment_kb_ry_per_kelvin, 1.0e-12_rp)
      argument = (energy - fermi_level)/kt
      if (argument >= 50.0_rp) then
         value = 0.0_rp
      else if (argument <= -50.0_rp) then
         value = 1.0_rp
      else
         value = 1.0_rp/(exp(argument) + 1.0_rp)
      end if
   end function fermi_occupation

   pure real(rp) function fermi_derivative(energy, fermi_level, temperature) result(value)
      real(rp), intent(in) :: energy, fermi_level, temperature
      real(rp) :: occupation, kt

      kt = max(temperature*density_moment_kb_ry_per_kelvin, 1.0e-12_rp)
      occupation = fermi_occupation(energy, fermi_level, temperature)
      value = -occupation*(1.0_rp - occupation)/kt
   end function fermi_derivative

   pure real(rp) function gp_value(power, energy, fermi_level, temperature) result(value)
      integer, intent(in) :: power
      real(rp), intent(in) :: energy, fermi_level, temperature
      value = energy**power*fermi_occupation(energy, fermi_level, temperature)
   end function gp_value

   pure real(rp) function gp_derivative(power, energy, fermi_level, temperature) result(value)
      integer, intent(in) :: power
      real(rp), intent(in) :: energy, fermi_level, temperature
      real(rp) :: occupation, derivative

      occupation = fermi_occupation(energy, fermi_level, temperature)
      derivative = fermi_derivative(energy, fermi_level, temperature)
      select case (power)
      case (0)
         value = derivative
      case (1)
         value = occupation + energy*derivative
      case (2)
         value = 2.0_rp*energy*occupation + energy*energy*derivative
      case default
         error stop 'DRESP-09V Frechet derivative: unsupported moment power'
      end select
   end function gp_derivative

end module lr_lmto_density_moment_tangent_mod
