! DRESP-10R endpoint-power matrix-action oracle.
!
! The fixture uses a dense Hermitian H with nontrivial eigenvectors and checks
! every retained branch against the band-state definition.  The 20/02 checks
! are explicit so a regression to one H multiplication cannot pass through a
! diagonal or commuting fixture.
program test_lr_lmto_endpoint_branch_action
   use precision_mod, only: rp
   use lr_lmto_endpoint_branches_mod, only: lmto_product_nbranch, lmto_product_branch_powers, &
      lmto_product_apply_branch_action
   implicit none

   integer, parameter :: n = 5
   real(rp), parameter :: tolerance = 2.0e-12_rp
   real(rp) :: eigenvalues(n), pi, error(lmto_product_nbranch), stored_dual_error(lmto_product_nbranch)
   complex(rp) :: unitary(n,n), diagonal(n,n), hamiltonian(n,n), vertex(n,n), transformed(n,n), action(n,n)
   complex(rp) :: expected(n,n), stored_expected(n,n), stored_action(n,n)
   integer :: i, j, branch, p, q

   pi = acos(-1.0_rp)
   eigenvalues = [-0.83_rp, -0.21_rp, 0.17_rp, 0.64_rp, 1.13_rp]
   do j = 1, n
      do i = 1, n
         unitary(i,j) = cmplx(cos(2.0_rp*pi*real((i-1)*(j-1),rp)/real(n,rp)), &
            sin(2.0_rp*pi*real((i-1)*(j-1),rp)/real(n,rp)), rp)/sqrt(real(n,rp))
      end do
   end do
   diagonal = cmplx(0.0_rp, 0.0_rp, rp)
   do i = 1, n
      diagonal(i,i) = cmplx(eigenvalues(i), 0.0_rp, rp)
   end do
   hamiltonian = matmul(unitary, matmul(diagonal, conjg(transpose(unitary))))

   do j = 1, n
      do i = 1, n
         vertex(i,j) = cmplx(0.017_rp*real(3*i+j,rp), -0.011_rp*real(i+2*j,rp), rp)
      end do
   end do
   transformed = matmul(conjg(transpose(unitary)), matmul(vertex, unitary))

   error = 0.0_rp
   stored_dual_error = 0.0_rp
   do branch = 1, lmto_product_nbranch
      call lmto_product_branch_powers(branch, p, q)
      call lmto_product_apply_branch_action(hamiltonian, vertex, branch, action)
      expected = cmplx(0.0_rp, 0.0_rp, rp)
      do j = 1, n
         do i = 1, n
            expected(i,j) = eigenvalues(i)**p*eigenvalues(j)**q*transformed(i,j)
         end do
      end do
      expected = matmul(unitary, matmul(expected, conjg(transpose(unitary))))
      error(branch) = frobenius_relative(action, expected)

      call lmto_product_apply_branch_action(hamiltonian, vertex, branch, stored_action, .true.)
      stored_expected = cmplx(0.0_rp, 0.0_rp, rp)
      do j = 1, n
         do i = 1, n
            stored_expected(i,j) = eigenvalues(i)**q*eigenvalues(j)**p*transformed(i,j)
         end do
      end do
      stored_expected = matmul(unitary, matmul(stored_expected, conjg(transpose(unitary))))
      stored_dual_error(branch) = frobenius_relative(stored_action, stored_expected)
   end do

   write (*, '(a,6(es12.4,1x))') '  direct_00_10_01_11_20_02=', error
   write (*, '(a,6(es12.4,1x))') '  stored_dual_00_10_01_11_20_02=', stored_dual_error
   if (maxval(error) > tolerance .or. maxval(stored_dual_error) > tolerance .or. &
       error(5) > tolerance .or. error(6) > tolerance) then
      error stop 'UnitLrLmtoEndpointBranchAction: FAIL'
   end if
   write (*, '(a)') 'UnitLrLmtoEndpointBranchAction: PASS (all six powers; explicit H2 20/02; stored-dual orientation)'

contains

   pure real(rp) function frobenius_relative(left, right) result(value)
      complex(rp), intent(in) :: left(:, :), right(:, :)
      value = sqrt(sum(abs(left-right)**2))/max(sqrt(sum(abs(right)**2)), tiny(1.0_rp))
   end function frobenius_relative

end program test_lr_lmto_endpoint_branch_action
