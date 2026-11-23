!> Reciprocal basis of lattice%a, where lattice%a(:, i) is primitive vector i.
!> Definitional check: b_i . a_j = 2 pi delta_ij for hcp, bcc, fcc and a triclinic cell.
!> Independent check: hcp b_i against the textbook hexagonal basis written out below.
program test_reciprocal_vectors
   use precision_mod, only: rp
   use reciprocal_mod, only: reciprocal
   use lattice_mod, only: lattice
   implicit none
   real(rp), parameter :: pi = acos(-1.0_rp), s3 = sqrt(3.0_rp), c_a = 1.633_rp
   real(rp), parameter :: hcp(3, 3) = reshape([1.0_rp, 0.0_rp, 0.0_rp, -0.5_rp, 0.5_rp*s3, 0.0_rp, 0.0_rp, 0.0_rp, c_a], [3, 3])
   real(rp), parameter :: bcc(3, 3) = reshape([-0.5_rp, 0.5_rp, 0.5_rp, 0.5_rp, -0.5_rp, 0.5_rp, 0.5_rp, 0.5_rp, -0.5_rp], [3, 3])
   real(rp), parameter :: fcc(3, 3) = reshape([0.0_rp, 0.5_rp, 0.5_rp, 0.5_rp, 0.0_rp, 0.5_rp, 0.5_rp, 0.5_rp, 0.0_rp], [3, 3])
   real(rp), parameter :: tri(3, 3) = reshape([1.1_rp, 0.2_rp, -0.1_rp, 0.3_rp, 0.9_rp, 0.15_rp, -0.2_rp, 0.4_rp, 1.3_rp], [3, 3])
   real(rp) :: b(3, 3), bhex(3, 3)

   ! Bug caught: a read by rows instead of columns (b1 = (2 pi, 0, 0) for hcp; b_i . a_j /= 2 pi delta_ij)
   call dual_check(bcc, 'bcc', b)
   call dual_check(fcc, 'fcc', b)
   call dual_check(tri, 'triclinic', b)

   ! Textbook hexagonal basis for a1 = (1, 0, 0), a2 = (-1/2, sqrt(3)/2, 0), a3 = (0, 0, c/a)
   bhex(:, 1) = 2.0_rp*pi*[1.0_rp, 1.0_rp/s3, 0.0_rp]
   bhex(:, 2) = 2.0_rp*pi*[0.0_rp, 2.0_rp/s3, 0.0_rp]
   bhex(:, 3) = 2.0_rp*pi*[0.0_rp, 0.0_rp, 1.0_rp/c_a]
   call dual_check(hcp, 'hcp', b)
   print '(A,ES10.2)', 'hcp max |b - textbook| = ', maxval(abs(b - bhex))
   if (maxval(abs(b - bhex)) > 1.0e-12_rp) error stop 'hcp reciprocal basis differs from textbook'
   print '(A)', 'test_reciprocal_vectors: all checks passed'

contains

   subroutine dual_check(a, name, bout)
      real(rp), intent(in) :: a(3, 3)
      character(*), intent(in) :: name
      real(rp), intent(out) :: bout(3, 3)
      type(reciprocal) :: recip
      type(lattice), target :: lat
      integer :: i, j
      real(rp) :: err

      call recip%restore_to_default()
      lat%a = a
      recip%lattice => lat
      call recip%generate_reciprocal_vectors()
      bout = recip%reciprocal_vectors
      err = 0.0_rp
      do i = 1, 3
         do j = 1, 3
            err = max(err, abs(dot_product(bout(:, i), a(:, j)) - merge(2.0_rp*pi, 0.0_rp, i == j)))
         end do
      end do
      print '(A,A,ES10.2)', name, ' max |b_i.a_j - 2 pi delta_ij| = ', err
      if (err > 1.0e-12_rp) error stop 'b_i . a_j /= 2 pi delta_ij'
   end subroutine dual_check
end program test_reciprocal_vectors
