!> Space group found by symmetry%initialize from lattice%a (columns = primitive vectors) and lattice%cr (Cartesian, units of alat).
!> Expected groups are fixed by crystallography: hcp P6_3/mmc #194 (24 operations), the same cell with two species
!> P-6m2 #187 (12), diamond Si Fd-3m #227 (48).
program test_symmetry_init
   use lattice_mod, only: lattice
   use symmetry_mod, only: symmetry
   implicit none
   integer, parameter :: rp = kind(1.0d0)
   real(rp), parameter :: s3 = sqrt(3.0_rp), ca = sqrt(8.0_rp/3.0_rp)
   real(rp) :: hex(3, 3), fcc(3, 3)

   hex = reshape([1.0_rp, 0.0_rp, 0.0_rp, -0.5_rp, 0.5_rp*s3, 0.0_rp, 0.0_rp, 0.0_rp, ca], [3, 3])
   fcc = reshape([0.0_rp, 0.5_rp, 0.5_rp, 0.5_rp, 0.0_rp, 0.5_rp, 0.5_rp, 0.5_rp, 0.0_rp], [3, 3])

   ! Bug caught: Cartesian cr passed to spglib as fractional (hcp found as P1, 1 operation) or a transposed
   call check(hex, reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp/s3, 0.5_rp*ca], [3, 2]), [26, 26], 194, 24)
   call check(hex, reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp/s3, 0.5_rp*ca], [3, 2]), [1, 2], 187, 12)
   call check(fcc, reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.25_rp, 0.25_rp, 0.25_rp], [3, 2]), [14, 14], 227, 48)
   print '(A)', 'test_symmetry_init: all checks passed'

contains

   subroutine check(a, cr, iz, nexpected, nops_expected)
      real(rp), intent(in) :: a(3, 3), cr(:, :)
      integer, intent(in) :: iz(:), nexpected, nops_expected
      type(lattice), target :: lat
      type(symmetry) :: sym

      lat%alat = 5.0_rp
      lat%a = a
      lat%nrec = size(iz)
      allocate(lat%cr(3, lat%nrec), lat%iz(lat%nrec))
      lat%cr = cr
      lat%iz = iz
      call sym%initialize(lat)
      print '(A,I4,A,I4,A,I4,A,I4)', 'space group', sym%space_group_number, ' expected', nexpected, &
         '  operations', sym%num_symmetry_operations, ' expected', nops_expected
      if (sym%space_group_number /= nexpected .or. sym%num_symmetry_operations /= nops_expected) &
         error stop 'space group or operation count differs from crystallography'
   end subroutine check

end program test_symmetry_init
