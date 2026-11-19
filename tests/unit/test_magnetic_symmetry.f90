!> Magnetic symmetry groups and the reduced k meshes built from them (spglib_interface).
!> Expected orders and mesh counts are fixed by the point groups (4/mm'm', -3m', m-3m) and by a brute-force
!> orbit count below, not by the code under test.
program test_magnetic_symmetry
   use spglib_interface_mod, only: spglib_interface
   implicit none
   integer, parameter :: rp = kind(1.0d0)
   real(rp), parameter :: s3 = 1.0_rp/sqrt(3.0_rp), s14 = 1.0_rp/sqrt(14.0_rp)
   real(rp), parameter :: bcc(3, 3) = reshape([-0.5_rp, 0.5_rp, 0.5_rp, 0.5_rp, -0.5_rp, 0.5_rp, 0.5_rp, 0.5_rp, -0.5_rp], [3, 3])
   real(rp), parameter :: fcc(3, 3) = reshape([0.0_rp, 0.5_rp, 0.5_rp, 0.5_rp, 0.0_rp, 0.5_rp, 0.5_rp, 0.5_rp, 0.0_rp], [3, 3])
   real(rp), parameter :: sc(3, 3) = reshape([1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp], [3, 3])
   type(spglib_interface) :: spg
   integer, allocatable :: rot(:, :, :), flips(:)
   real(rp), allocatable :: tr(:, :)
   real(rp) :: z0(3, 1) = 0.0_rp

   ! Bug caught: group order wrong for the moment direction (moments ignored or not treated as axial)
   call group(bcc, reshape([0.0_rp, 0.0_rp, 0.0_rp], [3, 1]), [26], reshape([0.0_rp, 0.0_rp, 1.0_rp], [3, 1]), .false., 48, 0)
   call group(bcc, reshape([0.0_rp, 0.0_rp, 0.0_rp], [3, 1]), [26], reshape([0.0_rp, 0.0_rp, 1.0_rp], [3, 1]), .true., 8, 8)
   call group(bcc, reshape([0.0_rp, 0.0_rp, 0.0_rp], [3, 1]), [26], reshape([s3, s3, s3], [3, 1]), .true., 6, 6)
   ! Bug caught: moment vector used as a species label wrongly (here one site: all 48 remain)
   call group(bcc, reshape([0.0_rp, 0.0_rp, 0.0_rp], [3, 1]), [26], reshape([1.0_rp, 2.0_rp, 3.0_rp]*s14, [3, 1]), .false., 48, 0)
   ! Bug caught: moment labels ignored in the no-SOC group (FM and AFM two-site sc cells must differ: 96 vs 48)
   call group(sc, reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.5_rp, 0.5_rp, 0.5_rp], [3, 2]), [26, 26], &
              reshape([0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 1.0_rp], [3, 2]), .false., 96, 0)
   call group(sc, reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.5_rp, 0.5_rp, 0.5_rp], [3, 2]), [26, 26], &
              reshape([0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, -1.0_rp], [3, 2]), .false., 48, 0)
   ! Diamond Si, zero moments: 48 unitary operations (spglib adds time-reversal copies, counted separately)
   call group(fcc, reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.25_rp, 0.25_rp, 0.25_rp], [3, 2]), [14, 14], &
              reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], [3, 2]), .true., 48, 48)

   ! Bug caught: time-reversal-combined operations not entered as -R (mesh too large), or R used in place of R^-T
   call mesh_count([0.0_rp, 0.0_rp, 1.0_rp], 4, 13)
   call mesh_count([0.0_rp, 0.0_rp, 1.0_rp], 10, 102)
   call mesh_count([s3, s3, s3], 4, 13)
   call mesh_count([s3, s3, s3], 10, 116)
   print '(A)', 'test_magnetic_symmetry: all checks passed'

contains

   !> Group order: unitary operations (flip = +1) and time-reversal-combined ones (flip = -1); totals 16 and 12 for the SOC bcc Fe cases.
   subroutine group(lat, pos, types, mom, soc, nu_expected, tr_expected)
      real(rp), intent(in) :: lat(3, 3), pos(:, :), mom(:, :)
      integer, intent(in) :: types(:), nu_expected, tr_expected
      logical, intent(in) :: soc
      integer :: n, ntr

      call spg%initialize(lat, pos, types)
      n = spg%get_magnetic_symmetry(pos, types, mom, soc, rot, tr, flips)
      ntr = count(flips < 0)
      print '(A,L2,A,I4,A,I4)', 'soc', soc, ' unitary', n - ntr, ' time-reversal-combined', ntr
      if (n - ntr /= nu_expected .or. ntr /= tr_expected) error stop 'magnetic group order differs from expected'
      call spg%finalize()
   end subroutine group

   !> Irreducible count of an nk^3 Gamma-centred mesh for bcc Fe SOC with moment mom, vs integer orbit count.
   subroutine mesh_count(mom, nk, expected)
      real(rp), intent(in) :: mom(3)
      integer, intent(in) :: nk, expected
      real(rp), allocatable :: kp(:, :), wt(:)
      integer, allocatable :: signed(:, :, :)
      integer :: n, nir, norb, ik, i, j, k, op, kv(3), kn(3), best
      real(rp) :: p(3, 1) = 0.0_rp

      call spg%initialize(bcc, p, [26])
      n = spg%get_magnetic_symmetry(p, [26], reshape(mom, [3, 1]), .true., rot, tr, flips)
      allocate(signed(3, 3, n))
      do op = 1, n
         signed(:, :, op) = flips(op)*rot(:, :, op)
      end do
      nir = spg%get_reduced_kpoint_mesh_with_points([nk, nk, nk], [0, 0, 0], kp, wt, .false., rotations=signed)

      ! orbits under k -> flip R^T k (mod nk), counted by the smallest linear index in each orbit
      norb = 0
      do ik = 0, nk**3 - 1
         kv = [ik/(nk*nk), mod(ik/nk, nk), mod(ik, nk)]
         best = ik
         do op = 1, n
            kn = modulo(matmul(transpose(signed(:, :, op)), kv), nk)
            best = min(best, kn(1)*nk*nk + kn(2)*nk + kn(3))
         end do
         if (best == ik) norb = norb + 1
      end do
      i = nk; j = nir; k = norb
      print '(A,I3,A,I5,A,I5,A,I5)', 'mesh', i, ' spglib', j, ' brute force', k, ' expected', expected
      if (nir /= expected .or. norb /= expected) error stop 'irreducible mesh count differs from expected'
      if (abs(sum(wt) - 1.0_rp) > 1.0e-12_rp) error stop 'weights do not sum to 1'
      call spg%finalize()
   end subroutine mesh_count

end program test_magnetic_symmetry
