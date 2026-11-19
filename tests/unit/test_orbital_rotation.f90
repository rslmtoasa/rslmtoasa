!> Orbital and spin rotation matrices (orbital_rotation, spin_rotation) and the site-matrix group average (group_average_dm).
!> Expected values: real harmonics written out as polynomials below (code order s; x,y,z; xy,yz,zx,x2-y2,3z2-r2), Pauli algebra,
!> and closed-form averages. The f block is not tested: init_math_operators gives its L matrices in the spherical basis.
program test_orbital_rotation
   use basis_mod, only: basis_init, norb
   use math_mod, only: inverse_3x3, init_math_operators, orbital_rotation, spin_rotation, group_average_dm
   use spglib_interface_mod, only: spglib_interface
   implicit none
   integer, parameter :: rp = kind(1.0d0)
   real(rp), parameter :: bcc(3, 3) = reshape([-0.5_rp, 0.5_rp, 0.5_rp, 0.5_rp, -0.5_rp, 0.5_rp, 0.5_rp, 0.5_rp, -0.5_rp], [3, 3])
   real(rp), parameter :: fcc(3, 3) = reshape([0.0_rp, 0.5_rp, 0.5_rp, 0.5_rp, 0.0_rp, 0.5_rp, 0.5_rp, 0.5_rp, 0.0_rp], [3, 3])
   type(spglib_interface) :: spg
   integer, allocatable :: rot(:, :, :), flips(:)
   real(rp), allocatable :: tr(:, :), rc(:, :, :)
   real(rp) :: z0(3, 1) = 0.0_rp
   integer :: n

   call basis_init(2)
   call init_math_operators()

   call spg%initialize(bcc, z0, [26])
   n = spg%get_magnetic_symmetry(z0, [26], z0, .false., rot, tr, flips)
   rc = cartesian(bcc, rot)
   print '(A,I3,A)', 'bcc Fe:', n, ' operations'
   call covariance(rc)
   call representation(rc)
   call spin(rc)
   call spg%finalize()
   call site_swap_diamond()
   call site_swap_hand()
   call time_reversal()
   print '(A)', 'test_orbital_rotation: all checks passed'

contains

   !> Cartesian matrices A R A^-1 of fractional rotations (lattice vectors are the columns of a).
   function cartesian(a, rf) result(r)
      real(rp), intent(in) :: a(3, 3)
      integer, intent(in) :: rf(:, :, :)
      real(rp) :: r(3, 3, size(rf, 3)), ai(3, 3)
      integer :: k

      ai = inverse_3x3(a)
      do k = 1, size(rf, 3)
         r(:, :, k) = matmul(a, matmul(real(rf(:, :, k), rp), ai))
         if (maxval(abs(matmul(transpose(r(:, :, k)), r(:, :, k)) - reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3]))) > 1.0e-12_rp) &
            error stop 'cartesian operation is not orthogonal'
      end do
   end function cartesian

   !> Real harmonics s; x,y,z; xy,yz,zx,x2-y2,3z2-r2 of a unit vector, orthonormal relative weights.
   function harmonics(v) result(y)
      real(rp), intent(in) :: v(3)
      real(rp) :: y(9), s3
      s3 = sqrt(3.0_rp)
      y = [1.0_rp, v(1), v(2), v(3), s3*v(1)*v(2), s3*v(2)*v(3), s3*v(3)*v(1), 0.5_rp*s3*(v(1)**2 - v(2)**2), &
           0.5_rp*(3.0_rp*v(3)**2 - 1.0_rp)]
   end function harmonics

   !> Bug caught: wrong rotation sense (R^-1), block ordering or sign of d; missing (-1)^l for improper operations.
   subroutine covariance(r)
      real(rp), intent(in) :: r(:, :, :)
      complex(rp) :: d(norb, norb)
      real(rp) :: v(3), err, e
      integer :: k, i

      err = 0.0_rp
      do k = 1, size(r, 3)
         call orbital_rotation(r(:, :, k), d)
         do i = 1, 6
            v = [sin(0.7_rp*i + 0.3_rp)*cos(1.9_rp*i), sin(0.7_rp*i + 0.3_rp)*sin(1.9_rp*i), cos(0.7_rp*i + 0.3_rp)]
            e = maxval(abs(matmul(d, cmplx(harmonics(v), 0.0_rp, rp)) - harmonics(matmul(r(:, :, k), v))))
            err = max(err, e)
         end do
      end do
      print '(A,ES10.2)', 'covariance D Y(r) = Y(R r), max error', err
      if (err > 1.0e-12_rp) error stop 'D Y(r) /= Y(R r)'
   end subroutine covariance

   !> Supporting only: D(R1 R2) = D(R1) D(R2) (a representation, but not tied to the code's harmonics).
   subroutine representation(r)
      real(rp), intent(in) :: r(:, :, :)
      complex(rp) :: d1(norb, norb), d2(norb, norb), d12(norb, norb)
      real(rp) :: err
      integer :: k

      err = 0.0_rp
      do k = 1, size(r, 3) - 1
         call orbital_rotation(r(:, :, k), d1)
         call orbital_rotation(r(:, :, k + 1), d2)
         call orbital_rotation(matmul(r(:, :, k), r(:, :, k + 1)), d12)
         err = max(err, maxval(abs(matmul(d1, d2) - d12)))
      end do
      print '(A,ES10.2)', 'representation D(R1 R2) = D(R1) D(R2), max error', err
      if (err > 1.0e-12_rp) error stop 'D is not a representation'
   end subroutine representation

   !> Bug caught: SU(2) matrix rotating spin vectors by R^-1 or without det(R): S sigma.m S^dagger = sigma.(Rp m), Rp = det(R) R.
   subroutine spin(r)
      real(rp), intent(in) :: r(:, :, :)
      complex(rp) :: s(2, 2), sm(2, 2), sr(2, 2)
      complex(rp), parameter :: sx(2, 2) = reshape([(0, 0), (1, 0), (1, 0), (0, 0)], [2, 2]), &
                                sy(2, 2) = reshape([(0, 0), (0, 1), (0, -1), (0, 0)], [2, 2]), &
                                sz(2, 2) = reshape([(1, 0), (0, 0), (0, 0), (-1, 0)], [2, 2])
      real(rp) :: m(3), mr(3), err, det
      integer :: k

      m = [0.3_rp, -0.5_rp, 0.8_rp]
      err = 0.0_rp
      do k = 1, size(r, 3)
         det = r(1, 1, k)*(r(2, 2, k)*r(3, 3, k) - r(2, 3, k)*r(3, 2, k)) - r(1, 2, k)*(r(2, 1, k)*r(3, 3, k) - r(2, 3, k)*r(3, 1, k)) &
               + r(1, 3, k)*(r(2, 1, k)*r(3, 2, k) - r(2, 2, k)*r(3, 1, k))
         call spin_rotation(r(:, :, k), s)
         sm = m(1)*sx + m(2)*sy + m(3)*sz
         mr = det*matmul(r(:, :, k), m)
         sr = mr(1)*sx + mr(2)*sy + mr(3)*sz
         err = max(err, maxval(abs(matmul(s, matmul(sm, conjg(transpose(s)))) - sr)))
      end do
      print '(A,ES10.2)', 'spin S sigma.m S^+ = sigma.(det R R m), max error', err
      if (err > 1.0e-12_rp) error stop 'spin rotation differs from the Pauli-algebra result'
   end subroutine spin

   !> Diamond Si, scalar site densities q_i: the average is (q1 + q2)/2 on both sites for the 48 operations.
   !> Bug caught: site permutation index ignored (identity) or applied to the wrong side.
   subroutine site_swap_diamond()
      real(rp), parameter :: pos(3, 2) = reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.25_rp, 0.25_rp, 0.25_rp], [3, 2])
      real(rp), allocatable :: rcd(:, :, :)
      complex(rp), allocatable :: u(:, :, :), rho(:, :, :)
      integer, allocatable :: perm(:, :)
      real(rp) :: x(3), dd(3)
      integer :: nops, k, i, j, a

      call spg%initialize(fcc, pos, [14, 14])
      nops = spg%get_magnetic_symmetry(pos, [14, 14], reshape([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], [3, 2]), &
                                       .false., rot, tr, flips)
      rcd = cartesian(fcc, rot)
      allocate(perm(2, nops), u(2*norb, 2*norb, nops), rho(2*norb, 2*norb, 2))
      u = (0.0_rp, 0.0_rp)
      do k = 1, nops
         do i = 1, 2
            x = matmul(real(rot(:, :, k), rp), pos(:, i)) + tr(:, k)
            perm(i, k) = 0
            do j = 1, 2
               dd = x - pos(:, j)
               if (maxval(abs(dd - nint(dd))) < 1.0e-6_rp) perm(i, k) = j
            end do
         end do
         call orbital_rotation(rcd(:, :, k), u(1:norb, 1:norb, k))
         u(norb + 1:, norb + 1:, k) = u(1:norb, 1:norb, k)
      end do
      rho = (0.0_rp, 0.0_rp)
      do a = 1, 2*norb
         rho(a, a, 1) = 1.0_rp
         rho(a, a, 2) = 3.0_rp
      end do
      call group_average_dm(rho, perm, u, flips < 0)
      print '(A,I3,A,I3,A,2ES10.2)', 'diamond Si:', nops, ' operations,', count(perm(1, :) == 2), ' swap the sites, max error', &
         maxval(abs(rho(:, :, 1) - eye(2*norb)*2.0_rp)), maxval(abs(rho(:, :, 2) - eye(2*norb)*2.0_rp))
      if (maxval(abs(rho(:, :, 1) - eye(2*norb)*2.0_rp)) > 1.0e-12_rp .or. maxval(abs(rho(:, :, 2) - eye(2*norb)*2.0_rp)) > 1.0e-12_rp) &
         error stop 'swap average of scalar site densities differs from (q1 + q2)/2'
      call spg%finalize()
   end subroutine site_swap_diamond

   function eye(n) result(e)
      integer, intent(in) :: n
      complex(rp) :: e(n, n)
      integer :: i
      e = (0.0_rp, 0.0_rp)
      do i = 1, n
         e(i, i) = 1.0_rp
      end do
   end function eye

   !> Two sites, {E, swap with U = cyclic permutation 1 -> 2 -> 3 -> 1 of the first three components (U^3 = 1, U /= U^dagger)}, unequal non-scalar input. U has one entry s_i in row i,
   !> column c_i, so (U B U^dagger)_ij = s_i s_j B_(c_i c_j) and the average is closed-form elementwise.
   !> Bug caught: permutation index, U versus U^dagger placement, matrix not conjugated by U on the correct side.
   subroutine site_swap_hand()
      complex(rp) :: rho(4, 4, 2), a(4, 4), b(4, 4), u(4, 4, 2), want1(4, 4), want2(4, 4)
      real(rp) :: sg(4) = [1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp]
      integer :: c(4) = [2, 3, 1, 4]
      integer :: i, j

      do i = 1, 4
         do j = 1, 4
            a(i, j) = cmplx(1.0_rp*i + 0.5_rp*j, 0.25_rp*(i - j), rp)
            b(i, j) = cmplx(2.0_rp*j - 0.3_rp*i, -0.4_rp*(i - j), rp)
         end do
      end do
      u = (0.0_rp, 0.0_rp)
      do i = 1, 4
         u(i, i, 1) = 1.0_rp
         u(i, c(i), 2) = sg(i)
      end do
      rho(:, :, 1) = a
      rho(:, :, 2) = b
      do i = 1, 4
         do j = 1, 4
            want1(i, j) = 0.5_rp*(a(i, j) + sg(i)*sg(j)*b(c(i), c(j)))
            want2(i, j) = 0.5_rp*(b(i, j) + sg(i)*sg(j)*a(c(i), c(j)))
         end do
      end do
      call group_average_dm(rho, reshape([1, 2, 2, 1], [2, 2]), u, [.false., .false.])
      print '(A,ES10.2)', 'two-site swap, closed form, max error', max(maxval(abs(rho(:, :, 1) - want1)), maxval(abs(rho(:, :, 2) - want2)))
      if (maxval(abs(rho(:, :, 1) - want1)) > 1.0e-13_rp .or. maxval(abs(rho(:, :, 2) - want2)) > 1.0e-13_rp) &
         error stop 'two-site swap average differs from the closed form'
   end subroutine site_swap_hand

   !> One site, {E, T} with T = i sigma_y K, 2 orbitals per spin: rho_uu -> (rho_uu + conj(rho_dd))/2, rho_ud -> (rho_ud - conj(rho_du))/2.
   !> Bug caught: missing conjugation, or missing sign of the spin-flip blocks.
   subroutine time_reversal()
      complex(rp) :: rho(4, 4, 1), want(4, 4), u(4, 4, 2)
      integer :: i, j

      do i = 1, 4
         do j = 1, 4
            rho(i, j, 1) = cmplx(1.0_rp*i + 0.5_rp*j, 0.25_rp*(i - 2*j), rp)
         end do
      end do
      u = (0.0_rp, 0.0_rp)
      u(:, :, 1) = eye(4)
      u(:, :, 2) = eye(4)
      do i = 1, 2
         do j = 1, 2
            want(i, j) = 0.5_rp*(rho(i, j, 1) + conjg(rho(i + 2, j + 2, 1)))
            want(i + 2, j + 2) = 0.5_rp*(rho(i + 2, j + 2, 1) + conjg(rho(i, j, 1)))
            want(i, j + 2) = 0.5_rp*(rho(i, j + 2, 1) - conjg(rho(i + 2, j, 1)))
            want(i + 2, j) = 0.5_rp*(rho(i + 2, j, 1) - conjg(rho(i, j + 2, 1)))
         end do
      end do
      call group_average_dm(rho, reshape([1, 1], [1, 2]), u, [.false., .true.])
      print '(A,ES10.2)', 'time reversal, closed form, max error', maxval(abs(rho(:, :, 1) - want))
      if (maxval(abs(rho(:, :, 1) - want)) > 1.0e-13_rp) error stop 'time-reversal average differs from the closed form'
   end subroutine time_reversal

end program test_orbital_rotation
