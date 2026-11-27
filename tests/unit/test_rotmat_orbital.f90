!> ROTMAT orbital block. Expected values: R = Rz(alpha) Ry(beta) written below, the code's L operators after hcpx('cart2sph'),
!> and the invariance of L.sigma under a joint rotation of orbital and spin space.
!> Catches: orbital block of ROTMAT equal to conj(D) of the spherical basis, visible only for alpha /= 0.
program test_rotmat_orbital
   use basis_mod, only: basis_init
   use math_mod, only: init_math_operators, hcpx, L_x, L_y, L_z, ROTMAT
   implicit none
   integer, parameter :: rp = kind(1.0d0)
   real(rp), parameter :: tol = 1.0e-12_rp
   real(rp) :: angs(2, 4), a, b, rz(3, 3), ry(3, 3), r(3, 3), worst
   complex(rp) :: lc(9, 9, 3), dm(18, 18), d(9, 9), h(18, 18), t(9, 9), rl(9, 9)
   complex(rp), parameter :: sig(2, 2, 3) = reshape([(0.0_rp, 0.0_rp), (1.0_rp, 0.0_rp), (1.0_rp, 0.0_rp), (0.0_rp, 0.0_rp), &
                                                    (0.0_rp, 0.0_rp), (0.0_rp, 1.0_rp), (0.0_rp, -1.0_rp), (0.0_rp, 0.0_rp), &
                                                    (1.0_rp, 0.0_rp), (0.0_rp, 0.0_rp), (0.0_rp, 0.0_rp), (-1.0_rp, 0.0_rp)], [2, 2, 3])
   integer :: ia, i, j, k, p, q

   angs(:, 1) = [0.0_rp, 1.3_rp]
   angs(:, 2) = [1.5707963267948966_rp, 1.5707963267948966_rp]
   angs(:, 3) = [0.7853981633974483_rp, 0.9553166181245093_rp]
   angs(:, 4) = [1.1071487177940904_rp, 0.6405223126794245_rp]
   call basis_init(2)
   call init_math_operators()
   lc(:, :, 1) = L_x; lc(:, :, 2) = L_y; lc(:, :, 3) = L_z
   do k = 1, 3
      call hcpx(lc(:, :, k), 'cart2sph')
   end do

   worst = 0.0_rp
   do ia = 1, 4
      a = angs(1, ia); b = angs(2, ia)
      rz = reshape([cos(a), sin(a), 0.0_rp, -sin(a), cos(a), 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp], [3, 3])
      ry = reshape([cos(b), 0.0_rp, -sin(b), 0.0_rp, 1.0_rp, 0.0_rp, sin(b), 0.0_rp, cos(b)], [3, 3])
      r = matmul(rz, ry)
      call ROTMAT(dm, a, b, 0.0_rp)
      d = dm(1:9, 1:9)/dm(1, 1)   ! orbital block: dm(1,1) is the s-s element times the spin phase
      do i = 1, 3
         t = matmul(d, matmul(lc(:, :, i), conjg(transpose(d))))
         rl = 0.0_rp
         do j = 1, 3
            rl = rl + r(j, i)*lc(:, :, j)
         end do
         worst = max(worst, maxval(abs(t - rl)))
      end do
      h = 0.0_rp
      do i = 1, 3
         do p = 1, 2
            do q = 1, 2
               h((p - 1)*9 + 1:p*9, (q - 1)*9 + 1:q*9) = h((p - 1)*9 + 1:p*9, (q - 1)*9 + 1:q*9) + sig(p, q, i)*lc(:, :, i)
            end do
         end do
      end do
      worst = max(worst, maxval(abs(matmul(conjg(transpose(dm)), matmul(h, dm)) - h)))
      print '(a,2f8.4,a,es10.2)', 'alpha beta', a, b, ' running max residual', worst
   end do
   if (worst > tol) then
      print '(a,es10.2)', 'FAIL: residual', worst
      error stop 1
   end if
   print '(a)', 'test_rotmat_orbital: all checks passed'
end program test_rotmat_orbital
