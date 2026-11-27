!> Displacement operator of calculate_jijk in the code basis. Expected values: the selection rule of a displacement along z,
!> Delta m = 0, with m read from the code's L_z after hcpx('cart2sph') (diagonal 2..-2 within the d shell).
!> Catches: disp_matrix (real cubic order) used as if it were in the spherical basis (the raw operator couples s-z, x-zx, ...).
program test_jijk_disp_basis
   use precision_mod, only: rp
   use basis_mod, only: basis_init
   use math_mod, only: init_math_operators, hcpx, L_z
   use symbolic_atom_mod, only: symbolic_atom
   use exchange_mod, only: spin_blocks_hcpx
   implicit none

   integer, parameter :: n = 9
   type(symbolic_atom) :: atom
   complex(rp), allocatable :: pmat(:, :, :), pcub(:, :, :), u(:, :, :), u_raw(:, :, :)
   complex(rp) :: lz(n, n)
   real(rp) :: viol, viol_raw, big
   integer :: i, j

   call basis_init(2)
   call init_math_operators()
   lz = L_z
   call hcpx(lz, 'cart2sph')
   allocate (pmat(2*n, 2*n, 1), pcub(2*n, 2*n, 1), u(2*n, 2*n, 1), u_raw(2*n, 2*n, 1))
   pmat = 0.0_rp
   do i = 1, 2*n
      pmat(i, i, 1) = 1.0_rp
   end do
   call atom%restore_to_default()

   ! Same call sequence as calculate_jijk.
   pcub = pmat
   call spin_blocks_hcpx(pcub, 'sph2cart')
   call atom%udisp_matrix(u, pcub, [0.0_rp, 0.0_rp, 1.0_rp], 2, 1.0_rp, 1)
   call spin_blocks_hcpx(u, 'cart2sph')
   ! Unconverted operator, as before the fix.
   call atom%udisp_matrix(u_raw, pmat, [0.0_rp, 0.0_rp, 1.0_rp], 2, 1.0_rp, 1)

   viol = 0.0_rp
   viol_raw = 0.0_rp
   big = maxval(abs(u))
   do i = 1, n
      do j = 1, n
         if (abs(real(lz(i, i) - lz(j, j), rp)) > 0.5_rp) then
            viol = max(viol, abs(u(i, j, 1)), abs(u(n + i, n + j, 1)))
            viol_raw = max(viol_raw, abs(u_raw(i, j, 1)), abs(u_raw(n + i, n + j, 1)))
         end if
      end do
   end do
   print '(a,es10.2,a,es10.2,a,es10.2)', 'max |U| ', big, '  Delta m /= 0 elements: converted', viol, ' raw', viol_raw
   if (big < 0.1_rp .or. viol > 1.0e-12_rp .or. viol_raw < 0.1_rp) then
      print '(a)', 'FAIL'
      error stop 1
   end if
   print '(a)', 'test_jijk_disp_basis: all checks passed'
end program test_jijk_disp_basis
