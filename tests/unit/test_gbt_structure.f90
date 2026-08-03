program test_gbt_structure
   use precision_mod, only: rp
   use gbt_structure_mod, only: gbt_bond_phase, gbt_endpoint_link, gbt_lift_orbital_block, gbt_contract_collinear
   implicit none

   real(rp), parameter :: tol = 1.0e-12_rp
   real(rp) :: max_error
   logical :: failed

   failed = .false.
   max_error = 0.0_rp
   call test_phase_units()
   call test_dense_link_and_reverse()
   call test_lift_and_contract()
   call test_q0_and_non_noop()
   call test_stoner_shift()

   write (*, '(a,es12.4)') 'GBT S-level maximum error: ', max_error
   if (failed .or. max_error >= tol) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine update_error(value, label)
      real(rp), intent(in) :: value
      character(len=*), intent(in) :: label
      max_error = max(max_error, value)
      if (value >= tol) then
         write (*, '(a,a,a,es12.4)') 'FAIL ', trim(label), ': ', value
         failed = .true.
      end if
   end subroutine update_error

   subroutine test_phase_units()
      real(rp) :: q(3), d(3), alpha, expected
      q = [0.17_rp, -0.08_rp, 0.21_rp]
      d = [2.1_rp, -0.7_rp, 1.4_rp]
      expected = 2.0_rp*acos(-1.0_rp)*(q(1)*d(1) + q(2)*d(2) + q(3)*d(3))/3.7_rp
      alpha = gbt_bond_phase(q, d, 3.7_rp)
      call update_error(abs(alpha - expected), 'phase units')
   end subroutine test_phase_units

   subroutine test_dense_link_and_reverse()
      real(rp) :: q(3), d(3), ta, pa, tb, pb, alpha
      complex(rp) :: link(2, 2), reverse(2, 2), ua(2, 2), ub(2, 2), rz(2, 2), oracle(2, 2)

      q = [0.13_rp, -0.04_rp, 0.09_rp]
      d = [1.2_rp, -2.3_rp, 0.8_rp]
      ta = 0.73_rp; pa = -0.41_rp
      tb = 1.12_rp; pb = 0.66_rp
      alpha = 2.0_rp*acos(-1.0_rp)*(q(1)*d(1) + q(2)*d(2) + q(3)*d(3))/2.9_rp
      call oracle_frame(ta, pa, ua)
      call oracle_frame(tb, pb, ub)
      rz = cmplx(0.0_rp, 0.0_rp, rp)
      rz(1, 1) = exp(cmplx(0.0_rp, -0.5_rp*alpha, rp))
      rz(2, 2) = exp(cmplx(0.0_rp,  0.5_rp*alpha, rp))
      oracle = matmul(transpose(conjg(ua)), matmul(rz, ub))

      call gbt_endpoint_link(q, d, 2.9_rp, ta, pa, tb, pb, link)
      call update_error(maxval(abs(link - oracle)), 'dense endpoint link')
      call gbt_endpoint_link(q, -d, 2.9_rp, tb, pb, ta, pa, reverse)
      call update_error(maxval(abs(reverse - transpose(conjg(link)))), 'reverse directed link')
   end subroutine test_dense_link_and_reverse

   subroutine test_lift_and_contract()
      complex(rp) :: s(2, 2), saved(2, 2), link(2, 2), lifted(4, 4), h(4, 4), oracle(4, 4)
      complex(rp) :: w0a(2), w1a(2), w0b(2), w1b(2)
      complex(rp) :: wa(4), wb(4)
      integer :: i, j

      s(1, 1) = cmplx(0.7_rp, -0.2_rp, rp); s(1, 2) = cmplx(-0.1_rp, 0.4_rp, rp)
      s(2, 1) = cmplx(0.3_rp, 0.1_rp, rp); s(2, 2) = cmplx(-0.5_rp, -0.6_rp, rp)
      saved = s
      call gbt_endpoint_link([0.11_rp, 0.0_rp, 0.0_rp], [1.7_rp, 0.0_rp, 0.0_rp], 3.1_rp, &
                             0.4_rp, 0.2_rp, 1.0_rp, -0.3_rp, link)
      call gbt_lift_orbital_block(s, link, lifted)
      call update_error(maxval(abs(s - saved)), 'raw S preservation')

      w0a = [cmplx(0.8_rp, 0.1_rp, rp), cmplx(1.1_rp, -0.2_rp, rp)]
      w1a = [cmplx(0.2_rp, -0.1_rp, rp), cmplx(-0.3_rp, 0.05_rp, rp)]
      w0b = [cmplx(0.9_rp, -0.15_rp, rp), cmplx(0.6_rp, 0.25_rp, rp)]
      w1b = [cmplx(-0.1_rp, 0.08_rp, rp), cmplx(0.4_rp, -0.12_rp, rp)]
      wa = [w0a + w1a, w0a - w1a]
      wb = [w0b - w1b, w0b + w1b]
      do i = 1, 4
         do j = 1, 4
            oracle(i, j) = wa(i)*lifted(i, j)*wb(j)
         end do
      end do
      call gbt_contract_collinear(lifted, w0a, w1a, w0b, w1b, 1.0_rp, -1.0_rp, h)
      call update_error(maxval(abs(h - oracle)), 'unequal endpoint collinear contraction')
   end subroutine test_lift_and_contract

   subroutine test_q0_and_non_noop()
      complex(rp) :: identity_link(2, 2), finite_link(2, 2), relative_link(2, 2), eye(2, 2)
      eye = cmplx(0.0_rp, 0.0_rp, rp)
      eye(1, 1) = cmplx(1.0_rp, 0.0_rp, rp); eye(2, 2) = eye(1, 1)
      call gbt_endpoint_link([0.0_rp, 0.0_rp, 0.0_rp], [1.0_rp, 0.0_rp, 0.0_rp], 2.0_rp, &
                             0.63_rp, -0.2_rp, 0.63_rp, -0.2_rp, identity_link)
      call update_error(maxval(abs(identity_link - eye)), 'q=0 common-frame identity')
      call gbt_endpoint_link([0.19_rp, 0.0_rp, 0.0_rp], [1.0_rp, 0.0_rp, 0.0_rp], 2.0_rp, &
                             0.63_rp, -0.2_rp, 0.63_rp, -0.2_rp, finite_link)
      if (maxval(abs(finite_link - eye)) < 1.0e-3_rp) then
         write (*, '(a)') 'FAIL finite-q link is a no-op'
         failed = .true.
      end if
      call gbt_endpoint_link([0.0_rp, 0.0_rp, 0.0_rp], [1.0_rp, 0.0_rp, 0.0_rp], 2.0_rp, &
                             0.2_rp, 0.1_rp, 1.1_rp, -0.4_rp, relative_link)
      if (maxval(abs(relative_link - eye)) < 1.0e-3_rp) then
         write (*, '(a)') 'FAIL relative-sublattice q=0 link is a no-op'
         failed = .true.
      end if
   end subroutine test_q0_and_non_noop

   subroutine test_stoner_shift()
      real(rp) :: q(3), d(3), alat, kval, alpha
      complex(rp) :: gp(2, 2), gm(2, 2), hp(2, 2), hm(2, 2), hk(2, 2)
      complex(rp) :: s1(1, 1), lifted(2, 2), w0(1), w1(1), phase_k
      real(rp) :: expected_up, expected_dn

      q = [0.16_rp, 0.0_rp, 0.0_rp]; d = [1.3_rp, 0.0_rp, 0.0_rp]; alat = 2.7_rp
      kval = 0.47_rp
      s1(1, 1) = cmplx(0.8_rp, 0.0_rp, rp)
      w0(1) = cmplx(1.05_rp, 0.0_rp, rp); w1(1) = cmplx(0.21_rp, 0.0_rp, rp)
      call gbt_endpoint_link(q, d, alat, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, gp)
      call gbt_lift_orbital_block(s1, gp, lifted)
      call gbt_contract_collinear(lifted, w0, w1, w0, w1, 1.0_rp, 1.0_rp, hp)
      call gbt_endpoint_link(q, -d, alat, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, gm)
      call gbt_lift_orbital_block(s1, gm, lifted)
      call gbt_contract_collinear(lifted, w0, w1, w0, w1, 1.0_rp, 1.0_rp, hm)
      phase_k = exp(cmplx(0.0_rp, kval, rp))
      hk = hp*phase_k + hm*conjg(phase_k)
      alpha = gbt_bond_phase(q, d, alat)
      expected_up = 2.0_rp*real(s1(1, 1), rp)*real((w0(1) + w1(1))**2, rp)*cos(kval - 0.5_rp*alpha)
      expected_dn = 2.0_rp*real(s1(1, 1), rp)*real((w0(1) - w1(1))**2, rp)*cos(kval + 0.5_rp*alpha)
      call update_error(abs(real(hk(1, 1), rp) - expected_up), 'Stoner k-q/2 up shift')
      call update_error(abs(real(hk(2, 2), rp) - expected_dn), 'Stoner k+q/2 down shift')
      call update_error(abs(hk(1, 2)) + abs(hk(2, 1)), 'Stoner shifted spin diagonality')
   end subroutine test_stoner_shift

   subroutine oracle_frame(theta, phi, frame)
      real(rp), intent(in) :: theta, phi
      complex(rp), intent(out) :: frame(2, 2)
      complex(rp) :: ep, em
      ep = exp(cmplx(0.0_rp, phi, rp)); em = conjg(ep)
      frame(1, 1) = cmplx(cos(theta/2.0_rp), 0.0_rp, rp)
      frame(1, 2) = -em*sin(theta/2.0_rp)
      frame(2, 1) = ep*sin(theta/2.0_rp)
      frame(2, 2) = frame(1, 1)
   end subroutine oracle_frame

end program test_gbt_structure
