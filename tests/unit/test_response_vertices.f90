!------------------------------------------------------------------------------
! TDDFT-02 -- generalized response vertices
!------------------------------------------------------------------------------
program test_response_vertices
   use precision_mod, only: rp
   use math_mod, only: i_unit, pi
   use response_components_mod, only: RESPONSE_CHARGE, RESPONSE_MX, RESPONSE_MY, RESPONSE_MZ, &
      RESPONSE_PLUS, RESPONSE_MINUS
   use response_vertices_mod, only: response_channel, site_projected_operator, &
      response_transition_vertex, response_transition_vectors
   implicit none

   real(rp), parameter :: tol = 1.0e-13_rp
   logical :: failed

   failed = .false.
   call test_dense_reference_and_batch()
   call test_collinear_selection_rules()
   call test_global_spin_rotation_covariance()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   end if
   write (*, '(a)') 'RESULT: PASS'

contains

   subroutine test_dense_reference_and_batch()
      integer, parameter :: counts(2) = [4, 2]
      type(response_channel) :: channels(6)
      complex(rp) :: bra(12), ket(12), dense_value
      complex(rp) :: bra_batch(12, 2), ket_batch(12, 2), vertices(6, 2)
      complex(rp), allocatable :: op(:, :)
      integer :: ichannel

      channels(1) = response_channel(1, RESPONSE_CHARGE)
      channels(2) = response_channel(1, RESPONSE_MX, l_sector=1)
      channels(3) = response_channel(1, RESPONSE_MY, orbital=3)
      channels(4) = response_channel(2, RESPONSE_MZ)
      channels(5) = response_channel(2, RESPONSE_PLUS, orbital=2)
      channels(6) = response_channel(1, RESPONSE_MINUS, l_sector=0)

      bra = [(cmplx(0.07_rp*real(ichannel, rp), -0.03_rp*real(ichannel - 2, rp), rp), ichannel=1, 12)]
      ket = [(cmplx(-0.02_rp*real(ichannel - 4, rp), 0.05_rp*real(ichannel, rp), rp), ichannel=1, 12)]
      do ichannel = 1, size(channels)
         op = site_projected_operator(channels(ichannel), counts)
         dense_value = sum(conjg(bra)*matmul(op, ket))
         call check_complex('dense response operator', response_transition_vertex(channels(ichannel), counts, bra, ket), &
            dense_value)
      end do

      bra_batch(:, 1) = bra
      ket_batch(:, 1) = ket
      bra_batch(:, 2) = cmplx(0.4_rp, -0.2_rp, rp)*ket
      ket_batch(:, 2) = cmplx(-0.3_rp, 0.6_rp, rp)*bra
      call response_transition_vectors(channels, counts, bra_batch, ket_batch, vertices)
      do ichannel = 1, size(channels)
         call check_complex('batched transition vector first pair', vertices(ichannel, 1), &
            response_transition_vertex(channels(ichannel), counts, bra_batch(:, 1), ket_batch(:, 1)))
         call check_complex('batched transition vector second pair', vertices(ichannel, 2), &
            response_transition_vertex(channels(ichannel), counts, bra_batch(:, 2), ket_batch(:, 2)))
      end do

      ! Circular response variables are m+/- = mx +/- i my, including the
      ! TDDFT-00 factor-of-two convention rather than half-normalized ladders.
      call check_complex('plus transform', response_transition_vertex(response_channel(1, RESPONSE_PLUS), counts, bra, ket), &
         response_transition_vertex(response_channel(1, RESPONSE_MX), counts, bra, ket) + &
         i_unit*response_transition_vertex(response_channel(1, RESPONSE_MY), counts, bra, ket))
      call check_complex('minus transform', response_transition_vertex(response_channel(1, RESPONSE_MINUS), counts, bra, ket), &
         response_transition_vertex(response_channel(1, RESPONSE_MX), counts, bra, ket) - &
         i_unit*response_transition_vertex(response_channel(1, RESPONSE_MY), counts, bra, ket))
   end subroutine test_dense_reference_and_batch

   subroutine test_collinear_selection_rules()
      integer, parameter :: counts(1) = [2]
      type(response_channel) :: charge, mz, plus, minus
      complex(rp) :: bra(4), ket(4), vc, vz, vp, vm
      integer :: n, m

      charge = response_channel(1, RESPONSE_CHARGE)
      mz = response_channel(1, RESPONSE_MZ)
      plus = response_channel(1, RESPONSE_PLUS)
      minus = response_channel(1, RESPONSE_MINUS)
      do n = 1, 4
         bra = cmplx(0.0_rp, 0.0_rp, rp)
         bra(n) = cmplx(1.0_rp, 0.0_rp, rp)
         do m = 1, 4
            ket = cmplx(0.0_rp, 0.0_rp, rp)
            ket(m) = cmplx(1.0_rp, 0.0_rp, rp)
            vc = response_transition_vertex(charge, counts, bra, ket)
            vz = response_transition_vertex(mz, counts, bra, ket)
            vp = response_transition_vertex(plus, counts, bra, ket)
            vm = response_transition_vertex(minus, counts, bra, ket)
            if ((n <= 2 .and. m > 2) .or. (n > 2 .and. m <= 2)) then
               call check_real('charge is spin conserving', abs(vc), 0.0_rp)
               call check_real('mz is spin conserving', abs(vz), 0.0_rp)
            else
               call check_real('plus/minus are spin flipping', abs(vp) + abs(vm), 0.0_rp)
            end if
            call check_real('transverse and charge-longitudinal decouple', &
               abs(vc*conjg(vp)) + abs(vc*conjg(vm)) + abs(vz*conjg(vp)) + abs(vz*conjg(vm)), 0.0_rp)
         end do
      end do
   end subroutine test_collinear_selection_rules

   subroutine test_global_spin_rotation_covariance()
      integer, parameter :: counts(2) = [2, 1]
      type(response_channel) :: charge, mx, my, mz
      complex(rp) :: bra(6), ket(6), bra_rotated(6), ket_rotated(6), u(2, 2)
      complex(rp) :: vc, vx, vy, vz, vc_rotated, vx_rotated, vy_rotated, vz_rotated
      real(rp) :: theta
      integer :: i

      charge = response_channel(1, RESPONSE_CHARGE)
      mx = response_channel(1, RESPONSE_MX)
      my = response_channel(1, RESPONSE_MY)
      mz = response_channel(1, RESPONSE_MZ)
      bra = [(cmplx(0.11_rp*real(i, rp), 0.04_rp*real(3 - i, rp), rp), i=1, 6)]
      ket = [(cmplx(-0.06_rp*real(i, rp), 0.09_rp*real(i - 1, rp), rp), i=1, 6)]
      theta = 0.37_rp*pi
      u = reshape([cmplx(cos(0.5_rp*theta), 0.0_rp, rp), cmplx(sin(0.5_rp*theta), 0.0_rp, rp), &
                   cmplx(-sin(0.5_rp*theta), 0.0_rp, rp), cmplx(cos(0.5_rp*theta), 0.0_rp, rp)], [2, 2])
      call rotate_global_spinor(bra, counts, u, bra_rotated)
      call rotate_global_spinor(ket, counts, u, ket_rotated)

      vc = response_transition_vertex(charge, counts, bra, ket)
      vx = response_transition_vertex(mx, counts, bra, ket)
      vy = response_transition_vertex(my, counts, bra, ket)
      vz = response_transition_vertex(mz, counts, bra, ket)
      vc_rotated = response_transition_vertex(charge, counts, bra_rotated, ket_rotated)
      vx_rotated = response_transition_vertex(mx, counts, bra_rotated, ket_rotated)
      vy_rotated = response_transition_vertex(my, counts, bra_rotated, ket_rotated)
      vz_rotated = response_transition_vertex(mz, counts, bra_rotated, ket_rotated)

      call check_complex('charge rotation invariant', vc_rotated, vc)
      call check_complex('x rotation covariance', vx_rotated, cos(theta)*vx + sin(theta)*vz)
      call check_complex('y rotation covariance', vy_rotated, vy)
      call check_complex('z rotation covariance', vz_rotated, -sin(theta)*vx + cos(theta)*vz)
   end subroutine test_global_spin_rotation_covariance

   subroutine rotate_global_spinor(spinor, site_counts, u, rotated)
      complex(rp), intent(in) :: spinor(:), u(2, 2)
      integer, intent(in) :: site_counts(:)
      complex(rp), intent(out) :: rotated(:)
      integer :: site, orbital, offset, nlocal

      offset = 0
      do site = 1, size(site_counts)
         nlocal = site_counts(site)
         do orbital = 1, nlocal
            rotated(offset + orbital) = u(1, 1)*spinor(offset + orbital) + u(1, 2)*spinor(offset + nlocal + orbital)
            rotated(offset + nlocal + orbital) = u(2, 1)*spinor(offset + orbital) + u(2, 2)*spinor(offset + nlocal + orbital)
         end do
         offset = offset + 2*nlocal
      end do
   end subroutine rotate_global_spinor

   subroutine check_complex(label, actual, expected)
      character(len=*), intent(in) :: label
      complex(rp), intent(in) :: actual, expected
      if (abs(actual - expected) > tol) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, abs(actual - expected)
         failed = .true.
      end if
   end subroutine check_complex

   subroutine check_real(label, actual, expected)
      character(len=*), intent(in) :: label
      real(rp), intent(in) :: actual, expected
      if (abs(actual - expected) > tol) then
         write (*, '(a,1x,a,1x,es12.4)') 'FAIL', label, abs(actual - expected)
         failed = .true.
      end if
   end subroutine check_real

end program test_response_vertices
