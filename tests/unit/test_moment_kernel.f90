!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_moment_kernel
!
!> @brief Exact double-Chebyshev moment known-answer tests (B5.1), dependency-free.
!> @details Pins the eigenbasis moment identity `moment_onsite_block` relies on,
!>          T_p(H~) = U diag(T_p(eps~)) U^dagger, with NO LMTO machinery: the
!>          k-space Hamiltonian is CONSTRUCTED from known eigenpairs H = U D U^dag
!>          so the kernel (which consumes eps, U) can be checked against a direct
!>          MATRIX Chebyshev recurrence on H (which never touches the eigenbasis).
!>          If the two agree, the identity -- the whole basis of the exact moment
!>          generator -- holds.
!>
!>          Pins:
!>          1. 1-band chain: mu_nm reduces to the scalar BZ sum
!>             (1/Nk) sum_k T_m(eps~) va T_n(eps~) vb -- kernel vs closed sum.
!>          2. 2x2 eigenbasis identity: random unitary U, random Hermitian
!>             velocities; kernel(eps,U) vs direct matrix-Chebyshev block, to
!>             machine precision (this is the real known answer).
!>          3. Multi-site site-block extraction: 4x4 (two 2-orbital sites); the
!>             kernel's ioff site block equals the direct operator's site block.
!>
!>          Exits non-zero (error stop) on any failure so ctest registers a fail.
!------------------------------------------------------------------------------
program test_moment_kernel
   use precision_mod, only: rp
   use moment_kernel_mod, only: moment_onsite_block
   implicit none

   real(rp), parameter :: tol = 1.0e-11_rp
   complex(rp), parameter :: iu = (0.0_rp, 1.0_rp)
   logical :: failed

   failed = .false.

   call test_chain_scalar()
   call test_eigenbasis_identity()
   call test_site_block()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !> Build T_{p-1}(M) for p=1..ll as explicit matrices (Chebyshev recurrence).
   subroutine matrix_chebyshev(m, a, b, ll, tmat)
      complex(rp), intent(in) :: m(:, :)
      real(rp), intent(in) :: a, b
      integer, intent(in) :: ll
      complex(rp), intent(out) :: tmat(:, :, :)   ! (n, n, ll)
      integer :: n, p, i
      complex(rp), allocatable :: mt(:, :)

      n = size(m, 1)
      allocate (mt(n, n))
      mt = m/cmplx(a, 0.0_rp, rp)
      do i = 1, n
         mt(i, i) = mt(i, i) - cmplx(b/a, 0.0_rp, rp)
      end do
      tmat(:, :, 1) = (0.0_rp, 0.0_rp)
      do i = 1, n
         tmat(i, i, 1) = (1.0_rp, 0.0_rp)
      end do
      if (ll >= 2) tmat(:, :, 2) = mt
      do p = 3, ll
         tmat(:, :, p) = 2.0_rp*matmul(mt, tmat(:, :, p - 1)) - tmat(:, :, p - 2)
      end do
      deallocate (mt)
   end subroutine matrix_chebyshev

   !> 2x2 unitary from (theta, phi): columns orthonormal by construction.
   function unitary2(theta, phi) result(u)
      real(rp), intent(in) :: theta, phi
      complex(rp) :: u(2, 2)
      u(1, 1) = cmplx(cos(theta), 0.0_rp, rp)
      u(2, 1) = sin(theta)*exp(iu*phi)
      u(1, 2) = -sin(theta)*exp(-iu*phi)
      u(2, 2) = cmplx(cos(theta), 0.0_rp, rp)
   end function unitary2

   !> Test 1: 1-band chain, mu_nm vs the closed scalar BZ sum.
   subroutine test_chain_scalar()
      integer, parameter :: nk = 7, ll = 5
      real(rp), parameter :: a = 3.0_rp, b = 0.2_rp, t = 1.0_rp
      real(rp) :: ev(1, nk), et, kf
      complex(rp) :: u(1, 1, nk), va(1, 1, nk), vb(1, 1, nk)
      complex(rp) :: mu(1, 1, ll, ll), ref
      real(rp) :: tcheb(ll), maxerr
      integer :: ik, n, m, p

      do ik = 1, nk
         kf = real(ik - 1, rp)/real(nk, rp)
         ev(1, ik) = -2.0_rp*t*cos(6.283185307179586_rp*kf)
         u(1, 1, ik) = (1.0_rp, 0.0_rp)
         va(1, 1, ik) = cmplx(2.0_rp*t*sin(6.283185307179586_rp*kf), 0.0_rp, rp)
         vb(1, 1, ik) = cmplx(0.5_rp + kf, 0.0_rp, rp)   ! arbitrary real Hermitian scalar
      end do

      call moment_onsite_block(ev, u, va, vb, a, b, 0, 1, ll, mu)

      maxerr = 0.0_rp
      do n = 1, ll
         do m = 1, ll
            ref = (0.0_rp, 0.0_rp)
            do ik = 1, nk
               et = (ev(1, ik) - b)/a
               tcheb(1) = 1.0_rp
               if (ll >= 2) tcheb(2) = et
               do p = 3, ll
                  tcheb(p) = 2.0_rp*et*tcheb(p - 1) - tcheb(p - 2)
               end do
               ref = ref + tcheb(m)*va(1, 1, ik)*tcheb(n)*vb(1, 1, ik)
            end do
            ref = ref/real(nk, rp)
            maxerr = max(maxerr, abs(mu(1, 1, n, m) - ref))
         end do
      end do
      write (*, '(a,es12.4)') 'Test 1 (1-band chain scalar)        max_err = ', maxerr
      if (maxerr > tol) failed = .true.
   end subroutine test_chain_scalar

   !> Test 2: 2x2 eigenbasis identity vs direct matrix Chebyshev.
   subroutine test_eigenbasis_identity()
      integer, parameter :: nk = 5, n = 2, ll = 6
      real(rp), parameter :: a = 2.5_rp, b = -0.3_rp
      real(rp) :: ev(n, nk)
      complex(rp) :: u(n, n, nk), va(n, n, nk), vb(n, n, nk)
      complex(rp) :: mu(n, n, ll, ll)
      complex(rp) :: hk(n, n), tmat(n, n, ll), oref(n, n)
      real(rp) :: maxerr, th, ph
      integer :: ik, nn, m, i, j

      do ik = 1, nk
         th = 0.3_rp + 0.4_rp*real(ik, rp)
         ph = 0.7_rp*real(ik, rp)
         ev(1, ik) = -1.0_rp + 0.3_rp*real(ik, rp)
         ev(2, ik) = 0.8_rp - 0.2_rp*real(ik, rp)
         u(:, :, ik) = unitary2(th, ph)
         ! Random Hermitian velocity operators (distinct a and b).
         va(:, :, ik) = herm2(0.5_rp*real(ik, rp), -0.2_rp, 0.9_rp, 0.4_rp)
         vb(:, :, ik) = herm2(-0.3_rp, 0.6_rp*real(ik, rp), 0.1_rp, -0.5_rp)
      end do

      call moment_onsite_block(ev, u, va, vb, a, b, 0, n, ll, mu)

      maxerr = 0.0_rp
      do nn = 1, ll
         do m = 1, ll
            oref = (0.0_rp, 0.0_rp)
            do ik = 1, nk
               ! H(k) = U diag(ev) U^dagger, then direct matrix Chebyshev.
               hk = hbuild(u(:, :, ik), ev(:, ik))
               call matrix_chebyshev(hk, a, b, ll, tmat)
               oref = oref + matmul(tmat(:, :, m), matmul(va(:, :, ik), &
                              matmul(tmat(:, :, nn), vb(:, :, ik))))
            end do
            oref = oref/real(nk, rp)
            do i = 1, n
               do j = 1, n
                  maxerr = max(maxerr, abs(mu(i, j, nn, m) - oref(i, j)))
               end do
            end do
         end do
      end do
      write (*, '(a,es12.4)') 'Test 2 (2x2 eigenbasis identity)    max_err = ', maxerr
      if (maxerr > tol) failed = .true.
   end subroutine test_eigenbasis_identity

   !> Test 3: 4x4 (two 2-orbital sites); ioff site block vs direct block.
   subroutine test_site_block()
      integer, parameter :: nk = 4, n = 4, nb = 2, ll = 5
      real(rp), parameter :: a = 3.2_rp, b = 0.1_rp
      integer, parameter :: ioff = 2   ! second site block, rows/cols 3..4
      real(rp) :: ev(n, nk)
      complex(rp) :: u(n, n, nk), va(n, n, nk), vb(n, n, nk)
      complex(rp) :: mu(nb, nb, ll, ll)
      complex(rp) :: hk(n, n), tmat(n, n, ll), oref(n, n)
      real(rp) :: maxerr
      integer :: ik, nn, m, i, j

      do ik = 1, nk
         ev(:, ik) = [ -1.2_rp+0.1_rp*ik, -0.3_rp, 0.5_rp, 1.1_rp-0.05_rp*ik ]
         u(:, :, ik) = unitary4(0.2_rp*ik, 0.5_rp*ik, 0.9_rp)
         va(:, :, ik) = herm4(real(ik, rp))
         vb(:, :, ik) = herm4(-0.4_rp*real(ik, rp))
      end do

      call moment_onsite_block(ev, u, va, vb, a, b, ioff, nb, ll, mu)

      maxerr = 0.0_rp
      do nn = 1, ll
         do m = 1, ll
            oref = (0.0_rp, 0.0_rp)
            do ik = 1, nk
               hk = hbuild(u(:, :, ik), ev(:, ik))
               call matrix_chebyshev(hk, a, b, ll, tmat)
               oref = oref + matmul(tmat(:, :, m), matmul(va(:, :, ik), &
                              matmul(tmat(:, :, nn), vb(:, :, ik))))
            end do
            oref = oref/real(nk, rp)
            do i = 1, nb
               do j = 1, nb
                  maxerr = max(maxerr, abs(mu(i, j, nn, m) - oref(ioff + i, ioff + j)))
               end do
            end do
         end do
      end do
      write (*, '(a,es12.4)') 'Test 3 (multi-site block extract)   max_err = ', maxerr
      if (maxerr > tol) failed = .true.
   end subroutine test_site_block

   ! ---- small builders ------------------------------------------------------

   function hbuild(u, ev) result(h)
      complex(rp), intent(in) :: u(:, :)
      real(rp), intent(in) :: ev(:)
      complex(rp) :: h(size(u, 1), size(u, 1))
      complex(rp) :: ud(size(u, 1), size(u, 1))
      integer :: i
      ud = conjg(transpose(u))
      do i = 1, size(u, 1)
         ud(i, :) = ev(i)*ud(i, :)
      end do
      h = matmul(u, ud)
      ! symmetrize away round-off so H is exactly Hermitian
      h = 0.5_rp*(h + conjg(transpose(h)))
   end function hbuild

   function herm2(d1, d2, re, im) result(v)
      real(rp), intent(in) :: d1, d2, re, im
      complex(rp) :: v(2, 2)
      v(1, 1) = cmplx(d1, 0.0_rp, rp)
      v(2, 2) = cmplx(d2, 0.0_rp, rp)
      v(1, 2) = cmplx(re, im, rp)
      v(2, 1) = cmplx(re, -im, rp)
   end function herm2

   function herm4(s) result(v)
      real(rp), intent(in) :: s
      complex(rp) :: v(4, 4)
      integer :: i, j
      v = (0.0_rp, 0.0_rp)
      do i = 1, 4
         v(i, i) = cmplx(0.3_rp*i + 0.1_rp*s, 0.0_rp, rp)
         do j = i + 1, 4
            v(i, j) = cmplx(0.2_rp*(i - j) + 0.05_rp*s, 0.1_rp*(i + j) - 0.03_rp*s, rp)
            v(j, i) = conjg(v(i, j))
         end do
      end do
   end function herm4

   !> 4x4 unitary as a block-diagonal pair of 2x2 unitaries mixed by a swap-ish
   !> rotation, giving a genuinely dense unitary (columns orthonormal).
   function unitary4(t1, t2, t3) result(u)
      real(rp), intent(in) :: t1, t2, t3
      complex(rp) :: u(4, 4), a(4, 4), c(4, 4)
      integer :: i
      a = (0.0_rp, 0.0_rp)
      a(1:2, 1:2) = unitary2(t1, 0.3_rp)
      a(3:4, 3:4) = unitary2(t2, -0.6_rp)
      ! Mixing rotation between orbital 2 and 3 (Givens, real).
      c = (0.0_rp, 0.0_rp)
      do i = 1, 4
         c(i, i) = (1.0_rp, 0.0_rp)
      end do
      c(2, 2) = cmplx(cos(t3), 0.0_rp, rp); c(3, 3) = cmplx(cos(t3), 0.0_rp, rp)
      c(2, 3) = cmplx(-sin(t3), 0.0_rp, rp); c(3, 2) = cmplx(sin(t3), 0.0_rp, rp)
      u = matmul(a, c)
   end function unitary4

end program test_moment_kernel
