!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: spin_response
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Transverse spin-response kernels on plain arrays (stage 1a: C0 only).
!> Conventions are those of cleanup/spin_response_restart_spec.md:
!>
!>   chi0_ij(w) = sum_k w_k sum_{n,m} (f_nk - f_m,k+q) T_i T_j^* / (w + e_nk - e_m,k+q + i eta)
!>   T_i^{nm}   = sum_mu conj(a^{nk}_{i mu up}) a^{m,k+q}_{i mu dn}
!>   (I + chi0 U) chi = chi0,   chi0(0,0) U M = -M,   L = -(chi - chi^H)/(2 pi i)
!>
!> Spin index of the amplitudes: 1 = up (majority), 2 = down. The module has no
!> dependence on the LMTO objects.
!------------------------------------------------------------------------------
module spin_response_mod
   use precision_mod, only: rp
   use math_mod, only: pi, i_unit
   implicit none
   private

   public :: accumulate_chi0
   public :: u_juelich
   public :: u_mills
   public :: solve_dyson
   public :: spectral_trace
   public :: pole_estimates

   !> Energy difference (Ry) below which the static value uses the Fermi divided difference.
   real(rp), parameter :: degenerate_tol = 1.0e-10_rp

contains

   !> @brief Add the Lehmann sum of one k set to chi0 (the caller zeroes chi0 and may chunk k).
   !> @details eta = 0 is valid only at omega = 0, where a pair with |e_n - e_m| < 1e-10 Ry
   !>          takes the Fermi divided difference -f(1-f)/kt. Every other pair is pruned only
   !>          when f_n - f_m is exactly zero.
   !> @param[in]    kt     Temperature (Ry).
   !> @param[in]    w      k weights, shape (nk).
   !> @param[in]    e_k    Energies at k, shape (nk, nb); f_k occupations, same shape.
   !> @param[in]    a_k    d amplitudes at k, shape (nk, nb, norb, nsite, 2).
   !> @param[in]    e_kq   Energies at k+q, shape (nk, nb); f_kq, a_kq likewise.
   !> @param[in]    omega  Energies, shape (nw) (Ry).
   !> @param[in]    eta    Broadening (Ry).
   !> @param[inout] chi0   Bare response, shape (nsite, nsite, nw).
   pure subroutine accumulate_chi0(kt, w, e_k, f_k, a_k, e_kq, f_kq, a_kq, omega, eta, chi0)
      real(rp), intent(in) :: kt, w(:), e_k(:, :), f_k(:, :), e_kq(:, :), f_kq(:, :), omega(:), eta
      complex(rp), intent(in) :: a_k(:, :, :, :, :), a_kq(:, :, :, :, :)
      complex(rp), intent(inout) :: chi0(:, :, :)

      integer :: ik, n, m, i, j, iw
      real(rp) :: df, de
      complex(rp) :: t(size(a_k, 4)), c

      do ik = 1, size(w)
         do n = 1, size(e_k, 2)
            do m = 1, size(e_kq, 2)
               df = f_k(ik, n) - f_kq(ik, m)
               de = e_k(ik, n) - e_kq(ik, m)
               do i = 1, size(t)
                  t(i) = sum(conjg(a_k(ik, n, :, i, 1))*a_kq(ik, m, :, i, 2))
               end do
               do iw = 1, size(omega)
                  if (eta == 0.0_rp .and. omega(iw) == 0.0_rp .and. abs(de) < degenerate_tol) then
                     c = cmplx(-f_k(ik, n)*(1.0_rp - f_k(ik, n))/kt, 0.0_rp, rp)
                  else if (df == 0.0_rp) then
                     cycle
                  else
                     c = df/(omega(iw) + de + i_unit*eta)
                  end if
                  do j = 1, size(t)
                     do i = 1, size(t)
                        chi0(i, j, iw) = chi0(i, j, iw) + w(ik)*c*t(i)*conjg(t(j))
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine accumulate_chi0

   !> @brief Juelich U from the Goldstone condition sum_j chi0_ij(0,0) M_j U_j = -M_i.
   !> @param[in]  chi0s  Static bare response chi0(0,0), real, shape (nsite, nsite).
   !> @param[in]  m      Site moments, shape (nsite).
   !> @param[out] u      Interaction (Ry), shape (nsite).
   subroutine u_juelich(chi0s, m, u)
      real(rp), intent(in) :: chi0s(:, :), m(:)
      real(rp), intent(out) :: u(:)

      real(rp) :: amat(size(m), size(m))
      integer :: j, info, ipiv(size(m))

      do j = 1, size(m)
         amat(:, j) = chi0s(:, j)*m(j)
      end do
      u = -m
      call dgesv(size(m), 1, amat, size(m), ipiv, u, size(m), info)
      if (info /= 0) error stop 'u_juelich: singular Goldstone matrix'
   end subroutine u_juelich

   !> @brief Mills U = (C_dn - C_up)/M per site, from the d band centres.
   elemental function u_mills(c_up, c_dn, m) result(u)
      real(rp), intent(in) :: c_up, c_dn, m
      real(rp) :: u

      u = (c_dn - c_up)/m
   end function u_mills

   !> @brief Solve (I + chi0 U) chi = chi0 at each omega by LU, with U diagonal in sites.
   !> @details Never call at eta = 0 for q = 0: the Goldstone condition makes I + chi0 U singular there.
   !> @param[in]  chi0  Bare response, shape (nsite, nsite, nw).
   !> @param[in]  u     Interaction (Ry), shape (nsite).
   !> @param[out] chi   Enhanced response, shape (nsite, nsite, nw).
   subroutine solve_dyson(chi0, u, chi)
      complex(rp), intent(in) :: chi0(:, :, :)
      real(rp), intent(in) :: u(:)
      complex(rp), intent(out) :: chi(:, :, :)

      integer :: ns, iw, i, j, info, ipiv(size(u))
      complex(rp) :: amat(size(u), size(u))

      ns = size(u)
      do iw = 1, size(chi0, 3)
         do j = 1, ns
            amat(:, j) = chi0(:, j, iw)*u(j)
            amat(j, j) = amat(j, j) + 1.0_rp
         end do
         chi(:, :, iw) = chi0(:, :, iw)
         call zgesv(ns, ns, amat, ns, ipiv, chi(:, :, iw), ns, info)
         if (info /= 0) error stop 'solve_dyson: singular I + chi0 U'
      end do
   end subroutine solve_dyson

   !> @brief tr L(omega) with L = -(chi - chi^H)/(2 pi i), i.e. -Im tr chi / pi.
   pure function spectral_trace(chi) result(trl)
      complex(rp), intent(in) :: chi(:, :, :)
      real(rp) :: trl(size(chi, 3))

      integer :: iw, i

      do iw = 1, size(chi, 3)
         trl(iw) = 0.0_rp
         do i = 1, size(chi, 1)
            trl(iw) = trl(iw) - aimag(chi(i, i, iw))/pi
         end do
      end do
   end function spectral_trace

   !> @brief One-site pole estimates: parabolic peak of tr L and first zero of Re(1 + chi0 U) with interpolated omega > 0.
   !> @details A maximum on either grid end gives has_peak = .false.; no sign change gives
   !>          has_crossing = .false. peak and crossing are undefined unless their flag is set.
   !> @param[in]  omega  Increasing grid, shape (nw).
   !> @param[in]  trl    tr L on the grid, shape (nw).
   !> @param[in]  chi0   Bare response on the grid, shape (nw).
   !> @param[in]  u      Interaction (Ry).
   pure subroutine pole_estimates(omega, trl, chi0, u, peak, crossing, has_peak, has_crossing)
      real(rp), intent(in) :: omega(:), trl(:), u
      complex(rp), intent(in) :: chi0(:)
      real(rp), intent(out) :: peak, crossing
      logical, intent(out) :: has_peak, has_crossing

      integer :: imax, iw
      real(rp) :: h, g0, g1, root

      peak = 0.0_rp
      crossing = 0.0_rp
      imax = maxloc(trl, 1)
      has_peak = imax > 1 .and. imax < size(trl)
      if (has_peak) then
         h = 0.5_rp*(omega(imax + 1) - omega(imax - 1))
         peak = omega(imax) + h*(trl(imax - 1) - trl(imax + 1)) &
                /(2.0_rp*(trl(imax - 1) - 2.0_rp*trl(imax) + trl(imax + 1)))
      end if

      has_crossing = .false.
      do iw = 1, size(omega) - 1
         g0 = 1.0_rp + u*real(chi0(iw), rp)
         g1 = 1.0_rp + u*real(chi0(iw + 1), rp)
         if (g0*g1 < 0.0_rp) then
            root = omega(iw) - g0*(omega(iw + 1) - omega(iw))/(g1 - g0)
            if (root > 0.0_rp) then
               crossing = root
               has_crossing = .true.
               exit
            end if
         end if
      end do
   end subroutine pole_estimates

end module spin_response_mod
