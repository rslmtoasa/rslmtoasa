!> k-space generalized-Bloch-theorem (GBT) kernels on plain arrays.
!>
!> Relation R1 (Herring 1966; Sandratskii, Adv. Phys. 47, 91 (1998)): without SOC a spiral is invariant under a
!> lattice translation combined with a spin rotation by q.R, so in the local frames
!>   H_{a s,b s'}(k) = C_{a s} delta_ab delta_ss' + W_{a s} [U_a^dag diag(S_ab(k - q/2), S_ab(k + q/2)) U_b]_{s s'} W_{b s'},
!>   S_ab(k) = sum_d exp(2 pi i k.d) S_ab(d),   d = R + tau_b - tau_a (atomic gauge).
!> The lab spin-up component carries k - q/2 because D(alpha) = diag(e^{-i alpha/2}, e^{+i alpha/2}) gives it the
!> phase e^{-i alpha/2}. k - q/2 and k + q/2 are formed from the k given and never folded separately.
module gbt_kspace_mod
   use precision_mod, only: rp
   implicit none
   private

   public :: gbt_kspace_h1, gbt_kspace_h2, gbt_kspace_hcc, grand_potential, gbt_kspace_horth, gbt_newton_field

contains

   !> First-order GBT Hamiltonian at one k. Basis index of (site a, local spin s, orbital i) is
   !> (a - 1)*2*norb + (s - 1)*norb + i (site-major, spin-blocked). The result is in the orbital basis of `sblock`;
   !> a change of orbital basis is the caller's business.
   !> @param[in] bond_d      bond vectors d in direct (lattice) coordinates of the cell; onsite is the bond a = b, d = 0
   !> @param[in] sblock      raw structure-constant block of each bond, sblock(i, j) = sbar(j, i, m, ino)
   !> @param[in] frame       site frames U_a
   !> @param[in] w, c        W and C per orbital, local spin (1 majority, 2 minority) and site
   !> @param[in] k, q        k and q in direct coordinates
   pure subroutine gbt_kspace_h1(norb, nsite, nbond, bond_a, bond_b, bond_d, sblock, frame, w, c, k, q, h)
      integer, intent(in) :: norb, nsite, nbond
      integer, intent(in) :: bond_a(nbond), bond_b(nbond)
      real(rp), intent(in) :: bond_d(3, nbond)
      complex(rp), intent(in) :: sblock(norb, norb, nbond)
      complex(rp), intent(in) :: frame(2, 2, nsite)
      complex(rp), intent(in) :: w(norb, 2, nsite), c(norb, 2, nsite)
      real(rp), intent(in) :: k(3), q(3)
      complex(rp), intent(out) :: h(2*norb*nsite, 2*norb*nsite)

      integer :: ibond, a, b, s, t, i, j
      complex(rp) :: link(2, 2)

      h = cmplx(0.0_rp, 0.0_rp, rp)
      do a = 1, nsite
         do s = 1, 2
            do i = 1, norb
               h(ix(a, s, i), ix(a, s, i)) = c(i, s, a)
            end do
         end do
      end do

      do ibond = 1, nbond
         a = bond_a(ibond)
         b = bond_b(ibond)
         call gbt_bond_link(frame(:, :, a), frame(:, :, b), k, q, bond_d(:, ibond), link)
         do s = 1, 2
            do t = 1, 2
               do j = 1, norb
                  do i = 1, norb
                     h(ix(a, s, i), ix(b, t, j)) = h(ix(a, s, i), ix(b, t, j)) + &
                        w(i, s, a)*link(s, t)*sblock(i, j, ibond)*w(j, t, b)
                  end do
               end do
            end do
         end do
      end do

   contains

      pure integer function ix(site, spin, orbital)
         integer, intent(in) :: site, spin, orbital

         ix = (site - 1)*2*norb + (spin - 1)*norb + orbital
      end function ix

   end subroutine gbt_kspace_h1

   !> Two-centre combined-correction block of the GBT Hamiltonian at one k, in the orbital basis and index order of `sblock`.
   !> The spin link multiplies the primitive S and Sdot bonds before they meet the endpoint factors, as build_ccor_pair_block_gbt:
   !>   H_{a s i, b t j} = lam link_st(d) w_{a s i} [ Sdot_ij + (ccd1_{a i} + ccd1_{b j}) S_ij ] w_{b t j}
   !> and the onsite term lam ccd0 w^2, with the scalar lam = VMT - E_lin. Hcc is added after the order step, never inside
   !> h - h o h.
   !> @param[in] sdotblock  sdotblock(i, j) = -avw^2 sdot(j, i, m, ino), the transposition and scaling of sblock
   !> @param[in] ccd        ccd(i, 0:1, site), the coefficients of build_ccor_coefficients
   !> @param[in] lam        VMT - E_lin
   pure subroutine gbt_kspace_hcc(norb, nsite, nbond, bond_a, bond_b, bond_d, sblock, sdotblock, frame, w, ccd, lam, k, q, h)
      integer, intent(in) :: norb, nsite, nbond
      integer, intent(in) :: bond_a(nbond), bond_b(nbond)
      real(rp), intent(in) :: bond_d(3, nbond)
      complex(rp), intent(in) :: sblock(norb, norb, nbond), sdotblock(norb, norb, nbond)
      complex(rp), intent(in) :: frame(2, 2, nsite)
      complex(rp), intent(in) :: w(norb, 2, nsite)
      real(rp), intent(in) :: ccd(norb, 0:1, nsite), lam
      real(rp), intent(in) :: k(3), q(3)
      complex(rp), intent(out) :: h(2*norb*nsite, 2*norb*nsite)

      integer :: ibond, a, b, s, t, i, j
      complex(rp) :: link(2, 2)

      h = cmplx(0.0_rp, 0.0_rp, rp)
      do a = 1, nsite
         do s = 1, 2
            do i = 1, norb
               h(ix(a, s, i), ix(a, s, i)) = lam*ccd(i, 0, a)*w(i, s, a)*w(i, s, a)
            end do
         end do
      end do

      do ibond = 1, nbond
         a = bond_a(ibond)
         b = bond_b(ibond)
         call gbt_bond_link(frame(:, :, a), frame(:, :, b), k, q, bond_d(:, ibond), link)
         do s = 1, 2
            do t = 1, 2
               do j = 1, norb
                  do i = 1, norb
                     h(ix(a, s, i), ix(b, t, j)) = h(ix(a, s, i), ix(b, t, j)) + link(s, t)*w(i, s, a)*w(j, t, b)* &
                        (lam*sdotblock(i, j, ibond) + (ccd(i, 1, a)*lam + ccd(j, 1, b)*lam)*sblock(i, j, ibond))
                  end do
               end do
            end do
         end do
      end do

   contains

      pure integer function ix(site, spin, orbital)
         integer, intent(in) :: site, spin, orbital

         ix = (site - 1)*2*norb + (spin - 1)*norb + orbital
      end function ix

   end subroutine gbt_kspace_hcc

   !> Spin link of one bond, link(s, t) = sum_sigma conj(U_a(sigma, s)) exp(2 pi i (k -+ q/2).d) U_b(sigma, t): sigma = 1 carries
   !> k - q/2 and sigma = 2 carries k + q/2. Shared by the S bond of the Hamiltonian and the S and Sdot bonds of the CCOR block.
   pure subroutine gbt_bond_link(frame_a, frame_b, k, q, d, link)
      complex(rp), intent(in) :: frame_a(2, 2), frame_b(2, 2)
      real(rp), intent(in) :: k(3), q(3), d(3)
      complex(rp), intent(out) :: link(2, 2)

      real(rp), parameter :: two_pi = 2.0_rp*acos(-1.0_rp)
      real(rp) :: angle
      complex(rp) :: phase(2)
      integer :: s, t, sigma

      angle = two_pi*dot_product(k - 0.5_rp*q, d)
      phase(1) = cmplx(cos(angle), sin(angle), rp)
      angle = two_pi*dot_product(k + 0.5_rp*q, d)
      phase(2) = cmplx(cos(angle), sin(angle), rp)
      do s = 1, 2
         do t = 1, 2
            link(s, t) = cmplx(0.0_rp, 0.0_rp, rp)
            do sigma = 1, 2
               link(s, t) = link(s, t) + conjg(frame_a(sigma, s))*phase(sigma)*frame_b(sigma, t)
            end do
         end do
      end do
   end subroutine gbt_bond_link

   !> Second-order Hamiltonian H = E_nu + h - h O h, the formula of the collinear assembler. O and E_nu are onsite blocks in
   !> the local frame of each site (O = obar, E_nu = enu), and h1, obar and enu must be in the same orbital basis. The
   !> onsite C inside h1 is the one the collinear route uses at second order (cex, not cx).
   pure subroutine gbt_kspace_h2(norb, nsite, h1, obar, enu, h2)
      integer, intent(in) :: norb, nsite
      complex(rp), intent(in) :: h1(2*norb*nsite, 2*norb*nsite)
      complex(rp), intent(in) :: obar(2*norb, 2*norb, nsite), enu(2*norb, 2*norb, nsite)
      complex(rp), intent(out) :: h2(2*norb*nsite, 2*norb*nsite)

      complex(rp), allocatable :: ho(:, :)
      integer :: a, lo, hi

      allocate (ho(2*norb*nsite, 2*norb*nsite))
      do a = 1, nsite
         lo = (a - 1)*2*norb + 1
         hi = a*2*norb
         ho(:, lo:hi) = matmul(h1(:, lo:hi), obar(:, :, a))
      end do
      h2 = h1 - matmul(ho, h1)
      do a = 1, nsite
         lo = (a - 1)*2*norb + 1
         hi = a*2*norb
         h2(lo:hi, lo:hi) = h2(lo:hi, lo:hi) + enu(:, :, a)
      end do
   end subroutine gbt_kspace_h2

   !> Mermin grand potential per cell (relation R11, Mermin 1965): sum_k w_k sum_n -kT ln(1 + exp(-(e - mu)/kT)), written
   !> as (e - mu) - kT ln(1 + exp((e - mu)/kT)) for e < mu so that the exponential never overflows. kT > 0 is
   !> required. One state holds one electron (spinor bands), so there is no spin factor.
   pure function grand_potential(nband, nk, eigenvalues, weights, mu, kT) result(omega)
      integer, intent(in) :: nband, nk
      real(rp), intent(in) :: eigenvalues(nband, nk), weights(nk), mu, kT
      real(rp) :: omega

      real(rp) :: x
      integer :: ik, ib

      omega = 0.0_rp
      do ik = 1, nk
         do ib = 1, nband
            x = eigenvalues(ib, ik) - mu
            if (x < 0.0_rp) then
               omega = omega + weights(ik)*(x - kT*log(1.0_rp + exp(x/kT)))
            else
               omega = omega - weights(ik)*kT*log(1.0_rp + exp(-x/kT))
            end if
         end do
      end do
   end function grand_potential

   !> Exact-ASA Hamiltonian of a rigid spiral in the orthogonal representation (relation C1: Andersen & Jepsen, Phys. Rev.
   !> Lett. 53, 2571 (1984); Turek et al., Electronic Structure of Disordered Alloys, Surfaces and Interfaces (1997)):
   !>   S^gamma = S^alpha [1 - Gamma S^alpha]^-1 = [1 - S^alpha Gamma]^-1 S^alpha,   H = C + dele S^gamma dele.
   !> srot is S^alpha_rot, which gbt_kspace_h1 returns with w = 1 and c = 0; dele and gam (= gamma - alpha) are the local-frame
   !> diagonals in the index order of h1, cc is C as a matrix in the same order (diagonal, plus an onsite spin block per site when a
   !> constraining field is added). S^gamma comes from one LU solve, not an explicit inverse. h is returned as computed, not
   !> symmetrised: its anti-Hermitian part is the caller's to inspect. The eigenvalues are the roots of det(P^alpha(e) - S^alpha) = 0.
   subroutine gbt_kspace_horth(n, srot, cc, dele, gam, h)
      integer, intent(in) :: n
      complex(rp), intent(in) :: srot(n, n)
      complex(rp), intent(in) :: cc(n, n)
      real(rp), intent(in) :: dele(n), gam(n)
      complex(rp), intent(out) :: h(n, n)

      complex(rp) :: a(n, n)
      integer :: ipiv(n), i, j, info
      external :: zgesv

      do j = 1, n
         a(:, j) = -srot(:, j)*gam(j)
         a(j, j) = a(j, j) + 1.0_rp
      end do
      h = srot
      call zgesv(n, n, a, n, ipiv, h, n, info)
      if (info /= 0) error stop 'gbt_kspace_horth: 1 - S Gamma is singular'
      do j = 1, n
         do i = 1, n
            h(i, j) = dele(i)*h(i, j)*dele(j) + cc(i, j)
         end do
      end do
   end subroutine gbt_kspace_horth

   !> Constant-Jacobian Newton step of the constraining field of one site: b_new = b_old - alpha chi^-1 m_perp. chi = dm_perp/dB
   !> (2x2, local frame, negative definite for +B.sigma), b and m_perp are the local x, y components.
   pure subroutine gbt_newton_field(chi, b_old, m_perp, alpha, b_new)
      real(rp), intent(in) :: chi(2, 2), b_old(2), m_perp(2), alpha
      real(rp), intent(out) :: b_new(2)
      real(rp) :: det

      det = chi(1, 1)*chi(2, 2) - chi(1, 2)*chi(2, 1)
      b_new(1) = b_old(1) - alpha*(chi(2, 2)*m_perp(1) - chi(1, 2)*m_perp(2))/det
      b_new(2) = b_old(2) - alpha*(chi(1, 1)*m_perp(2) - chi(2, 1)*m_perp(1))/det
   end subroutine gbt_newton_field

end module gbt_kspace_mod
