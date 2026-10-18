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

   public :: gbt_kspace_h1

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

      real(rp), parameter :: two_pi = 2.0_rp*acos(-1.0_rp)
      integer :: ibond, a, b, s, t, sigma, i, j
      real(rp) :: angle
      complex(rp) :: phase(2), link

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
         angle = two_pi*dot_product(k - 0.5_rp*q, bond_d(:, ibond))
         phase(1) = cmplx(cos(angle), sin(angle), rp)
         angle = two_pi*dot_product(k + 0.5_rp*q, bond_d(:, ibond))
         phase(2) = cmplx(cos(angle), sin(angle), rp)
         do s = 1, 2
            do t = 1, 2
               link = cmplx(0.0_rp, 0.0_rp, rp)
               do sigma = 1, 2
                  link = link + conjg(frame(sigma, s, a))*phase(sigma)*frame(sigma, t, b)
               end do
               do j = 1, norb
                  do i = 1, norb
                     h(ix(a, s, i), ix(b, t, j)) = h(ix(a, s, i), ix(b, t, j)) + &
                        w(i, s, a)*link*sblock(i, j, ibond)*w(j, t, b)
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

end module gbt_kspace_mod
