!> Full-spinor Liechtenstein + FLL functional for one d shell (10x10, not yet used by the SCF).
!>
!> Index: i = (sigma-1)*5 + (m+3), spin-major with m=-2..2. In the Hamiltonian this is
!> row (sigma-1)*spin_off + l*l + (m+l+1) with l=2. Operator orientation as in the
!> collinear code: P(a,b) = sum_n f_n psi_n(a) conj(psi_n(b)), V is inserted into H as is,
!> and dE = Re Tr(V dP) for Hermitian dP. The tensor is T(a,b,c,d) = <a c|b d> from
!> Coulomb_mat_spherical. The exchange term keeps all spin contractions.
module hubbard_spinor_mod
   use precision_mod, only: rp
   use math_mod, only: Coulomb_mat_spherical
   implicit none
   private
   public :: hubbard_d_tensor, hubbard_spinor_functional

   integer, parameter :: nd = 5

contains

   !> d-shell tensor with F0=U, F2=14J/1.625, F4=0.625 F2, as in the collinear routines.
   subroutine hubbard_d_tensor(u, jbar, t)
      real(rp), intent(in) :: u, jbar
      complex(rp), intent(out) :: t(nd, nd, nd, nd)
      real(rp) :: f2, f4
      integer :: a, b, c, d

      f2 = 14.0_rp*jbar/1.625_rp
      f4 = 0.625_rp*f2
      do d = 1, nd
         do c = 1, nd
            do b = 1, nd
               do a = 1, nd
                  t(a, b, c, d) = Coulomb_mat_spherical(2, a - 3, b - 3, c - 3, d - 3, u, f2, f4, 0.0_rp)
               end do
            end do
         end do
      end do
   end subroutine hubbard_d_tensor

   !> E_int = 1/2 sum T_abcd sum_{s,s'} (P_{b s,a s} P_{d s',c s'} - P_{d s',a s} P_{b s,c s'}),
   !> E_DC = U/2 N(N-1) - J/2 [Tr(D^2) - N] with D_{s s'} = sum_a P_{a s,a s'}
   !> (= U/2 N(N-1) - J/2 [N(N/2-1) + |m|^2/2]), E_Hub = E_int - E_DC, V = dE_Hub/dP.
   subroutine hubbard_spinor_functional(t, u, jbar, p, e_int, e_dc, e_hub, v)
      complex(rp), intent(in) :: t(nd, nd, nd, nd), p(2*nd, 2*nd)
      real(rp), intent(in) :: u, jbar
      real(rp), intent(out) :: e_int, e_dc, e_hub
      complex(rp), intent(out) :: v(2*nd, 2*nd)
      complex(rp) :: n(nd, nd), dmat(2, 2), direct, exchange, w
      real(rp) :: ntot, trd2, vdc
      integer :: a, b, c, d, s, r, ia, ib, ic, id

      n = (0.0_rp, 0.0_rp)
      dmat = (0.0_rp, 0.0_rp)
      do s = 1, 2
         n = n + p((s - 1)*nd + 1:s*nd, (s - 1)*nd + 1:s*nd)
      end do
      do s = 1, 2
         do r = 1, 2
            do a = 1, nd
               dmat(s, r) = dmat(s, r) + p((s - 1)*nd + a, (r - 1)*nd + a)
            end do
         end do
      end do
      ntot = real(dmat(1, 1) + dmat(2, 2), rp)

      direct = (0.0_rp, 0.0_rp)
      exchange = (0.0_rp, 0.0_rp)
      v = (0.0_rp, 0.0_rp)
      do d = 1, nd
         do c = 1, nd
            do b = 1, nd
               do a = 1, nd
                  w = t(a, b, c, d)
                  direct = direct + w*n(b, a)*n(d, c)
                  do s = 1, 2
                     do r = 1, 2
                        ia = (s - 1)*nd + a
                        ib = (s - 1)*nd + b
                        ic = (r - 1)*nd + c
                        id = (r - 1)*nd + d
                        exchange = exchange + w*p(id, ia)*p(ib, ic)
                        ! V_exchange(a r, d s) = -sum_bc T_abcd P(b r, c s)
                        v((r - 1)*nd + a, (s - 1)*nd + d) = v((r - 1)*nd + a, (s - 1)*nd + d) &
                           - w*p((r - 1)*nd + b, (s - 1)*nd + c)
                     end do
                     ! V_direct(a s, b s) = sum_cd T_abcd n(d, c)
                     v((s - 1)*nd + a, (s - 1)*nd + b) = v((s - 1)*nd + a, (s - 1)*nd + b) + w*n(d, c)
                  end do
               end do
            end do
         end do
      end do
      e_int = 0.5_rp*real(direct - exchange, rp)

      trd2 = real(dmat(1, 1)*dmat(1, 1) + dmat(1, 2)*dmat(2, 1) + dmat(2, 1)*dmat(1, 2) + dmat(2, 2)*dmat(2, 2), rp)
      e_dc = 0.5_rp*u*ntot*(ntot - 1.0_rp) - 0.5_rp*jbar*(trd2 - ntot)
      e_hub = e_int - e_dc

      ! V_DC(a s, a r) = [U(N-1/2) + J/2] delta_sr - J D_sr
      vdc = u*(ntot - 0.5_rp) + 0.5_rp*jbar
      do a = 1, nd
         do s = 1, 2
            do r = 1, 2
               v((s - 1)*nd + a, (r - 1)*nd + a) = v((s - 1)*nd + a, (r - 1)*nd + a) + jbar*dmat(s, r)
               if (s == r) v((s - 1)*nd + a, (s - 1)*nd + a) = v((s - 1)*nd + a, (s - 1)*nd + a) - vdc
            end do
         end do
      end do
   end subroutine hubbard_spinor_functional

end module hubbard_spinor_mod
