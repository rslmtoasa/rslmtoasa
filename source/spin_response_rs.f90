!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: spin_response_rs
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Real-space transverse spin response. Static bare susceptibility of the d channel per pair (i, j),
!> with the coefficient vertex sum_mu c+_{i mu up} c_{i mu dn}:
!>
!>   chi0_ij = -(X - conj(Y))/(2 pi i),  X = int f Tr_d[G^upup_ji G^dndn_ij] dE,  Y = int f Tr_d[G^dndn_ji G^upup_ij] dE
!>
!> with retarded G. Only the spin-diagonal blocks of the spinor G enter, also with SOC; with time reversal X = Y and
!> chi0_ij = -Im X/pi. Sign and units (1/Ry) are those of the k-space spin_response with method = 'mills',
!> channels = 'd'. The kernels take plain arrays; the exchange pair loop (calculate_pairwise) calls them.
!>
!> Goldstone diagnostics per atom type a, from the first self pair (i, i) with type a (the seed of a): its d moment
!> m_a = -(1/pi) Im int f Tr_d[G^upup_ii - G^dndn_ii] dE, and chi0(q = 0)_ab = sum of Re chi0_ij over the listed pairs
!> of that seed with type(j) = b (truncated to the pair list). As in spin_response: U from the Goldstone condition
!> sum_b chi0_ab U_b m_b = -m_a (by construction), U_Mills = (C_dn - C_up)_d/m and the Mills residual.
!------------------------------------------------------------------------------
module spin_response_rs_mod
   use precision_mod, only: rp
   use math_mod, only: pi, i_unit
   use basis_mod, only: spin_off
   use spin_response_mod, only: u_juelich, u_mills, mills_residual
   implicit none
   private

   public :: chi0_traces
   public :: chi0_from_traces
   public :: write_chi0_pairs
   public :: d_moment_trace
   public :: type_sums
   public :: write_chi0_types

   !> Angular momentum of the response channel.
   integer, parameter :: l_d = 2
   !> |m| (mu_B) below which a type is treated as non-magnetic and the U columns are n/a.
   real(rp), parameter :: m_tol = 1.0e-6_rp

contains

   !> @brief Tr_d[G^upup_ji G^dndn_ij] (column 1) and Tr_d[G^dndn_ji G^upup_ij] (column 2) for each index n of the third
   !>        dimension (energy point or moment).
   !> @param[in] gij, gji  Spinor blocks G_ij, G_ji, shape (nb, nb, n), spin-down block from spin_off + 1.
   pure function chi0_traces(gij, gji) result(t)
      complex(rp), intent(in) :: gij(:, :, :), gji(:, :, :)
      complex(rp) :: t(size(gij, 3), 2)

      integer :: n, u0, u1, d0, d1

      u0 = l_d**2 + 1
      u1 = (l_d + 1)**2
      d0 = u0 + spin_off
      d1 = u1 + spin_off
      do n = 1, size(gij, 3)
         t(n, 1) = sum(transpose(gji(u0:u1, u0:u1, n))*gij(d0:d1, d0:d1, n))
         t(n, 2) = sum(transpose(gji(d0:d1, d0:d1, n))*gij(u0:u1, u0:u1, n))
      end do
   end function chi0_traces

   !> @brief chi0_ij = -(x - conj(y))/(2 pi i) from the energy integrals x, y of the two columns of chi0_traces.
   elemental function chi0_from_traces(x, y) result(chi0)
      complex(rp), intent(in) :: x, y
      complex(rp) :: chi0

      chi0 = -(x - conjg(y))/(2.0_rp*pi*i_unit)
   end function chi0_from_traces

   !> @brief Write one row per pair: atoms, types, r_j - r_i and |r_j - r_i| (cluster coordinates), Re and Im chi0_ij.
   !> @param[in] ijpair  Pair list, shape (npair, 2).
   !> @param[in] iz      Type of each cluster atom.
   !> @param[in] cr      Cluster coordinates, shape (3, natom).
   !> @param[in] chi0    chi0_ij (1/Ry), shape (npair).
   subroutine write_chi0_pairs(fname, ijpair, iz, cr, chi0)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: ijpair(:, :), iz(:)
      real(rp), intent(in) :: cr(:, :)
      complex(rp), intent(in) :: chi0(:)

      integer :: funit, p, i, j

      open (newunit=funit, file=fname, action='write', status='replace')
      write (funit, '(a)') '# i j type_i type_j r_ij(3, cluster coordinates) |r_ij| re_chi0(1/Ry) im_chi0(1/Ry)'
      do p = 1, size(chi0)
         i = ijpair(p, 1)
         j = ijpair(p, 2)
         write (funit, '(4i8,4f12.6,2es24.14)') i, j, iz(i), iz(j), cr(:, j) - cr(:, i), norm2(cr(:, j) - cr(:, i)), &
            real(chi0(p), rp), aimag(chi0(p))
      end do
      close (funit)
   end subroutine write_chi0_pairs

   !> @brief Tr_d[G^upup - G^dndn] of a self pair for each index n of the third dimension (energy point or moment).
   pure function d_moment_trace(g) result(t)
      complex(rp), intent(in) :: g(:, :, :)
      complex(rp) :: t(size(g, 3))

      integer :: n, mu

      do n = 1, size(g, 3)
         t(n) = (0.0_rp, 0.0_rp)
         do mu = l_d**2 + 1, (l_d + 1)**2
            t(n) = t(n) + g(mu, mu, n) - g(mu + spin_off, mu + spin_off, n)
         end do
      end do
   end function d_moment_trace

   !> @brief Type sums of the pair list. seed(a) is the first atom i of a self pair (i, i) of type a (0: none), m(a) its
   !>        d moment md, npair(a) the number of pairs with first atom seed(a) and chi0s(a, b) the sum of Re chi0 over
   !>        those pairs with iz(j) = b. Rows of types without a seed stay 0.
   !> @param[in] ijpair  Pair list, shape (npair, 2); iz the type of each cluster atom.
   !> @param[in] chi0    chi0_ij (1/Ry) and md the d moment of each pair (used for self pairs only), shape (npair).
   pure subroutine type_sums(ijpair, iz, ntype, chi0, md, seed, npair, m, chi0s)
      integer, intent(in) :: ijpair(:, :), iz(:), ntype
      complex(rp), intent(in) :: chi0(:)
      real(rp), intent(in) :: md(:)
      integer, intent(out) :: seed(ntype), npair(ntype)
      real(rp), intent(out) :: m(ntype), chi0s(ntype, ntype)

      integer :: p, a

      seed = 0
      npair = 0
      m = 0.0_rp
      chi0s = 0.0_rp
      do p = 1, size(chi0)
         a = iz(ijpair(p, 1))
         if (ijpair(p, 1) == ijpair(p, 2) .and. seed(a) == 0) then
            seed(a) = ijpair(p, 1)
            m(a) = md(p)
         end if
      end do
      do p = 1, size(chi0)
         a = iz(ijpair(p, 1))
         if (seed(a) == 0 .or. ijpair(p, 1) /= seed(a)) cycle
         npair(a) = npair(a) + 1
         chi0s(a, iz(ijpair(p, 2))) = chi0s(a, iz(ijpair(p, 2))) + real(chi0(p), rp)
      end do
   end subroutine type_sums

   !> @brief Write one row per type with a seed: type, seed atom, pair count, m_d, U_Goldstone, U_Mills, Mills residual and
   !>        chi0(q = 0)(a, b) for every b. The three U columns are n/a unless every type has a seed and |m| >= m_tol.
   !> @param[in] c_up, c_dn  d band centres (Ry) of each type, shape (ntype).
   subroutine write_chi0_types(fname, seed, npair, m, chi0s, c_up, c_dn)
      character(len=*), intent(in) :: fname
      integer, intent(in) :: seed(:), npair(:)
      real(rp), intent(in) :: m(:), chi0s(:, :), c_up(:), c_dn(:)

      integer :: funit, a
      real(rp) :: ug(size(m)), um(size(m)), delta(size(m))
      character(len=72) :: ucols
      logical :: have_u

      have_u = all(seed > 0) .and. all(abs(m) >= m_tol)
      if (have_u) then
         call u_juelich(chi0s, m, ug)
         um = u_mills(c_up, c_dn, m)
         delta = mills_residual(chi0s, um, m)
      end if
      open (newunit=funit, file=fname, action='write', status='replace')
      write (funit, '(a)') '# type seed_atom n_pairs m_d(mu_B) u_goldstone(Ry, by construction) u_mills(Ry) mills_residual '// &
         'chi0_q0(type, b)(1/Ry) for b = 1..ntype'
      do a = 1, size(m)
         if (seed(a) == 0) cycle
         if (have_u) then
            write (ucols, '(3es24.14)') ug(a), um(a), delta(a)
         else
            write (ucols, '(3a24)') 'n/a', 'n/a', 'n/a'
         end if
         write (funit, '(3i8,es24.14,a,*(es24.14))') a, seed(a), npair(a), m(a), ucols, chi0s(a, :)
      end do
      close (funit)
   end subroutine write_chi0_types

end module spin_response_rs_mod
