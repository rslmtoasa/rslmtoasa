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
!------------------------------------------------------------------------------
module spin_response_rs_mod
   use precision_mod, only: rp
   use math_mod, only: pi, i_unit
   use basis_mod, only: spin_off
   implicit none
   private

   public :: chi0_traces
   public :: chi0_from_traces
   public :: write_chi0_pairs

   !> Angular momentum of the response channel.
   integer, parameter :: l_d = 2

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

end module spin_response_rs_mod
