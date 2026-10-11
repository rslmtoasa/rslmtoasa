!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: alsda_radial
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Radial data of the ALSDA transverse kernel for spin_response. A response run reads only the potential file, so the
!> atomic sphere of each site is solved again with the production lmtst on a copy of the site (the potential parameters
!> of the lattice are not touched), and the radial arrays are returned as plain arrays.
!------------------------------------------------------------------------------
module alsda_radial_mod
   use precision_mod, only: rp
   use lattice_mod, only: lattice
   use symbolic_atom_mod, only: symbolic_atom
   use self_mod, only: self
   use radial_ground_state_mod, only: radial_simpson_weight
   implicit none
   private

   public :: site_radial_data

contains

   !> @brief Radial arrays of one site from a fresh lmtst on a copy of the site, in physical spin.
   !> @details Radial function of channel l and spin s: phi at E_nu if juelich is false, else the frozen projector
   !>          (phi + x phidot)/sqrt(1 + x^2 p), x = e_fermi - E_nu, p = 1/ppar^2 (gfac kept at E_nu).
   !> @param[out] r, wq       Mesh (bohr) and quadrature weights including dr/di.
   !> @param[out] m, dv       n_up - n_dn (core included) and V_xc,up - V_xc,dn.
   !> @param[out] enu         Linearization energies of the re-run (Ry), shape (nch, 2).
   !> @param[out] g_up .. gfac_dn  Radial functions per channel, shape (nr, nch).
   subroutine site_radial_data(lattice_obj, isite, channel_l, e_fermi, juelich, r, wq, m, dv, enu, g_up, g_dn, s_up, s_dn, gfac_up, gfac_dn)
      type(lattice), target, intent(inout) :: lattice_obj
      integer, intent(in) :: isite, channel_l(:)
      real(rp), intent(in) :: e_fermi
      logical, intent(in) :: juelich
      real(rp), allocatable, intent(out) :: r(:), wq(:), m(:), dv(:), enu(:, :), g_up(:, :), g_dn(:, :), s_up(:, :), s_dn(:, :), &
                                            gfac_up(:, :), gfac_dn(:, :)

      type(self) :: scf
      type(symbolic_atom) :: atom
      real(rp) :: qsl(6), x, p
      real(rp), allocatable :: g(:, :, :), s(:, :, :)
      integer :: i, ic, l, sp, nr

      scf%lattice => lattice_obj
      scf%control => lattice_obj%control
      call scf%restore_to_default()
      atom = lattice_obj%symbolic_atoms(lattice_obj%nbulk + isite)
      qsl = scf%lmtst(atom, isite)

      associate (rg => atom%radial_ground_state)
         nr = size(rg%r)
         r = rg%r
         wq = [(radial_simpson_weight(i, nr)*rg%a*(rg%r(i) + rg%b), i=1, nr)]
         m = rg%n_up - rg%n_down
         dv = rg%delta_vxc
         allocate (g(nr, size(channel_l), 2), s(nr, size(channel_l), 2), enu(size(channel_l), 2), gfac_up(nr, size(channel_l)), &
                   gfac_dn(nr, size(channel_l)))
         do ic = 1, size(channel_l)
            l = channel_l(ic)
            do sp = 1, 2
               enu(ic, sp) = rg%sr_enu(l + 1, sp)
               x = merge(e_fermi - enu(ic, sp), 0.0_rp, juelich)
               p = 1.0_rp/atom%potential%ppar(l, sp)**2
               g(:, ic, sp) = (rg%pauli_large(:, l + 1, sp) + x*rg%sr_large_dot(:, l + 1, sp))/sqrt(1.0_rp + x*x*p)
               s(:, ic, sp) = (rg%pauli_small(:, l + 1, sp) + x*rg%sr_small_dot(:, l + 1, sp))/sqrt(1.0_rp + x*x*p)
            end do
            gfac_up(:, ic) = rg%sr_gfac(:, l + 1, 1)
            gfac_dn(:, ic) = rg%sr_gfac(:, l + 1, 2)
         end do
      end associate
      g_up = g(:, :, 1)
      g_dn = g(:, :, 2)
      s_up = s(:, :, 1)
      s_dn = s(:, :, 2)
   end subroutine site_radial_data

end module alsda_radial_mod
