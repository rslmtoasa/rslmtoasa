!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: bdg
!
!> DESCRIPTION:
!> Bogoliubov-de Gennes kernels on plain arrays (k-space, on-site singlet pairing).
!> Conventions (oracles: tests/bdg/oracles/README.md):
!>
!>   H_BdG(k) = [[H(k) - mu, D], [D^dag, -H(-k)^* + mu]],  D(iL up, iL dn) = +Delta_iL, D(iL dn, iL up) = -Delta_iL
!>   chi_iL   = sum_k w_k sum_n u_{iL up,n} conj(v_{iL dn,n}) (1 - f(E_n)),  Delta_iL = sum_L' K_LL' chi_iL'
!>
!> Matrix index a runs over nmat states; site_of, orb_of and spin_of (1 = up, 2 = dn) map it to (i, L, sigma).
!------------------------------------------------------------------------------
module bdg_mod
   use precision_mod, only: rp
   use math_mod, only: hermitian_eigensystem
   use reciprocal_mod, only: fermi_dirac_occupation
   implicit none
   private

   public :: bdg_hamiltonian_source, bdg_assemble, bdg_accumulate, bdg_gap, bdg_iterate

   !> @brief Supplies the normal-state Hamiltonian at k and -k for one k index.
   type, abstract :: bdg_hamiltonian_source
   contains
      procedure(hk_interface), deferred :: hk
   end type bdg_hamiltonian_source

   abstract interface
      subroutine hk_interface(this, ik, hk, hmk)
         import :: bdg_hamiltonian_source, rp
         class(bdg_hamiltonian_source), intent(in) :: this
         integer, intent(in) :: ik
         complex(rp), intent(out) :: hk(:, :), hmk(:, :)
      end subroutine hk_interface
   end interface

contains

   !> @brief Index of the spin-down partner of an up state with the same site and orbital.
   pure integer function partner(a, site_of, orb_of, spin_of)
      integer, intent(in) :: a, site_of(:), orb_of(:), spin_of(:)
      do partner = 1, size(site_of)
         if (spin_of(partner) == 2 .and. site_of(partner) == site_of(a) .and. orb_of(partner) == orb_of(a)) return
      end do
   end function partner

   !> @brief Builds H_BdG(k) (2 nmat x 2 nmat) from H(k), H(-k), mu and Delta(norb, nsite).
   subroutine bdg_assemble(hk, hmk, mu, delta, site_of, orb_of, spin_of, hbdg)
      complex(rp), intent(in) :: hk(:, :), hmk(:, :), delta(:, :)
      real(rp), intent(in) :: mu
      integer, intent(in) :: site_of(:), orb_of(:), spin_of(:)
      complex(rp), intent(out) :: hbdg(:, :)
      integer :: n, a, b

      n = size(hk, 1)
      hbdg = (0.0_rp, 0.0_rp)
      hbdg(1:n, 1:n) = hk
      hbdg(n + 1:2*n, n + 1:2*n) = -conjg(hmk)
      do a = 1, n
         hbdg(a, a) = hbdg(a, a) - mu
         hbdg(n + a, n + a) = hbdg(n + a, n + a) + mu
         if (spin_of(a) /= 1) cycle
         b = partner(a, site_of, orb_of, spin_of)
         hbdg(a, n + b) = delta(orb_of(a), site_of(a))
         hbdg(b, n + a) = -delta(orb_of(a), site_of(a))
         hbdg(n + b, a) = conjg(delta(orb_of(a), site_of(a)))
         hbdg(n + a, b) = -conjg(delta(orb_of(a), site_of(a)))
      end do
   end subroutine bdg_assemble

   !> @brief Adds one k point (weight wk) to chi(norb, nsite) and to the electron counts nel(2) = (up, dn).
   subroutine bdg_accumulate(evals, evecs, wk, kT, site_of, orb_of, spin_of, chi, nel)
      real(rp), intent(in) :: evals(:), wk, kT
      complex(rp), intent(in) :: evecs(:, :)
      integer, intent(in) :: site_of(:), orb_of(:), spin_of(:)
      complex(rp), intent(inout) :: chi(:, :)
      real(rp), intent(inout) :: nel(2)
      integer :: n, a, b, m
      real(rp) :: f

      n = size(site_of)
      do m = 1, size(evals)
         f = fermi_dirac_occupation(evals(m), 0.0_rp, kT)
         do a = 1, n
            nel(spin_of(a)) = nel(spin_of(a)) + wk*f*abs(evecs(a, m))**2
            if (spin_of(a) /= 1) cycle
            b = partner(a, site_of, orb_of, spin_of)
            chi(orb_of(a), site_of(a)) = chi(orb_of(a), site_of(a)) &
                                         + wk*evecs(a, m)*conjg(evecs(n + b, m))*(1.0_rp - f)
         end do
      end do
   end subroutine bdg_accumulate

   !> @brief Gap equation: delta(:, i) = K(type_of(i)) chi(:, i), with kmat(norb, norb, ntype) in Ry.
   subroutine bdg_gap(kmat, type_of, chi, delta)
      real(rp), intent(in) :: kmat(:, :, :)
      integer, intent(in) :: type_of(:)
      complex(rp), intent(in) :: chi(:, :)
      complex(rp), intent(out) :: delta(:, :)
      integer :: i

      do i = 1, size(type_of)
         delta(:, i) = matmul(kmat(:, :, type_of(i)), chi(:, i))
      end do
   end subroutine bdg_gap

   !> @brief Self-consistent loop in Delta: pass over all k, gap equation, linear mixing.
   !> @details max_iter = 0 does one pass at the given delta and leaves delta unchanged. chi and nel are those of
   !>          the last pass; dmax is max|delta_new - delta_old| of the last iteration.
   subroutine bdg_iterate(src, wk, mu, kT, kmat, type_of, site_of, orb_of, spin_of, mix, tol, max_iter, &
                          delta, chi, nel, niter, dmax)
      class(bdg_hamiltonian_source), intent(in) :: src
      real(rp), intent(in) :: wk(:), mu, kT, kmat(:, :, :), mix, tol
      integer, intent(in) :: type_of(:), site_of(:), orb_of(:), spin_of(:), max_iter
      complex(rp), intent(inout) :: delta(:, :)
      complex(rp), intent(out) :: chi(:, :)
      real(rp), intent(out) :: nel(2), dmax
      integer, intent(out) :: niter
      complex(rp), allocatable :: hk(:, :), hmk(:, :), hbdg(:, :), delta_new(:, :)
      real(rp), allocatable :: evals(:)
      integer :: n, ik

      n = size(site_of)
      allocate (hk(n, n), hmk(n, n), hbdg(2*n, 2*n), evals(2*n), delta_new(size(delta, 1), size(delta, 2)))
      niter = 0
      dmax = 0.0_rp
      do
         chi = (0.0_rp, 0.0_rp)
         nel = 0.0_rp
         do ik = 1, size(wk)
            call src%hk(ik, hk, hmk)
            call bdg_assemble(hk, hmk, mu, delta, site_of, orb_of, spin_of, hbdg)
            call hermitian_eigensystem(hbdg, evals)
            call bdg_accumulate(evals, hbdg, wk(ik), kT, site_of, orb_of, spin_of, chi, nel)
         end do
         if (max_iter == 0) exit
         call bdg_gap(kmat, type_of, chi, delta_new)
         dmax = maxval(abs(delta_new - delta))
         delta = delta + mix*(delta_new - delta)
         niter = niter + 1
         if (dmax < tol .or. niter >= max_iter) exit
      end do
   end subroutine bdg_iterate

end module bdg_mod
