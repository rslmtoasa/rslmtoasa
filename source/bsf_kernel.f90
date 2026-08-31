!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: bsf_kernel
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Bloch spectral function (BSF) numerical core (milestone B3). Dependency-free
!> (only `precision_mod` + `math_mod`) so the sum-rule known-answer test needs no
!> LMTO machinery. The Bloch spectral function is the partial trace of the
!> retarded k-space Green’s function,
!>
!>   A(k,E) = -(1/pi) Im Tr_sub G(k, E + i*eta),
!>
!> where `G(k,z)` is delivered by `reciprocal_green` (either backend: with Sigma=0
!> backend E/D both give the clean-crystal resolvent; with Sigma!=0, backend D
!> gives the disorder/correlation-broadened A(k,E) for CPA/DMFT through the SAME
!> code path -- that is the B3 design point). `Tr_sub` is a partial trace over a
!> chosen set of diagonal orbital channels: the full set (total A), the spin-up or
!> spin-down sub-blocks, a single site’s block, or one orbital -- all expressed as
!> an index list into the nmat x nmat resolvent.
!------------------------------------------------------------------------------
module bsf_kernel_mod
   use precision_mod, only: rp
   use math_mod, only: pi
   implicit none
   private

   public :: bsf_spectral_trace

contains

   !> @brief Partial-trace Bloch spectral weight A = -(1/pi) sum_{i in idx} Im G_ii.
   !> @details The single BSF convention pin (B3.3 review): the retarded spectral
   !>          function is the negative imaginary part of the diagonal resolvent,
   !>          normalized by 1/pi, summed over the selected diagonal channels. The
   !>          caller chooses `idx` to select total (all channels), spin-up /
   !>          spin-down (the per-site up/down orbital ranges), a site block, or a
   !>          single orbital. Units: states / Ry / (selected channels).
   !> @param[in] gk   nmat x nmat retarded resolvent G(k, E+i*eta).
   !> @param[in] idx  Diagonal channel indices to include in the partial trace.
   !> @return A(k,E) spectral weight for the selected channels.
   pure function bsf_spectral_trace(gk, idx) result(a)
      complex(rp), intent(in) :: gk(:, :)
      integer, intent(in) :: idx(:)
      real(rp) :: a
      integer :: m

      a = 0.0_rp
      do m = 1, size(idx)
         a = a - aimag(gk(idx(m), idx(m)))
      end do
      a = a/pi
   end function bsf_spectral_trace

end module bsf_kernel_mod
