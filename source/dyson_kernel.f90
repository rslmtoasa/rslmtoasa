!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: dyson_kernel
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Pure-math kernel for the direct-Dyson (backend D) k-space Green's function
!> (milestone B2, `reciprocal_green`). Given the k-space Hamiltonian block
!> H(k) and a self-energy block Sigma(z), it forms and inverts one Dyson matrix
!>
!>   G(k,z) = [ z*S - H(k) - Sigma(z) ]^{-1}      per (k, z)
!>
!> Overlap convention (pinned, B2.4). Backend E (strict Lehmann) diagonalizes
!> the STANDARD Hermitian eigenproblem (`reciprocal%diagonalize_hamiltonian`
!> calls `zheev`, not `zhegv`), i.e. it implicitly assumes an orthonormal basis
!> S = I. The permanent CI invariant "backend D with Sigma = 0 == backend E"
!> therefore holds ONLY when backend D also uses S = I. This kernel hard-codes
!> S = I (the Dyson matrix is z*I - H - Sigma); the generalized-overlap case
!> (`reciprocal_mode = 'generalized_overlap_proxy'`, `sk_overlap`) is deliberately
!> out of scope for B2.4 because it would break the D == E equivalence unless
!> backend E is first re-cast as a generalized eigenproblem. Do not add S(k)
!> here without also generalizing backend E and re-deriving the invariant.
!>
!> Like `lehmann_kernel_mod`, this routine deliberately carries NO dependence on
!> the `reciprocal`, `lattice`, or `green` objects: it is the isolated numerical
!> core so the equivalence known-answer test
!> (`tests/unit/test_dyson_equivalence.f90`) can exercise it against the Lehmann
!> kernel with no LMTO machinery. `reciprocal_green::fill_green_dyson` wires the
!> H(k) tiles, the block-diagonal Sigma, the pair->site map, the bond vectors,
!> and the 1/N_k inverse Bloch transform around it -- inverting ONCE per (k, z)
!> and distributing the sub-blocks to every pair (design B2 §1.4: stream (k, z),
!> never materialize all (k, z) or re-invert per pair).
!------------------------------------------------------------------------------
module dyson_kernel_mod
   use precision_mod, only: rp
   implicit none

   private
   public :: dyson_kspace_inverse

contains

   !> @brief Invert one Dyson matrix G(k,z) = [z*I - H(k) - Sigma(z)]^{-1}.
   !> @details Forms A = z*I - hk - sigma_full (S = I, see module header) and
   !>          returns its full inverse via LAPACK `zgetrf`/`zgetri` -- the same
   !>          factor/inverse pair the real-space route uses (`green.f90`). The
   !>          caller extracts the site sub-blocks it needs; inverting the full
   !>          nmat x nmat matrix once per (k, z) amortizes over all pairs.
   !> @param[in]  hk         H(k) block, shape (nmat, nmat).
   !> @param[in]  z          Complex energy on the retarded contour (z = E + i*eta).
   !> @param[in]  sigma_full Self-energy Sigma(z), shape (nmat, nmat) (block
   !>                        diagonal over sites; zero for backend E).
   !> @param[out] gk         Inverse [z*I - hk - sigma_full]^{-1}, shape (nmat, nmat).
   subroutine dyson_kspace_inverse(hk, z, sigma_full, gk)
      complex(rp), intent(in) :: hk(:, :)
      complex(rp), intent(in) :: z
      complex(rp), intent(in) :: sigma_full(:, :)
      complex(rp), intent(out) :: gk(:, :)

      integer :: nmat, i, info, lwork
      integer, allocatable :: ipiv(:)
      complex(rp), allocatable :: work(:)

      nmat = size(hk, 1)

      ! A = z*I - H(k) - Sigma  (S = I; see module header).
      gk = -hk - sigma_full
      do i = 1, nmat
         gk(i, i) = gk(i, i) + z
      end do

      allocate (ipiv(nmat))
      call zgetrf(nmat, nmat, gk, nmat, ipiv, info)
      if (info /= 0) then
         gk = cmplx(huge(0.0_rp), 0.0_rp, rp)
         deallocate (ipiv)
         return
      end if
      lwork = nmat*nmat
      allocate (work(lwork))
      call zgetri(nmat, gk, nmat, ipiv, work, lwork, info)
      deallocate (ipiv, work)
   end subroutine dyson_kspace_inverse

end module dyson_kernel_mod
