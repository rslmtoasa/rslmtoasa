!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: lehmann_kernel
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Pure-math kernel for the strict Lehmann (backend E) k-space Green's function
!> (milestone B2, `reciprocal_green`). Given the Hermitian eigenpairs of the
!> k-space Hamiltonian H(k), it accumulates the intersite resolvent block
!>
!>   G_ij(z) = (1/N_k) sum_{k,n} e^{i k.dR_ij} psi_{i,n}(k) psi^dagger_{j,n}(k)
!>                                              / (z - eps_{nk})
!>
!> where psi_{i,n}(k) is the sub-block of eigenvector n restricted to the
!> orbitals of unit-cell site i, and dR_ij = R_i - R_j is the real-space bond
!> vector between the two atoms of the pair (in fractional/direct coordinates).
!>
!> This routine deliberately carries NO dependence on the `reciprocal`,
!> `lattice`, or `green` objects: it is the isolated numerical core so the
!> 1-band-chain known-answer test (`tests/unit/test_lehmann_chain.f90`) can
!> exercise it against a closed form with no LMTO machinery. `reciprocal_green`
!> wires the LMTO eigenpairs, the pair->site map, and the bond vectors around it.
!>
!> Sign/phase convention (gate G-B2-1, C3). H(k) is assembled in the "atomic"
!> Bloch gauge of `reciprocal_fourier`: H_ij(k) = sum_R h_{i0,jR} e^{i k.R} with
!> R the true bond vector r_j - r_i (basis offset included), phase e^{+i 2pi
!> k_frac.R_frac} through `ham_vec_type_direct`. The inverse transform of
!> G_ij(k) = sum_n c_{i,n} c*_{j,n} / (z - eps) back to the real-space pair is
!>   G_{i0,jR}(z) = (1/N_k) sum_k e^{-i k.R} G_ij(k),   R = r_j - r_i,
!> i.e. e^{+i k.(r_i - r_j)}. Hence the caller passes dR_ij = R_i - R_j in
!> fractional coordinates and the kernel applies e^{+i 2pi k_frac.dR_ij}. This
!> reuses the B1 neighbor-vector table and sign verbatim (C3). For the on-site
!> block (i = j) dR_ij = 0 and the phase is unity.
!------------------------------------------------------------------------------
module lehmann_kernel_mod
   use precision_mod, only: rp
   use math_mod, only: two_pi
   implicit none

   private
   public :: lehmann_pair_block

contains

   !> @brief Accumulate one strict-Lehmann intersite Green's-function block.
   !> @details Streams over all k-points and bands, amortizing the eigenpairs
   !>          over every energy on the contour (design B2 §1.4). The result is
   !>          normalized by 1/N_k. See the module header for the phase
   !>          convention. `eigenvectors(:, n, ik)` is band n at k-point ik with
   !>          the site-i orbitals occupying rows `ioff+1 .. ioff+nblk`.
   !> @param[in]  eigenvalues  Real eigenvalues, shape (nmat, nk).
   !> @param[in]  eigenvectors Column eigenvectors, shape (nmat, nmat, nk).
   !> @param[in]  kfrac        k-points in fractional coordinates, shape (3, nk).
   !> @param[in]  z_contour    Complex energy contour, shape (ne).
   !> @param[in]  dr_direct    Bond vector R_i - R_j in fractional coordinates.
   !> @param[in]  ioff         Zero-based row offset of site i's block in nmat.
   !> @param[in]  joff         Zero-based row offset of site j's block in nmat.
   !> @param[in]  nblk         Block size (nb = 2*norb).
   !> @param[out] gblk         Accumulated block, shape (nblk, nblk, ne).
   subroutine lehmann_pair_block(eigenvalues, eigenvectors, kfrac, z_contour, &
                                 dr_direct, ioff, joff, nblk, gblk)
      real(rp), intent(in) :: eigenvalues(:, :)
      complex(rp), intent(in) :: eigenvectors(:, :, :)
      real(rp), intent(in) :: kfrac(:, :)
      complex(rp), intent(in) :: z_contour(:)
      real(rp), intent(in) :: dr_direct(3)
      integer, intent(in) :: ioff, joff, nblk
      complex(rp), intent(out) :: gblk(:, :, :)

      integer :: nmat, nk, ne, ik, n, ie, a, b
      real(rp) :: kdotr
      complex(rp) :: phase, w, cjb
      complex(rp) :: m_ab(nblk, nblk)

      nmat = size(eigenvalues, 1)
      nk = size(eigenvalues, 2)
      ne = size(z_contour)

      gblk = (0.0_rp, 0.0_rp)

      do ik = 1, nk
         kdotr = two_pi*dot_product(kfrac(:, ik), dr_direct)
         phase = cmplx(cos(kdotr), sin(kdotr), rp)
         do n = 1, nmat
            ! Outer product psi_i psi_j^dagger: m_ab = c_{i,a,n} * conjg(c_{j,b,n}).
            do b = 1, nblk
               cjb = conjg(eigenvectors(joff + b, n, ik))
               do a = 1, nblk
                  m_ab(a, b) = eigenvectors(ioff + a, n, ik)*cjb
               end do
            end do
            ! Amortize the eigenpair over the whole energy contour.
            do ie = 1, ne
               w = phase/(z_contour(ie) - eigenvalues(n, ik))
               gblk(:, :, ie) = gblk(:, :, ie) + w*m_ab
            end do
         end do
      end do

      gblk = gblk/real(nk, rp)
   end subroutine lehmann_pair_block

end module lehmann_kernel_mod
