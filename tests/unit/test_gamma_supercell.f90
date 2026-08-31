!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_gamma_supercell
!
!> @brief Gamma-only supercell identity + two-route DOS cross-validation (B2.5).
!> @details The ″beautiful self-contained identity″ of the B2 blueprint
!>          (validation item 2), with NO LMTO machinery: a periodic 1-band
!>          tight-binding ring of Nsc sites is described two equivalent ways, and
!>          the strict-Lehmann kernel (`lehmann_pair_block`, backend E) must give
!>          the same Green’s function both ways, to machine precision.
!>
!>          Route A -- PRIMITIVE cell (1 site), full BZ sampled at Nsc k-points
!>            k_n = 2*pi*n/Nsc, eps(k) = -2t cos k. G^A_{0,m}(z) is the k-sum
!>            (1/Nsc) sum_n e^{i k_n m} / (z - eps(k_n)) -- `lehmann_pair_block`
!>            with a 1-orbital eigensystem.
!>          Route L -- SUPERCELL (Nsc sites), a SINGLE k-point at Gamma. The
!>            supercell Hamiltonian H_sc is the Nsc x Nsc circulant hopping ring;
!>            its Gamma-only Lehmann sum sum_n psi_i psi_j^dagger / (z - E_n) over
!>            the H_sc eigenpairs is `lehmann_pair_block` with nk=1.
!>          Route C -- the real-space CLUSTER resolvent [z I - H_sc]^{-1} computed
!>            directly by `dyson_kspace_inverse` (backend D at Gamma).
!>
!>          Pins:
!>          1. Gamma-only supercell Lehmann  ==  cluster resolvent  (spectral
!>             theorem on H_sc: L route == C route), machine precision.
!>          2. Gamma-only supercell  ==  primitive-cell full-BZ Lehmann (the
!>             folding / supercell identity: L route == A route), machine
!>             precision. The eigenvalues of the circulant fold onto eps(k_n) and
!>             the eigenvectors are the Bloch plane waves, so the two k-samplings
!>             are literally the same sum.
!>          3. Two-route DOS: rho(E) = -1/pi Im Tr G(E + i eta) built from the
!>             Lehmann route (production kernel, site-summed on-site blocks) vs the
!>             Dyson route (trace of the direct inverse) agree to machine precision
!>             over an energy grid -- the self-contained k-space two-route DOS
!>             cross-validation. Also reports the integrated spectral weight, which
!>             tends to Nsc (C4 sum rule) as the grid widens / eta narrows.
!>
!>          Exits non-zero (error stop) on any failure so ctest registers a fail.
!------------------------------------------------------------------------------
program test_gamma_supercell
   use precision_mod, only: rp
   use math_mod, only: pi, two_pi, i_unit
   use lehmann_kernel_mod, only: lehmann_pair_block
   use dyson_kernel_mod, only: dyson_kspace_inverse
   implicit none

   integer, parameter :: nsc = 8          ! supercell / BZ-sampling size
   real(rp), parameter :: t = 1.0_rp
   real(rp), parameter :: tol_identity = 1.0e-12_rp
   real(rp), parameter :: tol_dos = 1.0e-12_rp
   logical :: failed

   failed = .false.

   call test_supercell_identity()
   call test_two_route_dos()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !> Build the Nsc x Nsc circulant nearest-neighbour hopping ring (periodic
   !> 1-band chain), H_sc(j,l) = -t for |j-l| = 1 (mod Nsc), else 0.
   subroutine build_supercell(hsc)
      complex(rp), intent(out) :: hsc(nsc, nsc)
      integer :: j, jp, jm

      hsc = (0.0_rp, 0.0_rp)
      do j = 1, nsc
         jp = merge(1, j + 1, j == nsc)
         jm = merge(nsc, j - 1, j == 1)
         hsc(j, jp) = cmplx(-t, 0.0_rp, rp)
         hsc(j, jm) = cmplx(-t, 0.0_rp, rp)
      end do
   end subroutine build_supercell

   !> Primitive-cell eigensystem: 1 site, Nsc k-points, eps(k) = -2t cos k.
   subroutine build_primitive(eigenvalues, eigenvectors, kfrac)
      real(rp), intent(out) :: eigenvalues(1, nsc)
      complex(rp), intent(out) :: eigenvectors(1, 1, nsc)
      real(rp), intent(out) :: kfrac(3, nsc)
      integer :: ik

      do ik = 1, nsc
         kfrac(:, ik) = 0.0_rp
         kfrac(1, ik) = real(ik - 1, rp)/real(nsc, rp)
         eigenvalues(1, ik) = -2.0_rp*t*cos(two_pi*kfrac(1, ik))
         eigenvectors(1, 1, ik) = (1.0_rp, 0.0_rp)
      end do
   end subroutine build_primitive

   !> Tests 1 & 2: Gamma-only supercell Lehmann == cluster resolvent == primitive
   !> full-BZ Lehmann, for the full first-row G_{0,m}, m = 0 .. Nsc-1.
   subroutine test_supercell_identity()
      complex(rp) :: hsc(nsc, nsc), vecs(nsc, nsc)
      real(rp) :: evals(nsc)
      real(rp) :: eig_sc(nsc, 1), kfrac_sc(3, 1)
      complex(rp) :: vec_sc(nsc, nsc, 1)
      real(rp) :: eig_pr(1, nsc), kfrac_pr(3, nsc)
      complex(rp) :: vec_pr(1, 1, nsc)
      complex(rp), allocatable :: z_contour(:)
      complex(rp) :: gL(1, 1, 3), gA(1, 1, 3), gC(nsc, nsc), sig0(nsc, nsc)
      real(rp) :: dr0(3), dr_pr(3)
      real(rp) :: err_LC, err_LA
      integer :: ie, ne, m

      call build_supercell(hsc)
      call hermitian_eig(hsc, nsc, evals, vecs)

      ! Gamma-only supercell eigensystem for lehmann_pair_block (nk = 1).
      kfrac_sc(:, 1) = 0.0_rp
      eig_sc(:, 1) = evals
      vec_sc(:, :, 1) = vecs

      call build_primitive(eig_pr, vec_pr, kfrac_pr)

      ne = 3
      allocate (z_contour(ne))
      z_contour(1) = cmplx(-1.3_rp, 0.30_rp, rp)
      z_contour(2) = cmplx(0.4_rp, 0.20_rp, rp)
      z_contour(3) = cmplx(1.1_rp, 0.45_rp, rp)
      sig0 = (0.0_rp, 0.0_rp)
      dr0 = 0.0_rp

      err_LC = 0.0_rp
      err_LA = 0.0_rp
      do m = 0, nsc - 1
         ! Route L: Gamma-only supercell Lehmann for the (0, m) block.
         call lehmann_pair_block(eig_sc, vec_sc, kfrac_sc, z_contour, dr0, 0, m, 1, gL)
         ! Route A: primitive full-BZ Lehmann, dR_{0m} = R_0 - R_m = -m.
         dr_pr = 0.0_rp
         dr_pr(1) = -real(m, rp)
         call lehmann_pair_block(eig_pr, vec_pr, kfrac_pr, z_contour, dr_pr, 0, 0, 1, gA)
         do ie = 1, ne
            call dyson_kspace_inverse(hsc, z_contour(ie), sig0, gC)
            err_LC = max(err_LC, abs(gL(1, 1, ie) - gC(1, m + 1)))
            err_LA = max(err_LA, abs(gL(1, 1, ie) - gA(1, 1, ie)))
         end do
      end do

      write (*, '(a,es12.4)') 'Test 1 (Gamma supercell == cluster)   max_err = ', err_LC
      write (*, '(a,es12.4)') 'Test 2 (Gamma supercell == primitive) max_err = ', err_LA
      if (err_LC > tol_identity) failed = .true.
      if (err_LA > tol_identity) failed = .true.

      deallocate (z_contour)
   end subroutine test_supercell_identity

   !> Test 3: two-route DOS rho(E) = -1/pi Im Tr G(E + i eta). Lehmann route
   !> (site-summed on-site blocks from the production kernel) vs Dyson route
   !> (trace of the direct inverse) -- machine-precision agreement + weight report.
   subroutine test_two_route_dos()
      complex(rp) :: hsc(nsc, nsc), vecs(nsc, nsc)
      real(rp) :: evals(nsc)
      real(rp) :: eig_sc(nsc, 1), kfrac_sc(3, 1)
      complex(rp) :: vec_sc(nsc, nsc, 1)
      complex(rp) :: z, g_on(1, 1, 1), gC(nsc, nsc), sig0(nsc, nsc)
      real(rp) :: dr0(3)
      real(rp), parameter :: eta = 0.20_rp
      real(rp), parameter :: emin = -3.5_rp, emax = 3.5_rp
      integer, parameter :: ne = 351
      real(rp) :: e, de, tr_le, tr_dy, dos_le, dos_dy
      real(rp) :: max_dos_err, weight_le, weight_dy
      real(rp) :: dos_le_prev, dos_dy_prev
      complex(rp) :: zc(1)
      integer :: ie, s

      call build_supercell(hsc)
      call hermitian_eig(hsc, nsc, evals, vecs)
      kfrac_sc(:, 1) = 0.0_rp
      eig_sc(:, 1) = evals
      vec_sc(:, :, 1) = vecs
      sig0 = (0.0_rp, 0.0_rp)
      dr0 = 0.0_rp

      de = (emax - emin)/real(ne - 1, rp)
      max_dos_err = 0.0_rp
      weight_le = 0.0_rp; weight_dy = 0.0_rp
      dos_le_prev = 0.0_rp; dos_dy_prev = 0.0_rp
      do ie = 1, ne
         e = emin + de*real(ie - 1, rp)
         z = cmplx(e, eta, rp)
         zc(1) = z

         ! Lehmann route: Tr G = sum over sites of the on-site (s,s) block.
         tr_le = 0.0_rp
         do s = 0, nsc - 1
            call lehmann_pair_block(eig_sc, vec_sc, kfrac_sc, zc, dr0, s, s, 1, g_on)
            tr_le = tr_le + aimag(g_on(1, 1, 1))
         end do
         dos_le = -tr_le/pi

         ! Dyson route: Tr of the direct inverse [z I - H_sc]^{-1}.
         call dyson_kspace_inverse(hsc, z, sig0, gC)
         tr_dy = 0.0_rp
         do s = 1, nsc
            tr_dy = tr_dy + aimag(gC(s, s))
         end do
         dos_dy = -tr_dy/pi

         max_dos_err = max(max_dos_err, abs(dos_le - dos_dy))
         if (ie > 1) then
            weight_le = weight_le + 0.5_rp*de*(dos_le + dos_le_prev)
            weight_dy = weight_dy + 0.5_rp*de*(dos_dy + dos_dy_prev)
         end if
         dos_le_prev = dos_le; dos_dy_prev = dos_dy
      end do

      write (*, '(a,es12.4)') 'Test 3 (two-route DOS agree)   max_err = ', max_dos_err
      write (*, '(a,f10.5,a,f10.5,a,i0,a)') 'Test 3 spectral weight  Lehmann = ', weight_le, &
         '  Dyson = ', weight_dy, '  (C4 sum rule -> ', nsc, ')'
      if (max_dos_err > tol_dos) failed = .true.
   end subroutine test_two_route_dos

   !> Hermitian eigensolve wrapper (zheev, jobz=’V’): columns of vecs are bands.
   subroutine hermitian_eig(a_in, n, w, vecs)
      integer, intent(in) :: n
      complex(rp), intent(in) :: a_in(n, n)
      real(rp), intent(out) :: w(n)
      complex(rp), intent(out) :: vecs(n, n)
      complex(rp) :: work_q(1)
      complex(rp), allocatable :: work(:)
      real(rp) :: rwork(max(1, 3*n - 2))
      integer :: info, lwork

      vecs = a_in
      call zheev('V', 'U', n, vecs, n, w, work_q, -1, rwork, info)
      lwork = max(1, nint(real(work_q(1), rp)))
      allocate (work(lwork))
      call zheev('V', 'U', n, vecs, n, w, work, lwork, rwork, info)
      deallocate (work)
      if (info /= 0) then
         write (*, '(a,i0)') 'hermitian_eig: zheev info = ', info
         failed = .true.
      end if
   end subroutine hermitian_eig

end program test_gamma_supercell
