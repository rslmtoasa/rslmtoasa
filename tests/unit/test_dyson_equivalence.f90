!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_dyson_equivalence
!
!> @brief Permanent CI invariant: backend D (Dyson) with Sigma = 0 == backend E
!>        (Lehmann), plus the Sigma sign pin (B2.4).
!> @details Exercises `dyson_kspace_inverse` against `lehmann_pair_block` on
!>          k-space Hamiltonians with NO LMTO machinery, reproducing the exact
!>          per-(k,z) inverse -> sub-block -> e^{i k.dR} phase -> 1/N_k
!>          accumulation `fill_green_dyson`/`fill_green_lehmann` perform:
!>
!>          1. 1-band chain H(k) = -2t cos k: Dyson route vs the closed-form
!>             lattice Green's function (on-site and intersite m=2) -- the same
!>             fixture test_lehmann_chain uses for backend E.
!>          2. Small multiband H(k) (2 sites x 2 orbitals, nmat=4): the intersite
!>             (ioff=0, joff=2) block with a nonzero bond phase, Dyson route vs
!>             Lehmann route (eigenpairs from zheev), to 1e-9 elementwise -- the
!>             D == E invariant on a nontrivial H(k) (the "small bcc-Fe H(k)"
!>             fixture of the B2.4 spec).
!>          3. Sigma sign pin: a constant real shift Sigma = s0*I must give the
!>             Dyson block at z equal to the Lehmann block at z - s0 (A = zI - H
!>             - Sigma). A sign error here would silently corrupt every provider.
!>
!>          Exits non-zero (error stop) on any failure so ctest registers a fail.
!------------------------------------------------------------------------------
program test_dyson_equivalence
   use precision_mod, only: rp
   use math_mod, only: pi, two_pi, i_unit
   use lehmann_kernel_mod, only: lehmann_pair_block
   use dyson_kernel_mod, only: dyson_kspace_inverse
   use, intrinsic :: ieee_exceptions, only: ieee_divide_by_zero, ieee_get_halting_mode, ieee_set_flag, ieee_set_halting_mode
   implicit none

   integer, parameter :: nk = 256
   real(rp), parameter :: t = 1.0_rp
   real(rp), parameter :: tol_chain = 1.0e-9_rp
   real(rp), parameter :: tol_equiv = 1.0e-9_rp
   logical :: failed

   failed = .false.

   call test_chain()
   call test_multiband()
   call test_sigma_sign()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !> Test 1: 1-band chain, Dyson route vs closed form (on-site + intersite).
   subroutine test_chain()
      real(rp) :: eigenvalues(1, nk), kfrac(3, nk), kang
      complex(rp) :: sig0(1, 1)
      complex(rp), allocatable :: z_contour(:)
      integer :: ik, ie, ne

      do ik = 1, nk
         kfrac(:, ik) = 0.0_rp
         kfrac(1, ik) = real(ik - 1, rp)/real(nk, rp)
         kang = two_pi*kfrac(1, ik)
         eigenvalues(1, ik) = -2.0_rp*t*cos(kang)
      end do

      ne = 9
      allocate (z_contour(ne))
      do ie = 1, ne
         z_contour(ie) = cmplx(-2.5_rp + 0.5_rp*real(ie - 1, rp), 0.3_rp, rp)
      end do
      sig0 = (0.0_rp, 0.0_rp)

      ! On-site (m=0) and intersite (m=2): Dyson-route accumulation vs closed form.
      call check_chain_block(0, z_contour, ne, eigenvalues, kfrac, sig0, 'Test 1a (chain on-site   Dyson)')
      call check_chain_block(2, z_contour, ne, eigenvalues, kfrac, sig0, 'Test 1b (chain intersite Dyson)')

      deallocate (z_contour)
   end subroutine test_chain

   !> Accumulate the Dyson-route chain block for bond m and compare to closed form.
   subroutine check_chain_block(m, z_contour, ne, eigenvalues, kfrac, sig0, label)
      integer, intent(in) :: m, ne
      complex(rp), intent(in) :: z_contour(:), sig0(1, 1)
      real(rp), intent(in) :: eigenvalues(1, nk), kfrac(3, nk)
      character(*), intent(in) :: label
      complex(rp) :: hk(1, 1), gk(1, 1), acc, phase, g_ref
      real(rp) :: dr(3), kdotr, err, max_err
      integer :: ik, ie

      dr = 0.0_rp
      dr(1) = -real(m, rp)   ! dR_ij = R_i - R_j = 0 - m
      max_err = 0.0_rp
      do ie = 1, ne
         acc = (0.0_rp, 0.0_rp)
         do ik = 1, nk
            hk(1, 1) = cmplx(eigenvalues(1, ik), 0.0_rp, rp)
            call dyson_kspace_inverse(hk, z_contour(ie), sig0, gk)
            kdotr = two_pi*dot_product(kfrac(:, ik), dr)
            phase = cmplx(cos(kdotr), sin(kdotr), rp)
            acc = acc + phase*gk(1, 1)
         end do
         acc = acc/real(nk, rp)
         g_ref = chain_gf_closed(z_contour(ie), m)
         err = abs(acc - g_ref)
         max_err = max(max_err, err)
      end do
      write (*, '(a,es12.4)') label//'   max_err = ', max_err
      if (max_err > tol_chain) failed = .true.
   end subroutine check_chain_block

   !> Test 2: multiband H(k), Dyson intersite block vs Lehmann (eigenpair) block.
   subroutine test_multiband()
      integer, parameter :: n = 4, nblk = 2, ioff = 0, joff = 2
      real(rp) :: kfrac(3, nk), theta, dr(3)
      complex(rp) :: hk_all(n, n, nk), d0(n, n), m0(n, n)
      real(rp) :: eigenvalues(n, nk)
      complex(rp) :: eigenvectors(n, n, nk)
      complex(rp) :: hk(n, n), gk(n, n), sig0(n, n)
      complex(rp), allocatable :: z_contour(:)
      complex(rp) :: gd(nblk, nblk), phase
      complex(rp), allocatable :: gle(:, :, :)
      real(rp) :: kdotr, max_err, err
      integer :: ik, ie, ne, a, b

      ! Fixed Hermitian building blocks: D0 diagonal (site energies), M0 hopping.
      d0 = (0.0_rp, 0.0_rp)
      d0(1, 1) = (-0.6_rp, 0.0_rp); d0(2, 2) = (0.2_rp, 0.0_rp)
      d0(3, 3) = (0.5_rp, 0.0_rp); d0(4, 4) = (-0.1_rp, 0.0_rp)
      m0 = (0.0_rp, 0.0_rp)
      m0(1, 3) = cmplx(0.30_rp, 0.10_rp, rp); m0(1, 4) = cmplx(-0.20_rp, 0.05_rp, rp)
      m0(2, 3) = cmplx(0.15_rp, -0.25_rp, rp); m0(2, 4) = cmplx(0.40_rp, 0.00_rp, rp)
      m0(1, 2) = cmplx(0.10_rp, 0.20_rp, rp); m0(3, 4) = cmplx(-0.30_rp, 0.10_rp, rp)

      do ik = 1, nk
         kfrac(:, ik) = 0.0_rp
         kfrac(1, ik) = real(ik - 1, rp)/real(nk, rp)
         theta = two_pi*kfrac(1, ik)
         ! H(k) = D0 + M0 e^{i theta} + M0^dagger e^{-i theta}  (Hermitian).
         hk_all(:, :, ik) = d0 + m0*exp(i_unit*theta) + &
                            conjg(transpose(m0))*exp(-i_unit*theta)
      end do

      ! Eigenpairs for the Lehmann route (one zheev per k).
      do ik = 1, nk
         hk = hk_all(:, :, ik)
         call hermitian_eig(hk, n, eigenvalues(:, ik), eigenvectors(:, :, ik))
      end do

      ne = 7
      allocate (z_contour(ne), gle(nblk, nblk, ne))
      do ie = 1, ne
         z_contour(ie) = cmplx(-1.2_rp + 0.4_rp*real(ie - 1, rp), 0.25_rp, rp)
      end do
      sig0 = (0.0_rp, 0.0_rp)
      dr = 0.0_rp; dr(1) = 0.37_rp

      ! Lehmann route: intersite (ioff,joff) block over the eigenpairs.
      call lehmann_pair_block(eigenvalues, eigenvectors, kfrac, z_contour, &
                              dr, ioff, joff, nblk, gle)

      ! Dyson route: same block by inverting [zI - H(k)] per (k,z).
      max_err = 0.0_rp
      do ie = 1, ne
         gd = (0.0_rp, 0.0_rp)
         do ik = 1, nk
            call dyson_kspace_inverse(hk_all(:, :, ik), z_contour(ie), sig0, gk)
            kdotr = two_pi*dot_product(kfrac(:, ik), dr)
            phase = cmplx(cos(kdotr), sin(kdotr), rp)
            gd = gd + phase*gk(ioff + 1:ioff + nblk, joff + 1:joff + nblk)
         end do
         gd = gd/real(nk, rp)
         do b = 1, nblk
            do a = 1, nblk
               err = abs(gd(a, b) - gle(a, b, ie))
               max_err = max(max_err, err)
            end do
         end do
      end do
      write (*, '(a,es12.4)') 'Test 2 (multiband D==E block)  max_err = ', max_err
      if (max_err > tol_equiv) failed = .true.

      deallocate (z_contour, gle)
   end subroutine test_multiband

   !> Test 3: Sigma = s0*I shifts the Dyson block by exactly s0 in z (sign pin).
   subroutine test_sigma_sign()
      integer, parameter :: n = 4, nblk = 2, ioff = 0, joff = 2
      real(rp) :: kfrac(3, nk), theta, dr(3)
      complex(rp) :: hk_all(n, n, nk), d0(n, n), m0(n, n)
      real(rp) :: eigenvalues(n, nk)
      complex(rp) :: eigenvectors(n, n, nk)
      complex(rp) :: gk(n, n), sig_full(n, n)
      complex(rp), allocatable :: z_contour(:), z_shift(:), gle(:, :, :)
      complex(rp) :: gd(nblk, nblk), phase, s0
      real(rp) :: kdotr, max_err, err
      integer :: ik, ie, ne, a, b

      d0 = (0.0_rp, 0.0_rp)
      d0(1, 1) = (-0.6_rp, 0.0_rp); d0(2, 2) = (0.2_rp, 0.0_rp)
      d0(3, 3) = (0.5_rp, 0.0_rp); d0(4, 4) = (-0.1_rp, 0.0_rp)
      m0 = (0.0_rp, 0.0_rp)
      m0(1, 3) = cmplx(0.30_rp, 0.10_rp, rp); m0(2, 4) = cmplx(0.40_rp, 0.00_rp, rp)
      m0(1, 2) = cmplx(0.10_rp, 0.20_rp, rp); m0(3, 4) = cmplx(-0.30_rp, 0.10_rp, rp)

      do ik = 1, nk
         kfrac(:, ik) = 0.0_rp
         kfrac(1, ik) = real(ik - 1, rp)/real(nk, rp)
         theta = two_pi*kfrac(1, ik)
         hk_all(:, :, ik) = d0 + m0*exp(i_unit*theta) + &
                            conjg(transpose(m0))*exp(-i_unit*theta)
         call hermitian_eig(hk_all(:, :, ik), n, eigenvalues(:, ik), eigenvectors(:, :, ik))
      end do

      s0 = cmplx(0.35_rp, 0.0_rp, rp)   ! real constant self-energy shift
      sig_full = (0.0_rp, 0.0_rp)
      do a = 1, n
         sig_full(a, a) = s0
      end do

      ne = 5
      allocate (z_contour(ne), z_shift(ne), gle(nblk, nblk, ne))
      do ie = 1, ne
         z_contour(ie) = cmplx(-1.0_rp + 0.4_rp*real(ie - 1, rp), 0.25_rp, rp)
         z_shift(ie) = z_contour(ie) - s0   ! Lehmann reference at z - Sigma
      end do
      dr = 0.0_rp; dr(1) = 0.37_rp

      call lehmann_pair_block(eigenvalues, eigenvectors, kfrac, z_shift, &
                              dr, ioff, joff, nblk, gle)

      max_err = 0.0_rp
      do ie = 1, ne
         gd = (0.0_rp, 0.0_rp)
         do ik = 1, nk
            call dyson_kspace_inverse(hk_all(:, :, ik), z_contour(ie), sig_full, gk)
            kdotr = two_pi*dot_product(kfrac(:, ik), dr)
            phase = cmplx(cos(kdotr), sin(kdotr), rp)
            gd = gd + phase*gk(ioff + 1:ioff + nblk, joff + 1:joff + nblk)
         end do
         gd = gd/real(nk, rp)
         do b = 1, nblk
            do a = 1, nblk
               err = abs(gd(a, b) - gle(a, b, ie))
               max_err = max(max_err, err)
            end do
         end do
      end do
      write (*, '(a,es12.4)') 'Test 3 (Sigma sign shift)      max_err = ', max_err
      if (max_err > tol_equiv) failed = .true.

      deallocate (z_contour, z_shift, gle)
   end subroutine test_sigma_sign

   !> Hermitian eigensolve wrapper (zheev, jobz='V'): columns of vecs are bands.
   subroutine hermitian_eig(a_in, n, w, vecs)
      integer, intent(in) :: n
      complex(rp), intent(in) :: a_in(n, n)
      real(rp), intent(out) :: w(n)
      complex(rp), intent(out) :: vecs(n, n)
      complex(rp) :: work_q(1)
      complex(rp), allocatable :: work(:)
      real(rp) :: rwork(max(1, 3*n - 2))
      integer :: info, lwork
      logical :: halt_divide_by_zero

      vecs = a_in
      call zheev('V', 'U', n, vecs, n, w, work_q, -1, rwork, info)
      lwork = max(1, nint(real(work_q(1), rp)))
      allocate (work(lwork))

      ! oneMKL's zheev implementation deliberately evaluates an intermediate
      ! divide-by-zero in its scaling path.  Keep the test's FPE diagnostics
      ! enabled everywhere else, but mask that external-library exception only
      ! for the LAPACK call and clear its status flag before restoring the mode.
      call ieee_get_halting_mode(ieee_divide_by_zero, halt_divide_by_zero)
      call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
      call zheev('V', 'U', n, vecs, n, w, work, lwork, rwork, info)
      call ieee_set_flag(ieee_divide_by_zero, .false.)
      call ieee_set_halting_mode(ieee_divide_by_zero, halt_divide_by_zero)
      deallocate (work)
      if (info /= 0) then
         write (*, '(a,i0)') 'hermitian_eig: zheev info = ', info
         failed = .true.
      end if
   end subroutine hermitian_eig

   !> Closed-form retarded lattice Green's function of the 1-band chain,
   !> G_0m(z) = x_in^|m| / (t (x_in - x_out)); for m=0, 1/sqrt(z^2 - 4t^2).
   function chain_gf_closed(z, m) result(g)
      complex(rp), intent(in) :: z
      integer, intent(in) :: m
      complex(rp) :: g, disc, xa, xb, x_in, x_out

      disc = sqrt((z/t)**2 - 4.0_rp)
      xa = 0.5_rp*(-(z/t) + disc)
      xb = 0.5_rp*(-(z/t) - disc)
      if (abs(xa) < abs(xb)) then
         x_in = xa; x_out = xb
      else
         x_in = xb; x_out = xa
      end if
      g = (x_in**abs(m))/(t*(x_in - x_out))
   end function chain_gf_closed

end program test_dyson_equivalence
