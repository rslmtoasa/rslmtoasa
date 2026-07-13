!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_lehmann_chain
!
!> @brief Known-answer test for the strict-Lehmann kernel (backend E, B2.2).
!> @details Exercises `lehmann_pair_block` on a 1-band tight-binding chain
!>          H(k) = -2t cos k, whose lattice Green's function has a closed form,
!>          with NO LMTO machinery involved. Pins three things:
!>
!>          1. On-site block G_00(z) vs the closed form 1/sqrt(z^2 - 4t^2)
!>             (retarded branch) to 1e-10 -- the core Lehmann sum + broadening.
!>          2. Intersite block G_0m(z) vs the residue closed form
!>             x_in^|m| / (t (x_in - x_out)) -- pins the e^{i k.dR} phase (C3).
!>          3. Normalization z*G_00(z) -> 1 as |z| -> infinity -- pins the
!>             1/N_k factor (C4, pure-math form).
!>
!>          Exits with a non-zero status (error stop) on any failure so ctest
!>          registers a fail.
!------------------------------------------------------------------------------
program test_lehmann_chain
   use precision_mod, only: rp
   use math_mod, only: pi, i_unit
   use lehmann_kernel_mod, only: lehmann_pair_block, pauli_decompose_block
   implicit none

   integer, parameter :: nk = 512
   real(rp), parameter :: t = 1.0_rp
   real(rp), parameter :: tol_onsite = 1.0e-10_rp
   real(rp), parameter :: tol_inter = 1.0e-9_rp
   real(rp), parameter :: tol_norm = 1.0e-6_rp

   real(rp) :: eigenvalues(1, nk)
   complex(rp) :: eigenvectors(1, 1, nk)
   real(rp) :: kfrac(3, nk)
   complex(rp), allocatable :: z_contour(:)
   complex(rp), allocatable :: gblk(:, :, :)
   real(rp) :: dr(3), kang
   integer :: ik, ie, ne, m
   real(rp) :: err, max_err
   complex(rp) :: g_ref
   logical :: failed

   failed = .false.

   ! --- Build the 1-band chain eigensystem: H(k) = -2t cos k, psi = 1. -------
   do ik = 1, nk
      kfrac(:, ik) = 0.0_rp
      kfrac(1, ik) = real(ik - 1, rp)/real(nk, rp)
      kang = 2.0_rp*pi*kfrac(1, ik)
      eigenvalues(1, ik) = -2.0_rp*t*cos(kang)
      eigenvectors(1, 1, ik) = (1.0_rp, 0.0_rp)
   end do

   ! --- Retarded contour: a spread of real energies with fixed broadening. ---
   ne = 9
   allocate (z_contour(ne), gblk(1, 1, ne))
   do ie = 1, ne
      z_contour(ie) = cmplx(-2.5_rp + 0.5_rp*real(ie - 1, rp), 0.3_rp, rp)
   end do

   ! --- Test 1: on-site block vs closed form. --------------------------------
   dr = 0.0_rp
   call lehmann_pair_block(eigenvalues, eigenvectors, kfrac, z_contour, dr, 0, 0, 1, gblk)
   max_err = 0.0_rp
   do ie = 1, ne
      g_ref = chain_gf_closed(z_contour(ie), 0)
      err = abs(gblk(1, 1, ie) - g_ref)
      max_err = max(max_err, err)
   end do
   write (*, '(a,es12.4)') 'Test 1 (on-site closed form)   max_err = ', max_err
   if (max_err > tol_onsite) failed = .true.

   ! --- Test 2: intersite block (m = 2) vs residue closed form (phase). ------
   m = 2
   dr = 0.0_rp
   dr(1) = -real(m, rp)   ! dR_ij = R_i - R_j = 0 - m
   call lehmann_pair_block(eigenvalues, eigenvectors, kfrac, z_contour, dr, 0, 0, 1, gblk)
   max_err = 0.0_rp
   do ie = 1, ne
      g_ref = chain_gf_closed(z_contour(ie), m)
      err = abs(gblk(1, 1, ie) - g_ref)
      max_err = max(max_err, err)
   end do
   write (*, '(a,es12.4)') 'Test 2 (intersite m=2 phase)   max_err = ', max_err
   if (max_err > tol_inter) failed = .true.

   ! --- Test 3: normalization z*G_00 -> 1 for large |z|. ---------------------
   deallocate (z_contour, gblk)
   allocate (z_contour(1), gblk(1, 1, 1))
   z_contour(1) = cmplx(0.0_rp, 1.0e6_rp, rp)
   dr = 0.0_rp
   call lehmann_pair_block(eigenvalues, eigenvectors, kfrac, z_contour, dr, 0, 0, 1, gblk)
   err = abs(z_contour(1)*gblk(1, 1, 1) - (1.0_rp, 0.0_rp))
   write (*, '(a,es12.4)') 'Test 3 (normalization z*G->1)  err     = ', err
   if (err > tol_norm) failed = .true.

   ! --- Test 4: Pauli (charge,x,y,z) decomposition transcription guard. -------
   ! Build a one-orbital spin block from known (c,x,y,z) via the standard
   ! assembly G = [[c+z, x-i y],[x+i y, c-z]] (up=orb1, dn=orb1+spin_off) and
   ! check pauli_decompose_block recovers each component -- pins the signs and
   ! the i factor of the torque algebra shared by gij and the eta ladder.
   block
      complex(rp) :: blk(2, 2, 2), gnmag(1, 1, 2), gz(1, 1, 2), gy(1, 1, 2), gx(1, 1, 2)
      complex(rp) :: c(2), sx(2), sy(2), sz(2)
      integer :: s

      c = [cmplx(0.7_rp, -0.2_rp, rp), cmplx(-1.1_rp, 0.4_rp, rp)]
      sx = [cmplx(0.3_rp, 0.5_rp, rp), cmplx(1.2_rp, -0.9_rp, rp)]
      sy = [cmplx(-0.6_rp, 0.1_rp, rp), cmplx(0.4_rp, 0.7_rp, rp)]
      sz = [cmplx(0.9_rp, -0.3_rp, rp), cmplx(-0.5_rp, 0.2_rp, rp)]
      do s = 1, 2
         blk(1, 1, s) = c(s) + sz(s)
         blk(2, 2, s) = c(s) - sz(s)
         blk(1, 2, s) = sx(s) - i_unit*sy(s)
         blk(2, 1, s) = sx(s) + i_unit*sy(s)
      end do
      call pauli_decompose_block(blk, 1, gnmag, gz, gy, gx)
      max_err = 0.0_rp
      do s = 1, 2
         max_err = max(max_err, abs(gnmag(1, 1, s) - c(s)))
         max_err = max(max_err, abs(gx(1, 1, s) - sx(s)))
         max_err = max(max_err, abs(gy(1, 1, s) - sy(s)))
         max_err = max(max_err, abs(gz(1, 1, s) - sz(s)))
      end do
      write (*, '(a,es12.4)') 'Test 4 (Pauli decomposition)   max_err = ', max_err
      if (max_err > tol_onsite) failed = .true.
   end block

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !> Closed-form retarded lattice Green's function of the 1-band chain,
   !> G_0m(z) = x_in^|m| / (t (x_in - x_out)), where x_in, x_out are the roots
   !> of x^2 + (z/t) x + 1 = 0 with |x_in| < 1 (residue inside the unit circle).
   !> For m = 0 this reduces to 1/sqrt(z^2 - 4t^2) on the retarded branch.
   function chain_gf_closed(z, m) result(g)
      complex(rp), intent(in) :: z
      integer, intent(in) :: m
      complex(rp) :: g
      complex(rp) :: disc, root, xa, xb, x_in, x_out

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

end program test_lehmann_chain
