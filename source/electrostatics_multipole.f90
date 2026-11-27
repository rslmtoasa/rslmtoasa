!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief
!> Higher-order (l=1 dipole) multipole moments for surface/interface
!> electrostatics — RS-LMTO-ASA audit item B6.
!>
!> The ASA Madelung problem is extended from monopoles (q_R) to include the
!> l=1 (z) dipole moment Q_{10} of each sphere’s charge. The generalized
!> Madelung matrices (dipole-monopole dsz, dipole-dipole dzz) are already
!> built by charge%madl2d following Skriver-Rosengaard, Phys. Rev. B 43, 9538;
!> this module supplies the missing charge-side moment Q_{10}, which
!> charge%surfpot then couples into the local potential shift via dsz.
!>
!> This is the *new* code that lives outside the legacy self.f90 fence
!> (CLAUDE.md rule 5): it reads the on-site Green function (bands%green%g0)
!> and the exported radial partial-wave amplitudes (self%phi_amp, passed in as
!> an argument to avoid a self_mod <-> this-module use cycle), and writes only
!> the derived scalar potential%q10.
!>
!> Q_{10} vanishes in the bulk limit (no broken symmetry -> no s-p / p-d
!> cross-orbital density), so with the feature disabled or in bulk the
!> electrostatics reduce bit-for-bit to the existing monopole result.
!---------------------------------------------------------------------------
module electrostatics_multipole_mod
   use mpi_mod, only: start_atom, end_atom, g2l_map, rank
#ifdef USE_MPI
   use mpi
#endif
   use bands_mod, only: bands
   use precision_mod, only: rp
   use math_mod, only: pi, gaunt
   use basis_mod, only: spin_off, lmax_basis
   use logger_mod, only: g_logger
   use string_mod, only: fmt
   implicit none

   private


   !> Orbital indices in the SPHERICAL-harmonic basis that g0 lives in. The
   !> Hamiltonian is transformed cartesian -> spherical via hcpx(...,’cart2sph’)
   !> on every on-site block (hamiltonian_build.f90 build_obarm/build_enim/hhmag),
   !> so g0 is ordered by m, descending within each shell, functions (-1)^m Y_lm (math.f90 hcpx, n=9):
   !>   (00)(11)(10)(1-1)(22)(21)(20)(2-1)(2-2)
   !> => 1=s, 2=p(m=+1), 3=p(m=0), 4=p(m=-1),
   !>    5=d(m=+2), 6=d(m=+1), 7=d(m=0), 8=d(m=-1), 9=d(m=-2).
   !> The l=1,m=0 (z) dipole couples s(1)<->p(m=0)=3 and p(m=0)=3<->d(m=0)=7.
   integer, parameter :: idx_s   = 1
   integer, parameter :: idx_pz  = 3
   integer, parameter :: idx_dz2 = 7

   !> Angular factors for the physical z-dipole moment Q10 = int n(r) z d^3r,
   !> with z = r*sqrt(4pi/3)*Y10 and D_{L1,L2} symmetrized over both orderings:
   !>   c = 2 * sqrt(4pi/3) * <Y_{l1,0}|Y_{1,0}|Y_{l2,0}>
   !> Gaunt (verified vs Wigner-3j): <Y00 Y10 Y10>=1/(2 sqrt pi),
   !> <Y10 Y10 Y20>=1/(sqrt5 sqrt pi), giving c_sp = 2/sqrt(3), c_pd = 4/sqrt(15).
   !> The factor 2 is folded into the symmetrized d_spz/d_pzd below, so the
   !> per-ordering coefficients used in the sum are:
   real(rp), parameter :: c_sp = 1.0_rp/sqrt(3.0_rp)
   real(rp), parameter :: c_pd = 2.0_rp/sqrt(15.0_rp)

   public :: compute_dipole_moments, compute_multipoles, real_gaunt_table, add_multipole_block

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Compute the l=1 (z) dipole charge moment Q_{10} for every recursion atom
   !> and store it in potential%q10.
   !>
   !> For each rank-local atom:
   !>   1. Extract the spin-resolved on-site cross-l density-matrix elements
   !>      D_{s,pz} and D_{pz,dz2} by integrating -Im G(L1,L2)/pi to the Fermi
   !>      level (same simpson_m machinery as bands%calculate_moments).
   !>   2. Build the radial dipole matrix elements R_sp, R_pd from the exported
   !>      partial-wave amplitudes u_l(r) = phi_amp over the WS mesh.
   !>   3. Assemble Q_{10} = sum_sigma [ C_sp D_{s,pz} R_sp + C_pd D_{pz,dz2} R_pd ].
   !> The per-atom results are reduced across MPI ranks so rank 0 (which runs
   !> charge%surfpot) sees the full set.
   !>
   !> @param[inout] bnd    Bands object (Green function, energy mesh, atoms).
   !> @param[in]    phi_amp Exported signed large-component radial amplitudes
   !>                       u_l(r)=r*phi_l(r), indexed (mesh, l+1, spin).
   !> @param[in]    mix    Linear mixing factor for the dipole moment across SCF
   !>                      steps: q10 <- mix*q10_new + (1-mix)*q10_old. The dipole
   !>                      feeds back into the potential, so an undamped update
   !>                      (mix=1) destabilizes the SCF; use a small value.
   !---------------------------------------------------------------------------
   subroutine compute_dipole_moments(bnd, phi_amp, mix)
      class(bands), intent(inout) :: bnd
      real(rp), dimension(:, :, :), intent(in) :: phi_amp
      real(rp), intent(in) :: mix
      integer :: na_glob, na_loc, plusbulk, isp, soff, nrec
      real(rp) :: d_spz, d_pzd, r_sp, r_pd, q10, q10_old
      real(rp), dimension(:), allocatable :: q10_all
      real(rp), dimension(size(bnd%en%ene)) :: y
#ifdef USE_MPI
      integer :: ierr
#endif

      nrec = bnd%lattice%nrec
      allocate (q10_all(nrec))
      q10_all(:) = 0.0_rp

      do na_glob = start_atom, end_atom
         na_loc = g2l_map(na_glob)
         plusbulk = bnd%lattice%nbulk + na_glob

         ! Radial dipole matrix elements from the current SCF partial waves.
         call radial_dipole_elements(bnd, plusbulk, phi_amp, r_sp, r_pd)

         q10 = 0.0_rp
         do isp = 1, 2
            soff = (isp - 1)*spin_off

            ! Symmetrized cross-l density-matrix elements (factor 2 accounts for
            ! the (L1,L2) and (L2,L1) orderings of the Hermitian on-site block).
            y(:) = -aimag(bnd%green%g0(idx_s + soff, idx_pz + soff, 1:size(bnd%en%ene), na_loc) + &
                          bnd%green%g0(idx_pz + soff, idx_s + soff, 1:size(bnd%en%ene), na_loc))/pi
            d_spz = bnd%integrate(y, 0)

            y(:) = -aimag(bnd%green%g0(idx_pz + soff, idx_dz2 + soff, 1:size(bnd%en%ene), na_loc) + &
                          bnd%green%g0(idx_dz2 + soff, idx_pz + soff, 1:size(bnd%en%ene), na_loc))/pi
            d_pzd = bnd%integrate(y, 0)

            q10 = q10 + c_sp*d_spz*r_sp + c_pd*d_pzd*r_pd
         end do

         q10_all(na_glob) = q10
      end do

#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, q10_all, nrec, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      do na_glob = 1, nrec
         ! Linear mixing against the previous SCF step to damp the dipole feedback.
         q10_old = bnd%symbolic_atom(bnd%lattice%nbulk + na_glob)%potential%q10
         q10_all(na_glob) = mix*q10_all(na_glob) + (1.0_rp - mix)*q10_old
         bnd%symbolic_atom(bnd%lattice%nbulk + na_glob)%potential%q10 = q10_all(na_glob)
      end do

      deallocate (q10_all)
   end subroutine compute_dipole_moments

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Radial dipole matrix elements over the Wigner-Seitz sphere:
   !>   R_sp = int u_s(r) u_pz(r) r dr,   R_pd = int u_pz(r) u_dz2(r) r dr,
   !> where u_l(r) = r*phi_l(r) is the exported large-component amplitude
   !> phi_amp(:, l+1, spin). The mesh is reconstructed exactly as in self.f90
   !> (logarithmic rofi, DRDI = a*(rofi+b)) and integrated with the code’s
   !> standard Simpson weights. The spin-averaged amplitude is used for the
   !> radial factor (spin dependence enters through the density matrix).
   !---------------------------------------------------------------------------
   subroutine radial_dipole_elements(bnd, plusbulk, phi_amp, r_sp, r_pd)
      class(bands), intent(inout) :: bnd
      integer, intent(in) :: plusbulk
      real(rp), dimension(:, :, :), intent(in) :: phi_amp
      real(rp), intent(out) :: r_sp, r_pd
      integer :: nr, ir
      real(rp) :: ws_r
      real(rp) :: phi_s, phi_p, phi_d
      real(rp), allocatable :: r(:), wd(:)

      r_sp = 0.0_rp
      r_pd = 0.0_rp

      nr = bnd%symbolic_atom(plusbulk)%mesh_grid_size()
      if (nr < 3) return
      ws_r = bnd%symbolic_atom(plusbulk)%potential%ws_r
      call radial_mesh(bnd, plusbulk, r, wd)

      do ir = 2, nr
         ! Spin-averaged large-component amplitudes u_l(r); channel index l+1 (s=1, p=2, d=3).
         phi_s = 0.5_rp*(phi_amp(ir, 1, 1) + phi_amp(ir, 1, 2))
         phi_p = 0.5_rp*(phi_amp(ir, 2, 1) + phi_amp(ir, 2, 2))
         phi_d = 0.5_rp*(phi_amp(ir, 3, 1) + phi_amp(ir, 3, 2))

         ! ASA+M multipole normalization: the l=1 moment uses (r/S)^l with l=1,
         ! i.e. the dipole length is measured in Wigner-Seitz radii (Vitos/Skriver;
         ! matches dsz carrying FACP = S*sqrt(3)). Hence the r/ws_r weight.
         r_sp = r_sp + wd(ir)*phi_s*phi_p*(r(ir)/ws_r)
         r_pd = r_pd + wd(ir)*phi_p*phi_d*(r(ir)/ws_r)
      end do
   end subroutine radial_dipole_elements

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Radial mesh of a site as in self.f90 (logarithmic rofi, DRDI = a*(rofi+b)):
   !> r(ir) and the integration weights wd(ir) = Simpson weight * DRDI, ir = 2..nr
   !> (r(1) = 0 and wd(1) = 0: the origin does not contribute).
   !---------------------------------------------------------------------------
   subroutine radial_mesh(bnd, plusbulk, r, wd)
      class(bands), intent(in) :: bnd
      integer, intent(in) :: plusbulk
      real(rp), allocatable, intent(out) :: r(:), wd(:)
      integer :: nr, ir
      real(rp) :: a, b, ea, rpb, wgt

      a = bnd%symbolic_atom(plusbulk)%a
      nr = bnd%symbolic_atom(plusbulk)%mesh_grid_size()
      b = bnd%symbolic_atom(plusbulk)%potential%ws_r/(exp(a*(nr - 1)) - 1.0_rp)
      allocate (r(nr), wd(nr))
      r(:) = 0.0_rp
      wd(:) = 0.0_rp
      ea = exp(a)
      rpb = b*ea
      do ir = 2, nr
         r(ir) = rpb - b
         ! Standard Simpson weights used throughout the atomic code.
         wgt = 2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp
         if (ir == 2 .or. ir == nr) wgt = 1.0_rp/3.0_rp
         wd(ir) = wgt*(a*rpb)
         rpb = rpb*ea
      end do
   end subroutine radial_mesh

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> General on-site multipoles Q_LM = int rho(r) r^L R_LM(r_hat) d^3r, L = 0..2*lmax, of every
   !> recursion atom from the raw energy moments of the on-site density matrix, written to multipoles.out.
   !>
   !> Derivation. With psi_l(E) = phi + (E - E_l) phidot + (E - E_l)^2/2 phiddot (radial functions at the
   !> expansion energy E_l, the EVAL of self.f90 shifted by vmad to the Hamiltonian energy scale) and the density
   !> operator rho = sum_ab |a> n_ab <b|, n_ab = -(1/2 pi i) int (G_ab - G_ba^*) dE, the pair term of rho(r) is, to
   !> second order in E - E_l,
   !>    n00 phi phi' + n10 phidot phi' + n01 phi phidot' + n11 phidot phidot' + (1/2)(n20 phiddot phi' + n02 phi phiddot'),
   !> with n^(ij)_ab = int (E - E_a)^i (E - E_b)^j n_ab(E) dE (E_a = E_l of the shell of a), from the raw moments
   !> M_p = int E^p n dE:  n00 = M0, n10 = M1 - E_a M0, n01 = M1 - E_b M0, n11 = M2 - (E_a + E_b) M1 + E_a E_b M0,
   !> n20 = M2 - 2 E_a M1 + E_a^2 M0, n02 = M2 - 2 E_b M1 + E_b^2 M0. For a = b (same shell) this is
   !> Q0 phi^2 + 2 Q1 phi phidot + Q2 (phidot^2 + phi phiddot), the spherical density of self.f90 NEWRHO. The radial
   !> products use its weights: sqrt(GFAC_a GFAC_b) f_a f_b (large component) + f_a f_b (small component, JR = IR + NR).
   !> rho(r) = sum_ab psi_a(r) n_ab conj(psi_b(r)), psi_a = u_l(r)/r b_a(r_hat), so Q_LM = sum_ab n_ab (radial_ab,L)
   !> int b_a conj(b_b) R_LM dOmega (real_gaunt_table). Only the spin-diagonal blocks enter (the charge density).
   !>
   !> @param[inout] bnd      Bands object (on-site Chebyshev moments, atoms).
   !> @param[in]    rad_fun  Radial functions (mesh, l+1, spin, 1:7, local atom): phi, phidot, phiddot (large), the same
   !>                        three (small component), GFAC; exported from NEWRHO.
   !> @param[in]    rad_enu  Expansion energy (l+1, spin, local atom) of rad_fun (radial energy scale).
   !> @param[in]    set_q10  If true, potential%q10 = (z-moment)/S = sqrt(4 pi/3) Q_10/S, mixed with mix as compute_dipole_moments.
   !> @param[in]    mix      Linear mixing factor of q10.
   !---------------------------------------------------------------------------
   subroutine compute_multipoles(bnd, rad_fun, rad_enu, set_q10, mix)
      class(bands), intent(inout) :: bnd
      real(rp), dimension(:, :, :, :, :), intent(in) :: rad_fun
      real(rp), dimension(:, :, :), intent(in) :: rad_enu
      logical, intent(in) :: set_q10
      real(rp), intent(in) :: mix
      integer :: nsh, lm2, ncomp, na_glob, na_loc, plusbulk, isp, so, ia, ib, t, l, m, u
      real(rp) :: herm, q10_old, vmad
      real(rp), dimension(:), allocatable :: r, wd, e_shell
      real(rp), dimension(:, :), allocatable :: q_all
      real(rp), dimension(:, :), allocatable :: qs
      real(rp), dimension(:, :, :, :), allocatable :: w
      complex(rp), dimension(:, :, :, :), allocatable :: a, mp
      complex(rp), dimension(:, :, :), allocatable :: n
      real(rp), dimension(:, :), allocatable :: ea, eb
      complex(rp), dimension(:, :), allocatable :: m0, m1, m2
#ifdef USE_MPI
      integer :: ierr
#endif

      nsh = spin_off
      lm2 = 2*lmax_basis
      ncomp = (lm2 + 1)**2
      if (size(rad_fun, 2) < lmax_basis + 1) call g_logger%fatal('compute_multipoles: radial functions miss l <= lmax', __FILE__, __LINE__)
      call real_gaunt_table(lmax_basis, a)
      call bnd%green%chebyshev_onsite_energy_moments(mp, herm)
#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, herm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
      if (rank == 0) call g_logger%info('multipoles: max |M_p - M_p^H| = '//fmt('es12.4', herm), __FILE__, __LINE__)

      allocate (q_all(ncomp, bnd%lattice%nrec), qs(0:lm2, -lm2:lm2), n(nsh, nsh, 6), e_shell(nsh))
      allocate (w(0:lmax_basis, 0:lmax_basis, 0:lm2, 6), ea(nsh, nsh), eb(nsh, nsh), m0(nsh, nsh), m1(nsh, nsh), m2(nsh, nsh))
      q_all(:, :) = 0.0_rp
      do na_glob = start_atom, end_atom
         na_loc = g2l_map(na_glob)
         plusbulk = bnd%lattice%nbulk + na_glob
         vmad = bnd%symbolic_atom(plusbulk)%potential%vmad
         call radial_mesh(bnd, plusbulk, r, wd)
         qs(:, :) = 0.0_rp
         do isp = 1, 2
            so = (isp - 1)*spin_off
            call radial_pair_weights(rad_fun(1:size(r), :, isp, :, na_loc), r, wd, lmax_basis, w)
            do ia = 1, nsh
               e_shell(ia) = rad_enu(shell_l(ia) + 1, isp, na_loc) + vmad
            end do
            ea = spread(e_shell, 2, nsh)
            eb = spread(e_shell, 1, nsh)
            m0 = mp(so + 1:so + nsh, so + 1:so + nsh, 0, na_loc)
            m1 = mp(so + 1:so + nsh, so + 1:so + nsh, 1, na_loc)
            m2 = mp(so + 1:so + nsh, so + 1:so + nsh, 2, na_loc)
            n(:, :, 1) = m0
            n(:, :, 2) = m1 - ea*m0
            n(:, :, 3) = m1 - eb*m0
            n(:, :, 4) = m2 - (ea + eb)*m1 + ea*eb*m0
            n(:, :, 5) = m2 - 2.0_rp*ea*m1 + ea*ea*m0
            n(:, :, 6) = m2 - 2.0_rp*eb*m1 + eb*eb*m0
            do t = 1, 6
               call add_multipole_block(lmax_basis, a, n(:, :, t), w(:, :, :, t), qs)
            end do
         end do
         do l = 0, lm2
            do m = -l, l
               q_all(l*l + l + m + 1, na_glob) = qs(l, m)
            end do
         end do
      end do
#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, q_all, size(q_all), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      if (rank == 0) then
         open (newunit=u, file='multipoles.out', status='replace', action='write')
         write (u, '(a)') '# site L M Q_LM [electrons*bohr^L]; real harmonics R_LM (M>0 cos, M<0 sin, Condon-Shortley), L = 0..2*lmax'
         do ia = 1, bnd%lattice%nrec
            do l = 0, lm2
               do m = -l, l
                  write (u, '(i6,2i4,es24.15)') ia, l, m, q_all(l*l + l + m + 1, ia)
               end do
            end do
         end do
         close (u)
      end if

      if (set_q10) then
         do ia = 1, bnd%lattice%nrec
            plusbulk = bnd%lattice%nbulk + ia
            q10_old = bnd%symbolic_atom(plusbulk)%potential%q10
            bnd%symbolic_atom(plusbulk)%potential%q10 = mix*sqrt(4.0_rp*pi/3.0_rp)*q_all(3, ia)/bnd%symbolic_atom(plusbulk)%potential%ws_r + &
                                                         (1.0_rp - mix)*q10_old
         end do
      end if
      deallocate (q_all, qs, n, e_shell, w, ea, eb, m0, m1, m2, a, mp)
   end subroutine compute_multipoles

   !> Shell l of the in-shell-ordered orbital index i (shell l holds l*l+1 .. (l+1)**2).
   pure integer function shell_l(i)
      integer, intent(in) :: i
      shell_l = int(sqrt(real(i - 1, rp)))
   end function shell_l

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Radial weights w(l1, l2, L, t) = int F_t(r) r^L dr for the six pair terms t = 1..6 of compute_multipoles
   !> (00, 10, 01, 11, 20, 02), F_t = c_t [sqrt(GFAC_1 GFAC_2) f1_large f2_large + f1_small f2_small] with
   !> (f1, f2) = (phi, phi), (phidot, phi), (phi, phidot), (phidot, phidot), (phiddot, phi), (phi, phiddot) and c_t = 1, 1, 1, 1, 1/2, 1/2.
   !> rf(mesh, l+1, 1:7) as rad_fun for one atom and spin; r, wd from radial_mesh.
   !---------------------------------------------------------------------------
   subroutine radial_pair_weights(rf, r, wd, lmax, w)
      real(rp), dimension(:, :, :), intent(in) :: rf
      real(rp), dimension(:), intent(in) :: r, wd
      integer, intent(in) :: lmax
      real(rp), dimension(0:, 0:, 0:, :), intent(out) :: w
      integer, dimension(6), parameter :: k1 = [1, 2, 1, 2, 3, 1], k2 = [1, 1, 2, 2, 1, 3]
      real(rp), dimension(6), parameter :: c = [1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp, 0.5_rp, 0.5_rp]
      integer :: l1, l2, ll, t
      real(rp), dimension(size(r)) :: f, g12

      do l2 = 0, lmax
         do l1 = 0, lmax
            g12 = sqrt(rf(:, l1 + 1, 7)*rf(:, l2 + 1, 7))
            do t = 1, 6
               f = c(t)*(g12*rf(:, l1 + 1, k1(t))*rf(:, l2 + 1, k2(t)) + rf(:, l1 + 1, 3 + k1(t))*rf(:, l2 + 1, 3 + k2(t)))
               do ll = 0, 2*lmax
                  w(l1, l2, ll, t) = sum(wd*r**ll*f)
               end do
            end do
         end do
      end do
   end subroutine radial_pair_weights

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Angular table a(ia, ib, L, M) = int b_ia conj(b_ib) R_LM dOmega of the orbital basis of the code
   !> (b_lm = (-1)^m Y_lm Condon-Shortley, m = +l..-l within each shell, index l*l + (l - m) + 1) and the real
   !> harmonics R_LM = Y_L0 (M = 0), sqrt2 (-1)^M Re Y_LM (M > 0), sqrt2 (-1)^M Im Y_L|M| (M < 0). With
   !> b_ia conj(b_ib) = (-1)^m1 Y_{l1 m1} Y_{l2 -m2}, the integral is the complex Gaunt coefficient of math_mod.
   !---------------------------------------------------------------------------
   subroutine real_gaunt_table(lmax, a)
      integer, intent(in) :: lmax
      complex(rp), dimension(:, :, :, :), allocatable, intent(out) :: a
      integer :: nsh, l1, l2, m1, m2, ll, mm, ia, ib, mu
      complex(rp) :: c, gp, gm
      complex(rp), parameter :: i_u = (0.0_rp, 1.0_rp)

      nsh = (lmax + 1)**2
      allocate (a(nsh, nsh, 0:2*lmax, -2*lmax:2*lmax))
      a(:, :, :, :) = (0.0_rp, 0.0_rp)
      do l1 = 0, lmax
         do m1 = l1, -l1, -1
            ia = l1*l1 + (l1 - m1) + 1
            do l2 = 0, lmax
               do m2 = l2, -l2, -1
                  ib = l2*l2 + (l2 - m2) + 1
                  do ll = abs(l1 - l2), l1 + l2
                     if (mod(l1 + l2 + ll, 2) /= 0) cycle
                     do mm = -ll, ll
                        mu = abs(mm)
                        gp = gaunt(l1, l2, ll, m1, -m2, mu)
                        gm = gaunt(l1, l2, ll, m1, -m2, -mu)
                        if (mm == 0) then
                           c = gp
                        else if (mm > 0) then
                           c = (real((-1)**mu, rp)*gp + gm)/sqrt(2.0_rp)
                        else
                           c = real((-1)**mu, rp)/(i_u*sqrt(2.0_rp))*(gp - real((-1)**mu, rp)*gm)
                        end if
                        a(ia, ib, ll, mm) = real((-1)**m1, rp)*c
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine real_gaunt_table

   !> q(L, M) += Re sum_ab n(a, b) a(a, b, L, M) w(l_a, l_b, L): one spin block n of the density matrix, a from real_gaunt_table, radial weights w.
   subroutine add_multipole_block(lmax, a, n, w, q)
      integer, intent(in) :: lmax
      complex(rp), dimension(:, :, 0:, -2*lmax:), intent(in) :: a
      complex(rp), dimension(:, :), intent(in) :: n
      real(rp), dimension(0:, 0:, 0:), intent(in) :: w
      real(rp), dimension(0:, -2*lmax:), intent(inout) :: q
      integer :: ia, ib, ll, mm

      do mm = -2*lmax, 2*lmax
         do ll = abs(mm), 2*lmax
            do ib = 1, size(n, 2)
               do ia = 1, size(n, 1)
                  q(ll, mm) = q(ll, mm) + real(n(ia, ib)*a(ia, ib, ll, mm), rp)*w(shell_l(ia), shell_l(ib), ll)
               end do
            end do
         end do
      end do
   end subroutine add_multipole_block

end module electrostatics_multipole_mod
