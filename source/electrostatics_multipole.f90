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
   use mpi_mod, only: rank, start_atom, end_atom, g2l_map
#ifdef USE_MPI
   use mpi
#endif
   use bands_mod, only: bands
   use precision_mod, only: rp
   use string_mod, only: fmt
   use logger_mod, only: g_logger
   use math_mod, only: simpson_m
   use basis_mod, only: spin_off
   implicit none

   private

   real(rp), parameter :: pi = 3.14159265358979323846_rp

   !> Orbital indices in the SPHERICAL-harmonic basis that g0 lives in. The
   !> Hamiltonian is transformed cartesian -> spherical via hcpx(...,’cart2sph’)
   !> on every on-site block (hamiltonian_build.f90 build_obarm/build_enim/hhmag),
   !> so g0 is ordered by m as (math.f90 hcpx, n=9):
   !>   (00)(1-1)(10)(11)(2-2)(2-1)(20)(21)(22)
   !> => 1=s, 2=p(m=-1), 3=p(m=0), 4=p(m=+1),
   !>    5=d(m=-2), 6=d(m=-1), 7=d(m=0), 8=d(m=+1), 9=d(m=+2).
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

   public :: compute_dipole_moments

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
      real(rp), dimension(bnd%en%channels_ldos + 10) :: y
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
            y(:) = -aimag(bnd%green%g0(idx_s + soff, idx_pz + soff, :, na_loc) + &
                          bnd%green%g0(idx_pz + soff, idx_s + soff, :, na_loc))/pi
            call simpson_m(d_spz, bnd%en%edel, bnd%en%fermi, bnd%nv1, y, bnd%e1, 0, bnd%en%ene)

            y(:) = -aimag(bnd%green%g0(idx_pz + soff, idx_dz2 + soff, :, na_loc) + &
                          bnd%green%g0(idx_dz2 + soff, idx_pz + soff, :, na_loc))/pi
            call simpson_m(d_pzd, bnd%en%edel, bnd%en%fermi, bnd%nv1, y, bnd%e1, 0, bnd%en%ene)

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
         ! if (rank == 0) then
         !    call g_logger%info(’B6 dipole moment Q10 of atom’//fmt(’i4’, na_glob)// &
         !                       ’ = ’//fmt(’f12.8’, q10_all(na_glob)), __FILE__, __LINE__)
         ! end if
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
      real(rp) :: a, b, ea, rpb, r, drdi, wgt, ws_r
      real(rp) :: phi_s, phi_p, phi_d

      r_sp = 0.0_rp
      r_pd = 0.0_rp

      a = bnd%symbolic_atom(plusbulk)%a
      nr = bnd%symbolic_atom(plusbulk)%mesh_grid_size()
      if (nr < 3) return
      ws_r = bnd%symbolic_atom(plusbulk)%potential%ws_r
      b = ws_r/(exp(a*(nr - 1)) - 1.0_rp)

      ea = exp(a)
      ! ir = 1 is the origin (rofi = 0, zero contribution); start at ir = 2.
      rpb = b*ea
      do ir = 2, nr
         r = rpb - b
         drdi = a*rpb                    ! DRDI = a*(rofi + b) = a*rpb
         ! Standard Simpson weights used throughout the atomic code.
         wgt = 2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp
         if (ir == 2 .or. ir == nr) wgt = 1.0_rp/3.0_rp

         ! Spin-averaged large-component amplitudes u_l(r).
         phi_s = 0.5_rp*(phi_amp(ir, idx_l(idx_s),   1) + phi_amp(ir, idx_l(idx_s),   2))
         phi_p = 0.5_rp*(phi_amp(ir, idx_l(idx_pz),  1) + phi_amp(ir, idx_l(idx_pz),  2))
         phi_d = 0.5_rp*(phi_amp(ir, idx_l(idx_dz2), 1) + phi_amp(ir, idx_l(idx_dz2), 2))

         ! ASA+M multipole normalization: the l=1 moment uses (r/S)^l with l=1,
         ! i.e. the dipole length is measured in Wigner-Seitz radii (Vitos/Skriver;
         ! matches dsz carrying FACP = S*sqrt(3)). Hence the r/ws_r weight.
         r_sp = r_sp + wgt*drdi*phi_s*phi_p*(r/ws_r)
         r_pd = r_pd + wgt*drdi*phi_p*phi_d*(r/ws_r)

         rpb = rpb*ea
      end do
   end subroutine radial_dipole_elements

   !> Map an orbital index (real-harmonics spd basis) to its (l+1) partial-wave
   !> channel used to index phi_amp / fun2: s->1, p->2, d->3.
   pure integer function idx_l(iorb)
      integer, intent(in) :: iorb
      if (iorb == idx_s) then
         idx_l = 1
      else if (iorb == idx_pz) then
         idx_l = 2
      else
         idx_l = 3
      end if
   end function idx_l

end module electrostatics_multipole_mod
