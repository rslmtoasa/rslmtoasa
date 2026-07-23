!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief
!> Higher-order (l=1 dipole) multipole moments for surface/interface
!> electrostatics — RS-LMTO-ASA audit item B6.
!>
!> The ASA Madelung problem is extended from monopoles (q_R) to include the
!> l=1 (z) dipole moment Q_{10} of each sphere's charge. The generalized
!> Madelung matrices (dipole-monopole dsz, dipole-dipole dzz) are already
!> built by charge%madl2d following Skriver-Rosengaard, Phys. Rev. B 43, 9538;
!> this module supplies the missing charge-side moment Q_{10}, which
!> charge%surfpot then couples into the local potential shift via dsz.
!>
!> This is the *new* code that lives outside the legacy self.f90 fence
!> (CLAUDE.md rule 5): it reads the on-site Green function (bands%green%g0)
!> and the exported radial partial-wave amplitudes (self%phi_amp), and writes
!> only the derived scalar potential%q10.
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
   use self_mod, only: self
   use bands_mod, only: bands
   use symbolic_atom_mod, only: symbolic_atom
   use lattice_mod, only: lattice
   use energy_mod, only: energy
   use math_mod, only: simpson_m, sqrt_three, pi
   use precision_mod, only: rp
   use string_mod, only: fmt
   use logger_mod, only: g_logger
   use basis_mod, only: spin_off
   implicit none

   private

   !> Orbital indices in the real-harmonics spd basis used throughout the code
   !> (see math_mod: s(1), p_x(2),p_y(3),p_z(4), d_xy(5),d_yz(6),d_zx(7),
   !> d_x2y2(8), d_3z2r2(9)). The l=1 (z) dipole couples s<->p_z and p_z<->d_z2.
   integer, parameter :: idx_s    = 1
   integer, parameter :: idx_pz   = 4
   integer, parameter :: idx_dz2  = 9

   !> Angular (Gaunt) factors  C = <Y_{L1}| cos(theta) |Y_{L2}>  for real
   !> spherical harmonics: <s|cos|p_z> = 1/sqrt(3), <p_z|cos|d_z2> = 2/sqrt(15).
   real(rp), parameter :: c_sp = 1.0_rp/sqrt(3.0_rp)
   real(rp), parameter :: c_pd = 2.0_rp/sqrt(15.0_rp)

   type, public :: electrostatics_multipole
      !> Parent self object (gives access to bands, lattice, energy, atoms and
      !> the exported radial amplitudes phi_amp).
      class(self), pointer :: self => null()
      class(bands), pointer :: bands => null()
      class(lattice), pointer :: lattice => null()
      class(energy), pointer :: en => null()
   contains
      final :: destructor
      procedure :: restore_to_default
      procedure :: compute_dipole_moments
   end type electrostatics_multipole

   interface electrostatics_multipole
      procedure :: constructor
   end interface

contains

   !---------------------------------------------------------------------------
   !> @brief Constructor. Wires up pointers from the parent self object.
   !---------------------------------------------------------------------------
   function constructor(self_obj) result(obj)
      type(electrostatics_multipole) :: obj
      class(self), target, intent(in) :: self_obj

      call obj%restore_to_default()
      obj%self => self_obj
      obj%bands => self_obj%bands
      obj%lattice => self_obj%lattice
      obj%en => self_obj%en
   end function constructor

   subroutine destructor(this)
      type(electrostatics_multipole) :: this
      this%self => null()
      this%bands => null()
      this%lattice => null()
      this%en => null()
   end subroutine destructor

   subroutine restore_to_default(this)
      class(electrostatics_multipole), intent(inout) :: this
      this%self => null()
      this%bands => null()
      this%lattice => null()
      this%en => null()
   end subroutine restore_to_default

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
   !>      partial-wave amplitudes u_l(r) = self%phi_amp over the WS mesh.
   !>   3. Assemble Q_{10} = sum_sigma [ C_sp D_{s,pz} R_sp + C_pd D_{pz,dz2} R_pd ].
   !> The per-atom results are reduced across MPI ranks so rank 0 (which runs
   !> charge%surfpot) sees the full set.
   !---------------------------------------------------------------------------
   subroutine compute_dipole_moments(this)
      class(electrostatics_multipole), intent(inout) :: this
      integer :: na_glob, na_loc, plusbulk, isp, soff
      integer :: nrec
      real(rp) :: d_spz, d_pzd, r_sp, r_pd, q10
      real(rp), dimension(:), allocatable :: q10_all
      real(rp), dimension(this%en%channels_ldos + 10) :: y
#ifdef USE_MPI
      integer :: ierr
#endif

      nrec = this%lattice%nrec
      allocate (q10_all(nrec))
      q10_all(:) = 0.0_rp

      do na_glob = start_atom, end_atom
         na_loc = g2l_map(na_glob)
         plusbulk = this%lattice%nbulk + na_glob

         ! Radial dipole matrix elements from the current SCF partial waves.
         call radial_dipole_elements(this, plusbulk, r_sp, r_pd)

         q10 = 0.0_rp
         do isp = 1, 2
            soff = (isp - 1)*spin_off

            ! Symmetrized cross-l density-matrix elements (factor 2 accounts for
            ! the (L1,L2) and (L2,L1) orderings of the Hermitian on-site block).
            y(:) = -aimag(this%bands%green%g0(idx_s + soff, idx_pz + soff, :, na_loc) + &
                          this%bands%green%g0(idx_pz + soff, idx_s + soff, :, na_loc))/pi
            call simpson_m(d_spz, this%en%edel, this%en%fermi, this%bands%nv1, y, &
                           this%bands%e1, 0, this%en%ene)

            y(:) = -aimag(this%bands%green%g0(idx_pz + soff, idx_dz2 + soff, :, na_loc) + &
                          this%bands%green%g0(idx_dz2 + soff, idx_pz + soff, :, na_loc))/pi
            call simpson_m(d_pzd, this%en%edel, this%en%fermi, this%bands%nv1, y, &
                           this%bands%e1, 0, this%en%ene)

            q10 = q10 + c_sp*d_spz*r_sp + c_pd*d_pzd*r_pd
         end do

         q10_all(na_glob) = q10
      end do

#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, q10_all, nrec, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      do na_glob = 1, nrec
         this%self%symbolic_atom(this%lattice%nbulk + na_glob)%potential%q10 = q10_all(na_glob)
         if (rank == 0) then
            call g_logger%info('B6 dipole moment Q10 of atom'//fmt('i4', na_glob)// &
                               ' = '//fmt('f12.8', q10_all(na_glob)), __FILE__, __LINE__)
         end if
      end do

      deallocate (q10_all)
   end subroutine compute_dipole_moments

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Radial dipole matrix elements over the Wigner-Seitz sphere:
   !>   R_sp = int u_s(r) u_pz(r) r dr,   R_pd = int u_pz(r) u_dz2(r) r dr,
   !> where u_l(r) = r*phi_l(r) is the exported large-component amplitude
   !> self%phi_amp(:, l+1, spin). The mesh is reconstructed exactly as in
   !> self.f90 (logarithmic rofi, DRDI = a*(rofi+b)) and integrated with the
   !> code's standard Simpson weights. The spin-averaged amplitude is used for
   !> the radial factor (spin dependence enters through the density matrix).
   !---------------------------------------------------------------------------
   subroutine radial_dipole_elements(this, plusbulk, r_sp, r_pd)
      class(electrostatics_multipole), intent(inout) :: this
      integer, intent(in) :: plusbulk
      real(rp), intent(out) :: r_sp, r_pd
      integer :: nr, ir
      real(rp) :: a, b, ea, rpb, r, drdi, wgt
      real(rp) :: phi_s, phi_p, phi_d

      r_sp = 0.0_rp
      r_pd = 0.0_rp

      a = this%self%symbolic_atom(plusbulk)%a
      nr = this%self%symbolic_atom(plusbulk)%mesh_grid_size()
      if (nr < 3) return
      b = this%self%symbolic_atom(plusbulk)%potential%ws_r/(exp(a*(nr - 1)) - 1.0_rp)

      ea = exp(a)
      rpb = b
      ! ir = 1 is the origin (rofi = 0): zero contribution, start at ir = 2.
      rpb = b*ea
      do ir = 2, nr
         r = rpb - b
         drdi = a*(rpb)                 ! DRDI = a*(rofi + b) = a*rpb
         ! Standard Simpson weights used throughout the atomic code.
         wgt = 2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp
         if (ir == 2 .or. ir == nr) wgt = 1.0_rp/3.0_rp

         ! Spin-averaged large-component amplitudes u_l(r).
         phi_s = 0.5_rp*(this%self%phi_amp(ir, idx_l(idx_s),   1) + this%self%phi_amp(ir, idx_l(idx_s),   2))
         phi_p = 0.5_rp*(this%self%phi_amp(ir, idx_l(idx_pz),  1) + this%self%phi_amp(ir, idx_l(idx_pz),  2))
         phi_d = 0.5_rp*(this%self%phi_amp(ir, idx_l(idx_dz2), 1) + this%self%phi_amp(ir, idx_l(idx_dz2), 2))

         r_sp = r_sp + wgt*drdi*phi_s*phi_p*r
         r_pd = r_pd + wgt*drdi*phi_p*phi_d*r

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
