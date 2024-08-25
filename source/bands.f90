!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Bands
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module to handle calculation of energy bands and parameters related
!with it
!------------------------------------------------------------------------------

module bands_mod
   use mpi_mod
   use control_mod
   use energy_mod
   use green_mod
   use lattice_mod
   use symbolic_atom_mod
   use density_of_states_mod
   use recursion_mod
   use precision_mod, only: rp
   use math_mod
   use logger_mod, only: g_logger
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   use string_mod
#ifdef USE_MPI
   use mpi
#endif
   implicit none

   private

   !> Module´s main structure
   type, public :: bands
      !> Green
      class(green), pointer :: green
      !> Lattice
      class(lattice), pointer :: lattice
      !> Symbolic atom
      class(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Density of states
      class(dos), pointer :: dos
      !> Energy
      class(energy), pointer :: en
      !> Control
      class(control), pointer :: control
      !> Recursion
      class(recursion), pointer :: recursion

      ! General variables
      !> Energy variable
      real(rp) :: e1, e1cheb
      !> Total number of valence electrons
      real(rp) :: qqv
      !> Band energy mesh
      integer :: nv1, ik1, nv1cheb, ik1cheb
      !> Total density of states
      real(rp), dimension(:), allocatable :: dtot, dtotcheb
      !> Projected DOS
      real(rp), dimension(:, :), allocatable :: dx, dy, dz, dnmag, dup, ddw
      real(rp), dimension(:, :, :, :, :), allocatable :: d_orb
      !> Projected Green Function
      complex(rp), dimension(:, :, :, :), allocatable :: g0_x, g0_y, g0_z
      !> Energy bands (?)
      real(rp), dimension(:, :, :), allocatable :: dspd
      real(rp) :: eband
      !> Magnetic force
      real(rp), dimension(:, :), allocatable :: mag_for
      !> Calculated effictive Hubbard U (= U - J)
      real(rp), dimension(:,:), allocatable :: hubbard_u_eff_old
      !> Checks if the calculated effective Hubbard U has converged
      logical :: hubbard_u_converged
   contains
      procedure :: calculate_projected_green
      procedure :: calculate_projected_dos
      procedure :: calculate_orbital_dos
      procedure :: calculate_moments
      procedure :: calculate_pl
      procedure :: calculate_moments_gauss_legendre
      procedure :: calculate_moments_chebgauss
      procedure :: calculate_band_energy
      procedure :: calculate_magnetic_moments
      procedure :: calculate_orbital_moments
      procedure :: calculate_magnetic_torques
      procedure :: calculate_fermi
      procedure :: calculate_fermi_gauss
      procedure :: calculate_occupation_gauss_legendre
      procedure :: Hubbard_U_Potential ! LDA+U+J method
      procedure :: Hubbard_V_Potential
      procedure :: calc_hubbard_U
      procedure :: fermi
      procedure :: restore_to_default
      final :: destructor
   end type

   interface bands
      procedure :: constructor
   end interface bands

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] fname Namelist file
   !> @return type(bands)
   !---------------------------------------------------------------------------
   function constructor(green_obj) result(obj)
      type(bands) :: obj
      type(green), target, intent(in) :: green_obj

      obj%green => green_obj
      obj%lattice => green_obj%dos%recursion%lattice
      obj%symbolic_atom => green_obj%dos%recursion%hamiltonian%charge%lattice%symbolic_atoms
      obj%dos => green_obj%dos
      obj%en => green_obj%dos%en
      obj%control => green_obj%dos%recursion%lattice%control
      obj%recursion => green_obj%dos%recursion

      call obj%restore_to_default()
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(bands) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%dtot)) call g_safe_alloc%deallocate('bands.dtot', this%dtot)
      if (allocated(this%dtotcheb)) call g_safe_alloc%deallocate('bands.dtotcheb', this%dtotcheb)
      if (allocated(this%dx)) call g_safe_alloc%deallocate('bands.dx', this%dx)
      if (allocated(this%dy)) call g_safe_alloc%deallocate('bands.dy', this%dy)
      if (allocated(this%dz)) call g_safe_alloc%deallocate('bands.dz', this%dz)
      if (allocated(this%dnmag)) call g_safe_alloc%deallocate('bands.dnmag', this%dnmag)
      if (allocated(this%dup)) call g_safe_alloc%deallocate('bands.dup', this%dup)
      if (allocated(this%ddw)) call g_safe_alloc%deallocate('bands.ddw', this%ddw)
      if (allocated(this%g0_x)) call g_safe_alloc%deallocate('bands.g0_x', this%g0_x)
      if (allocated(this%g0_y)) call g_safe_alloc%deallocate('bands.g0_y', this%g0_y)
      if (allocated(this%g0_z)) call g_safe_alloc%deallocate('bands.g0_z', this%g0_z)
      if (allocated(this%dspd)) call g_safe_alloc%deallocate('bands.dspd', this%dspd)
      if (allocated(this%d_orb)) call g_safe_alloc%deallocate('bands.ddw', this%d_orb)
      if (allocated(this%mag_for)) call g_safe_alloc%deallocate('bands.mag_for', this%mag_for)
#else
      if (allocated(this%dtot)) deallocate (this%dtot)
      if (allocated(this%dtotcheb)) deallocate (this%dtotcheb)
      if (allocated(this%dx)) deallocate (this%dx)
      if (allocated(this%dy)) deallocate (this%dy)
      if (allocated(this%dz)) deallocate (this%dz)
      if (allocated(this%dnmag)) deallocate (this%dnmag)
      if (allocated(this%dup)) deallocate (this%dup)
      if (allocated(this%ddw)) deallocate (this%ddw)
      if (allocated(this%g0_x)) deallocate (this%g0_x)
      if (allocated(this%g0_y)) deallocate (this%g0_y)
      if (allocated(this%g0_z)) deallocate (this%g0_z)
      if (allocated(this%dspd)) deallocate (this%dspd)
      if (allocated(this%d_orb)) deallocate (this%d_orb)
      if (allocated(this%mag_for)) deallocate (this%mag_for)
      if (allocated(this%hubbard_u_eff_old)) deallocate (this%hubbard_u_eff_old)
#endif
   end subroutine destructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      use mpi_mod
      class(bands) :: this

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('bands.dtot', this%dtot, (/this%en%channels_ldos + 10/))
      call g_safe_alloc%allocate('bands.dtotcheb', this%dtotcheb, (/this%en%channels_ldos + 10/))
      call g_safe_alloc%allocate('bands.dx', this%dx, (/this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('bands.dy', this%dy, (/this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('bands.dz', this%dz, (/this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('bands.g0_x', this%g0_x, (/9, 9, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('bands.g0_y', this%g0_y, (/9, 9, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('bands.g0_z', this%g0_z, (/9, 9, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('bands.dspd', this%dspd, (/6, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('bands.d_orb', this%d_orb, (/9, 9, 3, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('bands.mag_for', this%mag_for, (/3, atoms_per_process/))
#else
      allocate (this%dtot(this%en%channels_ldos + 10))
      allocate (this%dtotcheb(this%en%channels_ldos + 10))
      allocate (this%dx(this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%dy(this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%dz(this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%g0_x(9, 9, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%g0_y(9, 9, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%g0_z(9, 9, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%dspd(6, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%d_orb(9, 9, 3, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%mag_for(3, atoms_per_process))
      if ( this%recursion%hamiltonian%hubbardU_sc_check ) then
         allocate (this%hubbard_u_eff_old(this%lattice%nrec, 4))
      end if
#endif
      this%dtot(:) = 0.0d0
      this%dtotcheb(:) = 0.0d0
      this%dx(:, :) = 0.0d0
      this%dy(:, :) = 0.0d0
      this%dz(:, :) = 0.0d0
      this%g0_x(:, :, :, :) = 0.0d0
      this%g0_y(:, :, :, :) = 0.0d0
      this%g0_z(:, :, :, :) = 0.0d0
      this%dspd(:, :, :) = 0.0d0
      this%d_orb(:, :, :, :, :) = 0.0d0
      this%mag_for(:, :) = 0.0d0
      if ( this%recursion%hamiltonian%hubbardU_sc_check) then
         this%hubbard_u_eff_old(:,:) = 0.0d0
      end if
      this%hubbard_u_converged = .false.
   end subroutine restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the Fermi energy
   !---------------------------------------------------------------------------
   subroutine calculate_fermi(this)
      use mpi_mod
      class(bands) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, ia, ik1_mag, ik1, nv1, ifail
      real(rp) :: e1_mag, ef_mag, e1
      integer :: ia_glob
      real(rp), dimension(this%lattice%nrec, this%en%channels_ldos + 10) :: dosia
      real(rp), dimension(this%lattice%nrec, 18, this%en%channels_ldos + 10) :: dosial
      
      e1_mag = 0.0d0
      ef_mag = 0.0d0
      ik1_mag = 0
      ik1 = 0
      this%dtot(:) = 0.0d0
      dosial = 0.0d0
      dosia = 0.0d0

      ik1 = this%en%ik1
      nv1 = this%en%nv1

      ! Define the valence electrons from the bulk parameters
      this%qqv = real(sum(this%symbolic_atom(1:this%lattice%nbulk_bulk)%element%valence))
      if (rank == 0) call g_logger%info('Valence is:'//fmt('f16.6', this%qqv), __FILE__, __LINE__)
      ! Calculate the total density of states
      !do ia=1, this%lattice%nrec
      !call get_mpi_variables(rank,this%lattice%nrec)
      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)

         do i = 1, this%en%channels_ldos + 10
            do j = 1, 9
               this%dtot(i) = this%dtot(i) - aimag(this%green%g0(j, j, i, ia) + this%green%g0(j + 9, j + 9, i, ia))/pi
               dosia(ia_glob, i) = dosia(ia_glob, i) - aimag(this%green%g0(j, j, i, ia) + this%green%g0(j + 9, j + 9, i, ia))/pi
               dosial(ia_glob, j, i) = -aimag(this%green%g0(j, j, i, ia))/pi
               dosial(ia_glob, j + 9, i) = -aimag(this%green%g0(j + 9, j + 9, i, ia))/pi
            end do
            !write(250+ia,*) this%en%ene(i), dosia(ia,i)
            !write(450+ia,´(19f10.6)´) this%en%ene(i), dosial(ia,1,i), &
            !dosial(ia,2,i), dosial(ia,3,i), dosial(ia,4,i), &
            !dosial(ia,5,i), dosial(ia,6,i), dosial(ia,7,i), dosial(ia,8,i), dosial(ia,9,i), &
            !dosial(ia,10,i), &
            !dosial(ia,11,i), dosial(ia,12,i), dosial(ia,13,i), &
            !dosial(ia,14,i), dosial(ia,15,i), dosial(ia,16,i), dosial(ia,17,i), dosial(ia,18,i)
         end do
         !rewind(250+ia)
         !rewind(450+ia)
      end do
      ! Transfer total DOS across MPI ranks for Fermi surface determination.
#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%dtot, this%en%channels_ldos + 10, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, dosia, product(shape(dosia)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, dosial, product(shape(dosial)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      if (rank == 0) then
         do i = 1, this%en%channels_ldos + 10
            write (125, '(2f16.5)') this%en%ene(i), this%dtot(i)
         end do
         rewind (125)
      end if
      if (rank == 0) then
         do ia = 1, this%lattice%nrec
            do i = 1, this%en%channels_ldos + 10
               write (250 + ia, *) this%en%ene(i), dosia(ia, i)
            end do
            rewind (250 + ia) 
         end do
      end if
      if (rank == 0) then
         do ia = 1, this%lattice%nrec
            do i = 1, this%en%channels_ldos + 10
               write (450 + ia, '(19f10.6)') this%en%ene(i), dosial(ia, 1:18, i)
            end do
            rewind (450 + ia)
         end do
      end if

      ! Calculate the Fermi enery
      ef_mag = this%en%fermi
      this%en%chebfermi = this%en%fermi
      if (.not. (this%en%fix_fermi) .and. this%control%calctype == 'B') then
         e1_mag = ef_mag
         call this%fermi(ef_mag, this%en%edel, ik1_mag, this%en%energy_min, this%en%channels_ldos + 10, this%dtot, ifail, this%qqv, e1_mag)
         e1 = this%en%fermi
         call this%fermi(this%en%fermi, this%en%edel, ik1, this%en%energy_min, this%en%channels_ldos + 10, this%dtot, ifail, this%qqv, e1_mag)
         this%nv1 = ik1
         this%e1 = e1_mag
         if (rank == 0) call g_logger%info('Free Fermi energy:'//fmt('f10.6', this%en%fermi), __FILE__, __LINE__)
      else if (this%en%fix_fermi) then
         ik1 = nint((this%en%fermi - this%en%energy_min)/this%en%edel)
         e1 = this%en%energy_min + (ik1 - 1)*this%en%edel
         this%nv1 = ik1
         this%e1 = e1
         if (rank == 0) call g_logger%info('Fixed Fermi energy:'//fmt('f10.6', this%en%fermi), __FILE__, __LINE__)
      end if
   end subroutine calculate_fermi

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the band energy
   !---------------------------------------------------------------------------
   subroutine calculate_band_energy(this)
      class(bands) :: this
      ! Local variable

      call simpson_m(this%eband, this%en%edel, this%en%fermi, this%nv1, this%dtot, this%e1, 1, this%en%ene)
   end subroutine

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Auxiliar routine to calculate the Fermi enery
   !---------------------------------------------------------------------------
   subroutine fermi(this, EF, H, IK1, AINF, NPTS, Y, IFAIL, QQV, E1)
      class(bands) :: this
      ! Input
      integer, intent(in) :: NPTS
      integer, intent(inout) :: IK1
      integer, intent(out) :: IFAIL
      real(rp), intent(in) :: AINF, H, QQV
      real(rp), intent(inout) :: EF
      real(rp), intent(out) :: E1
      real(rp), dimension(NPTS), intent(in) :: Y
      ! Local variables
      integer :: i
      real(rp) :: AINT, AINT0, ALPHA

      IFAIL = 1
      AINT = 0.d0
      AINT0 = 0.d0
      do I = 2, NPTS - 1, 2
         AINT = AINT + H*(Y(I - 1) + 4.d0*Y(I) + Y(I + 1))/3.d0 !*fermifun(AINF+(I)*H, Ef, kBT)
         if (AINT >= QQV) goto 1000
         AINT0 = AINT
      end do
      return
1000  IFAIL = 0
      !print ´(3x, a, 2f12.6)´, ´Fermil´, QQV, AINT
      if (AINT == QQV) then
         IK1 = I + 1
         EF = AINF + H*I
         E1 = EF
      else
         ALPHA = (AINT - AINT0)/2.d0/H
         IK1 = I - 1
         E1 = AINF + H*(I - 2)
         EF = ((QQV - AINT0)/ALPHA) + E1
      end if
      !print ´(3x, a, 4f12.6)´, ´Fermil´, QQV, AINT, EF, E1
   end subroutine fermi

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the moments m^(q), q = 0, 1 and 2
   !---------------------------------------------------------------------------
   subroutine calculate_moments(this)
      use mpi_mod
      class(bands) :: this
      ! Local variables
      integer :: i, j, l, m, o ! Orbital index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index
      integer :: isp, soff, jo, nsp, plusbulk
      real(rp) :: sgef, pmef, smef, isgn, sums, sump, sumd, mnorm
      real(rp), dimension(this%en%channels_ldos + 10) :: y
      real(rp), dimension(18) :: chebmom(18), chebmom1(18), chebmom2(18)
      real(rp), dimension(this%lattice%nrec, 6) :: occ

      real(rp), dimension(:), allocatable :: pot_arr
      integer :: na_glob, pot_size
      real(rp), dimension(:, :), allocatable :: T_comm
      real(rp), dimension(3, atoms_per_process) :: mom_prev

      ! Store previous moments if needed for reverse rotation later
      do na_glob = start_atom, end_atom
         na = g2l_map(na_glob)
         plusbulk = this%lattice%nbulk + na_glob
         mom_prev(:, na) = this%symbolic_atom(plusbulk)%potential%mom
      end do

      !call this%calculate_magnetic_moments()
      if (this%control%nsp == 2 .or. this%control%nsp == 4) call this%calculate_orbital_moments()

      this%dspd(:, :, :) = 0.0d0
      !do na=1, this%lattice%nrec
      do na_glob = start_atom, end_atom
         na = g2l_map(na_glob)
         plusbulk = this%lattice%nbulk + na_glob
         do isp = 1, 2
            isgn = (-1.0d0)**(isp - 1)
            soff = 3*(isp - 1)
            do l = 1, 3
               do m = 1, 2*l - 1
                  o = (l - 1)**2 + m
                  do ie = 1, this%en%channels_ldos
                     this%dspd(l + soff, ie, na) = this%dspd(l + soff, ie, na) - aimag(this%green%g0(o, o, ie, na) + this%green%g0(o + 9, o + 9, ie, na)) - &
                                                   isgn*this%symbolic_atom(plusbulk)%potential%mom(3)*aimag(this%green%g0(o, o, ie, na) - this%green%g0(o + 9, o + 9, ie, na)) &
                                                   - isgn*this%symbolic_atom(plusbulk)%potential%mom(2)*aimag(i_unit*this%green%g0(o, o + 9, ie, na) - i_unit*this%green%g0(o + 9, o, ie, na)) &
                                                   - isgn*this%symbolic_atom(plusbulk)%potential%mom(1)*aimag(this%green%g0(o, o + 9, ie, na) + this%green%g0(o + 9, o, ie, na))
                  end do
               end do
            end do
         end do
         ! Rotate moments back to global frame if rotated earlier
         if (this%recursion%hamiltonian%local_axis) then
            call updatrotmom_single(this%symbolic_atom(plusbulk)%potential%mom, mom_prev(:, na))
            call this%symbolic_atom(plusbulk)%potential%copy_mom_to_scal()

            mnorm = this%symbolic_atom(plusbulk)%potential%mtot
            call g_logger%info('Global spin moment projections of atom'//fmt('i4', na)//' is '// &
                               fmt('f10.6', this%symbolic_atom(plusbulk)%potential%mom(1)*mnorm)//' '// &
                               fmt('f10.6', this%symbolic_atom(plusbulk)%potential%mom(2)*mnorm)//' '// &
                               fmt('f10.6', this%symbolic_atom(plusbulk)%potential%mom(3)*mnorm), __FILE__, __LINE__)
         end if
      end do

      this%dspd(:, :, :) = this%dspd(:, :, :)*0.5d0/pi

      occ = 0.0d0; sums = 0.0d0; sump = 0.0d0; sumd = 0.0d0
      !do na=1, this%lattice%nrec
      do na_glob = start_atom, end_atom
         na = g2l_map(na_glob)
         do i = 1, 6
            y(:) = 0.0d0
            if (i > 3) then
               nsp = 2
            else
               nsp = 1
            end if
            soff = 3*(nsp - 1)
            y(:) = this%dspd(i, :, na)
            sgef = 0.0d0; pmef = 0.0d0; smef = 0.0d0
            call simpson_m(sgef, this%en%edel, this%en%fermi, this%nv1, y, this%e1, 0, this%en%ene)
            call simpson_m(pmef, this%en%edel, this%en%fermi, this%nv1, y, this%e1, 1, this%en%ene)
            call simpson_m(smef, this%en%edel, this%en%fermi, this%nv1, y, this%e1, 2, this%en%ene)

            occ(na_glob, i) = sgef

            this%symbolic_atom(this%lattice%nbulk + na_glob)%potential%gravity_center(i - soff, nsp) = (pmef/sgef) - this%symbolic_atom(this%lattice%nbulk + na_glob)%potential%vmad
            this%symbolic_atom(this%lattice%nbulk + na_glob)%potential%ql(1, i - soff - 1, nsp) = sgef
            this%symbolic_atom(this%lattice%nbulk + na_glob)%potential%ql(2, i - soff - 1, nsp) = 0.0d0
            this%symbolic_atom(this%lattice%nbulk + na_glob)%potential%ql(3, i - soff - 1, nsp) = smef - 2.0d0*(pmef/sgef)*pmef + ((pmef/sgef)**2)*sgef
         end do
      end do

      call this%calculate_pl()

      ! Transfer calculated moments across MPI ranks
#ifdef USE_MPI
      pot_size = this%symbolic_atom(start_atom)%potential%sizeof_potential_lite()
      allocate (T_comm(pot_size, this%lattice%nrec))
      T_comm = 0.0_rp
      do na_glob = start_atom, end_atom
         call this%symbolic_atom(this%lattice%nbulk + na_glob)%potential%flatten_potential_lite(T_comm(:, na_glob))
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE, T_comm, product(shape(T_comm)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      do na_glob = 1, this%lattice%nrec
         call this%symbolic_atom(this%lattice%nbulk + na_glob)%potential%expand_potential_lite(T_comm(:, na_glob))
      end do
      deallocate (T_comm)
#endif

      !do na=1,this%lattice%nrec
      !write(*,*) ´Magnetic Moment for atom ´, na,´: ´, sum(occ(na,1:3)) - sum(occ(na,4:6))
      !write(*,*) ´Total Charge for atom ´,na, ´: ´, sum(occ(na,1:6)), ´s = ´, occ(na,1)+occ(na,4), ´p = ´, (occ(na,2))+(occ(na,5)), ´d = ´, (occ(na,3))+(occ(na,6))
      !sums = sums + occ(na,1)+occ(na,4)
      !sump = sump + occ(na,2)+occ(na,5)
      !sumd = sumd + occ(na,3)+occ(na,6)
      !end do
      !write(*,*) ´Total electrons s, p and d:´, sums, sump, sumd
   end subroutine calculate_moments

   subroutine calculate_moments_gauss_legendre(this)
      use mpi_mod
      class(bands) :: this
      ! Local variables
      integer :: i, j, l, m, o ! Orbital index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index
      integer :: isp, soff, jo, nsp
      integer :: fermi_point
      !complex(rp), dimension(18,18,this%lattice%nrec) :: g0_ef
      complex(rp), dimension(18, 18, atoms_per_process) :: g0_ef
      complex(rp) :: eta
      real(rp) :: sgef, pmef, smef, isgn, sumocc
      !real(rp), dimension(64,18,18,this%lattice%nrec) :: y
      real(rp), dimension(64, 18, 18, atoms_per_process) :: y
      real(rp), dimension(64) :: x, w
      real(rp) :: res, t
      real(rp), dimension(this%lattice%nrec, 18) :: occ

      integer :: m_glob

      ! Find the Gauss Legendre roots and weights
      call gauss_legendre(64, 0.0_rp, 1.0_rp, x, w)

      !call this%calculate_fermi_gauss(724.0d0,0.00001d0, 100, this%en%fermi)

      call this%en%e_mesh()

      do i = 1, this%en%channels_ldos + 10
         if ((this%en%ene(i) - this%en%fermi) .le. 0.0001d0) fermi_point = i
      end do
      if (rank == 0) write (*, *) this%en%fermi, fermi_point
      ! Calculate the Green Function at the Fermi energy plus eta
      do i = 1, 64
         eta = (0.0_rp, 0.0_rp)
         g0_ef = (0.0_rp, 0.0_rp)
         res = (1 - x(i))/x(i)
         eta = cmplx(0.0_rp, res)

         select case (this%control%recur)
         case ('block')
            call this%green%block_green_eta(eta, fermi_point, g0_ef)
         case ('chebyshev')
            call this%green%chebyshev_green_eta(eta, fermi_point, g0_ef)
         end select

         y(i, :, :, :) = real(g0_ef(:, :, :))
         y(i, :, :, :) = (y(i, :, :, :)*w(i))/(x(i)*x(i))
      end do

      occ(:, :) = 0.0d0
      !do m=1,this%lattice%nrec
      do m_glob = start_atom, end_atom
         m = g2l_map(m_glob)
         do j = 1, 64
            do i = 1, 18
               occ(m_glob, i) = occ(m_glob, i) + (y(j, i, i, m))/pi
            end do
         end do
         occ(m_glob, :) = occ(m_glob, :) + 0.5_rp
      end do

      ! Transfer calculated occupations across MPI
#ifdef USE_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, occ, product(shape(occ)), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      sumocc = 0.0d0
      do m = 1, this%lattice%nrec
         if (rank == 0) call g_logger%info('Spin moment of atom'//fmt('i4', m)//' is '//fmt('f10.6', sum(occ(m, 1:9)) - sum(occ(m, 10:18))), __FILE__, __LINE__)
         if (rank == 0) call g_logger%info('Total charge for atom'//fmt('i4', m)//' is '// &
                                           'total= '//fmt('f10.6', sum(occ(m, 1:18)))// &
                                           ' s= '//fmt('f10.6', occ(m, 1) + occ(m, 10))// &
                                           ' p= '//fmt('f10.6', sum(occ(m, 2:4)) + sum(occ(m, 11:13)))// &
                                           ' d= '//fmt('f10.6', sum(occ(m, 5:9)) + sum(occ(m, 14:18))), __FILE__, __LINE__)
         sumocc = sumocc + sum(occ(m, :))
      end do
      if (rank == 0) call g_logger%info('Total number of electrons is '//fmt('f16.6', sumocc), __FILE__, __LINE__)
   end subroutine calculate_moments_gauss_legendre

   subroutine calculate_occupation_gauss_legendre(this, fermi_energy, sumocc_out)
      class(bands) :: this
      ! Local variables
      real(rp), intent(in)  :: fermi_energy
      real(rp), intent(out) :: sumocc_out
      integer :: i, j, l, m, o ! Orbital index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index
      integer :: fermi_point
      complex(rp), dimension(18, 18, this%lattice%nrec) :: g0_ef
      complex(rp) :: eta
      real(rp) :: sgef, pmef, smef, isgn, sumocc
      real(rp), dimension(64, 18, 18, this%lattice%nrec) :: y
      real(rp), dimension(64) :: x, w
      real(rp) :: res, t
      real(rp), dimension(this%lattice%nrec, 18) :: occ

      ! Find the Gauss Legendre roots and weights
      call gauss_legendre(64, 0.0_rp, 1.0_rp, x, w)

      do i = 1, this%en%channels_ldos + 10
         if ((this%en%ene(i) - fermi_energy) .le. 0.001d0) fermi_point = i
      end do
      write (*, *) fermi_energy, fermi_point
      ! Calculate the Green Function at the Fermi energy plus eta
      do i = 1, 64
         eta = (0.0_rp, 0.0_rp)
         g0_ef = (0.0_rp, 0.0_rp)
         res = (1 - x(i))/x(i)
         eta = cmplx(0.0_rp, res)
         call this%green%block_green_eta(eta, fermi_point, g0_ef)
         y(i, :, :, :) = real(g0_ef(:, :, :))
         y(i, :, :, :) = (y(i, :, :, :)*w(i))/(x(i)*x(i))
      end do

      occ(:, :) = 0.0d0
      do m = 1, this%lattice%nrec
         do j = 1, 64
            do i = 1, 18
               occ(m, i) = occ(m, i) + (y(j, i, i, m))/pi
            end do
         end do
      end do

      occ(:, :) = occ(:, :) + 0.5_rp
      sumocc_out = 0.0d0
      do m = 1, this%lattice%nrec
         sumocc_out = sumocc_out + sum(occ(m, :))
      end do
      write (*, *) 'Total electrons:', sumocc_out
   end subroutine calculate_occupation_gauss_legendre

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the Fermi energy using the bisction method
   !> Uses the Gauss-Legendre implementation
   !---------------------------------------------------------------------------
   subroutine calculate_fermi_gauss(this, desired_electrons, tol, max_iter, fermi_energy)
      class(bands) :: this
      real(rp), intent(in)  :: desired_electrons, tol
      integer, intent(in)   :: max_iter
      real(rp), intent(out) :: fermi_energy
      real(rp) :: E_low, E_high, E_mid, error_low, error_high, error_mid
      integer :: iter
      real(rp) :: sumocc

      ! Initial guesses for Fermi energy
      E_low = -0.05 ! Lower bound guess
      E_high = 0.05 ! Upper bound guess

      ! Calculate the error for the initial guesses
      sumocc = 0.0d0
      call this%calculate_occupation_gauss_legendre(E_low, sumocc)
      error_low = desired_electrons - sumocc

      sumocc = 0.0d0
      call this%calculate_occupation_gauss_legendre(E_high, sumocc)
      error_high = desired_electrons - sumocc

      ! Ensure that the errors have different signs
      if (error_low*error_high >= 0.0) then
         print *, "Error: Initial guesses for Fermi energy do not bracket the root."
         return
      end if

      iter = 0
      do while (abs(E_high - E_low) > tol .and. iter < max_iter)
         iter = iter + 1
         ! Calculate mid-point
         E_mid = 0.5*(E_low + E_high)

         ! Evaluate the error at the mid-point
         sumocc = 0.0d0
         call this%calculate_occupation_gauss_legendre(E_mid, sumocc)
         error_mid = desired_electrons - sumocc

         ! Decide which interval to continue with
         if (error_mid*error_low < 0.0) then
            E_high = E_mid
            error_high = error_mid
         else
            E_low = E_mid
            error_low = error_mid
         end if
         this%en%fermi = E_mid
         call this%en%e_mesh()
      end do

      ! Return the Fermi energy
      fermi_energy = E_mid

   end subroutine calculate_fermi_gauss

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the band moments using the Chebyshev-Gauss quadrature method
   !---------------------------------------------------------------------------
   subroutine calculate_moments_chebgauss(this)
      class(bands), intent(inout) :: this
      ! Local variables
      real(rp), dimension(this%control%lld) :: x_i, w_i, q_i
      real(rp), dimension(this%control%lld) :: w, wscale
      real(rp), dimension(this%control%lld, 0:(this%control%lld)) :: polycheb
      real(rp), dimension(18, this%control%lld, this%lattice%ntype) :: doscheb_i
      real(rp), dimension(18, this%lattice%ntype) :: q0l, q1l, q2l
      real(rp) :: a, b, mom0, mom1, mom2, occ
      integer :: i, l, llplusone, n, ll

      ll = this%control%lld
      x_i(:) = 0.0d0
      w_i(:) = 0.0d0
      q_i(:) = 0.0d0

      call chebyshev_gauss_quadrature(ll, x_i(:), q_i(:))

      write (*, *) this%en%fermi
      w_i(:) = ((0.5*(x_i(:) + 1))*(this%en%energy_min))
      w_i(:) = w_i(:) + this%en%fermi
      q_i(:) = -((q_i(:)*(this%en%energy_min))/2)
      ! Defining rescaling coeficients
      a = (this%en%energy_max - this%en%energy_min)/(2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min)/2

      wscale(:) = ((w_i(:) - b)/a)
      ! Calculate the Chebyshev polynomials
      call t_polynomial(ll, ll, wscale(:), polycheb)

      doscheb_i(:, :, :) = 0.0d0

      occ = 0.0d0
      do n = 1, this%lattice%nrec
         ! Calculate the density of states
         do l = 1, 18
            do i = 1, ll
               doscheb_i(l, :, n) = doscheb_i(l, :, n) + real(this%recursion%mu_ng(l, l, i, n))*polycheb(:, i - 1)
            end do
         end do
         do l = 1, 18
            doscheb_i(l, :, n) = doscheb_i(l, :, n)/((sqrt((a**2) - ((w_i(:) - b)**2)))*pi)
         end do

         do l = 1, 18
            do i = 1, ll
               if (isnan(doscheb_i(l, i, n))) doscheb_i(l, i, n) = 0.0d0
            end do
         end do

         do l = 1, 18
            q0l(n, l) = sum(doscheb_i(l, :, n)*q_i(:))
            q1l(n, l) = sum(doscheb_i(l, :, n)*w_i(:)*q_i(:))
            q2l(n, l) = sum(doscheb_i(l, :, n)*(w_i(:)**2)*q_i(:))
            occ = occ + q0l(n, l)
         end do
         write (*, *) 'Charge of atom:', n, sum(q0l(n, :))
      end do
      write (*, *) 'Total number of electrons:', occ
   end subroutine calculate_moments_chebgauss

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the magnetic moments mx, my, mz and mtot
   !---------------------------------------------------------------------------
   subroutine calculate_magnetic_moments(this)
      use mpi_mod
      class(bands) :: this
      ! Local variables
      real(rp) :: mx, my, mz
      integer :: i, j ! Orbital index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index

      integer :: na_loc

      call this%calculate_projected_dos()

      do na = start_atom, end_atom
         na_loc = g2l_map(na)
         call simpson_m(this%symbolic_atom(this%lattice%nbulk + na)%potential%mx, this%en%edel, this%en%fermi, this%nv1, this%dx(:, na_loc), this%e1, 0, this%en%ene)
         call simpson_m(this%symbolic_atom(this%lattice%nbulk + na)%potential%my, this%en%edel, this%en%fermi, this%nv1, this%dy(:, na_loc), this%e1, 0, this%en%ene)
         call simpson_m(this%symbolic_atom(this%lattice%nbulk + na)%potential%mz, this%en%edel, this%en%fermi, this%nv1, this%dz(:, na_loc), this%e1, 0, this%en%ene)

         this%symbolic_atom(this%lattice%nbulk + na)%potential%mom0(1) = this%symbolic_atom(this%lattice%nbulk + na)%potential%mx
         this%symbolic_atom(this%lattice%nbulk + na)%potential%mom0(2) = this%symbolic_atom(this%lattice%nbulk + na)%potential%my
         this%symbolic_atom(this%lattice%nbulk + na)%potential%mom0(3) = this%symbolic_atom(this%lattice%nbulk + na)%potential%mz

         call simpson_m(this%symbolic_atom(this%lattice%nbulk + na)%potential%mom1(1), this%en%edel, this%en%fermi, this%nv1, this%dx(:, na_loc), this%e1, 1, this%en%ene)
         call simpson_m(this%symbolic_atom(this%lattice%nbulk + na)%potential%mom1(2), this%en%edel, this%en%fermi, this%nv1, this%dy(:, na_loc), this%e1, 1, this%en%ene)
         call simpson_m(this%symbolic_atom(this%lattice%nbulk + na)%potential%mom1(3), this%en%edel, this%en%fermi, this%nv1, this%dz(:, na_loc), this%e1, 1, this%en%ene)

         this%symbolic_atom(this%lattice%nbulk + na)%potential%mtot = sqrt((this%symbolic_atom(this%lattice%nbulk + na)%potential%mx**2) + &
                                                                           (this%symbolic_atom(this%lattice%nbulk + na)%potential%my**2) + &
                                                                           (this%symbolic_atom(this%lattice%nbulk + na)%potential%mz**2)) + 1.0d-15

         this%symbolic_atom(this%lattice%nbulk + na)%potential%mom(1) = this%symbolic_atom(this%lattice%nbulk + na)%potential%mx/this%symbolic_atom(this%lattice%nbulk + na)%potential%mtot
         this%symbolic_atom(this%lattice%nbulk + na)%potential%mom(2) = this%symbolic_atom(this%lattice%nbulk + na)%potential%my/this%symbolic_atom(this%lattice%nbulk + na)%potential%mtot
         this%symbolic_atom(this%lattice%nbulk + na)%potential%mom(3) = this%symbolic_atom(this%lattice%nbulk + na)%potential%mz/this%symbolic_atom(this%lattice%nbulk + na)%potential%mtot

         !if(rank==0) call g_logger%info(´Spin moment of atom´//fmt(´i4´,na)//´ is ´//fmt(´f10.6´,this%symbolic_atom(this%lattice%nbulk+na)%potential%mtot),__FILE__,__LINE__)
         call g_logger%info('Spin moment of atom'//fmt('i4', na)//' is '//fmt('f10.6', this%symbolic_atom(this%lattice%nbulk + na)%potential%mtot), __FILE__, __LINE__)
         mx = this%symbolic_atom(this%lattice%nbulk + na)%potential%mx
         my = this%symbolic_atom(this%lattice%nbulk + na)%potential%my
         mz = this%symbolic_atom(this%lattice%nbulk + na)%potential%mz
         !if(rank==0) call g_logger%info(´Spin moment projections of atom´//fmt(´i4´,na)//´ is ´//fmt(´f10.6´,mx)//´ ´//fmt(´f10.6´,my)//´ ´//fmt(´f10.6´,mz),__FILE__,__LINE__)

         if (this%recursion%hamiltonian%local_axis) then
            call g_logger%info('Local spin moment projections of atom'//fmt('i4', na)//' is '//fmt('f10.6', mx)//' '//fmt('f10.6', my)//' '//fmt('f10.6', mz), __FILE__, __LINE__)
         else
            call g_logger%info('Spin moment projections of atom'//fmt('i4', na)//' is '//fmt('f10.6', mx)//' '//fmt('f10.6', my)//' '//fmt('f10.6', mz), __FILE__, __LINE__)
         end if

         if (this%control%nsp < 3) this%symbolic_atom(this%lattice%nbulk + na)%potential%mom(:) = [0.0d0, 0.0d0, 1.00d0]

      end do
   end subroutine calculate_magnetic_moments

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the orbital moments lx, ly, lz and ltot
   !---------------------------------------------------------------------------
   subroutine calculate_orbital_moments(this)
      use mpi_mod
      class(bands) :: this
      ! Local variables
      real(rp) :: lx, ly, lz
      integer :: i, j ! Orbital index
      integer :: mdir ! Magnetic index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index
      real(rp), dimension(9, 9, 3) :: l_orb
      complex(rp), dimension(9, 9) :: mLx, mLy, mLz
      !

      integer :: na_loc

      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
      mLx(:, :) = L_x(:, :)
      mLy(:, :) = L_y(:, :)
      mLz(:, :) = L_z(:, :)

      ! Transforming them into the spherical harmonics coordinates
      call hcpx(mLx, 'cart2sph')
      call hcpx(mLy, 'cart2sph')
      call hcpx(mLz, 'cart2sph')

      call this%calculate_orbital_dos()

      !do na=1, this%lattice%nrec
      do na = start_atom, end_atom
         na_loc = g2l_map(na)

         do mdir = 1, 3
            do i = 1, 9
               do j = 1, 9
                  call simpson_m(l_orb(i, j, mdir), this%en%edel, this%en%fermi, this%nv1, this%d_orb(j, i, mdir, :, na_loc), this%e1, 0, this%en%ene)
               end do
            end do
         end do
         lx = 0.0_rp !-0.5_rp * rtrace9(matmul(mLx,l_orb(:,:,1)))
         ly = 0.0_rp !-0.5_rp * rtrace9(matmul(mLy,l_orb(:,:,2)))
         lz = -0.5_rp*rtrace9(matmul(mLz, l_orb(:, :, 3)))
         call g_logger%info('Orbital moment of atom'//fmt('i4', na)//' is '//fmt('f10.6', lz), __FILE__, __LINE__)
         this%symbolic_atom(this%lattice%nbulk + na)%potential%lmom(1) = lx
         this%symbolic_atom(this%lattice%nbulk + na)%potential%lmom(2) = ly
         this%symbolic_atom(this%lattice%nbulk + na)%potential%lmom(3) = lz

      end do
   end subroutine calculate_orbital_moments

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Builds the effective single-particle potentials for the LDA+U+J correction for s,p and d orbitals.
   !> (Only the slater integrals F0, F2, F4 and F6 are defined. Machinery works for f orbitals
   !> too, but the corresponding (32x32) basis needs implementation.)
   !> Created and improved by Viktor Frilén and Emil Beiersdorf 17.07.2024
   !---------------------------------------------------------------------------
   subroutine Hubbard_U_Potential(this)
      class(bands) :: this

      ! Local variables
      integer :: i, j ! Orbital indices
      integer :: na ! Atom index
      integer :: ie ! Energy channel index
      integer :: ispin ! Spin index
      integer :: l ! Orbital number, 0,1,2,3 (1,2,3,4) = s,p,d,f (index postions in arrays)
      integer :: cntr ! Counts number of orbitals with a U per atom.  
      real(rp), dimension(:,:), allocatable :: hub_u, hub_j ! Hubbard U and J parameters
      real(rp) :: result
      integer :: l_index ! Index for l 1,2,3,4 = s,p,d,f
      integer :: m, m1, m2, m3, m4, m_max, m1_val, m2_val, m3_val, m4_val ! Magnetic quantum numbers
      real(rp) :: f0, f2, f4, f6 ! Slater integrals
      real(rp), dimension(this%lattice%nrec,4,4) :: f ! Slater integrals
      real(rp), dimension(18, 18, this%en%channels_ldos + 10, this%lattice%nrec) :: im_g0
      ! Local variable array for local density matrix
      real(rp), dimension(:,:,:,:,:), allocatable :: LDM ! Local density matrix (LDM), works for spdf orbitals
      real(rp) :: dc ! Double counting term
      real(rp), dimension(4) :: U_energy, dc_energy 
      real(rp), dimension(this%lattice%nrec, 4, 2) :: n_spin ! LDM with m traced out
      real(rp), dimension(this%lattice%nrec, 4) :: n_tot  ! LDM with traced out spin  
      ! Temporary potential that will be put into this%hubbard_u_pot
      real(rp), dimension(this%lattice%nrec, 4, 2, 7, 7) :: hub_pot    
      
      type :: ArrayType
         integer, allocatable :: val(:)
      end type ArrayType
      type(ArrayType), dimension(4) :: ms
      type(ArrayType), dimension(this%lattice%nrec) :: l_arr
  
      ms(1)%val = [0]
      ms(2)%val = [-1, 0, 1]
      ms(3)%val = [-2, -1, 0, 1, 2]
      ms(4)%val = [-3, -2, -1, 0, 1, 2, 3]

      ! Checks which orbital to add U onto. Only works for spd-orbitals
      if (this%recursion%hamiltonian%hubbardU_check) then
         allocate(LDM(this%lattice%nrec, 4, 2, 7, 7)) ! (number of atoms, number of orbitals (spdf), m-value indexing for spdf)
         allocate(hub_u(this%lattice%nrec, 4))
         allocate(hub_j(this%lattice%nrec, 4))
         LDM(:,:,:,:,:) = 0.0d0
         im_g0(:,:,:,:) = 0.0d0
         this%recursion%hamiltonian%hubbard_u_pot = 0.0d0
         hub_pot = 0.0d0
         n_tot(:,:) = 0.0d0
         n_spin(:,:,:) = 0.0d0
         hub_u(:,:) = 0.0d0
         hub_j(:,:) = 0.0d0
         f = this%recursion%hamiltonian%f
         hub_u = this%recursion%hamiltonian%hub_u_sort 
         hub_j = this%recursion%hamiltonian%hub_j_sort

         
         !> Creates an array with each orbital for each atom
         do i = 1, this%lattice%nrec
            cntr = count(hub_u(i,:) > 1.0E-10) ! Counts orbitals with Hub U for each atom for allocation purposes
            allocate(l_arr(i)%val(cntr))
         end do

         do i = 1, this%lattice%nrec
            cntr = 0
            do j = 1, 4
               if (hub_u(i,j) > 1.0E-10) then
                  cntr = cntr + 1
                  l_arr(i)%val(cntr) = j ! Fills the l_arr list with the orbitals that have Hub U
               end if
            end do
         end do       
               
         !> Sets up the imaginary Green's function (only for spd orbitals)
         do na = 1, this%lattice%nrec
            do i = 1, 18
               do j = 1, 18
                  do ie  = 1, this%en%channels_ldos + 10
                     im_g0(i,j,ie,na) = im_g0(i,j,ie,na) - aimag(this%green%g0(i,j,ie,na))/pi
                  end do
               end do
            end do
         end do

         !> Calculates the local density matrix by integrating the Green's function over the energy channels.
         !> Only implemented for spd orbitals.
         do na = 1, this%lattice%nrec
            do l = 0, 2
               do ispin = 1, 2
                  do i = 1, 2*l + 1 !m3
                     do j = 1, 2*l + 1 !m4
                        call simpson_m(result, this%en%edel, this%en%fermi, this%nv1, im_g0(l**2 + i + (ispin-1)*9, l**2 + j + (ispin-1)*9, :, na), this%e1, 0, this%en%ene)
                        LDM(na, l + 1, ispin, i, j) = result
                     end do
                  end do
               end do
            end do
         end do
         
         !> Builds the Hubbard U+J potential 
         print *, ''
         print *, '-----------------------------------------------------------------------------------------------------------------'
         print *, 'Calculate Hubbard U potential for:'
         do na = 1, this%lattice%nrec
            print *, ' Atom ', na
            print *, ''
            !> Calculates traces of the local density matrix, n_spin is the trace in m, n_tot is trace in spin of n_spin.
            do l = 1, 4 !0 to 3 but indexing starts on 1
               do ispin = 1, 2
                  n_spin(na, l, ispin) = trace(LDM(na, l, ispin, :, :)) 
                  n_tot(na, l) = n_tot(na, l) + n_spin(na, l, ispin)
               end do
            end do
            U_energy = 0.0d0

            !> l_arr(na)%val(l) gives the orbital index for each atom, for s,p,d its value is 1,2,3
            !> l_arr(na)%val(l)-1 gives the actual l-value for each atom, for s,p,d its value is 0,1,2   
            do l = 1, size(l_arr(na)%val)
               l_index = l_arr(na)%val(l)
               m_max = 2*l_index-1
               do ispin = 1, 2
                  do m1 = 1, m_max
                     do m2 = 1, m_max
                        do m3 = 1, m_max
                           do m4 = 1, m_max
                              m1_val = ms(l_index)%val(m1)
                              m2_val = ms(l_index)%val(m2)
                              m3_val = ms(l_index)%val(m3)
                              m4_val = ms(l_index)%val(m4)
                              f0 = f(na,l_index,1)
                              f2 = f(na,l_index,2)
                              f4 = f(na,l_index,3)
                              f6 = f(na,l_index,4)
                              hub_pot(na, l_index, ispin, m1, m2) = hub_pot(na, l_index, ispin, m1, m2) &
                              + Coulomb_mat(l_index-1,m1_val,m3_val,m2_val,m4_val,f0,f2,f4,f6)*LDM(na,l_index,3-ispin,m3,m4) &
                              + (Coulomb_mat(l_index-1,m1_val,m3_val,m2_val,m4_val,f0,f2,f4,f6) &
                              - Coulomb_mat(l_index-1,m1_val,m3_val,m4_val,m2_val,f0,f2,f4,f6))*LDM(na,l_index,ispin,m3,m4)
                              U_energy(l_index) = U_energy(l_index) + 0.5_rp *( Coulomb_mat(l_index-1,m1_val,m3_val,m2_val,m4_val,f0,f2,f4,f6) &
                              * LDM(na,l_index,ispin,m1,m2) * LDM(na,l_index,3-ispin,m3,m4) &
                              + ( Coulomb_mat(l_index-1,m1_val,m3_val,m2_val,m4_val,f0,f2,f4,f6) &
                              - Coulomb_mat(l_index-1,m1_val,m3_val,m4_val,m2_val,f0,f2,f4,f6) ) &
                              * LDM(na,l_index,ispin,m1,m2) * LDM(na,l_index,ispin,m3,m4) )
                           end do
                        end do
                        ! Double counting terms only added to V_mm'^spin diagonal in m and m' (local density matrix is on-site)
                        if (m1 == m2) then
                           hub_pot(na, l_index, ispin, m1, m2) = hub_pot(na, l_index, ispin, m1, m2) &
                           - hub_u(na,l_index)*(n_tot(na,l_index) - 0.5_rp) + hub_j(na,l_index)*(n_spin(na,l_index,ispin) - 0.5_rp)
                        end if
                     end do
                  end do
               end do
            end do  

         !> Puts hub_pot into global hubbard_u_pot (only done for spd-orbitals)
            do l = 0, 2
               do i = 1, 2*l + 1
                  do j = 1, 2*l + 1
                     this%recursion%hamiltonian%hubbard_u_pot(l**2+i, l**2+j, na) = hub_pot(na, l+1, 1, i, j)
                     this%recursion%hamiltonian%hubbard_u_pot(l**2+i+9, l**2+j+9, na) = hub_pot(na, l+1, 2, i, j)
                  end do
               end do
               dc_energy(l+1) = 0.5_rp*(hub_u(na,l+1)*n_tot(na,l+1)*(n_tot(na,l+1) - 1.0_rp) &
               - hub_j(na,l+1)*(n_spin(na,l+1,1)*(n_spin(na,l+1,1) - 1.0_rp) + n_spin(na,l+1,2)*(n_spin(na,l+1,2) - 1.0_rp)))

               do i = 1, size(l_arr(na)%val)
                  if (l+1 == l_arr(na)%val(i)) then ! Only prints U_energy if that orbital has a specified U
                     print *, ' Orbital ', this%recursion%hamiltonian%orb_conv(l+1), ' :'
                     print *, '  (U, dc, total) [eV] = ', U_energy(l+1)*ry2ev, dc_energy(l+1)*ry2ev, U_energy(l+1)*ry2ev - dc_energy(l+1)*ry2ev
                     print *, ''
                  end if
               end do
            end do
            
         end do 
         print *, '-----------------------------------------------------------------------------------------------------------------'
         print *,''

      if (allocated(LDM)) deallocate(LDM)
      if (allocated(hub_u)) deallocate(hub_u)
      if (allocated(hub_j)) deallocate(hub_j)
      end if
      
   end subroutine Hubbard_U_Potential


   subroutine Hubbard_V_Potential(this)
      class(bands) :: this

   ! Local variables
      integer :: na
      integer :: ie
      integer :: nn 
      integer :: i, j, k, l1, l2, lmax, s1, s2, m1_max, m2_max, l1_ref, l2_ref, m1, m2, ia, nr
      integer :: ij, kl ! Pair index
      integer :: atom1, atom2 ! Atom index
      integer :: atom1_type, atom2_type ! Atom type
      integer :: tot_l_pair ! Total number of orbital pairs to be considered for an atom pair
      real(rp) :: result, result2, nn_dist, dist
      real(rp) :: vij ! Hubbard V correction between atom+orbital i and j
      logical :: calc_pair
      !this variable is calculated in lattice (1 = nearest neighbour)
      integer, dimension(this%lattice%njij, 2) :: ij_pair ! Change this later
      complex(rp), dimension(18,18,this%en%channels_ldos + 10,this%lattice%njij) :: gij, gji ! intersite greens, calculated elsewhere
      real(rp), dimension(18,18,this%en%channels_ldos + 10,this%lattice%njij) :: im_gij, im_gji ! -Im(gij)/2 (temp variable)
      real(rp), dimension(18,18,this%lattice%njij) :: nij, nji ! integrated im_gij
      real(rp), dimension(18,18) :: nji_temp
      real(rp), dimension(18,18,this%lattice%nrec) :: DM ! density matrix
      real(rp), dimension(:,:,:,:), allocatable :: hub_V_pot ! +V correction potential. nn should be changed to an appropriate value
      
      ! Restructuring of density matrix from nij(lms,l'm's',ij) to nij(ij,s,s',l,l',m,m')
      type :: ArrayType
         integer, allocatable :: mn(:,:)
      end type ArrayType

      type(ArrayType), dimension(2,2,3,3) :: nji_sorted


      ! nn = size(this%lattice%ijpair_sorted,3)
      nn = size(this%recursion%hamiltonian%ee, 3)
      na = this%lattice%nrec
      if (allocated(hub_V_pot)) deallocate (hub_V_pot)
      allocate(hub_V_pot(18,18,nn,this%lattice%nrec))
      im_gij = 0.0d0
      im_gji = 0.0d0
      hub_V_pot(:,:,:,:) = 0.0d0
      

      !> Calculates -im(gij)/pi
      do ij = 1, this%lattice%njij
         do i = 1, 18
            do j = 1, 18
               do ie  = 1, this%en%channels_ldos + 10
                  im_gij(i,j,ie,ij) = im_gij(i,j,ie,ij) - aimag(this%green%gij(i,j,ie,ij))/pi
                  im_gji(i,j,ie,ij) = im_gji(i,j,ie,ij) - aimag(this%green%gji(i,j,ie,ij))/pi
               end do
               call simpson_m(result, this%en%edel, this%en%fermi, this%nv1, im_gij(i,j,:,ij), this%e1, 0, this%en%ene)
               call simpson_m(result2, this%en%edel, this%en%fermi, this%nv1, im_gji(i,j,:,ij), this%e1, 0, this%en%ene)
               nij(i,j,ij) = result
               nji(i,j,ij) = result2
            end do
         end do
      end do


      print *, ''
      print *, '-----------------------------------------------------------------------------------------------------------------'
      do na = 1, this%lattice%nrec
         !> Loop over nearest neighbours
         ia = this%lattice%atlist(na) ! Atom number in clust
         nr = this%lattice%nn(ia, 1) ! Number of neighbours considered
         nn_dist = this%lattice%nn_dist(na) ! Distance of nearest neighbours for atom na
         do ij = 2, nr
            atom1 = ia
            atom1_type = this%lattice%iz(atom1)
            atom2 = this%lattice%nn(atom1, ij)
            atom2_type = this%lattice%iz(atom2)
            !> Checks if atom1 and atom2 are nearest neighbours.
            dist = sqrt((this%lattice%cr(1,atom1) - this%lattice%cr(1,atom2))**2 + (this%lattice%cr(2,atom1) &
                  - this%lattice%cr(2,atom2))**2 + (this%lattice%cr(3,atom1) - this%lattice%cr(3,atom2))**2)
            if ( dist - 0.001d0 > nn_dist ) then
               print *, 'Atom ', atom1, 'and atom ', atom2, 'are not nearest neighbours. Do not add a +V correction.'
            else
               print *, 'Atom ', atom1, 'and atom ', atom2, 'are nearest neighbours. Check to see if +V correction should be added.'
               !> Set correct density matrix for the given pair
               nji_temp = 0.0d0
               calc_pair = .false.
               do kl = 1, this%lattice%njij
                  if (atom1_type == this%lattice%iz(this%lattice%ijpair(kl,1)) .and. atom2_type == this%lattice%iz(this%lattice%ijpair(kl,2))) then
                     nji_temp(:,:) = TRANSPOSE(nji(:,:,kl))
                     calc_pair = .true.
                     exit
                  else if (atom1_type == this%lattice%iz(this%lattice%ijpair(kl,2)) .and. atom2_type == this%lattice%iz(this%lattice%ijpair(kl,1))) then
                     nji_temp(:,:) = TRANSPOSE(nij(:,:,kl))
                     calc_pair = .true.
                     exit
                  end if
               end do
               !> If there is a +V correction for this pair
               if (calc_pair) then
                  print *, 'Adding +V correction.'
                  print *, ''
                  !> Check between which orbitals there is a +V correction
                  !> Number of orbital pairs between the two atoms
                  tot_l_pair = this%recursion%hamiltonian%orbs_v_num(atom1_type, atom2_type)
                  do kl = 1, tot_l_pair
                     l1_ref = this%recursion%hamiltonian%orbs_v(atom1_type,atom2_type)%val(kl,1)
                     l2_ref = this%recursion%hamiltonian%orbs_v(atom1_type,atom2_type)%val(kl,2)
                     vij = this%recursion%hamiltonian%hubbard_v(atom1_type, atom2_type, l1_ref, l2_ref)
                     do s1 = 1, 2
                        do l1 = 1, 3
                           do l2 = 1, 3
                              if ((l1 == l1_ref) .and. (l2 == l2_ref)) then
                                 do m1 = 1, l1*2-1
                                    do m2 = 1, l2*2-1
                                       hub_V_pot((l1-1)**2 + (s1-1)*9 + m1, (l2-1)**2 + (s1-1)*9 + m2, ij, na) = &
                                                   -vij*nji_temp((l1-1)**2 + (s1-1)*9 + m1, (l2-1)**2 + (s1-1)*9 + m2)
                                    end do
                                 end do
                              end if
                           end do
                        end do
                     end do
                  end do
               else
                  print *, 'No Hubbard V correction for atom pair', atom1_type, atom2_type
               end if
            end if
         end do
      end do
      print *, ''
      print *, '-----------------------------------------------------------------------------------------------------------------'
      print *, ''

      this%recursion%hamiltonian%hubbard_v_pot = hub_V_pot

   end subroutine Hubbard_V_Potential

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the Hubbard U and J parameters and creates an effective U = U - J
   !> Based on the article
   !> Luis A. Agapito, Stefano Curtarolo, and Marco Buongiorno Nardelli - Phys. Rev. X 5, 011006 – Published 28 January 2015
   !> Created and improved by Viktor Frilén 16.08.2024
   !---------------------------------------------------------------------------
   subroutine calc_hubbard_U(this)
      class(bands) :: this

      ! Local variables
      integer :: i, j ! Orbital indices
      integer :: na ! Atom index
      integer :: ie ! Energy channel index
      integer :: ispin, ispin2 ! Spin index
      integer :: l ! Orbital number, 0,1,2,3 (1,2,3,4) = s,p,d,f (index postions in arrays)
      integer :: l_index ! Index for l 1,2,3,4 = s,p,d,f
      integer :: m, m1, m2, m3, m4, m_max, m1_val, m2_val, m3_val, m4_val ! Magnetic quantum numbers
      real(rp) :: f0, f2, f4, f6 ! Slater integrals
      real(rp), dimension(this%lattice%nrec, 4) :: hub_u_temp, hub_j_temp, hub_u_eff ! Variables used for calculating U, J and U_eff ab initio
      real(rp) :: hub_u_denominator, hub_j_denominator ! Variables to store the denominators of U and J
      real(rp), dimension(:,:), allocatable :: hub_u, hub_j ! Hubbard U and J parameters
      real(rp) :: result, temp, temp2
      real(rp), dimension(this%lattice%nrec,4,4) :: f ! Slater integrals
      real(rp), dimension(18, 18, this%en%channels_ldos + 10, this%lattice%nrec) :: im_g0
      ! Local variable array for local density matrix
      real(rp), dimension(:,:,:,:,:), allocatable :: LDM ! Local density matrix (LDM), works for spdf orbitals
      real(rp), dimension(:,:,:,:,:), allocatable :: LDM_screened ! Local density matrix renormalized to include screening according to L. A. Agapito
      real(rp) :: dc ! Double counting term
      real(rp), dimension(4) :: U_energy, dc_energy 
      real(rp), dimension(this%lattice%nrec, 4, 2) :: n_spin ! LDM with m traced out
      real(rp), dimension(this%lattice%nrec, 4) :: n_tot  ! LDM with traced out spin  
      ! Temporary potential that will be put into this%hubbard_u_pot in hamiltonian
      real(rp), dimension(this%lattice%nrec, 4, 2, 7, 7) :: hub_pot  
      integer :: cntr ! Counts number of orbitals with a U per atom.    
      
      type :: ArrayType
         integer, allocatable :: val(:)
      end type ArrayType
      type(ArrayType), dimension(4) :: ms
      type(ArrayType), dimension(this%lattice%nrec) :: l_arr
  
      ms(1)%val = [0]
      ms(2)%val = [-1, 0, 1]
      ms(3)%val = [-2, -1, 0, 1, 2]
      ms(4)%val = [-3, -2, -1, 0, 1, 2, 3]

      !> Sets up the DFT+U calculation usually done in hamiltonian.f90
      !> Only needed for the first iteration but this does not cause problems
      this%recursion%hamiltonian%hubbardU_check = .true.
      if (allocated(this%recursion%hamiltonian%F)) deallocate (this%recursion%hamiltonian%F)
      allocate(this%recursion%hamiltonian%F(this%lattice%nrec, 4, 4))
      this%recursion%hamiltonian%F = 0.0d0
      if (allocated(this%recursion%hamiltonian%hub_u_sort)) deallocate (this%recursion%hamiltonian%hub_u_sort)
      allocate(this%recursion%hamiltonian%hub_u_sort(this%lattice%nrec, 4))
      this%recursion%hamiltonian%hub_u_sort = 0.0d0
      if (allocated(this%recursion%hamiltonian%hub_j_sort)) deallocate (this%recursion%hamiltonian%hub_j_sort)
      allocate(this%recursion%hamiltonian%hub_j_sort(this%lattice%nrec, 4))
      this%recursion%hamiltonian%hub_j_sort = 0.0d0
      if (allocated(this%recursion%hamiltonian%hubbard_u_pot)) deallocate (this%recursion%hamiltonian%hubbard_u_pot)
      allocate(this%recursion%hamiltonian%hubbard_u_pot(18,18,this%lattice%nrec))
      this%recursion%hamiltonian%hubbard_u_pot = 0.0d0

      !> Reads of which atoms and orbitals should have +U correction. 
      !> Sets an arbitrary value to %hub_u_sort for the algorithm copied from Hubbard_U_Potential to work
      do na = 1, this%lattice%nrec
         do l = 1, 4
            if ( this%recursion%hamiltonian%hubbard_u_sc(na,l) == 1 ) then
               this%recursion%hamiltonian%hub_u_sort(na,l) = 1.0_rp
            end if
         end do
      end do

      !> Checks which orbital to add U onto. Only works for spd-orbitals
      allocate(LDM(this%lattice%nrec, 4, 2, 7, 7)) ! (number of atoms, number of orbitals (spdf), m-value indexing for spdf)
      allocate(LDM_screened(this%lattice%nrec, 4, 2, 7, 7))
      allocate(hub_u(this%lattice%nrec, 4))
      allocate(hub_j(this%lattice%nrec, 4))
      LDM(:,:,:,:,:) = 0.0d0
      LDM_screened(:,:,:,:,:) = 0.0d0
      im_g0(:,:,:,:) = 0.0d0
      n_tot(:,:) = 0.0d0
      n_spin(:,:,:) = 0.0d0
      hub_u(:,:) = 0.0d0
      hub_j(:,:) = 0.0d0
      f = this%recursion%hamiltonian%f
      hub_u = this%recursion%hamiltonian%hub_u_sort
      hub_j = this%recursion%hamiltonian%hub_j_sort
      hub_u_temp = 0.0d0
      hub_j_temp = 0.0d0
      hub_u_eff = 0.0d0
      
      !> Creates an array with each orbital for each atom
      do i = 1, this%lattice%nrec
         cntr = count(hub_u(i,:) > 1.0E-10) ! Counts orbitals with Hub U for each atom for allocation purposes
         allocate(l_arr(i)%val(cntr))
      end do

      do i = 1, this%lattice%nrec
         cntr = 0
         do j = 1, 4
            if (hub_u(i,j) > 1.0E-10) then
               cntr = cntr + 1
               l_arr(i)%val(cntr) = j ! Fills the l_arr list with the orbitals that have Hub U
            end if
         end do
      end do       
            
      !> Sets up the imaginary Green's function (only for spd orbitals)
      do na = 1, this%lattice%nrec
         do i = 1, 18
            do j = 1, 18
               do ie  = 1, this%en%channels_ldos + 10
                  im_g0(i,j,ie,na) = im_g0(i,j,ie,na) - aimag(this%green%g0(i,j,ie,na))/pi
               end do
            end do
         end do
      end do

      !> Calculates the local density matrix by integrating the Green's function over the energy channels.
      !> Only implemented for spd orbitals.
      do na = 1, this%lattice%nrec
         do l = 0, 2
            do ispin = 1, 2
               do i = 1, 2*l + 1 !m3
                  do j = 1, 2*l + 1 !m4
                     call simpson_m(result, this%en%edel, this%en%fermi, this%nv1, im_g0(l**2 + i + (ispin-1)*9, l**2 + j + (ispin-1)*9, :, na), this%e1, 0, this%en%ene)
                     LDM(na, l + 1, ispin, i, j) = result
                  end do
               end do
            end do
         end do
      end do

      !> Calculates the renormalized (or screened local density matrix)
      do na = 1, this%lattice%nrec
         do l = 0, 2
            do ispin = 1, 2
               do i = 1, 2*l + 1 !m3
                  do j = 1, 2*l + 1 !m4
                     LDM_screened(na,l+1,ispin,i,j) = LDM(na,l+1,ispin,i,j)*(LDM(na,l+1,ispin,i,i) + LDM(na,l+1,ispin,j,j))
                  end do
               end do
            end do
         end do
      end do
      
      !> calculates U and J
      print *, ''
      print *, 'Calculate Hubbard U parameter for:'
      do na = 1, this%lattice%nrec
         print *, 'Atom ', na
         U_energy = 0.0d0

         do l = 1, size(l_arr(na)%val)
            l_index = l_arr(na)%val(l)
            m_max = 2*l_index-1 
            do ispin = 1, 2  
               do ispin2 = 1, 2
                  do m1 = 1, m_max
                     do m2 = 1, m_max
                        do m3 = 1, m_max
                           do m4 = 1, m_max
                              m1_val = ms(l_index)%val(m1)
                              m2_val = ms(l_index)%val(m2)
                              m3_val = ms(l_index)%val(m3)
                              m4_val = ms(l_index)%val(m4)
                              
                              f0 = tabulated_slater_integrals(1,l_index,l_index,l_index,l_index)
                              f2 = tabulated_slater_integrals(2,l_index,l_index,l_index,l_index)
                              f4 = tabulated_slater_integrals(3,l_index,l_index,l_index,l_index)
                              f6 = tabulated_slater_integrals(4,l_index,l_index,l_index,l_index)
                              
                              ! Add onto Hubbard U
                              hub_u_temp(na, l_index) = hub_u_temp(na, l_index) &
                              + Coulomb_mat(l_index-1,m1_val,m3_val,m2_val,m4_val,f0,f2,f4,f6)*LDM_screened(na,l_index,ispin,m1,m2) &
                              *LDM_screened(na,l_index,ispin2,m3,m4)
                              
                              ! Add onto Hubbard J
                              if (ispin == ispin2) then
                                 hub_j_temp(na, l_index) = hub_j_temp(na, l_index) &
                                 + Coulomb_mat(l_index-1,m1_val,m3_val,m4_val,m2_val,f0,f2,f4,f6)*LDM_screened(na,l_index,ispin,m1,m2) &
                                 *LDM_screened(na,l_index,ispin2,m3,m4)
                              end if
                           end do
                        end do

                     end do
                  end do
               end do
            end do
            ! Calculate denominator of U and J
            hub_u_denominator = 0.0d0
            hub_j_denominator = 0.0d0
            do m1 = 1, m_max
               do m2 = 1, m_max
                  hub_u_denominator = hub_u_denominator + LDM(na,l_index,1,m1,m1)*LDM(na,l_index,2,m2,m2) &
                  + LDM(na,l_index,2,m1,m1)*LDM(na,l_index,1,m2,m2)
                  if (m1 /= m2) then
                     hub_u_denominator = hub_u_denominator + LDM(na,l_index,1,m1,m1)*LDM(na,l_index,1,m2,m2) &
                     + LDM(na,l_index,2,m1,m1)*LDM(na,l_index,2,m2,m2)

                     hub_j_denominator = hub_j_denominator + LDM(na,l_index,1,m1,m1)*LDM(na,l_index,1,m2,m2) &
                     + LDM(na,l_index,2,m1,m1)*LDM(na,l_index,2,m2,m2)

                  end if   
               end do
            end do

            temp = hub_u_temp(na,l_index)/hub_u_denominator
            hub_u_temp(na,l_index) = temp
            ! Check for s-electron since they don't have an exchange parameter
            if (l_index == 1) then
               temp = 0.0d0
            else
               temp = hub_j_temp(na,l_index)/hub_j_denominator
            end if
            hub_j_temp(na,l_index) = temp
            hub_u_eff(na,l_index) = hub_u_temp(na,l_index) - hub_j_temp(na,l_index)

            print *, '----------------------------------------------------------------------------------------------'
            print *, 'Calculated values of U, J and U_eff (= U-J)'
            print *, '----------------------------------------------------------------------------------------------'
            print *, 'Atom ', na, 'Orbital (spd) = (012)', l_index-1
            print *, 'U = ', hub_u_temp(na,l_index)
            print *, 'J = ', hub_j_temp(na,l_index)
            print *, 'U_eff = ', hub_u_eff(na,l_index)
            ! Stores the new value of Hubbard U
            this%recursion%hamiltonian%hub_u_sort(na,l_index) = hub_u_eff(na,l_index)
            this%recursion%hamiltonian%f(na,l_index,1) = hub_u_eff(na,l_index)
         end do
      end do 

      !> Check if calculated Hubbard U has converged
      !> Same criteria as Luis A. Agapito et al - Phys. Rev. X 5, 011006
      j = 0
      do na = 1, this%lattice%nrec
         do l = 1, 4
            if ( abs(this%hubbard_u_eff_old(na,l) - hub_u_eff(na,l)) < 3e-5 ) then
               j = j + 1
            end if
         end do
      end do
      if ( j == this%lattice%nrec*4 ) then
         this%hubbard_u_converged = .true.
      end if

      !> Store Hubbard U for later check
      this%hubbard_u_eff_old(:,:) = hub_u_eff(:,:)


      if (allocated(LDM)) deallocate(LDM)
      if (allocated(LDM_screened)) deallocate(LDM_screened)
      if (allocated(hub_u)) deallocate(hub_u)
      if (allocated(hub_j)) deallocate(hub_j)

   end subroutine calc_hubbard_U



   subroutine calculate_projected_dos(this)
      use mpi_mod
      class(bands) :: this
      ! Local variables
      integer :: i, j ! Orbital index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index

      integer :: na_glob

      this%dz = 0.0d0; this%dy = 0.0d0; this%dx = 0.0d0

      !do na = 1, this%lattice%nrec
      do na_glob = start_atom, end_atom
         na = g2l_map(na_glob)
         do ie = 1, this%en%channels_ldos + 10
            do i = 1, 9
               this%dz(ie, na) = this%dz(ie, na) - aimag(this%green%g0(i, i, ie, na) - this%green%g0(i + 9, i + 9, ie, na))/pi
               this%dy(ie, na) = this%dy(ie, na) - aimag(i_unit*this%green%g0(i, i + 9, ie, na) - i_unit*this%green%g0(i + 9, i, ie, na))/pi
               this%dx(ie, na) = this%dx(ie, na) - aimag(this%green%g0(i, i + 9, ie, na) + this%green%g0(i + 9, i, ie, na))/pi
            end do
         end do
      end do
   end subroutine calculate_projected_dos

   subroutine calculate_orbital_dos(this)
      use mpi_mod
      class(bands) :: this
      ! Local variables
      integer :: i, j ! Orbital index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index

      integer :: na_glob

      call this%calculate_projected_green()

      this%d_orb = 0.0_rp

      !do na = 1, this%lattice%nrec
      do na_glob = start_atom, end_atom
         na = g2l_map(na_glob)
         do ie = 1, this%en%channels_ldos + 10
            do i = 1, 9
               do j = 1, 9
                  this%d_orb(j, i, 1, ie, na) = this%d_orb(j, i, 1, ie, na) &
                                                - aimag(this%green%g0(j + 9, i, ie, na) + this%green%g0(j, i + 9, ie, na))/pi
                  this%d_orb(j, i, 2, ie, na) = this%d_orb(j, i, 1, ie, na) &
                                                - aimag(i_unit*this%green%g0(j + 9, i, ie, na) - i_unit*this%green%g0(j, i + 9, ie, na))/pi
                  this%d_orb(j, i, 3, ie, na) = this%d_orb(j, i, 1, ie, na) &
                                                - aimag(this%green%g0(j, i, ie, na) - this%green%g0(j + 9, i + 9, ie, na))/pi
               end do
            end do
         end do
      end do
   end subroutine calculate_orbital_dos

   subroutine calculate_projected_green(this)
      class(bands) :: this
      ! Local variables
      integer :: i, j ! Orbital index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index

      integer :: na_glob

      this%g0_z = 0.0d0; this%g0_y = 0.0d0; this%g0_x = 0.0d0

      !do na = 1, this%lattice%nrec
      do na_glob = start_atom, end_atom
         na = g2l_map(na_glob)
         do ie = 1, this%en%channels_ldos + 10
            do i = 1, 9
               do j = 1, 9
                  this%g0_z(i, j, ie, na) = this%g0_z(i, j, ie, na) + (this%green%g0(i, i, ie, na) - this%green%g0(i + 9, i + 9, ie, na))
                  this%g0_y(i, j, ie, na) = this%g0_y(i, j, ie, na) + (i_unit*this%green%g0(i, i + 9, ie, na) - i_unit*this%green%g0(i + 9, i, ie, na))
                  this%g0_x(i, j, ie, na) = this%g0_z(i, j, ie, na) + (this%green%g0(i, i + 9, ie, na) + this%green%g0(i + 9, i, ie, na))
               end do
            end do
         end do
      end do
   end subroutine calculate_projected_green

   subroutine calculate_pl(this)
      use mpi_mod
      class(bands) :: this
      ! Local variables
      integer :: i ! orbital index
      integer :: is ! spin channel index
      integer :: ia ! atom index
      integer :: plusbulk
      real(rp) :: rq, dnu, pli, delta2 ! Local variables

      !integer :: ia_glob

      !do ia=1, this%lattice%nrec
      do ia = start_atom, end_atom
         !ia = g2l_map(ia_glob)
         plusbulk = this%lattice%nbulk + ia
         do is = 1, 2
            do i = 1, 3
               rq = 1/this%symbolic_atom(plusbulk)%potential%qpar(i - 1, is)
               delta2 = this%symbolic_atom(plusbulk)%potential%srdel(i - 1, is)*this%symbolic_atom(plusbulk)%potential%srdel(i - 1, is)
               dnu = (i - 1.) + &
                     (2.*(i - 1) + 1.)/ &
                     (rq*(this%symbolic_atom(plusbulk)%potential%c(i - 1, is) - this%symbolic_atom(plusbulk)%potential%gravity_center(i, is))/2./(2*(i - 1) + 1.)/ &
                      (this%symbolic_atom(plusbulk)%potential%c(i - 1, is) - this%symbolic_atom(plusbulk)%potential%gravity_center(i, is) - delta2*rq) - 1.)
               pli = -atan(dnu)/pi + 0.5d0 + INT(this%symbolic_atom(plusbulk)%potential%pl(i - 1, is))
               this%symbolic_atom(plusbulk)%potential%pl(i - 1, is) = pli
            end do
         end do
      end do
   end subroutine calculate_pl

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the magnetic moments mx, my, mz and mtot
   !---------------------------------------------------------------------------
   subroutine calculate_magnetic_torques(this)
      use mpi_mod
      class(bands) :: this
      ! Local variables
      integer :: na ! Atom index

      real(rp), dimension(3) :: I_loc
      real(rp), dimension(3) :: tau_loc
      real(rp) :: pref_0, pref_1, fx, fy, fz, tx, ty, tz
      integer :: up = 1
      integer :: dw = 2
      integer :: d = 2
      integer :: na_loc
      integer :: plusbulk

      call this%calculate_projected_dos()

      !do na=1, this%lattice%nrec
      do na = start_atom, end_atom
         plusbulk = this%lattice%nbulk + na
         na_loc = g2l_map(na)

         pref_0 = this%symbolic_atom(plusbulk)%potential%c(d, up)*this%symbolic_atom(plusbulk)%potential%srdel(d, dw) &
                  /this%symbolic_atom(plusbulk)%potential%srdel(d, up) - this%symbolic_atom(plusbulk)%potential%c(d, dw) &
                  *this%symbolic_atom(plusbulk)%potential%srdel(d, up)/this%symbolic_atom(plusbulk)%potential%srdel(d, dw)
         pref_1 = this%symbolic_atom(plusbulk)%potential%srdel(d, dw)/this%symbolic_atom(plusbulk)%potential%srdel(d, up) &
                  - this%symbolic_atom(plusbulk)%potential%srdel(d, up)/this%symbolic_atom(plusbulk)%potential%srdel(d, dw)

         I_loc = pref_0*this%symbolic_atom(plusbulk)%potential%mom0 &
                 - pref_1*this%symbolic_atom(plusbulk)%potential%mom1

         tau_loc = cross_product(this%symbolic_atom(plusbulk)%potential%mom, I_loc)*ry2tesla

         this%mag_for(:, na_loc) = I_loc(:)*ry2tesla

         !this%mag_for(1:2,na_loc) = 0.0d0

         fx = this%mag_for(1, na_loc)
         fy = this%mag_for(2, na_loc)
         fz = this%mag_for(3, na_loc)

         tx = tau_loc(1)
         ty = tau_loc(2)
         tz = tau_loc(3)

         call g_logger%info('Magnetic field on atom'//fmt('i4', na)//' is '//fmt('f16.6', fx)//' '//fmt('f16.6', fy)//' '//fmt('f16.6', fz), __FILE__, __LINE__)
         call g_logger%info('Magnetic torque on atom'//fmt('i4', na)//' is '//fmt('f16.6', tx)//' '//fmt('f16.6', ty)//' '//fmt('f16.6', tz), __FILE__, __LINE__)

         !print ´(a,i4,a, 3f12.6)´ , "Magnetic mom0 for atom ", na, "=",this%symbolic_atom(plusbulk)%potential%mom0
         !print ´(a,i4,a, 3f12.6)´ , "Magnetic mom1 for atom ", na, "=", this%symbolic_atom(plusbulk)%potential%mom1
         !print ´(a,i4,a, 3f12.6)´ , "Field prefactors for atom ", na, "=", pref_0, pref_1
         !print ´(a,i4,a, 3f12.6)´ , "Magnetic force for atom ", na, "=", I_loc
         !print ´(a,i4,a, 3f12.6)´ , "Magnetic torque for atom ", na, "=", tau_loc
         !print ´(a,i4,a, 3f12.6)´ , "Field prefactors for atom ", na, "=", pref_0 * this%symbolic_atom(plusbulk)%potential%mom0
         !print ´(a,i4,a, 3f12.6)´ , "Field prefactors for atom ", na, "=", pref_1 * this%symbolic_atom(plusbulk)%potential%mom1
      end do
   end subroutine calculate_magnetic_torques
end module bands_mod
