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
      procedure :: calculate_angles
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
      real(rp), dimension(:, :), allocatable :: dosia
      real(rp), dimension(:, :, :), allocatable :: dosial

      allocate(dosia(this%lattice%nrec, this%en%channels_ldos + 10), dosial(this%lattice%nrec, 18, this%en%channels_ldos + 10))

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

      deallocate(dosia)
      deallocate(dosial)
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
      call this%calculate_orbital_moments()

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

         call g_logger%info('Spin moment of atom'//fmt('i4', na)//' is '//fmt('f10.6', this%symbolic_atom(this%lattice%nbulk + na)%potential%mtot), __FILE__, __LINE__)
         mx = this%symbolic_atom(this%lattice%nbulk + na)%potential%mx
         my = this%symbolic_atom(this%lattice%nbulk + na)%potential%my
         mz = this%symbolic_atom(this%lattice%nbulk + na)%potential%mz

         if (this%control%nsp < 3) this%symbolic_atom(this%lattice%nbulk + na)%potential%mom(:) = [0.0d0, 0.0d0, 1.00d0]

         if (this%recursion%hamiltonian%local_axis) then
            call g_logger%info('Local spin moment projections of atom'//fmt('i4', na)//' is '//fmt('f10.6', mx)//' '//fmt('f10.6', my)//' '//fmt('f10.6', mz), __FILE__, __LINE__)
         else
            call g_logger%info('Spin moment projections of atom'//fmt('i4', na)//' is '//fmt('f10.6', mx)//' '//fmt('f10.6', my)//' '//fmt('f10.6', mz), __FILE__, __LINE__)
         end if
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
      real(rp), dimension(this%en%channels_ldos+10) :: lxi, lyi, lzi
      integer :: i, j ! Orbital index
      integer :: mdir ! Magnetic index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index
      real(rp), dimension(9, 9, 3) :: l_orb
      complex(rp), dimension(9, 9) :: mLx, mLy, mLz
      complex(rp), dimension(18, 18) :: mLx_ext, mLy_ext, mLz_ext
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

      mLx_ext = 0.0d0
      mLx_ext(1:9, 1:9) = mLx(:, :)
      mLx_ext(10:18, 10:18) = mLx(:, :)
      mLy_ext = 0.d00
      mLy_ext(1:9, 1:9) = mLy(:, :)
      mLy_ext(10:18, 10:18) = mLy(:, :)
      mLz_ext = 0.0d0
      mLz_ext(1:9, 1:9) = mLz(:, :)
      mLz_ext(10:18, 10:18) = mLz(:, :)

      call this%calculate_orbital_dos()

      !do na=1, this%lattice%nrec
      do na = start_atom, end_atom
         na_loc = g2l_map(na)

         l_orb = 0.0d0
         lx = 0.0d0; ly = 0.0d0; lz = 0.d0


         do ie = 1, this%en%channels_ldos+10
            lxi(ie) = imtrace(matmul(mLx_ext, this%green%g0(:, :, ie, na_loc)))
            lyi(ie) = imtrace(matmul(mLy_ext, this%green%g0(:, :, ie, na_loc)))
            lzi(ie) = imtrace(matmul(mLz_ext, this%green%g0(:, :, ie, na_loc)))
         end do

         call simpson_m(lx, this%en%edel, this%en%fermi, this%nv1, lxi, this%e1, 0, this%en%ene)
         call simpson_m(ly, this%en%edel, this%en%fermi, this%nv1, lyi, this%e1, 0, this%en%ene)
         call simpson_m(lz, this%en%edel, this%en%fermi, this%nv1, lzi, this%e1, 0, this%en%ene)
         !do mdir = 1, 3
         !   do i = 1, 9
         !      do j = 1, 9
         !         call simpson_m(l_orb(i, j, mdir), this%en%edel, this%en%fermi, this%nv1, this%d_orb(i, j, mdir, :, na_loc), this%e1, 0, this%en%ene)
         !      end do
         !   end do
         !end do

         ! p contribution 
         ! lz
         !lz = lz + (l_orb(4, 4, 3) - l_orb(2, 2, 3)) 
         !write(*,'(2f10.6)') l_orb(4, 4, 3), l_orb(2, 2, 3)
         ! lx
         !lx = lx + (l_orb(4, 4, 1) - l_orb(2, 2, 1)) 
         !write(*,'(2f10.6)') l_orb(4, 4, 1), l_orb(2, 2, 1)
         ! ly
         !ly = ly + (l_orb(4, 4, 2) - l_orb(2, 2, 2)) 
         !write(*,'(2f10.6)') l_orb(4, 4, 2), l_orb(2, 2, 2)

         ! d contribution (up +down)
         ! lz
         !lz = lz + 2_rp * (l_orb(9, 9, 3) - l_orb(5, 5, 3)) + (l_orb(8, 8, 3) - l_orb(6, 6, 3))
         !write(*,'(4f10.6)') l_orb(9, 9, 3), l_orb(5, 5, 3), l_orb(8, 8, 3), l_orb(6, 6, 3)
         ! lx
         !lx = lx + 2_rp * (l_orb(9, 9, 1) - l_orb(5, 5, 1)) + (l_orb(8, 8, 1) - l_orb(6, 6, 1))
         !write(*,'(4f10.6)') l_orb(9, 9, 1), l_orb(5, 5, 1), l_orb(8, 8, 1), l_orb(6, 6, 1)
         ! ly
         !ly = ly + 2_rp * (l_orb(9, 9, 2) - l_orb(5, 5, 2)) + (l_orb(8, 8, 2) - l_orb(6, 6, 2)) 
         !write(*,'(4f10.6)') l_orb(9, 9, 2), l_orb(5, 5, 2), l_orb(8, 8, 2), l_orb(6, 6, 2)

         lz =  - (lz / pi) 
         lx =  - (lx / pi) 
         ly =  - (ly / pi) 

         !lz = -0.5_rp * rtrace9(matmul(l_orb(:, :, 3), mLz))
         !lx = -0.5_rp * rtrace9(matmul(l_orb(:, :, 1), mLx))
         !ly = -0.5_rp * rtrace9(matmul(l_orb(:, :, 2), mLy))

         call g_logger%info('Orbital moment of atom'//fmt('i4', na)//' is '//fmt('f10.6', lx)//' '//fmt('f10.6', ly)//' '//fmt('f10.6', lz), __FILE__, __LINE__)
         this%symbolic_atom(this%lattice%nbulk + na)%potential%lmom(1) = lx
         this%symbolic_atom(this%lattice%nbulk + na)%potential%lmom(2) = ly
         this%symbolic_atom(this%lattice%nbulk + na)%potential%lmom(3) = lz
      end do
   end subroutine calculate_orbital_moments

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
                                                - aimag(this%green%g0(j, i + 9, ie, na) + this%green%g0(j + 9, i, ie, na))/pi
                  this%d_orb(j, i, 2, ie, na) = this%d_orb(j, i, 2, ie, na) &
                                                - aimag(i_unit*this%green%g0(j, i + 9, ie, na) - i_unit*this%green%g0(j + 9, i, ie, na))/pi
                  this%d_orb(j, i, 3, ie, na) = this%d_orb(j, i, 3, ie, na) &
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
                  this%g0_x(i, j, ie, na) = this%g0_x(i, j, ie, na) + (this%green%g0(i, i + 9, ie, na) + this%green%g0(i + 9, i, ie, na))
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

   !**************************************************************************
   !> @brief Calculate angles between magnetic and orbital moments for all atoms.
   !> 
   !> This subroutine computes the angles between the magnetic moments and
   !> orbital moments of all atoms in a given system. The angles are calculated
   !> using the dot product of the vectors, normalized by the magnitudes of the
   !> vectors. The results are stored in two matrices: one for the magnetic moments
   !> and one for the orbital moments.
   !> 
   !> @param[in]  this          A derived band type 
   !> @param[in]  magmom        A 2D real array (nrec x 3) containing the magnetic moments for each atom.
   !> @param[in]  lmom          A 2D real array (nrec x 3) containing the orbital moments for each atom.
   !> @param[out] angles_magmom A 2D real array (nrec x nrec) to store the angles between magnetic moments of all atom pairs.
   !> @param[out] angles_lmom   A 2D real array (nrec x nrec) to store the angles between orbital moments of all atom pairs.
   !>
   !> @note The angles are calculated using the following formula:
   !>       \f$ \text{angle} = \cos^{-1} \left( \frac{\mathbf{a} \cdot \mathbf{b}}{|\mathbf{a}| |\mathbf{b}|} \right) \f$
   !>       where \f$\mathbf{a}\f$ and \f$\mathbf{b}\f$ are the moment vectors.
   !**************************************************************************
   subroutine calculate_angles(this, magmom, lmom)
      implicit none
      ! Input
      class(bands), intent(inout) :: this
      real(rp), dimension(this%lattice%nrec, 3), intent(in) :: magmom, lmom
      ! Local variables
      integer :: i, j, cols_per_line
      real(rp) :: dot_prod, mag_a, mag_b
      real(rp), dimension(this%lattice%nrec, this%lattice%nrec) :: angles_magmom
      real(rp), dimension(this%lattice%nrec, this%lattice%nrec) :: angles_lmom
      character(len=sl) :: format_string   

      ! Loop over all pairs of atoms
      do i = 1, this%lattice%nrec
         do j = 1, this%lattice%nrec
            if (i /= j) then
               ! Calculate the angle between magnetic moments
               dot_prod = dot_product(magmom(i, :), magmom(j, :))
               mag_a = sqrt(sum(magmom(i, :)**2))
               mag_b = sqrt(sum(magmom(j, :)**2))
               angles_magmom(i, j) = acos(dot_prod / (mag_a * mag_b))
   
               ! Calculate the angle between orbital moments
               dot_prod = dot_product(lmom(i, :), lmom(j, :))
               mag_a = sqrt(sum(lmom(i, :)**2))
               mag_b = sqrt(sum(lmom(j, :)**2))
               angles_lmom(i, j) = acos(dot_prod / (mag_a * mag_b))
            else
               ! Set the diagonal to zero since the angle between the same vectors is undefined
               angles_magmom(i, j) = 0.0_rp
               angles_lmom(i, j) = 0.0_rp
            end if
         end do
      end do

      ! Determine the number of columns per line
      cols_per_line = this%lattice%nrec
      write(format_string, '(A,I4,A)') '(', cols_per_line, 'f16.6)'

      ! Write the results to output files
      open(unit=10, file='angles_magmom.out', status='replace')
      open(unit=20, file='angles_lmom.out', status='replace')

      do i = 1, this%lattice%nrec
         write(10, format_string) angles_magmom(i, 1:this%lattice%nrec)*rad2deg
         write(20, format_string) angles_lmom(i, 1:this%lattice%nrec)*rad2deg
      end do

      close(10)
      close(20)
   end subroutine calculate_angles

end module bands_mod
