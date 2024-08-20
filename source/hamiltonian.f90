!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Hamiltonian
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
!> Module to handle procedures related to the Hamiltonian
!------------------------------------------------------------------------------

module hamiltonian_mod

   use control_mod
   use symbolic_atom_mod
   use element_mod
   use potential_mod
   use lattice_mod
   use charge_mod
   use precision_mod, only: rp
   use math_mod
   use string_mod
   use logger_mod, only: g_logger
   use timer_mod, only: g_timer
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   type :: ArrayType
      integer, allocatable :: val(:,:)
   end type ArrayType

   !> Module´s main procedure
   type, public :: hamiltonian
      !> Charge
      class(charge), pointer :: charge
      !> Lattice
      class(lattice), pointer :: lattice
      !> Control
      class(control), pointer :: control

      !> Spin-orbit coupling Hamiltonian
      complex(rp), dimension(:, :, :), allocatable :: lsham
      !> Torque operator T=[o, Hso]
      complex(rp), dimension(:, :, :, :), allocatable :: tmat
      !> Bulk Hamiltonian
      complex(rp), dimension(:, :, :, :), allocatable :: ee, eeo, eeoee
      !> Local Hamiltonian
      complex(rp), dimension(:, :, :, :), allocatable :: hall, hallo
      !> Hamiltonian built in chbar_nc (description to be improved)
      complex(rp), dimension(:, :, :, :), allocatable :: hmag
      !> Hamiltonian built in ham0m_nc (description to be improved
      complex(rp), dimension(:, :, :), allocatable :: hhmag
      !> Overlap Hamiltonian
      complex(rp), dimension(:, :, :), allocatable :: obarm
      !> Gravity center Hamiltonian
      complex(rp), dimension(:, :, :), allocatable :: enim
      !> Logical variable to include hoh term
      logical :: hoh
      !> Rotate Hamiltonian to local spin axis
      logical :: local_axis
      !> Add orbital polarization to Hamiltonian
      logical :: orb_pol
      !> Bulk Hamiltonian (backup for rotation)
      complex(rp), dimension(:, :, :, :), allocatable :: ee_glob, eeo_glob
      !> Local Hamiltonian (backup for rotation)
      complex(rp), dimension(:, :, :, :), allocatable :: hall_glob, hallo_glob
    !!> Spin-orbit coupling Hamiltonian (backup for rotation)
      !complex(rp), dimension(:, :, :), allocatable :: lsham
      !> Gravity center Hamiltonian (backup for rotation)
      complex(rp), dimension(:, :, :), allocatable :: enim_glob

      !> On-site potential for LDA+U+J correction
      real(rp), dimension(:,:,:), allocatable :: hubbard_pot
      real(rp), dimension(:,:,:,:), allocatable :: hubbard_v_pot

      !> Hubbard parameters for LDA+U+J implementation
      real(rp), dimension(:,:), allocatable :: hubbard_u
      real(rp), dimension(:,:), allocatable :: hub_u_sort
      real(rp), dimension(:,:), allocatable :: hubbard_j
      real(rp), dimension(:,:), allocatable :: hub_j_sort
      real(rp), dimension(:,:,:,:), allocatable :: hubbard_v ! (atom type 1, atom type 2, l1, l2) with l = 1,2,3,4 : s,p,d,f
      !> Orbitals for Hubbard U+J+V
      character(len=10), dimension(:), allocatable :: uj_orb
      !> Array that stores orbitals for +V correction
      type(ArrayType), dimension(:,:), allocatable :: orbs_v ! (atom type 1, atom type 2, l1, l2) with l = 1,2,3,4 : s,p,d,f
      !> Array that counts orbital pairs for +V correction
      integer, dimension(:,:), allocatable :: orbs_v_num 
      !> Slater integral array
      real(rp), dimension(:,:,:), allocatable :: F
      ! To be able to get the character 's', 'p', 'd', 'f' from element l+1 = 1,2,3,4 for printing purposes
      character(len=1), dimension(4) :: orb_conv 
      !> Orbital configuration for Hubbard U. 1=s, 2=p, 3=d, 4=sp, 5=sd, 6=pd, 7=spd, 8=f, 9=sf, 10=df, 11=pf. Used for iterating purposes.
      integer, dimension(:), allocatable :: hubbard_orb_config
      !> Maximum orbital number to create arrays of appropriate sizes. 0=s, 1=p, 2=d, 3=f
      integer :: hubbard_orb_max
      !> Maximum number of orbitals over the atoms (1 = s or p or d or f, 3 = spd)
      integer :: hubbard_nmb_orb
      !> Logical variables to include Hubbard U & J
      logical :: hubbardU_check = .false.
      logical :: hubbardJ_check = .false.
      logical :: hubbardV_check = .false.
      !> Logical variable to include self-consistent calculation of Hubbard U.
      logical :: hubbardU_sc_check = .false.
      !> Which atoms and orbitals to include self-consistent Hubbard U correction. (0 no correction. 1 include correction)
      integer, dimension(:,:), allocatable :: hubbard_u_sc
   contains
      procedure :: build_lsham
      procedure :: build_bulkham
      procedure :: build_locham
      procedure :: build_obarm
      procedure :: build_enim
      procedure :: build_from_paoflow
      procedure :: build_from_paoflow_opt
      procedure :: torque_operator_collinear
      procedure :: rs2pao
      procedure :: chbar_nc
      procedure :: ham0m_nc
      procedure :: hmfind
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: rotate_to_local_axis
      procedure :: rotate_from_local_axis
      final     :: destructor
   end type hamiltonian

   interface hamiltonian
      procedure :: constructor
   end interface

   

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] charge_obj Variable that holds charge_mod properties
   !> @return type(hamiltonian)
   !---------------------------------------------------------------------------
   function constructor(charge_obj) result(obj)
      type(hamiltonian) :: obj
      type(charge), target, intent(in) :: charge_obj

      obj%charge => charge_obj
      obj%lattice => charge_obj%lattice
      obj%control => charge_obj%lattice%control

      call obj%restore_to_default()
      call obj%build_from_file()
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(hamiltonian) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%lsham)) call g_safe_alloc%deallocate('hamiltonian.lsham', this%lsham)
      if (allocated(this%tmat)) call g_safe_alloc%deallocate('hamiltonian.tmat', this%tmat)
      if (allocated(this%ee)) call g_safe_alloc%deallocate('hamiltonian.ee', this%ee)
      if (allocated(this%hmag)) call g_safe_alloc%deallocate('hamiltonian.hmag', this%hmag)
      if (allocated(this%hhmag)) call g_safe_alloc%deallocate('hamiltonian.hhmag', this%hhmag)
      if (allocated(this%hall)) call g_safe_alloc%deallocate('hamiltonian.hall', this%hall)
      if (allocated(this%eeo)) call g_safe_alloc%deallocate('hamiltonian.eeo', this%eeo)
      if (allocated(this%eeoee)) call g_safe_alloc%deallocate('hamiltonian.eeoee', this%eeoee)
      if (allocated(this%hallo)) call g_safe_alloc%deallocate('hamiltonian.hallo', this%hallo)
      if (allocated(this%obarm)) call g_safe_alloc%deallocate('hamiltonian.obarm', this%obarm)
      if (allocated(this%enim)) call g_safe_alloc%deallocate('hamiltonian.enim', this%enim)
      if (allocated(this%ee_glob)) call g_safe_alloc%deallocate('hamiltonian.ee_glob', this%ee_glob)
      if (allocated(this%eeo_glob)) call g_safe_alloc%deallocate('hamiltonian.eeo_glob', this%eeo_glob)
      if (allocated(this%enim_glob)) call g_safe_alloc%deallocate('hamiltonian.enim_glob', this%enim_glob)
#else
      if (allocated(this%lsham)) deallocate (this%lsham)
      if (allocated(this%tmat)) deallocate (this%tmat)
      if (allocated(this%ee)) deallocate (this%ee)
      if (allocated(this%eeo)) deallocate (this%eeo)
      if (allocated(this%eeoee)) deallocate (this%eeoee)
      if (allocated(this%hmag)) deallocate (this%hmag)
      if (allocated(this%hhmag)) deallocate (this%hhmag)
      if (allocated(this%hall)) deallocate (this%hall)
      if (allocated(this%hallo)) deallocate (this%hallo)
      if (allocated(this%obarm)) deallocate (this%obarm)
      if (allocated(this%enim)) deallocate (this%enim)
      if (allocated(this%ee_glob)) deallocate (this%ee_glob)
      if (allocated(this%eeo_glob)) deallocate (this%eeo_glob)
      if (allocated(this%enim_glob)) deallocate (this%enim_glob)
      if (allocated(this%hubbard_u)) deallocate (this%hubbard_u)
      if (allocated(this%hubbard_j)) deallocate (this%hubbard_j)
      if (allocated(this%hubbard_v)) deallocate (this%hubbard_v)
      if (allocated(this%hub_u_sort)) deallocate (this%hub_u_sort)
      if (allocated(this%hub_j_sort)) deallocate (this%hub_j_sort)
      if (allocated(this%uj_orb)) deallocate (this%uj_orb)
      if (allocated(this%orbs_v)) deallocate (this%orbs_v)
      if (allocated(this%orbs_v_num)) deallocate (this%orbs_v_num)
      if (allocated(this%F)) deallocate(this%F)
      if (allocated(this%hubbard_orb_config)) deallocate (this%hubbard_orb_config)
      if (allocated(this%hubbard_pot)) deallocate (this%hubbard_pot)
      if (allocated(this%hubbard_v_pot)) deallocate (this%hubbard_v_pot)
      if (allocated(this%hubbard_u_sc)) deallocate (this%hubbard_u_sc)
#endif
   end subroutine destructor

   ! Member functions
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read parameters from input file
   !> Edited to prompt Hubbard U+J correction by calculating Slater integrals,
   !> check if Hubbard U and/or J should be in use and if there are incorrect input data, 
   !> by Emil Beiersdorf and Viktor Frilén 17.07.2024
   !---------------------------------------------------------------------------
   subroutine build_from_file(this)
      class(hamiltonian), intent(inout) :: this

      ! variables associated with the reading processes
      integer :: iostatus, funit, i, j, li, lj, k, max_orbs, length, cntr, na, l
      ! integer, dimension(this%lattice%nrec, this%lattice%nrec) :: orbs_v_num ! (atom1, atom2) = number of mutual orbs e.g. hubbard_v(1, 2, [p,d], [s,p]) -> orbs_v_num(1,2)=2
      ! Logical variable to check if input data is good
      logical :: implem_check = .true.


      include 'include_codes/namelists/hamiltonian.f90'

      hoh = this%hoh
      local_axis = this%local_axis
      orb_pol = this%orb_pol

      call move_alloc(this%hubbard_u_sc, hubbard_u_sc)
      call move_alloc(this%hubbard_u, hubbard_u)
      call move_alloc(this%hubbard_j, hubbard_j)
      call move_alloc(this%hubbard_v, hubbard_v) 
      call move_alloc(this%uj_orb, uj_orb)

      ! Reading
      open (newunit=funit, file=this%control%fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(this%control%fname)//' not found', __FILE__, __LINE__)
      end if

      !> Trick Ramon talked about
      read (funit, nml=hamiltonian, iostat=iostatus)

      max_orbs = 1
      do i = 1, this%lattice%nrec   
         max_orbs = max(max_orbs, len_trim(uj_orb(i)))
      end do

      
      if (size(hubbard_u,2) .ne. max_orbs) then
         deallocate(hubbard_u)
         deallocate(hubbard_j)
         deallocate(this%F)
         allocate(hubbard_u(this%lattice%nrec, max_orbs))
         allocate(hubbard_j(this%lattice%nrec, max_orbs))
         allocate(this%F(this%lattice%nrec, 4, 4))
         hubbard_u(:,:) = 0.0d0
         hubbard_j(:,:) = 0.0d0
         this%F(:,:,:) = 0.0d0
      end if
      
      rewind (funit)
      read (funit, nml=hamiltonian, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      this%hoh = hoh
      this%local_axis = local_axis
      this%orb_pol = orb_pol

      call move_alloc(hubbard_u_sc, this%hubbard_u_sc)
      call move_alloc(hubbard_u, this%hubbard_u)
      call move_alloc(hubbard_j, this%hubbard_j)
      call move_alloc(hubbard_v, this%hubbard_v)
      call move_alloc(uj_orb, this%uj_orb)

      ! Checks if self-consistent 
      do na = 1, this%lattice%nrec
         do l = 1, 4
            if ( this%hubbard_u_sc(na,l) == 1 ) then
               if (l == 1) then
                  print *, ''
                  print *, 'Self-consistent U added for s-electrons on atom', na
                  this%hubbardU_sc_check = .true.
               else if (l == 2) then
                  print *, ''
                  print *, 'Self-consistent U added for p-electrons on atom', na
                  this%hubbardU_sc_check = .true.
               else if (l == 3) then
                  print *, ''
                  print *, 'Self-consistent U added for d-electrons on atom', na
                  this%hubbardU_sc_check = .true.
               else if ( l == 4 ) then
                  print *, ''
                  print *, 'ERROR'
                  print *, 'Self-consistent U calculation not implemented for f-electrons'
                  print *, 'STOPS PROGRAM'
                  stop
               end if
            else if ( this%hubbard_u_sc(na,l) .gt. 1 ) then
               print *, ''
               print *, 'WRONG INPUT. hubbard_u_sc input was larger than 1, but only 0 and 1 is allowed.'
               print *, 'STOPS PROGRAM'
               stop
            end if
         end do
      end do
      ! Check if correct input is given for hubbard_u_sc
      if ( this%hubbardU_sc_check ) then
         print *, ''
         print *, '----------------------------------------------------------------------------------------'
         print *, 'Self-consistent calculation of Hubbard U initiated.'
         print *, '----------------------------------------------------------------------------------------'
      end if

      !Saves the Hubbard U orbital configuration in a variable
      do i = 1, this%lattice%nrec
         if (adjustl(this%uj_orb(i)) == 's') then
            this%hubbard_orb_config(i) = 1
            this%hubbard_nmb_orb = max(this%hubbard_nmb_orb, 1)
         else if (adjustl(this%uj_orb(i)) == 'p') then
            this%hubbard_orb_config(i) = 2
            this%hubbard_orb_max = max(this%hubbard_orb_max, 1)
            this%hubbard_nmb_orb = max(this%hubbard_nmb_orb, 1)
         else if (adjustl(this%uj_orb(i)) == 'd') then
            this%hubbard_orb_config(i) = 3
            this%hubbard_orb_max = max(this%hubbard_orb_max, 2)
            this%hubbard_nmb_orb = max(this%hubbard_nmb_orb, 1)
         else if (adjustl(this%uj_orb(i)) == 'sp') then
            this%hubbard_orb_config(i) = 4
            this%hubbard_orb_max = max(this%hubbard_orb_max, 1)
            this%hubbard_nmb_orb = max(this%hubbard_nmb_orb, 2)
         else if (adjustl(this%uj_orb(i)) == 'sd') then
            this%hubbard_orb_config(i) = 5
            this%hubbard_orb_max = max(this%hubbard_orb_max, 2)
            this%hubbard_nmb_orb = max(this%hubbard_nmb_orb, 2)
         else if (adjustl(this%uj_orb(i)) == 'pd') then
            this%hubbard_orb_config(i) = 6
            this%hubbard_orb_max = max(this%hubbard_orb_max, 2)
            this%hubbard_nmb_orb = max(this%hubbard_nmb_orb, 2)
         else if (adjustl(this%uj_orb(i)) == 'spd') then
            this%hubbard_orb_config(i) = 7
            this%hubbard_orb_max = max(this%hubbard_orb_max, 2)
            this%hubbard_nmb_orb = max(this%hubbard_nmb_orb, 3)
         else if (adjustl(this%uj_orb(i)) == 'f') then
            this%hubbard_orb_config(i) = 8
            this%hubbard_orb_max = max(this%hubbard_orb_max, 3)
            this%hubbard_nmb_orb = max(this%hubbard_nmb_orb, 1)
         else if (adjustl(this%uj_orb(i)) == 'sf') then
            this%hubbard_orb_config(i) = 9
            this%hubbard_orb_max = max(this%hubbard_orb_max, 3)
            this%hubbard_nmb_orb = max(this%hubbard_nmb_orb, 2)
         else if (adjustl(this%uj_orb(i)) == 'df') then
            this%hubbard_orb_config(i) = 10
            this%hubbard_orb_max = max(this%hubbard_orb_max, 3)
            this%hubbard_nmb_orb = max(this%hubbard_nmb_orb, 2)
         else if (adjustl(this%uj_orb(i)) == 'pf') then
            this%hubbard_orb_config(i) = 11
            this%hubbard_orb_max = max(this%hubbard_orb_max, 3)
            this%hubbard_nmb_orb = max(this%hubbard_nmb_orb, 2)
         else
            this%hubbard_orb_config(i) = 0
         end if
      end do

      ! Establishes if Hubbard U should be implemented by checking if any Hubbard U is specified
      outer2 : do i = 1, this%lattice%nrec 
         do j = 1, len_trim(this%uj_orb(i)) 
            if (this%hubbard_u(i,j) > 1.0E-10) then ! If a nonzero Hubbard U parameter is specified, Hubbard U is initiated
               this%hubbardU_check = .true.
               exit outer2
            end if
         end do
      end do outer2

      ! Establishes if Hubbard J should be implemented by checking if any Hubbard J is specified
      outer3 : do i = 1, this%lattice%nrec 
         do j = 1, len_trim(this%uj_orb(i)) 
            if (this%hubbard_j(i,j) > 1.0E-10) then ! If a nonzero Hubbard J parameter is specified, Hubbard J is initiated
               this%hubbardJ_check = .true.
               exit outer3
            end if
         end do
      end do outer3

      ! Establishes if Hubbard V should be implemented by checking if any Hubbard V is specified
      outer4 : do i = 1, this%lattice%nrec 
         do j = 1, this%lattice%nrec
            do li = 1, size(this%hubbard_v, 3)
               do lj = 1, size(this%hubbard_v, 4)
                  if (this%hubbard_v(i,j,li,lj) > 1.0E-10) then
                     this%hubbardV_check = .true.
                  end if
               end do
            end do
         end do
      end do outer4

      ! Organizes the orbital Hubbard U & J values into the correct position
      if (this%hubbardU_check) then
         do i = 1, this%lattice%nrec
            do j = 1, count(this%hubbard_u(i,:) > 1.0E-10)
               if (this%uj_orb(i)(j:j) == 's') then
                  this%hub_u_sort(i,1) = this%hubbard_u(i,j)
                  this%hub_j_sort(i,1) = this%hubbard_j(i,j)
               else if (this%uj_orb(i)(j:j) == 'p') then
                  this%hub_u_sort(i,2) = this%hubbard_u(i,j)
                  this%hub_j_sort(i,2) = this%hubbard_j(i,j)
               else if (this%uj_orb(i)(j:j) == 'd') then
                  this%hub_u_sort(i,3) = this%hubbard_u(i,j)
                  this%hub_j_sort(i,3) = this%hubbard_j(i,j)
               else if (this%uj_orb(i)(j:j) == 'f') then
                  this%hub_u_sort(i,4) = this%hubbard_u(i,j)
                  this%hub_j_sort(i,4) = this%hubbard_j(i,j)
               end if
            end do
         end do
      end if
      
      if (this%hubbardV_check) then

      !> Checks if Vij = Vji in input.
         do i = 1, this%lattice%nrec
            do j = 1, this%lattice%nrec
               do li = 1, 4
                  do lj = 1, 4
                     if (abs(this%hubbard_v(i,j,li,lj)) > 1.0E-10 .and. abs(this%hubbard_v(j,i,lj,li)) > 1.0E-10 .and. this%hubbard_v(i,j,li,lj) /= this%hubbard_v(j,i,lj,li)) then
                        print *, 'Incorrect +V implementation for interaction between atom type ',i , 'and ',j , '.'
                        call g_logger%error('Value for V(i,j,li,lj) must be the same as V(j,i,lj,li). Either delete one of them, or set them equal.', __FILE__, __LINE__)
                        stop
                     end if
                  end do
               end do
            end do
         end do
      end if

      if (this%hubbardV_check) then
      !> Finds and stores all orbital pairs (l+1) for each interaction in +V correction

         do i = 1, this%lattice%nrec
            !> Diagonal part
            cntr = 0
            do li = 1, 4
               do lj = 1, 4
                  if (this%hubbard_v(i,i,li,lj) > 1.0E-10) then
                     cntr = cntr + 1
                  end if
               end do
            end do
            this%orbs_v_num(i,i) = cntr
            allocate(this%orbs_v(i,i)%val(cntr,2))
            cntr = 0
            do li = 1, 4
               do lj = 1, 4
                  if (this%hubbard_v(i,i,li,lj) > 1.0E-10) then
                     cntr = cntr + 1
                     this%orbs_v(i,i)%val(cntr,:) = [li, lj]
                  end if
               end do
            end do

            !> Off-diagonal part
            if (i+1 <= this%lattice%nrec) then
               do j = i+1, this%lattice%nrec
                  cntr = 0
                  do li = 1, 4
                     do lj = 1, 4
                        if (this%hubbard_v(i,j,li,lj) > 1.0E-10 .or. this%hubbard_v(j,i,lj,li) > 1.0E-10) then
                           cntr = cntr + 1
                        end if
                     end do
                  end do
                  this%orbs_v_num(i,j) = cntr
                  this%orbs_v_num(j,i) = cntr
                  allocate(this%orbs_v(i,j)%val(cntr,2))
                  allocate(this%orbs_v(j,i)%val(cntr,2))
                  ! this%orbs_v(i,j)%val(:,:) = 0
                  cntr = 0
                  do li = 1, 4
                     do lj = 1, 4
                        if (this%hubbard_v(i,j,li,lj) > 1.0E-10 .or. this%hubbard_v(j,i,lj,li) > 1.0E-10) then
                           cntr = cntr + 1
                           this%orbs_v(i,j)%val(cntr,:) = [li, lj]
                           this%orbs_v(j,i)%val(cntr,:) = [lj, li]
                        end if
                     end do
                  end do
               end do             
            end if
         end do
            

      !> Sets the Vji = Vij for the +V correction based on input and raises error if Vij and Vji disagree if both are specified.
         do i = 1, this%lattice%nrec
            do j = 1, this%lattice%nrec
               if (j /= i) then
                  do li = 1, 4
                     do lj = 1, 4
                        if (abs(this%hubbard_v(i,j,li,lj)) > 1.0E-10 .and. abs(this%hubbard_v(j,i,lj,li)) < 1.0E-10) then
                           this%hubbard_v(j,i,lj,li) = this%hubbard_v(i,j,li,lj)
                        end if
                     end do
                  end do
               end if
            end do
         end do
      end if

      if (this%hubbardU_check) then
      ! Calculates the Slater integrals, F(atom, l+1, k) with k = 1,2,3,4 = l+1 indexing F^(2(k-1))
         do i = 1, this%lattice%nrec ! For each atom
            do j = 1, count(this%hubbard_u(i,:) > 1.0E-10) ! For each nonzero hubbard U value
               if (this%uj_orb(i)(j:j) == 's') then ! If corresponding orbital is s
               ! For s, p, d and f orbitals, F0 = U
                  this%F(i,1,1) = this%hubbard_u(i,j)
               else if (this%uj_orb(i)(j:j) == 'p') then ! If corresponding orbital is p
               ! For p orbitals, J = (1/5)*F2 
                  this%F(i,2,1) = this%hubbard_u(i,j)
                  this%F(i,2,2) = this%hubbard_j(i,j)*5.0_rp
               else if (this%uj_orb(i)(j:j) == 'd') then ! If corresponding orbital is d
               ! For d orbitals, J = (F2 + F4)/14 with F4/F2 ~ 0.625
                  this%F(i,3,1) = this%hubbard_u(i,j)
                  this%F(i,3,2) = 14.0_rp*this%hubbard_j(i,j)/1.625_rp 
                  this%F(i,3,3) = 0.625_rp*this%F(i,3,2)
               else if (this%uj_orb(i)(j:j) == 'f') then ! If corresponding orbital is f
               ! For f orbitals, J = (286F2 + 195F4 + 250F6)/6435 with F4/F2 ~ 0.67 and F6/F2 ~ 0.49
                  this%F(i,4,1) = this%hubbard_u(i,j)
                  this%F(i,4,2) = 6435_rp*this%hubbard_j(i,j)/(286_rp + 195_rp*0.67_rp + 250_rp*0.49_rp)
                  this%F(i,4,3) = 0.67_rp*this%F(i,4,2)
                  this%F(i,4,4) = 0.49_rp*this%F(i,4,2)
               else 
                  call g_logger%error('Only spdf orbitals have implemented Hubbard U+J framework.', __FILE__, __LINE__)
                  error stop
               end if
            end do
         end do      
      end if

      outer : do i = 1, this%lattice%nrec ! Steps through the atoms

         ! ----------------------------- FOR DEBUGGING ---------------------------------
         ! print *, 'U(i) len for i = ', i, ' is ', count(this%hubbard_u(i,:) > 1.0E-10)
         ! print *, 'J(i) len for i = ', i, 'is ', count(this%hubbard_j(i,:) > 1.0E-10), ' while Orb(i) len is', len_trim(this%uj_orb(i))
         ! -----------------------------------------------------------------------------

         ! If no input values are given for some atom(s), default values is set given that Hubbard U and or J should be included
         if (len_trim(this%uj_orb(i)) == 0 .and. count(this%hubbard_u(i,:) > 1.0E-10) == 0 .and. count(this%hubbard_j(i,:) > 1.0E-10) == 0 .and. this%hubbardU_check .and. this%hubbardJ_check) then
            print *, ''
            print *,''
            print *, '----------------------------------------------------------------------------------------'
            print *, 'No input values for Hubbard U+J specified for atom ', i, ', setting default values.'
            print *, '----------------------------------------------------------------------------------------'
         else if (len_trim(this%uj_orb(i)) == 0 .and. count(this%hubbard_u(i,:) > 1.0E-10) == 0 .and. count(this%hubbard_j(i,:) > 1.0E-10) == 0 .and. this%hubbardU_check) then
            print *, ''
            print *,''
            print *, '----------------------------------------------------------------------------------------'
            print *, 'No input values for Hubbard U specified for atom ', i, ', setting default values of zero.'
            print *, '----------------------------------------------------------------------------------------'
         ! Also if only U or J value is unspecified when orbital and J or U is, default value is set.
         else if (len_trim(this%uj_orb(i)) == 0 .and. count(this%hubbard_u(i,:) > 1.0E-10) == 0 .and. count(this%hubbard_j(i,:) > 1.0E-10) == 0 .and. this%hubbardJ_check) then
            print *, ''
            print *,''
            print *, '----------------------------------------------------------------------------------------'
            print *, 'No input values for Hubbard J specified for atom ', i, ', setting default values of zero.'
            print *, '----------------------------------------------------------------------------------------'
         else if (len_trim(this%uj_orb(i)) /= 0 .and. count(this%hubbard_u(i,:) > 1.0E-10) /= 0 .and. count(this%hubbard_j(i,:) > 1.0E-10) == 0 .and. this%hubbardU_check) then
            if (len_trim(this%uj_orb(i)) == count(this%hubbard_u(i,:) > 1.0E-10)) then
               print *, ''
               print *,''
               print *, '----------------------------------------------------------------------------------------'
               print *, 'No J parameter specified for atom ', i, ', setting default value of zero.'
               print *, '----------------------------------------------------------------------------------------'
            else 
               implem_check = .false.
               print *, ''
               print *, '----------------------------------------------------------------------------------------'
               print *, 'Number of orbitals and Hubbard U parameters for atom', i, ' disagrees.'
               print *, '----------------------------------------------------------------------------------------'
            end if
         else if (len_trim(this%uj_orb(i)) /=0 .and. count(this%hubbard_u(i,:) > 1.0E-10) == 0 .and. count(this%hubbard_j(i,:) > 1.0E-10) /= 0 .and. this%hubbardJ_check) then
            if (len_trim(this%uj_orb(i)) == count(this%hubbard_j(i,:) > 1.0E-10)) then
               print *, ''
               print *, ''
               print *, '----------------------------------------------------------------------------------------'
               print *, 'No U parameter specified for atom ', i, ', setting default value of zero'
               print *, '----------------------------------------------------------------------------------------'
            else 
               implem_check = .false.
               print *, ''
               print *, '----------------------------------------------------------------------------------------'
               print *, 'Number of orbitals and Hubbard J parameters for atom', i, ' disagrees.'
               print *, '----------------------------------------------------------------------------------------'
            end if
         else
            ! If the number of orbitals and Hubbard U parameters per atom in input file disagree, then a raised error is prompted
            if (count(this%hubbard_u(i,:) > 1.0E-10) /= len_trim(this%uj_orb(i))) then
               print *, 'Im here'
               implem_check = .false.
               print *, ''
               print *, '----------------------------------------------------------------------------------------'
               print *, 'Number of orbitals and Hubbard U parameters for atom', i, ' disagrees.'
               print *, '----------------------------------------------------------------------------------------'

            ! If the number of orbitals and Hubbard J parameters per atom in input file disagree, then a raised error is prompted 
            else if (count(this%hubbard_j(i,:) > 1.0E-10) /= len_trim(this%uj_orb(i))) then 
               implem_check = .false.
               print *, ''
               print *, '----------------------------------------------------------------------------------------'
               print *, 'Number of orbitals and Hubbard J parameters for atom', i, ' disagrees.'
               print *, '----------------------------------------------------------------------------------------'
            end if
         end if
      end do outer

      !> If the +V (intersite Coulomb reaction) is added, we use the simplified +U correction which is spherically averaged.
      !> This corresponds to setting F^2 = F^4 = J = 0.
      if (this%hubbardV_check) then
         print *, 'HubbardV_check is True.'
         this%hubbard_j = 0.0d0
         this%F(:,:,2) = 0.0d0
         this%F(:,:,3) = 0.0d0
         this%F(:,:,4) = 0.0d0
         call this%lattice%neigh(1)
         print *, 'ijpair after : ', this%lattice%ijpair
         print *, 'ijpair_sorted after : ', this%lattice%ijpair_sorted !(na = num of diff atoms, nn=1 or snn=2, num of nn, atom index)
         print *, 'size(ijpair_sorted) after : ', size(this%lattice%ijpair_sorted)
         print *, 'Nmb of nn : ', size(this%lattice%ijpair_sorted, 3)
         print *, 'njij after : ', this%lattice%njij
         print *, ''
         print *, 'Nearest neighbours'
         do i = 1, this%lattice%nrec
            print *, 'Atom ', i
            do j = 1, size(this%lattice%ijpair_sorted, 3)
               print *, this%lattice%ijpair_sorted(i,1,j,:)
            end do
         end do
         ! stop
      else 
         print *, 'HubbardV_check is False.'
         print *, 'ijpair : ', this%lattice%ijpair
         print *, 'ijpair_sorted : ', this%lattice%ijpair_sorted
         print *, 'njij : ', this%lattice%njij
      end if 

      !> Print Hubbard U+J info
      if (this%hubbardU_check .and. this%hubbardJ_check) then 
         if (implem_check) then
            print *,''
            print *, '----------------------------------------------------------------------------------------'
            print *, 'Stored input values for LDA+U+J'
            print *, '----------------------------------------------------------------------------------------'
            do i = 1, this%lattice%nrec
               print *, 'Atom ', i, ':'
               do j = 1, max_orbs
                  print *, '  Orbital ', this%uj_orb(i)(j:j), ': Hubbard U = ', this%hubbard_u(i,j), ' eV.', ' Hubbard J = ', this%hubbard_j(i,j), 'eV'
               end do
            end do
            print *, ''
            print *, 'Hubbard U+J input data implemented successfully, proceeding...'
            if (this%hubbardV_check) then
               print *, '+V (intersite Coulomb) input data implemented successfully, proceeding...'
               print *, ''
               print *, 'NOTE: +U correction restored to simplified version (F^2 = F^4 = J = 0) since the'
               print *, '      +V correction is asked to be included.'
            end if
            print *, '----------------------------------------------------------------------------------------'

            print *,''
            print *, '----------------------------------------------------------------------------------------'
            print *, 'Slater Integrals.'
            print *, '----------------------------------------------------------------------------------------'
            do i = 1, this%lattice%nrec
               print *, 'Atom ', i, ':'
               do j = 1, 4 ! Orbitals spdf have implementation
                  if (count(this%F(i,j,:) > 1.0E-10) /= 0) then
                     print *, ' Orbital ', this%orb_conv(j), ' :'
                     if (j == 1) then
                        print *, '             F0 = ', this%F(i,j,1), ' eV'
                     else if (j == 2) then
                        print *, '             F0 = ', this%F(i,j,1), ' eV'
                        print *, '             F2 = ', this%F(i,j,2), ' eV'
                     else if (j == 3) then
                        print *, '             F0 = ', this%F(i,j,1), ' eV'
                        print *, '             F2 = ', this%F(i,j,2), ' eV'
                        print *, '             F4 = ', this%F(i,j,3), ' eV'
                     else if (j == 4) then
                        print *, '             F0 = ', this%F(i,j,1), ' eV'
                        print *, '             F2 = ', this%F(i,j,2), ' eV'
                        print *, '             F4 = ', this%F(i,j,3), ' eV'
                        print *, '             F6 = ', this%F(i,j,4), ' eV'
                     end if
                  end if
               end do
            end do

            print *, '----------------------------------------------------------------------------------------'
         else 
            call g_logger%error('Implementation error in input data.', __FILE__, __LINE__)
            error stop
         end if      
         
      else if (this%hubbardU_check) then
         if (implem_check) then
            print *,''
            print *, '----------------------------------------------------------------------------------------'
            print *, 'Stored input values for LDA+U'
            print *, '----------------------------------------------------------------------------------------'
               do i = 1, this%lattice%nrec
                  print *, 'Atom ', i, ':'
                  do j = 1, max_orbs
                     print *, '  Orbital ', this%uj_orb(i)(j:j), ': Hubbard U = ', this%hubbard_u(i,j), ' eV.'
                  end do
               end do
               print *, ''
               print *, 'Hubbard U (no J) input data implemented successfully, proceeding...'
               if (this%hubbardV_check) then
                  print *, '+V (intersite Coulomb) input data implemented successfully, proceeding...'
                  print *, ''
                  print *, 'NOTE: +U correction restored to simplified version (F^2 = F^4 = J = 0) since the'
                  print *, '      +V correction is asked to be included.'
               end if
               print *, '----------------------------------------------------------------------------------------'
         else 
            call g_logger%error('Implementation error in input data.', __FILE__, __LINE__)
               error stop
         end if       

      else if (this%hubbardJ_check) then 
         if (implem_check) then
            print *,''
            print *, '----------------------------------------------------------------------------------------'
            print *, 'Stored input values for LDA+J'
            print *, '----------------------------------------------------------------------------------------'
               do i = 1, this%lattice%nrec
                  print *, 'Atom ', i, ':'
                  do j = 1, max_orbs
                     print *, '  Orbital ', this%uj_orb(i)(j:j), ' Hubbard J = ', this%hubbard_j(i,j), 'eV'
                  end do
               end do
               print *, ''
               print *, 'Hubbard J (no U) input data implemented successfully, proceeding...'
               if (this%hubbardV_check) then
                  print *, '+V (intersite Coulomb) input data implemented successfully, proceeding...'
               end if
               print *, 'N.B. that with +V correction, F^2 = F^4 = J = 0.'
            print *, '----------------------------------------------------------------------------------------'
         else 
            call g_logger%error('Implementation error in input data.', __FILE__, __LINE__)
            error stop
         end if

      else
         print *, ''
         print *, ''
         print *, '----------------------------------------------------------------------------------------'
         print *, 'No Hubbard U+J data was given as input, proceeding without.'
         print *, '----------------------------------------------------------------------------------------'
      end if      
      
      !> Converting the Hubbard U and J values to Ry from eV
      this%hubbard_u = this%hubbard_u/ry2ev
      this%hubbard_j = this%hubbard_j/ry2ev
      this%hubbard_v = this%hubbard_v/ry2ev
      this%F = this%F/ry2ev
      this%hub_u_sort = this%hub_u_sort/ry2ev
      this%hub_j_sort = this%hub_j_sort/ry2ev

      ! Checks if self-consistent U flag and input values of U and J have both been provided. They would interfer with eachother. 
      if ( (this%hubbardU_sc_check) .and. (this%hubbardU_check) .and. (this%hubbardJ_check)) then
         print *, 'ERROR'
         print *, 'Both input values for hubbard_u_sc and hubbard_u + hubbard_j has been provided.'
         print *, 'Only one of them is allowed. Stops program!!'
         stop
      else if ( (this%hubbardU_sc_check) .and. (this%hubbardU_check) ) then
         print *, 'ERROR'
         print *, 'Both input values for hubbard_u_sc and hubbard_u has been provided.'
         print *, 'Only one of them is allowed. Stops program!!'
         stop
      else if ( (this%hubbardU_sc_check) .and. (this%hubbardJ_check)) then 
         print *, 'ERROR'
         print *, 'Both input values for hubbard_u_sc and hubbard_j has been provided.'
         print *, 'Only one of them is allowed. Stops program!!'
         stop
      end if

   end subroutine build_from_file


   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      class(hamiltonian) :: this

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('hamiltonian.lsham', this%lsham, (/18, 18, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.tmat', this%tmat, (/18, 18, 3, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hhmag', this%hhmag, (/9, 9, 4/))
      call g_safe_alloc%allocate('hamiltonian.hmag', this%hmag, (/9, 9, this%charge%lattice%kk, 4/))
      call g_safe_alloc%allocate('hamiltonian.ee', this%ee, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hall', this%hall, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.hall_glob', this%hall_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.ee_glob', this%ee_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      !if (hoh) then
      call g_safe_alloc%allocate('hamiltonian.eeo', this%eeo, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.eeoee', this%eeoee, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hallo', this%hallo, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.obarm', this%obarm, (/18, 18, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.enim', this%enim, (/18, 18, this%charge%lattice%ntype/))
      !end if
      !if (local_axis)  then
      call g_safe_alloc%allocate('hamiltonian.hall_glob', this%hall_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.ee_glob', this%ee_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      !if (hoh) then
      call g_safe_alloc%allocate('hamiltonian.ee0_glob', this%eeo_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hallo_glob', this%hallo_glob, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.enim_glob', this%enim_glob, (/18, 18, this%charge%lattice%ntype/))
      !end if
      !end if
#else
      allocate (this%lsham(18, 18, this%charge%lattice%ntype))
      allocate (this%tmat(18, 18, 3, this%charge%lattice%ntype))
      allocate (this%hhmag(9, 9, 4), this%hmag(9, 9, this%charge%lattice%kk, 4))
      allocate (this%ee(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      print *, 'test : ', maxval(this%charge%lattice%nn(:, 1)) + 1
      ! allocate (this%ee(18, 18, size(this%lattice%ijpair, 3), this%charge%lattice%ntype))
      allocate (this%hall(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      !if (this%hoh) then
      allocate (this%eeo(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%eeoee(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hallo(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%obarm(18, 18, this%charge%lattice%ntype))
      allocate (this%enim(18, 18, this%charge%lattice%ntype))
      !end if
      !if (this%local_axis) then
      allocate (this%ee_glob(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hall_glob(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      !if (this%hoh) then
      allocate (this%eeo_glob(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hallo_glob(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%enim_glob(18, 18, this%charge%lattice%ntype))
      allocate (this%hubbard_u(this%lattice%nrec, 1))
      allocate(this%hubbard_j(this%lattice%nrec, 1))
      allocate(this%hubbard_v(this%lattice%nrec, this%lattice%nrec, 4, 4))
      allocate (this%hub_u_sort(this%lattice%nrec, 4))
      allocate(this%hub_j_sort(this%lattice%nrec, 4))
      allocate(this%orbs_v(this%lattice%nrec, this%lattice%nrec))
      allocate(this%orbs_v_num(this%lattice%nrec, this%lattice%nrec))
      allocate(this%F(this%lattice%nrec, 4, 4))
      allocate (this%uj_orb(this%lattice%nrec))
      allocate (this%hubbard_orb_config(this%lattice%nrec))
      allocate (this%hubbard_pot(18, 18, this%lattice%nrec))
      allocate (this%hubbard_v_pot(18, 18, size(this%ee, 3) , this%lattice%nrec)) ! (lm, l'm', number of NN, number of atom types)
      allocate (this%hubbard_u_sc(this%lattice%nrec,4))
      !end if
      !end if
#endif

      this%lsham(:, :, :) = 0.0d0
      this%tmat(:, :, :, :) = 0.0d0
      this%hhmag(:, :, :) = 0.0d0
      this%hmag(:, :, :, :) = 0.0d0
      this%hall(:, :, :, :) = 0.0d0
      this%ee(:, :, :, :) = 0.0d0
      !if (this%hoh) then
      this%hallo(:, :, :, :) = 0.0d0
      this%eeo(:, :, :, :) = 0.0d0
      this%eeoee(:, :, :, :) = 0.0d0
      this%obarm(:, :, :) = 0.0d0
      this%enim(:, :, :) = 0.0d0
      !end if
      !if (this%local_axis) then
      this%hall_glob(:, :, :, :) = 0.0d0
      this%ee_glob(:, :, :, :) = 0.0d0
      !  if (this%hoh) then
      this%hallo_glob(:, :, :, :) = 0.0d0
      this%eeo_glob(:, :, :, :) = 0.0d0
      this%enim_glob(:, :, :) = 0.0d0
      !  end if
      !end if
      this%hoh = .false.
      this%local_axis = .false.
      this%orb_pol = .false.
      this%hubbard_u(:,:) = 0.0d0
      this%hubbard_j(:,:) = 0.0d0
      this%hubbard_v(:,:,:,:) = 0.0d0
      this%hub_u_sort(:,:) = 0.0d0
      this%hub_j_sort(:,:) = 0.0d0
      this%uj_orb(:) = ''
      this%orbs_v_num(:,:) = 0
      this%F(:,:,:) = 0.0d0
      this%hubbard_orb_config = 0
      this%hubbard_orb_max = 0
      this%hubbard_nmb_orb = 0
      this%hubbardU_check = .false.
      this%hubbardJ_check = .false.
      this%hubbardV_check = .false.
      this%hubbard_pot(:,:,:) = 0.0d0
      this%hubbard_v_pot(:,:,:,:) = 0.0d0
      this%orb_conv(1) = 's'
      this%orb_conv(2) = 'p'
      this%orb_conv(3) = 'd'
      this%orb_conv(4) = 'f'
      this%hubbard_u_sc(:,:) = 0
   end subroutine restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build the spin-orbit coupling hamiltonian according to Wu/Freeman PRB 54, 61 (1996)
   !---------------------------------------------------------------------------
   subroutine build_lsham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k
      complex(rp) :: prefac, sg
      real(rp) :: soc_p, soc_d
      complex(rp), dimension(2) :: rac
      complex(rp), dimension(9, 9) :: Lx, Ly, Lz
      real(rp) :: lz_loc
      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
      Lx(:, :) = L_x(:, :)
      Ly(:, :) = L_y(:, :)
      Lz(:, :) = L_z(:, :)

      ! Transforming them into the spherical harmonics coordinates
      call hcpx(Lx, 'cart2sph')
      call hcpx(Ly, 'cart2sph')
      call hcpx(Lz, 'cart2sph')

      ! Writing the L.S hamiltonian
      this%lsham(:, :, :) = cmplx(0.0d0, 0.0d0)
      do k = 1, this%charge%lattice%ntype
         sg = cmplx(0.5d0, 0.0d0)
         soc_p = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_p(1)*this%charge%lattice%symbolic_atoms(k)%potential%xi_p(2))
         soc_d = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1)*this%charge%lattice%symbolic_atoms(k)%potential%xi_d(2))
         ! Check if orbital polarization is enabled
         if (this%orb_pol) then
            rac = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1)*this%charge%lattice%symbolic_atoms(k)%potential%rac)
            lz_loc = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1)*this%charge%lattice%symbolic_atoms(k)%potential%lmom(3))
         else
            rac = 0.0_rp
            lz_loc = 0.0_rp
         end if

         prefac = 0.0_rp
         do i = 1, 9
            do j = 1, 9
               if (i >= 2 .and. i <= 4 .and. j >= 2 .and. j <= 4) prefac = sg*soc_p
               if (i >= 5 .and. i <= 9 .and. j >= 5 .and. j <= 9) prefac = sg*soc_d
               this%lsham(j, i, k) = this%lsham(j, i, k) + prefac*Lz(j, i) + Lz(j, i)*rac(1)*lz_loc ! H11
               this%lsham(j, i + 9, k) = this%lsham(j, i + 9, k) + prefac*(Lx(j, i) - i_unit*Ly(j, i)) ! H12
               this%lsham(j + 9, i, k) = this%lsham(j + 9, i, k) + prefac*(Lx(j, i) + i_unit*Ly(j, i)) ! H21
               this%lsham(j + 9, i + 9, k) = this%lsham(j + 9, i + 9, k) - prefac*Lz(j, i) - Lz(j, i)*rac(2)*lz_loc ! H22
               !write(50,*) ´ntype=´, k
               !write(51,*) ´ntype=´, k
               !write(50,´(18f10.6)´) real(this%lsham(:,:,k))
               !write(51,´(18f10.6)´) aimag(this%lsham(:,:,k))
            end do
         end do
      end do
   end subroutine build_lsham

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build the torque operator for the collinear case (SQA in z)
   !> based on the definition of PRM 2, 013801 (2018).
   !> Implemented by Ivan Miranda on 07/11/2023.
   !---------------------------------------------------------------------------
   subroutine torque_operator_collinear(this)
      !
      class(hamiltonian), intent(inout) :: this
      !
      ! Local variables
      integer :: i, j, k
      complex(rp) :: prefac, sg, soc_p, soc_d
      complex(rp), dimension(9, 9) :: Lx, Ly, Lz
      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
      Lx(:, :) = L_x(:, :)
      Ly(:, :) = L_y(:, :)
      Lz(:, :) = L_z(:, :)

      ! Transforming them into the spherical harmonics coordinates
      call hcpx(Lx, 'cart2sph')
      call hcpx(Ly, 'cart2sph')
      call hcpx(Lz, 'cart2sph')

      ! Now write the torque operator matrix
      this%tmat(:, :, :, :) = cmplx(0.0_rp, 0.0_rp)
      do k = 1, this%charge%lattice%ntype
         sg = cmplx(0.5_rp, 0.0_rp)
         soc_p = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_p(1)*this%charge%lattice%symbolic_atoms(k)%potential%xi_p(2))
         soc_d = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1)*this%charge%lattice%symbolic_atoms(k)%potential%xi_d(2))
         prefac = 0.0_rp
         do i = 1, 9
            do j = 1, 9
               if (i >= 2 .and. i <= 4 .and. j >= 2 .and. j <= 4) prefac = sg*soc_p
               if (i >= 5 .and. i <= 9 .and. j >= 5 .and. j <= 9) prefac = sg*soc_d
               ! build Tx
               this%tmat(j, i, 1, k) = this%tmat(j, i, 1, k) + prefac*i_unit*Ly(j, i)*2.0_rp ! Tx_11
               this%tmat(j, i + 9, 1, k) = this%tmat(j, i + 9, 1, k) - prefac*Lz(j, i)*2.0_rp*cone ! Tx_12
               this%tmat(j + 9, i, 1, k) = this%tmat(j + 9, i, 1, k) + prefac*Lz(j, i)*2.0_rp*cone ! Tx_21
               this%tmat(j + 9, i + 9, 1, k) = this%tmat(j + 9, i + 9, 1, k) - prefac*i_unit*Ly(j, i)*2.0_rp ! Tx_22
               ! build Ty
               this%tmat(j, i, 2, k) = this%tmat(j, i, 2, k) - prefac*i_unit*Lx(j, i)*2.0_rp ! Ty_11
               this%tmat(j, i + 9, 2, k) = this%tmat(j, i + 9, 2, k) + prefac*i_unit*Lz(j, i)*2.0_rp ! Ty_12
               this%tmat(j + 9, i, 2, k) = this%tmat(j + 9, i, 2, k) + prefac*i_unit*Lz(j, i)*2.0_rp ! Ty_21
               this%tmat(j + 9, i + 9, 2, k) = this%tmat(j + 9, i + 9, 2, k) + prefac*i_unit*Lx(j, i)*2.0_rp ! Ty_22
               ! build Tz
               this%tmat(j, i + 9, 3, k) = this%tmat(j, i + 9, 3, k) + prefac*(Lx(j, i) - i_unit*Ly(j, i))*2.0_rp*cone ! Tz_12
               this%tmat(j + 9, i, 3, k) = this%tmat(j + 9, i, 3, k) + prefac*(Lx(j, i) + i_unit*Ly(j, i))*2.0_rp*cmone ! Tz_21
            end do
         end do
      end do

   end subroutine torque_operator_collinear

   subroutine build_obarm(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      complex(rp), dimension(9, 9) :: obm0, obm1
      complex(rp), dimension(3) :: mom
      integer :: ntype ! Atom type index
      integer :: l, m ! Orbital index

      this%obarm = 0.d00

      do ntype = 1, this%lattice%ntype
         obm0 = cmplx(0.0d0); obm1 = cmplx(0.0d0)
         do m = 1, 9
            obm0(m, m) = this%lattice%symbolic_atoms(ntype)%potential%obx0(m)
            obm1(m, m) = this%lattice%symbolic_atoms(ntype)%potential%obx1(m)
         end do
         mom(:) = cmplx(this%lattice%symbolic_atoms(ntype)%potential%mom(:), 0.0d0)
         do m = 1, 9
            do l = 1, 9
               this%obarm(m, l, ntype) = obm0(m, l) + obm1(m, l)*mom(3)
               this%obarm(m + 9, l + 9, ntype) = obm0(m, l) - obm1(m, l)*mom(3)
               this%obarm(l, m + 9, ntype) = obm1(m, l)*mom(1) - i_unit*obm1(m, l)*mom(2)
               this%obarm(l + 9, m, ntype) = obm1(m, l)*mom(1) + i_unit*obm1(m, l)*mom(2)
            end do
         end do
         call hcpx(this%obarm(1:9, 1:9, ntype), 'cart2sph')
         call hcpx(this%obarm(10:18, 10:18, ntype), 'cart2sph')
         call hcpx(this%obarm(1:9, 10:18, ntype), 'cart2sph')
         call hcpx(this%obarm(10:18, 1:9, ntype), 'cart2sph')
      end do
   end subroutine build_obarm

   subroutine build_enim(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      complex(rp), dimension(9, 9) :: em0, em1
      complex(rp), dimension(9) :: ex0, ex1
      complex(rp), dimension(3) :: mom
      complex(rp) :: eu, ed
      integer :: ntype ! Atom type index
      integer :: l, m ! Orbital index

      this%enim = 0.0d0

      do ntype = 1, this%lattice%ntype
         em0 = cmplx(0.0d0); em1 = cmplx(0.0d0)
         do m = 1, 9
            eu = this%lattice%symbolic_atoms(ntype)%potential%cx(m, 1) - this%lattice%symbolic_atoms(ntype)%potential%cex(m, 1)
            ed = this%lattice%symbolic_atoms(ntype)%potential%cx(m, 2) - this%lattice%symbolic_atoms(ntype)%potential%cex(m, 2)
            ex0(m) = 0.5*(eu + ed)
            ex1(m) = 0.5*(eu - ed)
            em0(m, m) = ex0(m)
            em1(m, m) = ex1(m)
         end do
         mom(:) = cmplx(this%lattice%symbolic_atoms(ntype)%potential%mom(:), 0.0d0)
         do m = 1, 9
            do l = 1, 9
               this%enim(m, l, ntype) = em0(m, l) + em1(m, l)*mom(3)
               this%enim(m + 9, l + 9, ntype) = em0(m, l) - em1(m, l)*mom(3)
               this%enim(l, m + 9, ntype) = em1(m, l)*mom(1) - i_unit*em1(m, l)*mom(2)
               this%enim(l + 9, m, ntype) = em1(m, l)*mom(1) + i_unit*em1(m, l)*mom(2)
            end do
         end do
         call hcpx(this%enim(1:9, 1:9, ntype), 'cart2sph')
         call hcpx(this%enim(10:18, 10:18, ntype), 'cart2sph')
         call hcpx(this%enim(1:9, 10:18, ntype), 'cart2sph')
         call hcpx(this%enim(10:18, 1:9, ntype), 'cart2sph')

         if (this%local_axis) then
            this%enim_glob = this%enim
         end if
      end do
   end subroutine build_enim

   subroutine build_bulkham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ij, ji, nr, ia
      integer :: ntype

      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         !write(123, *)´bulkham´
         call this%chbar_nc(ia, nr, ino, ntype)
         do m = 1, nr
            do i = 1, 9
               do j = 1, 9
                  this%ee(j, i, m, ntype) = this%hmag(j, i, m, 4) + this%hmag(j, i, m, 3)        ! H0+Hz
                  this%ee(j + 9, i + 9, m, ntype) = this%hmag(j, i, m, 4) - this%hmag(j, i, m, 3)        ! H0-Hz
                  this%ee(j, i + 9, m, ntype) = this%hmag(j, i, m, 1) - i_unit*this%hmag(j, i, m, 2) ! Hx-iHy
                  this%ee(j + 9, i, m, ntype) = this%hmag(j, i, m, 1) + i_unit*this%hmag(j, i, m, 2) ! Hx+iHy
               end do ! end of orbital j loop
            end do ! end of orbital i loop
            !write(128, *) ´m=´, m, ´ntype= ´, ntype
            !write(128, ´(18f10.6)´) real(this%ee(:, :, m, ntype))
         end do ! end of neighbour number
         ! Hubbard U correction.
         ! Only implemented for spd-orbitals
         if (this%hubbardU_check .and. this%hubbardJ_check) then
            print *, 'Add Hubbard U+J correction onto on-site Hamiltonian.'
            do i = 1, 9
               do j = 1, 9
                  this%ee(i, j, 1, ntype) = this%ee(i, j, 1, ntype) + this%hubbard_pot(i, j, ntype)
                  this%ee(i + 9, j + 9, 1, ntype) = this%ee(i + 9, j + 9, 1, ntype) + this%hubbard_pot(i + 9, j + 9, ntype)
               end do
            end do
            ! print *, 'hub_v_pot test ', this%hubbard_v_pot
            if (this%hubbardV_check) then
               print *, 'Adding hubbard_v_pot onto hamiltonian'
               do i = 1, 9
                  do j = 1, 9
                     do m = 1, nr
                        this%ee(i,j,m,ntype) = this%ee(i,j,m,ntype) + this%hubbard_v_pot(i,j,m,ntype)
                        this%ee(i+9,j+9,m,ntype) = this%ee(i+9,j+9,m,ntype) + this%hubbard_v_pot(i+9,j+9,m,ntype)
                     end do
                  end do
               end do
            end if
         end if

         if (this%hoh) then
            call this%build_obarm()
            call this%build_enim()
            do m = 1, nr
               ji = 0
               if (m > 1) then
                  ja = this%charge%lattice%nn(ia, m)
                  if (ja .ne. 0) then
                     ji = this%charge%lattice%iz(ja)
                  end if
               else
                  ji = this%charge%lattice%iz(ia)
               end if
               ! Check if neighbour ´m´ exists for atom ´ntype´, otherwise fill HoH Hamiltonian with zeros.
               if (ji > 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%ee(:, :, m, ntype), 18, this%obarm(:, :, ji), 18, czero, this%eeo(:, :, m, ntype), 18)
                  call zgemm('n', 'c', 18, 18, 18, cone, this%eeo(:, :, m, ntype), 18, this%ee(:, :, m, ntype), 18, czero, this%eeoee(:, :, m, ntype), 18)
               else
                  this%eeo(:, :, m, ntype) = 0.0d0
               end if
               !write(*,*) ´m=´, m
               !write(*,´(18f10.6)´) real(this%eeo(:,:,m,ntype))
               !write(*,*) ´ee´, m
               !write(*,´(18f10.6)´) real(this%ee(:,:,m,ntype))
            end do
         end if
      end do ! end of atom type number
      if (this%local_axis) then
         this%ee_glob = this%ee
         if (this%hoh) this%eeo_glob = this%eeo
      end if
      close (128)
   end subroutine build_bulkham

   subroutine build_locham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: it, ino, nr, nlim, m, i, j, ja, ji

      call g_timer%start('build local hamiltonian')
    !!$omp parallel do private(nlim, nr, ino, m, i, j, ji, ja, this)
      do nlim = 1, this%charge%lattice%nmax
         nr = this%charge%lattice%nn(nlim, 1) ! Number of neighbours considered
         ino = this%charge%lattice%num(nlim)
         call this%chbar_nc(nlim, nr, ino, nlim)
         do m = 1, nr
            do i = 1, 9
               do j = 1, 9
                  this%hall(j, i, m, nlim) = this%hmag(j, i, m, 4) + this%hmag(j, i, m, 3) ! H0+Hz
                  this%hall(j + 9, i + 9, m, nlim) = this%hmag(j, i, m, 4) - this%hmag(j, i, m, 3) ! H0-Hz
                  this%hall(j, i + 9, m, nlim) = this%hmag(j, i, m, 1) - i_unit*this%hmag(j, i, m, 2) ! Hx-iHy
                  this%hall(j + 9, i, m, nlim) = this%hmag(j, i, m, 1) + i_unit*this%hmag(j, i, m, 2) ! Hx+iHy
               end do
            end do
         end do
         if (this%hoh) then
            call this%build_obarm()
            call this%build_enim()
            do m = 1, nr
               ji = 0
               if (m > 1) then
                  ja = this%charge%lattice%nn(nlim, m)
                  if (ja .ne. 0) then
                     ji = this%charge%lattice%iz(ja)
                  end if
               else
                  ji = this%charge%lattice%iz(nlim)
               end if
               ! Check if neighbour ´m´ exists for atom ´nlim´, otherwise fill HoH Hamiltonian with zeros.
               if (ji > 0) then
                  call zgemm('n', 'n', 18, 18, 18, cone, this%hall(1, 1, m, nlim), 18, this%obarm(1, 1, ji), 18, czero, this%hallo(1, 1, m, nlim), 18)
               else
                  this%hallo(:, :, m, nlim) = 0.0d0
               end if
            end do
         end if
      end do
    !!$omp end parallel do
      if (this%local_axis) then
         this%hall_glob = this%hall
         if (this%hoh) this%hallo_glob = this%hallo
      end if
      call g_timer%stop('build local hamiltonian')
   end subroutine build_locham

   subroutine rs2pao(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      real(rp), dimension(3) :: rij, rijtest
      integer :: i, j, k, l, idxi, idxj, idxk, itype, ino, ja, jo, ji, nr, ia, iia, jja, ipao, jpao
      integer :: jj, jt, max_orbital, n_atoms
      integer :: ntype, iostat1, iostat2, iostatus
      real(rp), dimension(3) :: vet, vetpao, idx
      real(rp), dimension(3, 3) :: a_inv
      complex(rp), dimension(18, 18) :: dum
      n_atoms = this%charge%lattice%ntype
      max_orbital = 9

      open (unit=92, file='rs2paoham.dat', action='write', iostat=iostatus, status='replace')
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         do k = 1, nr
            jj = this%charge%lattice%nn(ia, k)
            !write(123, *)´ia, ii´, ia, m, this%charge%lattice%nn(ia, m)
            if (k == 1) then
               jj = ia
            end if
            if (jj /= 0) then
               rij(:) = this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, jj)

               rijtest(:) = 0.0d0
               do idxi = -5, 5
                  do idxj = -5, 5
                     do idxk = -5, 5
                        rijtest(:) = this%charge%lattice%cr(:, ia) - (this%charge%lattice%cr(:, this%charge%lattice%iz(jj)) &
                                                                      + idxi*this%charge%lattice%a(:, 1) + idxj*this%charge%lattice%a(:, 2) + idxk*this%charge%lattice%a(:, 3))
                        if (norm2(rij(:) - rijtest(:)) < 1.0d-3) then
                           idx(:) = [idxi, idxj, idxk]
                        end if
                     end do
                  end do
               end do
               !a_inv(:,:) = inverse_3x3(this%charge%lattice%a)
               !idx(:) = matmul(a_inv(:,:),b(:))

               !write(*,*) b(:), ´b´
               !write(*,*) matmul(this%charge%lattice%a(:,:),idx(:))

               !call hcpx(this%ee(1:9,1:9,k,ntype), ´sph2cart´)
               !call hcpx(this%ee(1:9,10:18,k,ntype), ´sph2cart´)
               !call hcpx(this%ee(10:18,1:9,k,ntype), ´sph2cart´)
               !call hcpx(this%ee(10:18,10:18,k,ntype), ´sph2cart´)
               !if(this%hoh)then
               !  !call zgemm(´n´,´c´,18,18,18,cone,this%eeo(:,:,k,ntype),18,this%ee(:,:,k,ntype),18,cone,dum(:,:),18)
               !  dum(:,:) = 0.0d0
               !  do i=1,nr
               !    dum(:,:) = dum(:,:) + this%eeoee(:,:,i,ntype)
               !  end do
               !  this%ee(:,:,k,ntype) = this%ee(:,:,k,ntype) - dum(:,:)
               !end if
               if (k == 1) this%ee(:, :, k, ntype) = this%ee(:, :, k, ntype) + this%lsham(:, :, ntype) !+ this%enim(:,:,ntype)
               do i = 1, 18
                  do j = 1, 18
                     ipao = 0; jpao = 0
                     call site2orb(i, ia, ipao, n_atoms, max_orbital)
                     call site2orb(j, this%charge%lattice%iz(jj), jpao, n_atoms, max_orbital)
                     write (92, '(3I4,2I7,2F22.14)') int(idx(:)), ipao, jpao, real(this%ee(i, j, k, ntype))*ry2ev, aimag(this%ee(i, j, k, ntype))*ry2ev
                  end do
               end do
            end if
         end do
      end do
      close (92)
   end subroutine rs2pao

   subroutine build_from_paoflow_opt(this)
      implicit none
      type hamData
         integer :: idxi, idxj, idxk
         integer :: orbl, orbm
         real :: dumre, dumcmplx
      end type hamData
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia, iia, jja
      integer :: jj, jt, orbl, orbm, idxi, idxj, idxk
      integer :: ntype, iostat1, iostat2, iostatus, n_atoms, max_orbital, numLines
      real(rp), dimension(3) :: vet, vetpao
      real(rp) :: dumre, dumcmplx
      integer, dimension(maxval(this%charge%lattice%nn(:, 1)) + 1, 3) :: idxup, idxdw, idx
      type(hamData) :: ham
      type(hamData), allocatable :: hamArray(:)

      n_atoms = this%charge%lattice%ntype
      max_orbital = 9
      numLines = countLines('paoham.dat')
      allocate (hamArray(numLines))

      ! Reading Hamiltonian data once
      open (unit=92, file='paoham.dat', action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file containing the paoflow Hamiltonian not found', __FILE__, __LINE__)
      end if
      do i = 1, numLines
         read (92, *, iostat=iostat1) hamArray(i)%idxi, hamArray(i)%idxj, hamArray(i)%idxk, &
            hamArray(i)%orbl, hamArray(i)%orbm, hamArray(i)%dumre, hamArray(i)%dumcmplx
      end do
      close (92)
      idx = 0
      write (*, *) 'PAOFLOW Hamiltonian has been read'
      call g_timer%start('Hamiltonian allocation')
      !$omp parallel do private(ntype, ia, ino, nr, jj, k, l, ham, vet, vetpao, i, j, iia, jja, idx) shared(this, numLines, hamArray, n_atoms, max_orbital)
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         ino = this%charge%lattice%num(ia)
         nr = this%charge%lattice%nn(ia, 1)
         write (*, *) ia, nr, n_atoms, max_orbital
         do k = 1, nr
            jj = this%charge%lattice%nn(ia, k)
            if (k == 1) then
               jj = ia
            end if
            if (jj /= 0) then
               vet(:) = (this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, jj))!*this%charge%lattice%alat
               do l = 1, numLines
                  ham = hamArray(l)
                  call orb2site(ham%orbl, i, iia, n_atoms, max_orbital)
                  call orb2site(ham%orbm, j, jja, n_atoms, max_orbital)
                  if (iia == ia) then
                     vetpao(:) = this%charge%lattice%cr(:, iia) - (this%charge%lattice%cr(:, jja) &
                                                                   + ham%idxi*this%charge%lattice%a(:, 1) + ham%idxj*this%charge%lattice%a(:, 2) + ham%idxk*this%charge%lattice%a(:, 3))!*this%charge%lattice%alat
                     if (norm2(vet(:) - vetpao(:)) < 1.0d-3) then
                        idx(k, :) = [ham%idxi, ham%idxj, ham%idxk]
                        this%ee(i, j, k, ntype) = cmplx(ham%dumre, ham%dumcmplx)/13.605703976
                        ! if(ntype==1.or.ntype==2)then
                        !   if(i==j.and.k==1.and.i>=5.and.i<=9)then
                        !     this%ee(i, j, k, ntype) = cmplx(ham%dumre-1.113 , ham%dumcmplx)/13.605703976
                        !   else if(i==j.and.k==1.and.i>=14.and.i<=18)then
                        !     this%ee(i, j, k, ntype) = cmplx(ham%dumre+1.113 , ham%dumcmplx)/13.605703976
                        !   end if
                        ! end if
                     end if
                  end if
               end do
            end if
            write (128, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (129, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (128, '(18f10.6)') real(this%EE(1:18, 1:18, k, ntype))*13.605703976
            write (129, '(18f10.6)') aimag(this%EE(1:18, 1:18, k, ntype))*13.605703976
            write (128, *) sum(real(this%ee(:, :, k, ntype)))
            write (129, *) sum(real(this%ee(:, :, k, ntype)))
         end do
      end do
      !$omp end parallel do
      deallocate (hamArray)
      call g_timer%stop('Hamiltonian allocation')
   end subroutine build_from_paoflow_opt

   subroutine build_from_paoflow(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia, iia, jja
      integer :: jj, jt, orbl, orbm, idxi, idxj, idxk
      integer :: ntype, iostat1, iostat2, iostatus, n_atoms, max_orbital
      real(rp), dimension(3) :: vet, vetpao, cri_dir, crj_dir, cri_cart, crj_cart
      integer, dimension(maxval(this%charge%lattice%nn(:, 1)) + 1, 3) :: idxup, idxdw, idx
      real(rp) :: dumre, dumcmplx

      n_atoms = this%charge%lattice%ntype
      max_orbital = 9
      open (unit=90, file='paoup.dat', action='read', iostat=iostatus, status='old')
      open (unit=91, file='paodw.dat', action='read', iostat=iostatus, status='old')
      open (unit=92, file='paoham.dat', action='read', iostat=iostatus, status='old')

      if (iostatus /= 0) then
         call g_logger%fatal('file containing the paoflow Hamiltonian not found', __FILE__, __LINE__)
      end if
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         write (*, *) ia, nr, n_atoms, max_orbital
         do k = 1, nr
            jj = this%charge%lattice%nn(ia, k)
            !write(123, *)´ia, ii´, ia, m, this%charge%lattice%nn(ia, m)
            if (k == 1) then
               jj = ia
            end if
            if (jj /= 0) then
               !cri_cart(:) = this%charge%lattice%cr(:, ia)
               !crj_cart(:) = this%charge%lattice%cr(:, jj)

               !cri_dir(:) = cartesian_to_direct(this%charge%lattice%a,cri_cart)
               !crj_dir(:) = cartesian_to_direct(this%charge%lattice%a,crj_cart)

               !vet(:) = cri_dir(:) - crj_dir(:)
               vet(:) = (this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, jj))!*this%charge%lattice%alat
               do
                  read (92, *, iostat=iostat1) idxi, idxj, idxk, orbl, orbm, dumre, dumcmplx
                  if (iostat1 /= 0) then
                     exit
                  end if
                  if (orbl <= n_atoms*max_orbital) then
                     i = modulo(orbl - 1, max_orbital) + 1
                     iia = int((orbl - 1)/max_orbital) + 1
                  else
                     i = modulo(orbl - 1, max_orbital) + 10
                     iia = int((orbl - 1 - n_atoms*max_orbital)/max_orbital) + 1
                  end if

                  if (orbm <= n_atoms*max_orbital) then
                     j = modulo(orbm - 1, max_orbital) + 1
                     jja = int((orbm - 1)/max_orbital) + 1
                  else
                     j = modulo(orbm - 1, max_orbital) + 10
                     jja = int((orbm - 1 - n_atoms*max_orbital)/max_orbital) + 1
                  end if

                  !cri_cart(:) = this%charge%lattice%cr(:, iia)
                  !crj_cart(:) = this%charge%lattice%cr(:, jja)

                  !cri_dir(:) = cartesian_to_direct(this%charge%lattice%a,cri_cart)
                  !crj_dir(:) = cartesian_to_direct(this%charge%lattice%a,crj_cart)
                  if (iia == ia) then
                     !vetpao(:) = cri_dir(:) - (crj_dir(:) + [idxi,idxj,idxk])

                     vetpao(:) = this%charge%lattice%cr(:, iia) - (this%charge%lattice%cr(:, jja) &
                                                                   + idxi*this%charge%lattice%a(:, 1) + idxj*this%charge%lattice%a(:, 2) + idxk*this%charge%lattice%a(:, 3))!*this%charge%lattice%alat
                     if (norm2(vet(:) - vetpao(:)) < 1.0d-3) then
                        idx(k, :) = [idxi, idxj, idxk]
                        this%ee(i, j, k, ntype) = cmplx(dumre, dumcmplx)/13.605703976
                     end if
                  end if
               end do
!          do
!            read(91,*,iostat=iostat2) idxi, idxj, idxk, orbl, orbm, dumre, dumcmplx
!            if (iostat2 /= 0) then
!              exit
!            end if
!            iia = orb_to_site(orbl,9)
!            jja = orb_to_site(orbm,9)
!
!            cri_cart(:) = this%charge%lattice%cr(:, iia)
!            crj_cart(:) = this%charge%lattice%cr(:, jja)
!
!            cri_dir(:) = cartesian_to_direct(this%charge%lattice%a,cri_cart)
!            crj_dir(:) = cartesian_to_direct(this%charge%lattice%a,crj_cart)
!            if(iia==ia)then
!              vetpao(:) = cri_dir(:) - (crj_dir(:) + [idxi,idxj,idxk])
!
!!              vetpao(:) = this%charge%lattice%cr(:,iia) - this%charge%lattice%cr(:,jja) &
!!                         + (idxi*this%charge%lattice%a(:,1) + idxj*this%charge%lattice%a(:,2) + idxk*this%charge%lattice%a(:,3))!*this%charge%lattice%alat
!              if(norm2(vet(:)-vetpao(:))<1.0d-3)then
!                idxdw(k,:) = [idxi,idxj,idxk]
!                this%ee(orbl+9-(iia-1)*9, orbm+9-(jja-1)*9,k,ntype) = cmplx(dumre,dumcmplx)/13.605703976
!              end if
!            end if
!          end do
            end if
            write (128, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (128, '(18f10.6)') real(this%EE(1:18, 1:18, k, ntype))!*13.605703976
            write (128, *) sum(real(this%ee(:, :, k, ntype)))
            !write(129,*)´m=´,k, ´Atom=´, jj, ´Coordinates=´, this%charge%lattice%cr(:, jj), ´Ntype=´,ntype, ´Index=´, idx(k,:)
            !write(129,´(9f10.6)´) real(this%EE(10:18,1:9,k,ntype))*13.605703976
            !write(129,*) sum(real(this%ee(:,:,k,ntype)))
            rewind (90)
            rewind (91)
            rewind (92)
         end do
      end do
   end subroutine build_from_paoflow

   subroutine ham0m_nc(this, it, jt, vet, hhh)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: it, jt ! Type of atom i and j
      real(rp), dimension(3), intent(in) :: vet
      real(rp), dimension(9, 9), intent(in) :: hhh
      ! Local Variables
      integer :: i, j, ilm, jlm, m
      complex(rp), dimension(3) :: cross
      complex(rp), dimension(9, 9) :: hhhc
      complex(rp), dimension(this%charge%lattice%ntype, 3) :: momc
      complex(rp) :: dot
      real(rp) :: vv

      this%hhmag(:, :, :) = 0.0d0

      vv = norm2(vet)

      ! Real to complex
      dot = cmplx(dot_product(this%charge%lattice%symbolic_atoms(it)%potential%mom, this%charge%lattice%symbolic_atoms(jt)%potential%mom), kind=kind(0.0d0))
      do i = 1, this%charge%lattice%ntype
         do j = 1, 3
            momc(i, j) = cmplx(this%charge%lattice%symbolic_atoms(i)%potential%mom(j), kind=kind(0.0d0))
         end do
      end do
      cross = cmplx(cross_product(this%charge%lattice%symbolic_atoms(it)%potential%mom, this%charge%lattice%symbolic_atoms(jt)%potential%mom), kind=kind(0.0d0))
      hhhc(:, :) = cmplx(hhh(:, :), kind=kind(0.0d0))

      do ilm = 1, 9
         do jlm = 1, 9
            this%hhmag(ilm, jlm, 4) = &
               this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm) + &
               this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*dot
         end do
      end do

!    do ilm=1, 9
!      write(123, ´(9f10.6)´) (real(this%hhmag(ilm, jlm, 4)), jlm=1, 9)
!    end do

      if (vv <= 0.01d0) then
         do ilm = 1, 9
            if (this%hoh) then
               this%hhmag(ilm, ilm, 4) = this%hhmag(ilm, ilm, 4) + this%charge%lattice%symbolic_atoms(it)%potential%cex0(ilm)
            else
               this%hhmag(ilm, ilm, 4) = this%hhmag(ilm, ilm, 4) + this%charge%lattice%symbolic_atoms(it)%potential%cx0(ilm)
            end if
         end do
      end if

      do m = 1, 3
         do jlm = 1, 9
            do ilm = 1, 9
               this%hhmag(ilm, jlm, m) = &
                  (this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm))*momc(it, m) + &
                  (this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm))*momc(jt, m) + &
                  i_unit*this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*cross(m)
            end do
         end do
      end do

      if (vv > 0.01d0) return
      do m = 1, 3
         do ilm = 1, 9
            if (this%hoh) then
               this%hhmag(ilm, ilm, m) = this%hhmag(ilm, ilm, m) + this%charge%lattice%symbolic_atoms(it)%potential%cex1(ilm)*momc(it, m)
            else
               this%hhmag(ilm, ilm, m) = this%hhmag(ilm, ilm, m) + this%charge%lattice%symbolic_atoms(it)%potential%cx1(ilm)*momc(it, m)
            end if
         end do
      end do

      !do m=1, 3
      !  write(123, *)´m=´, m
      !  do ilm=1, 9
      !    write(123, ´(9f10.6)´) (real(this%hhmag(ilm, jlm, m)), jlm=1, 9)
      !  end do
      !end do
   end subroutine ham0m_nc

   subroutine chbar_nc(this, ia, nr, ino, ntype)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: ia ! Atom number in clust
      integer, intent(in) :: nr ! Number of neighbours considered
      integer, intent(in) :: ino ! Atom bravais type of ia
      integer, intent(in) :: ntype ! Atom type
      ! Local variables
      real(rp) :: r2
      real(rp), dimension(3, size(this%charge%lattice%cr(1, :))) :: cralat ! Clust position times the lattice constant
      real(rp), dimension(3) :: vet
      real(rp), dimension(9, 9) :: hhh
      integer :: i, j, k, l, m, n, it, jt, jj, dummy
      integer :: ni, mdir
      integer :: kk ! Clust size number

      this%hmag(:, :, :, :) = 0.0d0

      r2 = this%charge%lattice%r2
      cralat(:, :) = this%charge%lattice%cr(:, :)*this%charge%lattice%alat
      kk = this%charge%lattice%kk

      call this%charge%lattice%clusba(r2, cralat, ia, kk, kk, dummy)

      !do m=1, nr
      !  print ´(9f10.6)´, real(this%charge%lattice%sbar(:, :, m, ino))
      !end do
      it = this%charge%lattice%iz(ia)
      do m = 1, nr
         jj = this%charge%lattice%nn(ia, m)
         !write(123, *)´ia, ii´, ia, m, this%charge%lattice%nn(ia, m)
         if (m == 1) then
            jj = ia
         end if
         if (jj /= 0) then
            jt = this%charge%lattice%iz(jj)
            vet(:) = (this%charge%lattice%cr(:, jj) - this%charge%lattice%cr(:, ia))*this%charge%lattice%alat
            !write(123, ´(3f10.6)´) vet(:)
            !write(123, ´(3f10.6)´) this%charge%lattice%sbarvec(:, m)
            !write(123, ´(a, 3i4, 3f10.6)´) ´nn ´, IA, m, JJ, VET(:)
            call this%hmfind(vet, nr, hhh, m, ia, m, ni, ntype)
            if (ni == 0) then
               this%charge%lattice%nn(ia, m) = 0
            end if
            call this%ham0m_nc(it, jt, vet, hhh)
            do mdir = 1, 4
               call hcpx(this%hhmag(:, :, mdir), 'cart2sph')
               this%hmag(:, :, m, mdir) = this%hhmag(:, :, mdir)
            end do
         end if
      end do
      !do m=1, nr
      !  write(123, *)´m=´, m
      !  do mdir=1, 4
      !    write(123, *)´mdir=´, mdir
      !    do i=1, 9
      !      write(123, ´(9f10.4)´)(real(this%hmag(i, j, m, mdir)), j=1, 9)
      !    end do
      !  end do
      !end do
   end subroutine chbar_nc

   subroutine hmfind(this, vet, nr, hhh, m, ia, jn, ni, ntype)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: ntype ! Atom type
      integer, intent(in) :: m ! Number of the given neighbour
      integer, intent(in) :: ia ! Atom number in clust
      integer, intent(in) :: jn ! ?
      integer, intent(in) :: nr ! Number of neighbours
      real(rp), dimension(3), intent(in) :: vet
      ! Output
      integer, intent(out) :: ni
      real(rp), dimension(9, 9), intent(inout) :: hhh
      ! Local variables
      real(rp) :: a1, a2, a3, aaa, eps
      integer :: i, ilm, jlm

      eps = 0.0001d0
      ni = 1
      a1 = 0.0d0
      a2 = 0.0d0
      a3 = 0.0d0
      aaa = 0.0d0
      do i = 1, nr
         !write(123, ´(a, i4, 3f10.4)´)´i´, i, this%charge%lattice%sbarvec(:, i)
         a1 = (vet(1) - this%charge%lattice%sbarvec(1, i))
         a2 = (vet(2) - this%charge%lattice%sbarvec(2, i))
         a3 = (vet(3) - this%charge%lattice%sbarvec(3, i))
         aaa = a1**2 + a2**2 + a3**2
         if (aaa < eps) goto 1000
      end do
      write (*, '(1x, a, i4, a, i4, a, 3f10.6)') ' Error in hamiltonian%hmfind: Neighbour vector not found for atom', ia, &
         ' neighbour', jn, 'vector', vet(:)

      ni = 0
1000  continue
      do ilm = 1, 9
         do jlm = 1, 9
            hhh(ilm, jlm) = real(this%charge%lattice%sbar(jlm, ilm, m, this%charge%lattice%num(ia)))
         end do
      end do

      !do ilm = 1, 9
      !    write(123, ´(9f10.6)´)(hhh(ilm, jlm), jlm=1, 9)
      !end do
   end subroutine hmfind

   subroutine orb2site(orb, i_out, ia_out, n_atoms, max_orbital)
      integer, intent(in) :: orb, n_atoms, max_orbital
      integer, intent(out) :: i_out, ia_out

      if (orb <= n_atoms*max_orbital) then
         i_out = modulo(orb - 1, max_orbital) + 1
         ia_out = int((orb - 1)/max_orbital) + 1
      else
         i_out = modulo(orb - 1, max_orbital) + 10
         ia_out = int((orb - 1 - n_atoms*max_orbital)/max_orbital) + 1
      end if
   end subroutine orb2site

   subroutine site2orb(i_in, ia_in, orb_out, n_atoms, max_orbital)
      integer, intent(in) :: i_in, ia_in, n_atoms, max_orbital
      integer, intent(out) :: orb_out

      if (i_in <= max_orbital) then
         orb_out = (ia_in - 1)*max_orbital + i_in
      else
         orb_out = (ia_in - 1)*max_orbital + i_in + (n_atoms - 1)*max_orbital
      end if
   end subroutine site2orb

   ! Rotate Hamiltonian to local axis
   subroutine rotate_to_local_axis(this, m_loc)
      use math_mod, only: rotmag_loc
      class(hamiltonian), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: m_loc

      ! Local variables
      integer :: sdim
      ! Rotate Hamiltonian to local axis if wanted
      if (this%local_axis) then
         sdim = product(shape(this%hall))/18/18
         call rotmag_loc(this%hall, this%hall_glob, sdim, m_loc)
         sdim = product(shape(this%ee))/18/18
         call rotmag_loc(this%ee, this%ee_glob, sdim, m_loc)
         if (this%hoh) then
            sdim = product(shape(this%eeo))/18/18
            call rotmag_loc(this%eeo, this%eeo_glob, sdim, m_loc)
            sdim = product(shape(this%hallo))/18/18
            call rotmag_loc(this%hallo, this%hallo_glob, sdim, m_loc)
            sdim = product(shape(this%enim))/18/18
            call rotmag_loc(this%enim, this%enim_glob, sdim, m_loc)
         end if
      end if
   end subroutine rotate_to_local_axis

   ! Rotate Hamiltonian back from local axis
   subroutine rotate_from_local_axis(this, m_loc)
      use math_mod, only: rotmag_loc
      class(hamiltonian), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: m_loc

      ! Local variables
      integer :: sdim
      ! Rotate Hamiltonian to local axis if wanted
      if (this%local_axis) then
         this%hall = this%hall_glob
         this%ee = this%ee_glob
         if (this%hoh) then
            this%eeo = this%eeo_glob
            this%hallo = this%hallo_glob
            this%enim = this%enim_glob
         end if
      end if
   end subroutine rotate_from_local_axis
end module hamiltonian_mod
