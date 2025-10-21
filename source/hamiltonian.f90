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
      ! !> Hamiltonian built in chbar_nc (description to be improved)
      ! complex(rp), dimension(:, :, :, :), allocatable :: hmag
      ! !> Hamiltonian built in ham0m_nc (description to be improved
      ! complex(rp), dimension(:, :, :), allocatable :: hhmag
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
      !> Velocity operators
      complex(rp), dimension(:, :, :, :), allocatable :: v_a, v_b, js_a, jl_a, vo_a, vo_b, jso_a, jlo_a
      character(len=10) :: js_alpha, jl_alpha
      real(rp), dimension(3) :: v_alpha, v_beta
      real(rp), dimension(:), allocatable :: velocity_scale
      !> Sparse Real Space Hamiltonian
      complex(rp), dimension(:, :), allocatable :: h_sparse
      !> Spin-spiral wave vector q
      real(rp), dimension(3) :: q_ss
      !
      !!! Testing Gershgorin bounds for later implementation
      !!! !> Upper Gershgorin bound
      !!! real(rp) :: g_max
      !!! !> Lower Gershgorin bound
      !!! real(rp) :: g_min

      !> On-site potential for LDA+U+J correction
      real(rp), dimension(:,:,:), allocatable :: hubbard_u_pot
      real(rp), dimension(:,:,:,:), allocatable :: hubbard_v_pot

      !> New improved Hubbard parameters that should work for both bulk and impurity
      real(rp), dimension(:,:), allocatable :: hubbard_u_general
      real(rp), dimension(:,:), allocatable :: hubbard_j_general
      logical :: hubbard_u_general_check = .False.

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
      !> Logical variables to include Hubbard U & J
      logical :: hubbardU_check = .false.
      logical :: hubbardJ_check = .false.
      logical :: hubbardV_check = .false.
      !> Logical variable to include self-consistent calculation of Hubbard U.
      logical :: hubbardU_sc_check = .false.
      !> Which atoms and orbitals to include self-consistent Hubbard U correction. (0 no correction. 1 include correction)
      integer, dimension(:,:), allocatable :: hubbard_u_sc
      !> Logical variable for including Hubbard correction U on impurities
      logical :: hubbardU_impurity_check = .false.
      !> Hubbard parameters for impurities
      real(rp), dimension(:,:), allocatable :: hubbard_u_impurity, hubbard_j_impurity ! (impurity type, l) with l = 1,2,3,4 : s,p,d,f orbitals
      !> Slater integral array for impurities
      real(rp), dimension(:,:,:), allocatable :: F_impurity
      !> On-site potential for LDA+U correction for impurities
      real(rp), dimension(:,:,:), allocatable :: hubbard_pot_impurity
   contains
      procedure :: build_lsham
      procedure :: build_bulkham
      procedure :: build_locham
      procedure :: build_obarm
      procedure :: build_enim
      procedure :: build_from_paoflow
      procedure :: build_from_paoflow_opt
      procedure :: build_realspace_velocity_operators
      procedure :: build_realspace_spin_operators
      procedure :: build_realspace_spin_torque_operators
      procedure :: build_realspace_orbital_velocity_operators
      procedure :: build_realspace_orbital_torque_operators
      procedure :: block_to_sparse
      procedure :: torque_operator_collinear
      procedure :: rs2pao
      procedure :: rs2txt
      procedure :: chbar_nc
      procedure :: ham0m_nc
      procedure :: hmfind
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: rotate_to_local_axis
      procedure :: rotate_from_local_axis
      procedure :: calculate_hubbard_u_potential_general
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
      !if (allocated(this%hmag)) call g_safe_alloc%deallocate('hamiltonian.hmag', this%hmag)
      !if (allocated(this%hhmag)) call g_safe_alloc%deallocate('hamiltonian.hhmag', this%hhmag)
      if (allocated(this%hall)) call g_safe_alloc%deallocate('hamiltonian.hall', this%hall)
      if (allocated(this%eeo)) call g_safe_alloc%deallocate('hamiltonian.eeo', this%eeo)
      if (allocated(this%eeoee)) call g_safe_alloc%deallocate('hamiltonian.eeoee', this%eeoee)
      if (allocated(this%hallo)) call g_safe_alloc%deallocate('hamiltonian.hallo', this%hallo)
      if (allocated(this%obarm)) call g_safe_alloc%deallocate('hamiltonian.obarm', this%obarm)
      if (allocated(this%enim)) call g_safe_alloc%deallocate('hamiltonian.enim', this%enim)
      if (allocated(this%ee_glob)) call g_safe_alloc%deallocate('hamiltonian.ee_glob', this%ee_glob)
      if (allocated(this%eeo_glob)) call g_safe_alloc%deallocate('hamiltonian.eeo_glob', this%eeo_glob)
      if (allocated(this%enim_glob)) call g_safe_alloc%deallocate('hamiltonian.enim_glob', this%enim_glob)
      if (allocated(this%v_a)) call g_safe_alloc%deallocate('hamiltonian.v_a', this%v_a)
      if (allocated(this%v_b)) call g_safe_alloc%deallocate('hamiltonian.v_b', this%v_b)
      if (allocated(this%vo_a)) call g_safe_alloc%deallocate('hamiltonian.vo_a', this%vo_a)
      if (allocated(this%vo_b)) call g_safe_alloc%deallocate('hamiltonian.vo_b', this%vo_b)
      if (allocated(this%jso_a)) call g_safe_alloc%deallocate('hamiltonian.jso_a', this%jso_a)
      if (allocated(this%jlo_a)) call g_safe_alloc%deallocate('hamiltonian.jlo_a', this%jlo_a)
      if (allocated(this%h_sparse)) call g_safe_alloc%deallocate('hamiltonian.h_sparse', this%h_sparse)
      if (allocated(this%velocity_scale)) call g_safe_alloc%deallocate('hamiltonian.velocity_scale', this%velocity_scale)
#else
      if (allocated(this%lsham)) deallocate (this%lsham)
      if (allocated(this%tmat)) deallocate (this%tmat)
      if (allocated(this%ee)) deallocate (this%ee)
      if (allocated(this%eeo)) deallocate (this%eeo)
      if (allocated(this%eeoee)) deallocate (this%eeoee)
      ! if (allocated(this%hmag)) deallocate (this%hmag)
      ! if (allocated(this%hhmag)) deallocate (this%hhmag)
      if (allocated(this%hall)) deallocate (this%hall)
      if (allocated(this%hallo)) deallocate (this%hallo)
      if (allocated(this%obarm)) deallocate (this%obarm)
      if (allocated(this%enim)) deallocate (this%enim)
      if (allocated(this%ee_glob)) deallocate (this%ee_glob)
      if (allocated(this%eeo_glob)) deallocate (this%eeo_glob)
      if (allocated(this%enim_glob)) deallocate (this%enim_glob)
      if (allocated(this%v_a)) deallocate(this%v_a)
      if (allocated(this%v_b)) deallocate(this%v_b)
      if (allocated(this%vo_a)) deallocate(this%vo_a)
      if (allocated(this%vo_b)) deallocate(this%vo_b)
      if (allocated(this%jso_a)) deallocate(this%jso_a)
      if (allocated(this%jlo_a)) deallocate(this%jlo_a)
      if (allocated(this%h_sparse)) deallocate(this%h_sparse)
      if (allocated(this%velocity_scale)) deallocate(this%velocity_scale)
      if (allocated(this%hubbard_u)) deallocate (this%hubbard_u)
      if (allocated(this%hubbard_j)) deallocate (this%hubbard_j)
      if (allocated(this%hubbard_v)) deallocate (this%hubbard_v)
      if (allocated(this%hub_u_sort)) deallocate (this%hub_u_sort)
      if (allocated(this%hub_j_sort)) deallocate (this%hub_j_sort)
      if (allocated(this%uj_orb)) deallocate (this%uj_orb)
      if (allocated(this%orbs_v)) deallocate (this%orbs_v)
      if (allocated(this%orbs_v_num)) deallocate (this%orbs_v_num)
      if (allocated(this%F)) deallocate(this%F)
      if (allocated(this%hubbard_u_pot)) deallocate (this%hubbard_u_pot)
      if (allocated(this%hubbard_v_pot)) deallocate (this%hubbard_v_pot)
      if (allocated(this%hubbard_u_sc)) deallocate (this%hubbard_u_sc)
      if (allocated(this%hubbard_u_impurity)) deallocate (this%hubbard_u_impurity)
      if (allocated(this%hubbard_j_impurity)) deallocate (this%hubbard_j_impurity)
      if (allocated(this%F_impurity)) deallocate (this%F_impurity)
      if (allocated(this%hubbard_pot_impurity)) deallocate (this%hubbard_pot_impurity)
      if (allocated(this%hubbard_u_general)) deallocate (this%hubbard_u_general)
      if (allocated(this%hubbard_j_general)) deallocate (this%hubbard_j_general)
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
      integer :: iostatus, funit, i, j, li, lj, k, max_orbs, length, cntr, na, l, ios
      ! integer, dimension(this%lattice%nrec, this%lattice%nrec) :: orbs_v_num ! (atom1, atom2) = number of mutual orbs e.g. hubbard_v(1, 2, [p,d], [s,p]) -> orbs_v_num(1,2)=2
      ! Logical variable to check if input data is good
      logical :: implem_check = .true.
      logical :: check
      real(rp), dimension(:,:,:), allocatable :: array


      include 'include_codes/namelists/hamiltonian.f90'

      hoh = this%hoh
      local_axis = this%local_axis
      orb_pol = this%orb_pol
      v_alpha(:) = this%v_alpha(:)
      v_beta(:) = this%v_beta(:)
      js_alpha = this%js_alpha
      jl_alpha = this%jl_alpha
      call move_alloc(this%velocity_scale, velocity_scale)
      q_ss = this%q_ss

      call move_alloc(this%hubbard_u_sc, hubbard_u_sc)
      call move_alloc(this%hubbard_u, hubbard_u)
      call move_alloc(this%hubbard_j, hubbard_j)
      call move_alloc(this%hubbard_v, hubbard_v) 
      call move_alloc(this%uj_orb, uj_orb)

      call move_alloc(this%hubbard_u_impurity, hubbard_u_impurity)
      call move_alloc(this%hubbard_j_impurity, hubbard_j_impurity)

      call move_alloc(this%hubbard_u_general, hubbard_u_general)
      call move_alloc(this%hubbard_j_general, hubbard_j_general)

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
      this%q_ss = q_ss
      this%orb_pol = orb_pol
      this%v_alpha(:) = v_alpha(:)
      this%v_beta(:) = v_beta(:)
      this%js_alpha = js_alpha
      this%jl_alpha = jl_alpha
      call move_alloc(velocity_scale, this%velocity_scale)

      call move_alloc(hubbard_u_sc, this%hubbard_u_sc)
      call move_alloc(hubbard_u, this%hubbard_u)
      call move_alloc(hubbard_j, this%hubbard_j)
      call move_alloc(hubbard_v, this%hubbard_v)
      call move_alloc(uj_orb, this%uj_orb)
      call move_alloc(hubbard_u_impurity, this%hubbard_u_impurity)
      call move_alloc(hubbard_j_impurity, this%hubbard_j_impurity)

      call move_alloc(hubbard_u_general, this%hubbard_u_general)
      call move_alloc(hubbard_j_general, this%hubbard_j_general)

      ! New check to see if Hubbard U correction should be initiated
      ! This general check is the only check that is needed for the +U in bulk and impurity. 
      ! It is not yet compatible with the +V implementation. Thus the other checks have not 
      ! been removed yet.
      do na = 1, this%lattice%ntype
         do l = 1, 3
            if (abs(this%hubbard_u_general(na, l)) > 1e-8) then
               this%hubbard_u_general_check = .true.
            end if
         end do
      end do
      if ( this%hubbard_u_general_check ) then
         print *, '----------------------------------------------------------------------------------------'
         print *, 'Hubbard U correction initiated'
         print *, '----------------------------------------------------------------------------------------'
      end if
      
      
      ! Checks if Hubbard U correction for impurity should be initiated (check for j also??)
      if ( size(this%hubbard_u_impurity, 1) > 0 ) then
         do i = 1, size(this%hubbard_u_impurity, 1)
            do j = 1, size(this%hubbard_u_impurity, 2)
               if (abs(this%hubbard_u_impurity(i,j)) > 1e-8) then
                  this%hubbardU_impurity_check = .true.
               else if (abs(this%hubbard_j_impurity(i,j)) > 1e-8) then
                  this%hubbardU_impurity_check = .true.
               end if
            end do
         end do
      end if

      if (this%hubbardU_impurity_check) then
         ! Print inputs
         print *, "Size of hubbard_u_impurity: ", size(hubbard_u_impurity)
         print *, "lattice%nclu = ", this%lattice%nclu
         print *, "hubbard_u_impurity:"
         do i = 1, this%lattice%nrec
            write(*, *) 'Impurity', i, '=', this%hubbard_u_impurity(i, :)
         end do
         print *, "hubbard_j_impurity:"
         do i = 1, this%lattice%nrec
            write(*, *) 'Impurity', i, '=', this%hubbard_j_impurity(i, :)
         end do

         print *, "Hubbard U correction for impurity initiated."
         ! Calculates the Slater integrals F_impurity for impurities
         ! F_impurity(impurity, l+1, k) with k = 1,2,3,4 = l+1 indexing F^(2(k-1)) for impurities
         do i = 1, this%lattice%nrec ! For each impurity (same as size(hubbard_u_impurity,1) )
            ! For s, p, d and f orbitals, F0 = U
            ! s-orbital
            this%F_impurity(i,1,1) = this%hubbard_u_impurity(i,1)
            ! p-orbital J = (1/5)*F2 
            this%F_impurity(i,2,1) = this%hubbard_u_impurity(i,2)
            this%F_impurity(i,2,2) = this%hubbard_j_impurity(i,2)*5.0_rp
            ! d-orbital, J = (F2 + F4)/14 with F4/F2 ~ 0.625
            this%F_impurity(i,3,1) = this%hubbard_u_impurity(i,3)
            this%F_impurity(i,3,2) = 14.0_rp*this%hubbard_j_impurity(i,3)/1.625_rp 
            this%F_impurity(i,3,3) = 0.625_rp*this%F_impurity(i,3,2)
            !> For f orbitals, J = (286F2 + 195F4 + 250F6)/6435 with F4/F2 ~ 0.67 and F6/F2 ~ 0.49
            this%F_impurity(i,4,1) = this%hubbard_u_impurity(i,4)
            this%F_impurity(i,4,2) = 6435_rp*this%hubbard_j_impurity(i,4)/(286_rp + 195_rp*0.67_rp + 250_rp*0.49_rp)
            this%F_impurity(i,4,3) = 0.67_rp*this%F(i,4,2)
            this%F_impurity(i,4,4) = 0.49_rp*this%F(i,4,2)               
         end do 
         ! Print the slater intergrals
         print *,''
         print *, '----------------------------------------------------------------------------------------'
         print *, 'Slater Integrals.'
         print *, '----------------------------------------------------------------------------------------'
         do i = 1, this%lattice%nrec
            print *, 'Impurity ', i, ':'
            do j = 1, 4 ! Orbitals spdf have implementation
               if (count(this%F_impurity(i,j,:) > 1.0E-10) /= 0) then
                  print *, ' Orbital ', this%orb_conv(j), ' :'
                  if (j == 1) then
                     print *, '             F0 = ', this%F_impurity(i,j,1), ' eV'
                  else if (j == 2) then
                     print *, '             F0 = ', this%F_impurity(i,j,1), ' eV'
                     print *, '             F2 = ', this%F_impurity(i,j,2), ' eV'
                  else if (j == 3) then
                     print *, '             F0 = ', this%F_impurity(i,j,1), ' eV'
                     print *, '             F2 = ', this%F_impurity(i,j,2), ' eV'
                     print *, '             F4 = ', this%F_impurity(i,j,3), ' eV'
                  else if (j == 4) then
                     print *, '             F0 = ', this%F_impurity(i,j,1), ' eV'
                     print *, '             F2 = ', this%F_impurity(i,j,2), ' eV'
                     print *, '             F4 = ', this%F_impurity(i,j,3), ' eV'
                     print *, '             F6 = ', this%F_impurity(i,j,4), ' eV'
                  end if
               end if
            end do
         end do
         print *, '----------------------------------------------------------------------------------------'  
      end if
            

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
                  call g_logger%error('Self-consistent U calculation not implemented for f-electrons', __FILE__, __LINE__)
                  error stop
               end if
            else if ( this%hubbard_u_sc(na,l) .gt. 1 ) then
               print *, ''
               call g_logger%error('WRONG INPUT. hubbard_u_sc input was larger than 1, but only 0 and 1 is allowed.', __FILE__, __LINE__)
               error stop
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

      ! Establishes if Hubbard U should be implemented by checking if any Hubbard U is specified
      outer2 : do i = 1, this%lattice%nrec 
         do j = 1, len_trim(this%uj_orb(i)) 
            if (this%hubbard_u(i,j) > 1.0E-10) then ! If a nonzero Hubbard U parameter is specified, Hubbard U is initiated
               this%hubbardU_check = .true.
               print *, 'Hubbard U correction initiated.', __FILE__, __LINE__
               exit outer2
            end if
         end do
      end do outer2

      ! Establishes if Hubbard J should be implemented by checking if any Hubbard J is specified
      outer3 : do i = 1, this%lattice%nrec 
         do j = 1, len_trim(this%uj_orb(i)) 
            if (this%hubbard_j(i,j) > 1.0E-10) then ! If a nonzero Hubbard J parameter is specified, Hubbard J is initiated
               this%hubbardJ_check = .true.
               print *, 'Hubbard J correction initiated.', __FILE__, __LINE__
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
                     print *, 'Hubbard V correction initiated.', __FILE__, __LINE__
                  end if
               end do
            end do
         end do
      end do outer4

      ! Testing reading from atom.nml instead!
      ! ! Read Hubbard potential matrix for a +U calculation
      ! if ( this%hubbardU_check ) then
      !    open(unit=10, file='hubbard_potential', form='unformatted', status='old', iostat=ios)
      !    if (ios /= 0) then
      !       print *, 'No Hubbard potential was found in input.'
      !    else
      !       read(10) na
      !       if ( na == this%lattice%nrec ) then
      !          print *, 'Reads input Hubbard potential.'
      !          read(10) this%hubbard_u_pot
      !       end if
      !       close(10)
      !    end if
      ! end if

      ! ! Read Hubbard potential matrix for a +U impurity calculation
      ! if ( this%hubbardU_impurity_check ) then
      !    open(unit=10, file='hubbard_potential', form='unformatted', status='old', iostat=ios)
      !    if (ios /= 0) then
      !       print *, 'No hubbard potential was found for the bulk. Continues without.'
      !       if (allocated(this%hubbard_u_pot)) deallocate (this%hubbard_u_pot)
      !       allocate (this%hubbard_u_pot(18, 18, this%lattice%ntype))
      !       this%hubbard_u_pot = 0.0d0
      !    else
      !       read(10) na
      !       if ( na /= this%lattice%nbulk ) then
      !          print *, 'ERROR: Input Hubbard potential has different dimension than bulk,'
      !          print *, 'Stops program'
      !          stop
      !       else
      !          print *, 'Reads input Hubbard potential for impurity calculation.'
      !          if (allocated(array)) deallocate (array)
      !          allocate (array(18, 18, na))
      !          read(10) array

      !          if (allocated(this%hubbard_u_pot)) deallocate (this%hubbard_u_pot)
      !          allocate (this%hubbard_u_pot(18, 18, this%lattice%ntype))

      !          this%hubbard_u_pot = 0.0d0
      !          do k = 1, this%lattice%nbulk
      !             this%hubbard_u_pot(:, :, k) = array(:,:,k)
      !          end do
      !          deallocate (array)
      !       end if
      !       close(10)
      !    end if
      ! end if
      
      !> Raises an error if V correction is wanted when +U correction is not. Why not!?!
      if (this%hubbardV_check .and. .not. this%hubbardU_check) then
         print *, ''
         print *, '----------------------------------------------------------------------------------------'
         print *, '+V correction specified without +U correction!'
         print *, '----------------------------------------------------------------------------------------'
         call g_logger%error('Cannot add +V without +U!', __FILE__, __LINE__)
         error stop
      end if

      ! Organizes the orbital Hubbard U & J values into the correct position
      print *,'Hubbard U & J parameters.', this%hubbardU_check, this%hubbardJ_check
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
                     if (abs(this%hubbard_v(i,j,li,lj)) > 1.0E-10 .and. abs(this%hubbard_v(j,i,lj,li)) > 1.0E-10 .and. (this%hubbard_v(i,j,li,lj) - this%hubbard_v(j,i,lj,li) > 1.0E-10) ) then
                        print *, 'Incorrect +V implementation for interaction between atom type ',i , 'and ',j , '.'
                        call g_logger%error('Value for V(i,j,li,lj) must be the same as V(j,i,lj,li). Either delete one of them, or set them equal.', __FILE__, __LINE__)
                        error stop
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
      !> Calculates the Slater integrals, F(atom, l+1, k) with k = 1,2,3,4 = l+1 indexing F^(2(k-1))
         do i = 1, this%lattice%nrec ! For each atom
            do j = 1, count(this%hubbard_u(i,:) > 1.0E-10) ! For each nonzero hubbard U value
               if (this%uj_orb(i)(j:j) == 's') then ! If corresponding orbital is s
               !> For s, p, d and f orbitals, F0 = U
                  this%F(i,1,1) = this%hubbard_u(i,j)
               else if (this%uj_orb(i)(j:j) == 'p') then ! If corresponding orbital is p
               !> For p orbitals, J = (1/5)*F2 
                  this%F(i,2,1) = this%hubbard_u(i,j)
                  this%F(i,2,2) = this%hubbard_j(i,j)*5.0_rp
               else if (this%uj_orb(i)(j:j) == 'd') then ! If corresponding orbital is d
               !> For d orbitals, J = (F2 + F4)/14 with F4/F2 ~ 0.625
                  this%F(i,3,1) = this%hubbard_u(i,j)
                  this%F(i,3,2) = 14.0_rp*this%hubbard_j(i,j)/1.625_rp 
                  this%F(i,3,3) = 0.625_rp*this%F(i,3,2)
               else if (this%uj_orb(i)(j:j) == 'f') then ! If corresponding orbital is f
               !> For f orbitals, J = (286F2 + 195F4 + 250F6)/6435 with F4/F2 ~ 0.67 and F6/F2 ~ 0.49
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
         !> If no input values are given for some atom(s), default values is set given that Hubbard U and or J should be included
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
         !> Also if only U or J value is unspecified when orbital and J or U is, default value is set.
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
            !> If the number of orbitals and Hubbard U parameters per atom in input file disagree, then a raised error is prompted
            if (count(this%hubbard_u(i,:) > 1.0E-10) /= len_trim(this%uj_orb(i))) then
               implem_check = .false.
               print *, ''
               print *, '----------------------------------------------------------------------------------------'
               print *, 'Number of orbitals and Hubbard U parameters for atom', i, ' disagrees.'
               print *, '----------------------------------------------------------------------------------------'
            !> If the number of orbitals and Hubbard J parameters per atom in input file disagree, then a raised error is prompted 
            else if (count(this%hubbard_j(i,:) > 1.0E-10) /= len_trim(this%uj_orb(i))) then 
               implem_check = .false.
               print *, ''
               print *, '----------------------------------------------------------------------------------------'
               print *, 'Number of orbitals and Hubbard J parameters for atom', i, ' disagrees.'
               print *, '----------------------------------------------------------------------------------------'
            end if
         end if
         if (.not. this%hubbardU_check .and. .not. this%hubbardJ_check .and. len_trim(this%uj_orb(i)) > 0 ) then
            call g_logger%error('Hubbard orbitals specified without U parameters.', __FILE__, __LINE__)
            error stop
         end if
      end do outer

      !> If the +V (intersite Coulomb reaction) is added, we use the simplified +U correction which is spherically averaged.
      !> This corresponds to setting F^2 = F^4 = J = 0.
      if (this%hubbardV_check) then
         this%hubbard_j = 0.0d0
         this%F(:,:,2) = 0.0d0
         this%F(:,:,3) = 0.0d0
         this%F(:,:,4) = 0.0d0
         call this%lattice%neigh(1)
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
                  if (this%uj_orb(i)(j:j) /= '') then
                     print *, '  Orbital ', this%uj_orb(i)(j:j), ': Hubbard U = ', this%hubbard_u(i,j), ' eV.', ' Hubbard J = ', this%hubbard_j(i,j), 'eV'
                  end if
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
                     if (this%uj_orb(i)(j:j) /= '') then
                        print *, '  Orbital ', this%uj_orb(i)(j:j), ': Hubbard U = ', this%hubbard_u(i,j), ' eV.'
                     end if
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
               print *, ''
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
                     if (this%uj_orb(i)(j:j) /= '') then
                        print *, '  Orbital ', this%uj_orb(i)(j:j), ' Hubbard J = ', this%hubbard_j(i,j), 'eV'
                     end if
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
         ! print *, '----------------------------------------------------------------------------------------'
         ! print *, 'No Hubbard U+J data was given as input, proceeding without.'
         ! print *, '----------------------------------------------------------------------------------------'
      end if      
      
      !> Converts the quantities from eV to Ry
      this%hubbard_u = this%hubbard_u/ry2ev
      this%hubbard_j = this%hubbard_j/ry2ev
      this%hubbard_v = this%hubbard_v/ry2ev
      this%F = this%F/ry2ev
      this%hub_u_sort = this%hub_u_sort/ry2ev
      this%hub_j_sort = this%hub_j_sort/ry2ev

      this%hubbard_u_impurity = this%hubbard_u_impurity/ry2ev
      this%hubbard_j_impurity = this%hubbard_j_impurity/ry2ev
      this%F_impurity = this%F_impurity/ry2ev

      this%hubbard_u_general = this%hubbard_u_general/ry2ev
      this%hubbard_j_general = this%hubbard_j_general/ry2ev

      !> Checks if self-consistent U flag and input values of U and J have both been provided. They would interfer with eachother. 
      if ( (this%hubbardU_sc_check) .and. (this%hubbardU_check) .and. (this%hubbardJ_check)) then
         call g_logger%error('Both input values for hubbard_u_sc and hubbard_u + hubbard_j has been provided. Only one of them are allowed.', __FILE__, __LINE__)
         error stop
      else if ( (this%hubbardU_sc_check) .and. (this%hubbardU_check) ) then
         call g_logger%error('Both input values for hubbard_u_sc and hubbard_u has been provided. Only one of them are allowed.', __FILE__, __LINE__)
         error stop
      else if ( (this%hubbardU_sc_check) .and. (this%hubbardJ_check)) then 
         call g_logger%error('Both input values for hubbard_u_sc and hubbard_j has been provided. Only one of them is allowed.', __FILE__, __LINE__)
         error stop
      end if

      ! Check to see if hubbard U is in atom.nml input
      check = .False.
      do na = 1, this%lattice%ntype
         do l = 1, 3
            if ( abs(this%lattice%symbolic_atoms(na)%potential%hubbard_u(l)) > 1e-8 ) then
               check = .True.
            end if
         end do
      end do

      if ( check ) then
         print *, 'Hubbard was provided in the atom.nml files. Calculates the initial Hubbard potential matrix.'
         do na = 1, this%lattice%ntype
            call this%lattice%symbolic_atoms(na)%potential%expand_ldm()
         end do
         call this%calculate_hubbard_u_potential_general()
         if (this%hubbard_u_general_check) then
            print *, 'Hubbard parameters was provided in input.nml. Use these for the rest of the calculation.'
            do na = 1, this%lattice%ntype
               do l = 1, 3
                  this%lattice%symbolic_atoms(na)%potential%hubbard_u(l) = this%hubbard_u_general(na, l)
                  this%lattice%symbolic_atoms(na)%potential%hubbard_j(l) = this%hubbard_j_general(na, l)
               end do
            end do
         else
            print *, 'No Hubbard parameters was provided in the input.nml.'
            print *, 'Uses no Hubbard correction for the rest of the calculation.'
            do na = 1, this%lattice%ntype
               do l = 1, 3
                  this%lattice%symbolic_atoms(na)%potential%hubbard_u(l) = this%hubbard_u_general(na, l)
                  this%lattice%symbolic_atoms(na)%potential%hubbard_j(l) = this%hubbard_j_general(na, l)
               end do
            end do
         end if
      else
         if ( this%hubbard_u_general_check ) then
            print *, 'No Hubbard corrections in atom.nml, but Hubbard U provided in input.nml'
            print *, 'Stores these values in symbolic atoms.'
            do na = 1, this%lattice%ntype
               do l = 1, 3
                  this%lattice%symbolic_atoms(na)%potential%hubbard_u(l) = this%hubbard_u_general(na, l)
                  this%lattice%symbolic_atoms(na)%potential%hubbard_j(l) = this%hubbard_j_general(na, l)
               end do
            end do
         end if
      end if

      ! Check if bulk U was provided in input, but not in atom.nml for the impurity calculation
      if (this%lattice%nrec .ne. this%lattice%ntype) then
         do na = 1, this%lattice%nbulk
            do l = 1, 3
               if ( abs(this%lattice%symbolic_atoms(na)%potential%hubbard_u(l) - this%hubbard_u_general(na,l)) > 1e-8 ) then
                  print *, 'ERROR! Hubbard U in atom.nml is not equal to the one in input.nml for atom ', na, ' l = ', l - 1
                  print *, 'Hubbard U in atom.nml = ', this%lattice%symbolic_atoms(na)%potential%hubbard_u(l), 'Ry'
                  print *, 'Hubbard U in input.nml = ', this%hubbard_u_general(na,l), 'Ry'
                  print *, 'Hubbard parameters for the bulk atoms has to be equal in an impurity calculation.'
                  print *, 'STOP PROGRAM'
                  stop
               else if ( abs(this%lattice%symbolic_atoms(na)%potential%hubbard_j(l) - this%hubbard_j_general(na,l)) > 1e-8) then
                  print *, 'ERROR! Hubbard J in atom.nml is not equal to the one in input.nml for atom ', na, ' l = ', l - 1
                  print *, 'Hubbard J in atom.nml = ', this%lattice%symbolic_atoms(na)%potential%hubbard_u(l), 'Ry'
                  print *, 'Hubbard J in input.nml = ', this%hubbard_u_general(na,l), 'Ry'
                  print *, 'Hubbard parameters for the bulk atoms has to be equal in an impurity calculation.'
                  print *, 'STOP PROGRAM'
                  stop
               end if
            end do
         end do
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
      call g_safe_alloc%allocate('hamiltonian.v_a', this%v_a, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.v_b', this%v_b, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.vo_a', this%vo_a, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.vo_b', this%vo_b, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.js_a', this%js_a, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.jl_a', this%jl_a, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.jso_a', this%jso_a, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.jlo_a', this%jlo_a, (/18, 18, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.velocity_scale', this%velocity_scale, (/this%charge%lattice%ntype/))
      !end if
      !end if
#else
      allocate (this%lsham(18, 18, this%charge%lattice%ntype))
      allocate (this%tmat(18, 18, 3, this%charge%lattice%ntype))
      !allocate (this%hhmag(9, 9, 4), this%hmag(9, 9, this%charge%lattice%kk, 4))
      allocate (this%ee(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
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
      ! Velocity operators
      allocate (this%v_a(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%v_b(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%vo_a(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%vo_b(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%js_a(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%jl_a(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%jso_a(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%jlo_a(18, 18, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%velocity_scale(this%charge%lattice%ntype))
      allocate (this%hubbard_u(this%lattice%nrec, 1))
      allocate(this%hubbard_j(this%lattice%nrec, 1))
      allocate(this%hubbard_v(this%lattice%nrec, this%lattice%nrec, 4, 4))
      allocate (this%hub_u_sort(this%lattice%nrec, 4))
      allocate(this%hub_j_sort(this%lattice%nrec, 4))
      allocate(this%orbs_v(this%lattice%nrec, this%lattice%nrec))
      allocate(this%orbs_v_num(this%lattice%nrec, this%lattice%nrec))
      allocate(this%F(this%lattice%nrec, 4, 4))
      allocate (this%uj_orb(this%lattice%nrec))
      allocate (this%hubbard_u_pot(18, 18, this%lattice%ntype))
      allocate (this%hubbard_v_pot(18, 18, size(this%ee, 3) , this%lattice%nrec)) ! (lm, l'm', number of NN, number of atom types)
      allocate (this%hubbard_u_sc(this%lattice%nrec,4))
      allocate (this%hubbard_u_impurity(this%lattice%nrec,4)) ! Size is equal to the number of impurities (nclu)
      allocate (this%hubbard_j_impurity(this%lattice%nrec,4))
      allocate (this%F_impurity(this%lattice%nrec, 4, 4))
      allocate (this%hubbard_pot_impurity(18, 18, this%lattice%nrec))
      allocate (this%hubbard_u_general(this%lattice%ntype, 4)) ! These should work for both impurity and bulk!
      allocate (this%hubbard_j_general(this%lattice%ntype, 4))
      !end if
      !end if
#endif

      this%lsham(:, :, :) = 0.0d0
      this%tmat(:, :, :, :) = 0.0d0
      ! this%hhmag(:, :, :) = 0.0d0
      ! this%hmag(:, :, :, :) = 0.0d0
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
      this%v_a(:, :, :, :) = 0.0d0
      this%v_b(:, :, :, :) = 0.0d0
      this%vo_a(:, :, :, :) = 0.0d0
      this%vo_b(:, :, :, :) = 0.0d0
      this%js_a(:, :, :, :) = 0.0d0
      this%jl_a(:, :, :, :) = 0.0d0
      this%jso_a(:, :, :, :) = 0.0d0
      this%jlo_a(:, :, :, :) = 0.0d0
      this%velocity_scale(:) = 1.0d0
      this%hoh = .false.
      this%q_ss = [0.0d0, 0.0d0, 0.0d0]
      this%local_axis = .false.
      this%orb_pol = .false.
      this%v_alpha(:) = [1, 0, 0]
      this%v_beta(:) = [1, 0, 0] 
      this%js_alpha = 'z'
      this%jl_alpha = 'z'
      this%hubbard_u(:,:) = 0.0d0
      this%hubbard_j(:,:) = 0.0d0
      this%hubbard_v(:,:,:,:) = 0.0d0
      this%hub_u_sort(:,:) = 0.0d0
      this%hub_j_sort(:,:) = 0.0d0
      this%uj_orb(:) = ''
      this%orbs_v_num(:,:) = 0
      this%F(:,:,:) = 0.0d0
      this%hubbardU_check = .false.
      this%hubbardJ_check = .false.
      this%hubbardV_check = .false.
      this%hubbard_u_pot(:,:,:) = 0.0d0
      this%hubbard_v_pot(:,:,:,:) = 0.0d0
      this%orb_conv(1) = 's'
      this%orb_conv(2) = 'p'
      this%orb_conv(3) = 'd'
      this%orb_conv(4) = 'f'
      this%hubbard_u_sc(:,:) = 0
      this%hubbard_u_impurity(:,:) = 0.0d0
      this%hubbard_j_impurity(:,:) = 0.0d0
      this%F_impurity(:,:,:) = 0.0d0
      this%hubbard_pot_impurity(:,:,:) = 0.0d0
      this%hubbard_u_general(:,:) = 0.0d0
      this%hubbard_j_general(:,:) = 0.0d0
   end subroutine restore_to_default

   subroutine block_to_sparse(this)
       !*************************************************************************
       !> @brief Constructs the sparse Hamiltonian matrix for a real-space cluster.
       !>
       !> This subroutine loops through all cluster atoms and their neighbors to
       !> assemble the sparse Hamiltonian matrix in a block-wise manner. The matrix
       !> size is determined by the number of atoms in the cluster and their spd basis.
       !>
       !> @param[in] this       Hamiltonian type derived
       !> @param[out] H_sparse  Sparse matrix structure to store the Hamiltonian.
       !*************************************************************************
       class(hamiltonian), intent(inout) :: this
   
       ! Local variables
       integer :: kk, m, nr, i, j, i_start, j_start
       integer :: neighbor  ! Neighbor atom index
       real(rp), dimension(3) :: rij  ! Displacement vector (unused for now but available)
  
       if (allocated(this%h_sparse)) deallocate(this%h_sparse) 
       allocate(this%h_sparse(this%lattice%kk*18,this%lattice%kk*18))

       ! Loop over cluster atoms
       do kk = 1, this%lattice%kk  ! Loop over all atoms in the cluster
           ! Number of neighbors for the kk-th atom
           nr = this%charge%lattice%nn(kk, 1)
           ! Loop over neighbors (including onsite, m = 1)
           do m = 1, nr
               ! Compute block indices in the global matrix
               i_start = 18 * (kk - 1) + 1
               if (m == 1) then
                   j_start = i_start  ! Onsite term: diagonal block
               else
                   neighbor = this%charge%lattice%nn(kk, m)  ! Neighbor atom index
                   j_start = 18 * (neighbor - 1) + 1
               end if
               ! Add the block to the sparse matrix
               if (neighbor .ne. 0) then
                  this%h_sparse(i_start:i_start+17, j_start:j_start+17) = 0.0d0
                  do i = 1, 18
                      do j = 1, 18
                         this%h_sparse(i_start + i - 1, j_start + j - 1) = this%h_sparse(i_start + i - 1, j_start + j - 1) &
                                                                           + this%ee(i, j, m, 1)
                      end do
                  end do
               end if
           end do
       end do
       ! Placeholder: Incorporate atom type and neighbor type logic in future if required.
   end subroutine block_to_sparse


   !**************************************************************************
   !> @brief Build orbital-velocity operators (\f$\mathbf{j}^{L}\f$) by combining 
   !>        an orbital operator (\f$L_\gamma\f$) with an existing velocity operator.
   !>
   !> This subroutine constructs orbital-velocity operators that couple an orbital 
   !> operator (\f$L_\gamma\f$) and a velocity operator (\f$v_\alpha\f$) in the 
   !> same orbital (and possibly spin–orbital) basis. Similar to spin–velocity 
   !> constructs, one can use the anticommutator for symmetrization:
   !>
   !> \f[
   !>   j_\alpha^{L}(\gamma) 
   !>     = \tfrac12 
   !>       \Bigl(
   !>         L_\gamma \,v_\alpha \;+\; v_\alpha \,L_\gamma
   !>       \Bigr),
   !> \f]
   !>
   !> storing the resulting blocks in \f$j_{\alpha}^{L}(\gamma)\f$. This approach 
   !> may be used to evaluate orbital currents or to track orbital angular momentum 
   !> flow in a lattice system.
   !>
   !> @param[inout] this
   !>   A derived type representing the Hamiltonian system. It contains both:
   !>   - The orbital velocity operator \f$\hat{v}\f$ (e.g., \f$this\%\text{v\_a}\f$),
   !>   - The orbital operator(s) \f$L_x, L_y, L_z\f$,
   !>   - The data structure to hold the new orbital–velocity operator 
   !>     \f$j^L_\alpha(\gamma)\f$.
   !> 
   !> @details
   !> - The subroutine loops over the same block structure (\f$\mathrm{dim}\times\mathrm{dim}\f$) 
   !>   as in real-space velocity. Each block is multiplied twice: 
   !>   \f$L_\gamma\,v_\alpha\f$ and \f$v_\alpha\,L_\gamma\f$, then averaged 
   !>   (\f$\times 0.5\f$).
   !> - One can repeat this logic for each \f$L_\gamma\f$ (\f$L_x,L_y,L_z\f$) and each velocity 
   !>   operator (\f$v_a,v_b\f$) to construct the entire set of orbital–velocity arrays 
   !>   for subsequent Chebyshev or Kubo expansions.
   !>
   !> @warning
   !>   - Ensure that \f$L_\gamma\f$ is defined in the **same** dimension/basis as 
   !>     \f$v_\alpha\f$. A mismatch in matrix size or basis alignment will 
   !>     invalidate the block multiplications.
   !>   - If additional factors (e.g., \(\hbar\) or a sign convention) are required 
   !>     for \f$L_\gamma\f$, confirm it is consistent with how \f$v_\alpha\f$ 
   !>     was defined.
   !>
   !**************************************************************************
   subroutine build_realspace_orbital_velocity_operators(this)
      class(hamiltonian), intent(inout) :: this
   
      integer :: ntype, ia, nr, m
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:,:), tmp2(:,:), L_op(:,:)  ! Temp matrices
      complex(rp), dimension(9, 9) :: mLx, mLy, mLz

      hblocksize = size(this%v_a, 1)
   
      ! Allocate local matrices for partial products
      allocate(tmp1(hblocksize,hblocksize), tmp2(hblocksize,hblocksize), L_op(hblocksize,hblocksize))
   
      ! Initialize the orbital–velocity array
      this%jl_a(:, :, :, :) = (0.0_rp, 0.0_rp)
      this%jlo_a(:, :, :, :) = (0.0_rp, 0.0_rp) 

      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
      mLx(:, :) = L_x(:, :)
      mLy(:, :) = L_y(:, :)
      mLz(:, :) = L_z(:, :)
      
      ! Transforming them into the spherical harmonics coordinates
      call hcpx(mLx, 'cart2sph')
      call hcpx(mLy, 'cart2sph')
      call hcpx(mLz, 'cart2sph')
   
      ! Pick which orbital operator L_x, L_y, or L_z based on some user choice
      select case (this%jl_alpha)   ! or whichever variable holds 'x','y','z'
      case ('x')
         L_op(1:9, 1:9) = mLx(:, :)
         L_op(10:18, 10:18) = mLx(:, :)
      case ('y')
         L_op(1:9, 1:9) = mLy(:, :)
         L_op(10:18, 10:18) = mLy(:, :)
      case ('z')
         L_op(1:9, 1:9) = mLz(:, :)
         L_op(10:18, 10:18) = mLz(:, :)
      end select
   
      ! Loop over each atom type
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         nr = this%charge%lattice%nn(ia, 1)
   
         ! For each neighbor block
         do m = 2, nr
            ! tmp1 = L_op * v_a(:,:,m,ntype)
            tmp1 = matmul(L_op, this%v_a(:,:,m,ntype))
   
            ! tmp2 = v_a(:,:,m,ntype) * L_op
            tmp2 = matmul(this%v_a(:,:,m,ntype), L_op)
   
            ! jl_a(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
            this%jl_a(:,:,m,ntype) = 0.5_rp * ( tmp1 + tmp2 )
   
            if (this%hoh) then
               tmp1 = 0.0d0; tmp2 = 0.0d0
               ! tmp1 = js_a * v_a(:,:,m,ntype)
               tmp1 = matmul(L_op, this%vo_a(:, :, m, ntype))

               ! tmp2 = v_a(:,:,m,ntype) * js_a
               tmp2 = matmul(this%vo_a(:, :, m, ntype), L_op)

               ! v_sza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
               this%jlo_a(:,:,m,ntype) = 0.5_rp * ( tmp1 + tmp2 )
            end if

            ! Optional debugging output:
            ! write(*,*) 'm=', m
            ! write(*,'(18f10.6)') real(this%jo_a(:,:,m,ntype))
            ! write(*,*)
            ! write(*,'(18f10.6)') aimag(this%jo_a(:,:,m,ntype))
         end do
      end do
   
      deallocate(tmp1, tmp2, L_op)
   end subroutine build_realspace_orbital_velocity_operators


   !**************************************************************************
   !> @brief Build spin-velocity operators (\f$\mathbf{j}^{s}\f$) by combining 
   !>        the spin operator (\f$S_\beta\f$) with an existing velocity operator.
   !>
   !> This subroutine constructs spin-velocity operators that couple a spin operator 
   !> (\f$S_\beta\f$) and a velocity operator (\f$v_\alpha\f$) within the same spin–orbital 
   !> basis. In many spin Hall or related calculations, one uses the anticommutator:
   !>
   !> \f[
   !>   j_\alpha^{s}(\beta) 
   !>     = \tfrac12 
   !>       \Bigl(
   !>         S_\beta \,v_\alpha \;+\; v_\alpha \,S_\beta
   !>       \Bigr),
   !> \f]
   !>
   !> storing the resulting blocks in \f$j_{\alpha}^{s}(\beta)\f$. This approach 
   !> is commonly used to evaluate spin currents along a chosen direction \(\alpha\) 
   !> with spin polarization \(\beta\).
   !>
   !> @param[inout] this
   !>   A derived type representing the Hamiltonian system. It contains both:
   !>   - The orbital velocity operator \f$\hat{v}\f$ (e.g., \f$this\%\text{v\_a}\f$),
   !>   - The spin operator(s) \f$S_x, S_y, S_z\f$,
   !>   - The data structure to hold the new spin–velocity operator 
   !>     \f$j^s_\alpha(\beta)\f$.
   !> 
   !> @details
   !> - The subroutine loops over the same blocks (\f$\mathrm{dim}\times\mathrm{dim}\f$) 
   !>   as in real-space velocity. Each block is multiplied twice: once \f$S_\beta v_\alpha\f$ 
   !>   and once \f$v_\alpha S_\beta\f$, and the results are averaged (multiplied by 0.5).
   !> - One typically repeats this for each spin operator component (\f$S_x, S_y, S_z\f$) 
   !>   and/or each velocity operator (\f$v_a, v_b\f$) to build the needed spin–velocity 
   !>   arrays for Chebyshev or Kubo expansions.
   !>
   !> @warning
   !>   - The spin operators \f$S_\beta\f$ must be defined in the **same** spin–orbital 
   !>     dimension as the velocity operator \f$v_\alpha\f$. Inconsistencies in dimension 
   !>     or basis alignment will lead to incorrect matrix products.
   !>   - Ensure that the factor \(\frac{1}{i}\) or \(\hbar\equiv1\) is consistently 
   !>     applied to both velocity and spin definitions, if required.
   !>
   !**************************************************************************
   subroutine build_realspace_spin_operators(this)
      class(hamiltonian), intent(inout) :: this
   
      integer :: ntype, ia, nr, m, ji, ja, atom_neighbor
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:, :), tmp2(:, :), S_op(:, :)  ! Temp matrices for partial products
   
      ! Derive dimension from your velocity array:
      hblocksize = size(this%v_a, 1)  ! e.g. first dimension of v_a
   
      ! Allocate temporary matrices for local block multiplication
      allocate(tmp1(hblocksize, hblocksize), tmp2(hblocksize, hblocksize), S_op(hblocksize, hblocksize))
   
      ! Initialize the spin–velocity array to zero
      this%js_a(:, :, :, :) = (0.0_rp, 0.0_rp)
      this%jso_a(:, :, :, :) = (0.0_rp, 0.0_rp)
      
      select case(this%js_alpha)
      
      case('z')
         S_op = S_z
      case('x')
         S_op = S_x
      case('y')
         S_op = S_y
      end select

      !write(*,'(18f10.6)') real(S_op)
      ! Loop over each atom type
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         nr = this%charge%lattice%nn(ia, 1)
   
         ! For each neighbor block 
         do m = 2, nr

            tmp1 = 0.0d0; tmp2 = 0.0d0
            ! tmp1 = js_a * v_a(:,:,m,ntype)
            tmp1 = matmul(S_op, this%v_a(:, :, m, ntype))
   
            ! tmp2 = v_a(:,:,m,ntype) * js_a
            tmp2 = matmul(this%v_a(:, :, m, ntype), S_op)

            ! v_sza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
            this%js_a(:,:,m,ntype) = 0.5_rp * ( tmp1 + tmp2 )
            !write(*,*) 'm=', m
            !write(*,'(18f10.6)') real(this%js_a(:,:,m,ntype))
            !write(*,*)
            !write(*,'(18f10.6)') aimag(this%js_a(:,:,m,ntype)) 
            if (this%hoh) then
               tmp1 = 0.0d0; tmp2 = 0.0d0
               ! tmp1 = js_a * v_a(:,:,m,ntype)
               tmp1 = matmul(S_op, this%vo_a(:, :, m, ntype))
   
               ! tmp2 = v_a(:,:,m,ntype) * js_a
               tmp2 = matmul(this%vo_a(:, :, m, ntype), S_op)
   
               ! v_sza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
               this%jso_a(:,:,m,ntype) = 0.5_rp * ( tmp1 + tmp2 )
            end if

         end do  ! m
      end do  ! ntype
   
      deallocate(tmp1, tmp2)
   end subroutine build_realspace_spin_operators

   !**************************************************************************
   !> @brief Build **spin-torque operators** (\f$\boldsymbol{\tau}\f$) by taking the
   !>        commutator of the spin operator (\f$S_\beta\f$) with the Hamiltonian.
   !>
   !> This subroutine constructs layer- or orbital-resolved spin-torque operators
   !> that couple a spin operator (\f$S_\beta\f$) to the system Hamiltonian
   !> (\f$H\f$) within the same spin–orbital basis.  For torkance calculations one
   !> uses the commutator
   !>
   !> \f[
   !>   \boxed{\;
   !>     \tau_\beta
   !>        = \frac{1}{i\hbar}\bigl[\,S_\beta,\,H\,\bigr]
   !>        = \frac{1}{i\hbar}
   !>          \Bigl(
   !>            S_\beta\,H \;-\; H\,S_\beta
   !>          \Bigr)
   !>   \;}
   !> \f]
   !>
   !> The resulting blocks are stored in \f$\tau_\beta\f$.  These operators enter
   !> linear-response formulas for the (field-like / damping-like) torkance
   !> tensor.
   !>
   !> @param[inout] this
   !>   A derived type representing the Hamiltonian system.  It contains:
   !>   - The Hamiltonian blocks \f$\hat{H}\f$ (e.g.\ \f$this\%\text{H}\f$),
   !>   - The spin operator(s) \f$S_x, S_y, S_z\f$,
   !>   - The data structure to hold the new spin-torque operator \f$\tau_\beta\f$.
   !>
   !> @details
   !> - The subroutine loops over the same blocks (\f$\mathrm{dim}\!\times\!\mathrm{dim}\f$)
   !>   used for the Hamiltonian.  For each block it evaluates
   !>   \f$S_\beta H - H S_\beta\f$ and multiplies by \f$(i/\hbar)\f$.
   !> - One repeats this for every spin component (\f$S_x, S_y, S_z\f$) to build
   !>   the full set of spin-torque arrays required by Chebyshev or Kubo routines.
   !>
   !> @warning
   !>   - The spin operators \f$S_\beta\f$ and Hamiltonian blocks \f$H\f$ must share
   !>     the **same** spin–orbital dimension and ordering; any mismatch yields an
   !>     incorrect commutator.
   !>   - Check that the factor \f$1/\hbar\f$ is consistent with units used
   !>     elsewhere (set \f$\hbar\!\equiv\!1\f$ if that is your convention).
   !>
   !**************************************************************************
   subroutine build_realspace_spin_torque_operators(this)
      class(hamiltonian), intent(inout) :: this

      integer :: ntype, ia, nr, m, ji, ja, atom_neighbor, ino
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:, :), tmp2(:, :), S_op(:, :)  ! Temp matrices for partial products
      complex(rp), dimension(18, 18) :: locham

      ! Derive dimension from your velocity array:
      hblocksize = size(this%v_a, 1)  ! e.g. first dimension of v_a

      ! Allocate temporary matrices for local block multiplication
      allocate(tmp1(hblocksize, hblocksize), tmp2(hblocksize, hblocksize), S_op(hblocksize, hblocksize))

      ! Initialize the spin–velocity array to zero
      this%js_a(:, :, :, :) = (0.0_rp, 0.0_rp)
      this%jso_a(:, :, :, :) = (0.0_rp, 0.0_rp)

      select case(this%js_alpha)

      case('z')
         S_op = S_z
      case('x')
         S_op = S_x
      case('y')
         S_op = S_y
      end select

      !write(*,'(18f10.6)') real(S_op)
      ! Loop over each atom type
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         nr = this%charge%lattice%nn(ia, 1)
          
         ! For each neighbor block 
         do m = 1, nr

            if (m==1) then
              locham(:,:) = this%ee(:, :, m, ntype) + this%lsham(:, :, ntype) 
            else
              locham(:,:) = this%ee(:, :, m, ntype)
            end if

            tmp1 = 0.0d0; tmp2 = 0.0d0
            ! tmp1 = js_a * ee(:,:,m,ntype)
            tmp1 = matmul(S_op, locham(:, :))

            ! tmp2 = ee(:,:,m,ntype) * js_a
            tmp2 = matmul(locham(:, :), S_op)

            ! v_sza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
            this%js_a(:,:,m,ntype) = (1 / i_unit) * ( tmp1 - tmp2 )
            !write(*,*) 'm=', m
            !write(*,'(18f10.6)') real(this%js_a(:,:,m,ntype))
            !write(*,*)
            !write(*,'(18f10.6)') aimag(this%js_a(:,:,m,ntype)) 
            if (this%hoh) then

               if (m==1) then
                 locham(:,:) = this%eeo(:, :, 1, ntype) + this%lsham(:, :, ntype)
               else
                 locham(:,:) = this%eeo(:, :, m, ntype)
               end if

               tmp1 = 0.0d0; tmp2 = 0.0d0
               ! tmp1 = js_a * eeo(:,:,m,ntype)
               tmp1 = matmul(S_op, locham(:, :))

               ! tmp2 = ee(:,:,m,ntype) * js_a
               tmp2 = matmul(locham(:, :), S_op)

               ! v_sza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
               this%jso_a(:,:,m,ntype) = (1 / i_unit) * ( tmp1 - tmp2 )
            end if

         end do  ! m
      end do  ! ntype

      deallocate(tmp1, tmp2)
   end subroutine build_realspace_spin_torque_operators


   subroutine build_realspace_orbital_torque_operators(this)
      class(hamiltonian), intent(inout) :: this

      integer :: ntype, ia, nr, m, ji, ja, atom_neighbor, ino
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:, :), tmp2(:, :), L_op(:, :)  ! Temp matrices for partial products
      complex(rp), dimension(9, 9) :: mLx, mLy, mLz
      complex(rp), dimension(18, 18) :: locham

      ! Derive dimension from your velocity array:
      hblocksize = size(this%v_a, 1)  ! e.g. first dimension of v_a

      ! Allocate temporary matrices for local block multiplication
      allocate(tmp1(hblocksize, hblocksize), tmp2(hblocksize, hblocksize), L_op(hblocksize, hblocksize))

      ! Initialize the spin–velocity array to zero
      this%jl_a(:, :, :, :) = (0.0_rp, 0.0_rp)
      this%jlo_a(:, :, :, :) = (0.0_rp, 0.0_rp)

      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
      mLx(:, :) = L_x(:, :)
      mLy(:, :) = L_y(:, :)
      mLz(:, :) = L_z(:, :)

      ! Transforming them into the spherical harmonics coordinates
      call hcpx(mLx, 'cart2sph')
      call hcpx(mLy, 'cart2sph')
      call hcpx(mLz, 'cart2sph')

      ! Pick which orbital operator L_x, L_y, or L_z based on some user choice
      select case (this%jl_alpha)   ! or whichever variable holds 'x','y','z'
      case ('x')
         L_op(1:9, 1:9) = mLx(:, :)
         L_op(10:18, 10:18) = mLx(:, :)
      case ('y')
         L_op(1:9, 1:9) = mLy(:, :)
         L_op(10:18, 10:18) = mLy(:, :)
      case ('z')
         L_op(1:9, 1:9) = mLz(:, :)
         L_op(10:18, 10:18) = mLz(:, :)
      end select

      ! Loop over each atom type
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         nr = this%charge%lattice%nn(ia, 1)

         ! For each neighbor block 
         do m = 1, nr

            if (m==1) then
              locham(:,:) = this%ee(:, :, m, ntype) + this%lsham(:, :, ntype)
            else
              locham(:,:) = this%ee(:, :, m, ntype)
            end if

            tmp1 = 0.0d0; tmp2 = 0.0d0
            ! tmp1 = jl_a * ee(:,:,m,ntype)
            tmp1 = matmul(L_op, locham(:, :))

            ! tmp2 = ee(:,:,m,ntype) * jl_a
            tmp2 = matmul(locham(:, :), L_op)

            ! v_lza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
            this%jl_a(:,:,m,ntype) = (1 / i_unit) * ( tmp1 - tmp2 )
            if (this%hoh) then

               if (m==1) then
                 locham(:,:) = this%eeo(:, :, 1, ntype) + this%lsham(:, :, ntype)
               else
                 locham(:,:) = this%eeo(:, :, m, ntype)
               end if

               tmp1 = 0.0d0; tmp2 = 0.0d0
               ! tmp1 = jl_a * eeo(:,:,m,ntype)
               tmp1 = matmul(L_op, locham(:, :))

               ! tmp2 = ee(:,:,m,ntype) * jl_a
               tmp2 = matmul(locham(:, :), L_op)

               ! v_lza(:,:,m,ntype) = 0.5 * ( tmp1 + tmp2 )
               this%jlo_a(:,:,m,ntype) = (1 / i_unit) * ( tmp1 - tmp2 )
            end if

         end do  ! m
      end do  ! ntype

      deallocate(tmp1, tmp2)
   end subroutine build_realspace_orbital_torque_operators


   !**************************************************************************
   !> @brief Build real-space velocity operators.
   !>
   !> This subroutine constructs the velocity operators (\f$v_x\f$, \f$v_y\f$, \f$v_z\f$)
   !> in real space for a given Hamiltonian system, using the displacement vectors
   !> between atom pairs and the intersite Hamiltonian blocks. The velocity operators
   !> are computed based on the relationship:
   !> \f[
   !> v_{\alpha} = i \cdot (\mathbf{r}_i - \mathbf{r}_j)_{\alpha} \cdot H_{ij}
   !> \f]
   !> where \f$\mathbf{r}_i\f$ and \f$\mathbf{r}_j\f$ are the positions of atoms \f$i\f$ 
   !> and \f$j\f$, and \f$H_{ij}\f$ is the intersite Hamiltonian block.
   !>
   !> @param[inout] this        A derived type representing the Hamiltonian system.
   !>                           Contains Hamiltonian blocks, lattice structure, and
   !>                           other system parameters.
   !>
   !> @warning Ensure that the lattice positions (\f$\mathbf{r}\f$) and Hamiltonian
   !>          blocks (\f$H_{ij}\f$) are correctly initialized before calling this
   !>          subroutine.
   !>
   !**************************************************************************
   subroutine build_realspace_velocity_operators(this)
      ! Arguments
      class(hamiltonian), intent(inout) :: this
   
      ! Local variables
      integer :: ia, ntype, nr, m, i, j, velotype, ja, ji    ! Atom and neighbor indices
      integer :: atom_neighbor                               ! Neighbor atom index
      real(rp) :: veloscale
      real(rp), dimension(3) :: rij                          ! Displacement vector (x, y, z components)
      real(rp), dimension(3) :: dir_a, dir_b                 ! Velocity operator directions
      real(rp) :: norm_a, norm_b, dot_a, dot_b
      ! Initialize velocity operators to zero
      this%v_a(:, :, :, :) = 0.0_rp
      this%v_b(:, :, :, :) = 0.0_rp
   
      norm_a = norm2(this%v_alpha)
      norm_b = norm2(this%v_beta)

      dir_a(:) = this%v_alpha(:) / norm_a
      dir_b(:) = this%v_beta(:) / norm_b

      ! Loop over atom types
      do ntype = 1, this%charge%lattice%ntype

         ia = this%charge%lattice%atlist(ntype)  ! Atom number in the cluster
         nr = this%charge%lattice%nn(ia, 1)     ! Number of neighbors for this atom type
    
         ! Loop over neighbors
         do m = 2, nr   ! Start from 2 to exclude the onsite term

            atom_neighbor = this%charge%lattice%nn(ia, m)  ! Neighbor atom number
            if (atom_neighbor /= 0) then
               ! Compute displacement vector rij = r_i - r_j
               rij(:) = (this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, atom_neighbor)) * this%charge%lattice%alat
   
               dot_a = dot_product(dir_a, rij); dot_b = dot_product(dir_b, rij)

               ! Compute velocity operator blocks
               this%v_a(:, :, m, ntype) = ((1 / i_unit) * dot_a * this%ee(:, :, m, ntype)) 
            
               velotype = this%charge%lattice%iz(atom_neighbor)
               veloscale = max(this%velocity_scale(ntype), this%velocity_scale(velotype))
               write(*,*) ntype, velotype, veloscale
               this%v_b(:, :, m, ntype) = ((1 / i_unit) * dot_b * this%ee(:, :, m, ntype)) * veloscale 
               ! If hoh is true, multiply the velocity operator by the overlap matrix, similarly to whats done to the Hamiltonian
               if (this%hoh) then
                  ji = this%charge%lattice%iz(atom_neighbor) 
                  call zgemm('n', 'n', 18, 18, 18, cone, this%v_a(:, :, m, ntype), 18, this%obarm(:, :, ji), 18, czero, this%vo_a(:, :, m, ntype), 18)
                  call zgemm('n', 'n', 18, 18, 18, cone, this%v_b(:, :, m, ntype), 18, this%obarm(:, :, ji), 18, czero, this%vo_b(:, :, m, ntype), 18)
               end if
            end if
         end do
      end do
   end subroutine build_realspace_velocity_operators

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
      complex(rp), dimension(:,:,:,:), allocatable :: hmag

      allocate(hmag(9, 9, this%charge%lattice%nn_max, 4))
      hmag = (0.0d0, 0.0d0)

      if (this%hoh) then
         call this%build_obarm()
         call this%build_enim()
      end if

      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         !write(123, *)´bulkham´
         call this%chbar_nc(ia, nr, hmag)
         do m = 1, nr
            do i = 1, 9
               do j = 1, 9
                  this%ee(j, i, m, ntype) = hmag(j, i, m, 4) + hmag(j, i, m, 3)        ! H0+Hz
                  this%ee(j + 9, i + 9, m, ntype) = hmag(j, i, m, 4) - hmag(j, i, m, 3)        ! H0-Hz
                  this%ee(j, i + 9, m, ntype) = hmag(j, i, m, 1) - i_unit*hmag(j, i, m, 2) ! Hx-iHy
                  this%ee(j + 9, i, m, ntype) = hmag(j, i, m, 1) + i_unit*hmag(j, i, m, 2) ! Hx+iHy
               end do ! end of orbital j loop
            end do ! end of orbital i loop
            write(128, *) 'm=', m, 'ntype= ', ntype
            write(128, '(18f10.6)') real(this%ee(:, :, m, ntype))
            ! print '(9f6.2)', real(hmag(:, :, m, 4))
            !print '(18f6.2)', real(this%hall(:, :, m, nlim))
            ! print *,'========================'
         end do ! end of neighbour number
         ! Hubbard U correction.
         ! Only implemented for spd-orbitals
         ! if (this%hubbardU_check .or. this%hubbardU_impurity_check) then
         if (this%hubbard_u_general_check) then
            print *, 'Adding Hubbard U correction to the Hamiltonian', sum(abs(this%hubbard_u_pot))
            do i = 1, 9
               do j = 1, 9
                  this%ee(i, j, 1, ntype) = this%ee(i, j, 1, ntype) + this%hubbard_u_pot(i, j, ntype)
                  !this%ee(i + 9, j + 9, 1, ntype) = this%ee(i + 9, j + 9, 1, ntype) - this%hubbard_u_pot(i + 9, j + 9, ntype)
                  this%ee(i + 9, j + 9, 1, ntype) = this%ee(i + 9, j + 9, 1, ntype) + this%hubbard_u_pot(i + 9, j + 9, ntype)
               end do
            end do
            print '(i3, a, 9f8.4)', ino, ' Up', (this%hubbard_u_pot(i, i, ntype), i=1,9)
            print '(i3, a, 9f8.4)', ino, ' Dw', (this%hubbard_u_pot(i, i, ntype), i=10,18)
            ! print *, 'hub_v_pot test ', this%hubbard_v_pot
            if (this%hubbardV_check) then
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
            ! call this%build_obarm()
            ! call this%build_enim()
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
      deallocate(hmag)
      !!! AB 270125 Testing Gershgorin circles for later implementation
      !!! ! Quick check/hack for Gershgorin circles
      !!! g_max = -100.0_rp
      !!! g_min = 100.0_rp
      !!! g_mid = 0.0_rp
      !!! do ntype = 1, this%charge%lattice%ntype
      !!!    nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
      !!!    !write(123, *)´bulkham´
      !!!    do m = 1, nr
      !!!       do i = 1, 9
      !!!          g_mid = this%ee(i, i, m, ntype)
      !!!          g_rad = 0.0_rp
      !!!          do j = 1, 9
      !!!             if (i /= j) g_rad = g_rad + abs(this%ee(j, i, m, ntype))
      !!!          end do ! end of orbital j loop
      !!!          g_min = min(g_mid - g_rad, g_min)
      !!!          g_max = max(g_mid + g_rad, g_max)
      !!!       end do ! end of orbital i loop
      !!!    end do ! end of neighbour number
      !!! end do ! end of neighbour number
      !!! print *, 'Gershgorin circles: ', g_min, g_max
      !!! ! Additional factor for safety
      !!! this%g_max = g_max * sqrt(2.0_rp)
      !!! this%g_min = g_min * sqrt(2.0_rp)
   end subroutine build_bulkham

   subroutine build_locham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: it, ino, nr, nlim, m, i, j, ja, ji, iz
      complex(rp), dimension(:,:,:,:), allocatable :: hmag

      ! print *, 'Building local Hamiltonian', this%charge%lattice%nmax, ' atoms'
      call g_timer%start('Build local hamiltonian')
      allocate(hmag(9, 9, this%charge%lattice%nn_max, 4))

      !$omp parallel do private(nlim, nr, ino, m, i, j, ji, ja, hmag)
      do nlim = 1, this%charge%lattice%nmax
         ! print *, 'Building local Hamiltonian for atom ', nlim, ' of ', this%charge%lattice%nmax
         nr = this%charge%lattice%nn(nlim, 1) ! Number of neighbours considered
         ino = this%charge%lattice%num(nlim)
         iz = this%charge%lattice%iz(nlim) ! Atom type in the input.nml file
         call this%chbar_nc(nlim, nr, hmag)
         do m = 1, nr
            do i = 1, 9
               do j = 1, 9
                  this%hall(j, i, m, nlim) = hmag(j, i, m, 4) + hmag(j, i, m, 3) ! H0+Hz
                  this%hall(j + 9, i + 9, m, nlim) = hmag(j, i, m, 4) - hmag(j, i, m, 3) ! H0-Hz
                  this%hall(j, i + 9, m, nlim) = hmag(j, i, m, 1) - i_unit*hmag(j, i, m, 2) ! Hx-iHy
                  this%hall(j + 9, i, m, nlim) = hmag(j, i, m, 1) + i_unit*hmag(j, i, m, 2) ! Hx+iHy
               end do
            end do
            ! print '(9f6.2)', real(hmag(:, :, m, 4))
            !print '(18f6.2)', real(this%hall(:, :, m, nlim))
            ! print *,'---------------------'
         end do
         if (this%hoh) then
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
         ! if (this%hubbardU_impurity_check) then
         if (this%hubbard_u_general_check) then
            do i = 1, 9
               do j = 1, 9
                  this%hall(i, j, 1, nlim) = this%hall(i, j, 1, nlim) + this%hubbard_u_pot(i, j, iz)
                  this%hall(i + 9, j + 9, 1, nlim) = this%hall(i + 9, j + 9, 1, nlim) + this%hubbard_u_pot(i + 9, j + 9, iz)
               end do
            end do
         end if
      end do
      !$omp end parallel do
      if (this%local_axis) then
         this%hall_glob = this%hall
         if (this%hoh) this%hallo_glob = this%hallo
      end if
      deallocate(hmag)
      !print '(2g12.5)', this%hall
      ! this%hall = (0.0d0, 0.0d0) ! Free memory
      call g_timer%stop('Build local hamiltonian')
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

               if (k == 1) this%ee(:, :, k, ntype) = this%ee(:, :, k, ntype) + this%lsham(:, :, ntype) !+ this%enim(:,:,ntype)

               call hcpx(this%ee(1:9,1:9,k,ntype), 'sph2cart')
               call hcpx(this%ee(1:9,10:18,k,ntype), 'sph2cart')
               call hcpx(this%ee(10:18,1:9,k,ntype), 'sph2cart')
               call hcpx(this%ee(10:18,10:18,k,ntype), 'sph2cart')
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

   subroutine ham0m_nc(this, ia, ja, it, jt, vet, hhh, hhmag)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: ia, ja ! Atom sites i and j
      integer, intent(in) :: it, jt ! Type of atom i and j
      real(rp), dimension(3), intent(in) :: vet
      real(rp), dimension(9, 9), intent(in) :: hhh
      ! Local Variables
      integer :: i, j, ilm, jlm, m
      real(rp), dimension(3) :: mom_ia, mom_ja
      real(rp), dimension(3) :: r_ia, r_ja
      complex(rp), dimension(3) :: cross
      complex(rp), dimension(9, 9) :: hhhc
      ! complex(rp), dimension(this%charge%lattice%ntype, 3) :: momc
      complex(rp), dimension(3) :: momc_ia, momc_ja
      complex(rp) :: dot
      real(rp) :: vv
      complex(rp), dimension(9, 9, 4), intent(out) :: hhmag

      hhmag(:, :, :) = 0.0d0
      !this%hhmag(:, :, :) = 0.0d0

      vv = norm2(vet)
      mom_ia = this%charge%lattice%symbolic_atoms(it)%potential%mom(:)
      mom_ja = this%charge%lattice%symbolic_atoms(jt)%potential%mom(:)
      if (norm2(this%q_ss)>0.00001_rp) then
         !print *, 'q:', this%q_ss
         r_ia = this%charge%lattice%cr(:, ia)
         r_ja = this%charge%lattice%cr(:, ja)
         mom_ia(3) = 0.0d0
         mom_ia(2) = sin(2.0d0*pi*dot_product(r_ia, this%q_ss))
         mom_ia(1) = cos(2.0d0*pi*dot_product(r_ia, this%q_ss))
         mom_ja(3) = 0.0d0
         mom_ja(2) = sin(2.0d0*pi*dot_product(r_ja, this%q_ss))
         mom_ja(1) = cos(2.0d0*pi*dot_product(r_ja, this%q_ss))
      end if
      ! Real to complex
      dot = cmplx(dot_product(mom_ia, mom_ja), kind=kind(0.0d0))
      momc_ia = cmplx(mom_ia, kind=kind(0.0d0))
      momc_ja = cmplx(mom_ja, kind=kind(0.0d0))
      cross = cmplx(cross_product(mom_ia, mom_ja), kind=kind(0.0d0))
      hhhc(:, :) = cmplx(hhh(:, :), kind=kind(0.0d0))

      do ilm = 1, 9
         do jlm = 1, 9
            hhmag(ilm, jlm, 4) = &
               this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm) + &
               this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*dot
         end do
      end do

      if (vv <= 0.01d0) then
         do ilm = 1, 9
            if (this%hoh) then
               hhmag(ilm, ilm, 4) = hhmag(ilm, ilm, 4) + this%charge%lattice%symbolic_atoms(it)%potential%cex0(ilm)
            else
               hhmag(ilm, ilm, 4) = hhmag(ilm, ilm, 4) + this%charge%lattice%symbolic_atoms(it)%potential%cx0(ilm)
            end if
         end do
      end if

      do m = 1, 3
         do jlm = 1, 9
            do ilm = 1, 9
               hhmag(ilm, jlm, m) = &
                  (this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm))*momc_ia(m) + &
                  (this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm))*momc_ja(m) + &
                  i_unit*this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*cross(m)
            end do
         end do
      end do

      if (vv > 0.01d0) return
      do m = 1, 3
         do ilm = 1, 9
            if (this%hoh) then
               hhmag(ilm, ilm, m) = hhmag(ilm, ilm, m) + this%charge%lattice%symbolic_atoms(it)%potential%cex1(ilm)*momc_ia(m)
            else
               hhmag(ilm, ilm, m) = hhmag(ilm, ilm, m) + this%charge%lattice%symbolic_atoms(it)%potential%cx1(ilm)*momc_ia(m)
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

   subroutine chbar_nc(this, ia, nr, hmag)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: ia ! Atom number in clust
      integer, intent(in) :: nr ! Number of neighbours considered
      complex(rp), dimension(9, 9, this%charge%lattice%nn_max, 4), intent(out) :: hmag
      ! Local variables
      real(rp) :: r2
      real(rp), dimension(3, this%charge%lattice%kk ) :: cralat ! Clust position times the lattice constant
      ! real(rp), dimension(3, size(this%charge%lattice%cr(1, :))) :: cralat ! Clust position times the lattice constant
      real(rp), dimension(3) :: vet
      real(rp), dimension(9, 9) :: hhh
      integer :: i, j, k, l, m, n, it, jt, ja, nn_max_loc
      integer :: ni, mdir
      integer :: kk ! Clust size number
      real(rp), dimension(:, :), allocatable :: ham_vec
      complex(rp), dimension(9, 9, 4) :: hhmag

      hmag(:, :, :, :) = (0.0d0, 0.0d0)
      !this%hmag(:, :, :, :) = 0.0d0

      r2 = this%charge%lattice%r2
      kk = this%charge%lattice%kk
      !cralat(:, :) = this%charge%lattice%cr(:, :)*this%charge%lattice%alat
      cralat(1:3, 1:kk) = this%charge%lattice%cr(1:3, 1:kk)*this%charge%lattice%alat
      ! print *,'Shape of cralat in chbar_nc', shape(cralat)
      ! print *,'kk:', kk ,'r2:', r2
      allocate(ham_vec(3, nr))
      nn_max_loc = nr

      ! Use clusba directly - we only need local ham_vec, no lattice%sbarvec storage
      call this%charge%lattice%clusba(r2, cralat, ia, kk, kk, nn_max_loc, ham_vec)

      ! DEBUG: Print neighbor vectors from chbar_nc for comparison
      ! if (ia == 1) then  ! Only for first atom
      !    print *, '=== DEBUG chbar_nc: Neighbor vectors for ia=1, nr=', nr
      !    do m = 1, min(5, nr)
      !       print '(A,I3,A,3F12.6)', '  ham_vec[', m, '] = ', ham_vec(:, m)
      !    end do
      ! end if

      !do m=1, nr
      !  print ´(9f10.6)´, real(this%charge%lattice%sbar(:, :, m, ino))
      !end do
      it = this%charge%lattice%iz(ia)
      do m = 1, nr
         ja = this%charge%lattice%nn(ia, m)
         !write(123, *)´ia, ii´, ia, m, this%charge%lattice%nn(ia, m)
         if (m == 1) then
            ja = ia
         end if
         ! print *, "CR test", shape(this%charge%lattice%cr)
         if (ja /= 0) then
            jt = this%charge%lattice%iz(ja)
            if (this%lattice%pbc) then
               call this%lattice%f_wrap_coord_diff(this%lattice%kk,this%lattice%cr*this%lattice%alat,ia,ja,vet)
            else
               vet(:) = (this%charge%lattice%cr(:, ja) - this%charge%lattice%cr(:, ia))*this%charge%lattice%alat
            end if
            !write(123, ´(3f10.6)´) vet(:)
            !write(123, ´(3f10.6)´) this%charge%lattice%sbarvec(:, m)
            !write(123, ´(a, 3i4, 3f10.6)´) ´nn ´, IA, m, JJ, VET(:)
            call this%hmfind(vet, nr, hhh, m, ia, m, ni, ham_vec)
            if (ni == 0) then
               this%charge%lattice%nn(ia, m) = 0
            end if
            call this%ham0m_nc(ia, ja, it, jt, vet, hhh, hhmag)
            do mdir = 1, 4
               !call hcpx(this%hhmag(:, :, mdir), 'cart2sph')
               call hcpx(hhmag(:, :, mdir), 'cart2sph')
               !this%hmag(:, :, m, mdir) = this%hhmag(:, :, mdir)
               hmag(:, :, m, mdir) = hhmag(:, :, mdir)
               !this%hmag(:, :, m, mdir) = hhmag(:, :, mdir)
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

   subroutine hmfind(this, vet, nr, hhh, m, ia, jn, ni, ham_vec)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: m ! Number of the given neighbour
      integer, intent(in) :: ia ! Atom number in clust
      integer, intent(in) :: jn ! ?
      integer, intent(in) :: nr ! Number of neighbours
      real(rp), dimension(3), intent(in) :: vet
      ! Output
      integer, intent(out) :: ni
      real(rp), dimension(9, 9), intent(inout) :: hhh
      real(rp), dimension(3, this%lattice%nn_max), intent(in) :: ham_vec
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
         a1 = (vet(1) - ham_vec(1, i))
         a2 = (vet(2) - ham_vec(2, i))
         a3 = (vet(3) - ham_vec(3, i))
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

   subroutine rs2txt(this)
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

      open (unit=92, file='rs2txt.dat', action='write', iostat=iostatus, status='replace')
      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         do k = 1, nr
            !write(123, *)´ia, ii´, ia, m, this%charge%lattice%nn(ia, m)
            if (k == 1) then
               jj = ia
            else
               jj = this%charge%lattice%nn(ia, k)
            end if
            if (jj /= 0) then
               rij(:) = this%charge%lattice%cr(:, ia) - this%charge%lattice%cr(:, jj)

               if (k == 1) dum = this%ee(:, :, k, ntype) + this%lsham(:, :, ntype) !+ this%enim(:,:,ntype)

               call hcpx(dum(1:9,1:9), 'sph2cart')
               call hcpx(dum(1:9,10:18), 'sph2cart')
               call hcpx(dum(10:18,1:9), 'sph2cart')
               call hcpx(dum(10:18,10:18), 'sph2cart')
               ! Write each matrix element as two columns: real and imaginary part
               do i = 1, 18
                 do j = 1, 18
                  write (92, '(I4,I7, 3F12.8, 2F22.14)') ntype, k, rij , real(this%ee(i, j, k, ntype)), aimag(this%ee(i, j, k, ntype))
                 end do
               end do
            end if
         end do
      end do
      close (92)
   end subroutine rs2txt



   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the on site Hubbard U potential matrix
   !> It used the rotationally invariant formulation bt Liechtensteins (add ref.)
   !> If hubbard_j was set to zero or not provided in the input, then Liechteinsteins formulation reduced to Dudarevs (add ref.).
   !> In this case the provided U becomes U_eff = U-J in Dudarevs paper.
   !> Inputs are number of atoms "na", Slater integrals "f", Hubbard U's "hubbard_u" and Hubbard J's, "hubbard_j"
   !> A generalization of the spdf_Hubbard subroutine that can be used for both bulk and impurity corrections
   !> Builds the effective single-particle potentials for the LDA+U correction for s,p and d orbitals.
   !> (Only the slater integrals F0, F2, F4 and F6 are defined. Machinery works for f orbitals
   !> too, but the corresponding (32x32) basis needs implementation.)
   !> Created by Viktor Frilén and Emil Beiersdorf 27.11.2024
   !> Modified by Viktor Frilén 13.12.2024
   !---------------------------------------------------------------------------
   subroutine calculate_hubbard_u_potential_general(this)
      class(hamiltonian) :: this
      
      ! Local variables
      integer :: nrec ! Number of atoms to perform recursion. This is different between bulk and impurity calculation. Which makes this subroutine compatable with both.
      integer :: i, j ! Orbital indices
      integer :: na ! Atom index
      integer :: ie ! Energy channel index
      integer :: ispin ! Spin index
      integer :: l ! Orbital number, 0,1,2,3 (1,2,3,4) = s,p,d,f (index postions in arrays)
      integer :: cntr ! Counts number of orbitals with a U per atom. 
      real(rp) :: result
      integer :: l_index ! Index for l 1,2,3,4 = s,p,d,f
      integer :: m, m1, m2, m3, m4, m_max, m1_val, m2_val, m3_val, m4_val ! Magnetic quantum numbers
      real(rp), dimension(this%lattice%ntype, 4) :: hub_u, hub_j
      real(rp), dimension(this%lattice%ntype, 4, 4) :: f
      real(rp) :: f0, f2, f4, f6 ! Slater integrals
      ! Local variable array for local density matrix
      real(rp), dimension(this%lattice%ntype, 4, 2, 7, 7) :: LDM ! Local density matrix (LDM), works for spdf orbitals
      real(rp) :: dc ! Double counting term
      real(rp), dimension(4) :: U_energy, dc_energy 
      real(rp), dimension(this%lattice%ntype, 4, 2) :: n_spin ! LDM with m traced out
      real(rp), dimension(this%lattice%ntype, 4) :: n_tot  ! LDM with traced out spin  
      ! Temporary potential that will be put into this%hubbard_u_pot
      real(rp), dimension(this%lattice%ntype, 4, 2, 7, 7) :: hub_pot
      
      type :: ArrayType
         integer, allocatable :: val(:)
      end type ArrayType
      type(ArrayType), dimension(4) :: ms
      type(ArrayType), dimension(this%lattice%ntype) :: l_arr
  
      ms(1)%val = [0]
      ms(2)%val = [-1, 0, 1]
      ms(3)%val = [-2, -1, 0, 1, 2]
      ms(4)%val = [-3, -2, -1, 0, 1, 2, 3]

      this%hubbard_u_pot(:,:,:) = 0.0d0

      LDM(:,:,:,:,:) = 0.0d0
      hub_pot(:,:,:,:,:) = 0.0d0
      n_tot(:,:) = 0.0d0
      n_spin(:,:,:) = 0.0d0
      f(:,:,:) = 0.0d0

      hub_u(:,:) = 0.0d0
      hub_j(:,:) = 0.0d0
      do na = 1, this%lattice%ntype
         do l = 1, 3
            hub_u(na,l) = this%lattice%symbolic_atoms(na)%potential%hubbard_u(l)
            hub_j(na,l) = this%lattice%symbolic_atoms(na)%potential%hubbard_j(l)
         end do
      end do

      ! Maybe do a loop similar to this one (so only the nrec hubbard_u_pot is calculated):
      ! do k = 1, this%lattice%nrec
      !    this%recursion%hamiltonian%hubbard_u_pot(:,:,this%lattice%nbulk + k) = array(:,:,k)
      ! end do

      ! Calculates slater integrals
      do na = 1, this%lattice%ntype
         ! For s, p, d and f orbitals, F0 = U
         ! s-orbital
         f(na,1,1) = hub_u(na,1)
         ! p-orbital J = (1/5)*F2 
         f(na,2,1) = hub_u(na,2)
         f(na,2,2) = hub_j(na,2)*5.0_rp
         ! d-orbital, J = (F2 + F4)/14 with F4/F2 ~ 0.625
         f(na,3,1) = hub_u(na,3)
         f(na,3,2) = 14.0_rp*hub_j(na,3)/1.625_rp 
         f(na,3,3) = 0.625_rp*f(na,3,2)
         !> For f orbitals, J = (286F2 + 195F4 + 250F6)/6435 with F4/F2 ~ 0.67 and F6/F2 ~ 0.49
         f(na,4,1) = hub_u(na,4)
         f(na,4,2) = 6435_rp*hub_j(na,4)/(286_rp + 195_rp*0.67_rp + 250_rp*0.49_rp)
         f(na,4,3) = 0.67_rp*f(na,4,2)
         f(na,4,4) = 0.49_rp*f(na,4,2)               
      end do 
      
      !> Creates an array with each orbital for each atom
      do i = 1, this%lattice%ntype
         cntr = count(abs(hub_u(i,:)) > 1.0E-10) ! Counts orbitals with Hub U for each atom for allocation purposes
         allocate(l_arr(i)%val(cntr))
      end do

      do i = 1, this%lattice%ntype
         cntr = 0
         do j = 1, 4
            if (abs(hub_u(i,j)) > 1.0E-10) then
               cntr = cntr + 1
               l_arr(i)%val(cntr) = j ! Fills the l_arr list with the orbitals that have Hub U
            end if
         end do
      end do       
      
      ! Sets up the local density matrix
      do na = 1, this%lattice%ntype
         do l = 0, 2
            do ispin = 1, 2
               do i = 1, 2*l + 1 !m3
                  do j = 1, 2*l + 1 !m4
                     LDM(na, l + 1, ispin, i, j) = this%lattice%symbolic_atoms(na)%potential%ldm(l + 1, ispin, i, j)
                  end do
               end do
            end do
         end do
      end do
      
      ! Builds the Hubbard U+J potential 
      print *, ''
      print *, '-----------------------------------------------------------------------------------------------------------------'
      print *, 'Calculate Hubbard U potential for impurity and bulk:'
      do na = 1, this%lattice%ntype
         print *, ' Atom ', na
         print *, ''
         ! Calculates traces of the local density matrix, n_spin is the trace in m, n_tot is trace in spin of n_spin.
         do l = 1, 3 !0 to 3 but indexing starts on 1
            do ispin = 1, 2
               n_spin(na, l, ispin) = trace(LDM(na, l, ispin, :, :)) 
               n_tot(na, l) = n_tot(na, l) + n_spin(na, l, ispin)
            end do
         end do
         U_energy = 0.0d0

         ! l_arr(na)%val(l) gives the orbital index for each atom, for s,p,d its value is 1,2,3
         ! l_arr(na)%val(l)-1 gives the actual l-value for each atom, for s,p,d its value is 0,1,2   
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
                     ! 
                     ! AMF (Around Mean Field) double counting:
                     ! E_DC^AMF = (1/2)*U*n^2 - (1/2)*J*(n_up^2 + n_down^2)
                     ! V_DC^AMF = dE_DC/dn_sigma = U*n - J*n_sigma
                     !  if (m1 == m2) then
                     !     hub_pot(na, l_index, ispin, m1, m2) = hub_pot(na, l_index, ispin, m1, m2) &
                     !     - hub_u(na,l_index)*n_tot(na,l_index) + hub_j(na,l_index)*n_spin(na,l_index,ispin)
                     !  end if
                     ! 
                     ! FLL (Fully Localized Limit) double counting:
                     ! E_DC^FLL = (1/2)*U*n*(n-1) - (1/2)*J*[n_up*(n_up-1) + n_down*(n_down-1)]
                     ! V_DC^FLL = dE_DC/dn_sigma = U*(n - 1/2) - J*(n_sigma - 1/2)
                     if (m1 == m2) then
                        hub_pot(na,l_index,ispin,m1,m2) = hub_pot(na,l_index,ispin,m1,m2) &
                              - hub_u(na,l_index)*(n_tot(na,l_index) - 0.5_rp) &
                              + hub_j(na,l_index)*(n_spin(na,l_index,ispin) - 0.5_rp)  ! Corrected: minus sign for J term
                      end if
                  end do
               end do
            end do
            print *, 'n_tot for atom', na, ':', n_tot(na, :)
            print *, 'n_spin up:', n_spin(na, :, 1)
            print *, 'n_spin down:', n_spin(na, :, 2)
            print *, 'Sample diagonal V_U:', hub_pot(na, 3, 1, 3, 3)
            print *, 'V_U diagonal spin-up (m=0):', hub_pot(na, 3, 1, 3, 3)
            print *, 'V_U diagonal spin-down (m=0):', hub_pot(na, 3, 2, 3, 3)
         end do  

      !> Puts hub_pot into global hubbard_u_pot (only done for spd-orbitals)
         do l = 0, 2
            do i = 1, 2*l + 1
               do j = 1, 2*l + 1
                  this%hubbard_u_pot(l**2+i, l**2+j, na) = hub_pot(na, l+1, 1, i, j)
                  this%hubbard_u_pot(l**2+i+9, l**2+j+9, na) = hub_pot(na, l+1, 2, i, j)
               end do
            end do
            ! Double counting energy correction
            ! 
            ! AMF (Around Mean Field):
            ! E_DC^AMF = (1/2)*U*n^2 - (1/2)*J*(n_up^2 + n_down^2)
            ! Simplifies to: E_DC^AMF = (1/2)*U*n*(n-1) + (1/2)*U*n - (1/2)*J*[n_up*(n_up-1) + n_down*(n_down-1)] - (1/2)*J*(n_up + n_down)
            !                         = (1/2)*U*n*(n-1) - (1/2)*J*[n_up*(n_up-1) + n_down*(n_down-1)] + (1/2)*(U - J)*n
            !  dc_energy(l+1) = 0.5_rp*hub_u(na,l+1)*n_tot(na,l+1)**2 &
            !                 - 0.5_rp*hub_j(na,l+1)*(n_spin(na,l+1,1)**2 + n_spin(na,l+1,2)**2)
            ! 
            ! FLL (Fully Localized Limit):
            ! E_DC^FLL = (1/2)*U*n*(n-1) - (1/2)*J*[n_up*(n_up-1) + n_down*(n_down-1)]
             dc_energy(l+1) = 0.5_rp * ( &
               hub_u(na,l+1) * n_tot(na,l+1) * (n_tot(na,l+1) - 1.0_rp) &
             - hub_j(na,l+1) * ( n_spin(na,l+1,1)*(n_spin(na,l+1,1) - 1.0_rp) &
                               + n_spin(na,l+1,2)*(n_spin(na,l+1,2) - 1.0_rp) ) )

            
      ! This part could be modified to work for imputrities
            do i = 1, size(l_arr(na)%val)
               if (l+1 == l_arr(na)%val(i)) then ! Only prints U_energy if that orbital has a specified U
                  print *, ' Orbital ', this%orb_conv(l+1), ' :'
                  print *, '  (U, dc, total) [eV] = ', U_energy(l+1)*ry2ev, dc_energy(l+1)*ry2ev, U_energy(l+1)*ry2ev - dc_energy(l+1)*ry2ev
                  print *, ''
               end if
            end do
         end do
         
      end do 
      print *, '-----------------------------------------------------------------------------------------------------------------'
      print *,''
      
   end subroutine calculate_hubbard_u_potential_general
end module hamiltonian_mod
