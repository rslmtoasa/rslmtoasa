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

   use mpi_mod, only: rank
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
   use spectrum_bounds_mod, only: compute_spectrum_bounds, bounds, normalize_bounds_algorithm, select_bounds_interval, apply_bounds_scaling
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   use basis_mod, only: nb, norb, spin_off
   use block_sparse_operator_mod, only: block_sparse_operator
   implicit none

   private

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
      complex(rp), dimension(:, :, :, :), allocatable :: hmag, hxc
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
      !> Velocity operators
      complex(rp), dimension(:, :, :, :), allocatable :: v_a, v_b, js_a, jl_a, vo_a, vo_b, jso_a, jlo_a
      character(len=10) :: js_alpha, jl_alpha
      real(rp), dimension(3) :: v_alpha, v_beta, q_ss
      real(rp) :: theta_ss
      real(rp), dimension(:), allocatable :: velocity_scale
      character(len=16) :: hubbard_u_potential_form
      !> Sparse Real Space Hamiltonian
      complex(rp), dimension(:, :), allocatable :: h_sparse
      !> Exported block operators used by recursion backends
      type(block_sparse_operator) :: exported_h_plain, exported_h_full
      type(block_sparse_operator) :: exported_v_a, exported_v_b, exported_vo_a, exported_vo_b
      !> On-site potential correction for LDA+U (+J) in spin-orbital basis
      real(rp), dimension(:, :, :), allocatable :: hubbard_u_pot
      !> Enable +U correction when any U/J is provided on symbolic atoms
      logical :: hubbard_u_general_check = .false.
      !> Optional impurity-only Hubbard inputs (eV in input.nml, stored as Ry here)
      real(rp), dimension(:, :), allocatable :: hubbard_u_impurity, hubbard_j_impurity
      logical :: hubbard_u_impurity_check = .false.
      !> Optional self-consistent Hubbard-U mask (itype,l) with l=1..4 => s,p,d,f.
      !> 0 = disabled, 1 = enable self-consistent U update for this channel.
      integer, dimension(:, :), allocatable :: hubbard_u_sc
      logical :: hubbard_u_sc_check = .false.
      !> Optional intersite Hubbard-V input and corresponding matrix correction
      real(rp), dimension(:, :, :, :), allocatable :: hubbard_v, hubbard_v_pot
      logical :: hubbard_v_check = .false.
      !> Spectrum bounds for Chebyshev scaling
      type(bounds) :: bounds
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
      procedure :: export_block_operator
      procedure :: export_velocity_operator
      procedure :: torque_operator_collinear
      procedure :: rs2pao
      procedure :: chbar_nc
      procedure :: ham0m_nc
      procedure :: hmfind
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: rotate_to_local_axis
      procedure :: rotate_from_local_axis
      procedure :: calculate_hubbard_u_potential_general
      procedure :: calculate_hubbard_v_potential
      procedure :: compute_hamiltonian_bounds
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
      if (allocated(this%v_a)) call g_safe_alloc%deallocate('hamiltonian.v_a', this%v_a)
      if (allocated(this%v_b)) call g_safe_alloc%deallocate('hamiltonian.v_b', this%v_b)
      if (allocated(this%vo_a)) call g_safe_alloc%deallocate('hamiltonian.vo_a', this%vo_a)
      if (allocated(this%vo_b)) call g_safe_alloc%deallocate('hamiltonian.vo_b', this%vo_b)
      if (allocated(this%jso_a)) call g_safe_alloc%deallocate('hamiltonian.jso_a', this%jso_a)
      if (allocated(this%jlo_a)) call g_safe_alloc%deallocate('hamiltonian.jlo_a', this%jlo_a)
      if (allocated(this%h_sparse)) call g_safe_alloc%deallocate('hamiltonian.h_sparse', this%h_sparse)
      if (allocated(this%velocity_scale)) call g_safe_alloc%deallocate('hamiltonian.velocity_scale', this%velocity_scale)
      if (allocated(this%hxc)) call g_safe_alloc%deallocate('hamiltonian.hxc', this%hxc)
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
      if (allocated(this%v_a)) deallocate(this%v_a)
      if (allocated(this%v_b)) deallocate(this%v_b)
      if (allocated(this%vo_a)) deallocate(this%vo_a)
      if (allocated(this%vo_b)) deallocate(this%vo_b)
      if (allocated(this%jso_a)) deallocate(this%jso_a)
      if (allocated(this%jlo_a)) deallocate(this%jlo_a)
      if (allocated(this%h_sparse)) deallocate(this%h_sparse)
      if (allocated(this%velocity_scale)) deallocate(this%velocity_scale)
      if (allocated(this%hxc)) deallocate(this%hxc)
      if (allocated(this%hubbard_u_pot)) deallocate(this%hubbard_u_pot)
      if (allocated(this%hubbard_u_impurity)) deallocate(this%hubbard_u_impurity)
      if (allocated(this%hubbard_j_impurity)) deallocate(this%hubbard_j_impurity)
      if (allocated(this%hubbard_u_sc)) deallocate(this%hubbard_u_sc)
      if (allocated(this%hubbard_v)) deallocate(this%hubbard_v)
      if (allocated(this%hubbard_v_pot)) deallocate(this%hubbard_v_pot)
#endif
   end subroutine destructor

   ! Member functions
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read parameters from input file
   !---------------------------------------------------------------------------
   subroutine build_from_file(this)
      class(hamiltonian), intent(inout) :: this

      ! variables associated with the reading processes
      integer :: iostatus, funit, i, j, l, li, lj, n_l_in
      integer :: nimp_in
      integer :: na
      integer :: ntype_nml, nrec_nml
      logical :: legacy_uj_present
      logical :: has_hubbard_general, has_legacy_uj, has_impurity_uj, has_hubbard_v, has_hubbard_u_sc
      character(len=1) :: orbch
      real(rp), parameter :: hubbard_sentinel = -9.87654321e30_rp
      integer, parameter :: hubbard_sc_sentinel = -999

      include 'include_codes/namelists/hamiltonian.f90'

      hoh = this%hoh
      local_axis = this%local_axis
      orb_pol = this%orb_pol
      v_alpha(:) = this%v_alpha(:)
      v_beta(:) = this%v_beta(:)
      q_ss(:) = this%q_ss(:)
      theta_ss = this%theta_ss
      js_alpha = this%js_alpha
      jl_alpha = this%jl_alpha
      call move_alloc(this%velocity_scale, velocity_scale)
      hubbard_u_potential_form = this%hubbard_u_potential_form
      bounds_algorithm = this%bounds%algorithm
      bounds_scaling = this%bounds%scaling

      ! Robust namelist buffers: pre-allocate allocatable inputs so namelist
      ! read does not depend on compiler-specific auto-allocation behavior.
      ntype_nml = max(1, this%lattice%ntype)
      nrec_nml = max(1, this%lattice%nrec)

      if (.not. allocated(hubbard_u_general)) allocate(hubbard_u_general(ntype_nml, 4))
      if (.not. allocated(hubbard_j_general)) allocate(hubbard_j_general(ntype_nml, 4))
      if (.not. allocated(hubbard_u)) allocate(hubbard_u(nrec_nml, 4))
      if (.not. allocated(hubbard_j)) allocate(hubbard_j(nrec_nml, 4))
      if (.not. allocated(hubbard_u_impurity)) allocate(hubbard_u_impurity(nrec_nml, 4))
      if (.not. allocated(hubbard_j_impurity)) allocate(hubbard_j_impurity(nrec_nml, 4))
      if (.not. allocated(hubbard_u_sc)) allocate(hubbard_u_sc(ntype_nml, 4))
      if (.not. allocated(hubbard_v)) allocate(hubbard_v(ntype_nml, ntype_nml, 4, 4))
      if (.not. allocated(uj_orb)) allocate(uj_orb(nrec_nml))

      hubbard_u_general(:, :) = hubbard_sentinel
      hubbard_j_general(:, :) = hubbard_sentinel
      hubbard_u(:, :) = hubbard_sentinel
      hubbard_j(:, :) = hubbard_sentinel
      hubbard_u_impurity(:, :) = hubbard_sentinel
      hubbard_j_impurity(:, :) = hubbard_sentinel
      hubbard_v(:, :, :, :) = hubbard_sentinel
      hubbard_u_sc(:, :) = hubbard_sc_sentinel
      uj_orb(:) = ''

      ! Reading
      open (newunit=funit, file=this%control%fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(this%control%fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=hamiltonian, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      this%hoh = hoh
      this%local_axis = local_axis
      this%orb_pol = orb_pol
      this%v_alpha(:) = v_alpha(:)
      this%v_beta(:) = v_beta(:)
      this%q_ss(:) = q_ss(:)
      this%theta_ss = theta_ss
      this%js_alpha = js_alpha
      this%jl_alpha = jl_alpha
      this%hubbard_u_potential_form = lower(trim(hubbard_u_potential_form))
      if (this%hubbard_u_potential_form /= 'liechtenstein' .and. this%hubbard_u_potential_form /= 'acbn0') then
         call g_logger%fatal("Invalid hubbard_u_potential_form. Use 'liechtenstein' or 'acbn0'.", __FILE__, __LINE__)
      end if
      this%bounds%algorithm = lower(trim(bounds_algorithm))
      if (len_trim(this%bounds%algorithm) == 0) this%bounds%algorithm = 'none'
      if (bounds_scaling > 0.0_rp) then
         this%bounds%scaling = bounds_scaling
      else
         this%bounds%scaling = 1.05_rp
      end if
      call move_alloc(velocity_scale, this%velocity_scale)

      has_hubbard_general = any(hubbard_u_general /= hubbard_sentinel) .or. any(hubbard_j_general /= hubbard_sentinel)
      has_legacy_uj = any(uj_orb /= '') .and. (any(hubbard_u /= hubbard_sentinel) .or. any(hubbard_j /= hubbard_sentinel))
      has_impurity_uj = any(hubbard_u_impurity /= hubbard_sentinel) .or. any(hubbard_j_impurity /= hubbard_sentinel)
      has_hubbard_v = any(hubbard_v /= hubbard_sentinel)
      has_hubbard_u_sc = any(hubbard_u_sc /= hubbard_sc_sentinel)

      ! Optional Hubbard inputs from &hamiltonian are provided in eV.
      ! Converted below to internal Ry units:
      ! - hubbard_u_general / hubbard_j_general
      ! - legacy hubbard_u / hubbard_j via uj_orb mapping
      ! - hubbard_u_impurity / hubbard_j_impurity
      ! - hubbard_v (intersite)
      if (has_hubbard_general) then
         n_l_in = min(size(hubbard_u_general, 2), size(hubbard_j_general, 2))
         do i = 1, min(this%lattice%ntype, size(hubbard_u_general, 1), size(hubbard_j_general, 1))
            do l = 1, min(n_l_in, size(this%lattice%symbolic_atoms(i)%potential%hubbard_u))
               if (hubbard_u_general(i, l) == hubbard_sentinel .or. hubbard_j_general(i, l) == hubbard_sentinel) cycle
               this%lattice%symbolic_atoms(i)%potential%hubbard_u(l) = hubbard_u_general(i, l)/ry2ev
               this%lattice%symbolic_atoms(i)%potential%hubbard_j(l) = hubbard_j_general(i, l)/ry2ev
            end do
         end do
      end if

      ! Legacy LDA+U(+J) input path from lda_u branch:
      ! map (hubbard_u/hubbard_j + uj_orb) into per-atom per-l potential arrays.
      legacy_uj_present = .false.
      if (has_legacy_uj) then
         do i = 1, min(this%lattice%nrec, size(uj_orb), size(hubbard_u, 1), size(hubbard_j, 1))
            na = this%lattice%nbulk + i
            if (na < 1 .or. na > this%lattice%ntype) cycle
            do j = 1, len_trim(uj_orb(i))
               if (j > size(hubbard_u, 2) .or. j > size(hubbard_j, 2)) exit
               orbch = uj_orb(i)(j:j)
               select case (orbch)
               case ('s', 'S')
                  l = 1
               case ('p', 'P')
                  l = 2
               case ('d', 'D')
                  l = 3
               case ('f', 'F')
                  l = 4
               case default
                  cycle
               end select
               if (l <= size(this%lattice%symbolic_atoms(na)%potential%hubbard_u)) then
                  if (hubbard_u(i, j) == hubbard_sentinel .or. hubbard_j(i, j) == hubbard_sentinel) cycle
                  this%lattice%symbolic_atoms(na)%potential%hubbard_u(l) = hubbard_u(i, j)/ry2ev
                  this%lattice%symbolic_atoms(na)%potential%hubbard_j(l) = hubbard_j(i, j)/ry2ev
                  if (abs(hubbard_u(i, j)) > 1.0e-10_rp .or. abs(hubbard_j(i, j)) > 1.0e-10_rp) legacy_uj_present = .true.
               end if
            end do
         end do
      end if
      if (legacy_uj_present) then
         call g_logger%info('Legacy uj_orb + hubbard_u/hubbard_j input detected and mapped to symbolic-atom Hubbard channels.', __FILE__, __LINE__)
      end if

      this%hubbard_u_impurity_check = .false.
      if (has_impurity_uj) then
         nimp_in = min(size(hubbard_u_impurity, 1), size(hubbard_j_impurity, 1))
         if (nimp_in > 0) then
            this%hubbard_u_impurity(:, :) = 0.0_rp
            this%hubbard_j_impurity(:, :) = 0.0_rp
            do i = 1, min(this%lattice%nrec, nimp_in)
               do l = 1, min(4, size(hubbard_u_impurity, 2), size(hubbard_j_impurity, 2))
                  if (hubbard_u_impurity(i, l) == hubbard_sentinel .or. hubbard_j_impurity(i, l) == hubbard_sentinel) cycle
                  this%hubbard_u_impurity(i, l) = hubbard_u_impurity(i, l)/ry2ev
                  this%hubbard_j_impurity(i, l) = hubbard_j_impurity(i, l)/ry2ev
               end do
            end do
         end if
      end if

      ! Apply impurity-specific U/J to recursive atoms (nbulk+1:ntype)
      if (this%lattice%nrec > 0) then
         do i = 1, this%lattice%nrec
            if (maxval(abs(this%hubbard_u_impurity(i, :))) > 1.0e-10_rp .or. &
                maxval(abs(this%hubbard_j_impurity(i, :))) > 1.0e-10_rp) then
               this%hubbard_u_impurity_check = .true.
               do l = 1, min(4, size(this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%hubbard_u))
                  this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%hubbard_u(l) = this%hubbard_u_impurity(i, l)
                  this%lattice%symbolic_atoms(this%lattice%nbulk + i)%potential%hubbard_j(l) = this%hubbard_j_impurity(i, l)
               end do
            end if
         end do
      end if

      ! Optional self-consistent U selector (0/1 mask by atom type and l-channel).
      this%hubbard_u_sc_check = .false.
      this%hubbard_u_sc(:, :) = 0
      if (has_hubbard_u_sc) then
         do i = 1, min(size(this%hubbard_u_sc, 1), size(hubbard_u_sc, 1))
            do l = 1, min(size(this%hubbard_u_sc, 2), size(hubbard_u_sc, 2))
               if (hubbard_u_sc(i, l) == hubbard_sc_sentinel) cycle
               if (hubbard_u_sc(i, l) < 0 .or. hubbard_u_sc(i, l) > 1) then
                  call g_logger%fatal('Invalid hubbard_u_sc value. Only 0 or 1 is allowed.', __FILE__, __LINE__)
               end if
               this%hubbard_u_sc(i, l) = hubbard_u_sc(i, l)
               if (this%hubbard_u_sc(i, l) == 1) this%hubbard_u_sc_check = .true.
            end do
         end do
      end if

      this%hubbard_v_check = .false.
      if (has_hubbard_v) then
         this%hubbard_v(:, :, :, :) = 0.0_rp
         do i = 1, min(size(this%hubbard_v, 1), size(hubbard_v, 1))
            do j = 1, min(size(this%hubbard_v, 2), size(hubbard_v, 2))
               do li = 1, min(size(this%hubbard_v, 3), size(hubbard_v, 3))
                  do lj = 1, min(size(this%hubbard_v, 4), size(hubbard_v, 4))
                     if (hubbard_v(i, j, li, lj) == hubbard_sentinel) cycle
                     this%hubbard_v(i, j, li, lj) = hubbard_v(i, j, li, lj)/ry2ev
                  end do
               end do
            end do
         end do

      ! Enforce V_ji(lj,li) = V_ij(li,lj) when only one side is provided.
      ! this%hubbard_v is already in Ry at this point.
      do i = 1, size(this%hubbard_v, 1)
            do j = 1, size(this%hubbard_v, 2)
               do li = 1, size(this%hubbard_v, 3)
                  do lj = 1, size(this%hubbard_v, 4)
                     if (abs(this%hubbard_v(i, j, li, lj)) > 1.0e-10_rp) then
                        this%hubbard_v_check = .true.
                        if (j <= size(this%hubbard_v, 1) .and. i <= size(this%hubbard_v, 2)) then
                           if (abs(this%hubbard_v(j, i, lj, li)) <= 1.0e-10_rp) then
                              this%hubbard_v(j, i, lj, li) = this%hubbard_v(i, j, li, lj)
                           end if
                        end if
                     end if
                  end do
               end do
            end do
         end do
      end if

      this%hubbard_u_general_check = .false.
      do i = 1, this%lattice%ntype
         if (maxval(abs(this%lattice%symbolic_atoms(i)%potential%hubbard_u(:))) > 1.0e-10_rp) then
            this%hubbard_u_general_check = .true.
         end if
         if (maxval(abs(this%lattice%symbolic_atoms(i)%potential%hubbard_j(:))) > 1.0e-10_rp) then
            this%hubbard_u_general_check = .true.
         end if
      end do

      ! Self-consistent U and fixed U/J are mutually exclusive in this implementation.
      if (this%hubbard_u_sc_check .and. this%hubbard_u_general_check) then
         call g_logger%fatal('Both hubbard_u_sc and explicit hubbard_u/hubbard_j are set. Use only one mode.', __FILE__, __LINE__)
      end if

      if (rank == 0) then
         if (this%hubbard_u_general_check .or. this%hubbard_u_sc_check .or. this%hubbard_v_check) then
         call g_logger%info('HUBBARD summary: form='//trim(this%hubbard_u_potential_form)// &
                            ' fixed_UJ='//merge('T', 'F', this%hubbard_u_general_check)// &
                            ' sc_U='//merge('T', 'F', this%hubbard_u_sc_check)// &
                            ' V='//merge('T', 'F', this%hubbard_v_check), __FILE__, __LINE__)
         end if
      end if
      if (this%hubbard_u_general_check) then
         block
            integer :: lch, nch
            character(len=512) :: u_msg, j_msg
            do i = 1, this%lattice%ntype
               u_msg = ''
               j_msg = ''
               nch = size(this%lattice%symbolic_atoms(i)%potential%hubbard_u)
               do lch = 1, nch
                  u_msg = trim(u_msg)//' '//fmt('f10.6', this%lattice%symbolic_atoms(i)%potential%hubbard_u(lch))
               end do
               nch = size(this%lattice%symbolic_atoms(i)%potential%hubbard_j)
               do lch = 1, nch
                  j_msg = trim(j_msg)//' '//fmt('f10.6', this%lattice%symbolic_atoms(i)%potential%hubbard_j(lch))
               end do
               if (rank == 0) then
                  call g_logger%info('HUBBARD fixed U/J type='//fmt('i4', i)//' [Ry] U='//trim(u_msg)//' J='//trim(j_msg), __FILE__, __LINE__)
               end if
            end do
         end block
      end if
      if (this%hubbard_u_sc_check) then
         do i = 1, size(this%hubbard_u_sc, 1)
            if (rank == 0) then
               call g_logger%info('HUBBARD sc_U mask type='//fmt('i4', i)//' [s p d f]='// &
                                  fmt('i2', this%hubbard_u_sc(i, 1))//' '//fmt('i2', this%hubbard_u_sc(i, 2))//' '// &
                                  fmt('i2', this%hubbard_u_sc(i, 3))//' '//fmt('i2', this%hubbard_u_sc(i, 4)), __FILE__, __LINE__)
            end if
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

      call this%exported_h_plain%clear()
      call this%exported_h_full%clear()
      call this%exported_v_a%clear()
      call this%exported_v_b%clear()
      call this%exported_vo_a%clear()
      call this%exported_vo_b%clear()

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('hamiltonian.lsham', this%lsham, (/nb, nb, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.tmat', this%tmat, (/nb, nb, 3, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hhmag', this%hhmag, (/norb, norb, 4/))
      call g_safe_alloc%allocate('hamiltonian.hmag', this%hmag, (/norb, norb, this%charge%lattice%kk, 4/))
      call g_safe_alloc%allocate('hamiltonian.ee', this%ee, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hall', this%hall, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.hall_glob', this%hall_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.ee_glob', this%ee_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      !if (hoh) then
      call g_safe_alloc%allocate('hamiltonian.eeo', this%eeo, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.eeoee', this%eeoee, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hallo', this%hallo, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.obarm', this%obarm, (/nb, nb, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.enim', this%enim, (/nb, nb, this%charge%lattice%ntype/))
      !end if
      !if (local_axis)  then
      call g_safe_alloc%allocate('hamiltonian.hall_glob', this%hall_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.ee_glob', this%ee_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      !if (hoh) then
      call g_safe_alloc%allocate('hamiltonian.ee0_glob', this%eeo_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hallo_glob', this%hallo_glob, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%nmax/))
      call g_safe_alloc%allocate('hamiltonian.enim_glob', this%enim_glob, (/nb, nb, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.v_a', this%v_a, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.v_b', this%v_b, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.vo_a', this%vo_a, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.vo_b', this%vo_b, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.js_a', this%js_a, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.jl_a', this%jl_a, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.jso_a', this%jso_a, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.jlo_a', this%jlo_a, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.velocity_scale', this%velocity_scale, (/this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hxc', this%hxc, (/nb, nb, (this%charge%lattice%nn(1, 1) + 1), this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hubbard_u_pot', this%hubbard_u_pot, (/nb, nb, this%charge%lattice%ntype/))
      call g_safe_alloc%allocate('hamiltonian.hubbard_u_impurity', this%hubbard_u_impurity, (/max(1, this%charge%lattice%nrec), 4/))
      call g_safe_alloc%allocate('hamiltonian.hubbard_j_impurity', this%hubbard_j_impurity, (/max(1, this%charge%lattice%nrec), 4/))
      call g_safe_alloc%allocate('hamiltonian.hubbard_u_sc', this%hubbard_u_sc, (/max(1, this%charge%lattice%ntype), 4/))
      ! Hubbard-V is indexed by symbolic atom type pair (itype,jtype), not by impurity index.
      call g_safe_alloc%allocate('hamiltonian.hubbard_v', this%hubbard_v, (/max(1, this%charge%lattice%ntype), max(1, this%charge%lattice%ntype), 4, 4/))
      call g_safe_alloc%allocate('hamiltonian.hubbard_v_pot', this%hubbard_v_pot, (/nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype/))
      !end if
      !end if
#else
      allocate (this%lsham(nb, nb, this%charge%lattice%ntype))
      allocate (this%tmat(nb, nb, 3, this%charge%lattice%ntype))
      allocate (this%hhmag(norb, norb, 4), this%hmag(norb, norb, this%charge%lattice%kk, 4))
      allocate (this%ee(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hall(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%hxc(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      !if (this%hoh) then
      allocate (this%eeo(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%eeoee(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hallo(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%obarm(nb, nb, this%charge%lattice%ntype))
      allocate (this%enim(nb, nb, this%charge%lattice%ntype))
      !end if
      !if (this%local_axis) then
      allocate (this%ee_glob(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hall_glob(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      !if (this%hoh) then
      allocate (this%eeo_glob(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%hallo_glob(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%nmax))
      allocate (this%enim_glob(nb, nb, this%charge%lattice%ntype))
      ! Velocity operators
      allocate (this%v_a(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%v_b(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%vo_a(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%vo_b(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%js_a(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%jl_a(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%jso_a(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%jlo_a(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      allocate (this%velocity_scale(this%charge%lattice%ntype))
      allocate (this%hubbard_u_pot(nb, nb, this%charge%lattice%ntype))
      allocate (this%hubbard_u_impurity(max(1, this%charge%lattice%nrec), 4))
      allocate (this%hubbard_j_impurity(max(1, this%charge%lattice%nrec), 4))
      allocate (this%hubbard_u_sc(max(1, this%charge%lattice%ntype), 4))
      ! Hubbard-V is indexed by symbolic atom type pair (itype,jtype), not by impurity index.
      allocate (this%hubbard_v(max(1, this%charge%lattice%ntype), max(1, this%charge%lattice%ntype), 4, 4))
      allocate (this%hubbard_v_pot(nb, nb, (maxval(this%charge%lattice%nn(:, 1)) + 1), this%charge%lattice%ntype))
      !end if
      !end if
#endif

      this%lsham(:, :, :) = 0.0d0
      this%tmat(:, :, :, :) = 0.0d0
      this%hhmag(:, :, :) = 0.0d0
      this%hmag(:, :, :, :) = 0.0d0
      this%hall(:, :, :, :) = 0.0d0
      this%hxc(:, :, :, :) = 0.0d0
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
      this%hubbard_u_pot(:, :, :) = 0.0d0
      this%bounds%algorithm = 'none'
      this%bounds%scaling = 1.05_rp
      this%bounds%e_min = 0.0_rp
      this%bounds%e_max = 0.0_rp
      this%hubbard_u_impurity(:, :) = 0.0d0
      this%hubbard_j_impurity(:, :) = 0.0d0
      this%hubbard_u_sc(:, :) = 0
      this%hubbard_v(:, :, :, :) = 0.0d0
      this%hubbard_v_pot(:, :, :, :) = 0.0d0
      this%hoh = .false.
      this%local_axis = .false.
      this%orb_pol = .false.
      this%v_alpha(:) = [1, 0, 0]
      this%v_beta(:) = [1, 0, 0]
      this%q_ss(:) = [0.0_rp, 0.0_rp, 0.0_rp]
      this%theta_ss = 0.0_rp
      this%js_alpha = 'z'
      this%jl_alpha = 'z'
      this%hubbard_u_potential_form = 'liechtenstein'
      this%hubbard_u_general_check = .false.
      this%hubbard_u_impurity_check = .false.
      this%hubbard_u_sc_check = .false.
      this%hubbard_v_check = .false.
   end subroutine restore_to_default

   subroutine block_to_sparse(this)
      class(hamiltonian), intent(inout) :: this
      integer :: row, entry, col, i0, j0
      character(len=:), allocatable :: base

      if (this%hoh) then
         call this%export_block_operator('hamiltonian_hoh', this%exported_h_full)
         if (allocated(this%h_sparse)) deallocate (this%h_sparse)
         allocate (this%h_sparse(this%lattice%kk*nb, this%lattice%kk*nb))
         this%h_sparse(:, :) = (0.0_rp, 0.0_rp)
         do row = 1, this%exported_h_full%n_sites
            i0 = (row - 1)*nb + 1
            do entry = this%exported_h_full%row_ptr(row), this%exported_h_full%row_ptr(row + 1) - 1
               col = this%exported_h_full%col_ind(entry)
               j0 = (col - 1)*nb + 1
               this%h_sparse(i0:i0 + nb - 1, j0:j0 + nb - 1) = this%exported_h_full%blocks(:, :, entry)
            end do
         end do
      else
         call this%export_block_operator('hamiltonian', this%exported_h_plain)
         if (allocated(this%h_sparse)) deallocate (this%h_sparse)
         allocate (this%h_sparse(this%lattice%kk*nb, this%lattice%kk*nb))
         this%h_sparse(:, :) = (0.0_rp, 0.0_rp)
         do row = 1, this%exported_h_plain%n_sites
            i0 = (row - 1)*nb + 1
            do entry = this%exported_h_plain%row_ptr(row), this%exported_h_plain%row_ptr(row + 1) - 1
               col = this%exported_h_plain%col_ind(entry)
               j0 = (col - 1)*nb + 1
               this%h_sparse(i0:i0 + nb - 1, j0:j0 + nb - 1) = this%exported_h_plain%blocks(:, :, entry)
            end do
         end do
      end if

      if (this%control%export_hamiltonian) then
         base = trim(this%control%export_hamiltonian_path)
         if (len_trim(base) > 0) then
            if (this%hoh) then
               call this%exported_h_full%write_text(base)
            else
               call this%exported_h_plain%write_text(base)
            end if
         end if
      end if
   end subroutine block_to_sparse

   subroutine export_block_operator(this, operator_name, op)
      class(hamiltonian), intent(inout) :: this
      character(len=*), intent(in) :: operator_name
      type(block_sparse_operator), intent(inout) :: op
      integer :: site, nr, m, nnzb, entry, neighbor, ih
      character(len=:), allocatable :: op_name
      complex(rp) :: block(nb, nb)

      op_name = lower(trim(operator_name))
      call op%clear()

      op%block_dim = nb
      op%n_sites = this%lattice%kk
      op%kind = op_name
      op%includes_overlap = (op_name == 'hamiltonian_overlap')
      op%includes_soc = (op_name == 'hamiltonian' .or. op_name == 'hamiltonian_hoh' .or. op_name == 'hamiltonian_shift')
      op%includes_enim = (op_name == 'hamiltonian_hoh' .or. op_name == 'hamiltonian_shift')

      allocate(op%site_types(this%lattice%kk))
      op%site_types(:) = this%lattice%iz(1:this%lattice%kk)

      nnzb = 0
      do site = 1, this%lattice%kk
         nr = max(1, this%lattice%nn(site, 1))
         nnzb = nnzb + 1
         do m = 2, nr
            neighbor = this%lattice%nn(site, m)
            if (neighbor /= 0) nnzb = nnzb + 1
         end do
      end do

      op%nnzb = nnzb
      allocate(op%row_ptr(this%lattice%kk + 1))
      allocate(op%col_ind(nnzb))
      allocate(op%blocks(nb, nb, nnzb))

      entry = 1
      op%row_ptr(1) = 1
      do site = 1, this%lattice%kk
         nr = max(1, this%lattice%nn(site, 1))

         call build_operator_block(site, 1, site, op_name, block)
         op%col_ind(entry) = site
         op%blocks(:, :, entry) = block(:, :)
         entry = entry + 1

         do m = 2, nr
            neighbor = this%lattice%nn(site, m)
            if (neighbor == 0) cycle
            call build_operator_block(site, m, neighbor, op_name, block)
            op%col_ind(entry) = neighbor
            op%blocks(:, :, entry) = block(:, :)
            entry = entry + 1
         end do
         op%row_ptr(site + 1) = entry
      end do

   contains

      subroutine build_operator_block(site, neighbor_idx, neighbor_site, name, block_out)
         integer, intent(in) :: site, neighbor_idx, neighbor_site
         character(len=*), intent(in) :: name
         complex(rp), intent(out) :: block_out(nb, nb)
         integer :: site_type

         block_out(:, :) = (0.0_rp, 0.0_rp)
         site_type = this%lattice%iz(site)

         if (site <= this%lattice%nmax .and. this%lattice%nmax /= 0) then
            select case (name)
            case ('hamiltonian')
               block_out(:, :) = this%hall(:, :, neighbor_idx, site)
               if (neighbor_idx == 1) block_out(:, :) = block_out(:, :) + this%lsham(:, :, site_type)
            case ('hamiltonian_core')
               block_out(:, :) = this%hall(:, :, neighbor_idx, site)
            case ('hamiltonian_overlap')
               block_out(:, :) = this%hallo(:, :, neighbor_idx, site)
            case ('hamiltonian_shift')
               if (neighbor_idx == 1) block_out(:, :) = this%enim(:, :, site_type) + this%lsham(:, :, site_type)
            case ('hamiltonian_hoh')
               block_out(:, :) = this%hall(:, :, neighbor_idx, site) - this%hallo(:, :, neighbor_idx, site)
               if (neighbor_idx == 1) block_out(:, :) = block_out(:, :) + this%enim(:, :, site_type) + this%lsham(:, :, site_type)
            case default
               call g_logger%fatal('Unsupported operator export: '//trim(name), __FILE__, __LINE__)
            end select
         else
            select case (name)
            case ('hamiltonian')
               block_out(:, :) = this%ee(:, :, neighbor_idx, site_type)
               if (neighbor_idx == 1) block_out(:, :) = block_out(:, :) + this%lsham(:, :, site_type)
            case ('hamiltonian_core')
               block_out(:, :) = this%ee(:, :, neighbor_idx, site_type)
            case ('hamiltonian_overlap')
               block_out(:, :) = this%eeo(:, :, neighbor_idx, site_type)
            case ('hamiltonian_shift')
               if (neighbor_idx == 1) block_out(:, :) = this%enim(:, :, site_type) + this%lsham(:, :, site_type)
            case ('hamiltonian_hoh')
               block_out(:, :) = this%ee(:, :, neighbor_idx, site_type) - this%eeo(:, :, neighbor_idx, site_type)
               if (neighbor_idx == 1) block_out(:, :) = block_out(:, :) + this%enim(:, :, site_type) + this%lsham(:, :, site_type)
            case default
               call g_logger%fatal('Unsupported operator export: '//trim(name), __FILE__, __LINE__)
            end select
         end if
      end subroutine build_operator_block

   end subroutine export_block_operator

   subroutine export_velocity_operator(this, operator_name, v_op, op)
      class(hamiltonian), intent(inout) :: this
      character(len=*), intent(in) :: operator_name
      complex(rp), intent(in) :: v_op(:, :, :, :)
      type(block_sparse_operator), intent(inout) :: op
      integer :: site, nr, m, nnzb, entry, neighbor, site_type
      complex(rp) :: block(nb, nb)

      call op%clear()
      op%block_dim = nb
      op%n_sites = this%lattice%kk
      op%kind = lower(trim(operator_name))
      op%includes_overlap = index(op%kind, 'vo_') == 1
      op%includes_soc = .false.
      op%includes_enim = .false.

      allocate(op%site_types(this%lattice%kk))
      op%site_types(:) = this%lattice%iz(1:this%lattice%kk)

      nnzb = 0
      do site = 1, this%lattice%kk
         nr = max(1, this%lattice%nn(site, 1))
         nnzb = nnzb + 1
         do m = 2, nr
            neighbor = this%lattice%nn(site, m)
            if (neighbor /= 0) nnzb = nnzb + 1
         end do
      end do

      op%nnzb = nnzb
      allocate(op%row_ptr(this%lattice%kk + 1))
      allocate(op%col_ind(nnzb))
      allocate(op%blocks(nb, nb, nnzb))

      entry = 1
      op%row_ptr(1) = 1
      do site = 1, this%lattice%kk
         nr = max(1, this%lattice%nn(site, 1))
         site_type = this%lattice%iz(site)

         block(:, :) = (0.0_rp, 0.0_rp)
         if (site > this%lattice%nmax .or. this%lattice%nmax == 0) block(:, :) = v_op(:, :, 1, site_type)
         op%col_ind(entry) = site
         op%blocks(:, :, entry) = block(:, :)
         entry = entry + 1

         do m = 2, nr
            neighbor = this%lattice%nn(site, m)
            if (neighbor == 0) cycle
            block(:, :) = (0.0_rp, 0.0_rp)
            if (site > this%lattice%nmax .or. this%lattice%nmax == 0) block(:, :) = v_op(:, :, m, site_type)
            op%col_ind(entry) = neighbor
            op%blocks(:, :, entry) = block(:, :)
            entry = entry + 1
         end do
         op%row_ptr(site + 1) = entry
      end do
   end subroutine export_velocity_operator


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
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz

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
         L_op(1:norb, 1:norb) = mLx(:, :)
         L_op(norb+1:nb, norb+1:nb) = mLx(:, :)
      case ('y')
         L_op(1:norb, 1:norb) = mLy(:, :)
         L_op(norb+1:nb, norb+1:nb) = mLy(:, :)
      case ('z')
         L_op(1:norb, 1:norb) = mLz(:, :)
         L_op(norb+1:nb, norb+1:nb) = mLz(:, :)
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
      complex(rp), dimension(nb, nb) :: locham

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

            !if (m==1) then
            !  locham(:,:) = this%ee(:, :, m, ntype) + this%lsham(:, :, ntype) 
            !else
            !  locham(:,:) = this%ee(:, :, m, ntype)
            !end if

            locham(:,:) = this%hxc(:, :, m, ntype)

            tmp1 = 0.0d0; tmp2 = 0.0d0
            ! tmp1 = js_a * hxc(:,:,m,ntype)
            tmp1 = matmul(S_op, locham(:, :))

            ! tmp2 = hxc(:,:,m,ntype) * js_a
            tmp2 = matmul(locham(:, :), S_op)

            ! v_sza(:,:,m,ntype) = 0.5 * ( tmp1 - tmp2 )
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
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz
      complex(rp), dimension(nb, nb) :: locham

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
         L_op(1:norb, 1:norb) = mLx(:, :)
         L_op(norb+1:nb, norb+1:nb) = mLx(:, :)
      case ('y')
         L_op(1:norb, 1:norb) = mLy(:, :)
         L_op(norb+1:nb, norb+1:nb) = mLy(:, :)
      case ('z')
         L_op(1:norb, 1:norb) = mLz(:, :)
         L_op(norb+1:nb, norb+1:nb) = mLz(:, :)
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
                  call zgemm('n', 'n', nb, nb, nb, cone, this%v_a(:, :, m, ntype), nb, this%obarm(:, :, ji), nb, czero, this%vo_a(:, :, m, ntype), nb)
                  call zgemm('n', 'n', nb, nb, nb, cone, this%v_b(:, :, m, ntype), nb, this%obarm(:, :, ji), nb, czero, this%vo_b(:, :, m, ntype), nb)
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
      complex(rp), dimension(norb, norb) :: Lx, Ly, Lz
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
         ! For f-orbitals, use a scaling based on d if not explicitly available
         ! Set soc_f to small value (f-orbital spin-orbit is typically weak)
         ! real soc_f = 0.0_rp  ! f-orbital s-o coupling (for future enhancement)

         ! Check if orbital polarization is enabled
         if (this%orb_pol) then
            rac = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1)*this%charge%lattice%symbolic_atoms(k)%potential%rac)
            lz_loc = sqrt(this%charge%lattice%symbolic_atoms(k)%potential%xi_d(1)*this%charge%lattice%symbolic_atoms(k)%potential%lmom(3))
         else
            rac = 0.0_rp
            lz_loc = 0.0_rp
         end if

         prefac = 0.0_rp
         do i = 1, norb
            do j = 1, norb
               ! p-orbitals (indices 2-4)
               if (i >= 2 .and. i <= 4 .and. j >= 2 .and. j <= 4) prefac = sg*soc_p
               ! d-orbitals (indices 5-9)
               if (i >= 5 .and. i <= 9 .and. j >= 5 .and. j <= 9) prefac = sg*soc_d
               ! f-orbitals (indices 10-16) - currently set to zero (no f-orbital s.o. coupling yet)
               if (i >= 10 .and. i <= 16 .and. j >= 10 .and. j <= 16) prefac = cmplx(0.0_rp, 0.0_rp)

               this%lsham(j, i, k) = this%lsham(j, i, k) + prefac*Lz(j, i) + Lz(j, i)*rac(1)*lz_loc ! H11
               this%lsham(j, i +spin_off, k) = this%lsham(j, i +spin_off, k) + prefac*(Lx(j, i) - i_unit*Ly(j, i)) ! H12
               this%lsham(j +spin_off, i, k) = this%lsham(j +spin_off, i, k) + prefac*(Lx(j, i) + i_unit*Ly(j, i)) ! H21
               this%lsham(j +spin_off, i +spin_off, k) = this%lsham(j +spin_off, i +spin_off, k) - prefac*Lz(j, i) - Lz(j, i)*rac(2)*lz_loc ! H22
            end do
         end do

         ! Debug output: Print on-site Hamiltonian for lmax=3
         if (norb == 16) then
            open(unit=999, file='debug_hamiltonian_lsham.txt', action='write', status='replace')
            write(999, '(A)') 'On-site Hamiltonian (lsham) for lmax=3 (SPDF basis)'
            write(999, '(A, I0)') 'Atom type: ', k
            write(999, '(A)') 'Real part:'
            do i = 1, norb+spin_off
               write(999, '(16F12.6)') (real(this%lsham(i, j, k)), j=1, norb+spin_off)
            end do
            write(999, '(A)') 'Imaginary part:'
            do i = 1, norb+spin_off
               write(999, '(16F12.6)') (aimag(this%lsham(i, j, k)), j=1, norb+spin_off)
            end do
            close(999)
         end if
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
      complex(rp), dimension(norb, norb) :: Lx, Ly, Lz
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
         do i = 1, norb
            do j = 1, norb
               if (i >= 2 .and. i <= 4 .and. j >= 2 .and. j <= 4) prefac = sg*soc_p
               if (i >= 5 .and. i <= 9 .and. j >= 5 .and. j <= 9) prefac = sg*soc_d
               ! build Tx
               this%tmat(j, i, 1, k) = this%tmat(j, i, 1, k) + prefac*i_unit*Ly(j, i)*2.0_rp ! Tx_11
               this%tmat(j, i +spin_off, 1, k) = this%tmat(j, i +spin_off, 1, k) - prefac*Lz(j, i)*2.0_rp*cone ! Tx_12
               this%tmat(j +spin_off, i, 1, k) = this%tmat(j +spin_off, i, 1, k) + prefac*Lz(j, i)*2.0_rp*cone ! Tx_21
               this%tmat(j +spin_off, i +spin_off, 1, k) = this%tmat(j +spin_off, i +spin_off, 1, k) - prefac*i_unit*Ly(j, i)*2.0_rp ! Tx_22
               ! build Ty
               this%tmat(j, i, 2, k) = this%tmat(j, i, 2, k) - prefac*i_unit*Lx(j, i)*2.0_rp ! Ty_11
               this%tmat(j, i +spin_off, 2, k) = this%tmat(j, i +spin_off, 2, k) + prefac*i_unit*Lz(j, i)*2.0_rp ! Ty_12
               this%tmat(j +spin_off, i, 2, k) = this%tmat(j +spin_off, i, 2, k) + prefac*i_unit*Lz(j, i)*2.0_rp ! Ty_21
               this%tmat(j +spin_off, i +spin_off, 2, k) = this%tmat(j +spin_off, i +spin_off, 2, k) + prefac*i_unit*Lx(j, i)*2.0_rp ! Ty_22
               ! build Tz
               this%tmat(j, i +spin_off, 3, k) = this%tmat(j, i +spin_off, 3, k) + prefac*(Lx(j, i) - i_unit*Ly(j, i))*2.0_rp*cone ! Tz_12
               this%tmat(j +spin_off, i, 3, k) = this%tmat(j +spin_off, i, 3, k) + prefac*(Lx(j, i) + i_unit*Ly(j, i))*2.0_rp*cmone ! Tz_21
            end do
         end do
      end do

   end subroutine torque_operator_collinear

   subroutine build_obarm(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      complex(rp), dimension(norb, norb) :: obm0, obm1
      complex(rp), dimension(3) :: mom
      integer :: ntype ! Atom type index
      integer :: l, m ! Orbital index

      this%obarm = 0.d00

      do ntype = 1, this%lattice%ntype
         obm0 = cmplx(0.0d0); obm1 = cmplx(0.0d0)
         do m = 1, norb
            obm0(m, m) = this%lattice%symbolic_atoms(ntype)%potential%obx0(m)
            obm1(m, m) = this%lattice%symbolic_atoms(ntype)%potential%obx1(m)
         end do
         mom(:) = cmplx(this%lattice%symbolic_atoms(ntype)%potential%mom(:), 0.0d0)
         do m = 1, norb
            do l = 1, norb
               this%obarm(m, l, ntype) = obm0(m, l) + obm1(m, l)*mom(3)
               this%obarm(m +spin_off, l +spin_off, ntype) = obm0(m, l) - obm1(m, l)*mom(3)
               this%obarm(l, m +spin_off, ntype) = obm1(m, l)*mom(1) - i_unit*obm1(m, l)*mom(2)
               this%obarm(l +spin_off, m, ntype) = obm1(m, l)*mom(1) + i_unit*obm1(m, l)*mom(2)
            end do
         end do
         call hcpx(this%obarm(1:norb, 1:norb, ntype), 'cart2sph')
         call hcpx(this%obarm(norb+1:nb, norb+1:nb, ntype), 'cart2sph')
         call hcpx(this%obarm(1:norb, norb+1:nb, ntype), 'cart2sph')
         call hcpx(this%obarm(norb+1:nb, 1:norb, ntype), 'cart2sph')
      end do
   end subroutine build_obarm

   subroutine build_enim(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      complex(rp), dimension(norb, norb) :: em0, em1
      complex(rp), dimension(norb) :: ex0, ex1
      complex(rp), dimension(3) :: mom
      complex(rp) :: eu, ed
      integer :: ntype ! Atom type index
      integer :: l, m ! Orbital index

      this%enim = 0.0d0

      do ntype = 1, this%lattice%ntype
         em0 = cmplx(0.0d0); em1 = cmplx(0.0d0)
         do m = 1, norb
            eu = this%lattice%symbolic_atoms(ntype)%potential%cx(m, 1) - this%lattice%symbolic_atoms(ntype)%potential%cex(m, 1)
            ed = this%lattice%symbolic_atoms(ntype)%potential%cx(m, 2) - this%lattice%symbolic_atoms(ntype)%potential%cex(m, 2)
            ex0(m) = 0.5*(eu + ed)
            ex1(m) = 0.5*(eu - ed)
            em0(m, m) = ex0(m)
            em1(m, m) = ex1(m)
         end do
         mom(:) = cmplx(this%lattice%symbolic_atoms(ntype)%potential%mom(:), 0.0d0)
         do m = 1, norb
            do l = 1, norb
               this%enim(m, l, ntype) = em0(m, l) + em1(m, l)*mom(3)
               this%enim(m +spin_off, l +spin_off, ntype) = em0(m, l) - em1(m, l)*mom(3)
               this%enim(l, m +spin_off, ntype) = em1(m, l)*mom(1) - i_unit*em1(m, l)*mom(2)
               this%enim(l +spin_off, m, ntype) = em1(m, l)*mom(1) + i_unit*em1(m, l)*mom(2)
            end do
         end do
         call hcpx(this%enim(1:norb, 1:norb, ntype), 'cart2sph')
         call hcpx(this%enim(norb+1:nb, norb+1:nb, ntype), 'cart2sph')
         call hcpx(this%enim(1:norb, norb+1:nb, ntype), 'cart2sph')
         call hcpx(this%enim(norb+1:nb, 1:norb, ntype), 'cart2sph')

         if (this%local_axis) then
            this%enim_glob = this%enim
         end if
      end do
   end subroutine build_enim

   subroutine build_bulkham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia
      integer :: ntype

      if (this%hubbard_u_general_check) then
         if (rank == 0) call g_logger%info('HUBBARD applying on-site +U correction to bulk Hamiltonian', __FILE__, __LINE__)
         call this%calculate_hubbard_u_potential_general()
      end if
      if (this%hubbard_v_check) then
         if (rank == 0) call g_logger%info('HUBBARD applying inter-site +V correction to bulk Hamiltonian', __FILE__, __LINE__)
         call this%calculate_hubbard_v_potential()
      end if

      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
         ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
         nr = this%charge%lattice%nn(ia, 1) ! Number of neighbours considered
         !write(123, *)´bulkham´
         ! call g_logger%info('Building Hamiltonian for atom type '//fmt('i5', ntype)//' with '//fmt('i5', nr)//' neighbours', __FILE__, __LINE__)
         call this%chbar_nc(ia, nr, ino, ntype)
         do m = 1, nr
            do i = 1, norb
               do j = 1, norb
                  this%ee(j, i, m, ntype) = this%hmag(j, i, m, 4) + this%hmag(j, i, m, 3)        ! H0+Hz
                  this%ee(j +spin_off, i +spin_off, m, ntype) = this%hmag(j, i, m, 4) - this%hmag(j, i, m, 3)        ! H0-Hz
                  this%ee(j, i +spin_off, m, ntype) = this%hmag(j, i, m, 1) - i_unit*this%hmag(j, i, m, 2) ! Hx-iHy
                  this%ee(j +spin_off, i, m, ntype) = this%hmag(j, i, m, 1) + i_unit*this%hmag(j, i, m, 2) ! Hx+iHy
                  ! Builds the magnetic part of the Hamiltonian only
                  this%hxc(j, i, m, ntype) = this%hmag(j, i, m, 3)        ! Hz
                  this%hxc(j +spin_off, i +spin_off, m, ntype) = - this%hmag(j, i, m, 3)        ! - Hz
                  this%hxc(j, i +spin_off, m, ntype) = this%hmag(j, i, m, 1) - i_unit*this%hmag(j, i, m, 2) ! Hx-iHy
                  this%hxc(j +spin_off, i, m, ntype) = this%hmag(j, i, m, 1) + i_unit*this%hmag(j, i, m, 2) ! Hx+iHy
               end do ! end of orbital j loop
            end do ! end of orbital i loop
            write(131,*) 'm=', m
            write(131,'(18f10.6)') real(this%ee(:,:,m,ntype))
            write(132,*) 'm=', m
            write(132,'(18f10.6)') aimag(this%ee(:,:,m,ntype))
         end do ! end of neighbour number
         if (this%hubbard_u_general_check) then
            do i = 1, nb
               do j = 1, nb
                  this%ee(i, j, 1, ntype) = this%ee(i, j, 1, ntype) + cmplx(this%hubbard_u_pot(i, j, ntype), 0.0_rp, kind=rp)
               end do
            end do
         end if
         if (this%hubbard_v_check) then
            do m = 1, nr
               do i = 1, nb
                  do j = 1, nb
                     this%ee(i, j, m, ntype) = this%ee(i, j, m, ntype) + cmplx(this%hubbard_v_pot(i, j, m, ntype), 0.0_rp, kind=rp)
                  end do
               end do
            end do
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
                  call zgemm('n', 'n', nb, nb, nb, cone, this%ee(:, :, m, ntype), nb, this%obarm(:, :, ji), nb, czero, this%eeo(:, :, m, ntype), nb)
                  call zgemm('n', 'c', nb, nb, nb, cone, this%eeo(:, :, m, ntype), nb, this%ee(:, :, m, ntype), nb, czero, this%eeoee(:, :, m, ntype), nb)
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
      if (trim(this%control%recur) == 'chebyshev') then
         call this%compute_hamiltonian_bounds(verbose=.false.)
      end if
      close (128)
   end subroutine build_bulkham

   subroutine build_locham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: it, ino, nr, nlim, m, i, j, ja, ji

      if (this%hubbard_u_general_check) then
         if (rank == 0) call g_logger%info('HUBBARD applying on-site +U correction to local Hamiltonian', __FILE__, __LINE__)
         call this%calculate_hubbard_u_potential_general()
      end if
      if (this%hubbard_v_check) then
         if (rank == 0) call g_logger%info('HUBBARD applying inter-site +V correction to local Hamiltonian', __FILE__, __LINE__)
         call this%calculate_hubbard_v_potential()
      end if

      call g_timer%start('build local hamiltonian')
    !!$omp parallel do private(nlim, nr, ino, m, i, j, ji, ja, this)
      do nlim = 1, this%charge%lattice%nmax
         nr = this%charge%lattice%nn(nlim, 1) ! Number of neighbours considered
         ino = this%charge%lattice%num(nlim)
         call this%chbar_nc(nlim, nr, ino, nlim)
         do m = 1, nr
            do i = 1, norb
               do j = 1, norb
                  this%hall(j, i, m, nlim) = this%hmag(j, i, m, 4) + this%hmag(j, i, m, 3) ! H0+Hz
                  this%hall(j +spin_off, i +spin_off, m, nlim) = this%hmag(j, i, m, 4) - this%hmag(j, i, m, 3) ! H0-Hz
                  this%hall(j, i +spin_off, m, nlim) = this%hmag(j, i, m, 1) - i_unit*this%hmag(j, i, m, 2) ! Hx-iHy
                  this%hall(j +spin_off, i, m, nlim) = this%hmag(j, i, m, 1) + i_unit*this%hmag(j, i, m, 2) ! Hx+iHy
               end do
            end do
         end do
         if (this%hubbard_u_general_check) then
            do i = 1, nb
               do j = 1, nb
                  this%hall(i, j, 1, nlim) = this%hall(i, j, 1, nlim) + cmplx(this%hubbard_u_pot(i, j, ino), 0.0_rp, kind=rp)
               end do
            end do
         end if
         if (this%hubbard_v_check) then
            do m = 1, nr
               do i = 1, nb
                  do j = 1, nb
                     this%hall(i, j, m, nlim) = this%hall(i, j, m, nlim) + cmplx(this%hubbard_v_pot(i, j, m, ino), 0.0_rp, kind=rp)
                  end do
               end do
            end do
         end if
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
                  call zgemm('n', 'n', nb, nb, nb, cone, this%hall(1, 1, m, nlim), nb, this%obarm(1, 1, ji), nb, czero, this%hallo(1, 1, m, nlim), nb)
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
      complex(rp), dimension(nb, nb) :: dum
      n_atoms = this%charge%lattice%ntype
      max_orbital = norb

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

               call hcpx(this%ee(1:norb,1:norb,k,ntype), 'sph2cart')
               call hcpx(this%ee(1:norb,norb+1:nb,k,ntype), 'sph2cart')
               call hcpx(this%ee(norb+1:nb,1:norb,k,ntype), 'sph2cart')
               call hcpx(this%ee(norb+1:nb,norb+1:nb,k,ntype), 'sph2cart')
               do i = 1, nb
                  do j = 1, nb
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
      max_orbital = norb
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
            write (128, '(18f10.6)') real(this%EE(1:nb, 1:nb, k, ntype))*13.605703976
            write (129, '(18f10.6)') aimag(this%EE(1:nb, 1:nb, k, ntype))*13.605703976
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
      max_orbital = norb
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
!                this%ee(orbl+spin_off-(iia-1)*9, orbm+spin_off-(jja-1)*9,k,ntype) = cmplx(dumre,dumcmplx)/13.605703976
!              end if
!            end if
!          end do
            end if
            write (128, *) 'm=', k, 'Atom=', jj, 'Coordinates=', this%charge%lattice%cr(:, jj), 'Ntype=', ntype, 'Index=', idx(k, :)
            write (128, '(18f10.6)') real(this%EE(1:nb, 1:nb, k, ntype))!*13.605703976
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

   subroutine ham0m_nc(this, ia, ja, it, jt, vet, hhh)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: ia, ja ! Atom sites i and j
      integer, intent(in) :: it, jt ! Type of atom i and j
      real(rp), dimension(3), intent(in) :: vet
      real(rp), dimension(norb, norb), intent(in) :: hhh
      ! Local Variables
      integer :: i, j, ilm, jlm, m
      real(rp), dimension(3) :: mom_ia, mom_ja
      real(rp), dimension(3) :: vet_dir
      complex(rp), dimension(3) :: cross
      complex(rp), dimension(norb, norb) :: hhhc
      complex(rp), dimension(this%charge%lattice%ntype, 3) :: momc
      complex(rp) :: dot
      real(rp) :: vv
      real(rp), dimension(3, 3) :: cell
      complex(rp), dimension(nb, nb) :: h_spin
      real(rp), parameter :: qtol = 1.0e-10_rp

      this%hhmag(:, :, :) = 0.0d0

      vv = norm2(vet)
      mom_ia = this%charge%lattice%symbolic_atoms(it)%potential%mom(:)
      mom_ja = this%charge%lattice%symbolic_atoms(jt)%potential%mom(:)
      if (norm2(this%q_ss) > 1.0e-5_rp .or. abs(sin(this%theta_ss)) > 1.0e-8_rp) then
         mom_ia = [sin(this%theta_ss), 0.0_rp, cos(this%theta_ss)]
         mom_ja = mom_ia
      end if

      ! Real to complex
      dot = cmplx(dot_product(mom_ia, mom_ja), kind=kind(0.0d0))
      do i = 1, this%charge%lattice%ntype
         do j = 1, 3
            momc(i, j) = cmplx(this%charge%lattice%symbolic_atoms(i)%potential%mom(j), kind=kind(0.0d0))
         end do
      end do
      momc(it, :) = cmplx(mom_ia, kind=kind(0.0d0))
      momc(jt, :) = cmplx(mom_ja, kind=kind(0.0d0))
      cross = cmplx(cross_product(mom_ia, mom_ja), kind=kind(0.0d0))
      hhhc(:, :) = cmplx(hhh(:, :), kind=kind(0.0d0))

      if (size(this%charge%lattice%symbolic_atoms(it)%potential%wx0) < norb) then
         call g_logger%fatal('wx0(it) too small: type='//int2str(it)// &
                             ' size='//int2str(size(this%charge%lattice%symbolic_atoms(it)%potential%wx0))// &
                             ' norb='//int2str(norb)// &
                             ' lmax='//int2str(this%charge%lattice%symbolic_atoms(it)%potential%lmax), __FILE__, __LINE__)
      end if
      if (size(this%charge%lattice%symbolic_atoms(jt)%potential%wx0) < norb) then
         call g_logger%fatal('wx0(jt) too small: type='//int2str(jt)// &
                             ' size='//int2str(size(this%charge%lattice%symbolic_atoms(jt)%potential%wx0))// &
                             ' norb='//int2str(norb)// &
                             ' lmax='//int2str(this%charge%lattice%symbolic_atoms(jt)%potential%lmax), __FILE__, __LINE__)
      end if

      do ilm = 1, norb
         do jlm = 1, norb
            this%hhmag(ilm, jlm, 4) = &
               this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm) + &
               this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*dot
         end do
      end do

!    do ilm = 1, norb
!      write(123, ´(9f10.6)´) (real(this%hhmag(ilm, jlm, 4)), jlm=1, 9)
!    end do

      if (vv <= 0.01d0) then
         do ilm = 1, norb
            if (this%hoh) then
               this%hhmag(ilm, ilm, 4) = this%hhmag(ilm, ilm, 4) + this%charge%lattice%symbolic_atoms(it)%potential%cex0(ilm)
            else
               this%hhmag(ilm, ilm, 4) = this%hhmag(ilm, ilm, 4) + this%charge%lattice%symbolic_atoms(it)%potential%cx0(ilm)
            end if
         end do
      end if

      do m = 1, 3
         do jlm = 1, norb
            do ilm = 1, norb
               this%hhmag(ilm, jlm, m) = &
                  (this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm))*momc(it, m) + &
                  (this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm))*momc(jt, m) + &
                  i_unit*this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm, jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*cross(m)
            end do
         end do
      end do

      if (norm2(this%q_ss) > qtol .and. vv > 0.01d0) then
         call build_spinor_from_hhmag(this%hhmag, h_spin)
         cell(:, :) = this%charge%lattice%a(:, :)*this%charge%lattice%alat
         vet_dir(:) = cartesian_to_direct(cell, vet)
         call rotate_spin_spiral_bond(this%q_ss, vet_dir, h_spin)
         call project_spinor_to_hhmag(h_spin, this%hhmag)
      end if

      if (vv > 0.01d0) return
      do m = 1, 3
         do ilm = 1, norb
            if (this%hoh) then
               this%hhmag(ilm, ilm, m) = this%hhmag(ilm, ilm, m) + this%charge%lattice%symbolic_atoms(it)%potential%cex1(ilm)*momc(it, m)
            else
               this%hhmag(ilm, ilm, m) = this%hhmag(ilm, ilm, m) + this%charge%lattice%symbolic_atoms(it)%potential%cx1(ilm)*momc(it, m)
            end if
         end do
      end do

      !do m=1, 3
      !  write(123, *)´m=´, m
      !  do ilm = 1, norb
      !    write(123, ´(9f10.6)´) (real(this%hhmag(ilm, jlm, m)), jlm=1, 9)
      !  end do
      !end do
   end subroutine ham0m_nc

   subroutine build_spinor_from_hhmag(hhmag_in, h_spin)
      complex(rp), intent(in) :: hhmag_in(norb, norb, 4)
      complex(rp), intent(out) :: h_spin(nb, nb)
      integer :: iorb, jorb

      h_spin(:, :) = cmplx(0.0_rp, 0.0_rp, rp)
      do iorb = 1, norb
         do jorb = 1, norb
            h_spin(iorb, jorb) = hhmag_in(iorb, jorb, 4) + hhmag_in(iorb, jorb, 3)
            h_spin(iorb + spin_off, jorb + spin_off) = hhmag_in(iorb, jorb, 4) - hhmag_in(iorb, jorb, 3)
            h_spin(iorb, jorb + spin_off) = hhmag_in(iorb, jorb, 1) - i_unit*hhmag_in(iorb, jorb, 2)
            h_spin(iorb + spin_off, jorb) = hhmag_in(iorb, jorb, 1) + i_unit*hhmag_in(iorb, jorb, 2)
         end do
      end do
   end subroutine build_spinor_from_hhmag

   subroutine project_spinor_to_hhmag(h_spin, hhmag_out)
      complex(rp), intent(in) :: h_spin(nb, nb)
      complex(rp), intent(out) :: hhmag_out(norb, norb, 4)
      integer :: iorb, jorb
      complex(rp) :: huu, hdd, hud, hdu

      hhmag_out(:, :, :) = cmplx(0.0_rp, 0.0_rp, rp)
      do iorb = 1, norb
         do jorb = 1, norb
            huu = h_spin(iorb, jorb)
            hdd = h_spin(iorb + spin_off, jorb + spin_off)
            hud = h_spin(iorb, jorb + spin_off)
            hdu = h_spin(iorb + spin_off, jorb)
            hhmag_out(iorb, jorb, 4) = 0.5_rp*(huu + hdd)
            hhmag_out(iorb, jorb, 3) = 0.5_rp*(huu - hdd)
            hhmag_out(iorb, jorb, 1) = 0.5_rp*(hud + hdu)
            hhmag_out(iorb, jorb, 2) = 0.5_rp*i_unit*(hud - hdu)
         end do
      end do
   end subroutine project_spinor_to_hhmag

   subroutine rotate_spin_spiral_bond(q_vec, r_vec_direct, h_spin)
      real(rp), intent(in) :: q_vec(3), r_vec_direct(3)
      complex(rp), intent(inout) :: h_spin(nb, nb)
      integer :: idx
      real(rp) :: phi
      complex(rp) :: phase_up, phase_dn
      complex(rp), dimension(nb) :: left_phase

      phi = 2.0_rp*pi*dot_product(q_vec, r_vec_direct)
      phase_up = exp(-0.5_rp*i_unit*phi)
      phase_dn = exp(0.5_rp*i_unit*phi)

      left_phase(1:spin_off) = phase_up
      left_phase(spin_off + 1:nb) = phase_dn

      do idx = 1, nb
         h_spin(idx, :) = left_phase(idx)*h_spin(idx, :)
      end do
      do idx = 1, nb
         h_spin(:, idx) = h_spin(:, idx)*conjg(left_phase(idx))
      end do
   end subroutine rotate_spin_spiral_bond

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
      real(rp), dimension(norb, norb) :: hhh
      integer :: i, j, k, l, m, n, it, jt, jj, nn_max_loc
      integer :: ni, mdir
      integer :: kk ! Clust size number
      real(rp), dimension(:, :), allocatable :: ham_vec

      this%hmag(:, :, :, :) = 0.0d0

      r2 = this%charge%lattice%r2
      cralat(:, :) = this%charge%lattice%cr(:, :)*this%charge%lattice%alat
      kk = this%charge%lattice%kk

      allocate(ham_vec(3, nr))
      nn_max_loc = nr
      call this%charge%lattice%clusba(r2, cralat, ia, kk, kk, nn_max_loc, ham_vec)

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
            if (this%lattice%pbc) then
               call this%lattice%f_wrap_coord_diff(this%lattice%kk,this%lattice%cr*this%lattice%alat,ia,jj,vet)
            else
               vet(:) = (this%charge%lattice%cr(:, jj) - this%charge%lattice%cr(:, ia))*this%charge%lattice%alat
            end if
            !write(123, ´(3f10.6)´) vet(:)
            !write(123, ´(3f10.6)´) this%charge%lattice%sbarvec(:, m)
            !write(123, ´(a, 3i4, 3f10.6)´) ´nn ´, IA, m, JJ, VET(:)
            call this%hmfind(vet, nr, hhh, m, ia, m, ni, ham_vec)
            if (ni == 0) then
               this%charge%lattice%nn(ia, m) = 0
            end if
            call this%ham0m_nc(ia, jj, it, jt, vet, hhh)
            do mdir = 1, 4
               call hcpx(this%hhmag(:, :, mdir), 'cart2sph')
               this%hmag(:, :, m, mdir) = this%hhmag(:, :, mdir)
            end do
         end if
      end do
      if (allocated(ham_vec)) deallocate(ham_vec)
      !do m=1, nr
      !  write(123, *)´m=´, m
      !  do mdir=1, 4
      !    write(123, *)´mdir=´, mdir
      !    do i = 1, norb
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
      real(rp), dimension(norb, norb), intent(inout) :: hhh
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
         !write(123, ´(a, i4, 3f10.4)´)´i´, i, this%charge%lattice%sbarvec(:, i)
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
      do ilm = 1, norb
         do jlm = 1, norb
            hhh(ilm, jlm) = real(this%charge%lattice%sbar(jlm, ilm, m, this%charge%lattice%num(ia)))
         end do
      end do

      !do ilm = 1, norb
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
         i_out = modulo(orb - 1, max_orbital) + max_orbital + 1
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
         sdim = product(shape(this%hall))/nb/nb
         call rotmag_loc(this%hall, this%hall_glob, sdim, m_loc)
         sdim = product(shape(this%ee))/nb/nb
         call rotmag_loc(this%ee, this%ee_glob, sdim, m_loc)
         if (this%hoh) then
            sdim = product(shape(this%eeo))/nb/nb
            call rotmag_loc(this%eeo, this%eeo_glob, sdim, m_loc)
            sdim = product(shape(this%hallo))/nb/nb
            call rotmag_loc(this%hallo, this%hallo_glob, sdim, m_loc)
            sdim = product(shape(this%enim))/nb/nb
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

   subroutine calculate_hubbard_u_potential_general(this)
      class(hamiltonian), intent(inout) :: this

      integer :: na, l, ispin, i, j, m1, m2, m3, m4, l_index, m_max
      integer :: m1_val, m2_val, m3_val, m4_val
      real(rp) :: f0, f2, f4, f6
      real(rp) :: ubar, jbar, ueff, d1, d2, eps_den, tr_n1mn
      real(rp) :: num_u, num_j, sum_occ_opposite, sum_occ_same_excl
      real(rp) :: common_pref, sum_u_aux, sum_j_aux, dUdn, dJdn, dUeff_dn
      real(rp) :: vdiag_up_avg, vdiag_dn_avg, vdiag_split
      real(rp), dimension(2, 7) :: occ_m
      real(rp), dimension(2, 7, 7) :: pbar
      logical :: use_acbn0
      real(rp), dimension(this%lattice%ntype, 4) :: hub_u, hub_j
      real(rp), dimension(this%lattice%ntype, 4, 4) :: f
      real(rp), dimension(this%lattice%ntype, 4, 2, 7, 7) :: ldm, hub_pot
      real(rp), dimension(this%lattice%ntype, 4, 2) :: n_spin
      real(rp), dimension(this%lattice%ntype, 4) :: n_tot
      type :: int_array
         integer, allocatable :: val(:)
      end type int_array
      type(int_array), dimension(this%lattice%ntype) :: l_arr
      type(int_array), dimension(4) :: ms
      integer :: cntr
      real(rp) :: ldm_trace_up, ldm_trace_dn, ldm_max_abs, hup_max_abs
      real(rp) :: ql_occ_up, ql_occ_dn

      ms(1)%val = [0]
      ms(2)%val = [-1, 0, 1]
      ms(3)%val = [-2, -1, 0, 1, 2]
      ms(4)%val = [-3, -2, -1, 0, 1, 2, 3]

      this%hubbard_u_pot(:, :, :) = 0.0_rp
      hub_u(:, :) = 0.0_rp
      hub_j(:, :) = 0.0_rp
      f(:, :, :) = 0.0_rp
      ldm(:, :, :, :, :) = 0.0_rp
      hub_pot(:, :, :, :, :) = 0.0_rp
      n_spin(:, :, :) = 0.0_rp
      n_tot(:, :) = 0.0_rp

      do na = 1, this%lattice%ntype
         do l = 1, min(4, size(this%lattice%symbolic_atoms(na)%potential%hubbard_u))
            hub_u(na, l) = this%lattice%symbolic_atoms(na)%potential%hubbard_u(l)
         end do
         do l = 1, min(4, size(this%lattice%symbolic_atoms(na)%potential%hubbard_j))
            hub_j(na, l) = this%lattice%symbolic_atoms(na)%potential%hubbard_j(l)
         end do
      end do

      do na = 1, this%lattice%ntype
         f(na, 1, 1) = hub_u(na, 1)
         f(na, 2, 1) = hub_u(na, 2)
         f(na, 2, 2) = hub_j(na, 2)*5.0_rp
         f(na, 3, 1) = hub_u(na, 3)
         f(na, 3, 2) = 14.0_rp*hub_j(na, 3)/1.625_rp
         f(na, 3, 3) = 0.625_rp*f(na, 3, 2)
         f(na, 4, 1) = hub_u(na, 4)
         f(na, 4, 2) = 6435.0_rp*hub_j(na, 4)/(286.0_rp + 195.0_rp*0.67_rp + 250.0_rp*0.49_rp)
         f(na, 4, 3) = 0.67_rp*f(na, 4, 2)
         f(na, 4, 4) = 0.49_rp*f(na, 4, 2)
      end do

      do na = 1, this%lattice%ntype
         cntr = count((abs(hub_u(na, :)) > 1.0e-10_rp) .or. (abs(hub_j(na, :)) > 1.0e-10_rp))
         allocate(l_arr(na)%val(cntr))
         cntr = 0
         do l = 1, 4
            if ((abs(hub_u(na, l)) > 1.0e-10_rp) .or. (abs(hub_j(na, l)) > 1.0e-10_rp)) then
               cntr = cntr + 1
               l_arr(na)%val(cntr) = l
            end if
         end do
      end do

      do na = 1, this%lattice%ntype
         do l = 0, min(3, this%control%lmax)
            do ispin = 1, 2
               do i = 1, 2*l + 1
                  do j = 1, 2*l + 1
                     ldm(na, l + 1, ispin, i, j) = this%lattice%symbolic_atoms(na)%potential%ldm(l + 1, ispin, i, j)
                  end do
               end do
            end do
         end do
      end do

      do na = 1, this%lattice%ntype
         do l = 1, min(4, this%control%lmax + 1)
            if ((abs(hub_u(na, l)) <= 1.0e-10_rp) .and. (abs(hub_j(na, l)) <= 1.0e-10_rp)) cycle
            ldm_trace_up = 0.0_rp
            ldm_trace_dn = 0.0_rp
            ldm_max_abs = 0.0_rp
            do i = 1, 2*l - 1
               ldm_trace_up = ldm_trace_up + ldm(na, l, 1, i, i)
               ldm_trace_dn = ldm_trace_dn + ldm(na, l, 2, i, i)
               do j = 1, 2*l - 1
                  ldm_max_abs = max(ldm_max_abs, abs(ldm(na, l, 1, i, j)))
                  ldm_max_abs = max(ldm_max_abs, abs(ldm(na, l, 2, i, j)))
               end do
            end do
            if (rank == 0) then
               call g_logger%info('HUBBARD_LDM type='//fmt('i4', na)//' l='//fmt('i2', l - 1)// &
                  ' tr_up='//fmt('f10.6', ldm_trace_up)//' tr_dn='//fmt('f10.6', ldm_trace_dn)// &
                  ' maxabs='//fmt('es12.4', ldm_max_abs), __FILE__, __LINE__)
               call g_logger%info('HUBBARD_OCC type='//fmt('i4', na)//' l='//fmt('i2', l - 1)// &
                  ' U='//fmt('f10.6', hub_u(na, l))//' J='//fmt('f10.6', hub_j(na, l))// &
                  ' nup='//fmt('f10.6', ldm_trace_up)//' ndn='//fmt('f10.6', ldm_trace_dn)// &
                  ' ntot='//fmt('f10.6', ldm_trace_up + ldm_trace_dn), __FILE__, __LINE__)
            end if
            ql_occ_up = 0.0_rp
            ql_occ_dn = 0.0_rp
            if ((l - 1) <= this%lattice%symbolic_atoms(na)%potential%lmax) then
               ql_occ_up = this%lattice%symbolic_atoms(na)%potential%ql(1, l - 1, 1)
               ql_occ_dn = this%lattice%symbolic_atoms(na)%potential%ql(1, l - 1, 2)
            end if
            if (rank == 0) then
               call g_logger%info('HUBBARD_QLCMP type='//fmt('i4', na)//' l='//fmt('i2', l - 1)// &
                  ' ql_up='//fmt('f10.6', ql_occ_up)//' ql_dn='//fmt('f10.6', ql_occ_dn)// &
                  ' ldm_up='//fmt('f10.6', ldm_trace_up)//' ldm_dn='//fmt('f10.6', ldm_trace_dn)// &
                  ' d_up(ql-ldm)='//fmt('f10.6', ql_occ_up - ldm_trace_up)// &
                  ' d_dn(ql-ldm)='//fmt('f10.6', ql_occ_dn - ldm_trace_dn), __FILE__, __LINE__)
            end if
         end do
      end do

      do na = 1, this%lattice%ntype
         do l = 1, 4
            do ispin = 1, 2
               do i = 1, 2*l - 1
                  n_spin(na, l, ispin) = n_spin(na, l, ispin) + ldm(na, l, ispin, i, i)
               end do
            end do
            ! Total occupancy for this channel is n_up + n_dn.
            n_tot(na, l) = n_spin(na, l, 1) + n_spin(na, l, 2)
         end do
      end do

      use_acbn0 = (trim(lower(this%hubbard_u_potential_form)) == 'acbn0')
      eps_den = 1.0e-12_rp

      do na = 1, this%lattice%ntype
         do l = 1, size(l_arr(na)%val)
            l_index = l_arr(na)%val(l)
            m_max = 2*l_index - 1

            if (use_acbn0) then
               occ_m(:, :) = 0.0_rp
               pbar(:, :, :) = 0.0_rp
               do ispin = 1, 2
                  do m1 = 1, m_max
                     occ_m(ispin, m1) = ldm(na, l_index, ispin, m1, m1)
                  end do
               end do
               do ispin = 1, 2
                  do m1 = 1, m_max
                     do m2 = 1, m_max
                        pbar(ispin, m1, m2) = (occ_m(ispin, m1) + occ_m(ispin, m2))*ldm(na, l_index, ispin, m1, m2)
                     end do
                  end do
               end do

               d1 = 0.0_rp
               d2 = 0.0_rp
               do m1 = 1, m_max
                  do m2 = 1, m_max
                     d1 = d1 + occ_m(1, m1)*occ_m(2, m2) + occ_m(2, m1)*occ_m(1, m2)
                     if (m1 /= m2) then
                        d1 = d1 + occ_m(1, m1)*occ_m(1, m2) + occ_m(2, m1)*occ_m(2, m2)
                        d2 = d2 + occ_m(1, m1)*occ_m(1, m2) + occ_m(2, m1)*occ_m(2, m2)
                     end if
                  end do
               end do

               f0 = f(na, l_index, 1)
               f2 = f(na, l_index, 2)
               f4 = f(na, l_index, 3)
               f6 = f(na, l_index, 4)
               num_u = 0.0_rp
               num_j = 0.0_rp
               do ispin = 1, 2
                  do i = 1, 2
                     do m1 = 1, m_max
                        do m2 = 1, m_max
                           do m3 = 1, m_max
                              do m4 = 1, m_max
                                 m1_val = ms(l_index)%val(m1)
                                 m2_val = ms(l_index)%val(m2)
                                 m3_val = ms(l_index)%val(m3)
                                 m4_val = ms(l_index)%val(m4)
                                 num_u = num_u + Coulomb_mat(l_index - 1, m1_val, m3_val, m2_val, m4_val, f0, f2, f4, f6)* &
                                         pbar(ispin, m1, m2)*pbar(i, m3, m4)
                              end do
                           end do
                        end do
                     end do
                  end do
                  do m1 = 1, m_max
                     do m2 = 1, m_max
                        if (m1 == m2) cycle
                        do m3 = 1, m_max
                           do m4 = 1, m_max
                              m1_val = ms(l_index)%val(m1)
                              m2_val = ms(l_index)%val(m2)
                              m3_val = ms(l_index)%val(m3)
                              m4_val = ms(l_index)%val(m4)
                              num_j = num_j + Coulomb_mat(l_index - 1, m1_val, m3_val, m4_val, m2_val, f0, f2, f4, f6)* &
                                      pbar(ispin, m1, m2)*pbar(ispin, m3, m4)
                           end do
                        end do
                     end do
                  end do
               end do
               ubar = num_u/max(d1, eps_den)
               jbar = num_j/max(d2, eps_den)
               ueff = ubar - jbar
               tr_n1mn = 0.0_rp
               do ispin = 1, 2
                  do m1 = 1, m_max
                     tr_n1mn = tr_n1mn + occ_m(ispin, m1)*(1.0_rp - occ_m(ispin, m1))
                  end do
               end do
            end if

            do ispin = 1, 2
               do m1 = 1, m_max
                  do m2 = 1, m_max
                     do m3 = 1, m_max
                        do m4 = 1, m_max
                           m1_val = ms(l_index)%val(m1)
                           m2_val = ms(l_index)%val(m2)
                           m3_val = ms(l_index)%val(m3)
                           m4_val = ms(l_index)%val(m4)
                           f0 = f(na, l_index, 1)
                           f2 = f(na, l_index, 2)
                           f4 = f(na, l_index, 3)
                           f6 = f(na, l_index, 4)
                           hub_pot(na, l_index, ispin, m1, m2) = hub_pot(na, l_index, ispin, m1, m2) &
                                + Coulomb_mat(l_index - 1, m1_val, m3_val, m2_val, m4_val, f0, f2, f4, f6)*ldm(na, l_index, 3 - ispin, m3, m4) &
                                + (Coulomb_mat(l_index - 1, m1_val, m3_val, m2_val, m4_val, f0, f2, f4, f6) &
                                -  Coulomb_mat(l_index - 1, m1_val, m3_val, m4_val, m2_val, f0, f2, f4, f6))*ldm(na, l_index, ispin, m3, m4)
                        end do
                     end do
                     if (use_acbn0) then
                        sum_u_aux = 0.0_rp
                        sum_j_aux = 0.0_rp
                        m1_val = ms(l_index)%val(m1)
                        m2_val = ms(l_index)%val(m2)
                        do j = 1, 2
                           do m3 = 1, m_max
                              do m4 = 1, m_max
                                 m3_val = ms(l_index)%val(m3)
                                 m4_val = ms(l_index)%val(m4)
                                 sum_u_aux = sum_u_aux + Coulomb_mat(l_index - 1, m1_val, m3_val, m2_val, m4_val, f0, f2, f4, f6)*pbar(j, m3, m4)
                              end do
                           end do
                        end do
                        do m3 = 1, m_max
                           do m4 = 1, m_max
                              m3_val = ms(l_index)%val(m3)
                              m4_val = ms(l_index)%val(m4)
                              sum_j_aux = sum_j_aux + Coulomb_mat(l_index - 1, m1_val, m3_val, m4_val, m2_val, f0, f2, f4, f6)*pbar(ispin, m3, m4)
                           end do
                        end do

                        common_pref = occ_m(ispin, m1) + occ_m(ispin, m2)
                        if (m1 == m2) common_pref = common_pref + 2.0_rp*occ_m(ispin, m1)
                        sum_occ_opposite = sum(occ_m(3 - ispin, 1:m_max))
                        sum_occ_same_excl = sum(occ_m(ispin, 1:m_max)) - occ_m(ispin, m1)
                        if (m1 == m2) then
                           dUdn = common_pref*sum_u_aux/max(d1, eps_den) - &
                                  (2.0_rp*ubar/max(d1, eps_den))*(sum_occ_opposite - sum_occ_same_excl)
                           dJdn = common_pref*sum_j_aux/max(d2, eps_den) - &
                                  (2.0_rp*jbar/max(d2, eps_den))*sum_occ_same_excl
                        else
                           dUdn = common_pref*sum_u_aux/max(d1, eps_den)
                           dJdn = common_pref*sum_j_aux/max(d2, eps_den)
                        end if
                        dUeff_dn = dUdn - dJdn
                     end if
                     if (m1 == m2) then
                        if (use_acbn0) then
                           hub_pot(na, l_index, ispin, m1, m2) = 0.0_rp
                           hub_pot(na, l_index, ispin, m1, m2) = hub_pot(na, l_index, ispin, m1, m2) + &
                              ueff*(0.5_rp - ldm(na, l_index, ispin, m1, m2)) + 0.5_rp*dUeff_dn*tr_n1mn
                        else
                           hub_pot(na, l_index, ispin, m1, m2) = hub_pot(na, l_index, ispin, m1, m2) &
                              - hub_u(na, l_index)*(n_tot(na, l_index) - 0.5_rp) &
                              + hub_j(na, l_index)*(n_spin(na, l_index, ispin) - 0.5_rp)
                        end if
                     else if (use_acbn0) then
                        hub_pot(na, l_index, ispin, m1, m2) = -ueff*ldm(na, l_index, ispin, m1, m2) + 0.5_rp*dUeff_dn*tr_n1mn
                     end if
                  end do
               end do
            end do

            ! Diagnostic for expected splitting: average diagonal on-site +U potential per spin.
            vdiag_up_avg = 0.0_rp
            vdiag_dn_avg = 0.0_rp
            do m1 = 1, m_max
               vdiag_up_avg = vdiag_up_avg + hub_pot(na, l_index, 1, m1, m1)
               vdiag_dn_avg = vdiag_dn_avg + hub_pot(na, l_index, 2, m1, m1)
            end do
            vdiag_up_avg = vdiag_up_avg/real(max(1, m_max), rp)
            vdiag_dn_avg = vdiag_dn_avg/real(max(1, m_max), rp)
            vdiag_split = vdiag_up_avg - vdiag_dn_avg
            if (rank == 0) then
               call g_logger%info('HUBBARD_VDIAG type='//fmt('i4', na)//' l='//fmt('i2', l_index - 1)// &
                  ' avg_up='//fmt('f11.6', vdiag_up_avg)//' avg_dn='//fmt('f11.6', vdiag_dn_avg)// &
                  ' split(up-dn)='//fmt('f11.6', vdiag_split), __FILE__, __LINE__)
            end if
         end do
      end do

      do na = 1, this%lattice%ntype
         do l = 0, min(3, this%control%lmax)
            do i = 1, 2*l + 1
               do j = 1, 2*l + 1
                  this%hubbard_u_pot(l**2 + i, l**2 + j, na) = hub_pot(na, l + 1, 1, i, j)
                  this%hubbard_u_pot(l**2 + i + spin_off, l**2 + j + spin_off, na) = hub_pot(na, l + 1, 2, i, j)
               end do
            end do
         end do
      end do

      do na = 1, this%lattice%ntype
         hup_max_abs = 0.0_rp
         do i = 1, nb
            do j = 1, nb
               hup_max_abs = max(hup_max_abs, abs(this%hubbard_u_pot(i, j, na)))
            end do
         end do
         if (rank == 0) call g_logger%info('HUBBARD_POT type='//fmt('i4', na)//' maxabs='//fmt('es12.4', hup_max_abs), __FILE__, __LINE__)
      end do

      do na = 1, this%lattice%ntype
         if (allocated(l_arr(na)%val)) deallocate(l_arr(na)%val)
      end do
   end subroutine calculate_hubbard_u_potential_general

   subroutine calculate_hubbard_v_potential(this)
      class(hamiltonian), intent(inout) :: this

      integer :: itype, ia, nr, m, ja, jt, ispin, li, lj, ii, i0, i1, jdim
      real(rp) :: rmin, rcur, tol, occ_up, occ_dn, nn_shell_tol
      real(rp), dimension(3) :: dr
      real(rp), dimension(this%lattice%ntype, 4, 2) :: n_spin
      logical, save :: warned_proxy_once = .false.

      this%hubbard_v_pot(:, :, :, :) = 0.0_rp
      n_spin(:, :, :) = 0.0_rp
      tol = 1.0e-4_rp

      if (rank == 0 .and. .not. warned_proxy_once) then
         call g_logger%warning('HUBBARD +V currently uses a diagonal local-occupation proxy. Full inter-site n^{JI}_{mm''} from inter-site Green functions is not yet wired into calculate_hubbard_v_potential.', __FILE__, __LINE__)
         warned_proxy_once = .true.
      end if

      ! Spin- and l-resolved occupations from local density matrices.
      do itype = 1, this%lattice%ntype
         do li = 1, min(4, this%control%lmax + 1)
            jdim = 2*li - 1
            do ispin = 1, 2
               do ii = 1, jdim
                  n_spin(itype, li, ispin) = n_spin(itype, li, ispin) + &
                     this%lattice%symbolic_atoms(itype)%potential%ldm(li, ispin, ii, ii)
               end do
            end do
         end do
      end do

      do itype = 1, this%lattice%ntype
         ia = this%lattice%atlist(itype)
         nr = this%lattice%nn(ia, 1)

         ! Find nearest-neighbour shell radius for this atom type.
         rmin = huge(1.0_rp)
         do m = 2, nr
            ja = this%lattice%nn(ia, m)
            if (ja == 0) cycle
            if (this%lattice%pbc) then
               call this%lattice%f_wrap_coord_diff(this%lattice%kk, this%lattice%cr*this%lattice%alat, ia, ja, dr)
            else
               dr(:) = (this%lattice%cr(:, ja) - this%lattice%cr(:, ia))*this%lattice%alat
            end if
            rcur = norm2(dr)
            if (rcur > tol .and. rcur < rmin) rmin = rcur
         end do
         if (rmin >= huge(1.0_rp)*0.5_rp) cycle

         do m = 2, nr
            ja = this%lattice%nn(ia, m)
            if (ja == 0) cycle
            jt = this%lattice%iz(ja)
            if (jt <= 0 .or. jt > this%lattice%ntype) cycle

            if (this%lattice%pbc) then
               call this%lattice%f_wrap_coord_diff(this%lattice%kk, this%lattice%cr*this%lattice%alat, ia, ja, dr)
            else
               dr(:) = (this%lattice%cr(:, ja) - this%lattice%cr(:, ia))*this%lattice%alat
            end if
            rcur = norm2(dr)
            nn_shell_tol = max(5.0e-4_rp, 2.5e-3_rp*rmin)
            if (abs(rcur - rmin) > nn_shell_tol) cycle

            ! Approximate intersite +V contribution (current proxy implementation):
            ! V^I,sigma_mm' += - sum_{J in NN} V^{IJ}_{li,lj} * n^{J,sigma}_{lj} / (2*li-1) * delta_mm'
            !
            ! IMPORTANT LIMITATION:
            ! The full +V expression requires inter-site density matrices n^{JI,sigma}_{mprime,m}
            ! from inter-site Green functions. That machinery is not yet wired here, so only
            ! an orbital-diagonal proxy based on local occupations is applied.
            do li = 1, min(4, this%control%lmax + 1)
               if (li > size(this%hubbard_v, 3)) cycle
               if (li > size(this%hubbard_v, 4)) cycle
               jdim = 2*li - 1
               i0 = (li - 1)*(li - 1) + 1
               i1 = i0 + jdim - 1
               do lj = 1, min(4, size(this%hubbard_v, 4), this%control%lmax + 1)
                  if (lj > size(this%hubbard_v, 3)) cycle
                  if (abs(this%hubbard_v(itype, jt, li, lj)) <= 1.0e-12_rp) cycle
                  occ_up = n_spin(jt, lj, 1)
                  occ_dn = n_spin(jt, lj, 2)
                  do ii = i0, i1
                     this%hubbard_v_pot(ii, ii, m, itype) = this%hubbard_v_pot(ii, ii, m, itype) - &
                        this%hubbard_v(itype, jt, li, lj)*occ_up/max(1, jdim)
                     this%hubbard_v_pot(ii + spin_off, ii + spin_off, m, itype) = &
                        this%hubbard_v_pot(ii + spin_off, ii + spin_off, m, itype) - &
                        this%hubbard_v(itype, jt, li, lj)*occ_dn/max(1, jdim)
                  end do
               end do
            end do
         end do
      end do
   end subroutine calculate_hubbard_v_potential

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Compute bulk Hamiltonian spectrum bounds for Chebyshev scaling.
   !---------------------------------------------------------------------------
   subroutine compute_hamiltonian_bounds(this, verbose)
      class(hamiltonian), intent(inout) :: this
      logical, intent(in), optional :: verbose

      integer :: ntype, ia, nr, i, j, m, n_orb, n_sites
      integer :: isite, jsite, i_start, i_end, j_start, j_end, ineigh, ntype_i, ia_loc, ja
      real(rp) :: g_min, g_max, center, radius
      real(rp) :: hgamma_min, hgamma_max
      logical :: verb, have_gamma
      character(len=16) :: algo
      character(len=256) :: msg
      type(bounds) :: gamma_bounds
      complex(rp), allocatable :: h_gamma(:, :)

      verb = .false.
      if (present(verbose)) verb = verbose

      if (.not. allocated(this%ee)) then
         call g_logger%warning('compute_hamiltonian_bounds: ee is not allocated; skipping bounds', __FILE__, __LINE__)
         return
      end if

      call normalize_bounds_algorithm(this%bounds%algorithm, algo)
      if (algo == 'none') return

      g_min = huge(1.0_rp)
      g_max = -huge(1.0_rp)

      do ntype = 1, this%charge%lattice%ntype
         ia = this%charge%lattice%atlist(ntype)
         nr = this%charge%lattice%nn(ia, 1)
         do i = 1, nb
            center = real(this%ee(i, i, 1, ntype))
            if (allocated(this%lsham)) center = center + real(this%lsham(i, i, ntype))
            radius = 0.0_rp
            do m = 1, nr
               do j = 1, nb
                  if (m == 1 .and. i == j) cycle
                  radius = radius + abs(this%ee(i, j, m, ntype))
               end do
            end do
            g_min = min(g_min, center - radius)
            g_max = max(g_max, center + radius)
         end do
      end do

      have_gamma = .false.
      hgamma_min = g_min
      hgamma_max = g_max
      n_orb = nb
      n_sites = this%charge%lattice%nrec
      if (n_sites > 0) then
         allocate(h_gamma(n_orb*n_sites, n_orb*n_sites))
         h_gamma(:, :) = cmplx(0.0_rp, 0.0_rp, kind=rp)
         do isite = 1, n_sites
            ntype_i = this%charge%lattice%ib(isite)
            ia_loc = this%charge%lattice%atlist(ntype_i)
            i_start = (isite - 1)*n_orb + 1
            i_end = isite*n_orb
            do ineigh = 1, this%charge%lattice%nn(ia_loc, 1)
               if (ineigh == 1) then
                  jsite = isite
               else
                  ja = this%charge%lattice%nn(ia_loc, ineigh)
                  jsite = this%charge%lattice%iz(ja)
                  if (jsite < 1 .or. jsite > n_sites) cycle
               end if
               j_start = (jsite - 1)*n_orb + 1
               j_end = jsite*n_orb
               h_gamma(i_start:i_end, j_start:j_end) = h_gamma(i_start:i_end, j_start:j_end) + this%ee(:, :, ineigh, ntype_i)
            end do
            if (allocated(this%lsham)) then
               h_gamma(i_start:i_end, i_start:i_end) = h_gamma(i_start:i_end, i_start:i_end) + this%lsham(:, :, ntype_i)
            end if
         end do
         call compute_spectrum_bounds(h_gamma, gamma_bounds, 'sturm', verbose=.false.)
         if (gamma_bounds%sturm_available) then
            have_gamma = .true.
            hgamma_min = gamma_bounds%e_min_sturm
            hgamma_max = gamma_bounds%e_max_sturm
         else
            have_gamma = .true.
            hgamma_min = gamma_bounds%e_min_gershgorin
            hgamma_max = gamma_bounds%e_max_gershgorin
         end if
         deallocate(h_gamma)
      end if

      call select_bounds_interval(algo, g_min, g_max, have_gamma, hgamma_min, hgamma_max, this%bounds%e_min, this%bounds%e_max)
      call apply_bounds_scaling(this%bounds%e_min, this%bounds%e_max, this%bounds%scaling)

      msg = 'Chebyshev bounds algorithm='//trim(algo)//' Gershgorin=['// &
            trim(fmt('f12.6', g_min))//','//trim(fmt('f12.6', g_max))//']'
      call g_logger%info(msg, __FILE__, __LINE__)
      if (have_gamma) then
         msg = 'Chebyshev bounds H(Gamma)=['//trim(fmt('f12.6', hgamma_min))//','// &
               trim(fmt('f12.6', hgamma_max))//']'
      else
         msg = 'Chebyshev bounds H(Gamma)=unavailable, fallback to Gershgorin'
      end if
      call g_logger%info(msg, __FILE__, __LINE__)
      msg = 'Chebyshev final scaled bounds=['//trim(fmt('f12.6', this%bounds%e_min))//','// &
            trim(fmt('f12.6', this%bounds%e_max))//'] scale='//trim(fmt('f8.4', this%bounds%scaling))
      call g_logger%info(msg, __FILE__, __LINE__)

      if (verb) then
         write(msg, '(A,L1)') 'Chebyshev bounds detailed mode, H(Gamma) available=', have_gamma
         call g_logger%info(trim(msg), __FILE__, __LINE__)
      end if
   end subroutine compute_hamiltonian_bounds
end module hamiltonian_mod
