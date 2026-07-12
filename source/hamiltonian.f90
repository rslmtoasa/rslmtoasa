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
   use basis_mod, only: nb, norb, spin_off, lmax_basis
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
      complex(rp), dimension(:, :, :, :), allocatable :: ee, eeo, eeoee, eecc
      !> Local Hamiltonian
      complex(rp), dimension(:, :, :, :), allocatable :: hall, hallo, hallcc
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
      !> Spin-spiral handled in k-space (GBT). When true, ham0m_nc does NOT rotate
      !> the site moments by the q_ss spiral phase (that real-space, absolute-position
      !> construction is not translationally invariant and would be wrong for a Bloch
      !> sum); the spiral is instead applied as a twist in the reciprocal module.
      !> When false (default), q_ss rotates the moments -> real-space spin spiral.
      logical :: gbt_kspace
      !> Add orbital polarization to Hamiltonian
      logical :: orb_pol
      !> Optional two-centre combined correction
      logical :: ccor_2c
      !> Spectral linearization energy for two-centre CCOR (Ry)
      real(rp) :: ccor_elin
      !> VMT strategy for two-centre CCOR
      character(len=16) :: ccor_vmt_mode
      !> Print CCOR diagnostics
      logical :: ccor_debug
      !> Promote CCOR warnings to fatal errors where supported
      logical :: ccor_strict
      !> Bulk Hamiltonian (backup for rotation)
      complex(rp), dimension(:, :, :, :), allocatable :: ee_glob, eeo_glob, eecc_glob
      !> Local Hamiltonian (backup for rotation)
      complex(rp), dimension(:, :, :, :), allocatable :: hall_glob, hallo_glob, hallcc_glob
    !!> Spin-orbit coupling Hamiltonian (backup for rotation)
      !complex(rp), dimension(:, :, :), allocatable :: lsham
      !> Gravity center Hamiltonian (backup for rotation)
      complex(rp), dimension(:, :, :), allocatable :: enim_glob
      !> Velocity operators
      complex(rp), dimension(:, :, :, :), allocatable :: v_a, v_b, js_a, jl_a, vo_a, vo_b, jso_a, jlo_a
      character(len=10) :: js_alpha, jl_alpha
      real(rp), dimension(3) :: v_alpha, v_beta, q_ss
      real(rp) :: theta_ss
      real(rp), dimension(:), allocatable :: theta_ss_sublattice, phi_ss_sublattice
      real(rp), dimension(:), allocatable :: velocity_scale
      character(len=16) :: hubbard_u_potential_form
      !> Sparse Real Space Hamiltonian (dense legacy format)
      complex(rp), dimension(:, :), allocatable :: h_sparse
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
      !> Hamiltonian export format: 'none', 'rs2pao', 'python'
      character(len=16) :: export
   contains
      procedure :: build_lsham
      procedure :: build_bulkham
      procedure :: build_locham
      procedure :: build_ccor_bulk
      procedure :: build_ccor_local
      procedure :: build_obarm
      procedure :: build_enim
      procedure :: build_from_paoflow
      procedure :: build_from_paoflow_opt
      procedure :: build_realspace_velocity_operators
      procedure :: build_realspace_spin_operators
      procedure :: build_realspace_spin_torque_operators
      procedure :: build_realspace_orbital_velocity_operators
      procedure :: build_realspace_orbital_torque_operators
      procedure :: torque_operator_collinear
      procedure :: rs2pao
      procedure :: export_rs_tb_all
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


   interface
   !> @brief Construct a Hamiltonian object from an initialized charge object.
   !> @details Wires the Hamiltonian to charge, lattice, and control state, resets
   !>          all owned arrays/flags, then reads the &hamiltonian namelist. This
   !>          object supplies real-space hopping blocks to recursion and k-space
   !>          Fourier paths.
   !> @param[in] charge_obj Charge object containing lattice and potential state.
   !> @return Initialized Hamiltonian object.
   module function constructor(charge_obj) result(obj)
      type(hamiltonian) :: obj
      type(charge), target, intent(in) :: charge_obj

   end function constructor

   !> @brief Finalize a Hamiltonian object.
   !> @details Releases spin-orbit, hopping, overlap, velocity, Hubbard, CCOR,
   !>          export, and cache arrays owned by the object.
   !> @param[inout] this Hamiltonian object being finalized.
   module subroutine destructor(this)
      type(hamiltonian) :: this
   end subroutine destructor

   !> @brief Read the &hamiltonian namelist and install Hamiltonian options.
   !> @details Parses HOH, local-axis rotation, CCOR, velocity directions, spectral
   !>          bounds, export mode, and Hubbard U/J/V inputs. Hubbard inputs are
   !>          accepted in eV and converted to internal Ry units.
   !> @param[inout] this Hamiltonian object whose control%fname selects the input file.
   !> @note This is a true input boundary and may raise fatal diagnostics for invalid options.
   module subroutine build_from_file(this)
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

   end subroutine build_from_file

   !> @brief Reset Hamiltonian flags, arrays, and cached bounds to defaults.
   !> @details Restores namelist-controlled options to their baseline values and
   !>          clears allocatable storage so a constructor can rebuild the object
   !>          for the current lattice/control state.
   !> @param[inout] this Hamiltonian object to reset.
   module subroutine restore_to_default(this)
      class(hamiltonian) :: this

   end subroutine restore_to_default

   !> @brief Build real-space orbital-current operator blocks.
   !> @details Forms orbital velocity operators from angular-momentum matrices and
   !>          real-space hopping geometry for orbital transport and orbital torque
   !>          workflows.
   !> @param[inout] this Hamiltonian object; fills jl_a/jl_b-style operator arrays.
   !> @note Uses the same site/neighbor layout as ee/hall velocity blocks.
   module subroutine build_realspace_orbital_velocity_operators(this)
      class(hamiltonian), intent(inout) :: this
   
      integer :: ntype, ia, nr, m
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:,:), tmp2(:,:), L_op(:,:)  ! Temp matrices
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz

   end subroutine build_realspace_orbital_velocity_operators

   !> @brief Build real-space spin-current operator blocks.
   !> @details Combines spin matrices with the velocity-operator layout so
   !>          stochastic conductivity can evaluate spin-current correlations.
   !> @param[inout] this Hamiltonian object; fills js_a and related spin-current arrays.
   module subroutine build_realspace_spin_operators(this)
      class(hamiltonian), intent(inout) :: this
   
      integer :: ntype, ia, nr, m, ji, ja, atom_neighbor
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:, :), tmp2(:, :), S_op(:, :)  ! Temp matrices for partial products
   
      ! Derive dimension from your velocity array:
   end subroutine build_realspace_spin_operators

   !> @brief Build real-space spin-torque operator blocks.
   !> @details Forms torque-current operators from spin operators and local
   !>          Hamiltonian/SOC blocks for spin-torque response calculations.
   !> @param[inout] this Hamiltonian object; fills js_a/jso_a torque arrays.
   module subroutine build_realspace_spin_torque_operators(this)
      class(hamiltonian), intent(inout) :: this

      integer :: ntype, ia, nr, m, ji, ja, atom_neighbor, ino
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:, :), tmp2(:, :), S_op(:, :)  ! Temp matrices for partial products
      complex(rp), dimension(nb, nb) :: locham

      ! Derive dimension from your velocity array:
   end subroutine build_realspace_spin_torque_operators

   !> @brief Build real-space orbital-torque operator blocks.
   !> @details Forms orbital torque-current operators from angular-momentum
   !>          matrices and local Hamiltonian/SOC blocks for orbital response.
   !> @param[inout] this Hamiltonian object; fills jl_a/jlo_a torque arrays.
   module subroutine build_realspace_orbital_torque_operators(this)
      class(hamiltonian), intent(inout) :: this

      integer :: ntype, ia, nr, m, ji, ja, atom_neighbor, ino
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:, :), tmp2(:, :), L_op(:, :)  ! Temp matrices for partial products
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz
      complex(rp), dimension(nb, nb) :: locham

      ! Derive dimension from your velocity array:
   end subroutine build_realspace_orbital_torque_operators

   !> @brief Build charge velocity operators in real-space hopping layout.
   !> @details Projects bond displacement vectors onto v_alpha/v_beta directions
   !>          and weights hopping blocks to form v_a/v_b for Kubo conductivity.
   !> @param[inout] this Hamiltonian object; fills v_a, v_b, vo_a, and vo_b.
   !> @note Honors velocity_scale and the HOH companion-operator layout.
   module subroutine build_realspace_velocity_operators(this)
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
   end subroutine build_realspace_velocity_operators

   !> @brief Build onsite spin-orbit coupling Hamiltonian blocks.
   !> @details Converts angular-momentum operators to the LMTO basis and combines
   !>          them with per-type SOC strengths to populate lsham for real-space
   !>          and reciprocal Hamiltonian construction.
   !> @param[inout] this Hamiltonian object; fills lsham(:,:,itype).
   module subroutine build_lsham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k
      complex(rp) :: prefac, sg
      real(rp) :: soc_p, soc_d
      complex(rp), dimension(2) :: rac
      complex(rp), dimension(norb, norb) :: Lx, Ly, Lz
      real(rp) :: lz_loc
      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
   end subroutine build_lsham

   !> @brief Build the collinear spin-orbit torque operator.
   !> @details Evaluates the commutator-style torque operator T=[o,H_so] in the
   !>          spin-orbital basis for workflows that need magnetic torques.
   !> @param[inout] this Hamiltonian object; fills tmat.
   module subroutine torque_operator_collinear(this)
      !
      class(hamiltonian), intent(inout) :: this
      !
      ! Local variables
      integer :: i, j, k
      complex(rp) :: prefac, sg, soc_p, soc_d
      complex(rp), dimension(norb, norb) :: Lx, Ly, Lz
      !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
   end subroutine torque_operator_collinear

   !> @brief Build overlap-bar matrices for the orthogonal representation.
   !> @details Converts potential/overlap information from symbolic atoms into
   !>          onsite obarm blocks used by HOH and representation transforms.
   !> @param[inout] this Hamiltonian object; fills obarm(:,:,itype).
   module subroutine build_obarm(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      complex(rp), dimension(norb, norb) :: obm0, obm1
      complex(rp), dimension(3) :: mom
      integer :: ntype ! Atom type index
      integer :: l, m ! Orbital index

   end subroutine build_obarm

   !> @brief Build onsite e_nu center matrices.
   !> @details Assembles spin-resolved gravity-center/e_nu blocks from symbolic
   !>          atom potentials for the orthogonalized Hamiltonian correction.
   !> @param[inout] this Hamiltonian object; fills enim(:,:,itype).
   module subroutine build_enim(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      complex(rp), dimension(norb, norb) :: em0, em1
      complex(rp), dimension(norb) :: ex0, ex1
      complex(rp), dimension(3) :: mom
      complex(rp) :: eu, ed
      integer :: ntype ! Atom type index
      integer :: l, m ! Orbital index

   end subroutine build_enim

   !> @brief Build bulk real-space Hamiltonian hopping blocks.
   !> @details Assembles ee, eeo, and related per-type neighbor blocks from the
   !>          lattice structure constants and symbolic-atom potentials. Real-space
   !>          recursion and reciprocal Fourier transforms consume these blocks.
   !> @param[inout] this Hamiltonian object; fills bulk hopping arrays.
   !> @note Also prepares optional Hubbard/CCOR contributions when enabled.
   module subroutine build_bulkham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia
      integer :: ntype

   end subroutine build_bulkham

   !> @brief Build local-cluster Hamiltonian hopping blocks.
   !> @details Assembles hall/hallo blocks for impurity or surface local regions
   !>          where atom-specific local geometry replaces bulk type blocks.
   !> @param[inout] this Hamiltonian object; fills local hopping arrays.
   module subroutine build_locham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: it, ino, nr, nlim, m, i, j, ja, ji

   end subroutine build_locham

   !> @brief Add two-centre combined-correction blocks to the bulk Hamiltonian.
   !> @details Builds CCOR contributions for each bulk atom type and neighbor
   !>          shell, using scalar or noncollinear pair-block helpers as required
   !>          by the magnetic state.
   !> @param[inout] this Hamiltonian object; fills eecc.
   !> @note validate_ccor_inputs must accept the current CCOR mode before use.
   module subroutine build_ccor_bulk(this)
      class(hamiltonian), intent(inout) :: this
      integer :: ntype, ia, ino, nr, m, jj, it, jt
      complex(rp), dimension(nb, nb) :: hcc

   end subroutine build_ccor_bulk

   !> @brief Add two-centre combined-correction blocks to local Hamiltonian regions.
   !> @details Builds CCOR contributions for local/surface/impurity neighbor
   !>          blocks in the same layout as hall, preserving local atom indexing.
   !> @param[inout] this Hamiltonian object; fills hallcc.
   module subroutine build_ccor_local(this)
      class(hamiltonian), intent(inout) :: this
      integer :: nlim, ino, nr, m, jj, it, jt
      complex(rp), dimension(nb, nb) :: hcc

   end subroutine build_ccor_local

   !> @brief Validate the current CCOR namelist and lattice state.
   !> @details Checks that two-centre combined-correction parameters are meaningful
   !>          for the active calculation and emits warnings or fatal diagnostics
   !>          according to ccor_strict.
   !> @param[in] this Hamiltonian object with CCOR flags and lattice state.
   module subroutine validate_ccor_inputs(this)
      class(hamiltonian), intent(in) :: this

   end subroutine validate_ccor_inputs

   !> @brief Build one scalar-relativistic CCOR pair block.
   !> @details Computes the two-centre combined-correction hopping contribution
   !>          between atom sites ia/ja and types it/jt for neighbor index m.
   !> @param[in] this Hamiltonian object containing potentials and CCOR settings.
   !> @param[in] ia Site index of the source atom.
   !> @param[in] ja Site index of the neighbor atom.
   !> @param[in] it Atom type of ia.
   !> @param[in] jt Atom type of ja.
   !> @param[in] ino Bravais/type index used by the target Hamiltonian layout.
   !> @param[in] m Neighbor-shell index.
   !> @param[out] hcc CCOR block in spin-orbital basis.
   module subroutine build_ccor_pair_block_scalar(this, ia, ja, it, jt, ino, m, hcc)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, ja, it, jt, ino, m
      complex(rp), dimension(nb, nb), intent(out) :: hcc

   end subroutine build_ccor_pair_block_scalar

   !> @brief Build one noncollinear CCOR pair block.
   !> @details Computes spin-dependent two-centre combined-correction terms using
   !>          local moments, D components, and VMT coefficients for the atom pair.
   !> @param[in] this Hamiltonian object containing potentials and magnetic state.
   !> @param[in] ia Site index of the source atom.
   !> @param[in] ja Site index of the neighbor atom.
   !> @param[in] it Atom type of ia.
   !> @param[in] jt Atom type of ja.
   !> @param[in] ino Bravais/type index used by the target Hamiltonian layout.
   !> @param[in] m Neighbor-shell index.
   !> @param[out] hcc CCOR block in spin-orbital basis.
   module subroutine build_ccor_pair_block_noncollinear(this, ia, ja, it, jt, ino, m, hcc)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, ja, it, jt, ino, m
      complex(rp), dimension(nb, nb), intent(out) :: hcc

      complex(rp), dimension(norb, norb) :: s_block, sdot_block
      complex(rp), dimension(norb, norb, 4) :: dcomp, ddotcomp, kcomp
      real(rp), dimension(norb, 0:2) :: ccd_i, ccd_j
      real(rp), dimension(3) :: mom_i
      real(rp) :: lambda
      integer :: ilm, jlm, idir

   end subroutine build_ccor_pair_block_noncollinear

      !> @brief Build the surface-mode noncollinear CCOR pair block.
      !> @details Combines endpoint and pair VMT terms with D and D-dot components
      !>          for surfaces or local regions where global bulk assumptions do
      !>          not apply.
      !> @param[in] this Hamiltonian object containing CCOR mode and moments.
      !> @param[in] ia Site index of the source atom.
      !> @param[in] ja Site index of the neighbor atom.
      !> @param[in] it Atom type of ia.
      !> @param[in] jt Atom type of ja.
      !> @param[in] m Neighbor-shell index.
      !> @param[in] dcomp Spin decomposed D matrix components.
      !> @param[in] ddotcomp Spin decomposed D-dot matrix components.
      !> @param[in] ccd_i Combined-correction coefficients for type it.
      !> @param[in] ccd_j Combined-correction coefficients for type jt.
      !> @param[out] hcc CCOR block in spin-orbital basis.
      module subroutine build_ccor_pair_surface_block(this, ia, ja, it, jt, m, dcomp, ddotcomp, ccd_i, ccd_j, hcc)
         class(hamiltonian), intent(in) :: this
         integer, intent(in) :: ia, ja, it, jt, m
         complex(rp), dimension(norb, norb, 4), intent(in) :: dcomp, ddotcomp
         real(rp), dimension(norb, 0:2), intent(in) :: ccd_i, ccd_j
         complex(rp), dimension(nb, nb), intent(out) :: hcc

         complex(rp), dimension(norb, norb, 4) :: hcomp
         complex(rp), dimension(4) :: lambda_i, lambda_j, dloc, ddotloc, term_l, term_r, onsite
         real(rp), dimension(2) :: lambda_pair
         real(rp), dimension(3) :: mom_i, mom_j
         integer :: ilm, jlm, idir

   end subroutine build_ccor_pair_surface_block

   !> @brief Decompose a structure-constant block into CCOR D components.
   !> @details Projects the orbital structure block onto scalar and spin-vector
   !>          components used by noncollinear CCOR pair-block construction.
   !> @param[in] this Hamiltonian object containing moments and spin-spiral state.
   !> @param[in] ia Site index of the source atom.
   !> @param[in] ja Site index of the neighbor atom.
   !> @param[in] it Atom type of ia.
   !> @param[in] jt Atom type of ja.
   !> @param[in] s_block Structure-constant block in orbital basis.
   !> @param[out] dcomp Decomposed D components.
   module subroutine build_ccor_d_components(this, ia, ja, it, jt, s_block, dcomp)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, ja, it, jt
      complex(rp), dimension(norb, norb), intent(in) :: s_block
      complex(rp), dimension(norb, norb, 4), intent(out) :: dcomp

      real(rp), dimension(3) :: mom_i, mom_j
      complex(rp), dimension(3) :: cross
      complex(rp) :: dotp
      integer :: ilm, jlm, idir

   end subroutine build_ccor_d_components

   !> @brief Compute orbital-resolved CCOR coefficients for one atom type.
   !> @details Evaluates the energy-linearized combined-correction coefficients
   !>          from potential parameters for the specified site/type.
   !> @param[in] this Hamiltonian object containing symbolic-atom potentials.
   !> @param[in] ia Site index used for local moment/context.
   !> @param[in] itype Atom type whose coefficients are requested.
   !> @param[out] ccd Coefficients indexed by orbital and angular channel.
   module subroutine build_ccor_coefficients(this, ia, itype, ccd)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, itype
      real(rp), dimension(norb, 0:2), intent(out) :: ccd

      integer :: ilm, l, alpha_idx
      real(rp) :: sw, dl, a, adot, l2p3, l2p1, l2m1, t1, t2, t3, avw

   end subroutine build_ccor_coefficients

   !> @brief Normalize a CCOR spin-dot product.
   !> @details Applies the selected CCOR normalization convention to a raw
   !>          spin-product scalar before it enters pair-block construction.
   !> @param[in] this Hamiltonian object with CCOR normalization settings.
   !> @param[in] sdot_raw Raw complex spin-dot value.
   !> @return Normalized spin-dot value.
   module function normalize_ccor_sdot(this, sdot_raw) result(sdot_cc)
      class(hamiltonian), intent(in) :: this
      complex(rp), intent(in) :: sdot_raw
      complex(rp) :: sdot_cc
      real(rp) :: avw

   end function normalize_ccor_sdot

      !> @brief Return the scalar VMT value for CCOR.
      !> @details Selects the configured scalar VMT strategy and reduces any
      !>          spin-resolved surface estimate to the scalar value needed by
      !>          scalar CCOR blocks.
      !> @param[in] this Hamiltonian object with ccor_vmt_mode.
      !> @return Scalar VMT value in internal units.
      module function ccor_vmt_scalar(this) result(vmt)
         class(hamiltonian), intent(in) :: this
         real(rp) :: vmt
         real(rp), dimension(2) :: vmt_spin

   end function ccor_vmt_scalar

      !> @brief Compute a global spin-resolved surface VMT estimate.
      !> @details Averages endpoint surface VMT values over symbolic atom types and
      !>          multiplicities to provide one spin-resolved correction scale.
      !> @param[in] this Hamiltonian object containing symbolic atoms.
      !> @return Spin-resolved VMT pair.
      module function ccor_vmt_global_surface(this) result(vmt_spin)
         class(hamiltonian), intent(in) :: this
         real(rp), dimension(2) :: vmt_spin
         real(rp), dimension(2) :: endpoint
         real(rp) :: weight, weight_sum
         integer :: itype, nsite, multiplicity

   end function ccor_vmt_global_surface

      !> @brief Compute a pair-specific spin-resolved surface VMT estimate.
      !> @details Blends endpoint VMT values for two atom types with pair weights
      !>          for pair_surface CCOR mode.
      !> @param[in] this Hamiltonian object containing symbolic atoms.
      !> @param[in] itype First atom type.
      !> @param[in] jtype Second atom type.
      !> @return Spin-resolved VMT pair.
      module function ccor_vmt_pair_surface(this, itype, jtype) result(vmt_spin)
         class(hamiltonian), intent(in) :: this
         integer, intent(in) :: itype, jtype
         real(rp), dimension(2) :: vmt_spin, vi, vj
         real(rp) :: wi, wj

   end function ccor_vmt_pair_surface

      !> @brief Compute an endpoint spin-resolved surface VMT estimate.
      !> @details Extracts the local surface VMT proxy for one symbolic atom type
      !>          from the available potential data.
      !> @param[in] this Hamiltonian object containing symbolic atoms.
      !> @param[in] itype Atom type whose endpoint value is requested.
      !> @return Spin-resolved VMT pair.
      module function ccor_vmt_endpoint_surface(this, itype) result(vmt_spin)
         class(hamiltonian), intent(in) :: this
         integer, intent(in) :: itype
         real(rp), dimension(2) :: vmt_spin
         real(rp) :: avg, diff

   end function ccor_vmt_endpoint_surface

      !> @brief Compute a scalar CCOR VMT estimate from Madelung potentials.
      !> @details Averages available vmad-like potential information over atom
      !>          types for the vmad_scalar CCOR mode.
      !> @param[in] this Hamiltonian object containing symbolic atoms.
      !> @return Scalar VMT value in internal units.
      module function ccor_vmt_scalar_from_vmad(this) result(vmt)
         class(hamiltonian), intent(in) :: this
         real(rp) :: vmt
         real(rp) :: weight, weight_sum
         integer :: itype

   end function ccor_vmt_scalar_from_vmad

      !> @brief Expand spin-resolved lambda values into Pauli components.
      !> @details Converts the up/down lambda pair and local moment direction into
      !>          scalar/vector components used by noncollinear CCOR algebra.
      !> @param[in] lambda_pair Spin-resolved lambda values.
      !> @param[in] mom Local magnetic moment direction.
      !> @param[out] lambda_comp Scalar and vector lambda components.
      module subroutine ccor_lambda_components(lambda_pair, mom, lambda_comp)
         real(rp), dimension(2), intent(in) :: lambda_pair
         real(rp), dimension(3), intent(in) :: mom
         complex(rp), dimension(4), intent(out) :: lambda_comp
         real(rp) :: lambda0, lambda1

   end subroutine ccor_lambda_components

      !> @brief Multiply two scalar/vector spin-component tuples.
      !> @details Applies Pauli-matrix product algebra to compact four-component
      !>          spin decompositions used by noncollinear CCOR construction.
      !> @param[in] a Left scalar/vector spin tuple.
      !> @param[in] b Right scalar/vector spin tuple.
      !> @param[out] c Product tuple.
      module subroutine ccor_spin_product(a, b, c)
         complex(rp), dimension(4), intent(in) :: a, b
         complex(rp), dimension(4), intent(out) :: c

   end subroutine ccor_spin_product

      !> @brief Apply the spin-spiral rotation to a local moment for CCOR.
      !> @details Rotates the supplied magnetic moment according to the atom
      !>          position and spin-spiral q/theta parameters before pair-block use.
      !> @param[in] this Hamiltonian object containing lattice and spin-spiral state.
      !> @param[in] ia Site index whose position sets the spin-spiral phase.
      !> @param[inout] mom Moment vector to rotate in place.
      module subroutine ccor_apply_spin_spiral(this, ia, mom)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia
      real(rp), dimension(3), intent(inout) :: mom
      real(rp), dimension(3) :: r_ia

   end subroutine ccor_apply_spin_spiral

   !> @brief Return angular momentum l from a packed orbital index.
   !> @details Maps the LMTO orbital index ilm to its shell quantum number for
   !>          coefficient lookup and export helpers.
   !> @param[in] ilm Packed orbital index.
   !> @return Angular momentum shell l.
   module integer pure function orbital_l_from_index(ilm) result(l)
      integer, intent(in) :: ilm
      integer :: lp1

   end function orbital_l_from_index

   !> @brief Emit diagnostics for generated CCOR blocks.
   !> @details Summarizes CCOR coefficient and block magnitudes for debugging the
   !>          two-centre correction without changing Hamiltonian data.
   !> @param[in] this Hamiltonian object containing CCOR settings.
   !> @param[in] hcc CCOR block array to summarize.
   !> @param[in] label Human-readable label for the diagnostic message.
   module subroutine log_ccor_debug(this, hcc, label)
      class(hamiltonian), intent(in) :: this
      complex(rp), dimension(:, :, :, :), intent(in) :: hcc
      character(len=*), intent(in) :: label
      real(rp), dimension(norb, 0:2) :: ccd
         real(rp) :: vmt, lambda, rms_ccd2, max_ccd2, avw
         real(rp), dimension(2) :: vmt_spin
      integer :: itype, count_ccd
      character(len=256) :: msg

   end subroutine log_ccor_debug

   !> @brief Export the real-space Hamiltonian in legacy PAOFLOW layout.
   !> @details Writes hopping records with lattice-cell indices and global orbital
   !>          numbers so external PAOFLOW-style tooling can consume the RS-LMTO
   !>          real-space Hamiltonian.
   !> @param[inout] this Hamiltonian object containing built hopping blocks.
   !> @note This is the legacy export entry point; export_rs_tb_all is the newer dispatcher.
   module subroutine rs2pao(this)
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
   end subroutine rs2pao

   !> @brief Import a PAOFLOW-format real-space Hamiltonian from paoham.dat.
   !> @details Reads legacy seven-column hopping records and fills the bulk
   !>          Hamiltonian blocks for PAOFLOW-to-real-space post-processing routes.
   !> @param[inout] this Hamiltonian object; replaces bulk hopping arrays from file data.
   !> @note PAOFLOW import is used by paoflow2rs, exchange_p2rs, and conductivity_p2rs.
   module subroutine build_from_paoflow_opt(this)
      implicit none
      type hamData
         integer :: idxi, idxj, idxk
         integer :: orbl, orbm
         real :: dumre, dumcmplx
      end type hamData
      class(hamiltonian), intent(inout) :: this
   end subroutine build_from_paoflow_opt

   !> @brief Import a PAOFLOW-format real-space Hamiltonian with the legacy reader.
   !> @details Maps PAOFLOW orbital/cell records back onto RS-LMTO neighbor blocks.
   !>          Kept for compatibility with older import paths.
   !> @param[inout] this Hamiltonian object; fills bulk hopping arrays from file data.
   module subroutine build_from_paoflow(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia, iia, jja
      integer :: jj, jt, orbl, orbm, idxi, idxj, idxk
      integer :: ntype, iostat1, iostat2, iostatus, n_atoms, max_orbital
      real(rp), dimension(3) :: vet, vetpao, cri_dir, crj_dir, cri_cart, crj_cart
      integer, dimension(maxval(this%charge%lattice%nn(:, 1)) + 1, 3) :: idxup, idxdw, idx
      real(rp) :: dumre, dumcmplx

   end subroutine build_from_paoflow

   !> @brief Build a noncollinear spin-orbital hopping block from orbital data.
   !> @details Lifts an orbital hopping/structure block hhh into the spin-orbital
   !>          basis using the local moments of atom sites ia and ja, including
   !>          spin-spiral phase handling where enabled.
   !> @param[inout] this Hamiltonian object containing magnetic moment state.
   !> @param[in] ia Site index of the source atom.
   !> @param[in] ja Site index of the neighbor atom.
   !> @param[in] it Atom type of ia.
   !> @param[in] jt Atom type of ja.
   !> @param[in] vet Bond vector between the sites.
   !> @param[in] hhh Orbital hopping/structure block.
   module subroutine ham0m_nc(this, ia, ja, it, jt, vet, hhh)
      class(hamiltonian), intent(inout) :: this
      ! Input
      integer, intent(in) :: ia, ja ! Atom sites i and j
      integer, intent(in) :: it, jt ! Type of atom i and j
      real(rp), dimension(3), intent(in) :: vet
      real(rp), dimension(norb, norb), intent(in) :: hhh
      ! Local Variables
      integer :: i, j, ilm, jlm, m
      real(rp), dimension(3) :: mom_ia, mom_ja
      real(rp), dimension(3) :: r_ia, r_ja
      complex(rp), dimension(3) :: cross
      complex(rp), dimension(norb, norb) :: hhhc
      complex(rp), dimension(this%charge%lattice%ntype, 3) :: momc
      complex(rp) :: dot
      real(rp) :: vv

   end subroutine ham0m_nc

   !> @brief Build noncollinear local hopping data for one cluster atom.
   !> @details Walks the neighbor list around atom ia, obtains orbital hopping
   !>          blocks, converts them with ham0m_nc, and stores the resulting local
   !>          Hamiltonian/field data.
   !> @param[inout] this Hamiltonian object; updates local noncollinear arrays.
   !> @param[in] ia Cluster atom index.
   !> @param[in] nr Number of neighbors considered.
   !> @param[in] ino Bravais/type index for ia.
   !> @param[in] ntype Atom type index.
   module subroutine chbar_nc(this, ia, nr, ino, ntype)
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

   end subroutine chbar_nc

   !> @brief Find the orbital hopping block matching a neighbor vector.
   !> @details Searches the precomputed neighbor/hopping table for the vector vet
   !>          and returns the matching orbital block and neighbor index.
   !> @param[inout] this Hamiltonian object containing lattice tolerances.
   !> @param[in] vet Candidate neighbor vector.
   !> @param[in] nr Number of neighbors in the search list.
   !> @param[inout] hhh Orbital hopping block output.
   !> @param[in] m Neighbor counter being resolved.
   !> @param[in] ia Cluster atom index.
   !> @param[in] jn Neighbor-list index.
   !> @param[out] ni Matched neighbor index.
   !> @param[in] ham_vec Precomputed neighbor vectors.
   module subroutine hmfind(this, vet, nr, hhh, m, ia, jn, ni, ham_vec)
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

   end subroutine hmfind

   !> @brief Convert a global orbital index to site and local-orbital indices.
   !> @details Used by PAOFLOW import/export helpers to translate flat orbital
   !>          numbering into RS-LMTO site-major layout.
   !> @param[in] orb Global orbital index.
   !> @param[out] i_out Site index.
   !> @param[out] ia_out Local orbital index on the site.
   !> @param[in] n_atoms Number of atoms in the exported/imported cell.
   !> @param[in] max_orbital Number of orbitals per atom in the flat layout.
   module subroutine orb2site(orb, i_out, ia_out, n_atoms, max_orbital)
      integer, intent(in) :: orb, n_atoms, max_orbital
      integer, intent(out) :: i_out, ia_out

   end subroutine orb2site

   !> @brief Convert site and local-orbital indices to a global orbital index.
   !> @details Used by PAOFLOW import/export helpers to write flat orbital
   !>          numbering from RS-LMTO site-major data.
   !> @param[in] i_in Site index.
   !> @param[in] ia_in Local orbital index on the site.
   !> @param[out] orb_out Global orbital index.
   !> @param[in] n_atoms Number of atoms in the exported/imported cell.
   !> @param[in] max_orbital Number of orbitals per atom in the flat layout.
   module subroutine site2orb(i_in, ia_in, orb_out, n_atoms, max_orbital)
      integer, intent(in) :: i_in, ia_in, n_atoms, max_orbital
      integer, intent(out) :: orb_out

   end subroutine site2orb

   !> @brief Rotate Hamiltonian blocks into a local magnetic axis.
   !> @details Saves global hopping/SOC data and rotates spin blocks so a local
   !>          moment direction can be treated as the quantization axis in legacy
   !>          recursion paths.
   !> @param[inout] this Hamiltonian object; updates rotated block arrays.
   !> @param[in] m_loc Local magnetic moment direction.
   !> @note Call rotate_from_local_axis to restore the global representation.
   module subroutine rotate_to_local_axis(this, m_loc)
      use math_mod, only: rotmag_loc
      class(hamiltonian), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: m_loc

      ! Local variables
      integer :: sdim
      ! Rotate Hamiltonian to local axis if wanted
   end subroutine rotate_to_local_axis

   !> @brief Restore Hamiltonian blocks from a local-axis rotation.
   !> @details Rotates local-axis Hamiltonian data back to the global spin frame
   !>          after an atom-specific recursion calculation.
   !> @param[inout] this Hamiltonian object; restores global block arrays.
   !> @param[in] m_loc Local magnetic moment direction used for the rotation.
   module subroutine rotate_from_local_axis(this, m_loc)
      use math_mod, only: rotmag_loc
      class(hamiltonian), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: m_loc

      ! Local variables
      integer :: sdim
      ! Rotate Hamiltonian to local axis if wanted
   end subroutine rotate_from_local_axis

   !> @brief Compute onsite Hubbard-U/J potential corrections.
   !> @details Builds the spin-orbital onsite correction from symbolic-atom
   !>          Hubbard inputs, including Liechtenstein and ACBN0-style forms and
   !>          optional self-consistent U channels.
   !> @param[inout] this Hamiltonian object; fills hubbard_u_pot and related masks.
   !> @note Results are added by the Hamiltonian builders, not by this routine directly.
   module subroutine calculate_hubbard_u_potential_general(this)
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
   end subroutine calculate_hubbard_u_potential_general

   !> @brief Compute intersite Hubbard-V potential corrections.
   !> @details Builds pair-dependent correction blocks from hubbard_v inputs and
   !>          approximate local occupations for later inclusion in real-space
   !>          hopping construction.
   !> @param[inout] this Hamiltonian object; fills hubbard_v_pot.
   module subroutine calculate_hubbard_v_potential(this)
      class(hamiltonian), intent(inout) :: this

      integer :: itype, ia, nr, m, ja, jt, ispin, li, lj, ii, i0, i1, jdim
      real(rp) :: rmin, rcur, tol, occ_up, occ_dn, nn_shell_tol
      real(rp), dimension(3) :: dr
      real(rp), dimension(this%lattice%ntype, 4, 2) :: n_spin
      logical, save :: warned_proxy_once = .false.

   end subroutine calculate_hubbard_v_potential

   !> @brief Estimate spectral bounds for Chebyshev scaling.
   !> @details Computes or selects Hamiltonian energy bounds according to the
   !>          configured bounds algorithm and scaling factor. Chebyshev recursion
   !>          uses these bounds to map H into [-1,1].
   !> @param[inout] this Hamiltonian object; updates this%bounds.
   !> @param[in] verbose Optional flag enabling diagnostic logging.
   module subroutine compute_hamiltonian_bounds(this, verbose)
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

   end subroutine compute_hamiltonian_bounds

   !> @brief Export real-space tight-binding data in all selected formats.
   !> @details Dispatcher for metadata and hopping-record writers used by the
   !>          hamiltonian export namelist option. Supports PAOFLOW legacy and
   !>          Python-friendly real-space records.
   !> @param[in] this Hamiltonian object containing built hopping blocks.
   !> @param[in] basename Optional output basename.
   !> @param[in] tol Optional threshold for skipping tiny records.
   !> @param[in] include_lsham Optional flag to include onsite SOC blocks.
   !> @param[in] transform_sph2cart Optional flag to transform orbital basis.
   module subroutine export_rs_tb_all(this, basename, tol, include_lsham, transform_sph2cart)
      implicit none
      class(hamiltonian), intent(in) :: this
      character(len=*), intent(in), optional :: basename
      real(rp), intent(in), optional :: tol
      logical, intent(in), optional :: include_lsham, transform_sph2cart
   
      character(len=512) :: base
      real(rp) :: eps
      logical :: add_lsham, do_sph2cart
   
   end subroutine export_rs_tb_all

   !> @brief Write metadata for real-space tight-binding exports.
   !> @details Records atom/orbital layout information needed to interpret the
   !>          exported hopping records.
   !> @param[in] this Hamiltonian object containing lattice and basis metadata.
   !> @param[in] filename Metadata output path.
   module subroutine export_rs_tb_metadata(this, filename)
      implicit none
      class(hamiltonian), intent(in) :: this
      character(len=*), intent(in) :: filename
   
      integer :: u, i, ios, n_atoms, max_orbital
   
   end subroutine export_rs_tb_metadata

   !> @brief Write PAOFLOW legacy-format hopping records.
   !> @details Emits the seven-column record layout read by build_from_paoflow_opt:
   !>          cell index, global orbital indices, and complex hopping value.
   !> @param[in] this Hamiltonian object containing built hopping blocks.
   !> @param[in] filename Output path.
   !> @param[in] tol Threshold for skipping tiny records.
   !> @param[in] include_lsham Include onsite SOC blocks in the output.
   !> @param[in] transform_sph2cart Transform orbital blocks before writing.
   module subroutine export_rs_paoflow_legacy(this, filename, tol, include_lsham, transform_sph2cart)
      implicit none
      class(hamiltonian), intent(in) :: this
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: tol
      logical, intent(in) :: include_lsham, transform_sph2cart
   
      integer :: u, ios
   
   end subroutine export_rs_paoflow_legacy

   !> @brief Write Python-friendly real-space hopping records.
   !> @details Emits hopping records plus metadata conventions intended for direct
   !>          parsing by scripts and post-processing tools.
   !> @param[in] this Hamiltonian object containing built hopping blocks.
   !> @param[in] filename Output path.
   !> @param[in] tol Threshold for skipping tiny records.
   !> @param[in] include_lsham Include onsite SOC blocks in the output.
   !> @param[in] transform_sph2cart Transform orbital blocks before writing.
   module subroutine export_rs_tb_hoppings(this, filename, tol, include_lsham, transform_sph2cart)
      implicit none
      class(hamiltonian), intent(in) :: this
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: tol
      logical, intent(in) :: include_lsham, transform_sph2cart
   
      integer :: u, ios
   
   end subroutine export_rs_tb_hoppings

   !> @brief Write real-space hopping records to an open unit.
   !> @details Shared implementation for legacy PAOFLOW and Python export modes,
   !>          walking local/bulk neighbor blocks and optionally transforming basis
   !>          or adding onsite spin-orbit terms.
   !> @param[in] this Hamiltonian object containing built hopping blocks.
   !> @param[in] u Open output unit.
   !> @param[in] mode Record format selector.
   !> @param[in] tol Threshold for skipping tiny records.
   !> @param[in] include_lsham Include onsite SOC blocks in the output.
   !> @param[in] transform_sph2cart Transform orbital blocks before writing.
   module subroutine write_rs_tb_records(this, u, mode, tol, include_lsham, transform_sph2cart)
      implicit none
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: u
      character(len=*), intent(in) :: mode
      real(rp), intent(in) :: tol
      logical, intent(in) :: include_lsham, transform_sph2cart
   
      integer :: ntype, ia, nr, k, jj, jbasis
      integer :: i, j, ipao, jpao, n_atoms, max_orbital
      integer :: idx(3)
      logical :: found
      complex(rp) :: hblock(nb, nb)
   
   end subroutine write_rs_tb_records

   !> @brief Resolve a neighbor into integer lattice-cell indices.
   !> @details Matches a real-space neighbor displacement against lattice
   !>          translations so exported hopping records can carry PAOFLOW-style
   !>          integer cell offsets.
   !> @param[in] this Hamiltonian object containing lattice vectors and coordinates.
   !> @param[in] ia Source atom index.
   !> @param[in] jj Neighbor-list entry.
   !> @param[out] idx Integer lattice-cell offset.
   !> @param[out] found True when a matching offset was found within tolerance.
   !> @param[in] tol Matching tolerance.
   module subroutine rs_neighbor_lattice_index(this, ia, jj, idx, found, tol)
      implicit none
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, jj
      integer, intent(out) :: idx(3)
      logical, intent(out) :: found
      real(rp), intent(in) :: tol
   
      integer :: idxi, idxj, idxk, jbasis
      real(rp) :: rij(3), rijtest(3), err, best_err
   
   end subroutine rs_neighbor_lattice_index

   end interface

end module hamiltonian_mod
