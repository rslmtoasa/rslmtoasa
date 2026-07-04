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
   module function constructor(charge_obj) result(obj)
      type(hamiltonian) :: obj
      type(charge), target, intent(in) :: charge_obj

   end function constructor

   module subroutine destructor(this)
      type(hamiltonian) :: this
   end subroutine destructor

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

   module subroutine restore_to_default(this)
      class(hamiltonian) :: this

   end subroutine restore_to_default

   module subroutine build_realspace_orbital_velocity_operators(this)
      class(hamiltonian), intent(inout) :: this
   
      integer :: ntype, ia, nr, m
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:,:), tmp2(:,:), L_op(:,:)  ! Temp matrices
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz

   end subroutine build_realspace_orbital_velocity_operators

   module subroutine build_realspace_spin_operators(this)
      class(hamiltonian), intent(inout) :: this
   
      integer :: ntype, ia, nr, m, ji, ja, atom_neighbor
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:, :), tmp2(:, :), S_op(:, :)  ! Temp matrices for partial products
   
      ! Derive dimension from your velocity array:
   end subroutine build_realspace_spin_operators

   module subroutine build_realspace_spin_torque_operators(this)
      class(hamiltonian), intent(inout) :: this

      integer :: ntype, ia, nr, m, ji, ja, atom_neighbor, ino
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:, :), tmp2(:, :), S_op(:, :)  ! Temp matrices for partial products
      complex(rp), dimension(nb, nb) :: locham

      ! Derive dimension from your velocity array:
   end subroutine build_realspace_spin_torque_operators

   module subroutine build_realspace_orbital_torque_operators(this)
      class(hamiltonian), intent(inout) :: this

      integer :: ntype, ia, nr, m, ji, ja, atom_neighbor, ino
      integer :: hblocksize
      complex(rp), allocatable :: tmp1(:, :), tmp2(:, :), L_op(:, :)  ! Temp matrices for partial products
      complex(rp), dimension(norb, norb) :: mLx, mLy, mLz
      complex(rp), dimension(nb, nb) :: locham

      ! Derive dimension from your velocity array:
   end subroutine build_realspace_orbital_torque_operators

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

   module subroutine build_obarm(this)
      implicit none
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      complex(rp), dimension(norb, norb) :: obm0, obm1
      complex(rp), dimension(3) :: mom
      integer :: ntype ! Atom type index
      integer :: l, m ! Orbital index

   end subroutine build_obarm

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

   module subroutine build_bulkham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia
      integer :: ntype

   end subroutine build_bulkham

   module subroutine build_locham(this)
      class(hamiltonian), intent(inout) :: this
      ! Local variables
      integer :: it, ino, nr, nlim, m, i, j, ja, ji

   end subroutine build_locham

   module subroutine build_ccor_bulk(this)
      class(hamiltonian), intent(inout) :: this
      integer :: ntype, ia, ino, nr, m, jj, it, jt
      complex(rp), dimension(nb, nb) :: hcc

   end subroutine build_ccor_bulk

   module subroutine build_ccor_local(this)
      class(hamiltonian), intent(inout) :: this
      integer :: nlim, ino, nr, m, jj, it, jt
      complex(rp), dimension(nb, nb) :: hcc

   end subroutine build_ccor_local

   module subroutine validate_ccor_inputs(this)
      class(hamiltonian), intent(in) :: this

   end subroutine validate_ccor_inputs

   module subroutine build_ccor_pair_block_scalar(this, ia, ja, it, jt, ino, m, hcc)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, ja, it, jt, ino, m
      complex(rp), dimension(nb, nb), intent(out) :: hcc

   end subroutine build_ccor_pair_block_scalar

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

   module subroutine build_ccor_coefficients(this, ia, itype, ccd)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia, itype
      real(rp), dimension(norb, 0:2), intent(out) :: ccd

      integer :: ilm, l, alpha_idx
      real(rp) :: sw, dl, a, adot, l2p3, l2p1, l2m1, t1, t2, t3, avw

   end subroutine build_ccor_coefficients

   module function normalize_ccor_sdot(this, sdot_raw) result(sdot_cc)
      class(hamiltonian), intent(in) :: this
      complex(rp), intent(in) :: sdot_raw
      complex(rp) :: sdot_cc
      real(rp) :: avw

   end function normalize_ccor_sdot

      module function ccor_vmt_scalar(this) result(vmt)
         class(hamiltonian), intent(in) :: this
         real(rp) :: vmt
         real(rp), dimension(2) :: vmt_spin

   end function ccor_vmt_scalar

      module function ccor_vmt_global_surface(this) result(vmt_spin)
         class(hamiltonian), intent(in) :: this
         real(rp), dimension(2) :: vmt_spin
         real(rp), dimension(2) :: endpoint
         real(rp) :: weight, weight_sum
         integer :: itype, nsite, multiplicity

   end function ccor_vmt_global_surface

      module function ccor_vmt_pair_surface(this, itype, jtype) result(vmt_spin)
         class(hamiltonian), intent(in) :: this
         integer, intent(in) :: itype, jtype
         real(rp), dimension(2) :: vmt_spin, vi, vj
         real(rp) :: wi, wj

   end function ccor_vmt_pair_surface

      module function ccor_vmt_endpoint_surface(this, itype) result(vmt_spin)
         class(hamiltonian), intent(in) :: this
         integer, intent(in) :: itype
         real(rp), dimension(2) :: vmt_spin
         real(rp) :: avg, diff

   end function ccor_vmt_endpoint_surface

      module function ccor_vmt_scalar_from_vmad(this) result(vmt)
         class(hamiltonian), intent(in) :: this
         real(rp) :: vmt
         real(rp) :: weight, weight_sum
         integer :: itype

   end function ccor_vmt_scalar_from_vmad

      module subroutine ccor_lambda_components(lambda_pair, mom, lambda_comp)
         real(rp), dimension(2), intent(in) :: lambda_pair
         real(rp), dimension(3), intent(in) :: mom
         complex(rp), dimension(4), intent(out) :: lambda_comp
         real(rp) :: lambda0, lambda1

   end subroutine ccor_lambda_components

      module subroutine ccor_spin_product(a, b, c)
         complex(rp), dimension(4), intent(in) :: a, b
         complex(rp), dimension(4), intent(out) :: c

   end subroutine ccor_spin_product

      module subroutine ccor_apply_spin_spiral(this, ia, mom)
      class(hamiltonian), intent(in) :: this
      integer, intent(in) :: ia
      real(rp), dimension(3), intent(inout) :: mom
      real(rp), dimension(3) :: r_ia

   end subroutine ccor_apply_spin_spiral

   module integer pure function orbital_l_from_index(ilm) result(l)
      integer, intent(in) :: ilm
      integer :: lp1

   end function orbital_l_from_index

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

   module subroutine build_from_paoflow_opt(this)
      implicit none
      type hamData
         integer :: idxi, idxj, idxk
         integer :: orbl, orbm
         real :: dumre, dumcmplx
      end type hamData
      class(hamiltonian), intent(inout) :: this
   end subroutine build_from_paoflow_opt

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

   module subroutine orb2site(orb, i_out, ia_out, n_atoms, max_orbital)
      integer, intent(in) :: orb, n_atoms, max_orbital
      integer, intent(out) :: i_out, ia_out

   end subroutine orb2site

   module subroutine site2orb(i_in, ia_in, orb_out, n_atoms, max_orbital)
      integer, intent(in) :: i_in, ia_in, n_atoms, max_orbital
      integer, intent(out) :: orb_out

   end subroutine site2orb

   module subroutine rotate_to_local_axis(this, m_loc)
      use math_mod, only: rotmag_loc
      class(hamiltonian), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: m_loc

      ! Local variables
      integer :: sdim
      ! Rotate Hamiltonian to local axis if wanted
   end subroutine rotate_to_local_axis

   module subroutine rotate_from_local_axis(this, m_loc)
      use math_mod, only: rotmag_loc
      class(hamiltonian), intent(inout) :: this
      real(rp), dimension(3), intent(in) :: m_loc

      ! Local variables
      integer :: sdim
      ! Rotate Hamiltonian to local axis if wanted
   end subroutine rotate_from_local_axis

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

   module subroutine calculate_hubbard_v_potential(this)
      class(hamiltonian), intent(inout) :: this

      integer :: itype, ia, nr, m, ja, jt, ispin, li, lj, ii, i0, i1, jdim
      real(rp) :: rmin, rcur, tol, occ_up, occ_dn, nn_shell_tol
      real(rp), dimension(3) :: dr
      real(rp), dimension(this%lattice%ntype, 4, 2) :: n_spin
      logical, save :: warned_proxy_once = .false.

   end subroutine calculate_hubbard_v_potential

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

   module subroutine export_rs_tb_metadata(this, filename)
      implicit none
      class(hamiltonian), intent(in) :: this
      character(len=*), intent(in) :: filename
   
      integer :: u, i, ios, n_atoms, max_orbital
   
   end subroutine export_rs_tb_metadata

   module subroutine export_rs_paoflow_legacy(this, filename, tol, include_lsham, transform_sph2cart)
      implicit none
      class(hamiltonian), intent(in) :: this
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: tol
      logical, intent(in) :: include_lsham, transform_sph2cart
   
      integer :: u, ios
   
   end subroutine export_rs_paoflow_legacy

   module subroutine export_rs_tb_hoppings(this, filename, tol, include_lsham, transform_sph2cart)
      implicit none
      class(hamiltonian), intent(in) :: this
      character(len=*), intent(in) :: filename
      real(rp), intent(in) :: tol
      logical, intent(in) :: include_lsham, transform_sph2cart
   
      integer :: u, ios
   
   end subroutine export_rs_tb_hoppings

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
