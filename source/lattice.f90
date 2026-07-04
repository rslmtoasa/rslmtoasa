!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Lattice
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
!> Module to handle structural properties
!------------------------------------------------------------------------------

module lattice_mod

#ifdef USE_MPI
   use mpi
#endif
   use mpi_mod
   use globals_mod
   use control_mod
   use string_mod
   use basis_mod, only: norb
   use math_mod
   use precision_mod, only: rp
   use symbolic_atom_mod, only: symbolic_atom, array_of_symbolic_atoms
   use namelist_generator_mod, only: namelist_generator
   use logger_mod, only: g_logger
   use timer_mod, only: g_timer
   use strux_lib, only: strux_compute, strux_lmto47_screening, strux_lmto47_autoalpha_screening, &
      strux_options, strux_result, STRUX_METHOD_LMTO47, STRUX_LMTO47_IALPHA_MANUAL, &
      STRUX_LMTO47_IALPHA_SIGMA, STRUX_LMTO47_IALPHA_FITD
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Module´s main structure
   type, public :: lattice
      !> Charge
      class(control), pointer :: control

      ! General variables

      !> Lattice parameter (in \f$Å\f$)
      !>
      !> Lattice parameter (in \f$Å\f$)
      real(rp) :: alat

      !> Sphere radius to cut cluster (in \f$Å\f$)
      !>
      !> Use \ref r2 = \ref ct\f$^2\f$. This radius (\ref ct) refers to the distance which includes all first neighbors,
      !> or all second nearest neighbors, etc. (use \ref ct and \ref r2 including 5 th neighs. to run newclu.x)
      !> Example: Pd fcc: let \f$R_1\f$ the distance to first neighbors, and \f$R_2\f$ the distance to second neighbors,
      !>
      !> \f$ R_1 = a \frac{\sqrt{2}}{2} = 2.7506 \, Å \f$,
      !>
      !> then,
      !>
      !> \f$ R_1^2 = 7.566 \, Å^2 \f$;
      !>
      !> and
      !>
      !> \f$R_2 = a = 3.89 \, Å\f$
      !>
      !> then,
      !>
      !> \f$R_2^2 = 15.13 \, Å^2\f$.
      !>
      !> Since \ref r2 \f$= 13\f$ is between \f$R_1^2\f$ and \f$R_2^2\f$ it will include all first neighs., but not second
      !> neighs.
      real(rp) :: r2

      !> TODO
      real(rp) :: celldm

      !> Wigner Seitz Radius (in \f$Å\f$)
      !>
      !> Wigner Seitz Radius (in \f$Å\f$)
      real(rp) :: wav

      !> Cell volume
      real(rp) :: vol

      !> Near neighbors cut radius (in \f$Å\f$)
      !>
      !> Near neighbors cut radius (in \f$Å\f$). See \ref r2 description for more details.
      real(rp), dimension(:), allocatable :: ct

      !> Auxiliar variables to set the correct ntype, nbulk and ntot
      integer :: nbulk_bulk

      !> Number of atoms of type bulk.
      !>
      !> Number of atoms of type bulk.
      !>
      !> Values:
      !>
      !> - \ref nbulk \f$ = 0\f$: for bulk
      !> - \ref nbulk \f$ = 1\f$: for impurities in bulk
      !> - \ref nbulk \f$ = 1\f$: for surface withou defects
      !> - \ref nbulk \f$ = 6\f$: for for a system with defects on surface, where this surface has been calculated with 5 layers plus bulk.
      integer :: nbulk

      !> \ref ntot \f$= 1\f$ for fcc and bcc, without relaxation; \ref ntot \f$= 2\f$ for hcp.
      !>
      !> - \ref ntot \f$= 1\f$ for fcc and bcc, without relaxation;
      !> - \ref ntot \f$= 2\f$ for hcp.
      !>
      !> TODO: poor description
      integer :: ntot

      !> Number of inequivalent atoms.
      !>
      !> Number of inequivalent atoms.
      !>
      !> Always equal to \ref nbulk + \ref nrec
      !>
      !> Usage:
      !>
      !> - \ref ntype \f$= 1\f$: bulk material
      !> - \ref ntype \f$= 2\f$: impurity embedded in bulk (single site calculation)
      !> - \ref ntype \f$= 3\f$: impurity (at1) ambedded in bulk (plus nearest neighs. at2)
      !> - \ref ntype \f$= 6\f$: free surface system with 5 layers (at1, at2, at3, at4, at5) (since \ref nbulk \f$= 1\f$)
      !> - \ref ntype \f$= 7\f$: adatom on surface (single site calculation), where the free surface were converged with 5 layers + bulk (since \ref nbulk \f$= 6\f$)
      integer :: ntype

      ! Atom type regarding the charge transfer - used for impurity mode
      integer, dimension(:), allocatable :: chargetrf_type
      !> Number of atoms to be considered in the recursion.
      !>
      !> Number of atoms to be considered in the recursion (for bulk \ref nrec\f$= 1\f$, for surface with 5 layers being calculated self-con. \ref nrec\f$= 5\f$, etc.)
      integer :: nrec

      !> TODO
      !> Clust size
      !>
      !> Clust size
      integer :: kk

      !> TODO
      !> Clust coordinates (expanded cluster)
      !>
      !> 'cr' holds the fully expanded cluster coordinates (cartesian coordinates)
      !> after the Bravais/cluster build. Shape is (3, kk) where 'kk' is the
      !> actual number of atoms in the constructed cluster. During construction
      !> a temporary local array with capacity (3, ndim) is often used and then
      !> moved into 'this%cr'. We shrink 'this%cr' to (3,kk) as soon as 'kk' is
      !> known to avoid holding a large buffer unnecessarily.
      real(rp), dimension(:, :), allocatable :: cr

      !> TODO
      !> Clust coordinates (primitive cell / basis)
      !>
      !> 'crd' stores the primitive cell (basis) coordinates (3, ntot). It is
      !> used as the starting point when expanding to the cluster (cr).
      real(rp), dimension(:, :), allocatable :: crd
      !> Clust atom number
      integer, dimension(:), allocatable :: ham_i
      !> TODO
      !> Atoms coordinates of the primitive cell
      !>
      !> Atoms coordinates of the primitive cell
      real(rp), dimension(:, :), allocatable :: primcell

      !> Number of shells/atoms to distribute charge
      !>
      !> Number of shells/atoms to distribute charge
      integer :: nbas

      !> Reduced \ref nbas
      !>
      !> Reduced \ref nbas
      integer :: reduced_nbas

      !>TODO
      !> Nmax
      integer :: nmax


      !> TODO
      !> Atom number for calculation.
      !>
      !> Vector of size \ref ntot.
      !>
      !> - \f$iu = 1\f$: for bulk material
      !> - Surface: this number IU is chosen in the clust file, looking for a site which can characterize the bulk layers, i.e. far from surface sites. The output on screen from buildsurf.x program gives this number.
      !> - Defects on surface: number given in newclu.x output (on screen).
      integer, dimension(:), allocatable :: iu

      !> TODO
      !> Atom number for calculation.
      !>
      !> Vector of size \ref ntot.
      !>
      !> - \f$ib = 1\f$: for bulk material
      !> - Surface: this number IB is chosen in the clust file, looking for a site which can characterize the bulk layers, i.e. far from surface sites. The output on screen from buildsurf.x program gives this number.
      !> - Defects on surface: number given in newclu.x output (on screen).
      integer, dimension(:), allocatable :: ib

      !> Refers to the atoms that will be calculated self-consistently, or to represent all equivalent atoms in the same neighboring shell.
      !>
      !> Refers to the atoms (in the clust file) that will be calculated self-consistently, or to represent all equivalent atoms in the same neighboring shell.
      !>
      !> - Those numbers are given by newclu.x output (on screen), for defects on surface and by buildsurf.x output (atomch.d file) for surface systems.
      integer, dimension(:), allocatable :: irec

      !> List containing the number of each atom inside the clust file
      integer, dimension(:), allocatable :: atlist

      !> Neighbouring map for each atom type
      integer, dimension(:, :), allocatable :: nn

      !> Maximum number of neighbours found in current structure map.
      integer :: nn_max

      !> Structure constant
      complex(rp), dimension(:, :, :, :), allocatable :: sbar
      !> Structure constant derivative
      complex(rp), dimension(:, :, :, :), allocatable :: sdot
      !> Vectors in the structure constant
      real(rp), dimension(:, :), allocatable :: sbarvec
      !> Screening constants for the local strux solve cluster
      real(rp), dimension(:, :, :), allocatable :: alpha
      !> Screening constant derivatives for the local strux solve cluster
      real(rp), dimension(:, :, :), allocatable :: alpha_dot
      !> Structure-constant backend selector
      character(len=16) :: strux_backend
      !> Screening mode for the strux backend
      character(len=16) :: screening
      !> Request derivative structure constants from strux
      logical :: strux_want_sdot
      !> Multiplier for the solve-cluster radius squared used by strux
      real(rp) :: strux_solve_scale
      !> Global l-resolved screening constants for manual mode
      real(rp), dimension(:), allocatable :: screening_alpha
      !> Global l-resolved normalized hard-core radii for sigma mode
      real(rp), dimension(:), allocatable :: screening_sigma
      ! Variables to build clust for bulk calculation

      !> TODO
      !> Cut radius
      !>
      !> Cut radius
      real(rp) :: rc

      !> Primitive vectors in units of lattice parameter \ref alat
      !>
      !> Primitive vectors in units of lattice parameter \ref alat
      real(rp), dimension(3, 3) :: a
      !> Inverse lattice matrix in Cartesian coordinates (if available).
      real(rp), dimension(3, 3) :: a_cart_inv
      !> Flag indicating whether `a_cart_inv` is valid.
      logical :: a_cart_inv_ready
      !> Variables to handle periodic boundary conditions
      !> 
      !> Variables to handle periodic boundary conditions
      logical :: pbc
      !> Flag to reorder the bulk cluster along a Morton (Z-order) space-filling
      !> curve to improve data locality in the recursion/Chebyshev kernels.
      !> Default .false. (preserves legacy lattice-scan ordering).
      logical :: morton_sfc
      logical :: b1, b2, b3
      integer :: n1, n2, n3  
      !> TODO
      !> TBD
      integer, dimension(:), allocatable :: izp

      !> TODO
      integer, dimension(:), allocatable :: iz

      !> TODO
      integer, dimension(:), allocatable :: num

      !> TODO
      integer, dimension(:), allocatable :: no

      !> TODO
      !> Paramemter that determines the clust size before cut. Should be large
      !enough, rarely changed. Default (ndim = 9900000, npe = 49)
      integer :: ndim, npe

      !> TODO
      !> Crystal symmetry. Options are ´bcc´, ´fcc´, ´hcp´ and ´nsy´
      !>
      !> Crystal symmetry. Options are ´bcc´, ´fcc´, ´hcp´ and ´nsy´
      character(len=4) :: crystal_sym

      ! Variables to build clust for surface calculation
      !> TODO
      !> Surface parameters
      real(rp) :: zmin, zmax, zstep

      !> TODO
      !> Layer coordinate
      real(rp), dimension(:), allocatable :: z

      !> TODO
      !> Surface symmetry. Options are ´111´, ´110´ and ´001´.
      !>
      !> Surface symmetry. Options are ´111´, ´110´ and ´001´.
      character(len=10) :: surftype

      !> TODO
      !> Number of layers
      !>
      !> Number of layers
      integer :: nlay

      !> TODO
      !> Surface indexes
      integer :: dx, dy, dz, dw
      !> Number of atoms per layer. Surface calculation only
      integer, dimension(:), allocatable :: natoms_layer
      !> TODO
      !> TDB
      integer, dimension(:), allocatable :: izpsurf, izsurf, nosurf

      !> TODO
      !> Clust variables from surface
      real(rp), dimension(:, :), allocatable :: crsurf

      ! Variables to build clust for impurity calculation
      !> New coordinates
      real(rp), dimension(:, :), allocatable :: inclu

      !> TDB
      integer, dimension(:), allocatable :: izpo

      !> Impurity number
      integer :: nclu

      !> Clust variables from impurity
      real(rp), dimension(:, :), allocatable :: acr

      !> TODO
      integer, dimension(:), allocatable :: reduced_acr

      !> Pair of atoms to calculate the exchange interactions
      integer, dimension(:, :), allocatable :: ijpair
      !> Number of pairs
      integer :: njij

      !> Trio of atoms to calculate the spin-lattice interactions
      real(rp), dimension(:, :), allocatable :: ijktrio
      !> Number of trios i,j,k
      integer :: njijk

      type(symbolic_atom), dimension(:), allocatable :: symbolic_atoms
   contains
      procedure :: build_from_file
      procedure :: build_from_lattice
      procedure :: restore_to_default
      procedure :: bravais
      procedure :: morton_reorder_bulk
      procedure :: build_data
      procedure :: build_clusup
      procedure :: build_surf
      procedure :: build_surf_full
      procedure :: newclu
      procedure :: structb
      procedure :: atomlist
      procedure :: check_atoms_in_volume
      procedure :: find_unique_struct
      procedure :: identify_unique_atoms
      procedure :: check_within_volume
      procedure :: f_wrap_coord_diff
      procedure :: remd
      procedure :: nncal
      procedure, private :: dbar1
      procedure, private :: structb_strux
      procedure, private :: init_strux_storage
      procedure, private :: build_strux_inputs
      procedure, private :: strux_mode
      procedure, private :: default_screening_alpha
      procedure, private :: build_hcr
      procedure, private :: build_rmt
      procedure, private :: load_symbolic_atoms_if_needed
      procedure :: clusba
      procedure :: calculate_nbas
      procedure :: print_state
      procedure :: print_state_full
      procedure :: print_state_formatted
      procedure, private :: check_all
      final :: destructor
   end type lattice

   interface lattice
      procedure :: constructor
   end interface lattice


   interface
   module function constructor(control_obj) result(obj)
      type(lattice) :: obj
      type(control), target, intent(in) :: control_obj

   end function constructor

   module subroutine destructor(this)
      type(lattice) :: this
   end subroutine destructor

   module subroutine build_from_file(this, fname)
      class(lattice), intent(inout) :: this
      character(len=*), intent(in), optional :: fname
      ! variables associated with the reading processes
      integer :: iostatus, funit, i
      character(len=sl) :: fname_
      real(rp), allocatable :: ct_first(:), screening_alpha_first(:), screening_sigma_first(:)

   end subroutine build_from_file

   module subroutine build_from_lattice(this)
      class(lattice), intent(inout) :: this
      integer :: nbulk_bulk, ntot, nbas, nrec, funit, iostatus, nsite_guess
      real(rp) :: r2, strux_solve_scale
      real(rp), dimension(3, 3) :: a
      real(rp), dimension(:), allocatable :: ct, screening_alpha, screening_sigma
      real(rp), dimension(:), allocatable :: ct_first, screening_alpha_first, screening_sigma_first
      integer, dimension(:), allocatable :: izp, no, iu, ib, irec
      real(rp), dimension(:, :), allocatable :: crd
      character(len=16) :: strux_backend, screening
      logical :: strux_want_sdot
   end subroutine build_from_lattice

   module subroutine build_data(this)
      class(lattice), intent(inout) :: this
      !> Local variables
      real(rp), dimension(3, 3) :: a
      integer :: i, j

   end subroutine build_data

   module subroutine restore_to_default(this, full)
      class(lattice) :: this
      logical, intent(in), optional :: full

   end subroutine restore_to_default

   module subroutine bravais(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      real(rp) :: rc, rs, lc, lcx, lcy, lcz
      integer, dimension(:), allocatable :: iz, num
      real(rp), dimension(:, :), allocatable :: cr, crbravais
   real(rp), allocatable :: tmp(:, :)
      integer :: npe, ndim, nx, ny, nz, npr, l, n, i, nl, k, kk
      logical :: isopen
      integer :: iostatus

   end subroutine bravais

   module subroutine morton_reorder_bulk(this)
      class(lattice), intent(inout) :: this
      integer :: kk, i, n
      integer, allocatable :: perm(:)
      integer(8), allocatable :: key(:)
      real(rp) :: cmin(3), cmax(3), span(3), inv_span(3)
      integer, parameter :: nbits = 21      ! 21 bits/axis -> 63-bit Morton key
      integer, parameter :: gmax = 2**nbits - 1
      real(rp), allocatable :: cr_new(:, :)
      integer, allocatable :: iz_new(:), num_new(:)
      integer :: q(3)

   end subroutine morton_reorder_bulk

   module pure function morton_encode3(x, y, z, nbits) result(code)
      integer, intent(in) :: x, y, z, nbits
      integer(8) :: code
      integer :: b
   end function morton_encode3

   module subroutine sort_by_key(key, perm, n)
      integer(8), intent(inout) :: key(:)
      integer, intent(inout) :: perm(:)
      integer, intent(in) :: n
      integer :: i, start, bottom
      integer(8) :: tkey
      integer :: tperm

      ! Build heap then sift down (standard heapsort, O(n log n)).
   end subroutine sort_by_key

   module subroutine build_clusup(this)
      class(lattice), intent(inout) :: this
      character(len=10) :: surftype
      ! Local variables
      real(rp), dimension(:), allocatable :: z
      integer, dimension(:), allocatable :: izpsurf
      real(rp) :: ds, ds2, new, one
      integer :: i, j, n

   end subroutine build_clusup

   module subroutine build_surf_full(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, natoms, nsurf, currentType, newType, newCrystalType
      integer :: nTypesTotal, nUnique, atomIdx, nTypesInLayer, ichoice, ichoicetype
      real(rp), dimension(:, :), allocatable :: crsurf
      real(rp), dimension(:), allocatable :: crh, crhd, z
      integer, dimension(:), allocatable :: atomType, crystalType, typesurf, crystalsurf, uniqueTypes, ichoicen, ichoicetypen
      integer, dimension(:), allocatable :: nTypesForCurrentLayer
      real(rp) :: dx, dy, dz, new, one, ds, ds2, disi, disi_min
      real(rp) :: zstep, zmin, zmax
      integer :: n, atomCount, maxType, nlay
      logical :: isUnique, isopen
      character(20) :: header
      ! Variables
      real(rp) :: rotated_cr(3, this%kk)
      real(rp) :: rotation_matrix(3, 3)
      real(rp) :: axis(3), theta
      real(rp) :: norm

      ! Initial definitions
   end subroutine build_surf_full

   module subroutine build_surf(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      integer :: i, j, n, k, kk, ichoice, icont, ichoicetype
      real(rp), dimension(this%kk) :: crh, crhd
      real(rp) :: disf, disi_min, disi
      real(rp), dimension(:, :), allocatable :: crsurf
      integer, dimension(:), allocatable :: izp, no, ichoicen, ichoicetypen
      logical :: isopen

   end subroutine build_surf

   module subroutine newclu(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      integer :: ndi, nnmx, ncnt, nmax
      logical :: isopen
      integer :: i, j, kk, k, ntypecount, ireccount, inclucheck
      integer, dimension(:), allocatable ::  ibulk
      integer, dimension(:), allocatable :: izpo, izp, no, nnmax, izimp, noimp
      integer, dimension(:, :), allocatable :: nn, nn2
      real(rp) :: nnscale
      real(rp), dimension(:, :), allocatable :: acr, crd, crimp
      real(rp), dimension(:), allocatable :: ctnew

      ! Set clust variables
   end subroutine newclu

   module subroutine structb(this, do_str)
      class(lattice), intent(inout) :: this
      logical, intent(in) :: do_str
      ! Local variables
      integer :: i, ia, nr, ii, j, nm, np, nlim, nomx, ncut, kk, nnmx
      integer :: sbar_dim, nm_store, nt_tmp
      integer, dimension(:, :), allocatable :: nn
      integer, dimension(:), allocatable :: idnn
      logical :: do_str_
      real(rp), dimension(:, :, :), allocatable :: set
      real(rp), dimension(3) :: ret
      real(rp) :: t_structb_start, t_cluster_ready, t_nnmap_start, t_nnmap_end
      real(rp) :: t_remd_start, t_remd_end, t_nm_store_start, t_nm_store_end
      real(rp) :: t_outmap_start, t_outmap_end, t_str_stage_end

   end subroutine structb

   module subroutine init_strux_storage(this, sbar_dim, nm)
      class(lattice), intent(inout) :: this
      integer, intent(in) :: sbar_dim, nm

   end subroutine init_strux_storage

   module pure function default_screening_alpha(this, nl) result(alpha_default)
      class(lattice), intent(in) :: this
      integer, intent(in) :: nl
      real(rp) :: alpha_default(nl)
      real(rp), parameter :: default_values(4) = [0.3485_rp, 0.0530_rp, 0.0107_rp, 0.00674_rp]
      integer :: i

   end function default_screening_alpha

   module integer function strux_mode(this)
      class(lattice), intent(in) :: this
      character(len=:), allocatable :: screening_mode

   end function strux_mode

   module subroutine load_symbolic_atoms_if_needed(this)
      class(lattice), intent(inout) :: this
      
      ! 1. Declare the temporary array. 
      type(symbolic_atom), allocatable :: temp_atoms(:)

   end subroutine load_symbolic_atoms_if_needed

   module subroutine build_rmt(this, nspec, species_labels, rmt)
      class(lattice), intent(inout) :: this
      integer, intent(in) :: nspec
      integer, intent(in) :: species_labels(nspec)
      real(rp), intent(out) :: rmt(nspec)
      integer :: is, label

   end subroutine build_rmt

   module subroutine build_hcr(this, nl, nspec, rmt, hcr)
      class(lattice), intent(in) :: this
      integer, intent(in) :: nl, nspec
      real(rp), intent(in) :: rmt(nspec)
      real(rp), intent(out) :: hcr(nl, nspec)
      integer :: is, l

   end subroutine build_hcr

   module subroutine build_strux_inputs(this, nspec, nl, species_labels, alpha_in, hcr, rmt)
      class(lattice), intent(inout) :: this
      integer, intent(in) :: nspec, nl
      integer, intent(in) :: species_labels(nspec)
      real(rp), intent(out) :: alpha_in(0:nl - 1, nspec)
      real(rp), intent(out) :: hcr(nl, nspec)
      real(rp), intent(out) :: rmt(nspec)
      integer :: is
      real(rp) :: alpha_global(nl)

   end subroutine build_strux_inputs

   module subroutine structb_strux(this)
      class(lattice), intent(inout) :: this

      integer, parameter :: max_orb = 16
      real(rp), parameter :: match_tol = 1.0e-5_rp
      integer :: ii, ia, ib, ja, jb, m, is, js, nspec, nl, nl2, sbar_dim, pair_idx, nt_store
      integer :: nbas, nttab
      integer :: label, species_idx
      integer, allocatable :: ips(:), lmxb(:), orb_map(:), species_labels(:)
      real(rp) :: pair_cutoff, solve_cutoff, pair_cutoff_bohr, solve_cutoff_bohr, alat_bohr, wav_bohr
      real(rp) :: t_total_start, t_total_end, t_map_start, t_map_end, t_inputs_start, t_inputs_end, t_store_start, t_store_end, t_remap_end
      real(rp) :: alpha_debug(4)
      character(len=16) :: effective_screening
      real(rp) :: vec_target(3)
      real(rp), allocatable :: pos(:,:), cralat(:,:), rmt(:), alpha_in(:,:), hcr(:,:), alpha_site(:,:), &
         adot_site(:,:), alpha_l_out(:,:), tral(:,:,:), trad(:,:,:)
      type(strux_options) :: opts
      type(strux_result) :: result

   end subroutine structb_strux

   module subroutine write_strux_block(this, nl2, m, ii, jclus)
      class(lattice), intent(in) :: this
      integer, intent(in) :: nl2, m, ii, jclus
      integer :: is, js

      ! write (*, '(" SBAR neighbor center=",i5," slot=",i5," iclus=",i5," vec=",3f12.6)') ii, m, iclus, &
      !    this%sbarvec(1, m), this%sbarvec(2, m), this%sbarvec(3, m)
   end subroutine write_strux_block

   module subroutine write_neighbor_vector_dump(iunit, sbarvec, nt)
      integer, intent(in) :: iunit, nt
      real(rp), intent(in) :: sbarvec(:, :)
      integer :: i

   end subroutine write_neighbor_vector_dump

   module integer function representative_atom_index(this, ii)
      class(lattice), intent(in) :: this
      integer, intent(in) :: ii

   end function representative_atom_index

   module integer function primitive_basis_label(this, ia)
      class(lattice), intent(in) :: this
      integer, intent(in) :: ia

   end function primitive_basis_label

   module integer function find_pair_by_vector(nttab, iax, plat, pos, alat, ib, jb, vec_target, tol)
      integer, intent(in) :: nttab, ib, jb
      integer, intent(in) :: iax(:,:)
      real(rp), intent(in) :: plat(3, 3), pos(:, :), alat, vec_target(3), tol
      integer :: i, n1, n2, n3
      real(rp) :: vec(3)
      real(rp) :: dmax, best_d
      integer :: best_i

   end function find_pair_by_vector

   module integer function find_neighbor_atom_by_vector(this, ia, vec_target, cralat, tol)
      class(lattice), intent(in) :: this
      integer, intent(in) :: ia
      real(rp), intent(in) :: vec_target(:), cralat(:, :)
      real(rp), intent(in) :: tol
      integer :: ja
      real(rp) :: dv(3)

   end function find_neighbor_atom_by_vector

   module subroutine build_orbital_map(norb, orb_map)
      integer, intent(in) :: norb
      integer, intent(out) :: orb_map(16)
      integer :: i
      integer, parameter :: full_map(16) = [ &
         1, 4, 2, 3, 5, 6, 8, 9, 7, 13, 14, 12, 15, 11, 16, 10 ]

   end subroutine build_orbital_map

   module subroutine atomlist(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      integer :: i, j, itype
      real(rp) :: mom_tmp(3)

   end subroutine atomlist

   module subroutine check_atoms_in_volume(this, cr, num, num_atoms, central_atom, a1, a2, a3, plane_constant, atoms_in_volume, atom_count)
       class(lattice), intent(inout) :: this
       real(rp), intent(in) :: cr(3, num_atoms), a1(3), a2(3), a3(3)
       integer, intent(in) :: num(num_atoms), num_atoms, central_atom
       real(rp), intent(in) :: plane_constant
       integer, allocatable, intent(out) :: atoms_in_volume(:)
       integer, intent(out) :: atom_count
       real(rp) :: relative_pos(3)
       logical :: inside
       integer :: i
   
       ! Initialize array for atoms inside the primitive cell volume
   end subroutine check_atoms_in_volume

   module subroutine check_within_volume(this, relative_pos, a1, a2, a3, inside)
       class(lattice), intent(inout) :: this
       real(rp), intent(in) :: relative_pos(3), a1(3), a2(3), a3(3)
       logical, intent(out) :: inside
       real(rp) :: dot11, dot12, dot13, dot22, dot23, dot33
       real(rp) :: dot1r, dot2r, dot3r, inv_denom, u, v, w
   
       ! Calculate dot products between the vectors
   end subroutine check_within_volume

   module subroutine find_unique_struct(this, num, num_atoms, unique_nums)
       class(lattice), intent(inout) :: this
       integer, intent(in) :: num(:), num_atoms
       integer, allocatable, intent(out) :: unique_nums(:)
       integer, allocatable :: temp_nums(:)
       integer :: i, j, num_unique
       logical :: found
   
       ! Allocate temporary array for unique numbers
   end subroutine find_unique_struct

   module subroutine identify_unique_atoms(this, cr, num_atoms, atoms_in_volume, atom_count, a1, a2, a3, unique_atoms, unique_atom_count)
       class(lattice), intent(inout) :: this
       real(rp), intent(in) :: cr(3, num_atoms), a1(3), a2(3), a3(3)
       integer, intent(in) :: num_atoms, atoms_in_volume(:), atom_count
       integer, allocatable, intent(out) :: unique_atoms(:)
       integer, intent(out) :: unique_atom_count
       logical :: found, is_transformed
       real(rp) :: trans_atom(3), delta(3)
       integer :: i, j, k, n, m, p
       integer, allocatable :: temp_unique_atoms(:)
   
       ! Initialize temporary array for unique atoms list
   end subroutine identify_unique_atoms

   module subroutine dbar1(this, ia, r2, wav, crd, nat, ndi, np, nr, ii)
      implicit none
      class(lattice), intent(inout) :: this
      ! Inputs
      integer, intent(in) :: ia, nat, ndi, np
      integer, intent(in) :: nr, ii
      real(rp), intent(in) :: r2, wav
      real(rp), dimension(3, ndi), intent(in) :: crd
      ! Local Scalars
      integer :: i, j, k, m, na, nrl, nt, jclus_dbg
      real(rp), dimension(:), allocatable :: bet, wk
      real(rp), dimension(:), allocatable :: a
      real(rp), dimension(:, :), allocatable :: cr
      real(rp), dimension(:, :), allocatable :: s
      real(rp), dimension(:, :, :), allocatable :: sbar
      real(rp), dimension(:, :), allocatable :: sbarvec
      !
      ! External Calls
      !external CLUSBA, MICHA

   end subroutine dbar1

   module subroutine clusba(this, r2, crd, ia, nat, ndi, n, sbarvec_out)
      implicit none
      class(lattice), intent(inout) :: this
      ! Inputs
      integer, intent(in) :: ia, nat, ndi
      real(rp), intent(in) :: r2
      real(rp), dimension(3, ndi), intent(in) :: crd
      ! Output
      integer, intent(inout) :: n
      real(rp), dimension(:, :), intent(inout), optional :: sbarvec_out
      ! Local variables
      integer :: i, ii, k, nn
      real(rp) :: s1
      real(rp), dimension(3) :: dum

   end subroutine clusba

   module subroutine micha(rws, r, nr, nlm, nrl, na, sbar, a, wk, bet, s, iclus, r2)
      implicit none
      ! Inputs              .
      integer, intent(in) :: iclus, na, nlm, nr, nrl
      real(rp), intent(in) :: rws, r2
      real(rp), dimension(3, nr), intent(in) :: r
      ! Outputs
      real(rp), dimension(na), intent(inout) :: a
      real(rp), dimension(nrl), intent(inout) :: bet, wk
      real(rp), dimension(nrl, nrl), intent(inout) :: s
      real(rp), dimension(nlm, nlm, nr), intent(inout) :: sbar
      ! Local variables
      real(rp) :: fak, pi
      real(rp), dimension(4) :: q
      ! External Calls
      !external SHLDCH, STREZE
      ! Intrinsic Functions
      intrinsic ATAN

   end subroutine micha

   module subroutine STREZE(w, r, nr, s, nrl, nlm)
      implicit none
      ! Input
      integer, intent(in) :: nlm, nr, nrl
      real(rp), intent(in) :: w
      real(rp), dimension(3, nr), intent(in) :: r
      ! Output
      real(rp), dimension(nrl, nrl), intent(inout) :: s
      ! Local variables
      integer :: ilm, ir, irl0, jlm, jr, jrl0
      real(rp) :: rr, w1
      real(rp), dimension(3) :: dr
      real(rp), dimension(16, 16) :: s0
      ! External calls
      !external CANSO
      ! Intrinsic Functions
      intrinsic SQRT
   end subroutine streze

   module subroutine shldch(r, nr, nlm, nrl, s, a, na, q, bet, wk, sbar, iclus, r2)
      implicit none
      !parameter for cutoff of sbar construction
      integer, parameter :: ncut = 9
      ! Input
      integer, intent(in) :: iclus, na, nlm, nr, nrl
      real(rp), intent(in) ::r2
      real(rp), dimension(4), intent(in) :: q
      real(rp), dimension(3, nr), intent(in) :: r
      ! Output
      real(rp), dimension(na), intent(inout) :: a
      real(rp), dimension(nrl), intent(inout) :: bet, wk
      real(rp), dimension(nrl, nrl), intent(inout) :: s
      real(rp), dimension(nlm, nlm, nr), intent(inout) :: sbar
      ! Local variables
      integer :: i, ia, ilm, ir, irl, irl0, isb, j, jlm, jsb, l, l2, lmax, m, ndef, hitc, info
      real(rp), dimension(:, :), allocatable :: s_temp
      ! External Calls
      !external chlr2f, chlr2s
      ! External Functions  .
      !integer, external :: LL
   end subroutine shldch

   module subroutine canso(w, dr, sc)
      implicit none
      ! Input
      real(rp), intent(in) :: w
      real(rp), dimension(3), intent(in) :: dr
      ! Output
      real(rp), dimension(16, 16), intent(out) :: sc
      ! Local variables
      integer :: i, j, l, ll
      real(rp) :: el, el2, elem, elen, em, em2, emen, en, en2, r1, r2, r3, rr, s2, s3, s4, s5, s6, s7, &
                  sbyr, sq3, sq5, sq7
      integer, dimension(16) :: ip
      real(rp), dimension(16, 16) :: s
      ! Intrinsic Functions
      intrinsic SQRT
      !.. Data Declarations ..
      ! original and correct
      !     S=1, X=2, Y=3, Z=4, XY=5, YZ=6, ZX=7, X**2-Y**2=8, 3Z*Z-R*R=9
      !     f-orbitals (10-16): fz3, fxz2, fyz2, fz(x2-y2), fxyz, fx3, fy3
   end subroutine canso

   module subroutine chlr2f(c, na, w, n, ndef)
      implicit none
      ! Input
      integer, intent(in) :: n, na
      ! Output
      integer, intent(out) :: ndef
      real(rp), dimension(n), intent(inout) :: w
      real(rp), dimension(na), intent(inout) :: c
      ! Local variables
      integer :: i, ic0, j, jc0, k
      real(rp) :: csum
      ! Intrinsic Functions
      intrinsic SQRT

   end subroutine chlr2f

   module subroutine chlr2s(c, na, v, n, m)
      implicit none
      ! Input
      integer, intent(in) :: m, n, na
      real(rp), dimension(na), intent(in) :: c
      ! Output
      real(rp), dimension(n, m), intent(inout) :: v
      ! Local Variables
      integer :: i, ic0, ind1, ind2, k, mm
      real(rp) :: csum

   end subroutine chlr2s

   module function LL(ilm)
      implicit none
      ! Input
      integer, intent(in) :: ilm
      ! Function Declaration
      integer :: LL
      ! Local variables
      integer, dimension(100) :: lla
      ! Data Declarations
   end function LL

   module subroutine remd(this, crd, no, iu, nn, nat, ntot, nomx, ndi, nnmx, set, idnn, ret)
      implicit none
      ! Inputs
      class(lattice), intent(inout) :: this
      integer, intent(in) :: nat, ndi, nnmx, nomx, ntot
      integer, dimension(ndi), intent(in) :: no
      real(rp), dimension(3, ndi), intent(in) :: crd
      integer, dimension(nomx), intent(in) :: iu
      ! Output
      integer, dimension(nnmx), intent(inout) :: idnn
      integer, dimension(ndi, nnmx), intent(inout) :: nn
      real(rp), dimension(3), intent(out) :: ret
      real(rp), dimension(3, nomx, nnmx), intent(inout) :: set
      ! Local variables
      integer :: i, ii, iii, imax, inn, ino, j, jj, jnn, jsz, k, la, lk, lm, m, n
      real(rp) :: a1, a2, a3, aaa, eps
      !-BUILDS VECTORS SET(3, NOMX, NNMX) CONNECTING NEIGHBORS OF EACH TYPE NO-
   end subroutine remd

   module subroutine outmap(IM, IZP, NN, NO, ND, NM, NTOT)
      implicit none
      ! Input
      integer, intent(in) :: IM, ND, NM, NTOT
      integer, dimension(ND), intent(in) :: NO
      integer, dimension(ND) :: IZP
      integer, dimension(ND, NM), intent(in) :: NN
      ! Local variables
      integer :: I, ID, IDM, J, K
      ! Intrinsic functions
      intrinsic MIN

   end subroutine outmap

   module integer function mapa(I, J, R2, DD, CT)
      implicit none
      ! Input
      integer, intent(in) :: I, J
      real(rp), intent(in) :: R2
      real(rp) :: DD
      real(rp), dimension(50), intent(in) :: CT
      ! Local variables
      real(rp) :: CTM, CTSM

   end function mapa

   module subroutine f_wrap_coord_diff(this,Natom,coord,i_atom,j_atom,cdiff)
      implicit none
      class(lattice), intent(inout) :: this
      integer, intent(in) :: Natom
      real(rp), dimension(3,Natom), intent(in) :: coord
      integer, intent(in) :: i_atom
      integer, intent(in) :: j_atom
      real(rp), dimension(3), intent(out) :: cdiff
      !
      real(rp), dimension(3) :: odiff, oshift, mdiff
      integer :: x,y,z
      integer :: xmin,xmax,ymin,ymax,zmin,zmax
      !
   end subroutine f_wrap_coord_diff

   module subroutine nncal(this,ct, crd, ndim, nat, izp, nn, nd, nm, ngbr, ntot)
      implicit none
      class(lattice), intent(inout) :: this
      ! Input
      integer, intent(in) :: NAT, ND, NDIM, NTOT
      integer, dimension(NAT), intent(in) :: IZP
      real(rp), dimension(NDIM, NAT), intent(in) :: CRD
      ! Output
      integer, intent(inout) :: NM
      integer, dimension(ND, NM), intent(inout) :: NN
      real(rp), dimension(50), intent(inout) :: CT
      ! External function
      integer, external :: NGBR
      ! Intrinsic function
      intrinsic ABS, MAX, MIN, FLOOR, NINT
      ! Local variables
      integer :: I, IADD, ID, II, IIP, ILJ, J, JJP, L, NNMAX
      integer :: NX, NY, NZ, NBIN, BX, BY, BZ, DBX, DBY, DBZ, BIN_ID
      integer :: CAP, CAND_COUNT, K, IX, IY, IZ, RX, RY, RZ
      real(rp) :: R2, RCUT, RCUT2, DETC
      real(rp), dimension(3) :: DDUM
      real(rp), dimension(3) :: MINC, MAXC, SPAN, BINW, DS
      real(rp), dimension(3, 3) :: CELL, CELL_INV
      real(rp), dimension(3, 2) :: CROSS_TMP
      real(rp), dimension(NM) :: DUM
      real(rp), allocatable :: FRAC(:, :)
      integer, allocatable :: HEAD(:), NEXT_ATOM(:), BIN_X(:), BIN_Y(:), BIN_Z(:), CANDIDATES(:)

   end subroutine

   module subroutine leia(alat, nndim, cr, iz, n, ip)
      implicit none
      ! Input
      integer, intent(in) :: ip, nndim
      real(rp), intent(in) :: alat
      ! Output
      integer, dimension(nndim + 10), intent(inout) :: iz, n
      real(rp), dimension(3, nndim + 10), intent(inout) :: cr
      ! Local variables
      integer :: i, j

   end subroutine

   module subroutine bubble(nl, ndim, m, nd, nt)
      implicit none
      ! Inputs
      integer, intent(in) :: ndim, nl, nd, nt
      ! Output
      real(rp), dimension(ndim, nd), intent(inout) :: m
      ! Local variables
      integer :: ind, inic, j, k
      real(rp) :: fim, z

   end subroutine bubble

   module subroutine cut(i, l, ndim, crd, cr, izp, iz, num, no, rs, ii)
      ! Inputs
      real(rp), intent(in) :: rs
      integer, intent(in) :: i, l, ndim
      real(rp), dimension(3, ndim), intent(in) :: crd
      integer, dimension(ndim), intent(in) :: izp, no
      ! Output
      integer, intent(out) :: ii
      real(rp), dimension(3, ndim), intent(out) :: cr
      integer, dimension(ndim), intent(out) :: iz, num
      ! Local variables
      real(rp) :: r2
      real(rp), dimension(3) :: dum
      integer :: na, j

   end subroutine cut

   module subroutine check_all(this)
      implicit none
      class(lattice) :: this

   end subroutine check_all

   module subroutine calculate_nbas(this)
      implicit none
      class(lattice) :: this
      integer :: size_iz
      integer, dimension(:) :: atype(this%nbas), amount(this%nbas)
      integer :: i, j

   end subroutine calculate_nbas

   module subroutine print_state_full(this, unit, file)
      class(lattice), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

   end subroutine print_state_full

   module subroutine print_state(this, unit, file)
      implicit none
      class(lattice), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

   end subroutine print_state

   module subroutine print_state_formatted(this, unit, file)
      class(lattice), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file

      type(namelist_generator) :: nml
      integer :: i

   end subroutine print_state_formatted

   end interface

end module lattice_mod
