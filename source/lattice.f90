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

   real(rp), parameter :: default_screening_alpha_values(4) = [0.3485_rp, 0.0530_rp, 0.0107_rp, 0.00674_rp]

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

      !> B7.5 (calctype='L'): number of FROZEN atomic layers of region A (the
      !> low-z semi-infinite reference) and of region B (the high-z one). The
      !> active zone is whatever lies between them, so the layer stack is
      !>     nlay_a  frozen A | nlay  active | nlay_b  frozen B
      !> and the total layer count is nlay_a + nlay + nlay_b.
      !>
      !> These are the LATTICE-side copies, read from the &lattice namelist.
      !> `charge%nlay_a`/`nlay_b` (from &charge) are the row counts the region
      !> registry partitions on and are deliberately kept separate: they count
      !> SITES (registry rows), while these count LAYERS. For a one-atom-per-
      !> layer cluster the two coincide; in general they do not, and conflating
      !> them would silently mis-place the A/active boundary.
      integer :: nlay_a, nlay_b

      !> B7.6 (calctype='L'): what region B physically IS. Two values:
      !>
      !>   'metal'  (default) region B is a second metallic reference, loaded
      !>            from its own converged parameter set through &atoms
      !>            label(:), exactly as region A is. This is the A | B
      !>            geometry and the pre-B7.6 behaviour.
      !>   'vacuum' region B is semi-infinite vacuum. Its frozen parameters
      !>            are NOT read from a label; they are generated per run by
      !>            `vacuum_lead` from the empty-sphere radius and the vacuum
      !>            level (B7 §1.6), and regenerated whenever the alignment
      !>            solver moves that level. This is the A | vacuum geometry.
      !>
      !> The distinction is confined to where region B's frozen parameters
      !> COME FROM and to the region kind the registry records. Layer binning,
      !> structure constants and the recursion are unaware of it -- B7 §0.2 is
      !> explicit that the Green-function machinery is region-agnostic and
      !> that a vacuum Green function is an ordinary one whose frozen
      !> parameters happen to describe an empty lattice.
      character(len=16) :: region_b_kind

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
      procedure :: build_interface_full
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
   !> @brief Construct a lattice object from an initialized control object.
   !> @details Wires the lattice to the run control state, restores default
   !>          geometry/structure-constant storage, and reads lattice input.
   !> @param[in] control_obj Control object that owns the input-file name and run options.
   !> @return Initialized lattice object.
   module function constructor(control_obj) result(obj)
      type(lattice) :: obj
      type(control), target, intent(in) :: control_obj

   end function constructor

   !> @brief Finalize a lattice object.
   !> @details Releases cluster, neighbor, surface, impurity, screening, and
   !>          structure-constant arrays owned by the object.
   !> @param[inout] this Lattice object being finalized.
   module subroutine destructor(this)
      type(lattice) :: this
   end subroutine destructor

   !> @brief Read the &lattice namelist and build derived lattice data.
   !> @details Parses bulk/surface/impurity geometry, cluster cutoffs, periodic
   !>          boundary flags, and strux-screening options, then prepares
   !>          primitive-cell and cluster state for later stages.
   !> @param[inout] this Lattice object to populate.
   !> @param[in] fname Optional input file; defaults to this%control%fname.
   !> @note This is an input boundary and may raise fatal diagnostics for invalid options.
   module subroutine build_from_file(this, fname)
      class(lattice), intent(inout) :: this
      character(len=*), intent(in), optional :: fname
      ! variables associated with the reading processes
      integer :: iostatus, funit, i
      character(len=sl) :: fname_
      real(rp), allocatable :: ct_first(:), screening_alpha_first(:), screening_sigma_first(:)

   end subroutine build_from_file

   !> @brief Rebuild lattice state from already installed lattice members.
   !> @details Reuses current geometry and namelist-derived settings to run the
   !>          same cluster, primitive-cell, and structure setup normally reached
   !>          from build_from_file.
   !> @param[inout] this Lattice object whose current members provide the input state.
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

   !> @brief Compute primitive-cell derived quantities.
   !> @details Converts lattice vectors to internal units, computes cell volume
   !>          and Wigner-Seitz radius, and prepares inverse Cartesian lattice
   !>          data used by periodic wrapping and neighbor searches.
   !> @param[inout] this Lattice object whose a/alat/celldm state is updated.
   module subroutine build_data(this)
      class(lattice), intent(inout) :: this
      !> Local variables
      real(rp), dimension(3, 3) :: a
      integer :: i, j

   end subroutine build_data

   !> @brief Reset lattice members to their default values.
   !> @details Clears allocatable storage and restores scalar defaults; with the
   !>          optional full flag it also clears persistent input-facing state.
   !> @param[inout] this Lattice object to reset.
   !> @param[in] full Optional flag requesting a full reset.
   module subroutine restore_to_default(this, full)
      class(lattice) :: this
      logical, intent(in), optional :: full

   end subroutine restore_to_default

   !> @brief Build the bulk Bravais cluster.
   !> @details Expands primitive-cell sites through lattice translations, cuts
   !>          them by the configured radius, and fills bulk cluster coordinates,
   !>          atomic labels, and optional Morton ordering.
   !> @param[inout] this Lattice object receiving cr/iz/num/no cluster state.
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

   !> @brief Reorder the bulk cluster by Morton space-filling-curve key.
   !> @details Sorts cluster atoms in normalized Cartesian space to improve
   !>          locality for real-space recursion kernels while preserving the
   !>          per-atom coordinate/type arrays as a consistent permutation.
   !> @param[inout] this Lattice object whose bulk cluster arrays are reordered.
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

   !> @brief Encode a 3D integer grid coordinate as a Morton key.
   !> @param[in] x Grid coordinate on the first axis.
   !> @param[in] y Grid coordinate on the second axis.
   !> @param[in] z Grid coordinate on the third axis.
   !> @param[in] nbits Number of bits to interleave from each coordinate.
   !> @return Interleaved 64-bit Morton code.
   module pure function morton_encode3(x, y, z, nbits) result(code)
      integer, intent(in) :: x, y, z, nbits
      integer(8) :: code
      integer :: b
   end function morton_encode3

   !> @brief Sort an integer permutation by associated Morton keys.
   !> @details Uses an in-place heapsort so key and permutation arrays are kept
   !>          synchronized without extra sorting dependencies.
   !> @param[inout] key Keys to sort in ascending order.
   !> @param[inout] perm Permutation entries carried with each key.
   !> @param[in] n Number of active entries in key and perm.
   module subroutine sort_by_key(key, perm, n)
      integer(8), intent(inout) :: key(:)
      integer, intent(inout) :: perm(:)
      integer, intent(in) :: n
      integer :: i, start, bottom
      integer(8) :: tkey
      integer :: tperm

      ! Build heap then sift down (standard heapsort, O(n log n)).
   end subroutine sort_by_key

   !> @brief Prepare surface-layer cluster metadata from a bulk cluster.
   !> @details Builds layer z positions and surface index arrays used by the
   !>          buildsurf path before selecting the active surface atoms.
   !> @param[inout] this Lattice object whose surface helper arrays are filled.
   module subroutine build_clusup(this)
      class(lattice), intent(inout) :: this
      character(len=10) :: surftype
      ! Local variables
      real(rp), dimension(:), allocatable :: z
      integer, dimension(:), allocatable :: izpsurf
      real(rp) :: ds, ds2, new, one
      integer :: i, j, n

   end subroutine build_clusup

   !> @brief Build the full surface cluster representation.
   !> @details Selects atoms by surface orientation and layer bounds, identifies
   !>          unique layer/type representatives, and installs the full surface
   !>          coordinate/type arrays for downstream surface calculations.
   !> @param[inout] this Lattice object receiving the full surface cluster.
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

   !> @brief Build the two-sided (region A | active | region B) interface cluster.
   !> @details The calctype='L' counterpart of build_surf_full: instead of one
   !>          vacuum boundary and one bulk boundary, it places two independent
   !>          FROZEN semi-infinite references around a central active zone.
   !>          Layers are selected by projecting onto the surface normal exactly
   !>          as build_surf_full does; the only structural difference is the
   !>          layer->type assignment, which now draws region-A types from the
   !>          leading block of the &atoms label list and region-B types from
   !>          the trailing block.
   !> @param[inout] this Lattice object receiving the interface cluster.
   module subroutine build_interface_full(this)
      class(lattice), intent(inout) :: this
   end subroutine build_interface_full

   !> @brief Build the legacy compact surface cluster.
   !> @details Chooses representative atoms from the bulk cluster for each
   !>          requested surface layer and writes the compact surface cluster
   !>          arrays consumed by newclu/buildsurf workflows.
   !> @param[inout] this Lattice object receiving compact surface state.
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

   !> @brief Build an impurity/local cluster from a bulk or surface host.
   !> @details Combines host cluster coordinates with impurity inclusions,
   !>          creates local atom/type maps, and builds neighbor tables for
   !>          impurity and defect calculations.
   !> @param[inout] this Lattice object receiving impurity cluster arrays.
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

   !> @brief Build neighbor maps and structure constants.
   !> @details Constructs the neighbor table around representative atoms, then
   !>          either computes screened structure constants or only prepares the
   !>          neighbor geometry depending on do_str.
   !> @param[inout] this Lattice object receiving nn/sbar/sdot/sbarvec state.
   !> @param[in] do_str If true, compute structure constants after neighbor mapping.
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

   !> @brief Allocate storage for strux structure constants.
   !> @param[inout] this Lattice object whose sbar/sdot/alpha storage is reset.
   !> @param[in] sbar_dim Orbital dimension of each structure-constant block.
   !> @param[in] nm Number of neighbor slots to store.
   module subroutine init_strux_storage(this, sbar_dim, nm)
      class(lattice), intent(inout) :: this
      integer, intent(in) :: sbar_dim, nm

   end subroutine init_strux_storage

   !> @brief Return default l-resolved screening constants.
   !> @param[in] this Lattice object providing context for the pure helper.
   !> @param[in] nl Number of angular-momentum channels requested.
   !> @return Default screening-alpha array of length nl.
   module pure function default_screening_alpha(this, nl) result(alpha_default)
      class(lattice), intent(in) :: this
      integer, intent(in) :: nl
      real(rp) :: alpha_default(nl)
      integer :: i

   end function default_screening_alpha

   !> @brief Convert the screening string to a strux-library mode code.
   !> @param[in] this Lattice object containing the screening option.
   !> @return STRUX_LMTO47_IALPHA_* selector used by strux_compute.
   module integer function strux_mode(this)
      class(lattice), intent(in) :: this
      character(len=:), allocatable :: screening_mode

   end function strux_mode

   !> @brief Load symbolic-atom data needed by strux helpers.
   !> @details Lazily populates this%symbolic_atoms from the control input when
   !>          muffin-tin radii or species labels are required.
   !> @param[inout] this Lattice object whose symbolic_atoms cache may be filled.
   module subroutine load_symbolic_atoms_if_needed(this)
      class(lattice), intent(inout) :: this
      
      ! 1. Declare the temporary array. 
      type(symbolic_atom), allocatable :: temp_atoms(:)

   end subroutine load_symbolic_atoms_if_needed

   !> @brief Build species muffin-tin radii for strux input.
   !> @param[inout] this Lattice object providing symbolic-atom data.
   !> @param[in] nspec Number of distinct species.
   !> @param[in] species_labels Atomic/species labels to map.
   !> @param[out] rmt Muffin-tin radius for each species.
   module subroutine build_rmt(this, nspec, species_labels, rmt)
      class(lattice), intent(inout) :: this
      integer, intent(in) :: nspec
      integer, intent(in) :: species_labels(nspec)
      real(rp), intent(out) :: rmt(nspec)
      integer :: is, label

   end subroutine build_rmt

   !> @brief Build hard-core radii for strux sigma screening.
   !> @param[in] this Lattice object containing screening_sigma settings.
   !> @param[in] nl Number of angular-momentum channels.
   !> @param[in] nspec Number of distinct species.
   !> @param[in] rmt Muffin-tin radii for each species.
   !> @param[out] hcr Hard-core radii indexed by l channel and species.
   module subroutine build_hcr(this, nl, nspec, rmt, hcr)
      class(lattice), intent(in) :: this
      integer, intent(in) :: nl, nspec
      real(rp), intent(in) :: rmt(nspec)
      real(rp), intent(out) :: hcr(nl, nspec)
      integer :: is, l

   end subroutine build_hcr

   !> @brief Assemble alpha, hard-core, and muffin-tin arrays for strux_compute.
   !> @details Combines manual/default screening inputs with per-species radius
   !>          data in the layout expected by the local strux library.
   !> @param[inout] this Lattice object providing screening and symbolic-atom data.
   !> @param[in] nspec Number of species in the strux solve cluster.
   !> @param[in] nl Number of angular-momentum channels.
   !> @param[in] species_labels Species labels in strux order.
   !> @param[out] alpha_in Screening constants in strux l/species layout.
   !> @param[out] hcr Hard-core radii in strux l/species layout.
   !> @param[out] rmt Muffin-tin radii by species.
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

   !> @brief Compute screened structure constants through the strux backend.
   !> @details Maps the RS-LMTO cluster and basis labels into strux input,
   !>          calls strux_compute, and remaps the returned screened constants
   !>          and optional derivatives back into lattice storage.
   !> @param[inout] this Lattice object receiving sbar/sdot/alpha data.
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

   !> @brief Write one strux structure-constant neighbor block for diagnostics.
   !> @param[in] this Lattice object containing sbar and neighbor vectors.
   !> @param[in] nl2 Block dimension to write.
   !> @param[in] m Neighbor slot.
   !> @param[in] ii Representative atom/type index.
   !> @param[in] jclus Cluster atom index associated with the block.
   module subroutine write_strux_block(this, nl2, m, ii, jclus)
      class(lattice), intent(in) :: this
      integer, intent(in) :: nl2, m, ii, jclus
      integer :: is, js

      ! write (*, '(" SBAR neighbor center=",i5," slot=",i5," iclus=",i5," vec=",3f12.6)') ii, m, iclus, &
      !    this%sbarvec(1, m), this%sbarvec(2, m), this%sbarvec(3, m)
   end subroutine write_strux_block

   !> @brief Write stored neighbor vectors to a diagnostic unit.
   !> @param[in] iunit Output unit.
   !> @param[in] sbarvec Neighbor vectors to dump.
   !> @param[in] nt Number of active vectors.
   module subroutine write_neighbor_vector_dump(iunit, sbarvec, nt)
      integer, intent(in) :: iunit, nt
      real(rp), intent(in) :: sbarvec(:, :)
      integer :: i

   end subroutine write_neighbor_vector_dump

   !> @brief Return the representative cluster atom for an inequivalent type.
   !> @param[in] this Lattice object containing irec/iu mappings.
   !> @param[in] ii Inequivalent atom/type index.
   !> @return Cluster atom index used as representative for ii.
   module integer function representative_atom_index(this, ii)
      class(lattice), intent(in) :: this
      integer, intent(in) :: ii

   end function representative_atom_index

   !> @brief Return the primitive-basis label associated with a cluster atom.
   !> @param[in] this Lattice object containing primitive and cluster labels.
   !> @param[in] ia Cluster atom index.
   !> @return Primitive-cell basis label for ia.
   module integer function primitive_basis_label(this, ia)
      class(lattice), intent(in) :: this
      integer, intent(in) :: ia

   end function primitive_basis_label

   !> @brief Find a strux pair table entry matching a target neighbor vector.
   !> @details Searches lattice translations and basis-pair labels for the
   !>          pair whose Cartesian displacement matches vec_target within tol.
   !> @param[in] nttab Number of pair-table entries.
   !> @param[in] iax Pair translation/index table from strux.
   !> @param[in] plat Primitive lattice vectors.
   !> @param[in] pos Basis positions.
   !> @param[in] alat Lattice parameter used to scale vectors.
   !> @param[in] ib Source basis label.
   !> @param[in] jb Target basis label.
   !> @param[in] vec_target Target displacement vector.
   !> @param[in] tol Matching tolerance.
   !> @return Matching pair-table index, or zero if no match is found.
   module integer function find_pair_by_vector(nttab, iax, plat, pos, alat, ib, jb, vec_target, tol)
      integer, intent(in) :: nttab, ib, jb
      integer, intent(in) :: iax(:,:)
      real(rp), intent(in) :: plat(3, 3), pos(:, :), alat, vec_target(3), tol
      integer :: i, n1, n2, n3
      real(rp) :: vec(3)
      real(rp) :: dmax, best_d
      integer :: best_i

   end function find_pair_by_vector

   !> @brief Find a cluster neighbor atom by displacement vector.
   !> @param[in] this Lattice object containing cluster metadata.
   !> @param[in] ia Source cluster atom.
   !> @param[in] vec_target Target displacement in lattice-scaled coordinates.
   !> @param[in] cralat Cluster coordinates scaled by the lattice parameter.
   !> @param[in] tol Matching tolerance.
   !> @return Cluster atom index matching the displacement, or zero if not found.
   module integer function find_neighbor_atom_by_vector(this, ia, vec_target, cralat, tol)
      class(lattice), intent(in) :: this
      integer, intent(in) :: ia
      real(rp), intent(in) :: vec_target(:), cralat(:, :)
      real(rp), intent(in) :: tol
      integer :: ja
      real(rp) :: dv(3)

   end function find_neighbor_atom_by_vector

   !> @brief Build the orbital-order map between RS-LMTO and strux layouts.
   !> @param[in] norb Number of orbitals in the active basis.
   !> @param[out] orb_map Mapping from local orbital index to strux orbital index.
   module subroutine build_orbital_map(norb, orb_map)
      integer, intent(in) :: norb
      integer, intent(out) :: orb_map(16)
      integer :: i
      integer, parameter :: full_map(16) = [ &
         1, 4, 2, 3, 5, 6, 8, 9, 7, 13, 14, 12, 15, 11, 16, 10 ]

   end subroutine build_orbital_map

   !> @brief Build atom-list metadata for self-consistent atom types.
   !> @details Fills atlist/chargetrf_type-style mappings used by impurity and
   !>          local-cluster charge transfer paths.
   !> @param[inout] this Lattice object whose atom-list arrays are updated.
   module subroutine atomlist(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      integer :: i, j, itype
      real(rp) :: mom_tmp(3)

   end subroutine atomlist

   !> @brief Select cluster atoms inside a primitive-cell volume.
   !> @details Tests positions relative to a central atom against the three
   !>          primitive vectors and returns the atom indices that lie inside.
   !> @param[inout] this Lattice object providing volume-test helper methods.
   !> @param[in] cr Cluster coordinates.
   !> @param[in] num Cluster atom labels.
   !> @param[in] num_atoms Number of atoms in cr/num.
   !> @param[in] central_atom Index of the atom used as the volume origin.
   !> @param[in] a1 First primitive vector.
   !> @param[in] a2 Second primitive vector.
   !> @param[in] a3 Third primitive vector.
   !> @param[in] plane_constant Boundary tolerance/plane constant.
   !> @param[out] atoms_in_volume Atom indices found inside the volume.
   !> @param[out] atom_count Number of returned atoms.
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

   !> @brief Test whether a relative position is inside a primitive parallelepiped.
   !> @param[inout] this Lattice object providing the type-bound helper context.
   !> @param[in] relative_pos Position measured from the cell origin.
   !> @param[in] a1 First primitive vector.
   !> @param[in] a2 Second primitive vector.
   !> @param[in] a3 Third primitive vector.
   !> @param[out] inside True if the point is inside the primitive volume.
   module subroutine check_within_volume(this, relative_pos, a1, a2, a3, inside)
       class(lattice), intent(inout) :: this
       real(rp), intent(in) :: relative_pos(3), a1(3), a2(3), a3(3)
       logical, intent(out) :: inside
       real(rp) :: dot11, dot12, dot13, dot22, dot23, dot33
       real(rp) :: dot1r, dot2r, dot3r, inv_denom, u, v, w
   
       ! Calculate dot products between the vectors
   end subroutine check_within_volume

   !> @brief Return unique integer structure/type labels.
   !> @param[inout] this Lattice object providing the type-bound helper context.
   !> @param[in] num Input labels.
   !> @param[in] num_atoms Number of active labels.
   !> @param[out] unique_nums Unique labels in first-seen order.
   module subroutine find_unique_struct(this, num, num_atoms, unique_nums)
       class(lattice), intent(inout) :: this
       integer, intent(in) :: num(:), num_atoms
       integer, allocatable, intent(out) :: unique_nums(:)
       integer, allocatable :: temp_nums(:)
       integer :: i, j, num_unique
       logical :: found
   
       ! Allocate temporary array for unique numbers
   end subroutine find_unique_struct

   !> @brief Identify symmetry-unique atoms inside a primitive volume.
   !> @details Removes atoms related by primitive translations from an input
   !>          volume selection and returns one representative per unique site.
   !> @param[inout] this Lattice object providing helper context.
   !> @param[in] cr Cluster coordinates.
   !> @param[in] num_atoms Number of atoms in cr.
   !> @param[in] atoms_in_volume Candidate atom indices.
   !> @param[in] atom_count Number of candidate atoms.
   !> @param[in] a1 First primitive vector.
   !> @param[in] a2 Second primitive vector.
   !> @param[in] a3 Third primitive vector.
   !> @param[out] unique_atoms Representative atom indices.
   !> @param[out] unique_atom_count Number of representatives.
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

   !> @brief Build one legacy screened-structure-constant block.
   !> @details Forms a local cluster around atom ia and drives the legacy
   !>          CLUSBA/MICHA screening path used by the non-strux backend.
   !> @param[inout] this Lattice object receiving legacy sbar data.
   !> @param[in] ia Center atom index.
   !> @param[in] r2 Cluster cutoff radius squared.
   !> @param[in] wav Wigner-Seitz radius.
   !> @param[in] crd Cluster coordinates.
   !> @param[in] nat Number of atoms in crd.
   !> @param[in] ndi Leading dimension/capacity of crd.
   !> @param[in] np Number of primitive atoms/neighbor entries.
   !> @param[in] nr Number of local-cluster vectors.
   !> @param[in] ii Inequivalent atom/type index.
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

   !> @brief Collect local cluster vectors within a cutoff.
   !> @details Selects neighbors around atom ia from crd, sorted by distance, for
   !>          use in the legacy structure-constant screening routines.
   !> @param[inout] this Lattice object providing helper context.
   !> @param[in] r2 Cutoff radius squared.
   !> @param[in] crd Candidate cluster coordinates.
   !> @param[in] ia Center atom index.
   !> @param[in] nat Number of candidate atoms.
   !> @param[in] ndi Leading dimension/capacity of crd.
   !> @param[inout] n Number of vectors found.
   !> @param[inout] sbarvec_out Optional output neighbor-vector storage.
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

   !> @brief Run the legacy MICHA screening construction.
   !> @details Builds raw structure constants, applies shell decomposition and
   !>          Cholesky-style screening, and accumulates screened sbar blocks.
   !> @param[in] rws Wigner-Seitz radius.
   !> @param[in] r Local cluster vectors.
   !> @param[in] nr Number of local cluster vectors.
   !> @param[in] nlm Orbital block dimension.
   !> @param[in] nrl Reduced screening dimension.
   !> @param[in] na Length of screening-coefficient storage.
   !> @param[inout] sbar Screened structure constants.
   !> @param[inout] a Screening coefficient work array.
   !> @param[inout] wk Work array.
   !> @param[inout] bet Screening beta work array.
   !> @param[inout] s Raw structure-constant work matrix.
   !> @param[in] iclus Center-cluster index.
   !> @param[in] r2 Cutoff radius squared.
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

   !> @brief Build raw canonical structure constants for a local cluster.
   !> @param[in] w Distance scaling factor.
   !> @param[in] r Local cluster vectors.
   !> @param[in] nr Number of local cluster vectors.
   !> @param[inout] s Raw structure-constant matrix.
   !> @param[in] nrl Reduced screening dimension.
   !> @param[in] nlm Orbital block dimension.
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

   !> @brief Apply shell decomposition and screening to legacy structure constants.
   !> @param[in] r Local cluster vectors.
   !> @param[in] nr Number of local cluster vectors.
   !> @param[in] nlm Orbital block dimension.
   !> @param[in] nrl Reduced screening dimension.
   !> @param[inout] s Raw/screened structure-constant work matrix.
   !> @param[inout] a Screening coefficient work array.
   !> @param[in] na Length of a.
   !> @param[in] q Angular-momentum shell weights.
   !> @param[inout] bet Screening beta work array.
   !> @param[inout] wk Work array.
   !> @param[inout] sbar Screened structure constants.
   !> @param[in] iclus Center-cluster index.
   !> @param[in] r2 Cutoff radius squared.
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

   !> @brief Evaluate canonical structure-constant angular blocks.
   !> @param[in] w Distance scaling factor.
   !> @param[in] dr Displacement vector.
   !> @param[out] sc 16-by-16 orbital structure-constant block.
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

   !> @brief Factor a packed real symmetric matrix for legacy screening.
   !> @param[inout] c Packed matrix/coefficient storage.
   !> @param[in] na Length of c.
   !> @param[inout] w Work/vector storage.
   !> @param[in] n Matrix order.
   !> @param[out] ndef Number of non-positive/deficient pivots detected.
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

   !> @brief Solve against a packed Cholesky factor for legacy screening.
   !> @param[in] c Packed factor storage.
   !> @param[in] na Length of c.
   !> @param[inout] v Right-hand sides overwritten by solved vectors.
   !> @param[in] n Matrix order.
   !> @param[in] m Number of right-hand sides.
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

   !> @brief Map an LM orbital index to angular-momentum quantum number l.
   !> @param[in] ilm One-based LM orbital index.
   !> @return Angular-momentum shell index for ilm.
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

   !> @brief Build representative neighbor vectors for each atom type.
   !> @details Compares cluster coordinates around representative atoms and
   !>          fills neighbor maps plus displacement sets used by structb.
   !> @param[inout] this Lattice object providing tolerance and helper context.
   !> @param[in] crd Cluster coordinates.
   !> @param[in] no Atom type labels for cluster atoms.
   !> @param[in] iu Representative atom indices.
   !> @param[inout] nn Neighbor table to fill.
   !> @param[in] nat Number of cluster atoms.
   !> @param[in] ntot Number of representative atoms.
   !> @param[in] nomx Number of representative/type slots.
   !> @param[in] ndi Leading dimension/capacity of crd.
   !> @param[in] nnmx Maximum neighbor slots.
   !> @param[inout] set Neighbor displacement vectors.
   !> @param[inout] idnn Neighbor identifier list.
   !> @param[out] ret Return/status vector used by legacy callers.
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

   !> @brief Write the neighbor-map table in legacy text format.
   !> @param[in] IM Output unit.
   !> @param[in] IZP Atomic numbers/labels for printed atoms.
   !> @param[in] NN Neighbor table.
   !> @param[in] NO Atom type labels.
   !> @param[in] ND Number of atoms in NO/IZP.
   !> @param[in] NM Number of neighbor slots.
   !> @param[in] NTOT Number of representative atoms.
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

   !> @brief Classify a neighbor-shell pair by distance cutoff.
   !> @param[in] I First shell/type index.
   !> @param[in] J Second shell/type index.
   !> @param[in] R2 Squared distance between sites.
   !> @param[in] DD Distance tolerance or shell spacing.
   !> @param[in] CT Shell cutoff radii.
   !> @return Legacy neighbor-shell map value.
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

   !> @brief Compute the minimum-image coordinate difference between two atoms.
   !> @param[inout] this Lattice object containing periodic-boundary settings.
   !> @param[in] Natom Number of atoms in coord.
   !> @param[in] coord Atomic coordinates.
   !> @param[in] i_atom First atom index.
   !> @param[in] j_atom Second atom index.
   !> @param[out] cdiff Minimum-image displacement from i_atom to j_atom.
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

   !> @brief Build the atom neighbor table using shell cutoffs.
   !> @details Uses spatial binning and optional periodic wrapping to find
   !>          neighbors within the configured cutoff shells and populate NN.
   !> @param[inout] this Lattice object containing cell/PBC state.
   !> @param[inout] ct Shell cutoff radii.
   !> @param[in] crd Atomic coordinates.
   !> @param[in] ndim Coordinate leading dimension.
   !> @param[in] nat Number of atoms.
   !> @param[in] izp Atomic numbers/labels.
   !> @param[inout] nn Neighbor table.
   !> @param[in] nd Number of atoms represented in nn.
   !> @param[inout] nm Maximum/actual neighbor slots.
   !> @param[in] ngbr Legacy neighbor-shell classification function.
   !> @param[in] ntot Number of representative atoms.
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

   !> @brief Read legacy cluster coordinate records from a unit.
   !> @param[in] alat Lattice parameter used for coordinate scaling.
   !> @param[in] nndim Maximum number of cluster atoms.
   !> @param[inout] cr Coordinate array to fill.
   !> @param[inout] iz Atomic numbers/labels.
   !> @param[inout] n Atom-number labels.
   !> @param[in] ip Input unit.
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

   !> @brief Sort a legacy matrix block by its first column.
   !> @param[in] nl First active row/starting index.
   !> @param[in] ndim Leading dimension of m.
   !> @param[inout] m Matrix block to sort.
   !> @param[in] nd Number of columns.
   !> @param[in] nt Number of active rows.
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

   !> @brief Cut translated primitive-cell atoms to a spherical cluster.
   !> @param[in] i Center atom index.
   !> @param[in] l Number of candidate atoms.
   !> @param[in] ndim Array capacity.
   !> @param[in] crd Candidate coordinates.
   !> @param[out] cr Coordinates that survive the cut.
   !> @param[in] izp Candidate atomic labels.
   !> @param[out] iz Atomic labels that survive the cut.
   !> @param[out] num Atom-number labels that survive the cut.
   !> @param[in] no Candidate atom-number labels.
   !> @param[in] rs Cut radius.
   !> @param[out] ii Number of atoms that survived the cut.
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

   !> @brief Check lattice option consistency after input parsing.
   !> @param[inout] this Lattice object whose options are validated.
   module subroutine check_all(this)
      implicit none
      class(lattice) :: this

   end subroutine check_all

   !> @brief Count reduced basis/species entries from the cluster atom labels.
   !> @details Computes nbas/reduced_nbas-style counts used by charge and
   !>          Hamiltonian setup after the cluster atom types are known.
   !> @param[inout] this Lattice object whose basis counters are updated.
   module subroutine calculate_nbas(this)
      implicit none
      class(lattice) :: this
      integer :: size_iz
      integer, dimension(:) :: atype(this%nbas), amount(this%nbas)
      integer :: i, j

   end subroutine calculate_nbas

   !> @brief Print the complete lattice state in the legacy diagnostic format.
   !> @param[in] this Lattice object to print.
   !> @param[in] unit Optional output unit.
   !> @param[in] file Optional output file path.
   module subroutine print_state_full(this, unit, file)
      class(lattice), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

   end subroutine print_state_full

   !> @brief Print the compact lattice state in the legacy diagnostic format.
   !> @param[in] this Lattice object to print.
   !> @param[in] unit Optional output unit.
   !> @param[in] file Optional output file path.
   module subroutine print_state(this, unit, file)
      implicit none
      class(lattice), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

   end subroutine print_state

   !> @brief Print lattice state as a generated namelist.
   !> @details Emits the current lattice configuration with namelist_generator
   !>          so diagnostics can be compared with input-facing settings.
   !> @param[in] this Lattice object to print.
   !> @param[in] unit Optional output unit.
   !> @param[in] file Optional output file path.
   module subroutine print_state_formatted(this, unit, file)
      class(lattice), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file

      type(namelist_generator) :: nml
      integer :: i

   end subroutine print_state_formatted

   end interface

end module lattice_mod
