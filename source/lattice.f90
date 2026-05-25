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

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] fname Namelist file
   !> @return type(lattice)
   !---------------------------------------------------------------------------
   function constructor(control_obj) result(obj)
      type(lattice) :: obj
      type(control), target, intent(in) :: control_obj

      obj%control => control_obj

      call obj%restore_to_default()
      call obj%build_from_file()
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(lattice) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%ct)) call g_safe_alloc%deallocate('lattice.ct', this%ct)
      if (allocated(this%cr)) call g_safe_alloc%deallocate('lattice.cr', this%cr)
      if (allocated(this%crd)) call g_safe_alloc%deallocate('lattice.crd', this%crd)
      if (allocated(this%iu)) call g_safe_alloc%deallocate('lattice.iu', this%iu)
      if (allocated(this%ib)) call g_safe_alloc%deallocate('lattice.ib', this%ib)
      if (allocated(this%irec)) call g_safe_alloc%deallocate('lattice.irec', this%irec)
      if (allocated(this%izp)) call g_safe_alloc%deallocate('lattice.izp', this%izp)
      if (allocated(this%iz)) call g_safe_alloc%deallocate('lattice.iz', this%iz)
      if (allocated(this%num)) call g_safe_alloc%deallocate('lattice.num', this%num)
      if (allocated(this%no)) call g_safe_alloc%deallocate('lattice.no', this%no)
      if (allocated(this%z)) call g_safe_alloc%deallocate('lattice.z', this%z)
      if (allocated(this%izpsurf)) call g_safe_alloc%deallocate('lattice.izpsurf', this%izpsurf)
      if (allocated(this%izsurf)) call g_safe_alloc%deallocate('lattice.izsurf', this%izsurf)
      if (allocated(this%nosurf)) call g_safe_alloc%deallocate('lattice.nosurf', this%nosurf)
      if (allocated(this%crsurf)) call g_safe_alloc%deallocate('lattice.crsurf', this%crsurf)
      if (allocated(this%inclu)) call g_safe_alloc%deallocate('lattice.inclu', this%inclu)
      if (allocated(this%izpo)) call g_safe_alloc%deallocate('lattice.izpo', this%izpo)
      if (allocated(this%acr)) call g_safe_alloc%deallocate('lattice.acr', this%acr)
      if (allocated(this%atlist)) call g_safe_alloc%deallocate('lattice.atlist', this%atlist)
      if (allocated(this%nn)) call g_safe_alloc%deallocate('lattice.nn', this%nn)
      if (allocated(this%sbar)) call g_safe_alloc%deallocate('lattice.sbar', this%sbar)
      if (allocated(this%sdot)) call g_safe_alloc%deallocate('lattice.sdot', this%sdot)
      if (allocated(this%sbarvec)) call g_safe_alloc%deallocate('lattice.sbarvec', this%sbarvec)
      if (allocated(this%alpha)) call g_safe_alloc%deallocate('lattice.alpha', this%alpha)
      if (allocated(this%alpha_dot)) call g_safe_alloc%deallocate('lattice.alpha_dot', this%alpha_dot)
      if (allocated(this%screening_alpha)) call g_safe_alloc%deallocate('lattice.screening_alpha', this%screening_alpha)
      if (allocated(this%screening_sigma)) call g_safe_alloc%deallocate('lattice.screening_sigma', this%screening_sigma)
      if (allocated(this%ijpair)) call g_safe_alloc%deallocate('lattice.ijpair', this%ijpair)
      if (allocated(this%ijktrio)) call g_safe_alloc%deallocate('lattice.ijktrio', this%ijktrio)
      if (allocated(this%natoms_layer)) call g_safe_alloc%deallocate('lattice.natoms_layer', this%natoms_layer)
#else
      if (allocated(this%ct)) deallocate (this%ct)
      if (allocated(this%cr)) deallocate (this%cr)
      if (allocated(this%crd)) deallocate (this%crd)
      if (allocated(this%iu)) deallocate (this%iu)
      if (allocated(this%ib)) deallocate (this%ib)
      if (allocated(this%irec)) deallocate (this%irec)
      if (allocated(this%izp)) deallocate (this%izp)
      if (allocated(this%iz)) deallocate (this%iz)
      if (allocated(this%num)) deallocate (this%num)
      if (allocated(this%no)) deallocate (this%no)
      if (allocated(this%z)) deallocate (this%z)
      if (allocated(this%izpsurf)) deallocate (this%izpsurf)
      if (allocated(this%izsurf)) deallocate (this%izsurf)
      if (allocated(this%nosurf)) deallocate (this%nosurf)
      if (allocated(this%crsurf)) deallocate (this%crsurf)
      if (allocated(this%inclu)) deallocate (this%inclu)
      if (allocated(this%izpo)) deallocate (this%izpo)
      if (allocated(this%acr)) deallocate (this%acr)
      if (allocated(this%atlist)) deallocate (this%atlist)
      if (allocated(this%nn)) deallocate (this%nn)
      if (allocated(this%sbar)) deallocate (this%sbar)
      if (allocated(this%sdot)) deallocate (this%sdot)
      if (allocated(this%sbarvec)) deallocate (this%sbarvec)
      if (allocated(this%alpha)) deallocate (this%alpha)
      if (allocated(this%alpha_dot)) deallocate (this%alpha_dot)
      if (allocated(this%screening_alpha)) deallocate (this%screening_alpha)
      if (allocated(this%screening_sigma)) deallocate (this%screening_sigma)
      if (allocated(this%ijpair)) deallocate (this%ijpair)
      if (allocated(this%ijktrio)) deallocate (this%ijktrio)
      if (allocated(this%natoms_layer)) deallocate (this%natoms_layer)
#endif
   end subroutine destructor

   ! Member functions
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read parameters from input file
   !
   !> @param[in] fname Namelist file
   !---------------------------------------------------------------------------
   subroutine build_from_file(this, fname)
      class(lattice), intent(inout) :: this
      character(len=*), intent(in), optional :: fname
      ! variables associated with the reading processes
      integer :: iostatus, funit, i
      character(len=sl) :: fname_
      real(rp), allocatable :: ct_first(:), screening_alpha_first(:), screening_sigma_first(:)

      include 'include_codes/namelists/lattice.f90'

      if (present(fname)) then
         fname_ = fname
         this%control%fname = fname
      else
         fname_ = this%control%fname
      end if

      ! Save previous values
      ! Bulk initialization
      ndim = this%ndim
      npe = this%npe
      rc = this%rc
      r2 = this%r2
      ntype = this%ntype
      crystal_sym = this%crystal_sym
      a = this%a
      wav = this%wav
      celldm = this%celldm
      strux_backend = this%strux_backend
      screening = this%screening
      strux_want_sdot = this%strux_want_sdot
      strux_solve_scale = this%strux_solve_scale
      b1 = this%b1
      b2 = this%b2
      b3 = this%b3
      n1 = this%n1 
      n2 = this%n2
      n3 = this%n3
      pbc = this%pbc

      if (.not. allocated(this%izp) .or. size(this%izp) .ne. this%ndim) then
#ifdef USE_SAFE_ALLOC
         if (allocated(this%izp)) call g_safe_alloc%deallocate('lattice.izp', this%izp)
         call g_safe_alloc%allocate('lattice.izp', this%izp, (/this%ndim/))
#else
         if (allocated(this%izp)) deallocate (this%izp)
         allocate (this%izp(this%ndim))
#endif
      end if
      if (.not. allocated(this%no) .or. size(this%no) .ne. this%ndim) then
#ifdef USE_SAFE_ALLOC
         if (allocated(this%no)) call g_safe_alloc%deallocate('lattice.no', this%no)
         call g_safe_alloc%allocate('lattice.no', this%no, (/this%ndim/))
#else
         if (allocated(this%no)) deallocate (this%no)
         allocate (this%no(this%ndim))
#endif
      end if
      if (.not. allocated(this%crd) .or. size(this%crd) .ne. 3*this%ndim) then
#ifdef USE_SAFE_ALLOC
         if (allocated(this%crd)) call g_safe_alloc%deallocate('lattice.crd', this%crd)
         call g_safe_alloc%allocate('lattice.crd', this%crd, (/3, this%ndim/))
#else
         if (allocated(this%crd)) deallocate (this%crd)
         allocate (this%crd(3, this%ndim))
#endif
      end if

      call move_alloc(this%izp, izp)
      call move_alloc(this%no, no)
      call move_alloc(this%crd, crd)

      ! Impurity initialization
      nclu = this%nclu
      if (.not. allocated(this%inclu) .or. size(this%inclu) .ne. 3*this%nclu) then
#ifdef USE_SAFE_ALLOC
         if (allocated(this%inclu)) call g_safe_alloc%deallocate('lattice.inclu', this%inclu)
         call g_safe_alloc%allocate('lattice.inclu', this%inclu, (/this%nclu, 3/))
#else
         if (allocated(this%inclu)) deallocate (this%inclu)
         allocate (this%inclu(this%nclu, 3))
#endif
      end if
      call move_alloc(this%inclu, inclu)

      ! Surface initialization
      surftype = this%surftype
      nlay = this%nlay
      ntype = this%ntype
      call move_alloc(this%ct, ct)
      call move_alloc(this%screening_alpha, screening_alpha)
      call move_alloc(this%screening_sigma, screening_sigma)
      njij = this%njij
      call move_alloc(this%ijpair, ijpair)
      njijk = this%njijk
      call move_alloc(this%ijktrio, ijktrio)

      ! Pre-size ct before the first read. ntype is unknown at this point
      ! (local ntype = 0 from restore_to_default), so ct has size 0. 1000 is
      ! a safe upper bound for atom types; the resize check below will shrink
      ! it to the actual ntype read from the file.
      if (.not. allocated(ct)) then
         allocate (ct(size(izp)))
      else if (size(ct) == 0) then
         deallocate (ct)
         allocate (ct(size(izp)))
      end if
      if (.not. allocated(screening_alpha)) then
         allocate (screening_alpha(4))
         screening_alpha = this%default_screening_alpha(size(screening_alpha))
      else if (size(screening_alpha) < 4) then
         deallocate (screening_alpha)
         allocate (screening_alpha(4))
         screening_alpha = this%default_screening_alpha(size(screening_alpha))
      end if
      if (.not. allocated(screening_sigma)) then
         allocate (screening_sigma(4))
         screening_sigma = 0.7_rp
      else if (size(screening_sigma) < 4) then
         deallocate (screening_sigma)
         allocate (screening_sigma(4))
         screening_sigma = 0.7_rp
      end if
      if (.not. allocated(z)) then
         allocate (z(1))
         z = 0.0_rp
      end if
      if (.not. allocated(primcell)) then
         allocate (primcell(1, 1))
         primcell = 0.0_rp
      end if
      if (.not. allocated(crsurf)) then
         allocate (crsurf(1, 1))
         crsurf = 0.0_rp
      end if
      if (.not. allocated(cr)) then
         allocate (cr(1, 1))
         cr = 0.0_rp
      end if
      if (.not. allocated(acr)) then
         allocate (acr(1, 1))
         acr = 0.0_rp
      end if
      if (.not. allocated(reduced_acr)) then
         allocate (reduced_acr(1))
         reduced_acr = 0
      end if
      if (.not. allocated(num)) then
         allocate (num(1))
         num = 0
      end if
      if (.not. allocated(izpsurf)) then
         allocate (izpsurf(1))
         izpsurf = 0
      end if
      if (.not. allocated(izsurf)) then
         allocate (izsurf(1))
         izsurf = 0
      end if
      if (.not. allocated(nosurf)) then
         allocate (nosurf(1))
         nosurf = 0
      end if
      if (.not. allocated(izpo)) then
         allocate (izpo(1))
         izpo = 0
      end if
      if (.not. allocated(iz)) then
         allocate (iz(1))
         iz = 0
      end if
      if (.not. allocated(iu)) then
         allocate (iu(1))
         iu = 0
      end if
      if (.not. allocated(irec)) then
         allocate (irec(1))
         irec = 0
      end if
      if (.not. allocated(ib)) then
         allocate (ib(1))
         ib = 0
      end if
      if (.not. allocated(ijpair)) then
         allocate (ijpair(1, 2))
         ijpair = 0
      end if
      if (.not. allocated(ijktrio)) then
         allocate (ijktrio(1, 6))
         ijktrio = 0.0_rp
      end if

      open (newunit=funit, file=fname_, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//fmt('A', trim(fname_))//'not found', __FILE__, __LINE__)
      end if

      read (funit, nml=lattice, iostat=iostatus)

      if (allocated(ct)) then
         allocate (ct_first(size(ct)))
         ct_first = ct
      end if
      if (allocated(screening_alpha)) then
         allocate (screening_alpha_first(size(screening_alpha)))
         screening_alpha_first = screening_alpha
      end if
      if (allocated(screening_sigma)) then
         allocate (screening_sigma_first(size(screening_sigma)))
         screening_sigma_first = screening_sigma
      end if

      if (size(izp) .ne. ndim) then
         deallocate (izp)
         allocate (izp(ndim))
      end if
      if (size(no) .ne. ndim) then
         deallocate (no)
         allocate (no(ndim))
      end if
      if (size(crd) .ne. 3*ndim) then
         deallocate (crd)
         allocate (crd(3, ndim))
      end if
      if (size(inclu) .ne. 3*nclu) then
         deallocate (inclu)
         allocate (inclu(nclu, 3))
      end if

      if (size(ijpair) .ne. 2*njij) then
         deallocate (ijpair)
         allocate (ijpair(njij, 2))
      end if

      if (size(ct) .ne. ntype) then
         deallocate (ct)
         allocate (ct(ntype))
      end if
      if (size(screening_alpha) < 4) then
         deallocate (screening_alpha)
         allocate (screening_alpha(4))
         screening_alpha = this%default_screening_alpha(size(screening_alpha))
      end if
      if (size(screening_sigma) < 4) then
         deallocate (screening_sigma)
         allocate (screening_sigma(4))
         screening_sigma = 0.7_rp
      end if

      if (size(ijktrio) .ne. 2*njijk) then
         deallocate (ijktrio)
         allocate (ijktrio(njijk, 6))
      end if

      rewind (funit)
      read (funit, nml=lattice, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//fmt('I0', iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      if (allocated(ct_first)) then
         if (allocated(ct)) then
            if (all(abs(ct) <= tiny(1.0_rp)) .and. any(abs(ct_first) > tiny(1.0_rp))) then
               ct(1:min(size(ct), size(ct_first))) = ct_first(1:min(size(ct), size(ct_first)))
            end if
         end if
         deallocate (ct_first)
      end if
      if (allocated(screening_alpha_first)) then
         if (allocated(screening_alpha)) then
            if (all(abs(screening_alpha) <= tiny(1.0_rp)) .and. any(abs(screening_alpha_first) > tiny(1.0_rp))) then
               screening_alpha(1:min(size(screening_alpha), size(screening_alpha_first))) = &
                  screening_alpha_first(1:min(size(screening_alpha), size(screening_alpha_first)))
            end if
         end if
         deallocate (screening_alpha_first)
      end if
      if (allocated(screening_sigma_first)) then
         if (allocated(screening_sigma)) then
            if (all(abs(screening_sigma) <= tiny(1.0_rp)) .and. any(abs(screening_sigma_first) > tiny(1.0_rp))) then
               screening_sigma(1:min(size(screening_sigma), size(screening_sigma_first))) = &
                  screening_sigma_first(1:min(size(screening_sigma), size(screening_sigma_first)))
            end if
         end if
         deallocate (screening_sigma_first)
      end if

      ! General intialization

      this%r2 = r2
      this%alat = alat
      this%celldm = celldm
      this%ntype = ntype
      this%strux_backend = trim(adjustl(strux_backend))
      this%screening = trim(adjustl(screening))
      this%strux_want_sdot = strux_want_sdot
      this%strux_solve_scale = strux_solve_scale
      if (len_trim(this%strux_backend) == 0) this%strux_backend = 'legacy'
      if (len_trim(this%screening) == 0) this%screening = 'default'
      if (this%strux_solve_scale <= 0.0_rp) this%strux_solve_scale = 9.0_rp
      this%pbc = pbc
      this%b1 = b1
      this%b2 = b2
      this%b3 = b3
      this%n1 = n1
      this%n2 = n2
      this%n3 = n3

      call move_alloc(ct, this%ct)
      if (allocated(this%ct)) then
         if (size(this%ct) > 0) then
            if (all(abs(this%ct) <= tiny(1.0_rp)) .and. this%r2 > 0.0_rp) then
               this%ct = sqrt(this%r2)
            end if
         end if
      end if
      call move_alloc(screening_alpha, this%screening_alpha)
      call move_alloc(screening_sigma, this%screening_sigma)
      ! Reads Wigner-Seitz radius if available
      this%wav = wav
      ! Bulk initialization
      this%ndim = ndim
      this%npe = npe
      this%rc = rc
      this%crystal_sym = crystal_sym
      this%a = a
      call move_alloc(izp, this%izp)
      call move_alloc(no, this%no)
      call move_alloc(crd, this%crd)

      ! Impurity initialization
      this%nclu = nclu
      call move_alloc(inclu, this%inclu)

      ! Surface initialization
      this%surftype = surftype
      this%nlay = nlay

      ! Exchange calculation initialization
      call move_alloc(ijpair, this%ijpair)
      this%njij = njij

      ! Spin-lattice calculation initialization
      call move_alloc(ijktrio, this%ijktrio)
      this%njijk = njijk

      ! Condition to store the trios in the original njij type

      if ((njijk .ne. 0) .and. (njij .eq. 0)) then
         this%njij = 3*this%njijk
#ifdef USE_SAFE_ALLOC
         if (allocated(this%ijpair)) call g_safe_alloc%deallocate('lattice.ijpair', this%ijpair)
         call g_safe_alloc%allocate('lattice.ijpair', this%ijpair, (/this%njij, 2/))
#else
         if (allocated(this%ijpair)) deallocate (this%ijpair)
         allocate (this%ijpair(this%njij, 2))
#endif
         do i = 1, this%njijk
            this%ijpair(3*(i - 1) + 1, 1) = int(this%ijktrio(i, 1))
            this%ijpair(3*(i - 1) + 1, 2) = int(this%ijktrio(i, 2))
            this%ijpair(3*(i - 1) + 2, 1) = int(this%ijktrio(i, 1))
            this%ijpair(3*(i - 1) + 2, 2) = int(this%ijktrio(i, 3))
            this%ijpair(3*(i - 1) + 3, 1) = int(this%ijktrio(i, 2))
            this%ijpair(3*(i - 1) + 3, 2) = int(this%ijktrio(i, 3))
         end do
      else if ((this%njijk .ne. 0) .and. (this%njij .ne. 0)) then
         call g_logger%fatal('not possible to calculate Jij and Jijk at the same time', __FILE__, __LINE__)
      end if

      ! checks
      call this%check_all()

   end subroutine build_from_file

   subroutine build_from_lattice(this)
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
      namelist /lattice/ r2, nbulk_bulk, ntot, nbas, nrec, &
         a, crd, &
         ct, izp, no, iu, ib, irec, strux_backend, screening, strux_want_sdot, &
         strux_solve_scale, screening_alpha, screening_sigma

      strux_backend = this%strux_backend
      screening = this%screening
      strux_want_sdot = this%strux_want_sdot
      strux_solve_scale = this%strux_solve_scale

      call move_alloc(this%crd, crd)
      call move_alloc(this%izp, izp)
      call move_alloc(this%no, no)
      call move_alloc(this%ct, ct)
      call move_alloc(this%screening_alpha, screening_alpha)
      call move_alloc(this%screening_sigma, screening_sigma)

      open (newunit=funit, file='lattice.nml', action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file lattice.nml not found', __FILE__, __LINE__)
      end if

      ! Pre-allocate ib, iu, irec before the sizing read to avoid undefined
      ! behaviour on compilers that write into unallocated storage when these
      ! variables appear in the namelist file. crd was moved from this%crd and
      ! has size (3, ntot), so its second dimension is a safe upper bound.
      nsite_guess = 1
      if (allocated(crd)) then
         nsite_guess = max(1, size(crd, 2))
         allocate (ib(size(crd, 2)), iu(size(crd, 2)), irec(size(crd, 2)))
      else if (allocated(izp)) then
         nsite_guess = max(1, size(izp))
         allocate (ib(nsite_guess), iu(nsite_guess), irec(nsite_guess))
      else
         allocate (ib(1), iu(1), irec(1))
      end if
      if (.not. allocated(ct)) then
         allocate (ct(nsite_guess))
         ct = 0.0_rp
      else if (size(ct) == 0) then
         deallocate (ct)
         allocate (ct(nsite_guess))
         ct = 0.0_rp
      end if
      if (.not. allocated(screening_alpha)) then
         allocate (screening_alpha(4))
         screening_alpha = this%default_screening_alpha(size(screening_alpha))
      else if (size(screening_alpha) < 4) then
         deallocate (screening_alpha)
         allocate (screening_alpha(4))
         screening_alpha = this%default_screening_alpha(size(screening_alpha))
      end if
      if (.not. allocated(screening_sigma)) then
         allocate (screening_sigma(4))
         screening_sigma = 0.7_rp
      else if (size(screening_sigma) < 4) then
         deallocate (screening_sigma)
         allocate (screening_sigma(4))
         screening_sigma = 0.7_rp
      end if
      read (funit, nml=lattice, iostat=iostatus)
      if (allocated(ct)) then
         allocate (ct_first(size(ct)))
         ct_first = ct
      end if
      if (allocated(screening_alpha)) then
         allocate (screening_alpha_first(size(screening_alpha)))
         screening_alpha_first = screening_alpha
      end if
      if (allocated(screening_sigma)) then
         allocate (screening_sigma_first(size(screening_sigma)))
         screening_sigma_first = screening_sigma
      end if
      deallocate (ib, iu, irec)
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('lattice.ib', ib, (/ntot/))
      call g_safe_alloc%allocate('lattice.iu', iu, (/ntot/))
      call g_safe_alloc%allocate('lattice.irec', irec, (/nrec/))
#else
      allocate (ib(ntot), iu(ntot), irec(nrec))
#endif
      if (allocated(izp)) nsite_guess = max(1, size(izp))
      if (size(ct) /= nsite_guess) then
         deallocate (ct)
         allocate (ct(nsite_guess))
         ct = 0.0_rp
      end if
      if (size(screening_alpha) < 4) then
         deallocate (screening_alpha)
         allocate (screening_alpha(4))
         screening_alpha = this%default_screening_alpha(size(screening_alpha))
      end if
      if (size(screening_sigma) < 4) then
         deallocate (screening_sigma)
         allocate (screening_sigma(4))
         screening_sigma = 0.7_rp
      end if
      rewind (funit)
      read (funit, nml=lattice, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//fmt('I0', iostatus), __FILE__, __LINE__)
      end if
      close (funit)
      if (allocated(ct_first)) then
         if (all(abs(ct) <= tiny(1.0_rp)) .and. any(abs(ct_first) > tiny(1.0_rp))) then
            ct(1:min(size(ct), size(ct_first))) = ct_first(1:min(size(ct), size(ct_first)))
         end if
         deallocate (ct_first)
      end if
      if (allocated(screening_alpha_first)) then
         if (all(abs(screening_alpha) <= tiny(1.0_rp)) .and. any(abs(screening_alpha_first) > tiny(1.0_rp))) then
            screening_alpha(1:min(size(screening_alpha), size(screening_alpha_first))) = &
               screening_alpha_first(1:min(size(screening_alpha), size(screening_alpha_first)))
         end if
         deallocate (screening_alpha_first)
      end if
      if (allocated(screening_sigma_first)) then
         if (all(abs(screening_sigma) <= tiny(1.0_rp)) .and. any(abs(screening_sigma_first) > tiny(1.0_rp))) then
            screening_sigma(1:min(size(screening_sigma), size(screening_sigma_first))) = &
               screening_sigma_first(1:min(size(screening_sigma), size(screening_sigma_first)))
         end if
         deallocate (screening_sigma_first)
      end if

      this%nbulk_bulk = nbulk_bulk
      this%ntot = ntot
      !this%ntype = ntype
      this%nbas = nbas
      this%nrec = nrec
      !this%r2 = r2
      this%nbulk = 0
      this%strux_backend = trim(adjustl(strux_backend))
      this%screening = trim(adjustl(screening))
      this%strux_want_sdot = strux_want_sdot
      this%strux_solve_scale = strux_solve_scale
      if (len_trim(this%strux_backend) == 0) this%strux_backend = 'legacy'
      if (len_trim(this%screening) == 0) this%screening = 'default'
      if (this%strux_solve_scale <= 0.0_rp) this%strux_solve_scale = 9.0_rp

      call move_alloc(ct, this%ct)
      if (allocated(this%ct)) then
         if (size(this%ct) > 0) then
            if (all(abs(this%ct) <= tiny(1.0_rp)) .and. this%r2 > 0.0_rp) then
               this%ct = sqrt(this%r2)
            end if
         end if
      end if
      call move_alloc(screening_alpha, this%screening_alpha)
      call move_alloc(screening_sigma, this%screening_sigma)
      call move_alloc(izp, this%izp)
      call move_alloc(no, this%no)
      call move_alloc(ib, this%ib)
      call move_alloc(irec, this%irec)
      call move_alloc(iu, this%iu)
      call move_alloc(crd, this%crd)
      this%a = a
   end subroutine build_from_lattice
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build structure according to crystal symetry provided
   !---------------------------------------------------------------------------
   subroutine build_data(this)
      class(lattice), intent(inout) :: this
      !> Local variables
      real(rp), dimension(3, 3) :: a
      integer :: i, j

      select case (this%crystal_sym)
      case ('bcc')
         a(:, 1) = [-0.50000000, 0.50000000, 0.50000000]
         a(:, 2) = [0.50000000, -0.50000000, 0.50000000]
         a(:, 3) = [0.50000000, 0.50000000, -0.50000000]
         this%crd(:, 1) = [0.0, 0.0, 0.0]
         this%izp(1) = 1
         this%no(1) = 1
         this%nbulk_bulk = 1
         this%ntot = 1
         if (this%control%calctype == 'B') then
            this%nrec = this%nbulk_bulk
            this%nbulk = 0
            this%ntype = 1
            this%nbas = this%ntot
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
!            call g_safe_alloc%allocate('lattice.ct', this%ct, (/this%ntype/))
#else
            allocate (this%ib(this%ntot), this%iu(this%ntot), this%irec(this%nrec))!, this%ct(this%ntype)) ! Now ct is defined at &lattice
#endif
            this%iu(1) = 1
            this%ib(1) = 1
            this%irec(1) = 1
!            this%ct(:) = this%alat + 0.1d0
         end if
         this%a = a
      case ('b2')
         a(:, 1) = [1.00000000, 0.00000000, 0.00000000]
         a(:, 2) = [0.00000000, 1.00000000, 0.00000000]
         a(:, 3) = [0.00000000, 0.00000000, 1.00000000]
         this%crd(:, 1) = [0.0, 0.0, 0.0]
         this%crd(:, 2) = [0.5, 0.5, 0.5]
         this%izp(1) = 1
         this%izp(2) = 2
         this%no(1) = 1
         this%no(2) = 2
         this%nbulk_bulk = 2
         this%ntot = 2
         if (this%control%calctype == 'B') then
            this%nrec = this%nbulk_bulk
            this%nbulk = 0
            this%ntype = 2
            this%nbas = this%ntot
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
!            call g_safe_alloc%allocate('lattice.ct', this%ct, (/this%ntype/))
#else
            allocate (this%ib(this%ntot), this%iu(this%ntot), this%irec(this%nrec)) !, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
            this%iu(1) = 1
            this%iu(2) = 2
            this%ib(1) = 1
            this%ib(2) = 2
            this%irec(1) = 1
            this%irec(2) = 2
!            this%ct(:) = this%alat + 0.3d0
         end if
         this%a = a
      case ('fcc')
         a(:, 1) = [0.00000000, 0.50000000, 0.50000000]
         a(:, 2) = [0.50000000, 0.00000000, 0.50000000]
         a(:, 3) = [0.50000000, 0.50000000, 0.00000000]
         this%crd(:, 1) = [0.0, 0.0, 0.0]
         this%izp(1) = 1
         this%no(1) = 1
         this%nbulk_bulk = 1
         this%ntot = 1
         if (this%control%calctype == 'B') then
            this%nrec = this%nbulk_bulk
            this%nbulk = 0
            this%ntype = 1
            this%nbas = this%ntot
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
!            call g_safe_alloc%allocate('lattice.ct', this%ct, (/this%ntype/))
#else
            allocate (this%ib(this%ntot), this%iu(this%ntot), this%irec(this%nrec))!, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
            this%iu(1) = 1
            this%ib(1) = 1
            this%irec(1) = 1
!            this%ct(:) = this%alat + 0.1d0
         end if
         this%a = a
      case ('fcc2')
         a(:, 1) = [0.00000000, 0.50000000, 0.50000000]
         a(:, 2) = [0.50000000, 0.00000000, 0.50000000]
         a(:, 3) = [0.50000000, 0.50000000, 0.00000000]
         this%crd(:, 1) = [0.00, 0.00, 0.00]
         this%crd(:, 2) = [0.50, 0.50, 0.50]
         this%izp(:) = [1, 2]
         this%no(:) = [1, 2]
         this%nbulk_bulk = 2
         this%ntot = 2
         if (this%control%calctype == 'B') then
            this%nrec = this%nbulk_bulk
            this%nbulk = 0
            this%ntype = 2
            this%nbas = this%ntot
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
!            call g_safe_alloc%allocate('lattice.ct', this%ct, (/this%ntype/))
#else
            allocate (this%ib(this%ntot), this%iu(this%ntot), this%irec(this%nrec))!, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
            this%iu(:) = [1, 2]
            this%ib(:) = [1, 2]
            this%irec(:) = [1, 2]
!            this%ct(:) = this%alat - 1.5d0
!            this%r2 = this%ct(1)**2
         end if
         this%a = a
      case ('fcc3')
         a(:, 1) = [0.50000000, 0.50000000, 0.00000000]
         a(:, 2) = [0.50000000, 0.00000000, 0.50000000]
         a(:, 3) = [0.00000000, 0.50000000, 0.50000000]
         this%crd(:, 1) = [0.25000000, 0.25000000, 0.25000000]
         this%crd(:, 2) = [0.00000000, 0.00000000, 0.00000000]
         this%crd(:, 3) = [0.50000000, 0.50000000, 0.50000000]
         this%crd(:, 4) = [-0.25000000, -0.25000000, -0.25000000]
         this%izp(1:4) = [1, 2, 3, 4]
         this%no(1:4) = [1, 2, 3, 4]
         this%nbulk_bulk = 4
         this%ntot = 4
         if (this%control%calctype == 'B') then
            this%nrec = this%nbulk_bulk
            this%nbulk = 0
            this%ntype = 4
            this%nbas = this%ntot
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
!            call g_safe_alloc%allocate('lattice.ct', this%ct, (/this%ntype/))
#else
            allocate (this%ib(this%ntot), this%iu(this%ntot), this%irec(this%nrec)) !, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
            this%iu(1:4) = [1, 2, 3, 4]
            this%ib(1:4) = [1, 2, 3, 4]
            this%irec(1:4) = [1, 2, 3, 4]
!            this%ct(1:4) = this%alat - 1.5d0
!            this%r2 = this%ct(1)**2
         end if
         this%a = a
      case ('hcp')
         if (this%celldm == 0.0d0) then
            print *, 'WARNING: hcp structure using as default c/a = 1.633'
            this%celldm = 1.633d0
         end if
         a(:, 1) = [1.00000000, 0.00000000, 0.00000000]
         a(:, 2) = [-.50000000, 0.86602500, 0.00000000]
         a(:, 3) = [0.00000000, 0.00000000, 0.00000000]
         a(3, 3) = this%celldm
         this%crd(:, 1) = [0.0, 0.0, 0.0]
         this%crd(:, 2) = [0.0, 0.577350000, 0.00000000]
         this%crd(3, 2) = (0.5d0)*this%celldm
         this%nbulk_bulk = 2
         this%ntot = 2
         this%izp(:) = [1, 2]
         this%no(:) = [1, 2]
         if (this%control%calctype == 'B') then
            this%nrec = this%nbulk_bulk
            this%nbulk = 0
            this%ntype = 2
            this%nbas = this%ntot
#ifdef USE_SAFE_ALLOC
            call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
            call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
!            call g_safe_alloc%allocate('lattice.ct', this%ct, (/this%ntype/))
#else
            allocate (this%ib(this%ntot), this%iu(this%ntot), this%irec(this%nrec)) !, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
            this%iu(:) = [1, 2]
            this%ib(:) = [1, 2]
            this%irec(:) = [1, 2]
!            this%ct(:) = this%alat + 0.1d0
         end if
         this%a = a
      case ('file')
         ! reads from input
         call this%build_from_lattice()
      end select

      ! Ensure neighbour-cutoff data is valid for built-in crystal presets.
      ! For legacy bulk/surface inputs, ct is often omitted and should default
      ! to include first + second shells (historical behaviour: alat + 0.1).
      if (this%crystal_sym /= 'file') then
         if (.not. allocated(this%ct) .or. size(this%ct) /= this%ntype) then
            if (allocated(this%ct)) deallocate(this%ct)
            allocate(this%ct(this%ntype))
            this%ct(:) = 0.0_rp
         end if
         if (all(abs(this%ct) <= tiny(1.0_rp))) then
            this%ct(:) = this%alat + 0.1_rp
         end if
         if (this%r2 <= 0.0_rp) then
            this%r2 = this%ct(1)**2
         end if
      end if

      ! Volume in cubic Angstroms
      this%vol = abs(dot_product(this%a(:, 3), cross_product(this%a(:, 1), this%a(:, 2))))*(this%alat**3)
      if (this%wav .eq. 0) then
         this%wav = (this%vol/((16.0d0/3.0d0)*atan(1.0d0)*this%ntot))**(1.0d0/3.0d0)
         write (*, *) 'wav', this%wav
      end if
      if (this%control%calctype == 'B' .or. this%control%calctype == 'S') this%nmax = 0
   end subroutine build_data

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this, full)
      class(lattice) :: this
      logical, intent(in), optional :: full

      this%ndim = 9900000
      this%npe = 49
      this%nclu = 0
      this%surftype = 'none'
      this%nlay = 0
      this%ntype = 0
      this%nbas = 0
      this%wav = 0
      this%celldm = 0.0d0
      this%njij = 0
      this%njijk = 0
      this%strux_backend = 'legacy'
      this%screening = 'default'
      this%strux_want_sdot = .false.
      this%strux_solve_scale = 9.0_rp
      this%nn_max = 0
      this%a_cart_inv = 0.0_rp
      this%a_cart_inv_ready = .false.
      this%b1 = .false.
      this%b2 = .false.
      this%b3 = .false.
      this%pbc = .false.
      this%n1 = 0
      this%n2 = 0
      this%n3 = 0
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('lattice.izp', this%izp, (/this%ndim/))
      call g_safe_alloc%allocate('lattice.no', this%no, (/this%ndim/))
      call g_safe_alloc%allocate('lattice.crd', this%crd, (/3, this%ndim/))
      call g_safe_alloc%allocate('lattice.inclu', this%inclu, (/this%nclu, 3/))
      call g_safe_alloc%allocate('lattice.ijpair', this%ijpair, (/this%njij, 2/))
      call g_safe_alloc%allocate('lattice.ijktrio', this%ijktrio, (/this%njijk, 6/))
      call g_safe_alloc%allocate('lattice.chargetrf_type', this%chargetrf_type, (/this%nbas/))
      call g_safe_alloc%allocate('lattice.ct', this%ct, (/this%ntype/))
      call g_safe_alloc%allocate('lattice.screening_alpha', this%screening_alpha, (/max(1, this%control%lmax + 1)/))
      call g_safe_alloc%allocate('lattice.screening_sigma', this%screening_sigma, (/max(1, this%control%lmax + 1)/))
#else
      allocate (this%izp(this%ndim), this%no(this%ndim))
      allocate (this%crd(3, this%ndim))
      allocate (this%inclu(this%nclu, 3))
      allocate (this%ijpair(this%njij, 2))
      allocate (this%ijktrio(this%njijk, 6))
      allocate (this%chargetrf_type(this%nbas))
      allocate (this%ct(this%ntype))
      allocate (this%screening_alpha(max(1, this%control%lmax + 1)))
      allocate (this%screening_sigma(max(1, this%control%lmax + 1)))
#endif

      this%izp = 0.0d0
      this%no = 0
      this%crd = 0.d0
      this%inclu = 0.0d0
      this%ijpair = 0.0d0
      this%ijktrio = 0.0d0
      this%chargetrf_type = 0.0d0
      this%ct = 0.0d0
      this%screening_alpha = this%default_screening_alpha(size(this%screening_alpha))
      this%screening_sigma = 0.7_rp
      if (associated(this%control)) then
         if (present(full)) then
            if (full) then
               call this%control%restore_to_default()
            end if
         end if
      end if

   end subroutine restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !---------------------------------------------------------------------------
   subroutine bravais(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      real(rp) :: rc, rs, lc, lcx, lcy, lcz
      integer, dimension(:), allocatable :: iz, num
      real(rp), dimension(:, :), allocatable :: cr, crbravais
   real(rp), allocatable :: tmp(:, :)
      integer :: npe, ndim, nx, ny, nz, npr, l, n, i, nl, k, kk
      logical :: isopen
      integer :: iostatus

      inquire (unit=10, opened=isopen)
      if (isopen) then
         call g_logger%fatal('lattice%bravais, file clust: Unit 10 is already open', __FILE__, __LINE__)
      else
         open (unit=10, file='clust')
      end if

      n = this%ntot
      ! Defining a size-like parameter for the clust before the cut
      npe = this%npe
      ndim = this%ndim
      rc = this%rc
      allocate (iz(ndim), num(ndim))
      allocate (cr(3, ndim), crbravais(3, ndim))
      ! Defining the cut radius
      rs = (0.8d0*int((npe/(1.0d0))/2.0d0))**2
      rs = dmin1(rs, rc)

      npr = int((ndim/(n*1.0d0))**(1.0d0/3.0d0))

      crbravais(:, :) = this%crd(:, :)
      if (this%pbc) then
         lcx = (this%n1 + 1)/2
         lcy = (this%n2 + 1)/2
         lcz = (this%n3 + 1)/2
      else
         lc = (npr + 1)/2
         lcx = lc
         lcy = lc
         lcz = lc
      end if
      l = n
      if (.not.this%pbc) then
         this%n1 = npr; this%n2 = npr; this%n3 = npr
      end if

      do i = 1, l
         do nx = 1, this%n1        
            do ny = 1, this%n2 
               do nz = 1, this%n3  
                  if (nx .eq. lcx .and. ny .eq. lcy .and. nz .eq. lcz) go to 13
                  n = n + 1
!  ...........Verifies dimension NDIM.........................
                  if (n .gt. ndim) then
                     write (6, 600) n
                     stop
                  end if
! ...............................................................
                  crbravais(1, n) = this%crd(1, i) + (nx - lcx)*this%a(1, 1) + (ny - lcy)*this%a(1, 2) + (nz - lcz)*this%a(1, 3)
                  crbravais(2, n) = this%crd(2, i) + (nx - lcx)*this%a(2, 1) + (ny - lcy)*this%a(2, 2) + (nz - lcz)*this%a(2, 3)
                  crbravais(3, n) = this%crd(3, i) + (nx - lcx)*this%a(3, 1) + (ny - lcy)*this%a(3, 2) + (nz - lcz)*this%a(3, 3)
                  this%no(n) = this%no(i)
                  this%izp(n) = this%izp(i)
13                continue
               end do
            end do
         end do
      end do
      nl = l*(npr**3)

      kk = 0
      if (rc == 0.0d0) rs = npr**3

      if (this%pbc) then
         ndim = this%n1*this%n2*this%n3*l 
         iz = this%izp
         num = this%no
         cr = crbravais
         kk = n
      else
          do i = 1, nl
             call cut(i, l, ndim, crbravais, cr, this%izp, iz, num, this%no, rs, kk)
          end do
      end if 

      if (int(kk/2) /= kk/2.d0) kk = kk - 1

      write (10, 300) kk
      do k = 1, kk, 2
         write (10, 200) (cr(1, k + i - 1), cr(2, k + i - 1), cr(3, k + i - 1), iz(k + i - 1), num(k + i - 1), i=1, 2)
      end do
      write (800, *) kk
      do k = 1, kk
         write (800, '(i5,3f14.8)') iz(k), cr(:, k)
      end do

200   format(3F14.8, 2I4, 3F14.8, 2I4)
300   format(3X, "II =", I7)
600   format(1X, 'K =', I7, 5X, 'redefine the dimension ndim')

      this%kk = kk
      call move_alloc(cr, this%cr)
      ! Shrink this%cr to the actual cluster size kk to avoid keeping the
      ! large temporary allocation of shape (3, ndim). We copy the first kk
      ! columns into a smaller array and move that allocation into this%cr.
      if (allocated(this%cr)) then
         if (size(this%cr, 2) > kk) then
            allocate(tmp(3, kk))
            tmp(:, :) = this%cr(:, 1:kk)
            call move_alloc(tmp, this%cr)
         end if
      end if
      call move_alloc(iz, this%iz)
      call move_alloc(num, this%num)
      close (10)

      ! Test to set nmax to the whole cluster
      ! this%nmax = kk
   end subroutine bravais

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !---------------------------------------------------------------------------
   subroutine build_clusup(this)
      class(lattice), intent(inout) :: this
      character(len=10) :: surftype
      ! Local variables
      real(rp), dimension(:), allocatable :: z
      integer, dimension(:), allocatable :: izpsurf
      real(rp) :: ds, ds2, new, one
      integer :: i, j, n

      select case (this%crystal_sym)
      case ('hcp')
         read (this%surftype, *) this%dx, this%dy, this%dz, this%dw
         if ((this%dz) .ne. (-1*(this%dx + this%dy))) then
            this%dz = -1*(this%dx + this%dy)
            call g_logger%fatal('Hexagonal Miller indices not right. Does it should be ['// &
                                fmt('I0', this%dx)//fmt('I0', this%dy)//fmt('I0', this%dz)//fmt('I0', this%dw)//']?', &
                                __FILE__, __LINE__)
         end if
         this%dx = (2*this%dx + this%dy)
         this%dy = (this%dx + 2*this%dy)
         this%dz = this%dw
      case default
         read (this%surftype, *) this%dx, this%dy, this%dz
      end select

!  ............Find zstep, zmin, zmax....................
      ds = 1000.0d0; ds2 = 1000.0d0
      select case (this%crystal_sym)
      case ('hcp')
         do i = 1, this%kk
            do j = 1, this%kk
               new = this%dx*this%cr(1, i) + this%dy*this%cr(2, i) + this%dz*this%cr(3, i)
               one = this%dx*this%cr(2, j) + this%dy*this%cr(2, j) + this%dz*this%cr(3, j)
               if ((abs(new - one) .gt. 1.d-6) .and. (abs(new - one) .le. ds)) then
                  ds = abs(new - one)
               end if
               if ((abs(new) .le. ds2)) then
                  ds2 = abs(new)
               end if
            end do
         end do
      case default
         do i = 1, this%kk
            new = this%dx*this%cr(1, i) + this%dy*this%cr(2, i) + this%dz*this%cr(3, i)
            if ((abs(new) .gt. 1.d-6) .and. (abs(new) .le. ds)) then
               ds = abs(new)
            end if
            if ((abs(new) .le. ds2)) then
               ds2 = abs(new)
            end if
         end do
      end select

      this%zstep = ds
      this%zmin = ds2 - this%zstep
      this%zmax = ds2 + 15*this%zstep

!  ......................................................

      n = int((this%zmax - this%zmin)/this%zstep) + 1

      allocate (z(n), izpsurf(n))
      do i = 1, n
         z(i) = this%zmin + ((i - 1)*this%zstep)
      end do
      call move_alloc(z, this%z)

      do i = 1, n
         izpsurf(i) = mod(i, this%ntot) + 1
      end do

      j = this%ntot
      do i = 1, this%nlay
         j = j + 1
         izpsurf(i) = j
      end do
      call move_alloc(izpsurf, this%izpsurf)
      if (this%control%calctype == 'S') then
         this%ntype = this%nbulk_bulk + this%nlay
         this%nbulk = this%nbulk_bulk
         this%nrec = this%ntype - this%nbulk
         this%nbas = 51
#ifdef USE_SAFE_ALLOC
         call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%nbulk/))
         call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
         call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
         call g_safe_alloc%allocate('lattice.ct', this%ct, (/this%ntype/))
#else
         allocate (this%ib(this%nbulk), this%irec(this%nrec), this%iu(this%ntot)) !, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
!         this%ct(:) = this%alat + 0.1d0
!         this%r2 = this%ct(1)**2
      end if

      this%surftype = clean_str(this%surftype)
   end subroutine build_clusup

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Cuts the cluster in the desirable direction and builds the surface cluster
   !---------------------------------------------------------------------------
   subroutine build_surf_full(this)
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
      natoms = this%kk

      ! Open the file and read header
      inquire (unit=10, opened=isopen)
      if (isopen) then
         call g_logger%fatal('lattice%build_surf, file clust: Unit 10 is already open', __FILE__, __LINE__)
      else
         open (unit=10, file='clust')
      end if
      ! Open the output for surface information
      inquire (unit=20, opened=isopen)
      if (isopen) then
         call g_logger%fatal('lattice%build_surf, file clust: Unit 10 is already open', __FILE__, __LINE__)
      else
         open (unit=20, file='surfclu.out')
      end if

      select case (this%crystal_sym)
      case ('hcp')
         read (this%surftype, *) this%dx, this%dy, this%dz, this%dw
         if ((this%dz) .ne. (-1*(this%dx + this%dy))) then
            this%dz = -1*(this%dx + this%dy)
            call g_logger%fatal('Hexagonal Miller indices not right. Does it should be ['// &
                                fmt('I0', this%dx)//fmt('I0', this%dy)//fmt('I0', this%dz)//fmt('I0', this%dw)//']?', &
                                __FILE__, __LINE__)
         end if
         this%dx = (2*this%dx + this%dy)
         this%dy = (this%dx + 2*this%dy)
         this%dz = this%dw
      case default
         read (this%surftype, *) this%dx, this%dy, this%dz
      end select

      ! Allocate arrays to store atom details
      allocate (crsurf(3, this%kk), crh(this%kk), typesurf(this%kk), crystalsurf(this%kk))
      allocate (uniqueTypes(this%kk), ichoicen(this%kk), ichoicetypen(this%kk), crhd(this%kk))

      ds = 1000.0d0
      ds2 = 1000.0d0
      do i = 1, this%kk
         do j = 1, this%kk
            new = this%dx*this%cr(1, i) + this%dy*this%cr(2, i) + this%dz*this%cr(3, i)
            one = this%dx*this%cr(1, j) + this%dy*this%cr(2, j) + this%dz*this%cr(3, j)
            if ((abs(new - one) .gt. 1.d-6) .and. (abs(new - one) .le. ds)) then
               ds = abs(new - one)
            end if
            if ((abs(new) .le. ds2)) then
               ds2 = abs(new)
            end if
         end do
      end do

      this%zstep = ds
      this%zmin = ds2 - this%zstep
      this%zmax = ds2 + 20*this%zstep

      n = int((this%zmax - this%zmin)/this%zstep) + 1

      ! Allocate the number of layers
      allocate (z(n))
      ! Allocate the number of atoms per layer variable
      allocate (nTypesForCurrentLayer(n))

      do i = 1, n
         z(i) = this%zmin + ((i - 1)*this%zstep)
      end do
      call move_alloc(z, this%z)

      ! Determining max atom type from input
      maxType = maxval(this%iz(:))
      call move_alloc(this%iz, atomType)
      call move_alloc(this%num, crystalType)

      nsurf = 0
      do i = 1, n
         nTypesForCurrentLayer(i) = 0
         disi_min = sqrt(this%z(i)**2) + 1.0d0
         do k = 1, natoms
            crh(k) = 0.0d0
            crh(k) = dot_product([this%dx, this%dy, this%dz], this%cr(:, k))
            if (abs(crh(k) - this%z(i)) < 1.0d-6) then
               nsurf = nsurf + 1
               crsurf(:, nsurf) = this%cr(:, k)
               crhd(nsurf) = dot_product([this%dx, this%dy, this%dz], crsurf(:, nsurf))
               if (i <= this%nlay) then
                  ! Determine the unique atom types and update typesurf and crystalsurf.
                  isUnique = .true.
                  do j = 1, nTypesForCurrentLayer(i)
                     if (atomType(k) == uniqueTypes(j)) then
                        isUnique = .false.
                        exit
                     end if
                  end do

                  if (isUnique) then
                     nTypesForCurrentLayer(i) = nTypesForCurrentLayer(i) + 1
                     uniqueTypes(nTypesForCurrentLayer(i)) = atomType(k)
                     maxType = maxType + 1
                     typesurf(nsurf) = maxType
                     crystalsurf(nsurf) = crystalType(k) !maxType
                     disi = norm2(crsurf(:, nsurf))
                     if (disi < disi_min) then
                        ichoicen(typesurf(nsurf)) = nsurf
                        ichoicetypen(typesurf(nsurf)) = typesurf(nsurf)
                     end if
                  else
                     typesurf(nsurf) = maxType - nTypesForCurrentLayer(i) + findloc(uniqueTypes, atomType(k), dim=1)
                     crystalsurf(nsurf) = crystalType(k)
                     disi = norm2(crsurf(:, nsurf))
                     if (disi < disi_min) then
                        ichoicen(typesurf(nsurf)) = nsurf
                        ichoicetypen(typesurf(nsurf)) = typesurf(nsurf)
                     end if
                  end if
               else
                  typesurf(nsurf) = atomType(k)
                  crystalsurf(nsurf) = crystalType(k)
                  disi = norm2(crsurf(:, nsurf))
                  if (i .le. (this%nlay + this%nbulk_bulk)) then
                     if (disi < disi_min) then
                        ichoicen(typesurf(nsurf)) = nsurf
                        ichoicetypen(typesurf(nsurf)) = typesurf(nsurf)
                     end if
                  end if
               end if
            end if
         end do
      end do

      ! Passing onto the global variable
      call move_alloc(nTypesForCurrentLayer, this%natoms_layer)
      write (20, *) 'Maxtype:', maxType
      write (20, *) 'Type of atoms chosen:', ichoicetypen(1:maxType)
      write (20, *) 'Atoms chosen:', ichoicen(1:maxType)
      write (20, *) 'Unique atoms type per layer:', this%natoms_layer(:)
      write (20, *) 'Layers are:', this%z
      write (20, *) 'Step is:', this%zstep, this%zmax, this%zmin
      if (this%control%calctype == 'S') then
         this%ntype = maxType
         this%nbulk = this%nbulk_bulk
         this%nrec = this%ntype - this%nbulk
         this%nbas = 49

         if (allocated(this%chargetrf_type)) deallocate (this%chargetrf_type)
         if (allocated(this%ib)) deallocate (this%ib)
         if (allocated(this%irec)) deallocate (this%irec)
         if (allocated(this%iu)) deallocate (this%iu)
!         if (allocated(this%ct)) deallocate (this%ct)

         allocate (this%ib(this%nbulk), this%irec(this%nrec), this%iu(this%ntot))!, this%ct(this%ntype)) Now ct is defined at &lattice
         allocate (this%chargetrf_type(this%nbas))

!         this%ct(:) = 4.0d0 !this%alat + 0.1d0
!         this%r2 = this%ct(1)**2

         do i = 1, this%nrec
            this%irec(i) = ichoicen(this%nbulk + i)
         end do
         do i = 1, this%ntot
            this%iu(ichoicetypen(i)) = ichoicen(i)
         end do
         do i = 1, this%nbulk
            this%ib(ichoicetypen(i)) = ichoicen(i)
         end do
      end if

      if (int(nsurf/2) /= nsurf/2.d0) nsurf = nsurf - 1

      ! Rotate the surface cluster so that z is perpendicular to the plane
      ! Calculate the unit vector normal to the crystallographic plane
      !axis(:) = [this%dx,this%dy,this%dz]

      !theta = acos(dot_product(axis,[0.0_rp,0.0_rp,1.0_rp]))

      !axis(:) = cross_product(axis(:),[0.0_rp,0.0_rp,1.0_rp])

      !norm = norm2(axis)

      !axis = axis / norm

    !! Setup the rotation matrix using Rodrigues´ formula
      !call setup_rotation_matrix(rotation_matrix, axis, theta)

      !do i = 1, nsurf
      !   rotated_cr(:,i) = matmul(rotation_matrix(:,:), crsurf(:,i))
      !   crsurf(:,i) = rotated_cr(:,i)
      !end do

      !rotated_cr(:,:) = 0.0_rp
    !! Rotate the primitive lattice vectors
      !do i=1,3
      !   rotated_cr(:,i) = matmul(rotation_matrix(:,:), this%a(:,i))
      !   write(*,*) rotated_cr(:,i)
      !end do

      write (10, 10004) nsurf
      do k = 1, nsurf - 1, 2
         write (10, 10002) &
            crsurf(1, k), crsurf(2, k), crsurf(3, k), typesurf(k), crystalsurf(k), crsurf(1, k + 1), crsurf(2, k + 1), &
            crsurf(3, k + 1), typesurf(k + 1), crystalsurf(k + 1)
      end do

      do i = 1, maxtype
         write (20, '(3f12.6, 2i5, f12.6)') crsurf(:, ichoicen(i)), ichoicetypen(i), crystalsurf(ichoicen(i)), &
                    & dot_product([this%dx, this%dy, this%dz], crsurf(:, ichoicen(i)))   
      end do

      ! Cleanup
      deallocate (crh, uniqueTypes)
      deallocate (this%cr)

      call move_alloc(crsurf, this%cr)
      call move_alloc(typesurf, this%iz)
      call move_alloc(crystalsurf, this%num)

      this%kk = nsurf
    !!!! FIX LATER
      if (this%control%calctype == 'I') this%nlay = maxtype - this%nbulk_bulk
    !!!!!
10002 format(3f14.6, 2i4, 3f14.6, 2i4)
10004 format(3x, "II =", i7)
      close (10)
      close (20)
   end subroutine build_surf_full
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Builds the surface cluster
   !---------------------------------------------------------------------------
   subroutine build_surf(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      integer :: i, j, n, k, kk, ichoice, icont, ichoicetype
      real(rp), dimension(this%kk) :: crh, crhd
      real(rp) :: disf, disi_min, disi
      real(rp), dimension(:, :), allocatable :: crsurf
      integer, dimension(:), allocatable :: izp, no, ichoicen, ichoicetypen
      logical :: isopen

      n = int((this%zmax - this%zmin)/this%zstep) + 1
      kk = this%kk
      allocate (crsurf(3, kk), izp(kk), no(kk), ichoicen(kk), ichoicetypen(kk))
      icont = 0
      crh = 0.0d0
      crhd = 0.0d0
      inquire (unit=10, opened=isopen)
      if (isopen) then
         call g_logger%fatal('lattice%build_surf, file clust: Unit 10 is already open', __FILE__, __LINE__)
      else
         open (unit=10, file='clust')
      end if

      do i = 1, n
         disf = 1.0d0
         disi_min = sqrt(this%z(i)**2) + 0.5d0
         do k = 1, kk
            crh(k) = 0.0d0
            crh(k) = this%dx*this%cr(1, k) + this%dy*this%cr(2, k) + this%dz*this%cr(3, k)
            if (abs(crh(k) - this%z(i)) < 1.0d-6) then
               icont = icont + 1
               do j = 1, 3
                  crsurf(j, icont) = this%cr(j, k)
               end do
               crhd(icont) = this%dx*crsurf(1, icont) + this%dy*crsurf(2, icont) + this%dz*crsurf(3, icont)
               if (abs(crhd(icont) - this%z(i)) < 1.0d-6) then
                  izp(icont) = this%izpsurf(i)
                  disi = norm2(crsurf(:, icont))
                  if (disi < disi_min) then
                     disi_min = disi
                     ichoice = icont
                     ichoicetype = izp(icont)
                  end if
                  no(icont) = this%num(k)
               else
               end if
            end if
         end do
         ichoicetypen(i) = ichoicetype
         ichoicen(i) = ichoice
      end do

      if (this%control%calctype == 'S') then
         do i = 1, this%nrec
            this%irec(i) = ichoicen(i)
         end do
         do i = 1, this%ntot
            this%iu(ichoicetypen(i + this%nrec)) = ichoicen(i + this%nrec)
         end do
         do i = 1, this%nbulk
            this%ib(ichoicetypen(i + this%nrec)) = ichoicen(i + this%nrec)
         end do
      end if

      if (int(icont/2) /= icont/2.d0) icont = icont - 1

      do k = 1, icont
         write (805, *) no(k), crsurf(:, k)
      end do
      write (10, 10004) icont
      do k = 1, icont - 1, 2
         write (10, 10002) &
            crsurf(1, k), crsurf(2, k), crsurf(3, k), izp(k), no(k), crsurf(1, k + 1), crsurf(2, k + 1), &
            crsurf(3, k + 1), izp(k + 1), no(k + 1)
      end do

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%deallocate('lattice.cr', this%cr)
      call g_safe_alloc%deallocate('lattice.iz', this%iz)
      call g_safe_alloc%deallocate('lattice.num', this%num)
#else
      deallocate (this%cr, this%iz, this%num)
#endif

      call move_alloc(crsurf, this%cr)
      call move_alloc(izp, this%iz)
      call move_alloc(no, this%num)

      this%kk = icont

10002 format(3f14.6, 2i4, 3f14.6, 2i4)
10004 format(3x, "II =", i7)
      close (10)
   end subroutine build_surf

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Adds the cluster/impurity atoms in ´inclu´ to ´clu0´ and outputs to ´clust´
   !---------------------------------------------------------------------------
   subroutine newclu(this)
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
      this%nbulk = this%nbulk_bulk + this%nlay
      this%ntype = this%nbulk_bulk + this%nlay + this%nclu
      this%nrec = this%ntype - this%nbulk
      write (*, *) this%nbulk, this%ntype, this%nrec
      ! Set the clust dimension
      kk = this%kk
      ndi = 150000
      nnmx = 200

      ! Allocating clust dimension
      if (allocated(this%ib)) deallocate (this%ib)
      if (allocated(this%irec)) deallocate (this%irec)
      if (allocated(this%iu)) deallocate (this%iu)
      !if (allocated(this%ct)) deallocate (this%ct)

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('lattice.ib', this%ib, (/this%nbulk/))
      call g_safe_alloc%allocate('lattice.iu', this%iu, (/this%ntot/))
      call g_safe_alloc%allocate('lattice.irec', this%irec, (/this%nrec/))
!      call g_safe_alloc%allocate('lattice.ct', this%ct, (/this%ntype/))
#else
      allocate (this%ib(this%nbulk), this%iu(this%ntot), this%irec(this%nrec))!, this%ct(this%ntype)) Now ct is defined at &lattice
#endif
      allocate (ctnew(this%ntype), ibulk(this%nbulk))
      allocate (nn(ndi, nnmx), nn2(ndi, nnmx))
      allocate (izpo(kk), izp(kk), no(kk), nnmax(kk), izimp(kk), noimp(kk))
      allocate (acr(kk, 7), crd(3, kk), crimp(3, kk))
      ! Setting ct values for impurity
      !this%ct(:) = this%alat+0.1d0
      ! Identify impurity atoms from ´inclu´
      do i = 1, kk
         izpo(i) = this%iz(i)
      end do
      ntypecount = this%nbulk
      inclucheck = 0
      do j = 1, this%nclu
         do i = 1, kk
            if (abs(this%cr(1, i) - this%inclu(j, 1)) < 1.0d-6 .and. &
                abs(this%cr(2, i) - this%inclu(j, 2)) < 1.0d-6 .and. &
                abs(this%cr(3, i) - this%inclu(j, 3)) < 1.0d-6) then
               inclucheck = inclucheck + 1
               ntypecount = ntypecount + 1
               this%iz(i) = ntypecount
            end if
         end do
      end do

      if (inclucheck /= this%nclu) then
         call g_logger%fatal('Atoms chosen for impurity were not found inside the clust. Please, check the inclu input', __FILE__, __LINE__)
      end if

      acr = 0.0d0
      inquire (unit=10, opened=isopen)
      if (isopen) then
         call g_logger%fatal('lattice%newclu, file clust: Unit 10 is already open', __FILE__, __LINE__)
      else
         open (unit=10, file='clust')
      end if

      open (unit=11, file='outnewclu')

      do k = 1, kk
         do i = 1, 3
            acr(k, i) = this%cr(i, k)
            acr(k, 6) = acr(k, 6) + (this%cr(i, k) - this%inclu(1, i))**2
         end do
         acr(k, 4) = (this%iz(k))
         acr(k, 5) = (this%num(k))
         acr(k, 7) = (izpo(k))
      end do
      call bubble(this%nclu, this%nclu, acr(1:this%nclu, 1:7), 7, 4)
      call bubble(kk - this%nclu, kk - this%nclu, acr(this%nclu + 1:kk, 1:7), 7, 6)

      write (10, 10004) kk
      do k = 1, kk - 1, 2
         write (10, 10002) &
            ( &
            acr(k + i - 1, 1), acr(k + i - 1, 2), acr(k + i - 1, 3), int(acr(k + i - 1, 4)) &
            , int(acr(k + i - 1, 5)), i=1, 2)
      end do

      rewind (10)

#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

      call leia(this%alat, kk, crd, izp, no, 10)

      close (10)

      nnscale = 0.95d0

      ctnew(:) = this%ct(:)*nnscale

      call this%nncal(ctnew, crd, 3, kk, izp, nn, ndi, nnmx, mapa, this%ntot)
      nnmx = 200
      call this%nncal(this%ct, crd, 3, kk, izp, nn2, ndi, nnmx, mapa, this%ntot)
      nnmx = 200

      do i = 1, kk
         if (acr(i, 4) > this%nbulk .and. acr(i, 4) <= this%ntype) then
            do j = 2, nn2(i, 1)
               if (acr(nn2(i, j), 4) <= this%nbulk) acr(nn2(i, j), 4) = 2000 + acr(i, 7)
            end do
         end if
      end do
      do i = 1, kk
         if (acr(i, 4) > this%nbulk .and. acr(I, 4) <= this%ntype) then
            do j = 2, nn(i, 1)
               if (acr(nn(i, j), 4) <= this%nbulk .or. acr(nn(i, j), 4) > 2000) acr(nn(i, j), 4) = 1000 + acr(i, 7)
            end do
         end if
      end do
      do i = 1, kk
         if (acr(i, 4) == 1) acr(i, 4) = 4000 + acr(I, 7)
      end do
      do i = 1, kk
         if (acr(i, 4) <= this%nbulk) acr(i, 4) = 3000 + acr(I, 7)
      end do
      call bubble(kk, kk, acr(1:kk, 1:7), 7, 4)
      ncnt = 0
      do I = 1, kk
         if (acr(i, 4) < 2000) ncnt = ncnt + 1
      end do
      do i = 1, kk
         if (acr(i, 4) > this%ntype) acr(i, 4) = acr(i, 7)
      end do
      call bubble(kk - ncnt, kk - ncnt, acr(ncnt + 1:(kk), 1:7), 7, 6)

      open (unit=10, file='clust')
      write (10, 10004) kk
      do k = 1, kk, 2
         write (10, 10002) &
            ( &
            acr(k + i - 1, 1), acr(k + i - 1, 2), acr(k + i - 1, 3), int(acr(k + i - 1, 4)) &
            , int(acr(k + i - 1, 5)), i=1, 2)
      end do
      rewind (10)
#ifdef USE_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
      call leia(this%alat, kk, crd, izp, no, 10)
      call this%nncal(this%ct, crd, 3, kk, izp, nn, ndi, nnmx, mapa, this%ntot)
      call outmap(11, izp, nn, no, ndi, nnmx, this%nrec)
      close (10)
      nmax = 0
      do i = 1, this%nrec
         do j = 2, nn(i, 1)
            if (nmax < nn(i, j)) nmax = nn(i, j)
         end do
      end do
      write (11, *) "--Info-for-bulcri------------------"
      !
      if (allocated(this%chargetrf_type)) deallocate (this%chargetrf_type)
      allocate (this%chargetrf_type(ncnt))
      !
      do k = 1, ncnt
         write (11, '(i5)') int(acr(k, 7))   !, int(ACR(K, 7))
         this%chargetrf_type(k) = int(acr(k, 7))
      end do
      this%nbas = ncnt
      this%nmax = nmax
      ! Temporary hack for full HALL
      !this%nmax = kk
      write (11, *) "--Info-for-control-----------------"
      write (11, '(a6, i4, a6, i6, a6, i6, a6, i6)') 'NTYPE=', this%ntype, 'NMAX=', this%nmax, 'NBAS=', this%nbas, 'NREC=', this%nrec

      nnmax = 0
      do i = nmax + 1, kk
         if (nnmax(izp(i)) < nn(i, 1)) then
            nnmax(izp(i)) = nn(i, 1)
            ibulk(izp(i)) = i
         end if
      end do
      write (11, '(50i7)') (nnmax(k), k=1, this%nbulk)
      write (11, '(50i7)') (ibulk(k), k=1, this%nbulk)
      do i = 1, this%nbulk
         do j = 2, nn(ibulk(i), 1)
            if (nn(ibulk(i), j) == 0) write (*, *) "Warning! Atom", ibulk(i), " is probably wrong."
         end do
      end do

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%deallocate('lattice.cr', this%cr)
      call g_safe_alloc%deallocate('lattice.cr', this%cr)
      call g_safe_alloc%deallocate('lattice.iz', this%iz)
      call g_safe_alloc%deallocate('lattice.num', this%num)
#else
      deallocate (this%cr, this%iz, this%num)
#endif

      do k = 1, kk
         crimp(1:3, k) = acr(k, 1:3)
      end do

      izimp(:) = int(acr(:, 4))
      noimp(:) = int(acr(:, 5))

      call move_alloc(izimp, this%iz)
      call move_alloc(noimp, this%num)
      call move_alloc(crimp, this%cr)

      do i = 1, this%nbulk
         this%ib(i) = ibulk(i)
      end do

      do i = 1, this%ntot
         this%iu(i) = ibulk(i)
      end do

      ireccount = 0
      do j = 1, this%nclu
         do i = 1, kk
            if (abs(this%cr(1, i) - this%inclu(j, 1)) < 1.0d-6 .and. &
                abs(this%cr(2, i) - this%inclu(j, 2)) < 1.0d-6 .and. &
                abs(this%cr(3, i) - this%inclu(j, 3)) < 1.0d-6) then
               ireccount = ireccount + 1
               this%irec(ireccount) = i
            end if
         end do
      end do

      close (11)
10002 format(3f14.6, 2i4, 3f14.6, 2i4)
10004 format(3x, "II =", i7)
   end subroutine newclu

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the tight-binding structure constant matrix elements from
   !> information in ´control´ and ´clust.
   !---------------------------------------------------------------------------
   subroutine structb(this, do_str)
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

      ! Open files
      open (12, file='map', form='unformatted')
      open (13, file="sbar", FORM="unformatted")
      open (16, file="view.sbar")
      if (this%strux_want_sdot) then
         open (14, file="sdot", form="unformatted")
         open (15, file="view.sdot")
      end if
      open (17, file='str.out')

      ! Clust parameters
      ncut = 9
      nnmx = 100 ! 5250
      nomx = this%ntot
      kk = this%kk
      allocate (set(3, nomx, nnmx)); set = 0.0d0
      allocate (nn(kk, nnmx))
      allocate (idnn(nnmx))
      nm = nnmx
      write (17, *) 'irec', this%nrec, this%irec
      write (17, *) 'irec type', this%iz(this%irec(:))
      write (17, *) 'ndi=', kk
      write (17, 10000) kk
      write (17, 10001)
      write (17, 10002) (i, (this%cr(j, i)*this%alat, j=1, 3), i=1, max(this%nmax, this%ntype))
      ! write (*, '(a,10f10.4)') 'NNCAL ct=', this%ct(1:min(this%ntype, size(this%ct)))
      call this%nncal(this%ct, this%cr*this%alat, 3, kk, this%iz, nn, kk, nm, mapa, this%ntype)

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('lattice.nn', this%nn, (/this%kk, nm + 1/))
#else
      allocate (this%nn(this%kk, nm + 1))
#endif
      do ii = 1, nm + 1
         this%nn(:, ii) = nn(:, ii)
      end do
      write (17, *) 'ndi=', kk
      write (17, *) 'remd'
      call this%remd(this%cr*this%alat, this%num, this%iu, this%nn, kk, this%ntot, nomx, kk, nnmx, set, idnn, ret)
      nm_store = max(1, maxval(this%nn(:, 1)))
      do ii = 1, this%ntot
         ia = this%iu(ii)
         if (ia <= 0) cycle
         nt_tmp = this%kk
         call this%clusba(this%r2, this%cr*this%alat, ia, kk, kk, nt_tmp)
         nm_store = max(nm_store, nt_tmp)
      end do
      this%nn_max = nm_store
      sbar_dim = max(norb, this%control%npold)
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('lattice.sbar', this%sbar, (/sbar_dim, sbar_dim, nm_store, this%ntot/))
#else
      allocate (this%sbar(sbar_dim, sbar_dim, nm_store, this%ntot))
#endif
      call this%init_strux_storage(sbar_dim, nm_store)
      write (17, *) 'outmap', this%nmax, maxval(this%irec)
      call outmap(17, this%iz, this%nn, this%num, kk, nnmx, max(this%nmax, maxval(this%irec)))
      write (17, 10003) kk, nm
      flush(17)
      if (do_str) then
         if (rank==0) call g_logger%info('STRUCTB structure-constants, backend='//trim(this%strux_backend), __FILE__, __LINE__)
         if (this%control%lmax > 2 .and. trim(lower(this%strux_backend)) /= 'strux_lib') then
            call g_logger%fatal('lmax='//int2str(this%control%lmax)//' requires strux_lib backend for structure constants; legacy backend is only valid up to spd (lmax=2).', __FILE__, __LINE__)
         end if
         if (trim(lower(this%strux_backend)) == 'strux_lib') then
            call this%structb_strux()
         else
            do ii = 1, this%ntot
               ia = this%iu(ii)
               nr = this%nn(ia, 1)
               write (17, '(1x, a, i5, a, i5)') 'Sbar atom no:', ii, ' Ntot:', this%ntot
               call this%dbar1(ia, ncut*this%r2, this%wav, this%cr*this%alat, kk, kk, norb, nr, ii)
            end do
         end if
      end if
10000 format(i5)
10001 format(" LATTICE COORDINATES")
10002 format(2(i5, 3f8.4))
10003 format(3i5)
10004 format(3i5)
10005 format(7x, i7)
   end subroutine structb

   subroutine init_strux_storage(this, sbar_dim, nm)
      class(lattice), intent(inout) :: this
      integer, intent(in) :: sbar_dim, nm

#ifdef USE_SAFE_ALLOC
      if (allocated(this%sdot)) call g_safe_alloc%deallocate('lattice.sdot', this%sdot)
      call g_safe_alloc%allocate('lattice.sdot', this%sdot, (/sbar_dim, sbar_dim, nm, this%ntot/))
      if (allocated(this%alpha)) call g_safe_alloc%deallocate('lattice.alpha', this%alpha)
      call g_safe_alloc%allocate('lattice.alpha', this%alpha, (/sbar_dim, this%kk, this%ntot/))
      if (allocated(this%alpha_dot)) call g_safe_alloc%deallocate('lattice.alpha_dot', this%alpha_dot)
      call g_safe_alloc%allocate('lattice.alpha_dot', this%alpha_dot, (/sbar_dim, this%kk, this%ntot/))
#else
      if (allocated(this%sdot)) deallocate (this%sdot)
      allocate (this%sdot(sbar_dim, sbar_dim, nm, this%ntot))
      if (allocated(this%alpha)) deallocate (this%alpha)
      allocate (this%alpha(sbar_dim, this%kk, this%ntot))
      if (allocated(this%alpha_dot)) deallocate (this%alpha_dot)
      allocate (this%alpha_dot(sbar_dim, this%kk, this%ntot))
#endif

      this%sbar(:, :, :, :) = czero
      this%sdot(:, :, :, :) = czero
      this%alpha(:, :, :) = 0.0_rp
      this%alpha_dot(:, :, :) = 0.0_rp
   end subroutine init_strux_storage

   pure function default_screening_alpha(this, nl) result(alpha_default)
      class(lattice), intent(in) :: this
      integer, intent(in) :: nl
      real(rp) :: alpha_default(nl)
      real(rp), parameter :: default_values(4) = [0.3485_rp, 0.0530_rp, 0.0107_rp, 0.00674_rp]
      integer :: i

      alpha_default = 0.0_rp
      do i = 1, nl
         alpha_default(i) = default_values(min(i, size(default_values)))
      end do
   end function default_screening_alpha

   integer function strux_mode(this)
      class(lattice), intent(in) :: this
      character(len=:), allocatable :: screening_mode

      screening_mode = trim(lower(this%screening))
      select case (screening_mode)
      case ('manual', 'default')
         strux_mode = STRUX_LMTO47_IALPHA_MANUAL
      case ('sigma')
         strux_mode = STRUX_LMTO47_IALPHA_SIGMA
      case ('fitted')
         strux_mode = STRUX_LMTO47_IALPHA_FITD
      case default
         call g_logger%fatal('Unsupported lattice.screening='//trim(this%screening), __FILE__, __LINE__)
      end select
   end function strux_mode

   subroutine load_symbolic_atoms_if_needed(this)
      class(lattice), intent(inout) :: this

      if (allocated(this%symbolic_atoms)) then
         if (size(this%symbolic_atoms) == this%ntype) return
         deallocate(this%symbolic_atoms)
      end if
      if (.not. allocated(this%symbolic_atoms)) then
         this%symbolic_atoms = array_of_symbolic_atoms(this%control%fname, this%ntype)
      end if
   end subroutine load_symbolic_atoms_if_needed

   subroutine build_rmt(this, nspec, species_labels, rmt)
      class(lattice), intent(inout) :: this
      integer, intent(in) :: nspec
      integer, intent(in) :: species_labels(nspec)
      real(rp), intent(out) :: rmt(nspec)
      integer :: is, label

      call this%load_symbolic_atoms_if_needed()

      do is = 1, nspec
         label = species_labels(is)
         if (label < 1 .or. label > size(this%symbolic_atoms)) then
            call g_logger%fatal('strux backend found invalid species label in lattice%no for WS radius lookup', __FILE__, __LINE__)
         end if
         if (this%symbolic_atoms(label)%potential%ws_r <= tiny(1.0_rp)) then
            call g_logger%fatal('strux backend requires positive potential%ws_r for sigma/fitted screening', __FILE__, __LINE__)
         end if
         ! NOTE: potential%ws_r is already stored in Bohr in this code path.
         ! Keep rmt in the same unit to remain consistent with avw_bohr.
         rmt(is) = this%symbolic_atoms(label)%potential%ws_r
      end do
   end subroutine build_rmt

   subroutine build_hcr(this, nl, nspec, rmt, hcr)
      class(lattice), intent(in) :: this
      integer, intent(in) :: nl, nspec
      real(rp), intent(in) :: rmt(nspec)
      real(rp), intent(out) :: hcr(nl, nspec)
      integer :: is, l

      do is = 1, nspec
         do l = 1, nl
            hcr(l, is) = this%screening_sigma(l)*rmt(is)
         end do
      end do
   end subroutine build_hcr

   subroutine build_strux_inputs(this, nspec, nl, species_labels, alpha_in, hcr, rmt)
      class(lattice), intent(inout) :: this
      integer, intent(in) :: nspec, nl
      integer, intent(in) :: species_labels(nspec)
      real(rp), intent(out) :: alpha_in(0:nl - 1, nspec)
      real(rp), intent(out) :: hcr(nl, nspec)
      real(rp), intent(out) :: rmt(nspec)
      integer :: is
      real(rp) :: alpha_global(nl)

      call this%build_rmt(nspec, species_labels, rmt)
      call this%build_hcr(nl, nspec, rmt, hcr)

      select case (trim(lower(this%screening)))
      case ('manual')
         alpha_global = this%screening_alpha(1:nl)
      case ('default')
         alpha_global = this%default_screening_alpha(nl)
      case default
         alpha_global = 0.0_rp
      end select

      alpha_in(:, :) = 0.0_rp
      do is = 1, nspec
         alpha_in(:, is) = alpha_global
      end do
   end subroutine build_strux_inputs

   subroutine structb_strux(this)
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

      nl = this%control%lmax + 1
      nl2 = nl*nl
      nbas = this%nbas
      if (nbas <= 0) nbas = max(1, maxval(this%no(1:min(this%kk, size(this%no)))))
      sbar_dim = size(this%sbar, 1)
      pair_cutoff = sqrt(this%r2)
      solve_cutoff = sqrt(max(this%strux_solve_scale, 1.0_rp)*this%r2)
      pair_cutoff_bohr = pair_cutoff*ang2au
      solve_cutoff_bohr = solve_cutoff*ang2au
      alat_bohr = this%alat*ang2au
      wav_bohr = this%wav*ang2au

      if (nl2 > sbar_dim) then
         call g_logger%fatal('strux basis size exceeds allocated sbar dimension', __FILE__, __LINE__)
      end if
      if (nbas <= 0) then
         call g_logger%fatal('strux backend requires lattice basis coordinates', __FILE__, __LINE__)
      end if
      if (.not. allocated(this%crd) .or. .not. allocated(this%no)) then
         call g_logger%fatal('strux backend requires primitive-cell coordinates and labels', __FILE__, __LINE__)
      end if
      if (size(this%crd, 2) < nbas .or. size(this%no) < nbas) then
         call g_logger%fatal('strux inferred basis size exceeds primitive-cell storage', __FILE__, __LINE__)
      end if
      if (any(this%no(1:nbas) <= 0)) then
         call g_logger%fatal('strux backend found non-positive primitive basis labels in lattice%no', __FILE__, __LINE__)
      end if

      allocate(orb_map(max_orb))
      call build_orbital_map(sbar_dim, orb_map)

      opts%method = STRUX_METHOD_LMTO47
      opts%auto_alpha = .false.
      opts%want_sdot = this%strux_want_sdot
      opts%lmaxw = -1
      opts%pair_cutoff = pair_cutoff_bohr
      opts%solve_cutoff = solve_cutoff_bohr
      effective_screening = trim(lower(this%screening))
      select case (effective_screening)
      case ('manual', 'default')
         opts%screening_mode = STRUX_LMTO47_IALPHA_MANUAL
      case ('sigma')
         opts%screening_mode = STRUX_LMTO47_IALPHA_SIGMA
      case ('fitted')
         opts%screening_mode = STRUX_LMTO47_IALPHA_FITD
      case default
         call g_logger%fatal('Unsupported effective screening mode='//trim(effective_screening), __FILE__, __LINE__)
      end select

      if (effective_screening == 'sigma' .or. effective_screening == 'fitted') then
         opts%auto_alpha = .true.
      end if

      allocate(species_labels(nbas))
      species_labels = 0
      nspec = 0
      allocate(pos(3, nbas), ips(nbas))
      allocate(cralat(3, this%kk))

      pos(:, :) = this%crd(:, 1:nbas)
      cralat(:, :) = this%cr(:, 1:this%kk)*this%alat
      do ib = 1, nbas
         label = this%no(ib)
         if (label <= 0) then
            call g_logger%fatal('strux backend found non-positive primitive basis labels in lattice%no', __FILE__, __LINE__)
         end if
         species_idx = 0
         do is = 1, nspec
            if (species_labels(is) == label) then
               species_idx = is
               exit
            end if
         end do
         if (species_idx == 0) then
            nspec = nspec + 1
            species_labels(nspec) = label
            species_idx = nspec
         end if
         ips(ib) = species_idx
      end do
      allocate(lmxb(nspec), rmt(nspec))
      allocate(alpha_in(0:nl - 1, nspec), hcr(nl, nspec))
      lmxb(:) = this%control%lmax
      call this%build_strux_inputs(nspec, nl, species_labels(1:nspec), alpha_in, hcr, rmt)
      write (17, '(a, i6, a, i6, a, f10.4, a, f10.4)') 'STRUX periodic solve nbas=', nbas, ' kk=', this%kk, ' pair_cutoff=', pair_cutoff, ' solve_cutoff=', solve_cutoff

      ! DEBUG_STRUX_DIAG_BEGIN (temporary instrumentation; safe to remove as one block)
      write (17, '(a)') 'DEBUG_STRUX_DIAG: begin'
      write (17, '(a,a)') 'DEBUG_STRUX_DIAG: effective_screening=', trim(effective_screening)
      write (17, '(a,i0)') 'DEBUG_STRUX_DIAG: opts%screening_mode=', opts%screening_mode
      write (17, '(a,l1)') 'DEBUG_STRUX_DIAG: opts%auto_alpha=', opts%auto_alpha
      write (17, '(a,f12.6,a,f12.6)') 'DEBUG_STRUX_DIAG: hcr[min,max]=', minval(hcr), ',', maxval(hcr)
      write (17, '(a,f12.6,a,f12.6)') 'DEBUG_STRUX_DIAG: alpha_in[min,max]=', minval(alpha_in), ',', maxval(alpha_in)
      write (17, '(a,f12.6,a,f12.6)') 'DEBUG_STRUX_DIAG: rmt[min,max]=', minval(rmt), ',', maxval(rmt)
      write (17, '(a,f12.6,a,f12.6)') 'DEBUG_STRUX_DIAG: avw_bohr, solve_cutoff_bohr=', wav_bohr, ',', solve_cutoff_bohr
      write (17, '(a)') 'DEBUG_STRUX_DIAG: end'
      flush(17)
      if (rank == 0) then
         call g_logger%info('DEBUG_STRUX_DIAG effective_screening='//trim(effective_screening)// &
                            ' mode='//int2str(opts%screening_mode)// &
                            ' auto='//merge('T','F',opts%auto_alpha), __FILE__, __LINE__)
      end if
      ! DEBUG_STRUX_DIAG_END

      if (effective_screening == 'manual' .or. effective_screening == 'default') then
         if (this%strux_want_sdot) then
            call g_logger%fatal('strux manual/default screening with sdot is unsupported; use sigma or fitted', __FILE__, __LINE__)
         end if
      end if

      if (allocated(result%iax)) deallocate(result%iax)
      if (allocated(result%alpha)) deallocate(result%alpha)
      if (allocated(result%alpha_l)) deallocate(result%alpha_l)
      if (allocated(result%s)) deallocate(result%s)
      if (allocated(result%sdot)) deallocate(result%sdot)

      select case (effective_screening)
      case ('manual', 'default')
         call strux_compute(opts, nbas, nspec, nl, alat_bohr, this%a, pos, ips, lmxb, wav_bohr, rmt, result, alpha_in=alpha_in)
         call strux_lmto47_screening(nbas, nspec, nl, wav_bohr, ips, lmxb, rmt, alpha_in, alpha_site, adot_site, &
            tral, trad, screening_mode=opts%screening_mode)
      case ('sigma')
         call strux_compute(opts, nbas, nspec, nl, alat_bohr, this%a, pos, ips, lmxb, wav_bohr, rmt, result, hcr=hcr)
         call strux_lmto47_autoalpha_screening(nbas, nspec, nl, wav_bohr, ips, lmxb, rmt, hcr, &
            alpha_l_out=alpha_l_out, alpha_out=alpha_site, adot=adot_site, tral=tral, trad=trad, &
            screening_mode=opts%screening_mode)
      case ('fitted')
         call strux_compute(opts, nbas, nspec, nl, alat_bohr, this%a, pos, ips, lmxb, wav_bohr, rmt, result)
         call strux_lmto47_autoalpha_screening(nbas, nspec, nl, wav_bohr, ips, lmxb, rmt, hcr, &
            alpha_l_out=alpha_l_out, alpha_out=alpha_site, adot=adot_site, tral=tral, trad=trad, &
            screening_mode=opts%screening_mode)
      end select
      do is = 1, nspec
         alpha_debug = 0.0_rp
         select case (effective_screening)
         case ('manual', 'default')
            alpha_debug(1:min(nl, 4)) = alpha_in(0:min(nl - 1, 3), is)
         case ('sigma', 'fitted')
            if (allocated(alpha_l_out)) then
               alpha_debug(1:min(nl, 4)) = alpha_l_out(0:min(nl - 1, 3), is)
            else if (allocated(result%alpha_l)) then
               alpha_debug(1:min(nl, 4)) = result%alpha_l(0:min(nl - 1, 3), is)
            end if
         end select
         write (17, '(a, a, a, i4, a, i6, a, 4f12.6)') 'STRUX screening=', trim(effective_screening), &
            ' species=', is, ' label=', species_labels(is), ' alpha=', alpha_debug
      end do
      flush(17)

      nttab = result%nttab
      write (17, '(a, i8, a, f10.3)') 'STRUX periodic result nttab=', nttab, ' cpu_s=', t_total_end - t_total_start

      do ii = 1, this%ntot
         this%alpha(:, :, ii) = 0.0_rp
         this%alpha_dot(:, :, ii) = 0.0_rp
         this%alpha(1:nl2, 1:nbas, ii) = alpha_site(:, 1:nbas)
         this%alpha_dot(1:nl2, 1:nbas, ii) = adot_site(:, 1:nbas)
      end do

      do is = 1, nspec
         label = species_labels(is)
         if (label < 1 .or. label > size(this%symbolic_atoms)) then
            call g_logger%fatal('strux backend found invalid species label for screening alpha storage', __FILE__, __LINE__)
         end if
         if (allocated(this%symbolic_atoms(label)%potential%screening_alpha)) then
            if (lbound(this%symbolic_atoms(label)%potential%screening_alpha, 1) /= 0 .or. &
                ubound(this%symbolic_atoms(label)%potential%screening_alpha, 1) /= nl - 1) then
               deallocate(this%symbolic_atoms(label)%potential%screening_alpha)
               allocate(this%symbolic_atoms(label)%potential%screening_alpha(0:nl - 1))
            end if
         else
            allocate(this%symbolic_atoms(label)%potential%screening_alpha(0:nl - 1))
         end if
         select case (effective_screening)
         case ('manual', 'default')
            this%symbolic_atoms(label)%potential%screening_alpha(0:nl - 1) = alpha_in(0:nl - 1, is)
         case ('sigma', 'fitted')
            if (allocated(alpha_l_out)) then
               this%symbolic_atoms(label)%potential%screening_alpha(0:nl - 1) = alpha_l_out(0:nl - 1, is)
            else if (allocated(result%alpha_l)) then
               this%symbolic_atoms(label)%potential%screening_alpha(0:nl - 1) = result%alpha_l(0:nl - 1, is)
            else
               call g_logger%fatal('strux backend did not produce screening alpha data for potential transform', __FILE__, __LINE__)
            end if
         end select
      end do
      do ii = 1, this%ntot
         ia = representative_atom_index(this, ii)
         ib = primitive_basis_label(this, ia)
         nt_store = this%kk
         call this%clusba(this%r2, cralat, ia, this%kk, this%kk, nt_store)
         call cpu_time(t_map_start)
         write (17, '(a, i5, a, i5, a, i5, a, i6)') 'STRUX map   center ', ii, ' atom ', ia, ' basis=', ib, ' nt=', nt_store
         call write_neighbor_vector_dump(17, this%sbarvec, nt_store)

         do m = 1, nt_store
            vec_target(:) = this%sbarvec(:, m)
            if (m == 1) then
               ja = ia
               jb = ib
            else
               ja = this%nn(ia, m)
               if (ja == 0) then
                  ja = find_neighbor_atom_by_vector(this, ia, vec_target, cralat, match_tol)
                  if (ja == 0) cycle
               end if
               jb = primitive_basis_label(this, ja)
            end if
            pair_idx = find_pair_by_vector(nttab, result%iax, this%a, pos, this%alat, ib, jb, vec_target, match_tol)
            if (pair_idx == 0) then
               write (17, '(a,3f14.8,a,2i6,a,2i6)') 'STRUX pair miss vec=', vec_target, ' center/neigh=', ia, ja, ' basis=', ib, jb
               call g_logger%fatal('Failed to map strux periodic pair to nn ordering', __FILE__, __LINE__)
            end if

            do is = 1, nl2
               if (orb_map(is) <= 0) cycle
               do js = 1, nl2
                  if (orb_map(js) <= 0) cycle
                  this%sbar(is, js, m, ii) = cmplx(result%s(orb_map(is), orb_map(js), pair_idx), 0.0_rp, kind=rp)
                  if (this%strux_want_sdot) then
                     this%sdot(is, js, m, ii) = cmplx(result%sdot(orb_map(is), orb_map(js), pair_idx), 0.0_rp, kind=rp)
                  end if
               end do
            end do
            call write_strux_block(this, nl2, m, ia, ja)
         end do
         call cpu_time(t_map_end)
         write (17, '(a, i5, a, f10.3, a, i6)') 'STRUX done  center ', ii, ' cpu_s=', t_map_end - t_map_start, ' stored_nn=', nt_store
         flush(17)
      end do

      if (allocated(alpha_l_out)) deallocate(alpha_l_out)
      if (allocated(tral)) deallocate(tral)
      if (allocated(trad)) deallocate(trad)
      if (allocated(alpha_site)) deallocate(alpha_site)
      if (allocated(adot_site)) deallocate(adot_site)
      deallocate(pos, cralat, ips, lmxb, rmt, alpha_in, hcr, orb_map, species_labels)
   end subroutine structb_strux

   subroutine write_strux_block(this, nl2, m, ii, jclus)
      class(lattice), intent(in) :: this
      integer, intent(in) :: nl2, m, ii, jclus
      integer :: is, js

      ! write (*, '(" SBAR neighbor center=",i5," slot=",i5," iclus=",i5," vec=",3f12.6)') ii, m, iclus, &
      !    this%sbarvec(1, m), this%sbarvec(2, m), this%sbarvec(3, m)
      write (16, '(" VECTOR=",3f12.6,"   ICLUS=",i5," JCLUS=",i5)') this%sbarvec(1, m), this%sbarvec(2, m), this%sbarvec(3, m), ii, jclus
      if (this%strux_want_sdot) then
         write (15, '(" VECTOR=",3f12.6,"   ICLUS=",i5," JCLUS=",i5)') this%sbarvec(1, m), this%sbarvec(2, m), this%sbarvec(3, m), ii, jclus
      end if

      do is = 1, nl2
         do js = 1, nl2
            write (13) this%sbar(is, js, m, ii)
            if (this%strux_want_sdot) write (14) this%sdot(is, js, m, ii)
         end do
         write (16, '(*(f10.4))') (real(this%sbar(is, js, m, ii), rp), js=1, nl2)
         if (this%strux_want_sdot) then
            write (15, '(*(f10.4))') (real(this%sdot(is, js, m, ii), rp), js=1, nl2)
         end if
      end do
      flush(16)
      flush(17)
      if (this%strux_want_sdot) flush(15)
   end subroutine write_strux_block

   subroutine write_neighbor_vector_dump(iunit, sbarvec, nt)
      integer, intent(in) :: iunit, nt
      real(rp), intent(in) :: sbarvec(:, :)
      integer :: i

      write (iunit, '(i5)') nt
      do i = 1, nt
         write (iunit, '(3f8.4)') sbarvec(1, i), sbarvec(2, i), sbarvec(3, i)
      end do
   end subroutine write_neighbor_vector_dump

   integer function representative_atom_index(this, ii)
      class(lattice), intent(in) :: this
      integer, intent(in) :: ii

      representative_atom_index = 0

      if (allocated(this%iu)) then
         if (ii <= size(this%iu)) then
            if (this%iu(ii) > 0) representative_atom_index = this%iu(ii)
         end if
      end if
      if (representative_atom_index == 0 .and. allocated(this%ib)) then
         if (ii <= size(this%ib)) then
            if (this%ib(ii) > 0) representative_atom_index = this%ib(ii)
         end if
      end if
      if (representative_atom_index == 0 .and. allocated(this%irec)) then
         if (ii <= size(this%irec)) then
            if (this%irec(ii) > 0) representative_atom_index = this%irec(ii)
         end if
      end if
      if (representative_atom_index == 0) then
         if (ii <= this%kk) representative_atom_index = ii
      end if
      if (representative_atom_index <= 0 .or. representative_atom_index > this%kk) then
         call g_logger%fatal('Failed to resolve representative atom index for strux mapping', __FILE__, __LINE__)
      end if
   end function representative_atom_index

   integer function primitive_basis_label(this, ia)
      class(lattice), intent(in) :: this
      integer, intent(in) :: ia

      primitive_basis_label = 0
      if (allocated(this%no)) then
         if (ia >= 1 .and. ia <= size(this%no)) then
            primitive_basis_label = this%no(ia)
         end if
      end if
      if (primitive_basis_label <= 0 .and. allocated(this%num)) then
         if (ia >= 1 .and. ia <= size(this%num)) then
            primitive_basis_label = this%num(ia)
         end if
      end if
      if (primitive_basis_label <= 0) then
         call g_logger%fatal('Failed to resolve primitive basis label for strux mapping', __FILE__, __LINE__)
      end if
   end function primitive_basis_label

   integer function find_pair_by_vector(nttab, iax, plat, pos, alat, ib, jb, vec_target, tol)
      integer, intent(in) :: nttab, ib, jb
      integer, intent(in) :: iax(:,:)
      real(rp), intent(in) :: plat(3, 3), pos(:, :), alat, vec_target(3), tol
      integer :: i, n1, n2, n3
      real(rp) :: vec(3)
      real(rp) :: dmax, best_d
      integer :: best_i

      find_pair_by_vector = 0
      best_d = huge(1.0_rp)
      best_i = 0
      do i = 1, nttab
         if (iax(1, i) /= ib) cycle
         if (iax(2, i) /= jb) cycle
         n1 = iax(3, i)
         n2 = iax(4, i)
         n3 = iax(5, i)
         vec(:) = alat*(pos(:, jb) - pos(:, ib) + real(n1, rp)*plat(:, 1) + real(n2, rp)*plat(:, 2) + real(n3, rp)*plat(:, 3))
         dmax = maxval(abs(vec - vec_target))
         if (dmax < best_d) then
            best_d = dmax
            best_i = i
         end if
         if (dmax <= tol) then
            find_pair_by_vector = i
            return
         end if
      end do
      if (best_i /= 0) find_pair_by_vector = best_i
   end function find_pair_by_vector

   integer function find_neighbor_atom_by_vector(this, ia, vec_target, cralat, tol)
      class(lattice), intent(in) :: this
      integer, intent(in) :: ia
      real(rp), intent(in) :: vec_target(:), cralat(:, :)
      real(rp), intent(in) :: tol
      integer :: ja
      real(rp) :: dv(3)

      find_neighbor_atom_by_vector = 0
      do ja = 1, this%kk
         if (ja == ia) cycle
         dv(:) = cralat(:, ja) - cralat(:, ia)
         if (maxval(abs(dv - vec_target)) <= tol) then
            find_neighbor_atom_by_vector = ja
            return
         end if
      end do
   end function find_neighbor_atom_by_vector

   subroutine build_orbital_map(norb, orb_map)
      integer, intent(in) :: norb
      integer, intent(out) :: orb_map(16)
      integer :: i
      integer, parameter :: full_map(16) = [ &
         1, 4, 2, 3, 5, 6, 8, 9, 7, 13, 14, 12, 15, 11, 16, 10 ]

      orb_map = 0
      do i = 1, min(norb, size(full_map))
         orb_map(i) = full_map(i)
      end do
   end subroutine build_orbital_map

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Makes a list of atoms with their respective number inside the clust file
   !---------------------------------------------------------------------------
   subroutine atomlist(this)
      class(lattice), intent(inout) :: this
      ! Local variables
      integer :: i, j, itype
      real(rp) :: mom_tmp(3)

#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('lattice.atlist', this%atlist, (/this%ntype/))
      call g_safe_alloc%allocate('lattice.ham_i', this%ham_i, (/this%kk/))
#else
      allocate (this%atlist(this%ntype))
      allocate (this%ham_i(this%kk))
#endif

      do i = 1, this%kk
         this%ham_i(i) = i
      end do

      j = 0
      do i = 1, this%nbulk
         j = j + 1
         this%atlist(j) = this%ib(i)
      end do
      do i = 1, this%nrec
         j = j + 1
         this%atlist(j) = this%irec(i)
      end do
      call this%load_symbolic_atoms_if_needed()

      
      write (805, *) this%kk
      write (805, *)
      do i = 1, this%kk
         itype = 1
         if (allocated(this%iz)) then
            if (i <= size(this%iz)) then
               if (this%iz(i) >= 1 .and. this%iz(i) <= size(this%symbolic_atoms)) itype = this%iz(i)
            end if
         end if
         mom_tmp(:) = [0.0_rp, 0.0_rp, 1.0_rp]
         if (allocated(this%symbolic_atoms(itype)%potential%mom)) then
            if (size(this%symbolic_atoms(itype)%potential%mom) >= 3) then
               mom_tmp(:) = this%symbolic_atoms(itype)%potential%mom(1:3)
            end if
         end if
         write (805, '(A6,6F16.6)') (elem_var(int(this%symbolic_atoms(itype)%element%atomic_number))), &
            this%cr(:, i), mom_tmp(:)
      end do
   end subroutine atomlist

   !**************************************************************************
   !> @brief Identifies atoms within the volume defined by vectors a1, a2, and a3.
   !> 
   !> This subroutine checks if atoms, relative to a central atom, are located 
   !> within the parallelepiped volume formed by translating vectors a1, a2, and a3.
   !> 
   !> @param[in]  cr              A 2D real(rp) array (3 x num_atoms) containing the coordinates of all atoms.
   !> @param[in]  num             An integer array (num_atoms) representing structure identifiers for all atoms.
   !> @param[in]  num_atoms       An integer representing the total number of atoms.
   !> @param[in]  central_atom    An integer representing the index of the central atom.
   !> @param[in]  a1, a2, a3      3-element real(rp) arrays representing the primitive vectors defining the volume.
   !> @param[in]  plane_constant  A real(rp) value representing the constant of the plane equation for the surface.
   !> @param[out] atoms_in_volume An allocatable integer array containing the indices of atoms found within the volume.
   !> @param[out] atom_count      An integer representing the number of atoms found within the volume.
   !>
   !> @note The volume is defined by the parallelepiped formed from vectors a1, a2, and a3.
   !**************************************************************************
   subroutine check_atoms_in_volume(this, cr, num, num_atoms, central_atom, a1, a2, a3, plane_constant, atoms_in_volume, atom_count)
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
       allocate(atoms_in_volume(num_atoms))
       atom_count = 0
       atoms_in_volume = 0
       ! Loop over all atoms to find those inside the primitive cell volume
       do i = 1, num_atoms
          ! Calculate the position relative to the central atom
          relative_pos = cr(:, i) - cr(:, central_atom)
          ! Check if the atom is within the parallelepiped defined by a1, a2, and a3
          call this%check_within_volume(relative_pos, a1, a2, a3, inside)
   
          if (inside) then
             atom_count = atom_count + 1
             atoms_in_volume(atom_count) = i
          end if
       end do
       !write(*,*) 'atoms_in_volume', atoms_in_volume
       ! Resize atoms_in_volume array to the actual number of atoms found
       !if (atom_count > 0) then
       !   deallocate(atoms_in_volume); allocate(atoms_in_volume(atom_count))
       !else
       !   deallocate(atoms_in_volume)
       !end if
   
   end subroutine check_atoms_in_volume

   !**************************************************************************
   !> @brief Checks if an atom is inside the volume defined by a1, a2, and a3.
   !> 
   !> This subroutine checks if an atom's position, relative to a central atom, is 
   !> inside the parallelepiped volume formed by the vectors a1, a2, and a3.
   !>
   !> The atom's position, relative to the central atom, is expressed as a linear 
   !> combination of the vectors a1, a2, and a3, i.e.:
   !>
   !> \f$ \mathbf{r}_{\text{rel}} = u \mathbf{a}_1 + v \mathbf{a}_2 + w \mathbf{a}_3 \f$
   !>
   !> Where \f$ u, v, w \f$ are the coordinates in the basis defined by \f$ \mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3 \f$.
   !>
   !> The atom is considered inside the volume if:
   !> \f$ 0 \leq u \leq 1 \f$, \f$ 0 \leq v \leq 1 \f$, and \f$ 0 \leq w \leq 1 \f$
   !>
   !> @param[in]  relative_pos  A 3-element real(rp) array representing the atom's position relative to the central atom.
   !> @param[in]  a1, a2, a3    3-element real(rp) arrays representing the vectors defining the volume.
   !> @param[out] inside        A logical value indicating whether the atom is inside the parallelepiped.
   !>
   !> @note The coordinates \( u, v, w \) are computed by solving the linear system:
   !>       \f[
   !>       \begin{aligned}
   !>       \mathbf{r}_{\text{rel}} &= u \mathbf{a}_1 + v \mathbf{a}_2 + w \mathbf{a}_3 \\
   !>       \end{aligned}
   !>       \f]
   !>       The inverse of the matrix of dot products is used to calculate \( u \), \( v \), and \( w \).
   !**************************************************************************
   subroutine check_within_volume(this, relative_pos, a1, a2, a3, inside)
       class(lattice), intent(inout) :: this
       real(rp), intent(in) :: relative_pos(3), a1(3), a2(3), a3(3)
       logical, intent(out) :: inside
       real(rp) :: dot11, dot12, dot13, dot22, dot23, dot33
       real(rp) :: dot1r, dot2r, dot3r, inv_denom, u, v, w
   
       ! Calculate dot products between the vectors
       dot11 = dot_product(a1, a1)
       dot12 = dot_product(a1, a2)
       dot13 = dot_product(a1, a3)
       dot22 = dot_product(a2, a2)
       dot23 = dot_product(a2, a3)
       dot33 = dot_product(a3, a3)
       dot1r = dot_product(a1, relative_pos)
       dot2r = dot_product(a2, relative_pos)
       dot3r = dot_product(a3, relative_pos)
   
       ! Calculate inverse of the denominator for the linear combination
       inv_denom = 1.0_rp / (dot11 * (dot22 * dot33 - dot23 * dot23) &
                           - dot12 * (dot12 * dot33 - dot23 * dot13) &
                           + dot13 * (dot12 * dot23 - dot22 * dot13))
   
       ! Calculate coordinates (u, v, w) in the basis defined by a1, a2, and a3
       u = ((dot22 * dot33 - dot23 * dot23) * dot1r + &
           (dot13 * dot23 - dot12 * dot33) * dot2r + &
           (dot12 * dot23 - dot13 * dot22) * dot3r) * inv_denom
   
       v = ((dot13 * dot23 - dot12 * dot33) * dot1r + &
           (dot11 * dot33 - dot13 * dot13) * dot2r + &
           (dot12 * dot13 - dot11 * dot23) * dot3r) * inv_denom
   
       w = ((dot12 * dot23 - dot13 * dot22) * dot1r + &
           (dot12 * dot13 - dot11 * dot23) * dot2r + &
           (dot11 * dot22 - dot12 * dot12) * dot3r) * inv_denom
   
       ! Check if the atom is inside the parallelepiped
       inside = (u >= 0.0_rp .and. u <= 1.0_rp .and. &
                 v >= 0.0_rp .and. v <= 1.0_rp .and. &
                 w >= 0.0_rp .and. w <= 1.0_rp)
   end subroutine check_within_volume

   !**************************************************************************
   !> @brief Finds all unique structure types within a given list of atoms.
   !> 
   !> This subroutine identifies all unique structure types (or identifiers) present 
   !> within a set of atoms and returns a list of these unique identifiers.
   !> 
   !> @param[in]  num          An integer array (num_atoms) containing structure identifiers for all atoms.
   !> @param[in]  num_atoms    An integer representing the total number of atoms.
   !> @param[out] unique_nums  An allocatable integer array containing the unique structure identifiers.
   !> 
   !> @note This subroutine loops through the structure identifiers and ensures only
   !>       unique ones are collected.
   !**************************************************************************
   subroutine find_unique_struct(this, num, num_atoms, unique_nums)
       class(lattice), intent(inout) :: this
       integer, intent(in) :: num(:), num_atoms
       integer, allocatable, intent(out) :: unique_nums(:)
       integer, allocatable :: temp_nums(:)
       integer :: i, j, num_unique
       logical :: found
   
       ! Allocate temporary array for unique numbers
       allocate(temp_nums(num_atoms))
       num_unique = 0
   
       ! Loop over all atoms to find unique structure types
       do i = 1, num_atoms
          found = .false.
          do j = 1, num_unique
             if (num(i) == temp_nums(j)) then
                found = .true.
                exit
             end if
          end do
          if (.not. found) then
             num_unique = num_unique + 1
             temp_nums(num_unique) = num(i)
          end if
       end do
   
       ! Resize array to the actual number of unique nums
       allocate(unique_nums(num_unique))
       unique_nums(:) = temp_nums(1:num_unique)
       ! Deallocate temporary array
       deallocate(temp_nums)
   end subroutine find_unique_struct

   !**************************************************************************
   !> @brief Identifies unique atomic positions within the primitive cell volume.
   !> 
   !> This subroutine identifies which atomic positions within the primitive cell 
   !> volume are unique. An atom is considered unique if it cannot be generated 
   !> from another atom using a combination of translations along a1, a2, and a3.
   !> 
   !> @param[in]  cr                A 2D real(rp) array (3 x num_atoms) containing coordinates of all atoms.
   !> @param[in]  num_atoms         An integer representing the total number of atoms.
   !> @param[in]  atoms_in_volume   An integer array containing indices of atoms within the volume.
   !> @param[in]  atom_count        An integer representing the number of atoms within the volume.
   !> @param[in]  a1, a2, a3        3-element real(rp) arrays representing the vectors defining the volume.
   !> @param[out] unique_atoms      An allocatable integer array containing indices of unique atoms within the volume.
   !> @param[out] unique_atom_count An integer representing the number of unique atoms within the volume.
   !> 
   !> @note The uniqueness of an atom is determined by comparing its position to all
   !>       other atoms and checking if it can be generated by a translation using 
   !>       a1, a2, and a3.
   !**************************************************************************
   subroutine identify_unique_atoms(this, cr, num_atoms, atoms_in_volume, atom_count, a1, a2, a3, unique_atoms, unique_atom_count)
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
       allocate(temp_unique_atoms(atom_count))
       unique_atom_count = 0
       temp_unique_atoms = -1
       ! First pass to identify all unique atoms within the primitive cell volume
       do i = 1, atom_count
           found = .false.
   
           ! Check if this atom can be generated by translating another atom
           do j = 1, unique_atom_count
               ! Translate atom by all combinations of a1, a2, and a3 within the cell
               do k = -1, 1
                   do n = -1, 1
                       do p = -1, 1
                           trans_atom = cr(:, temp_unique_atoms(j)) + k * a1 + n * a2 + p * a3
                           delta = cr(:, atoms_in_volume(i)) - trans_atom
                           if (norm2(delta) < 1.0d-6) then
                               found = .true.
                               exit
                           end if
                       end do
                       if (found) exit
                   end do
                   if (found) exit
               end do
           end do
   
           ! If the atom is not redundant, add it to the temporary list of unique atoms
           if (.not. found) then
               unique_atom_count = unique_atom_count + 1
               temp_unique_atoms(unique_atom_count) = atoms_in_volume(i)
           end if
       end do
   
       ! Allocate the final unique_atoms array to the correct size
       allocate(unique_atoms(unique_atom_count))
       unique_atoms(:) = temp_unique_atoms(1:unique_atom_count)
       ! Deallocate temporary array
       deallocate(temp_unique_atoms)
   end subroutine identify_unique_atoms

   ! Local library
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !---------------------------------------------------------------------------
   subroutine dbar1(this, ia, r2, wav, crd, nat, ndi, np, nr, ii)
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

      nt = this%kk ! Use cluster size instead of fixed 500 for neighbours
      ! allocate (cr(3, nt))
      allocate (sbar(np, np, nt))
      allocate(sbarvec(3, nt))
      call this%clusba(r2, crd, ia, nat, ndi, nt, sbarvec)
      write (17, 10000) nt
      write (17, 10002) ((sbarvec(j, i), j=1, 3), i=1, nt)
      !call write_neighbor_vector_dump(17, this%sbarvec, nt)
      !write (17, 10001) ((this%sbarvec(j, i), j=1, 3), i=1, nt)
      nrl = np*nt
      na = (nrl*(nrl + 1))/2
      allocate (a(na))
      allocate (bet(nrl))
      allocate (wk(nrl))
      allocate (s(nrl, nrl))
      call micha(wav, sbarvec, nt, np, nrl, na, sbar, a, wk, bet, s, ia, r2)
      !call micha(wav, this%sbarvec, nt, np, nrl, na, sbar, a, wk, bet, s, ia, r2)

      ! Saving parameters to be used in the Hamiltonian build
      ! Call clusba again for proper number of TB neighbours
      call this%clusba((r2/9.0d0), crd, ia, nat, ndi, nt, sbarvec)
      ! Store the number of neighbours
      this%nn_max = nt

      do m = 1, nt
         do i = 1, np
            do j = 1, np
               this%sbar(i, j, m, ii) = sbar(i, j, m)
            end do
         end do
         jclus_dbg = ia
         if (m > 1 .and. allocated(this%nn)) then
            if (ia >= 1 .and. ia <= size(this%nn, 1) .and. m <= size(this%nn, 2)) then
               if (this%nn(ia, m) > 0) jclus_dbg = this%nn(ia, m)
            end if
         end if
         write (16, '(" VECTOR=",3f12.6,"   ICLUS=",i5," JCLUS=",i5)') this%sbarvec(1, m), this%sbarvec(2, m), this%sbarvec(3, m), ia, jclus_dbg
         do i = 1, np
            write (16, '(*(f10.4))') (real(this%sbar(i, j, m, ii), rp), j=1, np)
         end do
      end do
      flush(17)

      !do m=1, nt
      !  write(*, *) ia, m
      !  do i=1, 9
      !    write(*, ´(9f10.4)´)((real(this%sbar(i, j, m, ia))), j=1, 9)
      !  end do
      !end do
      deallocate (a, bet, wk, s, sbar, sbarvec)
      return

10000 format(i5)
10001 format(3f8.4)
10002 format(3f8.4)
   end subroutine dbar1

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !
   !> @param[in] r2
   !> @param[in] crd
   !> @param[in] ia
   !> @param[in] nat
   !> @param[in] ndi
   !> @param[inout] n
   !> @param[inout] cr
   !> @return type(calculation)
   !---------------------------------------------------------------------------
   subroutine clusba(this, r2, crd, ia, nat, ndi, n, sbarvec_out)
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

#ifdef USE_SAFE_ALLOC
      if (allocated(this%sbarvec)) call g_safe_alloc%deallocate('lattice.sbarvec', this%sbarvec)
      call g_safe_alloc%allocate('lattice.sbarvec', this%sbarvec, (/3, this%kk/))
#else
      if (allocated(this%sbarvec)) deallocate (this%sbarvec)
      allocate (this%sbarvec(3, this%kk))
#endif

      this%sbarvec(:, :) = 0.0d0
      if (present(sbarvec_out)) sbarvec_out(:, :) = 0.0d0

      ii = 1
      do k = 1, 3
         this%sbarvec(k, 1) = 0.0d0
      end do
      if (present(sbarvec_out)) sbarvec_out(:, 1) = 0.0d0

      do nn = 1, nat
         s1 = 0.0
         do i = 1, 3
            dum(i) = (crd(i, nn) - crd(i, ia))**2
            s1 = s1 + dum(i)
         end do
         if (s1 < r2 .and. s1 > 0.0001) then
            ii = ii + 1
            this%sbarvec(1, ii) = crd(1, nn) - crd(1, ia)
            this%sbarvec(2, ii) = crd(2, nn) - crd(2, ia)
            this%sbarvec(3, ii) = crd(3, nn) - crd(3, ia)
            if (present(sbarvec_out)) then
               if (ii <= size(sbarvec_out, 2)) then
                  sbarvec_out(:, ii) = this%sbarvec(:, ii)
               end if
            end if
         end if
      end do
      n = ii
   end subroutine clusba

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Program to make screened structure constants
   !> S(ME)=-0.5*S(OLE)   J(ME)=2*J(OLE)  QBAR(ME)=2*QBAR(OLE)
   !
   !> @param[in] rws
   !> @param[in] r
   !> @param[in] nr
   !> @param[in] nlm
   !> @param[in] nrl
   !> @param[in] na
   !> @param[inout] sbar
   !> @param[inout] a
   !> @param[inout] wk
   !> @param[inout] bet
   !> @param[inout] s
   !> @param[in] iclus
   !> @param[in] r2
   !---------------------------------------------------------------------------
   subroutine micha(rws, r, nr, nlm, nrl, na, sbar, a, wk, bet, s, iclus, r2)
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

      pi = 4.d0*ATAN(1.d0)
      ! ----------------------------------
      call STREZE(rws, r, nr, s, nrl, nlm)
      ! -------------------------------------------
      fak = 2.d0
      !Original faktors
      q(1) = 0.3485d0*fak
      q(2) = 0.05303d0*fak
      q(3) = 0.010714d0*fak
      q(4) = 0.00337d0*fak   ! f-channel screening parameter
      ! Factors from LMTO47
      !q(1) = 0.33727d0 * fak
      !q(2) = 0.05115d0 * fak
      !q(3) = 0.01397d0 * fak
  !! ! New from LMTO47 (fcc Pt/Es)    data QM/0.3500000d0, 0.0571667d0, 0.0168070d0/
      !q(1) = 0.3500000d0* fak
      !q(2) = 0.0571667d0* fak
      !q(3) = 0.0168070d0* fak
      !     NA=(NRL*(NRL+1))/2
      call SHLDCH(r, nr, nlm, nrl, s, a, na, q, bet, wk, sbar, iclus, r2)
   end subroutine micha

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !>  Makes big matrix of all structure constants connecting sites R
   !
   !> @param[in] w
   !> @param[in] r
   !> @param[in] nr
   !> @param[inout] s
   !> @param[in] nrl
   !> @param[in] nlm
   !---------------------------------------------------------------------------
   subroutine STREZE(w, r, nr, s, nrl, nlm)
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
      do ir = 1, nr
         irl0 = (ir - 1)*nlm
         do jr = 1, nr
            jrl0 = (jr - 1)*nlm
            dr(1) = (r(1, jr) - r(1, ir))/w
            dr(2) = (r(2, jr) - r(2, ir))/w
            dr(3) = (r(3, jr) - r(3, ir))/w
            rr = SQRT(dr(1)**2 + dr(2)**2 + dr(3)**2)
            w1 = 1.d0
            call CANSO(w1, dr, s0)
            do jlm = 1, nlm
               do ilm = 1, nlm
                  s(ilm + irl0, jlm + jrl0) = s0(ilm, jlm)
               end do
            end do
         end do
      end do
   end subroutine streze

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !>  Solves for the screened structure constants
   !
   !> @param[in] r
   !> @param[in] nr
   !> @param[in] nlm
   !> @param[in] nrl
   !> @param[inout] s
   !> @param[inout] a
   !> @param[in] na
   !> @param[in] q
   !> @param[inout] bet
   !> @param[inout] wk
   !> @param[inout] sbar
   !> @param[in] iclus
   !> @param[in] r2
   !---------------------------------------------------------------------------
   subroutine shldch(r, nr, nlm, nrl, s, a, na, q, bet, wk, sbar, iclus, r2)
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
      ndef = 0
      lmax = LL(nlm)
      allocate (s_temp(nrl, nrl))
      write (17, 10000) lmax, q
      irl = 0
      do ir = 1, nr
         do l = 0, lmax
            do m = 1, 2*l + 1
               irl = irl + 1
               bet(irl) = 1.d0/q(l + 1)
            end do
         end do
      end do
      s_temp = s
      do i = 1, nrl
         s_temp(i, i) = s_temp(i, i) + bet(i)
      end do
!    ia = 0
!    do i = 1, nrl
!      do j = 1, i
!        ia = ia + 1
!        a(ia) = s(i, j)
!        if (i == j) then
!          a(ia) = a(ia) + bet(i)
!        end if
!      end do
!    end do
!    call chlr2f(a, na, wk, nrl, ndef)
      call DPOTRF('U', nrl, s_temp, nrl, info)
      write (17, 10001) ndef
!    call chlr2s(a, na, s, nrl, nlm)
      call DPOTRS('U', nrl, nlm, s_temp, nrl, s, nrl, INFO)
      deallocate (s_temp)
      do ilm = 1, nlm
         do irl = 1, nrl
            s(irl, ilm) = -bet(irl)*s(irl, ilm)
         end do
      end do
      ! --------------------------------
      ir = 0
      hitc = 0
      !print *, ´ nr = ´, nr
      do ir = 1, nr
         if (abs(R(1, IR)**2 + R(2, IR)**2 + R(3, IR)**2 - &
                 R(1, 1)**2 + R(2, 1)**2 + R(3, 1)**2) <= (R2/ncut)) then
            ! Legacy view.sbar printing is emitted in dbar1 where neighbour
            ! atom indices are available from nn(ia,m).
	    ! write (16, 10002) r(1, ir), r(2, ir), r(3, ir), iclus
            hitc = hitc + 1
            irl0 = (ir - 1)*nlm
            do ilm = 1, nlm
               do jlm = 1, nlm
                  sbar(ilm, jlm, hitc) = 2.*s(ilm + irl0, jlm)
               end do
            end do
          !!! Changed for ´symmetrizing´ Sbar
          !! NOT NEEDED SINCE WITH LARGER CUTOFF
          !! SBAR IS CALCULATED PROPERLY
            ! do l = 2, 9
            !    l2 = l - 1
            !    do j = 1, l2
            !       sbar(l, j, ir) = sbar(j, l, ir)
            !    end do
            ! end do
            ! do l = 1, 3
            !    sbar(l+1, 1, ir) = -sbar(l+1, 1, ir)
            ! end do
            ! do l = 5, 9
            !    do j = 2, 4
            !       sbar(l, j, ir) = -sbar(l, j, ir)
            !    end do
            ! end do
            !
            do isb = 1, nlm
               write (13) (sbar(isb, jsb, hitc), jsb=1, nlm)
            end do
            ! 208 FORMAT(5F12.6:/(12X, 4F12.6))
            do isb = 1, nlm
               ! Legacy view.sbar printing is emitted in dbar1.
            end do
         end if
      end do
      !print *, ´ hitc = ´, hitc
      return
      !
      ! ... Format Declarations ...
      !
10000 format(" LMAX=", i2, "   Q=", 4f10.6)
10001 format(" NDEF=", i10)
10002 format(" VECTOR=", 3f12.6, "   ICLUS=", i5)
10003 format(9f10.4)
   end subroutine shldch

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculation of the unscreened structure constants in real
   !> space using the slater_koster table.
   !
   !> @param[in] w
   !> @param[in] dr
   !> @param[inout] sc
   !---------------------------------------------------------------------------
   subroutine canso(w, dr, sc)
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
      data ip/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16/
      !     S=1, X=2, Y=3, Z=4, XY=5, YZ=6, ZX=7, X**2-Y**2=8, 3Z*Z-R*R=9
      !     (10-16 reserved for f-orbitals)
      r1 = dr(1)
      r2 = dr(2)
      r3 = dr(3)
      rr = SQRT(r1*r1 + r2*r2 + r3*r3)
      do i = 1, 16
         do j = 1, 16
            sc(i, j) = 0.0d0
         end do
      end do
      if (rr/w <= 0.30d0) return
      sbyr = w/rr
      s2 = sbyr*sbyr
      s3 = s2*sbyr
      s4 = s3*sbyr
      s5 = s4*sbyr
      s6 = s5*sbyr
      s7 = s6*sbyr
      sq3 = SQRT(3.d0)
      sq5 = SQRT(5.d0)
      sq7 = SQRT(7.d0)
      el = r1/rr
      em = r2/rr
      en = r3/rr
      el2 = el*el
      em2 = em*em
      en2 = en*en
      elem = el*em
      elen = el*en
      emen = em*en
      !
      sc(1, 1) = -2.0d0*sbyr
      sc(1, 2) = el*s2*2.d0*sq3
      sc(1, 3) = em*s2*2.d0*sq3
      sc(1, 4) = en*s2*2.d0*sq3
      sc(1, 5) = -2.d0*sq3*sq5*elem*s3
      sc(1, 6) = -2.d0*sq3*sq5*emen*s3
      sc(1, 7) = -2.d0*sq3*sq5*elen*s3
      sc(1, 8) = -sq3*sq5*s3*(el2 - em2)
      sc(1, 9) = sq5*s3*(1.d0 - 3.d0*en2)
      !-----------------------------------------------------------------------
      sc(2, 2) = (3.0d0*el2 - 1.0)*6.d0*s3
      sc(2, 3) = 18.0d0*s3*elem
      sc(2, 4) = 18.0d0*s3*elen
      sc(2, 5) = 6.d0*sq5*s4*em*(1.0d0 - 5.0d0*el2)
      sc(2, 6) = -30.0d0*sq5*s4*elem*en
      sc(2, 7) = 6.d0*sq5*s4*en*(1.0d0 - 5.0d0*el2)
      sc(2, 8) = 6.0d0*sq5*s4*el*(1.0d0 - 2.5d0*el2 + 2.5d0*em2)
      sc(2, 9) = 3.d0*sq3*sq5*s4*el*(1.0d0 - 5.0d0*en2)
      !-----------------------------------------------------------------------
      sc(3, 3) = 6.d0*s3*(3.0d0*em2 - 1.d0)
      sc(3, 4) = 18.d0*s3*emen
      sc(3, 5) = 6.0d0*sq5*s4*el*(1.d0 - 5.d0*em2)
      sc(3, 6) = 6.d0*sq5*s4*en*(1.d0 - 5.d0*em2)
      sc(3, 7) = sc(2, 6)
      sc(3, 8) = -6.d0*sq5*s4*em*(1.d0 - 2.5d0*em2 + 2.5d0*el2)
      sc(3, 9) = 3.d0*sq3*sq5*s4*em*(1.d0 - 5.d0*en2)
      !----------------------------------------------------------------------
      sc(4, 4) = 6.d0*s3*(3.d0*en2 - 1.d0)
      sc(4, 5) = sc(2, 6)
      sc(4, 6) = 6.d0*sq5*s4*em*(1.d0 - 5.d0*en2)
      sc(4, 7) = 6.d0*sq5*s4*el*(1.d0 - 5.d0*en2)
      sc(4, 8) = -15.d0*sq5*s4*en*(el*el - em2)
      sc(4, 9) = 3.d0*sq3*sq5*s4*en*(3.0d0 - 5.0d0*en2)
      !-----------------------------------------------------------------------
      sc(5, 5) = 10.d0*s5*(-35.d0*el2*em2 - 5.d0*en2 + 4.d0)
      sc(5, 6) = -50.d0*s5*elen*(7.d0*em2 - 1.d0)
      sc(5, 7) = -50.d0*s5*emen*(7.d0*el2 - 1.d0)
      sc(5, 8) = -175.d0*s5*elem*(el2 - em2)
      sc(5, 9) = -25.d0*sq3*s5*elem*(7.d0*en2 - 1.d0)
      !----------------------------------------------------------------------
      sc(6, 6) = 10.d0*s5*(-35.d0*em2*en2 - 5.d0*el2 + 4.d0)
      sc(6, 7) = -50.d0*s5*elem*(7.d0*en2 - 1.d0)
      sc(6, 8) = 50.d0*s5*emen*(3.5d0*em2 - 3.5d0*el2 - 1.d0)
      sc(6, 9) = -25.d0*sq3*s5*emen*(7.d0*en2 - 3.d0)
      !-----------------------------------------------------------------------
      sc(7, 7) = 10.d0*s5*(-35.d0*el2*en2 - 5.d0*em2 + 4.d0)
      sc(7, 8) = -50.d0*s5*elen*(3.5d0*el2 - 3.5d0*em2 - 1.0d0)
      sc(7, 9) = -25.d0*sq3*s5*elen*(7.d0*en2 - 3.d0)
      !-----------------------------------------------------------------------
      sc(8, 8) = 10.d0*s5*(-8.75d0*(el2 - em2)**2 - 5.d0*en2 + 4.d0)
      sc(8, 9) = -12.5d0*sq3*s5*(7.d0*en2 - 1.d0)*(el2 - em2)
      !-----------------------------------------------------------------------
      sc(9, 9) = -7.5d0*s5*(35.d0*en2*en2 - 30.d0*en2 + 3.d0)
      !-----------------------------------------------------------------------
      ! F-ORBITAL BLOCK (indices 10-16)
      ! From Methfessel, Mossner & Springborg, J. Phys. C 20, 1069 (1987), Table 1
      ! Orbital order: ilm 10=fz3, 11=fxz2, 12=fyz2, 13=fz(x2-y2), 14=fxyz, 15=fx3, 16=fy3
      !
      ! s-f interactions (row 1, cols 10-16)
      sc(1, 10) = sq7*s3*(5.d0*en2*en2 - 3.d0*en2)
      sc(1, 11) = -2.d0*sq7*sq5*elen*s4*(1.d0 - 5.d0*en2)
      sc(1, 12) = -2.d0*sq7*sq5*emen*s4*(1.d0 - 5.d0*en2)
      sc(1, 13) = -2.d0*sq7*sq5*(el2 - em2)*en*s4
      sc(1, 14) = -4.d0*sq7*sq5*elem*en*s4
      sc(1, 15) = sq7*s4*el*(5.d0*el2 - 3.d0)
      sc(1, 16) = sq7*s4*em*(5.d0*em2 - 3.d0)
      !
      ! p-f interactions (rows 2-4, cols 10-16)
      sc(2, 10) = 6.d0*sq7*sq3*s4*el*(5.d0*en2*en2 - 1.d0)
      sc(2, 11) = 3.d0*sq7*sq5*s5*em*(35.d0*en2*en2 - 30.d0*en2 + 3.d0)
      sc(2, 12) = -30.d0*sq7*sq5*s5*elem*en*(en2 - 1.d0)
      sc(2, 13) = 3.d0*sq7*sq5*s5*el*(7.d0*en2 - 1.d0)*(el2 - em2)
      sc(2, 14) = 6.d0*sq7*sq5*s5*elem*en*(7.d0*en2 - 3.d0)
      sc(2, 15) = 3.d0*sq7*s5*el*(35.d0*el2*em2 - 5.d0*el2 - 20.d0*em2 + 4.d0)
      sc(2, 16) = 3.d0*sq7*s5*el2*em*(35.d0*em2 - 15.d0)
      !
      sc(3, 10) = 6.d0*sq7*sq3*s4*em*(5.d0*en2*en2 - 1.d0)
      sc(3, 11) = -30.d0*sq7*sq5*s5*elem*en*(en2 - 1.d0)
      sc(3, 12) = 3.d0*sq7*sq5*s5*el*(35.d0*en2*en2 - 30.d0*en2 + 3.d0)
      sc(3, 13) = 3.d0*sq7*sq5*s5*em*(7.d0*en2 - 1.d0)*(el2 - em2)
      sc(3, 14) = 6.d0*sq7*sq5*s5*emen*en*(7.d0*en2 - 3.d0)
      sc(3, 15) = 3.d0*sq7*s5*el*em2*(35.d0*el2 - 15.d0)
      sc(3, 16) = 3.d0*sq7*s5*em*(35.d0*em2*el2 - 20.d0*el2 - 5.d0*em2 + 4.d0)
      !
      sc(4, 10) = 3.d0*sq7*sq5*s4*en*(35.d0*en2 - 28.d0)
      sc(4, 11) = 6.d0*sq7*sq5*s5*elen*en*(7.d0*en2 - 3.d0)
      sc(4, 12) = 6.d0*sq7*sq5*s5*emen*en*(7.d0*en2 - 3.d0)
      sc(4, 13) = -12.d0*sq7*sq5*s5*(el2 - em2)*en*(7.d0*en2 - 1.d0)
      sc(4, 14) = -12.d0*sq7*sq5*s5*elem*en*(7.d0*en2 - 1.d0)
      sc(4, 15) = 3.d0*sq7*s5*el*(5.d0*el2 - 1.d0)*(7.d0*en2 - 3.d0)
      sc(4, 16) = 3.d0*sq7*s5*em*(5.d0*em2 - 1.d0)*(7.d0*en2 - 3.d0)
      !
      ! d-f interactions (rows 5-9, cols 10-16)
      sc(5, 10) = -15.d0*sq7*sq5*s5*elem*(7.d0*en2*en2 - 6.d0*en2 + 1.d0)
      sc(5, 11) = -35.d0*sq7*s6*emen*elem*(1.d0 - 7.d0*en2)
      sc(5, 12) = -35.d0*sq7*s6*el2*en*elem*(1.d0 - 7.d0*en2)
      sc(5, 13) = -35.d0*sq7*s6*elem*(el2 - em2)*(7.d0*en2 - 1.d0)
      sc(5, 14) = -70.d0*sq7*s6*elem*(el2*en2 - em2*(1.d0 - en2))
      sc(5, 15) = -15.d0*sq7*s6*elem*(el2 - 3.d0*em2)*(7.d0*en2 - 1.d0)
      sc(5, 16) = -15.d0*sq7*s6*elem*(em2 - 3.d0*el2)*(7.d0*en2 - 1.d0)
      !
      sc(6, 10) = -15.d0*sq7*sq5*s5*emen*(7.d0*en2*en2 - 6.d0*en2 + 1.d0)
      sc(6, 11) = -35.d0*sq7*s6*el2*en*emen*(1.d0 - 7.d0*en2)
      sc(6, 12) = -35.d0*sq7*s6*em2*en*emen*(1.d0 - 7.d0*en2)
      sc(6, 13) = 35.d0*sq7*s6*emen*(el2 - em2)*(7.d0*en2 - 1.d0)
      sc(6, 14) = -70.d0*sq7*s6*emen*el2*(7.d0*en2 - 1.d0)
      sc(6, 15) = 15.d0*sq7*s6*emen*el*(em2 - 7.d0*el2)
      sc(6, 16) = -15.d0*sq7*s6*emen*em*(el2 - 7.d0*em2)
      !
      sc(7, 10) = -15.d0*sq7*sq5*s5*elen*(7.d0*en2*en2 - 6.d0*en2 + 1.d0)
      sc(7, 11) = -35.d0*sq7*s6*el2*en*elen*(1.d0 - 7.d0*en2)
      sc(7, 12) = 35.d0*sq7*s6*em2*elen*(1.d0 - 7.d0*en2)
      sc(7, 13) = 35.d0*sq7*s6*elen*(el2 - em2)*(7.d0*en2 - 1.d0)
      sc(7, 14) = -70.d0*sq7*s6*elen*em2*(7.d0*en2 - 1.d0)
      sc(7, 15) = 15.d0*sq7*s6*elen*el*(em2 - 7.d0*el2)
      sc(7, 16) = -15.d0*sq7*s6*elen*em*(el2 - 7.d0*em2)
      !
      sc(8, 10) = -15.d0*sq7*sq5*s5*(el2 - em2)*(7.d0*en2*en2 - 6.d0*en2 + 1.d0)
      sc(8, 11) = 70.d0*sq7*s6*en*(el2 - em2)*elen
      sc(8, 12) = 70.d0*sq7*s6*en*(el2 - em2)*emen
      sc(8, 13) = -70.d0*sq7*s6*(el2 - em2)*(el2*7.d0*en2 - el2 - 7.d0*em2*en2 + em2)
      sc(8, 14) = 140.d0*sq7*s6*elem*(el2 - em2)*(3.d0 - 7.d0*en2)
      sc(8, 15) = 35.d0*sq7*s6*en*el*(el2 - em2)*(7.d0*el2 - 9.d0)
      sc(8, 16) = -35.d0*sq7*s6*en*em*(el2 - em2)*(7.d0*em2 - 9.d0)
      !
      sc(9, 10) = 5.d0*sq7*sq5*s5*en*(21.d0*en2*en2 - 14.d0*en2 + 1.d0)
      sc(9, 11) = 35.d0*sq7*s6*elen*(21.d0*en2 - 5.d0)
      sc(9, 12) = 35.d0*sq7*s6*emen*(21.d0*en2 - 5.d0)
      sc(9, 13) = -70.d0*sq7*s6*en*(el2 - em2)*(21.d0*en2 - 5.d0)
      sc(9, 14) = -140.d0*sq7*s6*elem*en*(21.d0*en2 - 5.d0)
      sc(9, 15) = 35.d0*sq7*s6*el*(21.d0*en2 - 3.d0)*(7.d0*en2 - 1.d0)
      sc(9, 16) = 35.d0*sq7*s6*em*(21.d0*en2 - 3.d0)*(7.d0*en2 - 1.d0)
      !
      ! f-f interactions (rows 10-16, cols 10-16)
      ! Diagonal terms
      sc(10, 10) = -7.d0*s7*(99.d0*en2*en2*en2 - 135.d0*en2*en2 + 55.d0*en2 - 5.d0)
      sc(11, 11) = s7*(-385.d0*en2*en2*(el2 + em2) + 70.d0*en2*(el2 + em2) + 245.d0*el2*em2 + 5.d0)
      sc(12, 12) = sc(11, 11)  ! same for y
      sc(13, 13) = s7*(245.d0*en2*(el2 - em2)**2 - 70.d0*(el2 - em2)**2 - 385.d0*en2*en2*(el2 - em2)**2 + 5.d0)
      sc(14, 14) = s7*(-1225.d0*en2*en2*el2*em2 + 350.d0*en2*el2*em2 - 35.d0*el2*em2 + 5.d0)
      sc(15, 15) = s7*(-231.d0*el2*el2*em2 - 385.d0*el2*el2*en2*en2 + 70.d0*el2*el2*en2 - 55.d0*el2*el2 + &
                       245.d0*el2*em2*em2 + 5.d0)
      sc(16, 16) = s7*(-231.d0*em2*em2*el2 - 385.d0*em2*em2*en2*en2 + 70.d0*em2*em2*en2 - 55.d0*em2*em2 + &
                       245.d0*em2*el2*el2 + 5.d0)
      ! Off-diagonal f-f terms (selected important ones)
      sc(11, 12) = -35.d0*s7*elem*(35.d0*en2*en2 - 10.d0*en2 + 1.d0)
      sc(11, 13) = s7*el*(-245.d0*en2*en2*el2 + 245.d0*en2*en2*em2 + 70.d0*en2*el2 - 70.d0*en2*em2 - 5.d0*el2 + 5.d0*em2)
      sc(11, 14) = -70.d0*s7*elen*em*(35.d0*en2*en2 - 10.d0*en2 + 1.d0)
      sc(11, 15) = s7*el*em*(245.d0*en2*(el2 - em2) - 49.d0*(el2 - em2))
      sc(11, 16) = 70.d0*s7*em2*el*(35.d0*en2 - 5.d0)*(en2 - 1.d0)
      sc(12, 13) = s7*em*(-245.d0*en2*en2*el2 + 245.d0*en2*en2*em2 + 70.d0*en2*el2 - 70.d0*en2*em2 - 5.d0*el2 + 5.d0*em2)
      sc(12, 14) = -70.d0*s7*emen*el*(35.d0*en2*en2 - 10.d0*en2 + 1.d0)
      sc(12, 15) = 70.d0*s7*el2*em*(35.d0*en2 - 5.d0)*(en2 - 1.d0)
      sc(12, 16) = s7*el*em*(245.d0*en2*(el2 - em2) - 49.d0*(el2 - em2))
      sc(13, 14) = -140.d0*s7*elem*(el2 - em2)*(7.d0*en2*en2 - 5.d0*en2 + 1.d0)
      sc(13, 15) = -35.d0*s7*el*(el2 - em2)*(el2 - 7.d0*em2)*(7.d0*en2 - 1.d0)
      sc(13, 16) = -35.d0*s7*em*(el2 - em2)*(em2 - 7.d0*el2)*(7.d0*en2 - 1.d0)
      sc(14, 15) = -35.d0*s7*elem*em*(el2 - 3.d0*em2)*(7.d0*en2 - 1.d0)
      sc(14, 16) = -35.d0*s7*elem*el*(em2 - 3.d0*el2)*(7.d0*en2 - 1.d0)
      sc(15, 16) = 35.d0*s7*el*em*(el2 - em2)*(el2 + em2 - 5.d0*en2*(el2 + em2))
      !-----------------------------------------------------------------------
      do l = 2, 16
         ll = l - 1
         do j = 1, ll
            sc(l, j) = sc(j, l)
         end do
      end do
      do l = 1, 3
         sc(l + 1, 1) = -sc(l + 1, 1)
      end do
      do l = 5, 9
         do j = 2, 4
            sc(l, j) = -sc(l, j)
         end do
      end do
      do l = 10, 16
         do j = 2, 4
            sc(l, j) = -sc(l, j)
         end do
      end do
      ! ------ THIS PART CHANGES THE YLM ORDER AND MULTIPLIES BY -0.5 ----
      do i = 1, 16
         do j = 1, 16
            s(ip(j), ip(i)) = -0.5d0*sc(j, i)
         end do
      end do
      do i = 1, 16
         do j = 1, 16
            sc(j, i) = s(j, i)
         end do
      end do
   end subroutine canso

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Cholesky factorization of (N X N) REAL*8 matrix c.
   !> returns as soon as matrix becomes indefinite, rest is unchanged
   !>
   !
   !> @param[inout] C
   !> @param[in] NA
   !> @param[inout] W Workspace of dimension N
   !> @param[in] N
   !> @param[out] NDEF
   !---------------------------------------------------------------------------
   subroutine chlr2f(c, na, w, n, ndef)
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

      write (17, *) 'chol', na, n
      ic0 = 0
      do i = 1, n
         jc0 = 0
         do j = 1, i - 1
            csum = c(ic0 + j)
            do k = 1, j - 1
               csum = csum - w(k)*c(jc0 + k)
            end do
            jc0 = jc0 + j
            w(j) = csum/c(jc0)
         end do
         csum = c(ic0 + i)
         do k = 1, i - 1
            csum = csum - w(k)*w(k)
         end do
         if (csum <= 0.d0) then
            ndef = i - 1
            goto 1000
         else
            do j = 1, i - 1
               c(ic0 + j) = w(j)
            end do
            !
            ic0 = ic0 + i
            c(ic0) = SQRT(csum)
         end if
      end do
      ndef = n
1000  if (ndef < n) then
         write (17, 10000) n, ndef
      end if
      return
10000 format(" CHLR2F:    N=", i4, "    NDEF=", i4)
   end subroutine chlr2f

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Back substitution after cholesky factorization
   !> solution overwrites RHS in V
   !
   !> @param[in] C
   !> @param[in] NA
   !> @param[inout] V
   !> @param[in] N
   !> @param[in] M Number of RHS´s
   !---------------------------------------------------------------------------
   subroutine chlr2s(c, na, v, n, m)
      implicit none
      ! Input
      integer, intent(in) :: m, n, na
      real(rp), dimension(na), intent(in) :: c
      ! Output
      real(rp), dimension(n, m), intent(inout) :: v
      ! Local Variables
      integer :: i, ic0, ind1, ind2, k, mm
      real(rp) :: csum

      do mm = 1, m
         ic0 = 0
         do i = 1, n
            csum = v(i, mm)
            do k = 1, i - 1
               csum = csum - v(k, mm)*c(ic0 + k)
            end do
            ic0 = ic0 + i
            v(i, mm) = csum/c(ic0)
         end do
         !
         do i = n, 1, -1
            ind2 = (i*(i + 1))/2
            csum = v(i, mm)
            ind1 = ind2 + i
            do k = i + 1, n
               csum = csum - v(k, mm)*c(ind1)
               ind1 = ind1 + k
            end do
            v(i, mm) = csum/c(ind2)
         end do
      end do
   end subroutine chlr2s

   function LL(ilm)
      implicit none
      ! Input
      integer, intent(in) :: ilm
      ! Function Declaration
      integer :: LL
      ! Local variables
      integer, dimension(100) :: lla
      ! Data Declarations
      data lla/0, 3*1, 5*2, 7*3, 9*4, 11*5, 13*6, 15*7, 17*8, 19*9/

      LL = lla(ilm)
   end function LL

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reorders map according to eq vectors
   !
   !> @param[in] crd
   !> @param[in] no
   !> @param[in] iu
   !> @param[inout] nn
   !> @param[in] nat
   !> @param[in] ntot
   !> @param[in] nomx
   !> @param[in] ndi
   !> @param[in] nnmx
   !> @param[inout] set
   !> @param[inout] idnn
   !> @param[out] ret
   !---------------------------------------------------------------------------
   subroutine remd(this, crd, no, iu, nn, nat, ntot, nomx, ndi, nnmx, set, idnn, ret)
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
      do i = 1, ntot
         la = iu(i)
         jsz = nn(la, 1)
         do j = 2, jsz
            jj = nn(la, j)
            if (this%pbc) then 
               call this%f_wrap_coord_diff(nat,crd,la,jj,set(:,i,j))
            else
               do m = 1, 3
                  set(m, i, j) = crd(m, la) - crd(m, jj)
               end do
            end if
         end do
      end do
      do i = 1, nat
         n = no(i)
         do lk = 1, ntot
            ino = iu(lk)
            lm = no(ino)
            if (n == lm) goto 1000
         end do
         write (17, 10003)
1000     imax = nn(ino, 1)
         do iii = 1, imax
            idnn(iii) = 0
         end do
         jsz = nn(i, 1)
         do j = 2, jsz
            jj = nn(i, j)
            if (this%pbc) then
               call this%f_wrap_coord_diff(nat,crd,i,jj,ret)
            else
               do m = 1, 3
                  ret(m) = crd(m, i) - crd(m, jj)
               end do
            end if
            !----------FINDS EQUIVALENT VECTOR------------------------
            eps = .0001
            do ii = 2, imax
               a1 = ret(1) - set(1, n, ii)
               a2 = ret(2) - set(2, n, ii)
               a3 = ret(3) - set(3, n, ii)
               aaa = a1**2 + a2**2 + a3**2
               if (aaa < eps) goto 1100
            end do
            goto 1200
1100        k = ii
            idnn(k) = jj
         end do
         nn(i, 1) = imax
         do j = 2, imax
            nn(i, j) = idnn(j)
         end do
      end do
      do inn = 1, nat
         write (12) (nn(inn, jnn), jnn=1, nn(inn, 1))
      end do
      return
1200  write (17, 10002)
      print *, 'atom:', i
      stop
      !  9  WRITE(6, 223)I, SET(1, I, J), SET(2, I, J), SET(3, I, J)
10000 format(i5, 3f9.4)
      !--REORDERS NEIGHBORS FOR I LARGER THAN NTOT ACCORDING TO TYPICAL NO--
10001 format(3i5)
10002 format(" VECTOR   NOT FOUND ")
10003 format(" TYPE NO  NOT FOUND ")
   end subroutine remd

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> TODO
   !
   !> @param[in] IM
   !> @param[in] IZP
   !> @param[in] NN
   !> @param[in] NO
   !> @param[in] ND
   !> @param[in] NM
   !> @param[in] NTOT
   !---------------------------------------------------------------------------
   subroutine outmap(IM, IZP, NN, NO, ND, NM, NTOT)
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

      write (IM, 10000)
      do I = 1, NTOT
!!    K = NO(I)
         K = IZP(I)
         ID = NN(I, 1)
!!    IDM = MIN(ID, 21)
         IDM = MIN(ID, 16)
         write (IM, 10001) I, K, ID, (NN(I, J), J=2, IDM)
         if (IDM /= NN(I, 1)) then
            !     ID=NN(I, 1)
            write (IM, 10002) (NN(I, J), J=17, ID)
         end if
      end do
      return

10000 format( &
         " NEAREST NEIGHBOUR MAP"/, "      ATOM   TYPE  CONNECTIVITY", 5x, &
         "NEIGHBOURS")
10001 format(1x, 3i4, 2x, 20i5)
10002 format(29x, 16i5)
   end subroutine outmap

   integer function mapa(I, J, R2, DD, CT)
      implicit none
      ! Input
      integer, intent(in) :: I, J
      real(rp), intent(in) :: R2
      real(rp) :: DD
      real(rp), dimension(50), intent(in) :: CT
      ! Local variables
      real(rp) :: CTM, CTSM

      CTM = (CT(1) + CT(1))/2.
      CTSM = CTM**2
      if (R2 >= CTSM) then
         MAPA = 0
      else
         MAPA = 1
      end if
   end function mapa

   subroutine f_wrap_coord_diff(this,Natom,coord,i_atom,j_atom,cdiff)
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
      odiff=coord(:,j_atom) - coord(:,i_atom)
      !
      xmax=0;xmin=0;ymax=0;ymin=0;zmax=0;zmin=0
      if(this%b1)then
         xmax=1
         xmin=-1
      end if
      if(this%b2)then
         ymax=1
         ymin=-1
      end if
      if(this%b3)then
         zmax=1
         zmin=-1
      end if
      
      mdiff=odiff
      do z=zmin,zmax
         do y=ymin,ymax
            do x=xmin,xmax
               oshift = odiff + x*(this%n1)*this%a(:, 1)*this%alat& 
                              + y*(this%n2)*this%a(:, 2)*this%alat&
                              + z*(this%n3)*this%a(:, 3)*this%alat
               if(norm2(oshift)<norm2(mdiff))  mdiff = oshift
            end do
         end do
      end do
      cdiff=mdiff
      return
      !
   end subroutine f_wrap_coord_diff
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Recursion library
   !
   !> @param[inout] ct
   !> @param[in] crd
   !> @param[in] ndim
   !> @param[in] nat
   !> @param[in] izp
   !> @param[inout] nn
   !> @param[in] nd
   !> @param[inout] nm
   !> @param[in] ngbr
   !> @param[in] ntot
   !---------------------------------------------------------------------------
   subroutine nncal(this,ct, crd, ndim, nat, izp, nn, nd, nm, ngbr, ntot)
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

      CAP = NM
      NN(:, :) = 0
      NN(:, 1) = 1

      NNMAX = 0
      IADD = 1
      if (IZP(1) < 0) then
         NN(1, 1) = 1
      end if
      do II = 1, NAT
         if (IZP(II) < 0) goto 1000
      end do
      NN(1, 1) = 1
      do I = 1, NAT
         do J = 2, NM
            NN(I, J) = 0
         end do
      end do
      IADD = 0
      II = 2
1000  if (II <= 1) then
         II = 2
      end if

      RCUT = CT(1)
      RCUT2 = RCUT*RCUT
      if (RCUT2 <= 0.0_rp) then
         NM = 0
         return
      end if

      if (this%pbc .and. this%b1 .and. this%b2 .and. this%b3) then
         CELL(:, 1) = real(this%n1, rp)*this%a(:, 1)*this%alat
         CELL(:, 2) = real(this%n2, rp)*this%a(:, 2)*this%alat
         CELL(:, 3) = real(this%n3, rp)*this%a(:, 3)*this%alat
         CELL_INV = inverse_3x3(CELL)
         DETC = abs(determinant(CELL))

         CROSS_TMP(:, 1) = [ &
            CELL(2, 2)*CELL(3, 3) - CELL(3, 2)*CELL(2, 3), &
            CELL(3, 2)*CELL(1, 3) - CELL(1, 2)*CELL(3, 3), &
            CELL(1, 2)*CELL(2, 3) - CELL(2, 2)*CELL(1, 3) ]
         CROSS_TMP(:, 2) = [ &
            CELL(2, 3)*CELL(3, 1) - CELL(3, 3)*CELL(2, 1), &
            CELL(3, 3)*CELL(1, 1) - CELL(1, 3)*CELL(3, 1), &
            CELL(1, 3)*CELL(2, 1) - CELL(2, 3)*CELL(1, 1) ]
         DDUM = [ &
            CELL(2, 1)*CELL(3, 2) - CELL(3, 1)*CELL(2, 2), &
            CELL(3, 1)*CELL(1, 2) - CELL(1, 1)*CELL(3, 2), &
            CELL(1, 1)*CELL(2, 2) - CELL(2, 1)*CELL(1, 2) ]

         SPAN(1) = DETC/max(sqrt(sum(CROSS_TMP(:, 1)**2)), tiny(1.0_rp))
         SPAN(2) = DETC/max(sqrt(sum(CROSS_TMP(:, 2)**2)), tiny(1.0_rp))
         SPAN(3) = DETC/max(sqrt(sum(DDUM(:)**2)), tiny(1.0_rp))
         NX = max(1, int(SPAN(1)/RCUT))
         NY = max(1, int(SPAN(2)/RCUT))
         NZ = max(1, int(SPAN(3)/RCUT))
         NBIN = NX*NY*NZ

         allocate(FRAC(3, NAT), HEAD(NBIN), NEXT_ATOM(NAT), BIN_X(NAT), BIN_Y(NAT), BIN_Z(NAT), CANDIDATES(NAT))
         HEAD = 0
         NEXT_ATOM = 0
         do I = 1, NAT
            FRAC(:, I) = matmul(CELL_INV, CRD(:, I))
            do L = 1, 3
               FRAC(L, I) = FRAC(L, I) - floor(FRAC(L, I))
            end do
            BIN_X(I) = 1 + min(NX - 1, int(FRAC(1, I)*NX))
            BIN_Y(I) = 1 + min(NY - 1, int(FRAC(2, I)*NY))
            BIN_Z(I) = 1 + min(NZ - 1, int(FRAC(3, I)*NZ))
            BIN_ID = ((BIN_Z(I) - 1)*NY + (BIN_Y(I) - 1))*NX + BIN_X(I)
            NEXT_ATOM(I) = HEAD(BIN_ID)
            HEAD(BIN_ID) = I
         end do

         do I = II, NAT
            if (IZP(I)*IADD > 0) cycle
            CAND_COUNT = 0
            BX = BIN_X(I); BY = BIN_Y(I); BZ = BIN_Z(I)
            RX = merge(1, 0, NX > 1)
            RY = merge(1, 0, NY > 1)
            RZ = merge(1, 0, NZ > 1)
            do DBZ = -RZ, RZ
               IZ = modulo(BZ - 1 + DBZ, NZ) + 1
               do DBY = -RY, RY
                  IY = modulo(BY - 1 + DBY, NY) + 1
                  do DBX = -RX, RX
                     IX = modulo(BX - 1 + DBX, NX) + 1
                     BIN_ID = ((IZ - 1)*NY + (IY - 1))*NX + IX
                     J = HEAD(BIN_ID)
                     do while (J > 0)
                        if (J < I) then
                           DS(:) = FRAC(:, J) - FRAC(:, I)
                           DS(:) = DS(:) - real(nint(DS(:)), rp)
                           DDUM(:) = matmul(CELL, DS)
                           R2 = dot_product(DDUM, DDUM)
                           if (R2 < RCUT2) then
                              CAND_COUNT = CAND_COUNT + 1
                              CANDIDATES(CAND_COUNT) = J
                           end if
                        end if
                        J = NEXT_ATOM(J)
                     end do
                  end do
               end do
            end do

            call sort_integer_list(CANDIDATES, CAND_COUNT)
            do K = 1, CAND_COUNT
               J = CANDIDATES(K)
               ID = NN(I, 1) + 1
               NN(I, 1) = ID
               NN(I, ID) = J
               NNMAX = MAX(NNMAX, ID)
               ID = NN(J, 1) + 1
               NN(J, 1) = ID
               NN(J, ID) = I
               NNMAX = MAX(NNMAX, ID)
               if (NNMAX > CAP) then
                  write (6, '(" TOO MANY NEIGHBOURS")')
                  write (6, '(" NEIGHBOUR MAP AS FAR AS", i6, "TH SITE")') I
                  write (6, *) NNMAX, ID, CAP
                  stop
               end if
            end do
         end do

         deallocate(FRAC, HEAD, NEXT_ATOM, BIN_X, BIN_Y, BIN_Z, CANDIDATES)
      else if (.not. this%pbc) then
         do L = 1, 3
            MINC(L) = minval(CRD(L, 1:NAT))
            MAXC(L) = maxval(CRD(L, 1:NAT))
            SPAN(L) = max(MAXC(L) - MINC(L), RCUT)
         end do
         NX = max(1, int(SPAN(1)/RCUT))
         NY = max(1, int(SPAN(2)/RCUT))
         NZ = max(1, int(SPAN(3)/RCUT))
         BINW(1) = SPAN(1)/real(NX, rp)
         BINW(2) = SPAN(2)/real(NY, rp)
         BINW(3) = SPAN(3)/real(NZ, rp)
         NBIN = NX*NY*NZ

         allocate(HEAD(NBIN), NEXT_ATOM(NAT), BIN_X(NAT), BIN_Y(NAT), BIN_Z(NAT), CANDIDATES(NAT))
         HEAD = 0
         NEXT_ATOM = 0
         do I = 1, NAT
            BIN_X(I) = 1 + min(NX - 1, int((CRD(1, I) - MINC(1))/BINW(1)))
            BIN_Y(I) = 1 + min(NY - 1, int((CRD(2, I) - MINC(2))/BINW(2)))
            BIN_Z(I) = 1 + min(NZ - 1, int((CRD(3, I) - MINC(3))/BINW(3)))
            BIN_ID = ((BIN_Z(I) - 1)*NY + (BIN_Y(I) - 1))*NX + BIN_X(I)
            NEXT_ATOM(I) = HEAD(BIN_ID)
            HEAD(BIN_ID) = I
         end do

         do I = II, NAT
            if (IZP(I)*IADD > 0) cycle
            CAND_COUNT = 0
            BX = BIN_X(I); BY = BIN_Y(I); BZ = BIN_Z(I)
            do DBZ = max(-1, 1 - BZ), min(1, NZ - BZ)
               IZ = BZ + DBZ
               do DBY = max(-1, 1 - BY), min(1, NY - BY)
                  IY = BY + DBY
                  do DBX = max(-1, 1 - BX), min(1, NX - BX)
                     IX = BX + DBX
                     BIN_ID = ((IZ - 1)*NY + (IY - 1))*NX + IX
                     J = HEAD(BIN_ID)
                     do while (J > 0)
                        if (J < I) then
                           DDUM(:) = CRD(:, J) - CRD(:, I)
                           R2 = dot_product(DDUM, DDUM)
                           if (R2 < RCUT2) then
                              CAND_COUNT = CAND_COUNT + 1
                              CANDIDATES(CAND_COUNT) = J
                           end if
                        end if
                        J = NEXT_ATOM(J)
                     end do
                  end do
               end do
            end do

            call sort_integer_list(CANDIDATES, CAND_COUNT)
            do K = 1, CAND_COUNT
               J = CANDIDATES(K)
               ID = NN(I, 1) + 1
               NN(I, 1) = ID
               NN(I, ID) = J
               NNMAX = MAX(NNMAX, ID)
               ID = NN(J, 1) + 1
               NN(J, 1) = ID
               NN(J, ID) = I
               NNMAX = MAX(NNMAX, ID)
               if (NNMAX > CAP) then
                  write (6, '(" TOO MANY NEIGHBOURS")')
                  write (6, '(" NEIGHBOUR MAP AS FAR AS", i6, "TH SITE")') I
                  write (6, *) NNMAX, ID, CAP
                  stop
               end if
            end do
         end do

         deallocate(HEAD, NEXT_ATOM, BIN_X, BIN_Y, BIN_Z, CANDIDATES)
      else
         do I = II, NAT
            if (IZP(I)*IADD <= 0) then
               NN(I, 1) = 1
               IIP = IZP(I)
               IIP = ABS(IIP)
               ILJ = I - 1
               do J = 1, ILJ
                  JJP = IZP(J)
                  JJP = ABS(JJP)
                  R2 = 0.0
                  if (this%pbc) then
                     call this%f_wrap_coord_diff(nat, crd, i, j, ddum)
                     r2 = sum(ddum(:)**2)
                  else        
                     do L = 1, 3
                        DDUM(L) = CRD(L, I) - CRD(L, J)
                        R2 = R2 + DDUM(L)*DDUM(L)
                     end do
                  end if
                  ID = NGBR(IIP, JJP, R2, DUM, CT)
                  if (ID /= 0) then
                     ID = NN(I, 1) + 1
                     NN(I, 1) = ID
                     NN(I, ID) = J
                     NNMAX = MAX(NNMAX, ID)
                     ID = NN(J, 1) + 1
                     NN(J, 1) = ID
                     NN(J, ID) = I
                     NNMAX = MAX(NNMAX, ID)
                     if (NNMAX > NM) then
                        write (6, '(" TOO MANY NEIGHBOURS")')
                        write (6, '(" NEIGHBOUR MAP AS FAR AS", i6, "TH SITE")') I
                        write (6, *) NNMAX, ID, NM
                        stop
                     end if
                  end if
               end do
            end if
         end do
      end if
      nm = nnmax
      return

contains

      subroutine sort_integer_list(list, n)
         integer, intent(inout) :: list(:)
         integer, intent(in) :: n
         integer :: a, b, value

         do a = 2, n
            value = list(a)
            b = a - 1
            do while (b >= 1)
               if (list(b) <= value) exit
               list(b + 1) = list(b)
               b = b - 1
            end do
            list(b + 1) = value
         end do
      end subroutine sort_integer_list
   end subroutine
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read CR of clust file
   !
   !> @param[in] alat
   !> @param[in] nndim
   !> @param[inout] cr
   !> @param[inout] iz
   !> @param[inout] n
   !> @param[in] ip Unit of opened file
   !---------------------------------------------------------------------------
   subroutine leia(alat, nndim, cr, iz, n, ip)
      implicit none
      ! Input
      integer, intent(in) :: ip, nndim
      real(rp), intent(in) :: alat
      ! Output
      integer, dimension(nndim + 10), intent(inout) :: iz, n
      real(rp), dimension(3, nndim + 10), intent(inout) :: cr
      ! Local variables
      integer :: i, j

      read (IP, *)
      do I = 1, nndim, 2
         read (IP, *) &
            CR(1, I), CR(2, I), CR(3, I), IZ(I), N(I), CR(1, I + 1), CR(2, I + 1), &
            CR(3, I + 1), IZ(I + 1), N(I + 1)
      end do
      do I = 1, nndim
         do J = 1, 3
            CR(J, I) = ALAT*CR(J, I)
         end do
      end do
      return
   end subroutine
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Subrotina que ordena um determinado vetor ´M´ em ordem decrescente
   !>     segundo o metodo da bolha (BUBBLE)
   !> ...M(1) >= M(2)...  ; onde:
   !>     NL, ndim   -  input
   !>     M         -  input/output (na saida M contem os elementos dispostos
   !>     em ordem decrescente em modulo.
   !>
   !> @param[in] nl
   !> @param[in] ndim
   !> @param[inout] m
   !> @param[in] nd
   !> @param[in] nt
   !---------------------------------------------------------------------------
   subroutine bubble(nl, ndim, m, nd, nt)
      implicit none
      ! Inputs
      integer, intent(in) :: ndim, nl, nd, nt
      ! Output
      real(rp), dimension(ndim, nd), intent(inout) :: m
      ! Local variables
      integer :: ind, inic, j, k
      real(rp) :: fim, z

      IND = 1
      INIC = 2
      FIM = NL
      do while ((IND == 1) .and. (INIC <= FIM))
         IND = 0
         do J = INT(FIM), INIC, -1
            if (M(J, NT) < M(J - 1, NT)) then
            do K = 1, ND
               Z = M(J, K)
               M(J, K) = M(J - 1, K)
               M(J - 1, K) = Z
            end do
            end if
         end do
         FIM = FIM - 1
         do J = INIC, INT(FIM)
            if (M(J + 1, NT) < M(J, NT)) then
            do K = 1, ND
               Z = M(J + 1, K)
               M(J + 1, K) = M(J, K)
               M(J, K) = Z
            end do
            IND = 1
            end if
         end do
         INIC = INIC + 1
      end do
   end subroutine bubble

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> CUT:subprogram of the builds package
   !> Cuts a sphere of square raius RS around every atom in the unit
   !> cell. Best for cells with large numbers of atoms...  hould
   !> Carefull with monoatomic systems!!! In this case should cut
   !> around several sites, since spherical clusters should be avoided.
   !>
   !> @param[in] i
   !> @param[in] l
   !> @param[in] ndim
   !> @param[in] crd
   !> @param[inout] cr
   !> @param[in] izp
   !> @param[inout] iz
   !> @param[inout] num
   !> @param[in] no
   !> @param[in] rs
   !> @param[in] ii
   !---------------------------------------------------------------------------
   subroutine cut(i, l, ndim, crd, cr, izp, iz, num, no, rs, ii)
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

      do na = 1, l
         r2 = 0.0d0
         do j = 1, 3
            dum(j) = (crd(j, na) - crd(j, i))**2
            r2 = r2 + dum(j)
         end do
         if (r2 .le. rs) then
            ii = ii + 1
            cr(1, ii) = crd(1, i)
            cr(2, ii) = crd(2, i)
            cr(3, ii) = crd(3, i)
            iz(ii) = izp(i)
            num(ii) = no(i)
            return
         end if
      end do
      return
   end subroutine cut

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Check if object is well fulfilled
   !---------------------------------------------------------------------------
   subroutine check_all(this)
      implicit none
      class(lattice) :: this

      if (this%crystal_sym /= 'bcc' &
          .and. this%crystal_sym /= 'b2' &
          .and. this%crystal_sym /= 'fcc' &
          .and. this%crystal_sym /= 'hcp' &
          .and. this%crystal_sym /= 'fcc2' &
          .and. this%crystal_sym /= 'fcc3' &
          .and. this%crystal_sym /= 'file') then
         call g_logger%fatal('lattice%crystal_sym must be one of: ''bcc'', ''fcc'', ''hcp'' or  ''file''")', __FILE__, __LINE__)
      end if
   end subroutine check_all

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculate reduced_acr and reduced_nbas
   !---------------------------------------------------------------------------
   subroutine calculate_nbas(this)
      implicit none
      class(lattice) :: this
      integer :: size_iz
      integer, dimension(:) :: atype(this%nbas), amount(this%nbas)
      integer :: i, j

      amount(:) = 0
      atype(:) = 0

      j = 1

      amount(j) = 1
      atype(j) = this%iz(j)

      do i = 2, this%nbas
         if (this%iz(i) .eq. this%iz(i - 1)) then
            amount(j) = amount(j) + 1
         else
            j = j + 1
            atype(j) = this%iz(i)
            amount(j) = 1
         end if
      end do

      this%reduced_nbas = 0

      do i = 1, size(amount)
         if (amount(i) .eq. 0) exit
         this%reduced_nbas = this%reduced_nbas + 1
      end do

#ifdef USE_SAFE_ALLOC
      if (allocated(this%reduced_acr)) call g_safe_alloc%deallocate('lattice.reduced_acr', this%reduced_acr)
      call g_safe_alloc%allocate('lattice.reduced_acr', this%reduced_acr, (/this%reduced_nbas/))
#else
      if (allocated(this%reduced_acr)) deallocate (this%reduced_acr)
      allocate (this%reduced_acr(this%reduced_nbas))
#endif

      do i = 1, size(this%reduced_acr)
         this%reduced_acr(i) = atype(i)
      end do

   end subroutine calculate_nbas

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
   !> @param[in] unit File unit used to write namelist
   !> @param[in] file File name used to write namelist
   !---------------------------------------------------------------------------
   subroutine print_state_full(this, unit, file)
      class(lattice), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

      include 'include_codes/namelists/lattice.f90'

      ! scalar

      zmin = this%zmin
      zmax = this%zmax
      zstep = this%zstep
      wav = this%wav
      vol = this%vol
      rc = this%rc
      r2 = this%r2
      celldm = this%celldm
      alat = this%alat
      reduced_nbas = this%reduced_nbas
      ntype = this%ntype
      ntot = this%ntot
      nrec = this%nrec
      nmax = this%nmax
      nlay = this%nlay
      ndim = this%ndim
      npe = this%npe
      nclu = this%nclu
      nbulk_bulk = this%nbulk_bulk
      nbulk = this%nbulk
      nbas = this%nbas
      kk = this%kk
      dx = this%dx
      dy = this%dy
      dz = this%dz
      dw = this%dw
      crystal_sym = this%crystal_sym
      surftype = this%surftype
      strux_backend = this%strux_backend
      screening = this%screening
      strux_want_sdot = this%strux_want_sdot
      strux_solve_scale = this%strux_solve_scale
      a = this%a

      ! one dimensional allocatables

      if (allocated(this%z)) then
         allocate (z, mold=this%z)
         z = this%z
      else
         allocate (z(0))
      end if
      if (allocated(this%ct)) then
         allocate (ct, mold=this%ct)
         ct = this%ct
      else
         allocate (ct(0))
      end if
      if (allocated(this%screening_alpha)) then
         allocate (screening_alpha, mold=this%screening_alpha)
         screening_alpha = this%screening_alpha
      else
         allocate (screening_alpha(0))
      end if
      if (allocated(this%screening_sigma)) then
         allocate (screening_sigma, mold=this%screening_sigma)
         screening_sigma = this%screening_sigma
      else
         allocate (screening_sigma(0))
      end if
      if (allocated(this%reduced_acr)) then
         allocate (reduced_acr, mold=this%reduced_acr)
         reduced_acr = this%reduced_acr
      else
         allocate (reduced_acr(0))
      end if
      if (allocated(this%num)) then
         allocate (num, mold=this%num)
         num = this%num
      else
         allocate (num(0))
      end if
      if (allocated(this%no)) then
         allocate (no, mold=this%no)
         no = this%no
      else
         allocate (no(0))
      end if
      if (allocated(this%izpsurf)) then
         allocate (izpsurf, mold=this%izpsurf)
         izpsurf = this%izpsurf
      else
         allocate (izpsurf(0))
      end if
      if (allocated(this%izsurf)) then
         allocate (izsurf, mold=this%izsurf)
         izsurf = this%izsurf
      else
         allocate (izsurf(0))
      end if
      if (allocated(this%nosurf)) then
         allocate (nosurf, mold=this%nosurf)
         nosurf = this%nosurf
      else
         allocate (nosurf(0))
      end if
      if (allocated(this%izpo)) then
         allocate (izpo, mold=this%izpo)
         izpo = this%izpo
      else
         allocate (izpo(0))
      end if
      if (allocated(this%izp)) then
         allocate (izp, mold=this%izp)
         izp = this%izp
      else
         allocate (izp(0))
      end if
      if (allocated(this%iz)) then
         allocate (iz, mold=this%iz)
         iz = this%iz
      else
         allocate (iz(0))
      end if
      if (allocated(this%iu)) then
         allocate (iu, mold=this%iu)
         iu = this%iu
      else
         allocate (iu(0))
      end if
      if (allocated(this%irec)) then
         allocate (irec, mold=this%irec)
         irec = this%irec
      else
         allocate (irec(0))
      end if
      if (allocated(this%ib)) then
         allocate (ib, mold=this%ib)
         ib = this%ib
      else
         allocate (ib(0))
      end if

      ! two dimensional allocatables

      if (allocated(this%primcell)) then
         allocate (primcell, mold=this%primcell)
         primcell = this%primcell
      else
         allocate (primcell(0, 0))
      end if
      if (allocated(this%inclu)) then
         allocate (inclu, mold=this%inclu)
         inclu = this%inclu
      else
         allocate (inclu(0, 0))
      end if
      if (allocated(this%crsurf)) then
         allocate (crsurf, mold=this%crsurf)
         crsurf = this%crsurf
      else
         allocate (crsurf(0, 0))
      end if
      if (allocated(this%crd)) then
         allocate (crd, mold=this%crd)
         crd = this%crd
      else
         allocate (crd(0, 0))
      end if
      if (allocated(this%cr)) then
         allocate (cr, mold=this%cr)
         cr = this%cr
      else
         allocate (cr(0, 0))
      end if
      if (allocated(this%acr)) then
         allocate (acr, mold=this%acr)
         acr = this%acr
      else
         allocate (acr(0, 0))
      end if

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         write (unit, nml=lattice)
      else if (present(file)) then
         open (newunit=newunit, file=file)
         write (newunit, nml=lattice)
         close (newunit)
      else
         write (*, nml=lattice)
      end if

   end subroutine print_state_full

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
   !> @param[in] unit File unit used to write namelist
   !> @param[in] file File name used to write namelist
   !---------------------------------------------------------------------------
   subroutine print_state(this, unit, file)
      implicit none
      class(lattice), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

      include 'include_codes/namelists/lattice.f90'

      ! scalar

      wav = this%wav
      rc = this%rc
      celldm = this%celldm
      alat = this%alat
      nlay = this%nlay
      ndim = this%ndim
      npe = this%npe
      nclu = this%nclu
      crystal_sym = this%crystal_sym
      surftype = this%surftype
      a = this%a

      ! ! one dimensional allocatables

      if (allocated(this%no)) then
         allocate (no, mold=this%no)
         no = this%no
      else
         allocate (no(0))
      end if
      if (allocated(this%izp)) then
         allocate (izp, mold=this%izp)
         izp = this%izp
      else
         allocate (izp(0))
      end if

      ! ! two dimensional allocatables

      if (allocated(this%inclu)) then
         allocate (inclu, mold=this%inclu)
         inclu = this%inclu
      else
         allocate (inclu(0, 0))
      end if
      if (allocated(this%crd)) then
         allocate (crd, mold=this%crd)
         crd = this%crd
      else
         allocate (crd(0, 0))
      end if

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         call this%print_state_formatted(unit=unit)
      else if (present(file)) then
         call this%print_state_formatted(file=file)
      else
         call this%print_state_formatted()
      end if

   end subroutine print_state

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Print class members values in namelist format
   !>
   !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
   !> @param[in] unit File unit used to write namelist
   !> @param[in] file File name used to write namelist
   !---------------------------------------------------------------------------
   subroutine print_state_formatted(this, unit, file)
      class(lattice), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file

      type(namelist_generator) :: nml
      integer :: i

      nml = namelist_generator('lattice')

      ! scalar

      call nml%add('zmin', this%zmin)
      call nml%add('zmax', this%zmax)
      call nml%add('zstep', this%zstep)
      call nml%add('wav', this%wav)
      call nml%add('vol', this%vol)
      call nml%add('rc', this%rc)
      call nml%add('r2', this%r2)
      call nml%add('celldm', this%celldm)
      call nml%add('alat', this%alat)
      call nml%add('reduced_nbas', this%reduced_nbas)
      call nml%add('ntype', this%ntype)
      call nml%add('ntot', this%ntot)
      call nml%add('nrec', this%nrec)
      call nml%add('nmax', this%nmax)
      call nml%add('nlay', this%nlay)
      call nml%add('ndim', this%ndim)
      call nml%add('npe', this%npe)
      call nml%add('nclu', this%nclu)
      call nml%add('nbulk_bulk', this%nbulk_bulk)
      call nml%add('nbulk', this%nbulk)
      call nml%add('nbas', this%nbas)
      call nml%add('kk', this%kk)
      call nml%add('dx', this%dx)
      call nml%add('dy', this%dy)
      call nml%add('dz', this%dz)
      call nml%add('dw', this%dw)
      call nml%add('crystal_sym', this%crystal_sym)
      call nml%add('surftype', this%surftype)
      call nml%add('njij', this%njij)
      call nml%add('njijk', this%njijk)
      call nml%add('strux_backend', trim(this%strux_backend))
      call nml%add('screening', trim(this%screening))
      call nml%add('strux_want_sdot', this%strux_want_sdot)
      call nml%add('strux_solve_scale', this%strux_solve_scale)

      ! one dimensional allocatables
      ! TODO: implement test inside namelist_generator
      if (allocated(this%z)) call nml%add('z', this%z)
      if (allocated(this%ct)) call nml%add('ct', this%ct)
      if (allocated(this%screening_alpha)) call nml%add('screening_alpha', this%screening_alpha)
      if (allocated(this%screening_sigma)) call nml%add('screening_sigma', this%screening_sigma)
      if (allocated(this%reduced_acr)) call nml%add('reduced_acr', this%reduced_acr)
      if (allocated(this%num)) call nml%add('num', this%num)
      if (allocated(this%no)) call nml%add('no', this%no)
      if (allocated(this%izpsurf)) call nml%add('izpsurf', this%izpsurf)
      if (allocated(this%izsurf)) call nml%add('izsurf', this%izsurf)
      if (allocated(this%nosurf)) call nml%add('nosurf', this%nosurf)
      if (allocated(this%izpo)) call nml%add('izpo', this%izpo)
      if (allocated(this%izp)) call nml%add('izp', this%izp)
      if (allocated(this%iz)) call nml%add('iz', this%iz)
      if (allocated(this%iu)) call nml%add('iu', this%iu)
      if (allocated(this%irec)) call nml%add('irec', this%irec)
      if (allocated(this%ib)) call nml%add('ib', this%ib)

      ! ! two dimensional allocatables
      ! TODO: implement test inside namelist_generator
      call nml%add('a', this%a)
      if (allocated(this%primcell)) call nml%add('primcell', this%primcell)
      if (allocated(this%inclu)) call nml%add('inclu', this%inclu)
      if (allocated(this%crsurf)) call nml%add('crsurf', this%crsurf)
      if (allocated(this%crd)) call nml%add('crd', this%crd)
      if (allocated(this%cr)) call nml%add('cr', this%cr)
      if (allocated(this%acr)) call nml%add('acr', this%acr)

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         call nml%generate_namelist(unit=unit)
      else if (present(file)) then
         call nml%generate_namelist(file=file)
      else
         call nml%generate_namelist()
      end if
   end subroutine print_state_formatted

end module lattice_mod
