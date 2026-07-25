!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Region Registry
!
! DESCRIPTION:
!> @brief
!> Explicit per-site region bookkeeping for the embedding cluster (B7.1).
!>
!> Today's `charge%surfpot` locates the active/frozen boundary implicitly:
!> `init = 6` is a magic offset encoding the number of leading (vacuum-side)
!> frozen rows, and `iq = init + ibas` maps a layer index onto a matrix row.
!> That works only because the buildsurf layout is fixed (one vacuum side,
!> one bulk side, `init` hardcoded). B7 (interfaces and vacuum leads) needs
!> more than one frozen reference region to coexist in a single cluster, so
!> the offset arithmetic has to become data.
!>
!> This module supplies that data as an explicit registry, one entry per
!> cluster site (dimension `nbas`, matching `charge%qz`/`charge%wssurf`):
!>   - `region_id`      which frozen/active region the site belongs to
!>   - `active`         .true. for sites that relax self-consistently
!>   - `layer_index`    physical atomic layer (1..nlay for active sites)
!>   - `z`               the site's z coordinate, carried as DATA (not
!>                        derived from layer_index) so relaxed-z is later a
!>                        parameter change, not a rewrite (B7 §2.10)
!>   - `w`               per-site Wigner-Seitz radius
!>   - `reference_type`  which frozen parameter class a site reverts to when
!>                        frozen (index into the region's reference-type
!>                        table; mirrors lattice%chargetrf_type)
!>
!> `build_from_buildsurf` reproduces today's buildsurf index map exactly, as
!> a registry instance, so it can be validated against the surface regression
!> without touching `charge%surfpot` itself (which remains the permanent
!> regression oracle -- B7 explicitly forbids changing its behaviour). Later
!> work (B7.3/B7.4) is expected to consume a registry built by a genuinely
!> two-sided constructor; this module does not gate that on rewiring surfpot.
!------------------------------------------------------------------------------
module region_registry_mod
   use precision_mod, only: rp
   use logger_mod, only: g_logger
   use string_mod, only: fmt
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Region kinds. VACUUM and BULK are the two ends of today's buildsurf
   !> slab; LEAD_A/LEAD_B are reserved for the two-sided A|B interface
   !> geometry (B7.3+) so the registry does not need to grow a new field
   !> when that lands.
   integer, parameter, public :: region_kind_vacuum = 1
   integer, parameter, public :: region_kind_bulk = 2
   integer, parameter, public :: region_kind_active = 3
   integer, parameter, public :: region_kind_lead_a = 4
   integer, parameter, public :: region_kind_lead_b = 5

   !> One row per region: metadata shared by every site assigned to it.
   type, public :: region_descriptor
      integer :: kind = region_kind_active
      character(len=32) :: name = ''
      !> .true. if sites of this region are held frozen (never updated).
      logical :: frozen = .false.
   end type region_descriptor

   !> Explicit per-site region registry, dimension nbas.
   type, public :: region_registry
      !> Number of sites (== charge%lattice%nbas for a buildsurf-derived
      !> instance; the extended interaction zone, user-sized).
      integer :: nsite = 0
      !> Number of distinct regions.
      integer :: nregion = 0
      !> Region metadata, dimension(nregion).
      type(region_descriptor), dimension(:), allocatable :: region

      !> Per-site region id, 1..nregion.
      integer, dimension(:), allocatable :: region_id
      !> Per-site active/frozen mask. .true. = active (relaxes self-
      !> consistently); .false. = frozen (parameters imported, never
      !> updated). This is the registry-explicit form of "the charge loop
      !> fills only rows init+1..nbas": frozen sites are deliberately held
      !> at zero deviation charge, not accidentally skipped.
      logical, dimension(:), allocatable :: active
      !> Physical atomic layer index. Meaningful (1..nlay) for active
      !> sites; 0 for frozen sites, which are not part of the layer stack.
      integer, dimension(:), allocatable :: layer_index
      !> Site z coordinate, carried as DATA -- not derived from
      !> layer_index -- so relaxed-z is a parameter change later (B7 §2.10).
      real(rp), dimension(:), allocatable :: z
      !> Per-site Wigner-Seitz radius (mirrors charge%wssurf).
      real(rp), dimension(:), allocatable :: w
      !> Reference-type index a frozen site reverts to (mirrors
      !> lattice%chargetrf_type); 0 for active sites, which have no frozen
      !> reference.
      integer, dimension(:), allocatable :: reference_type
   contains
      procedure :: restore_to_default => region_registry_restore_to_default
      procedure :: build_from_buildsurf => region_registry_build_from_buildsurf
      procedure :: dump => region_registry_dump
      procedure :: is_active => region_registry_is_active
      procedure :: region_of => region_registry_region_of
      final :: region_registry_destructor
   end type region_registry

   interface region_registry
      procedure :: region_registry_constructor
   end interface region_registry

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor. Produces an empty (zero-site) registry; call
   !> build_from_buildsurf (or a future two-sided builder) to populate it.
   !---------------------------------------------------------------------------
   function region_registry_constructor() result(obj)
      type(region_registry) :: obj

      call obj%restore_to_default()
   end function region_registry_constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor.
   !---------------------------------------------------------------------------
   subroutine region_registry_destructor(this)
      type(region_registry) :: this

      if (allocated(this%region)) deallocate (this%region)
      if (allocated(this%region_id)) deallocate (this%region_id)
      if (allocated(this%active)) deallocate (this%active)
      if (allocated(this%layer_index)) deallocate (this%layer_index)
      if (allocated(this%z)) deallocate (this%z)
      if (allocated(this%w)) deallocate (this%w)
      if (allocated(this%reference_type)) deallocate (this%reference_type)
   end subroutine region_registry_destructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset to an empty, zero-site registry.
   !---------------------------------------------------------------------------
   subroutine region_registry_restore_to_default(this)
      class(region_registry), intent(inout) :: this

      if (allocated(this%region)) deallocate (this%region)
      if (allocated(this%region_id)) deallocate (this%region_id)
      if (allocated(this%active)) deallocate (this%active)
      if (allocated(this%layer_index)) deallocate (this%layer_index)
      if (allocated(this%z)) deallocate (this%z)
      if (allocated(this%w)) deallocate (this%w)
      if (allocated(this%reference_type)) deallocate (this%reference_type)

      this%nsite = 0
      this%nregion = 0
      allocate (this%region(0))
      allocate (this%region_id(0))
      allocate (this%active(0))
      allocate (this%layer_index(0))
      allocate (this%z(0))
      allocate (this%w(0))
      allocate (this%reference_type(0))
   end subroutine region_registry_restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Build a registry that reproduces today's buildsurf/surfpot layout
   !> exactly, as an explicit instance, so it can be validated against the
   !> surface regression without touching surfpot itself.
   !>
   !> Reproduces the index map documented in B7 §2.10 and verified against
   !> charge%surfpot (source/charge.f90):
   !>   - rows 1..init            vacuum-side frozen rows (init = 6, the
   !>                              historical magic offset)
   !>   - rows init+1..init+nlay  active layers (charge loop fills these)
   !>   - rows init+nlay+1..nbas  bulk-side frozen rows
   !> The charge loop in surfpot only ever touches rows init+1..nbas, so the
   !> vacuum rows 1..init exist in the geometry but are deliberately held at
   !> exactly zero deviation charge (B7 §1.5). That is encoded here directly:
   !> the vacuum region is frozen, not merely "not yet visited".
   !>
   !> @param[in] nbas   Total number of sites in the extended interaction
   !>                    zone (charge%lattice%nbas). User-sized; no width is
   !>                    assumed.
   !> @param[in] nlay   Number of active physical layers
   !>                    (charge%lattice%nlay).
   !> @param[in] z      Per-site z coordinate, dimension(nbas)
   !>                    (charge%qz, already fully z-general -- madl2d runs
   !>                    off QPPZ = QZ(IQ)-QZ(JQ) for arbitrary pairs).
   !> @param[in] w      Per-site Wigner-Seitz radius, dimension(nbas)
   !>                    (charge%wssurf, already per-site).
   !> @param[in] reference_type  Per-site frozen-reference class index,
   !>                    dimension(nbas) (lattice%chargetrf_type). Optional;
   !>                    if absent, reference_type is left 0 for every site
   !>                    (no reference table available in that context).
   !> @param[in] init    Vacuum-side frozen row count. Optional, default 6
   !>                    -- the value hardcoded in charge%surfpot today.
   !---------------------------------------------------------------------------
   subroutine region_registry_build_from_buildsurf(this, nbas, nlay, z, w, reference_type, init)
      class(region_registry), intent(inout) :: this
      integer, intent(in) :: nbas
      integer, intent(in) :: nlay
      real(rp), dimension(:), intent(in) :: z
      real(rp), dimension(:), intent(in) :: w
      integer, dimension(:), intent(in), optional :: reference_type
      integer, intent(in), optional :: init
      integer :: init_
      integer :: i, ilay
      integer, parameter :: reg_vacuum = 1
      integer, parameter :: reg_active = 2
      integer, parameter :: reg_bulk = 3

      init_ = 6
      if (present(init)) init_ = init

      call this%restore_to_default()

      this%nsite = nbas
      this%nregion = 3

      deallocate (this%region)
      allocate (this%region(this%nregion))
      this%region(reg_vacuum) = region_descriptor(kind=region_kind_vacuum, name='vacuum', frozen=.true.)
      this%region(reg_active) = region_descriptor(kind=region_kind_active, name='active', frozen=.false.)
      this%region(reg_bulk) = region_descriptor(kind=region_kind_bulk, name='bulk', frozen=.true.)

      deallocate (this%region_id, this%active, this%layer_index, this%z, this%w, this%reference_type)
      allocate (this%region_id(nbas))
      allocate (this%active(nbas))
      allocate (this%layer_index(nbas))
      allocate (this%z(nbas))
      allocate (this%w(nbas))
      allocate (this%reference_type(nbas))

      this%z(1:nbas) = z(1:nbas)
      this%w(1:nbas) = w(1:nbas)
      this%reference_type(:) = 0
      if (present(reference_type)) this%reference_type(1:nbas) = reference_type(1:nbas)

      ! Rows 1..init: vacuum-side frozen. Held at exactly zero deviation
      ! charge deliberately (B7 §1.5) -- this is the explicit statement of
      ! what today's surfpot does implicitly by never writing tdq() there.
      do i = 1, min(init_, nbas)
         this%region_id(i) = reg_vacuum
         this%active(i) = .false.
         this%layer_index(i) = 0
      end do

      ! Rows init+1..init+nlay: the active layers.
      ilay = 0
      do i = init_ + 1, min(init_ + nlay, nbas)
         ilay = ilay + 1
         this%region_id(i) = reg_active
         this%active(i) = .true.
         this%layer_index(i) = ilay
      end do

      ! Rows init+nlay+1..nbas: bulk-side frozen.
      do i = init_ + nlay + 1, nbas
         this%region_id(i) = reg_bulk
         this%active(i) = .false.
         this%layer_index(i) = 0
      end do
   end subroutine region_registry_build_from_buildsurf

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Whether site i is active (relaxes self-consistently).
   !---------------------------------------------------------------------------
   pure logical function region_registry_is_active(this, i) result(res)
      class(region_registry), intent(in) :: this
      integer, intent(in) :: i
      res = this%active(i)
   end function region_registry_is_active

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Region descriptor for site i.
   !---------------------------------------------------------------------------
   pure function region_registry_region_of(this, i) result(res)
      class(region_registry), intent(in) :: this
      integer, intent(in) :: i
      type(region_descriptor) :: res
      res = this%region(this%region_id(i))
   end function region_registry_region_of

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Debug dump: one line per site, plus a region summary. Writes to the
   !> given unit if open, otherwise to a fresh file 'region_registry.out'.
   !---------------------------------------------------------------------------
   subroutine region_registry_dump(this, unit, fname)
      class(region_registry), intent(in) :: this
      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: fname
      integer :: iu, i, ir
      character(len=64) :: fname_
      logical :: opened_here

      opened_here = .false.
      if (present(unit)) then
         iu = unit
      else
         fname_ = 'region_registry.out'
         if (present(fname)) fname_ = fname
         open (newunit=iu, file=trim(fname_), action='write', status='replace')
         opened_here = .true.
      end if

      write (iu, '(a)') '# Region registry dump (B7.1)'
      write (iu, '(a,i0,a,i0)') '# nsite= ', this%nsite, '  nregion= ', this%nregion
      write (iu, '(a)') '#'
      write (iu, '(a)') '# regions:'
      do ir = 1, this%nregion
         write (iu, '(a,i0,a,a,a,l1)') '#   id=', ir, ' name=', trim(this%region(ir)%name), &
            ' frozen=', this%region(ir)%frozen
      end do
      write (iu, '(a)') '#'
      write (iu, '(a)') '#  site  region_id  region_name  active  layer_index'// &
         '            z               w  reference_type'
      do i = 1, this%nsite
         ir = this%region_id(i)
         write (iu, '(2x,i5,4x,i5,4x,a12,4x,l1,6x,i5,2x,f14.8,2x,f14.8,4x,i5)') &
            i, ir, trim(this%region(ir)%name), this%active(i), this%layer_index(i), &
            this%z(i), this%w(i), this%reference_type(i)
      end do

      if (opened_here) close (iu)
   end subroutine region_registry_dump

end module region_registry_mod
