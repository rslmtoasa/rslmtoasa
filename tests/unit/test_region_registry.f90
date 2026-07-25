!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_region_registry
!
!> @brief Pins the B7.1 region registry (source/region_registry.f90) against
!>        today's charge%surfpot index map, documented in B7 §2.10 and
!>        verified against the sources during B7.0:
!>
!>          - rows 1..init            vacuum-side frozen (init = 6, the
!>                                     historical magic offset)
!>          - rows init+1..init+nlay  active layers
!>          - rows init+nlay+1..nbas  bulk-side frozen
!>
!>        and the deliberate-not-accidental fact that surfpot's charge loop
!>        only ever fills rows init+1..nbas, so the vacuum rows exist in the
!>        geometry but are held at exactly zero deviation charge (B7 §1.5).
!>
!>        Tests:
!>          1. Region assignment matches the index map exactly for a
!>             representative (nbas, nlay, init) matching the fccCu001
!>             surface regression case (nbas=49, nlay=6, init=6).
!>          2. Vacuum and bulk regions are frozen; the active region is not.
!>          3. z is carried as DATA -- the registry's z array is bit-identical
!>             to the input z array, not recomputed from layer_index (B7
!>             §2.10: this is what makes relaxed-z later a parameter change).
!>          4. w (Wigner-Seitz radius) is carried per-site, not collapsed to
!>             a single system-wide value.
!>          5. Dump routine runs without error and produces one line per site.
!>
!>        Exits non-zero (error stop) on any failure so ctest registers a fail.
!------------------------------------------------------------------------------
program test_region_registry
   use precision_mod, only: rp
   use region_registry_mod, only: region_registry, region_kind_vacuum, region_kind_bulk, region_kind_active
   implicit none

   logical :: failed

   failed = .false.

   call test_index_map_matches_surfpot()
   call test_frozen_active_mask()
   call test_z_is_data_not_derived()
   call test_per_site_w()
   call test_dump()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !> 1. Region id / active mask / layer_index reproduce surfpot's index map
   !> exactly for the fccCu001 surface regression geometry (nbas=49, nlay=6,
   !> init=6, per example/surface/fccCu001/Test/input.nml and
   !> lattice_cluster.f90's build_surf_full).
   subroutine test_index_map_matches_surfpot()
      type(region_registry) :: reg
      integer, parameter :: nbas = 49, nlay = 6, init = 6
      real(rp), dimension(nbas) :: z, w
      integer, dimension(nbas) :: reftype
      integer :: i, expect_layer

      do i = 1, nbas
         z(i) = real(i, rp)*0.7_rp   ! arbitrary but distinguishable per-site z
         w(i) = 2.66_rp
         reftype(i) = 0
      end do

      call reg%build_from_buildsurf(nbas, nlay, z, w, reftype, init)

      if (reg%nsite /= nbas) then
         write (*, '(a,i0)') 'FAIL nsite mismatch: ', reg%nsite
         failed = .true.
      end if

      do i = 1, nbas
         if (i <= init) then
            if (reg%region(reg%region_id(i))%kind /= region_kind_vacuum) then
               write (*, '(a,i0)') 'FAIL site expected vacuum, site=', i
               failed = .true.
            end if
            if (reg%active(i)) then
               write (*, '(a,i0)') 'FAIL vacuum site marked active, site=', i
               failed = .true.
            end if
            if (reg%layer_index(i) /= 0) then
               write (*, '(a,i0)') 'FAIL vacuum site has nonzero layer_index, site=', i
               failed = .true.
            end if
         else if (i <= init + nlay) then
            if (reg%region(reg%region_id(i))%kind /= region_kind_active) then
               write (*, '(a,i0)') 'FAIL site expected active region, site=', i
               failed = .true.
            end if
            if (.not. reg%active(i)) then
               write (*, '(a,i0)') 'FAIL active-layer site not marked active, site=', i
               failed = .true.
            end if
            expect_layer = i - init
            if (reg%layer_index(i) /= expect_layer) then
               write (*, '(a,i0,a,i0,a,i0)') 'FAIL layer_index site=', i, ' got=', reg%layer_index(i), &
                  ' expect=', expect_layer
               failed = .true.
            end if
         else
            if (reg%region(reg%region_id(i))%kind /= region_kind_bulk) then
               write (*, '(a,i0)') 'FAIL site expected bulk, site=', i
               failed = .true.
            end if
            if (reg%active(i)) then
               write (*, '(a,i0)') 'FAIL bulk site marked active, site=', i
               failed = .true.
            end if
         end if
      end do

      if (.not. failed) write (*, '(a)') 'ok  index map matches surfpot (init=6, nlay=6, nbas=49)'
   end subroutine test_index_map_matches_surfpot

   !> 2. Frozen/active mask: vacuum and bulk regions are frozen (parameters
   !> never updated); the active region is not. This is the registry-explicit
   !> statement of "the charge loop fills only rows init+1..nbas, so the
   !> vacuum rows exist in the geometry but are held at exactly zero
   !> deviation charge" (B7 §1.5) -- preserved DELIBERATELY, not by omission.
   subroutine test_frozen_active_mask()
      type(region_registry) :: reg
      integer, parameter :: nbas = 20, nlay = 4, init = 3
      real(rp), dimension(nbas) :: z, w
      integer :: i, ir

      z = 0.0_rp; w = 2.5_rp
      do i = 1, nbas
         z(i) = real(i, rp)
      end do

      call reg%build_from_buildsurf(nbas, nlay, z, w, init=init)

      do ir = 1, reg%nregion
         select case (reg%region(ir)%kind)
         case (region_kind_vacuum, region_kind_bulk)
            if (.not. reg%region(ir)%frozen) then
               write (*, '(a,i0)') 'FAIL vacuum/bulk region not marked frozen, region=', ir
               failed = .true.
            end if
         case (region_kind_active)
            if (reg%region(ir)%frozen) then
               write (*, '(a)') 'FAIL active region marked frozen'
               failed = .true.
            end if
         end select
      end do

      ! Every frozen site's active() must be .false., matching region%frozen.
      do i = 1, nbas
         if (reg%region(reg%region_id(i))%frozen .eqv. reg%active(i)) then
            write (*, '(a,i0)') 'FAIL active(i) inconsistent with region%frozen at site=', i
            failed = .true.
         end if
      end do

      if (.not. failed) write (*, '(a)') 'ok  vacuum/bulk frozen, active region not frozen'
   end subroutine test_frozen_active_mask

   !> 3. z is carried as DATA. Perturb z at a few sites (as a stand-in for a
   !> future relaxed-z geometry) and confirm the registry reproduces exactly
   !> the perturbed values, not a value derived from layer spacing/index --
   !> this is what B7 §2.10 means by "relaxed-z later is a parameter change".
   subroutine test_z_is_data_not_derived()
      type(region_registry) :: reg
      integer, parameter :: nbas = 12, nlay = 3, init = 4
      real(rp), dimension(nbas) :: z, w
      integer :: i

      w = 3.0_rp
      do i = 1, nbas
         z(i) = real(i, rp)*1.5_rp
      end do
      ! Perturb one active-layer site off its "regular" spacing, as relaxed-z
      ! would.
      z(init + 2) = z(init + 2) + 0.1234_rp

      call reg%build_from_buildsurf(nbas, nlay, z, w, init=init)

      do i = 1, nbas
         if (reg%z(i) /= z(i)) then
            write (*, '(a,i0,2es20.10)') 'FAIL registry z not bit-identical to input at site=', i, reg%z(i), z(i)
            failed = .true.
         end if
      end do

      if (.not. failed) write (*, '(a)') 'ok  z carried as data, bit-identical including a perturbed site'
   end subroutine test_z_is_data_not_derived

   !> 4. w (Wigner-Seitz radius) is per-site, not a single system-wide value
   !> collapsed across the cluster -- this is the ingredient B7.0/C5 flags as
   !> required for two regions at different w.
   subroutine test_per_site_w()
      type(region_registry) :: reg
      integer, parameter :: nbas = 8, nlay = 2, init = 3
      real(rp), dimension(nbas) :: z, w
      integer :: i

      do i = 1, nbas
         z(i) = real(i, rp)
      end do
      w = 2.66_rp
      w(1) = 3.05_rp   ! a different region's w, e.g. a vacuum empty sphere

      call reg%build_from_buildsurf(nbas, nlay, z, w, init=init)

      if (reg%w(1) /= 3.05_rp) then
         write (*, '(a)') 'FAIL per-site w not preserved at site 1'
         failed = .true.
      end if
      do i = 2, nbas
         if (reg%w(i) /= 2.66_rp) then
            write (*, '(a,i0)') 'FAIL per-site w not preserved at site=', i
            failed = .true.
         end if
      end do

      if (.not. failed) write (*, '(a)') 'ok  w carried per-site, not collapsed to a system-wide average'
   end subroutine test_per_site_w

   !> 5. Dump routine: runs without error, writes one data line per site plus
   !> a region summary header. Just exercises the debugging path end-to-end.
   subroutine test_dump()
      type(region_registry) :: reg
      integer, parameter :: nbas = 10, nlay = 3, init = 3
      real(rp), dimension(nbas) :: z, w
      character(len=*), parameter :: fname = 'test_region_registry_dump.out'
      integer :: i, unit, ios, nlines
      character(len=256) :: line
      logical :: exists

      do i = 1, nbas
         z(i) = real(i, rp)
      end do
      w = 2.66_rp

      call reg%build_from_buildsurf(nbas, nlay, z, w, init=init)
      call reg%dump(fname=fname)

      inquire (file=fname, exist=exists)
      if (.not. exists) then
         write (*, '(a)') 'FAIL dump did not create output file'
         failed = .true.
         return
      end if

      nlines = 0
      open (newunit=unit, file=fname, status='old', action='read')
      do
         read (unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         nlines = nlines + 1
      end do
      close (unit, status='delete')

      ! At least one line per site, plus header/region-summary lines.
      if (nlines < nbas) then
         write (*, '(a,i0)') 'FAIL dump produced too few lines: ', nlines
         failed = .true.
      end if

      if (.not. failed) write (*, '(a,i0,a)') 'ok  dump routine wrote ', nlines, ' lines'
   end subroutine test_dump

end program test_region_registry
