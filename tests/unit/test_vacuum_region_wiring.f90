!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
! PROGRAM: test_vacuum_region_wiring
!
!> @brief B7.6: pins the *wiring* that makes `source/vacuum_lead.f90` reachable
!>        from the `buildinterface` (calctype='L') path.
!>
!>        ## Why this is a separate file from test_vacuum_lead.f90
!>
!>        `test_vacuum_lead` checks the GENERATOR: that empty-lattice
!>        parameters produced by the code's own radial solver agree with the
!>        analytic spherical-Bessel oracle. That test passed for the entire
!>        period during which the generator had no caller anywhere in
!>        production code -- a correct component wired to nothing.
!>
!>        This file checks the opposite thing: that a vacuum region can be
!>        CONSTRUCTED on the interface path at all, and that the vacuum-aware
!>        behaviour already written into the consumers actually engages when it
!>        is. A component test cannot catch an integration gap; that gap is
!>        precisely what B7.6 found.
!>
!>        ## What is pinned
!>
!>          1. `build_from_interface` defaults to the metallic A | B layout
!>             (`region_kind_lead_b`), so the pre-B7.6 behaviour is unchanged
!>             when `kind_b` is not passed. Backwards compatibility is a
!>             standing B7 requirement (G-B7-2).
!>          2. Passing `kind_b = region_kind_vacuum` produces a vacuum region,
!>             and names it 'vacuum' rather than 'B' -- the name is what every
!>             alignment diagnostic prints and what `fix_fermi_to_region`
!>             matches on.
!>          3. The vacuum region is FROZEN. "Frozen at a shifted value" is a
!>             change of reference level, not a relaxation (B7 §1.3); vacuum
!>             sites must never acquire self-consistent charge.
!>          4. `gauge_anchor` skips the vacuum region. Vacuum carries no states
!>             at E_F, so "deep vacuum is neutral bulk vacuum" cannot drive a
!>             residual and must not be the gauge anchor (B7 §1.3). With region
!>             B vacuum the anchor must come out region A.
!>          5. The rigid-shift property the self-consistent vacuum level relies
!>             on: regenerating at V0 + dV shifts enu, C and VL by exactly dV
!>             and leaves srdel, qpar, ppar invariant. `refresh_vacuum_region`
!>             re-runs the generator at the solved level each iteration rather
!>             than shifting parameters in place, so this must hold in the
!>             generator for that to be equivalent.
!>
!>        Exits non-zero (error stop) on any failure so ctest registers a fail.
!------------------------------------------------------------------------------
program test_vacuum_region_wiring
   use precision_mod, only: rp
   use region_registry_mod, only: region_registry, region_kind_vacuum, &
                                  region_kind_lead_a, region_kind_lead_b, region_kind_active
   use vacuum_lead_mod, only: vacuum_lead
   use self_mod, only: self
   implicit none

   logical :: failed

   failed = .false.

   call test_default_is_metallic()
   call test_vacuum_kind_and_name()
   call test_vacuum_region_is_frozen()
   call test_gauge_anchor_skips_vacuum()
   call test_rigid_shift_property()

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   !> Build a representative interface registry. Geometry matches the
   !> fccCu111 interface examples: 49 Madelung rows, one frozen row each side.
   subroutine build(reg, kind_b)
      type(region_registry), intent(out) :: reg
      integer, intent(in), optional :: kind_b
      integer, parameter :: nbas = 49, nlay_a = 1, nlay_active = 47
      real(rp), dimension(nbas) :: z, w
      integer :: i

      do i = 1, nbas
         z(i) = real(i, rp)*0.7_rp
      end do
      w = 2.66899_rp

      if (present(kind_b)) then
         call reg%build_from_interface(nbas, nlay_a, nlay_active, z, w, kind_b=kind_b)
      else
         call reg%build_from_interface(nbas, nlay_a, nlay_active, z, w)
      end if
   end subroutine build

   !> 1. Omitting kind_b leaves the metallic A | B layout untouched.
   subroutine test_default_is_metallic()
      type(region_registry) :: reg

      call build(reg)

      if (reg%region(3)%kind /= region_kind_lead_b) then
         write (*, '(a,i0)') 'FAIL default region B kind is not lead_b, got ', reg%region(3)%kind
         failed = .true.
      end if
      if (trim(reg%region(3)%name) /= 'B') then
         write (*, '(a)') 'FAIL default region B name is not ''B'', got '//trim(reg%region(3)%name)
         failed = .true.
      end if

      if (.not. failed) write (*, '(a)') 'ok  default build_from_interface is still metallic A | B'
   end subroutine test_default_is_metallic

   !> 2. kind_b = vacuum produces a vacuum region, named 'vacuum'.
   subroutine test_vacuum_kind_and_name()
      type(region_registry) :: reg

      call build(reg, kind_b=region_kind_vacuum)

      if (reg%region(3)%kind /= region_kind_vacuum) then
         write (*, '(a,i0)') 'FAIL region B kind is not vacuum, got ', reg%region(3)%kind
         failed = .true.
      end if
      if (trim(reg%region(3)%name) /= 'vacuum') then
         write (*, '(a)') 'FAIL vacuum region name is not ''vacuum'', got '//trim(reg%region(3)%name)
         failed = .true.
      end if
      ! Region A and the active zone must be untouched by the kind switch.
      if (reg%region(1)%kind /= region_kind_lead_a) then
         write (*, '(a)') 'FAIL region A kind changed when region B became vacuum'
         failed = .true.
      end if
      if (reg%region(2)%kind /= region_kind_active) then
         write (*, '(a)') 'FAIL active region kind changed when region B became vacuum'
         failed = .true.
      end if

      if (.not. failed) write (*, '(a)') 'ok  kind_b = vacuum builds a vacuum region named ''vacuum'''
   end subroutine test_vacuum_kind_and_name

   !> 3. Vacuum sites are frozen -- never self-consistent (B7 §1.3).
   subroutine test_vacuum_region_is_frozen()
      type(region_registry) :: reg
      integer :: i

      call build(reg, kind_b=region_kind_vacuum)

      if (.not. reg%region(3)%frozen) then
         write (*, '(a)') 'FAIL vacuum region is not marked frozen'
         failed = .true.
      end if

      do i = 1, reg%nsite
         if (reg%region_id(i) == 3 .and. reg%active(i)) then
            write (*, '(a,i0)') 'FAIL vacuum site marked active at row ', i
            failed = .true.
         end if
      end do

      if (.not. failed) write (*, '(a)') 'ok  vacuum region and all its sites are frozen'
   end subroutine test_vacuum_region_is_frozen

   !> 4. The gauge anchor must be region A, not the vacuum region (B7 §1.3):
   !> vacuum has no states at E_F, so it cannot drive an alignment residual.
   subroutine test_gauge_anchor_skips_vacuum()
      type(region_registry) :: reg
      integer :: ianchor

      call build(reg, kind_b=region_kind_vacuum)
      ianchor = reg%gauge_anchor()

      if (ianchor /= 1) then
         write (*, '(a,i0)') 'FAIL gauge anchor is not region A (1), got ', ianchor
         failed = .true.
      end if
      if (ianchor >= 1) then
         if (reg%region(ianchor)%kind == region_kind_vacuum) then
            write (*, '(a)') 'FAIL gauge anchor landed on the vacuum region'
            failed = .true.
         end if
      end if

      if (.not. failed) write (*, '(a)') 'ok  gauge anchor skips vacuum and selects region A'
   end subroutine test_gauge_anchor_skips_vacuum

   !> 5. Rigid shift. `refresh_vacuum_region` regenerates at the solved vacuum
   !> level rather than shifting stored parameters, so the generator itself
   !> must carry the shift exactly: enu/C/VL move by dV, srdel/qpar/ppar do not.
   subroutine test_rigid_shift_property()
      type(vacuum_lead) :: v0set, v1set
      type(self) :: solver
      real(rp), parameter :: ws = 2.66899_rp, dv = 0.15_rp
      real(rp), parameter :: tol = 1.0e-9_rp
      integer :: l, isp

      v0set = vacuum_lead(2)
      v1set = vacuum_lead(2)

      call v0set%generate(solver, ws, 0.0_rp)
      call v1set%generate(solver, ws, dv)

      do isp = 1, 2
         do l = 0, 2
            if (abs((v1set%enu(l, isp) - v0set%enu(l, isp)) - dv) > tol) then
               write (*, '(a,2i3)') 'FAIL enu did not shift rigidly at l,isp=', l, isp
               failed = .true.
            end if
            if (abs((v1set%c(l, isp) - v0set%c(l, isp)) - dv) > tol) then
               write (*, '(a,2i3)') 'FAIL C did not shift rigidly at l,isp=', l, isp
               failed = .true.
            end if
            if (abs(v1set%srdel(l, isp) - v0set%srdel(l, isp)) > tol) then
               write (*, '(a,2i3)') 'FAIL srdel changed under a rigid shift at l,isp=', l, isp
               failed = .true.
            end if
            if (abs(v1set%qpar(l, isp) - v0set%qpar(l, isp)) > tol) then
               write (*, '(a,2i3)') 'FAIL qpar changed under a rigid shift at l,isp=', l, isp
               failed = .true.
            end if
            if (abs(v1set%ppar(l, isp) - v0set%ppar(l, isp)) > tol) then
               write (*, '(a,2i3)') 'FAIL ppar changed under a rigid shift at l,isp=', l, isp
               failed = .true.
            end if
         end do
      end do

      if (.not. failed) then
         write (*, '(a)') 'ok  rigid shift: enu/C move by dV, srdel/qpar/ppar invariant'
      end if
   end subroutine test_rigid_shift_property

end program test_vacuum_region_wiring
