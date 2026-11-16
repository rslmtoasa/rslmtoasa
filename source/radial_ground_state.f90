!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Radial ground-state contract
!
!> Accepted radial state of one site, captured at the end of ATOMSC/NEWRHO and
!> read by the ALSDA radial response (alsda_radial.f90).
!------------------------------------------------------------------------------

module radial_ground_state_mod
   use precision_mod, only: rp
   implicit none
   private

   real(rp), parameter, public :: RADIAL_PI = 3.1415926535897932384626433832795_rp

   type, public :: radial_ground_state
      ! Set by capture.
      logical :: valid = .false.

      ! Mesh metadata.  r is in bohr and A/B reproduce dr/di=A*(r+B).
      real(rp) :: a = 0.0_rp
      real(rp) :: b = 0.0_rp
      real(rp), allocatable :: r(:)
      ! Accepted SR number densities: n_up-n_down defines the ALSDA derivative.
      real(rp), allocatable :: n_up(:), n_down(:)
      ! Vxc_up - Vxc_down (Ry).
      real(rp), allocatable :: delta_vxc(:)

      ! Accepted scalar-relativistic radial solutions per (r, l+1, spin) at the
      ! accepted potential parameters, and their linearisation energies (Ry).
      real(rp), allocatable :: pauli_large(:, :, :), pauli_small(:, :, :)
      real(rp), allocatable :: sr_large_dot(:, :, :), sr_small_dot(:, :, :)
      real(rp), allocatable :: sr_gfac(:, :, :)
      real(rp), allocatable :: sr_enu(:, :)
      integer :: lmax = -1
   contains
      procedure :: clear
      procedure :: capture
      procedure :: begin_pauli_basis
      procedure :: set_basis_channel
      final :: finalize_radial_ground_state
   end type radial_ground_state

   public :: radial_simpson_weight

contains

   pure function radial_simpson_weight(ir, nr) result(weight)
      integer, intent(in) :: ir, nr
      real(rp) :: weight

      weight = 2.0_rp*(mod(ir + 1, 2) + 1)/3.0_rp
      if (ir == 1 .or. ir == nr) weight = 1.0_rp/3.0_rp
   end function radial_simpson_weight

   subroutine clear(this)
      class(radial_ground_state), intent(inout) :: this

      if (allocated(this%r)) deallocate(this%r)
      if (allocated(this%n_up)) deallocate(this%n_up)
      if (allocated(this%n_down)) deallocate(this%n_down)
      if (allocated(this%delta_vxc)) deallocate(this%delta_vxc)
      if (allocated(this%pauli_large)) deallocate(this%pauli_large)
      if (allocated(this%pauli_small)) deallocate(this%pauli_small)
      if (allocated(this%sr_large_dot)) deallocate(this%sr_large_dot)
      if (allocated(this%sr_small_dot)) deallocate(this%sr_small_dot)
      if (allocated(this%sr_gfac)) deallocate(this%sr_gfac)
      if (allocated(this%sr_enu)) deallocate(this%sr_enu)
      this%valid = .false.
      this%a = 0.0_rp
      this%b = 0.0_rp
      this%lmax = -1
   end subroutine clear

   subroutine capture(this, a, b, rofi, rho, density_origin, vxc_up, vxc_down)
      class(radial_ground_state), intent(inout) :: this
      real(rp), intent(in) :: a, b
      real(rp), intent(in) :: rofi(:)
      real(rp), intent(in) :: rho(:, :)
      real(rp), intent(in) :: density_origin(2)
      real(rp), intent(in) :: vxc_up(:), vxc_down(:)

      integer :: ir, nr
      real(rp) :: four_pi_r2

      nr = size(rofi)
      if (nr < 3 .or. size(rho, 1) /= nr .or. size(rho, 2) /= 2 .or. &
          size(vxc_up) /= nr .or. size(vxc_down) /= nr) then
         error stop 'radial_ground_state%capture: inconsistent radial dimensions'
      end if

      call this%clear()
      allocate(this%r(nr), this%n_up(nr), this%n_down(nr), this%delta_vxc(nr))

      this%a = a
      this%b = b
      this%r = rofi
      this%delta_vxc = vxc_up - vxc_down
      this%valid = .true.

      ! The first mesh point is r=0.  VXC0SP obtains the regular physical
      ! density there by extrapolation; every other point follows exactly
      ! n_sigma=RHO_sigma/(4*pi*r^2).
      this%n_up(1) = density_origin(1)
      this%n_down(1) = density_origin(2)
      do ir = 2, nr
         four_pi_r2 = 4.0_rp*RADIAL_PI*this%r(ir)**2
         this%n_up(ir) = rho(ir, 1)/four_pi_r2
         this%n_down(ir) = rho(ir, 2)/four_pi_r2
      end do
   end subroutine capture

   subroutine begin_pauli_basis(this, lmax)
      class(radial_ground_state), intent(inout) :: this
      integer, intent(in) :: lmax

      integer :: nr

      if (.not. this%valid .or. .not. allocated(this%r)) then
         error stop 'radial_ground_state%begin_pauli_basis: no captured radial state'
      end if
      if (lmax < 0) error stop 'radial_ground_state%begin_pauli_basis: invalid lmax'

      if (allocated(this%pauli_large)) deallocate(this%pauli_large, this%pauli_small)
      if (allocated(this%sr_large_dot)) deallocate(this%sr_large_dot, this%sr_small_dot, this%sr_gfac, this%sr_enu)
      nr = size(this%r)
      allocate(this%pauli_large(nr, lmax + 1, 2), this%pauli_small(nr, lmax + 1, 2))
      allocate(this%sr_large_dot(nr, lmax + 1, 2), this%sr_small_dot(nr, lmax + 1, 2), &
               this%sr_gfac(nr, lmax + 1, 2), this%sr_enu(lmax + 1, 2))
      this%pauli_large = 0.0_rp
      this%pauli_small = 0.0_rp
      this%sr_large_dot = 0.0_rp
      this%sr_small_dot = 0.0_rp
      this%sr_gfac = 1.0_rp
      this%sr_enu = 0.0_rp
      this%lmax = lmax
   end subroutine begin_pauli_basis

   !> Store one (l, spin) channel of the RSEQSR/PHDFSR solution used by NEWRHO.
   subroutine set_basis_channel(this, l, ispin, enu, large, small, large_dot, small_dot, gfac)
      class(radial_ground_state), intent(inout) :: this
      integer, intent(in) :: l, ispin
      real(rp), intent(in) :: enu
      real(rp), intent(in) :: large(:), small(:), large_dot(:), small_dot(:), gfac(:)

      if (.not. allocated(this%pauli_large) .or. l < 0 .or. l > this%lmax .or. ispin < 1 .or. ispin > 2 .or. &
          any([size(large), size(small), size(large_dot), size(small_dot), size(gfac)] /= size(this%r))) then
         error stop 'radial_ground_state%set_basis_channel: inconsistent basis channel'
      end if
      this%pauli_large(:, l + 1, ispin) = large
      this%pauli_small(:, l + 1, ispin) = small
      this%sr_large_dot(:, l + 1, ispin) = large_dot
      this%sr_small_dot(:, l + 1, ispin) = small_dot
      this%sr_gfac(:, l + 1, ispin) = gfac
      this%sr_enu(l + 1, ispin) = enu
   end subroutine set_basis_channel

   subroutine finalize_radial_ground_state(this)
      type(radial_ground_state) :: this

      call this%clear()
   end subroutine finalize_radial_ground_state

end module radial_ground_state_mod
