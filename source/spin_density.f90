!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: spin_density
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> The one rotating-frame density contract shared by both SCF solvers
!> (GBT completion blueprint WP7 / gate G7, §2.5 and §4).
!>
!> Before WP7 the two routes disagreed about what "the density" is. The real-
!> space route reconstructed spin information from the full on-site Green
!> function and projected it, per energy channel, onto whatever
!> `potential%mom` happened to hold; the k-space route built projected spin
!> channels from eigenvectors along the *same* implicit axis and then updated a
!> Cartesian moment separately, afterwards. Neither route stated an axis, so
!> neither could be compared with the other and both could feed inconsistent
!> spin channels into the radial SCF.
!>
!> This module defines the object both routes must fill. Per site `a` and
!> angular channel `l` it retains the COMPLETE rotating-frame spin matrix
!>
!>   rho_al = ( rho_uu  rho_ud )      equivalently  rho = (n*1 + m.sigma)/2
!>            ( rho_du  rho_dd )
!>
!> for the three energy moments the radial SCF consumes,
!>
!>   order 1: M0 = int rho(E) dE       (occupation)
!>   order 2: M1 = int E rho(E) dE     (band centre numerator)
!>   order 3: M2 = int E^2 rho(E) dE   (second moment numerator)
!>
!> so that the up/down radial channels are formed ONLY after accumulation, by
!> `project_radial`/`radial_band_moments`, against an axis the caller must
!> state explicitly (`set_axis`). There is no default projection axis and no
!> read of `potential%mom` anywhere in this module.
!>
!> Conventions (do-not-get-this-wrong):
!>   * rho_{s s'} = psi_s psi*_{s'}, i.e. the k-space producer accumulates
!>     rho_uu = |u|^2, rho_dd = |d|^2, rho_ud = u conjg(d), which reproduces
!>     m_x = 2 Re(conjg(u) d), m_y = 2 Im(conjg(u) d), m_z = |u|^2 - |d|^2 --
!>     the same signs the pre-WP7 spinor projections used.
!>   * From a retarded Green function the matching Hermitian construction is
!>     rho = (i/2pi)(G - G^dagger), whose Cartesian components are
!>     m_x = -Im(G_ud + G_du)/pi, m_y = -Re(G_ud - G_du)/pi,
!>     m_z = -Im(G_uu - G_dd)/pi, n = -Im(G_uu + G_dd)/pi -- again the signs
!>     the pre-WP7 real-space accumulation used.
!>   * The stored matrix is the ROTATING-frame density. Lab-frame moments are
!>     reconstructed only for output/comparison, by `lab_frame_moment`, using
!>     the gauge D(phi) = diag(e^{-i phi/2}, e^{+i phi/2}) with phi = q.R_na.
!>
!> Two SCF policies are supported and are deliberately different objects, not
!> two readings of one code path (`resolve_site_axis`):
!>   * `constrained_spiral`: the reference direction is imposed and fixed. Only
!>     charge and the LONGITUDINAL moment magnitude are mixed; the transverse
!>     density is reported as a constraint residual and torque.
!>   * `relaxed_reference`: the full rotating-frame Cartesian moment is mixed
!>     and the per-sublattice reference axis follows it, retaining the single-q
!>     ansatz. The transverse residual is zero by construction.
!------------------------------------------------------------------------------
module spin_density_mod
   use precision_mod, only: rp
   use gbt_structure_mod, only: gbt_frame_t, gbt_frame_from_phase, gbt_rotating_to_lab_vector
   use logger_mod, only: g_logger
   use string_mod, only: int2str, real2str
   implicit none
   private

   !> Energy-moment orders retained per (site, l): int rho, int E rho, int E^2 rho.
   integer, parameter, public :: sd_orders = 3

   integer, parameter, public :: sd_name_len = 24

   !> SCF density policies (see the module header).
   character(len=sd_name_len), parameter, public :: sd_constrained_spiral = 'constrained_spiral'
   character(len=sd_name_len), parameter, public :: sd_relaxed_reference = 'relaxed_reference'

   !> Producer tags -- which solver filled the object.
   character(len=sd_name_len), parameter, public :: sd_producer_none = 'none'
   character(len=sd_name_len), parameter, public :: sd_producer_rs = 'real_space_green'
   character(len=sd_name_len), parameter, public :: sd_producer_kspace = 'kspace_eigenvectors'

   type, public :: spin_density
      !> Number of sites carried by this object.
      integer :: nsite = 0
      !> Number of angular channels per site (l = 0 .. nl-1).
      integer :: nl = 0
      !> SCF policy governing how the axis and the mixed variables are chosen.
      character(len=sd_name_len) :: policy = sd_constrained_spiral
      !> Which producer filled the object.
      character(len=sd_name_len) :: producer = sd_producer_none
      !> rho(s, s', l, site, order) -- rotating-frame 2x2 spin density matrix.
      complex(rp), dimension(:, :, :, :, :), allocatable :: rho
      !> Explicit per-site projection axis (unit vector). No default meaning.
      real(rp), dimension(:, :), allocatable :: axis
      !> .true. once set_axis has been called for that site.
      logical, dimension(:), allocatable :: axis_set
   contains
      procedure :: restore_to_default
      procedure :: zero_density
      procedure :: set_axis
      procedure :: get_axis
      procedure :: accumulate_block
      procedure :: cartesian
      procedure :: site_cartesian
      procedure :: project_radial
      procedure :: radial_band_moments
      procedure :: electron_count
      procedure :: check_physicality
      procedure :: max_deviation
      procedure :: lab_frame_moment
      procedure :: resolve_site_axis
      final :: destructor
   end type spin_density

   interface spin_density
      procedure :: constructor
   end interface spin_density

   public :: sd_matrix_from_cartesian
   public :: sd_cartesian_from_matrix

contains

   !---------------------------------------------------------------------------
   !> @brief Allocate an empty rotating-frame density object.
   !> @param[in] nsite Number of sites.
   !> @param[in] nl    Number of angular channels per site (l = 0 .. nl-1).
   !> @return Zeroed object with every axis unset.
   !---------------------------------------------------------------------------
   function constructor(nsite, nl) result(obj)
      integer, intent(in) :: nsite, nl
      type(spin_density) :: obj

      if (nsite < 1 .or. nl < 1) then
         call g_logger%fatal('spin_density: nsite and nl must be positive (got nsite='// &
                             trim(int2str(nsite))//', nl='//trim(int2str(nl))//').', __FILE__, __LINE__)
      end if

      obj%nsite = nsite
      obj%nl = nl
      allocate (obj%rho(2, 2, nl, nsite, sd_orders))
      allocate (obj%axis(3, nsite))
      allocate (obj%axis_set(nsite))
      call obj%restore_to_default()
   end function constructor

   subroutine destructor(this)
      type(spin_density) :: this

      if (allocated(this%rho)) deallocate (this%rho)
      if (allocated(this%axis)) deallocate (this%axis)
      if (allocated(this%axis_set)) deallocate (this%axis_set)
   end subroutine destructor

   !> @brief Zero the density and forget every projection axis.
   subroutine restore_to_default(this)
      class(spin_density), intent(inout) :: this

      this%rho(:, :, :, :, :) = (0.0_rp, 0.0_rp)
      this%axis(:, :) = 0.0_rp
      this%axis_set(:) = .false.
      this%policy = sd_constrained_spiral
      this%producer = sd_producer_none
   end subroutine restore_to_default

   !> @brief Clear the accumulated density and every stated axis, keeping the
   !>        configured policy.
   !> @details What a producer calls before re-accumulating. The axes go with
   !>          the density on purpose: an axis stated for a previous density is
   !>          exactly the stale projection definition WP7 exists to remove.
   subroutine zero_density(this)
      class(spin_density), intent(inout) :: this

      this%rho(:, :, :, :, :) = (0.0_rp, 0.0_rp)
      this%axis(:, :) = 0.0_rp
      this%axis_set(:) = .false.
   end subroutine zero_density

   !---------------------------------------------------------------------------
   !> @brief State the explicit projection axis for one site.
   !> @details This is the ONLY way an axis enters the contract. A vanishing
   !>          direction is a caller bug, not something to paper over with a
   !>          silent z fallback, so it is fatal.
   !---------------------------------------------------------------------------
   subroutine set_axis(this, isite, axis)
      class(spin_density), intent(inout) :: this
      integer, intent(in) :: isite
      real(rp), dimension(3), intent(in) :: axis
      real(rp) :: anorm

      call check_site(this, isite, 'set_axis')
      anorm = sqrt(sum(axis(:)**2))
      if (anorm <= 1.0e-12_rp) then
         call g_logger%fatal('spin_density%set_axis: zero-length projection axis for site '// &
                             trim(int2str(isite))//'; the caller must state a direction.', __FILE__, __LINE__)
      end if
      this%axis(:, isite) = axis(:)/anorm
      this%axis_set(isite) = .true.
   end subroutine set_axis

   !> @brief Return the explicit axis of one site; fatal if it was never set.
   subroutine get_axis(this, isite, axis)
      class(spin_density), intent(in) :: this
      integer, intent(in) :: isite
      real(rp), dimension(3), intent(out) :: axis

      call check_site(this, isite, 'get_axis')
      if (.not. this%axis_set(isite)) then
         call g_logger%fatal('spin_density%get_axis: no explicit axis was set for site '// &
                             trim(int2str(isite))//'; radial projection has no defined meaning.', &
                             __FILE__, __LINE__)
      end if
      axis(:) = this%axis(:, isite)
   end subroutine get_axis

   !---------------------------------------------------------------------------
   !> @brief Add one 2x2 contribution to (site, l, order).
   !> @details Producers call this; nothing else writes `rho`.
   !---------------------------------------------------------------------------
   subroutine accumulate_block(this, isite, il, iorder, block)
      class(spin_density), intent(inout) :: this
      integer, intent(in) :: isite, il, iorder
      complex(rp), dimension(2, 2), intent(in) :: block

      call check_index(this, isite, il, iorder, 'accumulate_block')
      this%rho(:, :, il, isite, iorder) = this%rho(:, :, il, isite, iorder) + block(:, :)
   end subroutine accumulate_block

   !---------------------------------------------------------------------------
   !> @brief Cartesian (n, m) form of one stored spin matrix.
   !> @param[out] n Trace (charge-like component).
   !> @param[out] m Cartesian moment, m_i = Tr(sigma_i rho).
   !---------------------------------------------------------------------------
   subroutine cartesian(this, isite, il, iorder, n, m)
      class(spin_density), intent(in) :: this
      integer, intent(in) :: isite, il, iorder
      real(rp), intent(out) :: n
      real(rp), dimension(3), intent(out) :: m

      call check_index(this, isite, il, iorder, 'cartesian')
      call sd_cartesian_from_matrix(this%rho(:, :, il, isite, iorder), n, m)
   end subroutine cartesian

   !> @brief Site total (n, m) of one energy-moment order, summed over l.
   subroutine site_cartesian(this, isite, iorder, n, m)
      class(spin_density), intent(in) :: this
      integer, intent(in) :: isite, iorder
      real(rp), intent(out) :: n
      real(rp), dimension(3), intent(out) :: m
      integer :: il
      real(rp) :: nl_val, ml_val(3)

      call check_index(this, isite, 1, iorder, 'site_cartesian')
      n = 0.0_rp
      m(:) = 0.0_rp
      do il = 1, this%nl
         call sd_cartesian_from_matrix(this%rho(:, :, il, isite, iorder), nl_val, ml_val)
         n = n + nl_val
         m(:) = m(:) + ml_val(:)
      end do
   end subroutine site_cartesian

   !---------------------------------------------------------------------------
   !> @brief Project one accumulated channel onto radial up/down.
   !> @details The projection happens HERE, after accumulation, against the
   !>          explicit site axis -- never per energy point against a moment
   !>          array that may since have been mixed or rotated.
   !---------------------------------------------------------------------------
   subroutine project_radial(this, isite, il, iorder, up, down)
      class(spin_density), intent(in) :: this
      integer, intent(in) :: isite, il, iorder
      real(rp), intent(out) :: up, down
      real(rp) :: n, m(3), axis(3)

      call this%get_axis(isite, axis)
      call this%cartesian(isite, il, iorder, n, m)
      up = 0.5_rp*(n + dot_product(m, axis))
      down = 0.5_rp*(n - dot_product(m, axis))
   end subroutine project_radial

   !---------------------------------------------------------------------------
   !> @brief The three radial band moments the SCF consumes for one (site,l,spin).
   !> @details q0 = occupation, centre = <E>, q2 = int (E-<E>)^2 rho dE. These
   !>          are exactly `potential%ql(1,l,isp)`, `gravity_center` (before the
   !>          vmad shift) and `potential%ql(3,l,isp)`.
   !> @param[in] ispin 1 = up (along the axis), 2 = down.
   !---------------------------------------------------------------------------
   subroutine radial_band_moments(this, isite, il, ispin, q0, centre, q2)
      class(spin_density), intent(in) :: this
      integer, intent(in) :: isite, il, ispin
      real(rp), intent(out) :: q0, centre, q2
      real(rp) :: up, down, p(sd_orders)
      integer :: iorder

      if (ispin /= 1 .and. ispin /= 2) then
         call g_logger%fatal('spin_density%radial_band_moments: ispin must be 1 or 2 (got '// &
                             trim(int2str(ispin))//').', __FILE__, __LINE__)
      end if
      do iorder = 1, sd_orders
         call this%project_radial(isite, il, iorder, up, down)
         if (ispin == 1) then
            p(iorder) = up
         else
            p(iorder) = down
         end if
      end do

      q0 = p(1)
      if (abs(q0) > epsilon(1.0_rp)) then
         centre = p(2)/p(1)
         q2 = p(3) - 2.0_rp*centre*p(2) + centre*centre*p(1)
      else
         q0 = 0.0_rp
         centre = 0.0_rp
         q2 = 0.0_rp
      end if
   end subroutine radial_band_moments

   !> @brief Total electron count carried by the object (order 1, all sites/l).
   real(rp) function electron_count(this)
      class(spin_density), intent(in) :: this
      integer :: isite
      real(rp) :: n, m(3)

      electron_count = 0.0_rp
      do isite = 1, this%nsite
         call this%site_cartesian(isite, 1, n, m)
         electron_count = electron_count + n
      end do
   end function electron_count

   !---------------------------------------------------------------------------
   !> @brief Density physicality assertions (blueprint WP7 step 8).
   !> @details For the occupation order, every (site, l) block must be
   !>          Hermitian, have a non-negative trace, have non-negative
   !>          eigenvalues (n +- |m|)/2, and satisfy |m| <= n. Optionally the
   !>          total electron count is compared against an expected value.
   !> @param[in]  tol       Absolute tolerance.
   !> @param[in]  expected_electrons Optional expected electron count.
   !> @param[out] ok        .true. if every check passed.
   !> @param[out] message   Empty when ok, else the first violation found.
   !---------------------------------------------------------------------------
   subroutine check_physicality(this, tol, ok, message, expected_electrons)
      class(spin_density), intent(in) :: this
      real(rp), intent(in) :: tol
      logical, intent(out) :: ok
      character(len=*), intent(out) :: message
      real(rp), intent(in), optional :: expected_electrons

      integer :: isite, il
      real(rp) :: n, m(3), mabs, eig_lo, herm, total
      complex(rp) :: blk(2, 2)
      character(len=64) :: where_str

      ok = .true.
      message = ''

      do isite = 1, this%nsite
         do il = 1, this%nl
            blk(:, :) = this%rho(:, :, il, isite, 1)
            where_str = 'site '//trim(int2str(isite))//' l '//trim(int2str(il - 1))

            herm = max(abs(blk(1, 2) - conjg(blk(2, 1))), &
                       abs(aimag(blk(1, 1))), abs(aimag(blk(2, 2))))
            if (herm > tol) then
               ok = .false.
               message = 'non-Hermitian density at '//trim(where_str)// &
                         ' (residual '//trim(real2str(herm, '(ES12.5)'))//')'
               return
            end if

            call sd_cartesian_from_matrix(blk, n, m)
            if (n < -tol) then
               ok = .false.
               message = 'negative trace at '//trim(where_str)// &
                         ' (n = '//trim(real2str(n, '(ES12.5)'))//')'
               return
            end if

            mabs = sqrt(sum(m(:)**2))
            eig_lo = 0.5_rp*(n - mabs)
            if (eig_lo < -tol) then
               ok = .false.
               message = 'negative density eigenvalue at '//trim(where_str)// &
                         ' (lambda_min = '//trim(real2str(eig_lo, '(ES12.5)'))//')'
               return
            end if
            if (mabs > n + tol) then
               ok = .false.
               message = '|m| > n at '//trim(where_str)// &
                         ' (|m| = '//trim(real2str(mabs, '(ES12.5)'))// &
                         ', n = '//trim(real2str(n, '(ES12.5)'))//')'
               return
            end if
         end do
      end do

      if (present(expected_electrons)) then
         total = this%electron_count()
         if (abs(total - expected_electrons) > tol) then
            ok = .false.
            message = 'electron count '//trim(real2str(total, '(F14.8)'))// &
                      ' differs from expected '//trim(real2str(expected_electrons, '(F14.8)'))
            return
         end if
      end if
   end subroutine check_physicality

   !---------------------------------------------------------------------------
   !> @brief Largest absolute difference between two filled contracts.
   !> @details The producer-equivalence measure: RS and k-space must fill the
   !>          same object, so this is what a cross-solver test asserts on.
   !---------------------------------------------------------------------------
   real(rp) function max_deviation(this, other)
      class(spin_density), intent(in) :: this
      type(spin_density), intent(in) :: other

      if (this%nsite /= other%nsite .or. this%nl /= other%nl) then
         call g_logger%fatal('spin_density%max_deviation: shape mismatch ('// &
                             trim(int2str(this%nsite))//'x'//trim(int2str(this%nl))//' vs '// &
                             trim(int2str(other%nsite))//'x'//trim(int2str(other%nl))//').', &
                             __FILE__, __LINE__)
      end if
      max_deviation = maxval(abs(this%rho - other%rho))
   end function max_deviation

   !---------------------------------------------------------------------------
   !> @brief Lab-frame moment of one channel, for output/comparison only.
   !> @details rho_lab = D(phi) rho D^dagger(phi) with
   !>          D(phi) = diag(e^{-i phi/2}, e^{+i phi/2}) and phi = q.R_na, so
   !>          the transverse moment rotates by +phi about z and the
   !>          longitudinal component is untouched. The SCF must never consume
   !>          this; it exists so a single-q run can be compared with an
   !>          explicit supercell.
   !---------------------------------------------------------------------------
   subroutine lab_frame_moment(this, isite, il, iorder, phase, n, m_lab)
      class(spin_density), intent(in) :: this
      integer, intent(in) :: isite, il, iorder
      real(rp), intent(in) :: phase
      real(rp), intent(out) :: n
      real(rp), dimension(3), intent(out) :: m_lab
      real(rp) :: m(3)
      type(gbt_frame_t) :: frame

      call this%cartesian(isite, il, iorder, n, m)
      ! A density stored by the SCF contract is already in the rotating
      ! reference coordinates.  For output/comparison, use the same
      ! authoritative frame map as the Hamiltonian rather than duplicating a
      ! z-rotation here.
      call gbt_frame_from_phase(0.0_rp, 0.0_rp, phase, frame)
      call gbt_rotating_to_lab_vector(frame, m, m_lab)
   end subroutine lab_frame_moment

   !---------------------------------------------------------------------------
   !> @brief Apply the SCF policy to one site and return what may be mixed.
   !> @details The two policies are genuinely different reductions of the same
   !>          accumulated density:
   !>
   !>   constrained_spiral -- `reference` is imposed and returned unchanged as
   !>     the projection axis. `m_long` is the longitudinal magnitude (the only
   !>     magnetic variable that may be mixed), `m_transverse` is the transverse
   !>     residual and `torque` = reference x m is the constraining torque that
   !>     must be reported.
   !>
   !>   relaxed_reference -- the axis follows the full rotating-frame Cartesian
   !>     moment, `m_long` is |m|, and the transverse residual and torque vanish
   !>     identically. The single-q ansatz is untouched: only the per-sublattice
   !>     reference direction moves.
   !>
   !> @param[in]  reference Imposed reference direction (used by both policies;
   !>                       the relaxed policy only falls back to it when the
   !>                       accumulated moment is numerically zero).
   !> @param[out] axis_out  Axis the caller must hand to `set_axis`.
   !---------------------------------------------------------------------------
   subroutine resolve_site_axis(this, isite, reference, axis_out, charge, m_long, m_transverse, torque)
      class(spin_density), intent(in) :: this
      integer, intent(in) :: isite
      real(rp), dimension(3), intent(in) :: reference
      real(rp), dimension(3), intent(out) :: axis_out
      real(rp), intent(out) :: charge
      real(rp), intent(out) :: m_long
      real(rp), dimension(3), intent(out) :: m_transverse
      real(rp), dimension(3), intent(out) :: torque
      real(rp) :: m(3), ref(3), refnorm, mabs

      call check_site(this, isite, 'resolve_site_axis')
      refnorm = sqrt(sum(reference(:)**2))
      if (refnorm <= 1.0e-12_rp) then
         call g_logger%fatal('spin_density%resolve_site_axis: zero-length reference direction '// &
                             'for site '//trim(int2str(isite))//'.', __FILE__, __LINE__)
      end if
      ref(:) = reference(:)/refnorm

      call this%site_cartesian(isite, 1, charge, m)
      mabs = sqrt(sum(m(:)**2))

      select case (trim(this%policy))
      case (sd_constrained_spiral)
         axis_out(:) = ref(:)
         m_long = dot_product(m, ref)
         m_transverse(:) = m(:) - m_long*ref(:)
         torque(:) = cross(ref, m)
      case (sd_relaxed_reference)
         if (mabs > 1.0e-12_rp) then
            axis_out(:) = m(:)/mabs
         else
            axis_out(:) = ref(:)
         end if
         m_long = mabs
         m_transverse(:) = 0.0_rp
         torque(:) = 0.0_rp
      case default
         call g_logger%fatal('spin_density%resolve_site_axis: unknown policy "'// &
                             trim(this%policy)//'"; expected "'//trim(sd_constrained_spiral)// &
                             '" or "'//trim(sd_relaxed_reference)//'".', __FILE__, __LINE__)
      end select
   end subroutine resolve_site_axis

   !---------------------------------------------------------------------------
   !> @brief Build the 2x2 spin matrix from (n, m): rho = (n*1 + m.sigma)/2.
   !---------------------------------------------------------------------------
   pure subroutine sd_matrix_from_cartesian(n, m, block)
      real(rp), intent(in) :: n
      real(rp), dimension(3), intent(in) :: m
      complex(rp), dimension(2, 2), intent(out) :: block

      block(1, 1) = cmplx(0.5_rp*(n + m(3)), 0.0_rp, rp)
      block(2, 2) = cmplx(0.5_rp*(n - m(3)), 0.0_rp, rp)
      block(1, 2) = cmplx(0.5_rp*m(1), -0.5_rp*m(2), rp)
      block(2, 1) = conjg(block(1, 2))
   end subroutine sd_matrix_from_cartesian

   !---------------------------------------------------------------------------
   !> @brief Cartesian components of a 2x2 spin matrix: n = Tr rho,
   !>        m_i = Tr(sigma_i rho).
   !---------------------------------------------------------------------------
   pure subroutine sd_cartesian_from_matrix(block, n, m)
      complex(rp), dimension(2, 2), intent(in) :: block
      real(rp), intent(out) :: n
      real(rp), dimension(3), intent(out) :: m

      n = real(block(1, 1) + block(2, 2), rp)
      m(1) = real(block(1, 2) + block(2, 1), rp)
      m(2) = -aimag(block(1, 2) - block(2, 1))
      m(3) = real(block(1, 1) - block(2, 2), rp)
   end subroutine sd_cartesian_from_matrix

   pure function cross(a, b) result(c)
      real(rp), dimension(3), intent(in) :: a, b
      real(rp), dimension(3) :: c

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   end function cross

   subroutine check_site(this, isite, caller)
      class(spin_density), intent(in) :: this
      integer, intent(in) :: isite
      character(len=*), intent(in) :: caller

      if (isite < 1 .or. isite > this%nsite) then
         call g_logger%fatal('spin_density%'//caller//': site index '//trim(int2str(isite))// &
                             ' outside 1..'//trim(int2str(this%nsite))//'.', __FILE__, __LINE__)
      end if
   end subroutine check_site

   subroutine check_index(this, isite, il, iorder, caller)
      class(spin_density), intent(in) :: this
      integer, intent(in) :: isite, il, iorder
      character(len=*), intent(in) :: caller

      call check_site(this, isite, caller)
      if (il < 1 .or. il > this%nl) then
         call g_logger%fatal('spin_density%'//caller//': angular channel '//trim(int2str(il))// &
                             ' outside 1..'//trim(int2str(this%nl))//'.', __FILE__, __LINE__)
      end if
      if (iorder < 1 .or. iorder > sd_orders) then
         call g_logger%fatal('spin_density%'//caller//': energy-moment order '//trim(int2str(iorder))// &
                             ' outside 1..'//trim(int2str(sd_orders))//'.', __FILE__, __LINE__)
      end if
   end subroutine check_index

end module spin_density_mod
