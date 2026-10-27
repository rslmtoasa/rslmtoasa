!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Mix
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
!> lorem ipsum
!------------------------------------------------------------------------------

module mix_mod
   use control_mod
   use lattice_mod
   use charge_mod
   use symbolic_atom_mod, only: symbolic_atom
   use precision_mod, only: rp
   use mpi_mod
   use string_mod
   use logger_mod, only: g_logger
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   use basis_mod, only: lmax_basis
   use math_mod, only: hubbard_density_diagnostics
   use scf_guard_mod, only: guard_q0_capacity, guard_nonfinite_scalar
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_MPI
   use mpi, only: MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_MAX, MPI_MIN, MPI_IN_PLACE, MPI_ALLREDUCE
#endif
   implicit none

   private

   ! Frozen pre-campaign convergence targets: dimensionless full-block density
   ! norm and Ry potential norm. These are goals, not claimed integration-noise
   ! floors; production discretization refinement must verify attainability.
   real(rp), parameter, public :: fixed_u_density_tolerance = 1.0e-7_rp
   real(rp), parameter, public :: fixed_u_potential_tolerance = 1.0e-7_rp
   real(rp), parameter, public :: fixed_u_test_tight_density_tolerance = 1.0e-8_rp
   real(rp), parameter, public :: fixed_u_test_tight_potential_tolerance = 1.0e-8_rp
   integer, parameter, public :: fixed_u_required_stable_updates = 3
   public :: fixed_u_residual_metrics

   !> Module´s main structure
   type, public :: mix
      ! Lattice
      class(lattice), pointer :: lattice
      ! Charge
      class(charge), pointer :: charge
      ! Symbolic atom
      class(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Description
      !>
      !> Description
      integer :: var
      !> Variables to save parameters for mixing
      real(rp), dimension(:, :), allocatable :: qia, qia_new, qia_old, qiaprev
      !> Variables to save magnetic moments for mixing
      real(rp), dimension(:, :), allocatable :: mag_old, mag_new, mag_mix
      !> Variables to save LDA+U density matrices for optional linear mixing
      !> Dimensions: (atom, l-channel, spin, flattened m,m’)
      complex(rp), dimension(:, :, :, :), allocatable :: ldm_old, ldm_new, ldm_mix
      !> Mixing parameter
      real(rp) :: beta
      !> Linear mixing factor for LDA+U density matrix (ldm).
      !> 1.0 => no ldm mixing, 0.0 => keep previous ldm.
      real(rp) :: ldm_beta
      !> Magnetic mixing parameter
      real(rp), dimension(:), allocatable :: magbeta
      !> Difference between interactions
      real(rp) :: delta
      !> Type of mixing. Can be linear or broyden
      character(len=7) :: mixtype
      !> Variables for Broyden mixing
      real(rp), dimension(:), allocatable :: v_broy, u_broy, fo_broy, muo_broy
      real(rp) :: fsqo
      integer :: itr, nmix
      ! Number of consecutive fixed-U iterations meeting the full complex P/V
      ! residual targets. Invalid matrices reset this count.
      integer :: fixed_u_stable_updates = 0
      !> Variable to managed induced magnetic moments
      logical, dimension(:), allocatable :: is_induced
   contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: print_state
      procedure :: save_to
      procedure :: mixpq
      procedure :: save_ldm_to
      procedure :: mix_ldm_linear
      procedure :: advance_fixed_u_convergence
      procedure :: mix_magnetic_moments
      final :: destructor
   end type mix

   interface mix
      procedure :: constructor
   end interface mix

contains

   !> Compute full-complex fixed-U residuals from explicit, active-shell
   !> snapshots. Inputs are packed as (site,spin,row,column); the Frobenius
   !> norm combines both collinear spin blocks for each physical site.
   subroutine fixed_u_residual_metrics(p_in, p_out, v_in, v_out, rp_abs, rp_rel, rv_abs, rv_rel, &
                                       physical, antiherm_p, antiherm_v, worst_site)
      complex(rp), intent(in) :: p_in(:, :, :, :), p_out(:, :, :, :)
      complex(rp), intent(in) :: v_in(:, :, :, :), v_out(:, :, :, :)
      real(rp), intent(out) :: rp_abs, rp_rel, rv_abs, rv_rel
      logical, intent(out) :: physical
      real(rp), intent(out) :: antiherm_p, antiherm_v
      integer, intent(out) :: worst_site
      real(rp), allocatable :: rp_site(:), rp_rel_site(:), rv_site(:), rv_rel_site(:)
      real(rp) :: tr, emin, emax, ah, pin_norm, pout_norm, vin_norm, vout_norm
      integer :: nsite, nspin, norb, site, spin, info, local_physical
#ifdef USE_MPI
      integer :: mpierr
#endif

      nsite = size(p_in, 1)
      nspin = size(p_in, 2)
      norb = size(p_in, 3)
      rp_abs = huge(1.0_rp); rp_rel = huge(1.0_rp)
      rv_abs = huge(1.0_rp); rv_rel = huge(1.0_rp)
      antiherm_p = huge(1.0_rp); antiherm_v = huge(1.0_rp)
      physical = .false.; worst_site = 0
      ! nspin=2: collinear spin blocks; nspin=1: one full spinor block per site (nsp=2 spinor U)
      if (nsite < 1 .or. nspin < 1 .or. nspin > 2 .or. norb < 1) return
      if (size(p_out,1) /= nsite .or. size(v_in,1) /= nsite .or. size(v_out,1) /= nsite .or. &
          size(p_out,2) /= nspin .or. size(v_in,2) /= nspin .or. size(v_out,2) /= nspin .or. &
          size(p_out,3) /= norb .or. size(p_out,4) /= norb .or. &
          size(p_in,4) /= norb .or. size(v_in,3) /= norb .or. size(v_in,4) /= norb .or. &
          size(v_out,3) /= norb .or. size(v_out,4) /= norb) return
      if (.not. all(ieee_is_finite(real(p_in,rp))) .or. .not. all(ieee_is_finite(aimag(p_in))) .or. &
          .not. all(ieee_is_finite(real(p_out,rp))) .or. .not. all(ieee_is_finite(aimag(p_out))) .or. &
          .not. all(ieee_is_finite(real(v_in,rp))) .or. .not. all(ieee_is_finite(aimag(v_in))) .or. &
          .not. all(ieee_is_finite(real(v_out,rp))) .or. .not. all(ieee_is_finite(aimag(v_out)))) return

      allocate(rp_site(nsite), rp_rel_site(nsite), rv_site(nsite), rv_rel_site(nsite))
      local_physical = 1
      antiherm_p = 0.0_rp; antiherm_v = 0.0_rp
      do site = 1, nsite
         pin_norm = 0.0_rp; pout_norm = 0.0_rp; vin_norm = 0.0_rp; vout_norm = 0.0_rp
         rp_site(site) = 0.0_rp; rv_site(site) = 0.0_rp
         do spin = 1, nspin
            pin_norm = pin_norm + sum(abs(p_in(site,spin,:,:))**2)
            pout_norm = pout_norm + sum(abs(p_out(site,spin,:,:))**2)
            vin_norm = vin_norm + sum(abs(v_in(site,spin,:,:))**2)
            vout_norm = vout_norm + sum(abs(v_out(site,spin,:,:))**2)
            rp_site(site) = rp_site(site) + sum(abs(p_out(site,spin,:,:)-p_in(site,spin,:,:))**2)
            rv_site(site) = rv_site(site) + sum(abs(v_out(site,spin,:,:)-v_in(site,spin,:,:))**2)
            call hubbard_density_diagnostics(p_in(site,spin,:,:), tr, emin, emax, ah, info)
            antiherm_p = max(antiherm_p, ah)
            ! The 1e-6 anti-Hermitian and 1e-3 occupation-bound allowances are
            ! diagnostic tolerances inherited from the fixed-U occupation
            ! checks; no clipping/projection is performed.
            if (info /= 0 .or. ah > 1.0e-6_rp .or. emin < -1.0e-3_rp .or. emax > 1.0_rp + 1.0e-3_rp .or. &
                tr < -1.0e-3_rp .or. tr > real(norb,rp) + 1.0e-3_rp) local_physical = 0
            call hubbard_density_diagnostics(p_out(site,spin,:,:), tr, emin, emax, ah, info)
            antiherm_p = max(antiherm_p, ah)
            if (info /= 0 .or. ah > 1.0e-6_rp .or. emin < -1.0e-3_rp .or. emax > 1.0_rp + 1.0e-3_rp .or. &
                tr < -1.0e-3_rp .or. tr > real(norb,rp) + 1.0e-3_rp) local_physical = 0
            ah = maxval(abs(v_in(site,spin,:,:)-transpose(conjg(v_in(site,spin,:,:)))))
            antiherm_v = max(antiherm_v, ah)
            if (ah > 1.0e-6_rp) local_physical = 0
            ah = maxval(abs(v_out(site,spin,:,:)-transpose(conjg(v_out(site,spin,:,:)))))
            antiherm_v = max(antiherm_v, ah)
            if (ah > 1.0e-6_rp) local_physical = 0
         end do
         rp_site(site) = sqrt(rp_site(site))
         rv_site(site) = sqrt(rv_site(site))
         rp_rel_site(site) = rp_site(site)/max(1.0_rp,sqrt(pin_norm),sqrt(pout_norm))
         rv_rel_site(site) = rv_site(site)/max(1.0_rp,sqrt(vin_norm),sqrt(vout_norm))
      end do
#ifdef USE_MPI
      if (numprocs > 1) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, rp_site, nsite, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, rp_rel_site, nsite, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, rv_site, nsite, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, rv_rel_site, nsite, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, antiherm_p, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, antiherm_v, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, local_physical, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, mpierr)
      end if
#endif
      rp_abs = maxval(rp_site); rp_rel = maxval(rp_rel_site)
      rv_abs = maxval(rv_site); rv_rel = maxval(rv_rel_site)
      worst_site = maxloc(max(rp_rel_site,rv_rel_site), dim=1)
      physical = local_physical == 1
   end subroutine fixed_u_residual_metrics

   !> Advance a three-success stability gate. A bad physicality check or a
   !> residual above target clears the streak; callers skip this for inactive U.
   subroutine advance_fixed_u_convergence(this, rp_residual, rv_residual, physical, converged, &
                                          density_tolerance, potential_tolerance)
      class(mix), intent(inout) :: this
      real(rp), intent(in) :: rp_residual, rv_residual
      logical, intent(in) :: physical
      logical, intent(out) :: converged
      real(rp), intent(in), optional :: density_tolerance, potential_tolerance
      real(rp) :: rp_target, rv_target
      rp_target = fixed_u_density_tolerance
      rv_target = fixed_u_potential_tolerance
      if (present(density_tolerance)) rp_target = density_tolerance
      if (present(potential_tolerance)) rv_target = potential_tolerance
      if (.not. ieee_is_finite(rp_target) .or. .not. ieee_is_finite(rv_target) .or. &
          rp_target <= 0.0_rp .or. rv_target <= 0.0_rp) then
         call g_logger%fatal('Fixed-U convergence thresholds must be positive and finite.', __FILE__, __LINE__)
      end if
      if (.not. physical .or. .not. ieee_is_finite(rp_residual) .or. .not. ieee_is_finite(rv_residual)) then
         this%fixed_u_stable_updates = 0
      else if (rp_residual <= rp_target .and. rv_residual <= rv_target) then
         this%fixed_u_stable_updates = this%fixed_u_stable_updates + 1
      else
         this%fixed_u_stable_updates = 0
      end if
      converged = this%fixed_u_stable_updates >= fixed_u_required_stable_updates
   end subroutine advance_fixed_u_convergence

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] lattice_obj pointer
   !> @return type(mix)
   !---------------------------------------------------------------------------
   function constructor(lattice_obj, charge_obj) result(obj)
      type(mix) :: obj
      type(lattice), target, intent(in) :: lattice_obj
      type(charge), target, intent(in) :: charge_obj

      obj%lattice => lattice_obj
      obj%charge => charge_obj
      obj%symbolic_atom => lattice_obj%symbolic_atoms

      call obj%restore_to_default()
      call obj%build_from_file()
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(mix) :: this

#ifdef USE_SAFE_ALLOC
      if (allocated(this%qia)) call g_safe_alloc%deallocate('mix.qia', this%qia)
      if (allocated(this%qia_new)) call g_safe_alloc%deallocate('mix.qia_new', this%qia_new)
      if (allocated(this%qia_old)) call g_safe_alloc%deallocate('mix.qia_old', this%qia_old)
      if (allocated(this%qiaprev)) call g_safe_alloc%deallocate('mix.qiaprev', this%qiaprev)
      if (allocated(this%v_broy)) call g_safe_alloc%deallocate('mix.v_broy', this%v_broy)
      if (allocated(this%u_broy)) call g_safe_alloc%deallocate('mix.u_broy', this%u_broy)
      if (allocated(this%fo_broy)) call g_safe_alloc%deallocate('mix.fo_broy', this%fo_broy)
      if (allocated(this%muo_broy)) call g_safe_alloc%deallocate('mix.muo_broy', this%muo_broy)
      if (allocated(this%magbeta)) call g_safe_alloc%deallocate('mix.magbeta', this%magbeta)
      if (allocated(this%mag_new)) call g_safe_alloc%deallocate('mix.mag_new', this%mag_new)
      if (allocated(this%mag_old)) call g_safe_alloc%deallocate('mix.mag_old', this%mag_old)
      if (allocated(this%mag_mix)) call g_safe_alloc%deallocate('mix.mag_mix', this%mag_mix)
      if (allocated(this%ldm_old)) call g_safe_alloc%deallocate('mix.ldm_old', this%ldm_old)
      if (allocated(this%ldm_new)) call g_safe_alloc%deallocate('mix.ldm_new', this%ldm_new)
      if (allocated(this%ldm_mix)) call g_safe_alloc%deallocate('mix.ldm_mix', this%ldm_mix)
      if (allocated(this%is_induced)) call g_safe_alloc%deallocate('mix.is_induced', this%is_induced)
#else
      if (allocated(this%qia)) deallocate (this%qia)
      if (allocated(this%qia_new)) deallocate (this%qia_new)
      if (allocated(this%qia_old)) deallocate (this%qia_old)
      if (allocated(this%qiaprev)) deallocate (this%qiaprev)
      if (allocated(this%v_broy)) deallocate (this%v_broy)
      if (allocated(this%u_broy)) deallocate (this%u_broy)
      if (allocated(this%fo_broy)) deallocate (this%fo_broy)
      if (allocated(this%muo_broy)) deallocate (this%muo_broy)
      if (allocated(this%magbeta)) deallocate (this%magbeta)
      if (allocated(this%mag_new)) deallocate (this%mag_new)
      if (allocated(this%mag_old)) deallocate (this%mag_old)
      if (allocated(this%mag_mix)) deallocate (this%mag_mix)
      if (allocated(this%ldm_old)) deallocate (this%ldm_old)
      if (allocated(this%ldm_new)) deallocate (this%ldm_new)
      if (allocated(this%ldm_mix)) deallocate (this%ldm_mix)
      if (allocated(this%is_induced)) deallocate (this%is_induced)
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
   subroutine build_from_file(this)
      class(mix), intent(inout) :: this
      character(len=sl) :: fname
      ! variables associated with the reading processes
      integer :: iostatus, funit

      include 'include_codes/namelists/mix.f90'

      fname = this%lattice%control%fname
      var = this%var
      mixtype = this%mixtype
      beta = this%beta
      ldm_beta = this%ldm_beta
      if (allocated(this%magbeta)) then
         allocate(magbeta(size(this%magbeta)))
         magbeta = this%magbeta
      end if

      open (newunit=funit, file=fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//fmt('A', fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=mix, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus) .and. iostatus /= 5010) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//fmt('I0', iostatus), __FILE__, __LINE__)
      else if (iostatus == 5010) then
         call g_logger%warning('Non-fatal namelist read value mismatch in '//trim(fname)//' (iostatus=5010).', __FILE__, __LINE__)
      end if
      close(funit)
      ! PATCH_MIX_NML_FIX:
      ! Persist parsed namelist values back to object state.
      ! Without this, this%mixtype keeps restore_to_default() value (’linear’).
      this%var = var
      this%mixtype = trim(mixtype)
      this%beta = beta
      this%ldm_beta = min(1.0_rp, max(0.0_rp, ldm_beta))
      if (allocated(magbeta)) then
         call move_alloc(magbeta, this%magbeta)
      end if
   end subroutine build_from_file

      !--------------------------------------------------------------------------
      ! DESCRIPTION:
      !> @brief
      !> Restore defaults and allocate mixing arrays
      subroutine restore_to_default(this)
      class(mix), intent(inout) :: this
      integer :: qia_width
      integer :: ldm_flat

      qia_width = 6*(lmax_basis + 1)
      ldm_flat = (2*lmax_basis + 1)*(2*lmax_basis + 1)
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('mix.qia', this%qia, (/this%lattice%nrec, qia_width/))
      call g_safe_alloc%allocate('mix.qia_new', this%qia_new, (/this%lattice%nrec, qia_width/))
      call g_safe_alloc%allocate('mix.qia_old', this%qia_old, (/this%lattice%nrec, qia_width/))
      call g_safe_alloc%allocate('mix.qiaprev', this%qiaprev, (/this%lattice%nrec, qia_width/))
      call g_safe_alloc%allocate('mix.v_broy', this%v_broy, (/this%lattice%nrec*qia_width/))
      call g_safe_alloc%allocate('mix.u_broy', this%u_broy, (/this%lattice%nrec*qia_width/))
      call g_safe_alloc%allocate('mix.fo_broy', this%fo_broy, (/this%lattice%nrec*qia_width/))
      call g_safe_alloc%allocate('mix.muo_broy', this%muo_broy, (/this%lattice%nrec*qia_width/))
      call g_safe_alloc%allocate('mix.magbeta', this%magbeta, (/this%lattice%nrec/))
      call g_safe_alloc%allocate('mix.mag_old', this%mag_old, (/this%lattice%nrec, 3/))
      call g_safe_alloc%allocate('mix.mag_new', this%mag_new, (/this%lattice%nrec, 3/))
      call g_safe_alloc%allocate('mix.mag_mix', this%mag_mix, (/this%lattice%nrec, 3/))
      call g_safe_alloc%allocate('mix.ldm_old', this%ldm_old, (/this%lattice%nrec, lmax_basis + 1, 2, ldm_flat/))
      call g_safe_alloc%allocate('mix.ldm_new', this%ldm_new, (/this%lattice%nrec, lmax_basis + 1, 2, ldm_flat/))
      call g_safe_alloc%allocate('mix.ldm_mix', this%ldm_mix, (/this%lattice%nrec, lmax_basis + 1, 2, ldm_flat/))
      call g_safe_alloc%allocate('mix.is_induced', this%is_induced, (/this%lattice%nrec/))
#else
      allocate (this%qia(this%lattice%nrec, qia_width), this%qia_new(this%lattice%nrec, qia_width), this%qia_old(this%lattice%nrec, qia_width))
      allocate (this%qiaprev(this%lattice%nrec, qia_width))
      allocate (this%v_broy(this%lattice%nrec*qia_width), this%u_broy(this%lattice%nrec*qia_width), this%fo_broy(this%lattice%nrec*qia_width), this%muo_broy(this%lattice%nrec*qia_width))
      allocate (this%magbeta(this%lattice%nrec), this%mag_old(this%lattice%nrec, 3), this%mag_new(this%lattice%nrec, 3), this%mag_mix(this%lattice%nrec, 3))
      allocate (this%ldm_old(this%lattice%nrec, lmax_basis + 1, 2, ldm_flat))
      allocate (this%ldm_new(this%lattice%nrec, lmax_basis + 1, 2, ldm_flat))
      allocate (this%ldm_mix(this%lattice%nrec, lmax_basis + 1, 2, ldm_flat))
      allocate (this%is_induced(this%lattice%nrec))
#endif

      this%qia(:, :) = 0.0_rp
      this%qia_new(:, :) = 0.0_rp
      this%qia_old(:, :) = 0.0_rp
      this%qiaprev(:, :) = 0.0_rp
      this%ldm_old(:, :, :, :) = 0.0_rp
      this%ldm_new(:, :, :, :) = 0.0_rp
      this%ldm_mix(:, :, :, :) = 0.0_rp
      this%v_broy(:) = 0.0d0
      this%u_broy(:) = 0.0d0
      this%fo_broy(:) = 0.0d0
      this%muo_broy(:) = 0.0d0
      this%fsqo = 1.0d0
      this%itr = 0
      this%fixed_u_stable_updates = 0
      this%var = 0
      this%beta = 0.1d0
      this%ldm_beta = 1.0d0
      this%magbeta(:) = 1.0d0
      this%nmix = 2
      this%mixtype = 'linear'
      this%is_induced = .false.
      end subroutine restore_to_default

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
      class(mix), intent(in) :: this

      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: file
      integer :: newunit

      include 'include_codes/namelists/mix.f90'

      var = this%var
      beta = this%beta
      mixtype = this%mixtype
      ldm_beta = this%ldm_beta
      if (allocated(this%magbeta)) then
         if (allocated(magbeta)) deallocate(magbeta)
         allocate(magbeta(size(this%magbeta)))
         magbeta = this%magbeta
      end if

      if (present(unit) .and. present(file)) then
         call g_logger%fatal('Argument error: both unit and file are present', __FILE__, __LINE__)
      else if (present(unit)) then
         write (unit, nml=mix)
      else if (present(file)) then
         open (newunit=newunit, file=file)
         write (newunit, nml=mix)
         close(newunit)
      end if
   end subroutine print_state

   !--------------------------------------------------------------------------
   !> @brief
   !> Save or restore QL/PL parameters to/from `this%qia` arrays.
   subroutine save_to(this, whereto, iteration, guard_collinear)
      use mpi_mod
      class(mix), intent(inout) :: this
      character(len=*), intent(in) :: whereto
      integer, intent(in), optional :: iteration
      logical, intent(in), optional :: guard_collinear
      integer :: it, I, lcount, qcols, spin
      integer :: iteration_value
      logical :: guard_active

      iteration_value = 0
      guard_active = .false.
      if (present(iteration)) iteration_value = iteration
      if (present(guard_collinear)) guard_active = guard_collinear

      select case (trim(whereto))
      case ('new')
         this%qia_new = 0.0_rp
         do it = 1, this%lattice%nrec
            lcount = min(this%symbolic_atom(this%lattice%nbulk + it)%potential%lmax, lmax_basis) + 1
            qcols = size(this%qia_new, 2)
            do I = 1, lcount
               if (guard_active) then
                  do spin = 1, 2
                     call guard_nonfinite_scalar(this%symbolic_atom(this%lattice%nbulk + it)%element%symbol, &
                        it, I - 1, spin, 'new q1', iteration_value, &
                        this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(2, I - 1, spin))
                  end do
               end if
               if (I              <= qcols) this%qia_new(it, I             ) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, I - 1, 1)
               if (I +     lcount <= qcols) this%qia_new(it, I +     lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, I - 1, 2)
               if (I + 2 * lcount <= qcols) this%qia_new(it, I + 2 * lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, I - 1, 1)
               if (I + 3 * lcount <= qcols) this%qia_new(it, I + 3 * lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, I - 1, 2)
               if (I + 4 * lcount <= qcols) this%qia_new(it, I + 4 * lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(I - 1, 1)
               if (I + 5 * lcount <= qcols) this%qia_new(it, I + 5 * lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(I - 1, 2)
            end do
            if (rank == 0) then
                  if (any(this%qia_new(it, :) /= this%qia_new(it, :)) .or. maxval(abs(this%qia_new(it, :))) > 1.0e6_rp) then
                     call g_logger%warning('Suspicious values in qia_new for atom '//fmt('i4', it), __FILE__, __LINE__)
               end if
            end if
         end do
         if (guard_active) call validate_packed_q(this, this%qia_new, iteration_value, 'new')
      case ('old')
         this%qia_old = 0.0_rp
         do it = 1, this%lattice%nrec
            lcount = min(this%symbolic_atom(this%lattice%nbulk + it)%potential%lmax, lmax_basis) + 1
            qcols = size(this%qia_old, 2)
            do I = 1, lcount
               if (I              <= qcols) this%qia_old(it, I             ) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, I - 1, 1)
               if (I +     lcount <= qcols) this%qia_old(it, I +     lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, I - 1, 2)
               if (I + 2 * lcount <= qcols) this%qia_old(it, I + 2 * lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, I - 1, 1)
               if (I + 3 * lcount <= qcols) this%qia_old(it, I + 3 * lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, I - 1, 2)
               if (I + 4 * lcount <= qcols) this%qia_old(it, I + 4 * lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(I - 1, 1)
               if (I + 5 * lcount <= qcols) this%qia_old(it, I + 5 * lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(I - 1, 2)
            end do
            if (rank == 0) then
               if (any(this%qia_old(it, :) /= this%qia_old(it, :)) .or. maxval(abs(this%qia_old(it, :))) > 1.0e6_rp) then
                  call g_logger%warning('Suspicious values in qia_old for atom '//fmt('i4', it), __FILE__, __LINE__)
               end if
            end if
         end do
      case ('prev')
         this%qiaprev = 0.0_rp
         do it = 1, this%lattice%nrec
            lcount = min(this%symbolic_atom(this%lattice%nbulk + it)%potential%lmax, lmax_basis) + 1
            qcols = size(this%qiaprev, 2)
            do I = 1, lcount
               if (I              <= qcols) this%qiaprev(it, I             ) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, I - 1, 1)
               if (I +     lcount <= qcols) this%qiaprev(it, I +     lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, I - 1, 2)
               if (I + 2 * lcount <= qcols) this%qiaprev(it, I + 2 * lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, I - 1, 1)
               if (I + 3 * lcount <= qcols) this%qiaprev(it, I + 3 * lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, I - 1, 2)
               if (I + 4 * lcount <= qcols) this%qiaprev(it, I + 4 * lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(I - 1, 1)
               if (I + 5 * lcount <= qcols) this%qiaprev(it, I + 5 * lcount) = this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(I - 1, 2)
            end do
            if (rank == 0) then
               if (any(this%qiaprev(it, :) /= this%qiaprev(it, :)) .or. maxval(abs(this%qiaprev(it, :))) > 1.0e6_rp) then
                  call g_logger%warning('Suspicious values in qiaprev for atom '//fmt('i4', it), __FILE__, __LINE__)
               end if
            end if
         end do

      case ('current')
         do it = 1, this%lattice%nrec
            lcount = min(this%symbolic_atom(this%lattice%nbulk + it)%potential%lmax, lmax_basis) + 1
            qcols = size(this%qia, 2)
            do I = 1, lcount
               if (I              <= qcols) this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, I - 1, 1) = this%qia(it, I             )
               if (I +     lcount <= qcols) this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(1, I - 1, 2) = this%qia(it, I +     lcount)
               if (I + 2 * lcount <= qcols) this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, I - 1, 1) = this%qia(it, I + 2 * lcount)
               if (I + 3 * lcount <= qcols) this%symbolic_atom(this%lattice%nbulk + it)%potential%ql(3, I - 1, 2) = this%qia(it, I + 3 * lcount)
               if (I + 4 * lcount <= qcols) this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(I - 1, 1)    = this%qia(it, I + 4 * lcount)
               if (I + 5 * lcount <= qcols) this%symbolic_atom(this%lattice%nbulk + it)%potential%pl(I - 1, 2)    = this%qia(it, I + 5 * lcount)
            end do
         end do
      end select
   end subroutine save_to

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Mix the magnetic moments.
   !>
   !> Mix the magnetic moments.
   !> @param[in] magnetic moments read from the potential file
   !> @return type(symbolic_atom): new magnetic moments after mixing.
   !---------------------------------------------------------------------------
   subroutine mix_magnetic_moments(this, mag_old, mag_new, mag_mix, mtot)
      use mpi_mod
      implicit none
      class(mix), intent(inout) :: this
      ! Input variables
      real(rp), dimension(:, :), intent(in) :: mag_old, mag_new
      real(rp), dimension(:), intent(in) :: mtot
      ! Output variable
      real(rp), dimension(:, :), intent(out) :: mag_mix

      ! Local variable
      integer :: ia ! Self-consistent atom index

      ! Calculate the mixed magnetic moments
      do ia = 1, this%lattice%nrec 
         if (mtot(ia + this%lattice%nbulk) < 0.5d0) then
            this%is_induced(ia) = .true.
            !this%magbeta(ia) = 0.0d0
            if (rank == 0) call g_logger%info('Spin moment at atom '//fmt('i4', (ia))//' is being considered induced', __FILE__, __LINE__)
         end if
         ! Defensive checks on magbeta
         if (.not. allocated(this%magbeta)) then
            if (rank == 0) call g_logger%warning('magbeta not allocated — resetting to 1.0 for all atoms', __FILE__, __LINE__)
            allocate(this%magbeta(this%lattice%nrec))
            this%magbeta(:) = 1.0d0
         end if
         if (this%magbeta(ia) /= this%magbeta(ia) .or. this%magbeta(ia) < 0.0d0 .or. this%magbeta(ia) > 1.0e6_rp) then
            if (rank == 0) call g_logger%warning('Suspicious magbeta at atom '//fmt('i4', ia)//' -> resetting to 1.0', __FILE__, __LINE__)
            this%magbeta(ia) = 1.0d0
         end if
         ! Clamp to [0,1]
         this%magbeta(ia) = min(1.0d0, max(0.0d0, this%magbeta(ia)))

         mag_mix(ia, :) = (1.d0 - this%magbeta(ia))*mag_old(ia, :) + this%magbeta(ia)*mag_new(ia, :)
      end do

      ! Quick diagnostic printout for the first few atoms to help debug mixing
      ! if (rank == 0) then
      !    do ia = 1, min(3, this%lattice%nrec)
      !       call g_logger%info(’mix diagnostics atom ’//fmt(’i4’, ia)//″: magbeta=’″//fmt(’f6.3’, this%magbeta(ia)) &
      !            //’ mag_old=’//fmt(’f10.6’, mag_old(ia,1))//’ ’//fmt(’f10.6’, mag_old(ia,2))//’ ’//fmt(’f10.6’, mag_old(ia,3)) &
      !            , __FILE__, __LINE__)
      !       call g_logger%info(’mix diagnostics atom ’//fmt(’i4’, ia)//″: mag_new=’″//fmt(’f10.6’, mag_new(ia,1))//’ ’//fmt(’f10.6’, mag_new(ia,2))//’ ’//fmt(’f10.6’, mag_new(ia,3)) &
      !            , __FILE__, __LINE__)
      !       call g_logger%info(’mix diagnostics atom ’//fmt(’i4’, ia)//″: mag_mix=’″//fmt(’f10.6’, mag_mix(ia,1))//’ ’//fmt(’f10.6’, mag_mix(ia,2))//’ ’//fmt(’f10.6’, mag_mix(ia,3)) &
      !            , __FILE__, __LINE__)
      !    end do
      ! end if
   end subroutine mix_magnetic_moments

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Mix the parameters QL and PL. Can be linear, andersen or broyden.
   !>
   !> Mix the parameters QL and PL. Can be linear, andersen or broyden.
   !> @param[in] dummy variable of type mix
   !> @param[in] old and new potential variables qia_old and qia_new, respectively.
   !> @return type(mix): qia, the mix potential variable.
   !---------------------------------------------------------------------------
   subroutine mixpq(this, qia_old, qia_new, iteration, guard_collinear)
      use mpi_mod
      class(mix), intent(inout) :: this
      real(rp), dimension(:, :), intent(in) :: qia_old, qia_new
      integer, intent(in), optional :: iteration
      logical, intent(in), optional :: guard_collinear
      real(rp), allocatable :: qi_to(:, :), qi_tn(:, :) ! Local variables for broyden mixing
      real(rp) :: delta_atom
      integer :: ia ! Atom index
      logical :: reset
      real(rp) :: Bnorm
      integer :: nb_slice, lcount, occ_cols, packed_cols
      integer :: ncols
      real(rp) :: denom
      integer :: broy_width
      integer :: iteration_value
      logical :: guard_active

      iteration_value = 0
      guard_active = .false.
      if (present(iteration)) iteration_value = iteration
      if (present(guard_collinear)) guard_active = guard_collinear

      if (size(qia_old, 1) /= this%lattice%nrec .or. size(qia_new, 1) /= this%lattice%nrec) then
         call g_logger%fatal('mixpq received qia arrays with inconsistent atom count', __FILE__, __LINE__)
      end if
      if (size(qia_old, 2) /= size(qia_new, 2) .or. size(qia_old, 2) /= size(this%qia, 2)) then
         call g_logger%fatal('mixpq received qia arrays with inconsistent packed widths', __FILE__, __LINE__)
      end if
      broy_width = this%lattice%nrec*size(qia_old, 2)

      if (guard_active) then
         call validate_packed_q(this, qia_old, iteration_value, 'old')
         call validate_packed_q(this, qia_new, iteration_value, 'new')
      end if

      select case (trim(this%mixtype))
      case ('linear')
         this%qia(:, :) = (1.0d0 - this%beta)*qia_old(:, :) + this%beta*qia_new(:, :)
      case ('broyden')
         ! call g_logger%fatal(´Broyden mixing not implemented yet!´, __FILE__, __LINE__)
         reset = .false.
         allocate(qi_to(size(qia_old, 1), size(qia_old, 2)))
         allocate(qi_tn(size(qia_new, 1), size(qia_new, 2)))
         qi_to(:, :) = this%qia_old(:, :)
         qi_tn(:, :) = this%qia_new(:, :)
         call broydn(this%beta, this%beta, reset, qi_to, qi_tn, Bnorm, broy_width, this%itr, this%fsqo, this%u_broy, this%v_broy, this%muo_broy, this%fo_broy, this%nmix)
         this%qia(:, :) = qi_to(:, :)
         bnorm = bnorm**0.5d0
         deallocate(qi_to, qi_tn)
      end select
      ! Safety: detect NaNs produced by mixing (e.g., in Broyden) and fallback
      if (guard_active) then
         call validate_packed_q(this, this%qia, iteration_value, 'mixed')
      else if (any(this%qia /= this%qia)) then
         call g_logger%error('NaN detected in mixed qia — reverting to linear mix', __FILE__, __LINE__)
         this%qia(:, :) = (1.0d0 - this%beta)*qia_old(:, :) + this%beta*qia_new(:, :)
      end if

      ! Debug: print first atoms qia slices for inspection (up to 6 cols)
      ! if (rank == 0) then
      !    do ia = 1, min(3, this%lattice%nrec)
      !       nprint = min(6, size(qia_old, 2))
      !       msg = ’’
      !       do k = 1, nprint
      !          msg = msg // ’ ’ // fmt(’f10.6’, qia_old(ia, k))
      !       end do
      !       call g_logger%info(’mixpq diagnostics atom ’//fmt(’i4’, ia)//’: qia_old=’//trim(msg), __FILE__, __LINE__)

      !       nprint = min(6, size(qia_new, 2))
      !       msg = ’’
      !       do k = 1, nprint
      !          msg = msg // ’ ’ // fmt(’f10.6’, qia_new(ia, k))
      !       end do
      !       call g_logger%info(’mixpq diagnostics atom ’//fmt(’i4’, ia)//’: qia_new=’//trim(msg), __FILE__, __LINE__)

      !       nprint = min(6, size(this%qia, 2))
      !       msg = ’’
      !       do k = 1, nprint
      !          msg = msg // ’ ’ // fmt(’f10.6’, this%qia(ia, k))
      !       end do
      !       call g_logger%info(’mixpq diagnostics atom ’//fmt(’i4’, ia)//’: qia_mixed=’//trim(msg), __FILE__, __LINE__)
      !    end do
      ! end if
      this%charge%dq(:) = 0.0d0
      do ia = 1, this%lattice%nrec
         lcount = min(this%symbolic_atom(this%lattice%nbulk + ia)%potential%lmax, lmax_basis) + 1
         occ_cols = min(2*lcount, size(this%qia, 2))
         this%charge%trq = 0.0d0
         this%charge%cht = 0.0d0
         this%charge%cht = sum(this%qia_new(ia, 1:occ_cols)) - this%symbolic_atom(this%lattice%nbulk + ia)%element%valence
         this%charge%dq(ia) = sum(this%qia(ia, 1:occ_cols)) - this%symbolic_atom(this%lattice%nbulk + ia)%element%valence
         this%charge%trq = this%charge%dq(ia)
         if (rank == 0) call g_logger%info('Valence at atom'//fmt('i4', (ia))//' '//'is'//' '//fmt('f10.6', (this%symbolic_atom(this%lattice%nbulk + ia)%element%valence)), __FILE__, __LINE__)
         if (rank == 0) call g_logger%info('Charge transfer at atom '//fmt('i4', (ia))//': '// &
                                           fmt('f10.6', this%charge%cht)//' '//fmt('f10.6', this%charge%trq)//' '//fmt('f10.6', this%charge%cht - this%charge%trq) &
                                           , __FILE__, __LINE__)
         !if(rank==0) call g_logger%info(´Charge transfer at atom ´//fmt(´i4´,(ia)), __FILE__, __LINE__)
         !if(rank==0) call g_logger%info(fmt(´f10.6´,this%charge%cht)//´ ´//fmt(´f10.6´,this%charge%trq)//´ ´//fmt(´f10.6´,this%charge%cht-this%charge%trq), __FILE__, __LINE__)
      end do

      do ia = 1, this%lattice%nrec
         lcount = min(this%symbolic_atom(this%lattice%nbulk + ia)%potential%lmax, lmax_basis) + 1
         packed_cols = min(6*lcount, size(qia_old, 2))
         denom = max(1.0_rp, real(packed_cols)/2.0_rp)
         delta_atom = sqrt(sum((this%qia_old(ia, 1:packed_cols) - this%qia_new(ia, 1:packed_cols))**2))/denom
         if (rank == 0) call g_logger%info('Moment diff for atom:'//fmt('i4', (ia))//' '//'is'//' '//fmt('f10.6', delta_atom), __FILE__, __LINE__)
      end do
      ncols = size(qia_old, 2)
      nb_slice = min(6*(lmax_basis + 1), ncols)
      denom = max(1.0_rp, real(nb_slice)/2.0_rp)
      this%delta = sqrt(sum((this%qia_old(:, 1:nb_slice) - this%qia_new(:, 1:nb_slice))**2))/denom/this%lattice%nrec
   end subroutine mixpq

   subroutine validate_packed_q(this, values, iteration, stage)
      class(mix), intent(in) :: this
      real(rp), intent(in) :: values(:, :)
      integer, intent(in) :: iteration
      character(len=*), intent(in) :: stage
      integer :: it, lcount, l, spin, q0_col, q2_col, pl_col, qcols
      character(len=48) :: field

      do it = 1, this%lattice%nrec
         lcount = min(this%symbolic_atom(this%lattice%nbulk + it)%potential%lmax, lmax_basis) + 1
         qcols = size(values, 2)
         do spin = 1, 2
            do l = 0, lcount - 1
               q0_col = l + 1 + (spin - 1)*lcount
               q2_col = l + 1 + (spin + 1)*lcount
               pl_col = l + 1 + (spin + 3)*lcount
               if (q0_col <= qcols) then
                  call guard_q0_capacity(this%symbolic_atom(this%lattice%nbulk + it)%element%symbol, &
                     it, l, spin, iteration, values(it, q0_col))
               end if
               if (q2_col <= qcols) then
                  field = trim(stage)//' q2'
                  call guard_nonfinite_scalar(this%symbolic_atom(this%lattice%nbulk + it)%element%symbol, &
                     it, l, spin, trim(field), iteration, values(it, q2_col))
               end if
               if (pl_col <= qcols) then
                  field = trim(stage)//' PL'
                  call guard_nonfinite_scalar(this%symbolic_atom(this%lattice%nbulk + it)%element%symbol, &
                     it, l, spin, trim(field), iteration, values(it, pl_col))
               end if
            end do
         end do
      end do
   end subroutine validate_packed_q

   subroutine save_ldm_to(this, whereto)
      class(mix), intent(inout) :: this
      character(len=*), intent(in) :: whereto
      integer :: it, ispin, l, lcount, max_flat, local_flat

      max_flat = size(this%ldm_old, 4)

      select case (trim(whereto))
      case ('old', 'new')
         do it = 1, this%lattice%nrec
            call this%symbolic_atom(this%lattice%nbulk + it)%potential%flatten_ldm()
            lcount = min(this%symbolic_atom(this%lattice%nbulk + it)%potential%lmax, lmax_basis) + 1
            if (trim(whereto) == 'old') then
               this%ldm_old(it, :, :, :) = 0.0_rp
            else
               this%ldm_new(it, :, :, :) = 0.0_rp
            end if
            do l = 1, lcount
               local_flat = (2*l - 1)*(2*l - 1)
               do ispin = 1, 2
                  if (trim(whereto) == 'old') then
                     this%ldm_old(it, l, ispin, 1:min(local_flat, max_flat)) = cmplx( &
                        this%symbolic_atom(this%lattice%nbulk + it)%potential%ldm_flatten(l, ispin, 1:min(local_flat, max_flat)), &
                        this%symbolic_atom(this%lattice%nbulk + it)%potential%ldm_flatten_imag(l, ispin, 1:min(local_flat, max_flat)), rp)
                  else
                     this%ldm_new(it, l, ispin, 1:min(local_flat, max_flat)) = cmplx( &
                        this%symbolic_atom(this%lattice%nbulk + it)%potential%ldm_flatten(l, ispin, 1:min(local_flat, max_flat)), &
                        this%symbolic_atom(this%lattice%nbulk + it)%potential%ldm_flatten_imag(l, ispin, 1:min(local_flat, max_flat)), rp)
                  end if
               end do
            end do
         end do

      case ('current')
         do it = 1, this%lattice%nrec
            lcount = min(this%symbolic_atom(this%lattice%nbulk + it)%potential%lmax, lmax_basis) + 1
            this%symbolic_atom(this%lattice%nbulk + it)%potential%ldm_flatten(:, :, :) = 0.0_rp
            this%symbolic_atom(this%lattice%nbulk + it)%potential%ldm_flatten_imag(:, :, :) = 0.0_rp
            do l = 1, lcount
               local_flat = (2*l - 1)*(2*l - 1)
               do ispin = 1, 2
                  this%symbolic_atom(this%lattice%nbulk + it)%potential%ldm_flatten(l, ispin, 1:min(local_flat, max_flat)) = &
                     real(this%ldm_mix(it, l, ispin, 1:min(local_flat, max_flat)), rp)
                  this%symbolic_atom(this%lattice%nbulk + it)%potential%ldm_flatten_imag(l, ispin, 1:min(local_flat, max_flat)) = &
                     aimag(this%ldm_mix(it, l, ispin, 1:min(local_flat, max_flat)))
               end do
            end do
            call this%symbolic_atom(this%lattice%nbulk + it)%potential%expand_ldm()
         end do
      end select
   end subroutine save_ldm_to

   subroutine mix_ldm_linear(this)
      class(mix), intent(inout) :: this
      real(rp) :: beta_loc

      beta_loc = min(1.0_rp, max(0.0_rp, this%ldm_beta))
      this%ldm_mix(:, :, :, :) = (1.0_rp - beta_loc)*this%ldm_old(:, :, :, :) + beta_loc*this%ldm_new(:, :, :, :)
   end subroutine mix_ldm_linear
   subroutine broydn(pmix, amix, reset, mu, f, fsq, imu, itr, fsqo, u, v, muo, fo, nmix)
      ! mix potentials with broydens jacobian updating method as
      ! implemented by g.p.srivastava, j. phys a 17, l317(1984),
      ! changed 021015 by Lars Bergqvist in order to avoid any writing on files
      ! eqs. 12 through 16, and yesim n. darici 8/85.
      ! input  : mu =  input potential vector ( destroyed )
      !          f  = output potential vector ( destroyed )
      ! output : mu =  new input vector
      ! ierr :  0 normal processing
      !         1 input reset
      !
      !.. Implicit Declarations ..
      implicit none
      !
      !.. Formal Arguments ..
      logical, intent(inout) :: reset
      integer, intent(in) :: imu
      real(rp), intent(in) :: pmix
      real(rp), intent(in) :: amix
      real(rp), intent(out) :: fsq
      real(rp), dimension(imu), intent(inout) :: f, mu, u, v, muo, fo
      integer, intent(inout) :: itr
      real(rp), intent(inout) :: fsqo
      integer, intent(in) :: nmix
      !
      !.. Local Scalars ..
      logical :: wjac
      integer :: i, ierr, it, n, itrn, diverg, myid
      real(rp) :: df2, t
      !
      !.. Local Arrays ..
      real(rp), dimension(imu) :: &
         df, dmu, w1, w2
      !
      !.. External Functions ..
      real(rp), external :: ddot
      !
      ! ... Executable Statements ...
      !
      myid = 0
      !
      wjac = .true.
      diverg = 0
      ierr = 0
      !  n = 0
      n = imu
      !
      !  do i = 1, nvs
      !     n = n + nves(i)
      !  end do
      do i = 1, n
         f(i) = f(i) - mu(i)
      end do
      fsq = ddot(n, f, 1, f, 1)/n
      !print *, ´BROYDN:´, fsq, fsqo, fsq>fsqo
      !print *, ´BROYDN:´, itr, nmix, reset
      if (.not. reset) then
         if (itr == 0 .or. fsq > fsqo) then
            ! if (itr==0 .or. itr>nmix .or. fsq>fsqo) then
            ! if (itr==0 .or. itr>nmix) then
            ! if (itr==0 .or. itr>4 .or. fsq>fsqo) then
            ! if (itr==0 .or. itr>10 .or. fsq>fsqo) then
            reset = .true.
            if (fsq > fsqo) then
               diverg = 1
            end if
         end if
      end if
      if (reset) then
         itr = 0
         ierr = 1
      end if
      !print *, ´Bmix mid´, itr, reset, diverg
      ! do i = 1, nvs
      !     pmix(i) = amix
      !     rmix(i) = pmix(i)
      ! end do
      if (itr /= 0) then
         ! read mu and f from previous iteration
         dmu = muo
         df = fo
      end if
    !!! if (myid==0) then
    !!!    write (6, ″(/ a / a / 1x, i6, 1pe12.4, 2i12 /)″) &
    !!!         ″ Broydn:″, ″    itr         fsq      length        ierr″, itr, fsq, &
    !!!         n, ierr
    !!!    write (0, ″(/ a / a / 1x, i6, 1pe12.4, 2i12 /)″) &
    !!!         ″ Broydn:″, ″    itr         fsq      length        ierr″, itr, fsq, &
    !!!         n, ierr
    !!! end if
      !write (876, ″(/ a / a / 1x, i6, 1pe12.4, 2i12 /)″) &
      !      ″ Broydn:″, ″    itr         fsq      length        ierr″, itr, fsq, &
      !      n, ierr
      if (wjac) then
         !
         itrn = itr + 1
         muo = mu
         fo = f
      end if
      ! set mix vector
      ! k = 0
      ! do i = 1, nvs
      !     do j = 1, nves(i)
      !       k = k + 1
      !       mixv(k) = rmix(i)
      !     end do
      ! end do
      !
      !print *, ´BROYDEN´, itr, itrn, n, nmix
      if (itr == 0) then
         do i = 1, n
            ! mu(i) = mu(i) + mixv(i)*f(i)
            mu(i) = mu(i) + pmix*f(i)
            ! mu(i) = mu(i) + 0.001d0*f(i)
         end do
      elseif (itr == 1) then
         do i = 1, n
            ! u(i) = mu(i) - dmu(i) + mixv(i)*(f(i)-df(i))
            u(i) = mu(i) - dmu(i) + amix*(f(i) - df(i))
            v(i) = f(i) - df(i)
         end do
         !
         df2 = ddot(n, v, 1, v, 1)
         do i = 1, n
            v(i) = v(i)/df2
         end do
         !
         t = ddot(n, v, 1, f, 1)
         do i = 1, n
            ! mu(i) = mu(i) + mixv(i)*f(i) - u(i)*t
            mu(i) = mu(i) + amix*f(i) - u(i)*t
         end do
      else
         ! itr > 1
         do i = 1, n
            dmu(i) = mu(i) - dmu(i)
            df(i) = f(i) - df(i)
            w1(i) = 0
            w2(i) = 0
         end do
         !
         do it = 1, itr - 1
            !
            t = ddot(n, v, 1, f, 1)
            do i = 1, n
               w1(i) = w1(i) + u(i)*t
            end do
            !
            t = ddot(n, v, 1, df, 1)
            do i = 1, n
               w2(i) = w2(i) + u(i)*t
            end do
            !
         end do
         !
         do i = 1, n
            ! u(i) = dmu(i) + mixv(i)*df(i) - w2(i)
            u(i) = dmu(i) + amix*df(i) - w2(i)
            v(i) = df(i)
         end do
         !
         df2 = ddot(n, v, 1, v, 1)
         do i = 1, n
            v(i) = v(i)/df2
         end do
         !
         t = ddot(n, v, 1, f, 1)
         do i = 1, n
            w1(i) = w1(i) + u(i)*t
         end do
         !
         do i = 1, n
            ! mu(i) = mu(i) + mixv(i)*f(i) - w1(i)
            !write(877, ´(i4, 3f12.6)´) i, mu(i), amix*f(i), w1(i)
            mu(i) = mu(i) + amix*f(i) - w1(i)
         end do
      end if
      itr = itrn
      fsqo = fsq
      if (itr > nmix) itr = 1
      !print *, ´Bmix done´, itr, fsq
   end subroutine broydn
end module mix_mod
