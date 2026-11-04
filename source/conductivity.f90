!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Conductivity
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
!> Module to handle the conductivity related processes
!------------------------------------------------------------------------------

module conductivity_mod

   use mpi_mod
   use control_mod
   use self_mod
   use energy_mod
   use lattice_mod
   use charge_mod
   use symbolic_atom_mod
   use hamiltonian_mod
   use recursion_mod
   use green_mod
   use density_of_states_mod
   use bands_mod
   use exchange_mod
   use mix_mod
   use math_mod
   use precision_mod
   use string_mod
   use self_mod
   use timer_mod, only: g_timer
   use logger_mod, only: g_logger
   implicit none

   private
   public :: conductivity_integral, response_operator_units
   type, public :: conductivity
      !> Control
      class(control), pointer :: control
      !> Lattice
      class(lattice), pointer :: lattice
      !> Self
      class(self), pointer :: self
      !> Energy 
      class(energy), pointer :: en
      !> Recursion 
      class(recursion), pointer :: recursion

      !> Number of recursion steps for the conductivity tensor
      integer :: ll_cond
      ! Bulk volume (Angstrom^3), sheet area (Angstrom^2), or per-site response.
      ! A positive measure is the TOTAL sampled volume/area, never per atom.
      character(len=16) :: response_geometry = 'bulk'
      real(rp) :: response_measure = 0.0_rp
      !> Pre factor Gamma_nm for the conductivity tensor calculation
      complex(rp), dimension(:, :, :), allocatable :: gamma_nm

   contains
      ! Destructor
      final :: destructor
      ! Procedures
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: calculate_gamma_nm
      procedure :: calculate_conductivity_tensor
      procedure :: response_normalization
      procedure :: write_physical_response
   end type
  
   interface conductivity
      procedure :: constructor
   end interface

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] self_obj Pointer to system´s self object
   !> @return type(conductivity) 
   !---------------------------------------------------------------------------
   function constructor(self_obj) result(obj)
      type(conductivity) :: obj
      class(self), target, intent(in) :: self_obj

      obj%self => self_obj
      obj%control => self_obj%control
      obj%lattice => self_obj%lattice
      obj%en => self_obj%en
      obj%recursion => self_obj%recursion

      call obj%restore_to_default()
      call obj%build_from_file()
   end function

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(conductivity) :: this
   end subroutine destructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Read parameters from input file
   !---------------------------------------------------------------------------
   subroutine build_from_file(this)
      class(conductivity), intent(inout) :: this

      ! variables associated with the reading processes
      integer :: iostatus, funit, i

      include 'include_codes/namelists/conductivity.f90'

      ll_cond = this%ll_cond
      response_geometry = this%response_geometry
      response_measure = this%response_measure

      ! Reading
      open (newunit=funit, file=this%control%fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(this%control%fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=conductivity, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%fatal('Error while reading conductivity namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      this%ll_cond = ll_cond
      this%response_geometry = response_geometry
      this%response_measure = response_measure
   end subroutine build_from_file   

   subroutine restore_to_default(this)
      class(conductivity), intent(inout) :: this

      this%ll_cond = 200
      this%response_geometry = 'bulk'
      this%response_measure = 0.0_rp
   end subroutine restore_to_default

   !*************************************************************************
   !> @brief Computes the Gamma_nm function for the Chebyshev polynomial expansion.
   !> 
   !> This subroutine calculates the Gamma_nm function used in the conductivity
   !> formula. The function incorporates the Lorentz kernel for smoothing and uses
   !> Chebyshev polynomials up to the specified maximum order. It also incorporates
   !> the (1 - energy**2)**2, Eq. (4) -> PRL 114, 116602 (2015).
   !> 
   !> @param[out] gamma_nm      Gamma_nm array (dimension: energy_grid, recursion_level, recursion_level).
   !*************************************************************************
   subroutine calculate_gamma_nm(this)
      class(conductivity), intent(inout) :: this

      ! Local variables
      integer :: i, n, m
      real(rp) :: a, b
      real(rp), dimension(:), allocatable :: g_kernel(:)       ! Lorentz damping kernel
      real(rp), dimension(:), allocatable :: weights(:)        ! Weight factors
      real(rp), dimension(:), allocatable :: acos_x, sqrt_term, wscale
      real(rp), dimension(:, :), allocatable :: chebyshev_poly
      complex(rp), dimension(:, :), allocatable :: cn, cm
      
      ! Initialize global variable
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('recursion.gamma_nm', this%gamma_nm, (/this%en%channels_ldos + 10, &
#else
      allocate (this%gamma_nm(this%en%channels_ldos + 10, this%lattice%control%cond_ll, this%lattice%control%cond_ll))
#endif

      ! Initialize variables
      allocate(acos_x(this%en%channels_ldos + 10), sqrt_term(this%en%channels_ldos + 10), wscale(this%en%channels_ldos + 10))
      allocate(chebyshev_poly(this%en%channels_ldos + 10, this%control%cond_ll))
      allocate(cn(this%en%channels_ldos + 10, this%control%cond_ll), cm(this%en%channels_ldos + 10, this%control%cond_ll))
      allocate(g_kernel(this%control%cond_ll), weights(this%control%cond_ll))

      ! Precompute acos(x) and sqrt(1 - x^2) with scaled energy
      a = (this%en%energy_max - this%en%energy_min)/(2.0_rp - 0.3_rp)
      b = (this%en%energy_max + this%en%energy_min)/2

      wscale(:) = (this%en%ene(:) - b)/a
      if (this%control%cond_ll < 2 .or. any(abs(wscale) >= 1.0_rp)) then
         call g_logger%fatal('Conductivity requires cond_ll>=2 and all scaled energies inside (-1,1)', __FILE__, __LINE__)
      end if
      acos_x(:) = acos(wscale(:))
      sqrt_term(:) = sqrt(1.0_rp - wscale(:)**2)

      ! Preserve the existing Lorentz damping (lambda=6); no kernel change.
      !call jackson_kernel((this%control%cond_ll), g_kernel)
      call lorentz_kernel(this%control%cond_ll, g_kernel, 6.0d0)
      ! Calculate weights
      weights(:) = 1.0d0
      weights(1) = 0.5d0

      ! Compute Cn and Cm
      do n = 1, this%control%cond_ll
         cn(:, n) = (wscale(:) - i_unit * real(n-1, rp) * sqrt_term(:)) * exp(i_unit * real(n-1, rp) * acos_x(:))
         cm(:, n) = (wscale(:) + i_unit * real(n-1, rp) * sqrt_term(:)) * exp(-i_unit * real(n-1, rp) * acos_x(:))
      end do

      chebyshev_poly(:, 1) = 1.0_rp
      chebyshev_poly(:, 2) = wscale(:)
      do n = 3, this%control%cond_ll
         chebyshev_poly(:, n) = 2.0_rp * wscale(:) * chebyshev_poly(:, n - 1) - chebyshev_poly(:, n - 2)
      end do

      ! Initialize Gamma_nm
      this%gamma_nm(:, :, :) = 0.0_rp

      ! Compute Gamma_nm
      do n = 1, this%control%cond_ll
         do m = 1, this%control%cond_ll
            this%gamma_nm(:, n, m) = (cn(:, n) * chebyshev_poly(:, m) + cm(:, m) * chebyshev_poly(:, n))
            this%gamma_nm(:, n, m) = this%gamma_nm(:, n, m) / ((1.0_rp - wscale(:)**2)**2)
            this%gamma_nm(:, n, m) = this%gamma_nm(:, n, m) * g_kernel(n) * g_kernel(m) * weights(n) * weights(m)
         end do
      end do

      ! Clean up
      deallocate(acos_x, sqrt_term, chebyshev_poly, cn, cm, g_kernel, weights)
   end subroutine 

   subroutine calculate_conductivity_tensor(this)
      implicit none
      ! Input
      class(conductivity), intent(inout) :: this
      ! Local variables
      integer :: i, m, n, l1, l2, ntype, loop_over
      complex(rp), dimension(:, :, :), allocatable :: integrand
      complex(rp), dimension(:, :, :, :), allocatable :: integrand_at
      real(rp), dimension(:, :), allocatable :: integrand_l_im, integrand_l_real
      real(rp), dimension(:), allocatable :: integrand_tot_real, integrand_tot_im, fermi_f, wscale, real_part_l, im_part_l
      complex(rp), dimension(18) :: temp
      real(rp) :: a, b, real_part, im_part, factor
      ! Printing variables
      character(len=*), parameter :: fname_cond_total = "cond_total.out"
      character(len=*), parameter :: fname_cond_orb_real = "cond_total_orb_real.out"
      character(len=*), parameter :: fname_cond_orb_im   = "cond_total_orb_im.out"
      character(len=sl) :: fname_r, fname_i, fname_orb_r, fname_orb_i, symbol

      allocate(integrand(18, 18, this%en%channels_ldos + 10), real_part_l(18), im_part_l(18))
      allocate(integrand_tot_real(this%en%channels_ldos + 10), integrand_tot_im(this%en%channels_ldos + 10))
      allocate(wscale(this%en%channels_ldos + 10))
      allocate(integrand_l_real(18, this%en%channels_ldos + 10), integrand_l_im(18, this%en%channels_ldos + 10))
      allocate(integrand_at(18, 18, this%en%channels_ldos + 10, this%lattice%ntype))
      
      integrand(:, :, :) = (0.0d0, 0.0d0)
      real_part_l(:) = 0.0d0
      im_part_l(:) = 0.0d0
      integrand_tot_real(:) = 0.0d0
      integrand_tot_im(:) = 0.0d0
      integrand_l_real(:, :) = 0.0d0
      integrand_l_im(:, :) = 0.0d0
      integrand_at(:, :, :, :) = (0.0d0, 0.0d0)
      temp(:) = 0.0d0
    
      a = (this%en%energy_max - this%en%energy_min)/(2.0_rp - 0.3_rp)
      b = (this%en%energy_max + this%en%energy_min)/2

      wscale(:) = (this%en%ene(:) - b)/a

      select case(this%control%cond_calctype)
      case('per_type')
         loop_over = this%lattice%ntype
      case('random_vec')
         loop_over = this%control%random_vec_num
      end select  

      ! x=(E-b)/a uses a=DeltaE/1.7, not DeltaE/2. The Bastin energy
      ! Jacobian is 1/a^2. Constants and spatial normalization are restored
      ! only in response_*.out; cond_*.out retain the bare numerical kernels.
      factor = 4.0_rp / (pi*a*a)
      do ntype = 1, loop_over
         !$omp parallel do default(shared) private(i, n, m, l2) schedule(dynamic)
         do i = 1, this%en%channels_ldos + 10
            
            ! Accumulate the contributions over the Chebyshev recursion indices.
            do n = 1, this%control%cond_ll
               do m = 1, this%control%cond_ll
                  do l2 = 1, 18
                     integrand(l2, l2, i) = integrand(l2, l2, i) + factor * this%gamma_nm(i, n, m) * &
                                this%recursion%mu_nm_stochastic(l2, l2, n, m, ntype)
                     if (this%control%cond_calctype == 'per_type') then
                        integrand_at(l2, l2, i, ntype) = integrand_at(l2, l2, i, ntype) + factor * this%gamma_nm(i, n, m) * &
                                this%recursion%mu_nm_stochastic(l2, l2, n, m, ntype)
                     end if
                   end do
               end do
            end do
         end do
         !$omp end parallel do
      end do

      integrand_tot_real(:) = 0.0d0
      integrand_tot_im(:) = 0.0d0

      do l2 = 1, 18
         integrand_tot_real(:) = integrand_tot_real(:) + real(integrand(l2, l2, :))
         integrand_tot_im(:) = integrand_tot_im(:) + aimag(integrand(l2, l2, :))
         integrand_l_real(l2, :) = real(integrand(l2, l2, :))
         integrand_l_im(l2, :) = aimag(integrand(l2, l2, :))
      end do

      ! Starting writing statements
      open(unit=3, file=fname_cond_total, status='replace', action='write')
      open(unit=32, file=fname_cond_orb_real, status='replace', action='write')
      open(unit=33, file=fname_cond_orb_im,   status='replace', action='write')

      do i = 1, this%en%channels_ldos + 10
         real_part = 0.0d0; im_part = 0.0d0; real_part_l(:) = 0.0d0; im_part_l(:) = 0.0d0
         write(123,'(3es16.6)') (a*wscale(i)+b) - this%en%fermi, integrand_tot_real(i), integrand_tot_im(i) 
         real_part = conductivity_integral(integrand_tot_real(:), wscale, i)
         im_part = conductivity_integral(integrand_tot_im(:), wscale, i)
         write(3, '(3es16.6)') (a*wscale(i)+b) - this%en%fermi, real_part / real(loop_over),  im_part / real(loop_over)
         do l2 = 1, 18
            real_part_l(l2) = conductivity_integral(integrand_l_real(l2, :), wscale, i)
            im_part_l(l2) = conductivity_integral(integrand_l_im(l2, :), wscale, i)
         end do
         write(32,'(19es16.6)') (a*wscale(i)+b) - this%en%fermi, real_part_l(1:18) / real(loop_over)
         write(33,'(19es16.6)') (a*wscale(i)+b) - this%en%fermi, im_part_l(1:18) / real(loop_over)
      end do


      if (this%control%cond_calctype == 'per_type') then
         ! Loop over each atomic type
         do ntype = 1, loop_over

            integrand_tot_real(:) = 0.0d0
            integrand_tot_im(:)   = 0.0d0
            integrand_l_real(:, :) = 0.0d0
            integrand_l_im(:, :) = 0.0d0 

            do l2 = 1, 18
               integrand_tot_real(:) = integrand_tot_real(:) + real(integrand_at(l2, l2, :, ntype))
               integrand_tot_im(:)   = integrand_tot_im(:)   + aimag(integrand_at(l2, l2, :, ntype))
               integrand_l_real(l2, :) = real(integrand_at(l2, l2, :, ntype))
               integrand_l_im(l2, :) = aimag(integrand_at(l2, l2, :, ntype))
            end do
         
            fname_r = trim(this%lattice%symbolic_atoms(ntype)%element%symbol) // "_cond.out"
            fname_orb_r = trim(this%lattice%symbolic_atoms(ntype)%element%symbol) // "_cond_orb_real.out"
            fname_orb_i = trim(this%lattice%symbolic_atoms(ntype)%element%symbol) // "_cond_orb_im.out"

            open(unit=100+ntype, file=fname_r, status='replace', action='write')
            open(unit=300+ntype, file=fname_orb_r, status='replace', action='write')
            open(unit=400+ntype, file=fname_orb_i, status='replace', action='write')
         
            do i = 1, this%en%channels_ldos + 10
               real_part = 0.0d0; im_part = 0.0d0; real_part_l(:) = 0.0d0; im_part_l(:) = 0.0d0
               ! Integrate over wscale for real and imaginary
               real_part = conductivity_integral(integrand_tot_real(:), wscale, i)
               im_part = conductivity_integral(integrand_tot_im(:), wscale, i)
               write(100+ntype, '(3es16.6)') (a*wscale(i)+b) - this%en%fermi, real_part, im_part
               do l2 = 1, 18
                  real_part_l(l2) = conductivity_integral(integrand_l_real(l2, :), wscale, i)
                  im_part_l(l2) = conductivity_integral(integrand_l_im(l2, :), wscale, i)
               end do
               write(300+ntype,'(19es16.6)') (a*wscale(i)+b) - this%en%fermi, real_part_l(1:18) 
               write(400+ntype,'(19es16.6)') (a*wscale(i)+b) - this%en%fermi, im_part_l(1:18) 
               end do
            close(100+ntype)
            close(300+ntype)
            close(400+ntype)
         end do  ! end do over ntype
      end if
      ! Close raw files explicitly before producing normalized outputs.
      close(3)
      close(32)
      close(33)
      call this%write_physical_response(integrand, integrand_at, wscale, a, b)
      ! End writing statements

      deallocate(integrand, integrand_tot_real, integrand_tot_im, wscale, real_part_l, im_part_l, integrand_l_real, integrand_l_im, integrand_at)
   end subroutine calculate_conductivity_tensor

   !---------------------------------------------------------------------------
   ! Integrate the zero-temperature occupied interval x(1):x(upper).
   ! Unlike simpson_f, assumed-shape arguments cannot invent ten extra points.
   ! Integrate quadratic interpolants, including a partial panel for an odd
   ! number of intervals. This is exact for quadratics on nonuniform grids too.
   ! There is no extra energy Jacobian here: x is already dimensionless.
   !---------------------------------------------------------------------------
   function conductivity_integral(y, x, upper) result(value)
      real(rp), intent(in) :: y(:), x(:)
      integer, intent(in) :: upper
      real(rp) :: value, h1, d1, d2, ta, tb
      integer :: i, j, last
      value = 0.0_rp
      if (size(y) /= size(x) .or. upper < 1 .or. upper > size(x)) then
         call g_logger%fatal('Invalid conductivity integration bounds', __FILE__, __LINE__)
      end if
      if (any(x(2:)-x(:size(x)-1) <= 0.0_rp)) then
         call g_logger%fatal('Conductivity energy grid must increase strictly', __FILE__, __LINE__)
      end if
      if (upper == 1) return
      if (size(x) == 2) then
         value = (x(2)-x(1))*(y(1)+y(2))/2
         return
      end if
      i = 1
      do while (i < upper)
         j = i
         last = min(i+2, upper)
         if (i+2 > size(x)) j = i-1
         h1 = x(j+1)-x(j)
         d1 = (y(j+1)-y(j))/h1
         d2 = ((y(j+2)-y(j+1))/(x(j+2)-x(j+1))-d1)/(x(j+2)-x(j))
         ta = x(i)-x(j)
         tb = x(last)-x(j)
         value = value + y(j)*(tb-ta) + d1*(tb*tb-ta*ta)/2 &
                 + d2*((tb**3-ta**3)/3-h1*(tb*tb-ta*ta)/2)
         i = last
      end do
   end function conductivity_integral

   ! Bare matrices carry Ry**power_e * Angstrom**power_l. The rank labels
   ! the physical observable q0=hbar**rank, not an extra factor in the matrix.
   ! Charge is a velocity here; e is restored exactly once at output.
   subroutine response_operator_units(name, power_e, power_l, rank)
      character(len=*), intent(in) :: name
      integer, intent(out) :: power_e, power_l, rank
      power_e = 1
      power_l = 1
      rank = 0
      select case(trim(name))
      case('charge')
         return
      case('spin','orbital','spin_accumulation','orbital_accumulation', &
           'spin_torque','orbital_torque','spin_soc_torque','soc_spin_torque')
         rank = 1
      case('quadrupole','oap','quadrupole_accumulation','oap_accumulation')
         rank = 2
      case('octupole','octupole_current','octupole_accumulation')
         rank = 3
      case default
         call g_logger%fatal('Unknown response operator units: '//trim(name), __FILE__, __LINE__)
      end select
      if (index(name,'accumulation') > 0) then
         power_e = 0
         power_l = 0
      else if (index(name,'torque') > 0) then
         power_l = 0
      end if
   end subroutine response_operator_units

   !---------------------------------------------------------------------------
   ! Convert the trace estimator to a density or per-site average. For per_type
   ! the representative trace is multiplied by the ACTUAL population of each
   ! type. For random_vec each normalized seed estimates Tr/N, so restore N.
   ! n3=1 does not select sheet units: a one-cell periodic stack is still bulk.
   ! Automatic measures require a complete generated Bravais supercell. For a
   ! cut cluster/custom geometry supply its total measure explicitly.
   !---------------------------------------------------------------------------
   subroutine response_normalization(this, weights, conversion, unit_label)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      class(conductivity), intent(in) :: this
      real(rp), intent(out) :: weights(:), conversion
      character(len=*), intent(out) :: unit_label
      real(rp) :: measure, density_si, area
      integer :: it, pe, pl, rank, pie, pil, ri, ndim
      character(len=32) :: qlabel

      if (this%lattice%kk <= 0) call g_logger%fatal('Empty response sample', __FILE__, __LINE__)
      measure = this%response_measure
      if (.not. ieee_is_finite(measure) .or. measure < 0) then
         call g_logger%fatal('response_measure must be finite and nonnegative', __FILE__, __LINE__)
      end if
      select case(trim(this%response_geometry))
      case('bulk')
         ndim = 3
      case('sheet')
         ndim = 2
      case('per_atom')
         ndim = 0
         if (measure /= 0) call g_logger%fatal('per_atom does not accept response_measure', __FILE__, __LINE__)
         measure = real(this%lattice%kk,rp)
      case default
         call g_logger%fatal('response_geometry must be bulk, sheet or per_atom', __FILE__, __LINE__)
      end select
      if (measure == 0) then
         if (.not. this%lattice%pbc) then
            call g_logger%fatal('Supply total response_measure for an open/cut sample', __FILE__, __LINE__)
         end if
         if (this%lattice%kk /= this%lattice%ntot*this%lattice%n1*this%lattice%n2*this%lattice%n3) then
            call g_logger%fatal('Incomplete supercell: supply total response_measure', __FILE__, __LINE__)
         end if
         area = norm2(cross_product(this%lattice%a(:,1),this%lattice%a(:,2)))*this%lattice%alat**2
         if (ndim == 2) then
            measure = area*this%lattice%n1*this%lattice%n2
         else
            measure = abs(dot_product(this%lattice%a(:,3),cross_product(this%lattice%a(:,1), &
                          this%lattice%a(:,2))))*this%lattice%alat**3 &
                          *this%lattice%n1*this%lattice%n2*this%lattice%n3
         end if
      end if
      if (.not. ieee_is_finite(measure) .or. measure <= 0) then
         call g_logger%fatal('Response volume/area must be positive; select sheet for a planar cell', __FILE__, __LINE__)
      end if
      select case(trim(this%control%cond_calctype))
      case('per_type')
         if (size(weights) /= this%lattice%ntype) error stop 'Response type weight size'
         do it = 1, size(weights)
            weights(it) = real(count(this%lattice%iz(1:this%lattice%kk)==it),rp)/measure
         end do
      case('random_vec')
         if (size(weights) /= this%control%random_vec_num) error stop 'Response random weight size'
         weights = real(this%lattice%kk,rp)/(measure*size(weights))
      case default
         call g_logger%fatal('Unknown response trace estimator', __FILE__, __LINE__)
      end select

      call response_operator_units(this%control%linear_out, pe, pl, rank)
      call response_operator_units(this%control%linear_in, pie, pil, ri)
      ! An arbitrary input operator has no implicitly defined electric field.
      ! Preserve this useful generalized response, but do not label it S/cm.
      if (trim(this%control%linear_in) /= 'charge') then
         conversion = 1.0_rp
         write(unit_label,'(a,i0,a,i0,a,a)') 'raw Ry^',pe+pie-2,' Angstrom^',pl+pil, &
               ' / normalization_measure; geometry=',trim(this%response_geometry)
         return
      end if
      write(qlabel,'(a,i0)') 'hbar^',rank
      ! Physical factors are absorbed in the NAMED unit for angular powers.
      ! This avoids tiny hbar^2/hbar^3 numbers and preserves their dimensions.
      if (pl == 1) then
         select case(ndim)
         case(3)
            conversion = e_const**2/hbar_const*1.0e8_rp
            unit_label = 'S/cm'
            if (rank > 0) unit_label = '('//trim(qlabel)//'/e) S/cm'
         case(2)
            conversion = 2*pi
            unit_label = 'e^2/h'
            if (rank > 0) unit_label = '('//trim(qlabel)//') e/h'
         case(0)
            conversion = e_const**2/hbar_const*1.0e-20_rp
            unit_label = 'S m^2 per atom'
            if (rank > 0) unit_label = '('//trim(qlabel)//'/e) S m^2 per atom'
         end select
      else
         ! Accumulations keep 1/Ry; torques contain [H,Q]/hbar, so replace Ry
         ! by hbar. Density weights above use Angstrom^ndim, hence 10^(10*ndim).
         density_si = 10.0_rp**(10*ndim)
         conversion = e_const*1.0e-10_rp/ry2joule*density_si
         unit_label = trim(qlabel)//' / (V/m)'
         if (pe == 1) then
            conversion = e_const*1.0e-10_rp/hbar_const*density_si
            unit_label = '('//trim(qlabel)//'/s) / (V/m)'
         end if
         select case(ndim)
         case(3)
            unit_label = trim(unit_label)//' per m^3'
         case(2)
            unit_label = trim(unit_label)//' per m^2'
         case(0)
            unit_label = trim(unit_label)//' per atom'
         end select
      end if
   end subroutine response_normalization

   !---------------------------------------------------------------------------
   ! Write explicitly labeled responses separately from the bare cond_* files.
   ! Constants use the existing positive-e Bastin convention: no silent Hall
   ! sign change. Multipoles are angular powers (q0=hbar**rank), not electric
   ! position multipoles. All 18 spin-orbitals are traced: no degeneracy factor.
   !---------------------------------------------------------------------------
   subroutine write_physical_response(this, integrand, integrand_at, x, a, b)
      class(conductivity), intent(in) :: this
      complex(rp), intent(in) :: integrand(:,:,:), integrand_at(:,:,:,:)
      real(rp), intent(in) :: x(:), a, b
      real(rp), allocatable :: weights(:), re(:,:), im(:,:), type_re(:), type_im(:)
      real(rp) :: conversion, energy, values_re(18), values_im(18)
      integer :: ntrace, it, l, i, info, total, orbital_re, orbital_im, type_unit
      integer :: pe, pl, rank, pie, pil, ri
      character(len=160) :: units
      character(len=80) :: filename

      ntrace = this%lattice%ntype
      if (trim(this%control%cond_calctype) == 'random_vec') ntrace = this%control%random_vec_num
      if (ntrace < 1) call g_logger%fatal('Response trace count must be positive', __FILE__, __LINE__)
      allocate(weights(ntrace), re(18,size(x)), im(18,size(x)), type_re(size(x)), type_im(size(x)))
      call this%response_normalization(weights, conversion, units)
      re = 0; im = 0
      if (trim(this%control%cond_calctype) == 'per_type') then
         do it = 1, ntrace
            type_re = 0; type_im = 0
            do l = 1, 18
               re(l,:) = re(l,:) + weights(it)*real(integrand_at(l,l,:,it))
               im(l,:) = im(l,:) + weights(it)*aimag(integrand_at(l,l,:,it))
               type_re = type_re + weights(it)*real(integrand_at(l,l,:,it))
               type_im = type_im + weights(it)*aimag(integrand_at(l,l,:,it))
            end do
            write(filename,'(a,i0,a)') 'response_type_',it,'.out'
            open(newunit=type_unit,file=trim(filename),status='replace')
            write(type_unit,'(a)') '# E-Ef [eV], Re, Im; additive type contribution; units: '//trim(units)
            do i = 1, size(x)
               write(type_unit,'(3es25.16)') (a*x(i)+b-this%en%fermi)*ry2ev, &
                    conversion*conductivity_integral(type_re,x,i), conversion*conductivity_integral(type_im,x,i)
            end do
            close(type_unit)
         end do
      else
         do l = 1, 18
            re(l,:) = weights(1)*real(integrand(l,l,:))
            im(l,:) = weights(1)*aimag(integrand(l,l,:))
         end do
      end if
      open(newunit=info,file='response_units.txt',status='replace')
      call response_operator_units(this%control%linear_out, pe, pl, rank)
      call response_operator_units(this%control%linear_in, pie, pil, ri)
      write(info,'(a,i0,a,i0)') 'raw kernel units = Ry^',pe+pie-2,' Angstrom^',pl+pil
      write(info,'(a)') 'linear_in = '//trim(this%control%linear_in)
      write(info,'(a)') 'linear_out = '//trim(this%control%linear_out)
      write(info,'(a)') 'geometry = '//trim(this%response_geometry)
      write(info,'(a)') 'units = '//trim(units)
      write(info,'(a)') 'energy column = E-Ef in eV; raw cond_* energy remains Ry'
      write(info,'(a)') 'multipoles = angular powers q0=hbar^rank; spin matrix already sigma/2'
      write(info,'(a)') 'sign = existing positive-e Bastin convention; no additional sign flip'
      write(info,'(a)') 'per_type = population-weighted representative-site approximation'
      write(info,'(a)') 'random_vec = normalized seeds estimate Tr/N; N restored in weights'
      write(info,'(a)') 'raw cond_* = bare kernels; historical unweighted type average, no volume or constants'
      write(info,'(a,es25.16)') 'normalization_measure (Angstrom^3, Angstrom^2, or site count) = ', &
                                    real(this%lattice%kk,rp)/sum(weights)
      write(info,'(a,3es25.16)') 'output direction = ',this%recursion%transport_dir_a
      write(info,'(a,3es25.16)') 'input direction = ',this%recursion%transport_dir_b
      write(info,'(a)') 'pol_alpha = '//trim(this%recursion%hamiltonian%pol_alpha)
      write(info,'(a)') 'pol_beta = '//trim(this%recursion%hamiltonian%pol_beta)
      write(info,'(a,es25.16)') 'output conversion = ',conversion
      write(info,'(a,*(es25.16,1x))') 'trace weights = ',weights
      close(info)
      open(newunit=total,file='response_total.out',status='replace')
      open(newunit=orbital_re,file='response_orbitals_real.out',status='replace')
      open(newunit=orbital_im,file='response_orbitals_imag.out',status='replace')
      write(total,'(a)') '# E-Ef [eV], Re, Im; units: '//trim(units)
      write(orbital_re,'(a)') '# E-Ef [eV], 18 spin-orbital contributions; units: '//trim(units)
      write(orbital_im,'(a)') '# E-Ef [eV], 18 spin-orbital contributions; units: '//trim(units)
      do i = 1, size(x)
         energy = (a*x(i)+b-this%en%fermi)*ry2ev
         do l = 1, 18
            values_re(l) = conversion*conductivity_integral(re(l,:),x,i)
            values_im(l) = conversion*conductivity_integral(im(l,:),x,i)
         end do
         write(total,'(3es25.16)') energy, sum(values_re), sum(values_im)
         write(orbital_re,'(19es25.16)') energy, values_re
         write(orbital_im,'(19es25.16)') energy, values_im
      end do
      close(total)
      close(orbital_re)
      close(orbital_im)
   end subroutine write_physical_response

end module conductivity_mod
