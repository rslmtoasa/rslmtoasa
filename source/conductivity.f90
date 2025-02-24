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

      ! Reading
      open (newunit=funit, file=this%control%fname, action='read', iostat=iostatus, status='old')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//trim(this%control%fname)//' not found', __FILE__, __LINE__)
      end if

      read (funit, nml=conductivity, iostat=iostatus)
      if (iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
         call g_logger%error('Error while reading namelist', __FILE__, __LINE__)
         call g_logger%error('iostatus = '//int2str(iostatus), __FILE__, __LINE__)
      end if
      close (funit)

      this%ll_cond = ll_cond
   end subroutine build_from_file   

   subroutine restore_to_default(this)
      class(conductivity), intent(inout) :: this

      this%ll_cond = 200
   end subroutine restore_to_default

   !*************************************************************************
   !> @brief Computes the Gamma_nm function for the Chebyshev polynomial expansion.
   !> 
   !> This subroutine calculates the Gamma_nm function used in the conductivity
   !> formula. The function incorporates the Jackson kernel for smoothing and uses
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
      real(rp), dimension(:), allocatable :: g_kernel(:)       ! Jackson kernel
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
      a = (this%en%energy_max - this%en%energy_min)/(2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min)/2

      wscale(:) = (this%en%ene(:) - b)/a

      acos_x(:) = acos(wscale(:))
      sqrt_term(:) = sqrt(1.0_rp - wscale(:)**2)

      ! Calculating the Jackson Kernel
      call jackson_kernel((this%control%cond_ll), g_kernel)

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
      real(rp) :: a, b, real_part, im_part, factor, volume
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
    
      a = (this%en%energy_max - this%en%energy_min)/(2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min)/2

      wscale(:) = (this%en%ene(:) - b)/a

      select case(this%control%cond_calctype)
      case('per_type')
         loop_over = this%lattice%ntype
      case('random_vec')
         loop_over = this%control%random_vec_num
      end select  

      ! factor = 2 * (e / pi) * omega^-1 
      volume = dot_product(this%lattice%a(:, 1), (cross_product(this%lattice%a(:, 2), this%lattice%a(:, 3)))) * (this%lattice%alat * ang2cent) ** 3
      factor = (e_const**2) * (hbar_const / e_const) / (hbar_const * volume)

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
         call simpson_f(real_part, wscale, wscale(i), this%en%nv1, integrand_tot_real(:), .true., .false., 0.0d0)
         call simpson_f(im_part, wscale, wscale(i), this%en%nv1, integrand_tot_im(:), .true., .false., 0.0d0)
         write(3, '(3f16.6)') (a*wscale(i)+b) - this%en%fermi, real_part / real(loop_over),  im_part / real(loop_over)
         do l2 = 1, 18
            call simpson_f(real_part_l(l2), wscale, wscale(i), this%en%nv1, integrand_l_real(l2, :), .true., .false., 0.0d0)
            call simpson_f(im_part_l(l2), wscale, wscale(i), this%en%nv1, integrand_l_im(l2, :), .true., .false., 0.0d0)
         end do
         write(32,'(19f16.6)') (a*wscale(i)+b) - this%en%fermi, real_part_l(1:18) / real(loop_over)
         write(33,'(19f16.6)') (a*wscale(i)+b) - this%en%fermi, im_part_l(1:18) / real(loop_over)
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
               call simpson_f(real_part, wscale, wscale(i), this%en%nv1, integrand_tot_real(:), .true., .false., 0.0d0)
               call simpson_f(im_part,   wscale, wscale(i), this%en%nv1, integrand_tot_im(:), .true., .false., 0.0d0)
               write(100+ntype, '(3f16.6)') (a*wscale(i)+b) - this%en%fermi, real_part, im_part
               do l2 = 1, 18
                  call simpson_f(real_part_l(l2), wscale, wscale(i), this%en%nv1, integrand_l_real(l2, :), .true., .false., 0.0d0)
                  call simpson_f(im_part_l(l2), wscale, wscale(i), this%en%nv1, integrand_l_im(l2, :), .true., .false., 0.0d0)
               end do
               write(300+ntype,'(19f16.6)') (a*wscale(i)+b) - this%en%fermi, real_part_l(1:18) 
               write(400+ntype,'(19f16.6)') (a*wscale(i)+b) - this%en%fermi, im_part_l(1:18) 
               end do
            close(100+ntype)
            close(300+ntype)
            close(400+ntype)
         end do  ! end do over ntype
      end if
      ! End writing statements

      deallocate(integrand, integrand_tot_real, integrand_tot_im, wscale, real_part_l, im_part_l, integrand_l_real, integrand_l_im, integrand_at)
   end subroutine calculate_conductivity_tensor

end module conductivity_mod
