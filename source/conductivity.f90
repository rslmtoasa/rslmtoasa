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

   use iso_fortran_env, only: int64, real32
   use iso_c_binding, only: c_f_pointer, c_int, c_loc, c_long_long
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
   use kpm_profile_mod, only: g_kpm_profile
   use logger_mod, only: g_logger
   use cfd
   use basis_mod, only: nb, norb, spin_off
   implicit none

   private
   public :: gamma_mu_reference, gamma_mu_blas, gamma_mu_cblas, pack_gamma_mu_diagonal

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

   !> @brief Pack the diagonal orbital blocks consumed by Gamma*mu.
   !> @details The conductivity reconstruction only reads
   !> `mu_nm_stochastic(l,l,n,m,t)`.  The packed matrix uses the exact
   !> column-major pair mapping q = n + (m-1)*M, so it can be consumed by
   !> ZGEMM as U(q,l).  The workspace is supplied by the caller and is
   !> intentionally reusable across trace/type contributions.
   subroutine pack_gamma_mu_diagonal(mu_nm, u)
      complex(rp), intent(in) :: mu_nm(:, :, :, :)
      complex(rp), intent(out) :: u(:, :)

      integer :: l, m, n, q, nb_local, m_local

      nb_local = size(mu_nm, 1)
      m_local = size(mu_nm, 3)
      u = (0.0_rp, 0.0_rp)
      do l = 1, nb_local
         do m = 1, m_local
            do n = 1, m_local
               q = n + (m - 1)*m_local
               u(q, l) = mu_nm(l, l, n, m)
            end do
         end do
      end do
   end subroutine pack_gamma_mu_diagonal

   !> @brief Scalar reference for one stochastic trace/type contribution.
   !> @details This is the unoptimized production algebra retained for focused
   !> validation.  In particular, it has no conjugation and preserves the
   !> existing n-then-m summation order.
   subroutine gamma_mu_reference(gamma_nm, mu_nm, factor, c)
      complex(rp), intent(in) :: gamma_nm(:, :, :)
      complex(rp), intent(in) :: mu_nm(:, :, :, :)
      real(rp), intent(in) :: factor
      complex(rp), intent(out) :: c(:, :)

      integer :: i, l, m, n

      c = (0.0_rp, 0.0_rp)
      do i = 1, size(gamma_nm, 1)
         do n = 1, size(gamma_nm, 2)
            do m = 1, size(gamma_nm, 3)
               do l = 1, size(mu_nm, 1)
                  c(i, l) = c(i, l) + factor*gamma_nm(i, n, m)*mu_nm(l, l, n, m)
               end do
            end do
         end do
      end do
   end subroutine gamma_mu_reference

   !> @brief Reconstruct one trace/type contribution with a BLAS-3 contraction.
   !> @details `gamma_nm` is already allocated as (NE,M,M).  Its first
   !> dimension is contiguous in Fortran storage, so C_F_POINTER supplies a
   !> zero-copy view with shape (NE,M*M), where column q is the (n,m) pair
   !> q=n+(m-1)*M.  No conjugation or transpose is used: this computes
   !> C = factor * G * U exactly as the scalar Gamma*mu expression requires.
   !> `u` and `c` are caller-owned reusable workspaces.
   subroutine gamma_mu_blas(gamma_nm, u, factor, c)
      complex(rp), target, intent(in) :: gamma_nm(:, :, :)
      complex(rp), intent(in) :: u(:, :)
      real(rp), intent(in) :: factor
      complex(rp), intent(out) :: c(:, :)

      complex(rp), pointer, contiguous :: gamma_matrix(:, :)
      complex(rp) :: alpha, beta
      integer :: ne
      integer(c_int) :: shape(2)
      external :: zgemm

      ne = size(gamma_nm, 1)
      shape = [int(ne, c_int), int(size(gamma_nm, 2)*size(gamma_nm, 3), c_int)]
      call c_f_pointer(c_loc(gamma_nm), gamma_matrix, shape)

      alpha = cmplx(factor, 0.0_rp, kind=rp)
      beta = (0.0_rp, 0.0_rp)
      call zgemm('N', 'N', ne, size(u, 2), size(u, 1), alpha, gamma_matrix, ne, &
                 u, size(u, 1), beta, c, ne)
   end subroutine gamma_mu_blas

   !> @brief Single-precision BLAS analogue for precision-fair reconstruction.
   !> @details This is intentionally a small helper rather than a second
   !> transport path.  A caller that owns FP32 Gamma and mu workspaces can use
   !> CGEMM with the same column-major pair mapping as gamma_mu_blas.
   subroutine gamma_mu_cblas(gamma_matrix, u, factor, c)
      complex(real32), intent(in) :: gamma_matrix(:, :)
      complex(real32), intent(in) :: u(:, :)
      real(real32), intent(in) :: factor
      complex(real32), intent(out) :: c(:, :)

      complex(real32) :: alpha, beta
      integer :: ne
      external :: cgemm

      ne = size(gamma_matrix, 1)
      alpha = cmplx(factor, 0.0_real32, kind=real32)
      beta = (0.0_real32, 0.0_real32)
      call cgemm('N', 'N', ne, size(u, 2), size(u, 1), alpha, gamma_matrix, ne, &
                 u, size(u, 1), beta, c, ne)
   end subroutine gamma_mu_cblas

   !> @brief Batch the zero-temperature Fermi-weighted Simpson integral.
   !> @details The conductivity route always calls simpson_f with
   !>          fermi=.true., dfermi=.false., and T=0.  With the same energy
   !>          grid used for both Ene and EF, fermifun is exactly one below
   !>          the requested grid point, one half at that point, and zero
   !>          above it (the 1e-15 kBT floor is many orders below the grid
   !>          spacing).  Prefixing the fixed Simpson weights preserves that
   !>          calculation while reducing the repeated NE-by-NE work to
   !>          O(NE*ncolumn).
   subroutine integrate_fermi_batch(ene, y, result)
      real(rp), intent(in) :: ene(:), y(:, :)
      real(rp), intent(out) :: result(:, :)
      integer :: i, j, n, ncolumn, simpson_end
      real(rp) :: h, coefficient, running

      n = size(ene)
      ncolumn = size(y, 2)
      if (size(y, 1) /= n .or. size(result, 1) /= n .or. size(result, 2) /= ncolumn) then
         result = 0.0_rp
         return
      end if

      h = ene(2) - ene(1)
      ! simpson_f is called with NPTS = size(Ene)-9 and its last panel ends
      ! at NPTS+8 = size(Ene)-1.  The final grid point is retained for the
      ! Fermi-edge query but is not part of the Simpson domain.
      simpson_end = n - 1
      result = 0.0_rp
      do j = 1, ncolumn
         running = 0.0_rp
         do i = 1, n
            if (i == 1 .or. i == simpson_end) then
               coefficient = 1.0_rp
            else if (i < simpson_end .and. mod(i, 2) == 0) then
               coefficient = 4.0_rp
            else if (i < simpson_end) then
               ! Interior odd points occur as the upper endpoint of one
               ! Simpson panel and the lower endpoint of the next one.
               coefficient = 2.0_rp
            else
               coefficient = 0.0_rp
            end if
            ! Add the half-weight endpoint before advancing the prefix.  This
            ! retains the same ascending summation order as simpson_f for the
            ! zero-temperature grid-point Fermi edge.
            result(i, j) = h*(running + 0.5_rp*coefficient*y(i, j))/3.0_rp
            running = running + coefficient*y(i, j)
         end do
      end do
   end subroutine integrate_fermi_batch

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
      ! initialize constraining if requested for transport calculations
      if (associated(obj%control)) then
         if (obj%control%constraints_enable) then
            call initialize_cfd(obj%lattice%nrec, 1, obj%control%constraints_i_cons, obj%control%constraints_code_prefac)
         end if
      end if
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

      ! The optimized CUDA transport path retains only packed diagonal moments
      ! and constructs Gamma in energy tiles during reconstruction.  Keep the
      ! canonical host Gamma allocation available for CPU and diagnostic paths.
      if (associated(this%recursion%gpu_backend)) then
         if (this%recursion%gpu_backend%resident_moments_available()) then
            call g_logger%info('CUDA KPM reconstruction uses tiled resident Gamma blocks.', __FILE__, __LINE__)
            return
         end if
      end if

      ! The opt-in CPU FP32 route computes Gamma in single precision and
      ! leaves only the reduced conductivity result in the canonical FP64
      ! host representation.  The default route below remains unchanged.
      if (trim(this%control%cheb_backend) == 'fast' .and. &
          trim(this%control%cpu_reconstruction_precision) == 'fp32') then
         call calculate_gamma_nm_fp32(this)
         return
      end if

      call g_kpm_profile%start('P_gamma_basis_setup')
      
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

      call g_kpm_profile%start('D_gamma_basis')
      ! Precompute acos(x) and sqrt(1 - x^2) with scaled energy
      a = (this%en%energy_max - this%en%energy_min)/(2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min)/2

      wscale(:) = (this%en%ene(:) - b)/a
      acos_x(:) = acos(wscale(:))
      sqrt_term(:) = sqrt(1.0_rp - wscale(:)**2)

      ! Calculating the Jackson Kernel
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
      call g_kpm_profile%stop('D_gamma_basis')
      call g_kpm_profile%stop('P_gamma_basis_setup')

      ! Initialize Gamma_nm
      call g_kpm_profile%start('P_gamma_generation')
      call g_kpm_profile%start('D_gamma_fill')
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
      call g_kpm_profile%stop('D_gamma_fill')
      deallocate(acos_x, sqrt_term, chebyshev_poly, cn, cm, g_kernel, weights)
      call g_kpm_profile%stop('P_gamma_generation')
   end subroutine 

   !> @brief Build the CPU precision-matched FP32 Gamma tensor.
   !> @details This is deliberately local to the transport reconstruction
   !> route.  Gamma is formed with real32/complex(real32) arithmetic and
   !> widened only when stored in the existing canonical host tensor.
   subroutine calculate_gamma_nm_fp32(this)
      class(conductivity), intent(inout) :: this

      integer :: i, n, m, ne, m_cond
      real(real32) :: a_sp, b_sp, theta
      real(real32), allocatable :: g_kernel(:), weights(:), acos_x(:), sqrt_term(:), wscale(:)
      real(real32), allocatable :: chebyshev_poly(:, :)
      complex(real32), allocatable :: cn(:, :), cm(:, :), gamma_sp(:)

      ne = this%en%channels_ldos + 10
      m_cond = this%control%cond_ll
      call g_kpm_profile%start('P_gamma_basis_setup')
#ifdef USE_SAFE_ALLOC
      call g_safe_alloc%allocate('recursion.gamma_nm', this%gamma_nm, (/ne, m_cond, m_cond/))
#else
      allocate(this%gamma_nm(ne, m_cond, m_cond))
#endif
      allocate(acos_x(ne), sqrt_term(ne), wscale(ne), chebyshev_poly(ne, m_cond))
      allocate(cn(ne, m_cond), cm(ne, m_cond), g_kernel(m_cond), weights(m_cond), gamma_sp(ne))

      a_sp = real((this%en%energy_max - this%en%energy_min)/(2.0_rp - 0.3_rp), real32)
      b_sp = real((this%en%energy_max + this%en%energy_min)/2.0_rp, real32)
      wscale(:) = real((this%en%ene(:) - real(b_sp, rp))/real(a_sp, rp), real32)
      acos_x(:) = acos(wscale(:))
      sqrt_term(:) = sqrt(1.0_real32 - wscale(:)**2)
      do i = 1, m_cond
         theta = 6.0_real32 * (1.0_real32 - (real(i, real32) - 1.0_real32)/real(m_cond, real32))
         g_kernel(i) = sinh(theta)/sinh(6.0_real32)
         weights(i) = 1.0_real32
         if (i == 1) weights(i) = 0.5_real32
      end do
      do n = 1, m_cond
         cn(:, n) = (wscale(:) - cmplx(0.0_real32, real(n - 1, real32), real32)*sqrt_term(:)) * &
            exp(cmplx(0.0_real32, real(n - 1, real32), real32)*acos_x(:))
         cm(:, n) = (wscale(:) + cmplx(0.0_real32, real(n - 1, real32), real32)*sqrt_term(:)) * &
            exp(-cmplx(0.0_real32, real(n - 1, real32), real32)*acos_x(:))
      end do
      chebyshev_poly(:, 1) = 1.0_real32
      if (m_cond >= 2) chebyshev_poly(:, 2) = wscale(:)
      do n = 3, m_cond
         chebyshev_poly(:, n) = 2.0_real32*wscale(:)*chebyshev_poly(:, n - 1) - chebyshev_poly(:, n - 2)
      end do
      call g_kpm_profile%stop('P_gamma_basis_setup')

      call g_kpm_profile%start('P_gamma_generation')
      do n = 1, m_cond
         do m = 1, m_cond
            gamma_sp(:) = (cn(:, n)*chebyshev_poly(:, m) + cm(:, m)*chebyshev_poly(:, n)) / &
               ((1.0_real32 - wscale(:)**2)**2)
            gamma_sp(:) = gamma_sp(:)*g_kernel(n)*g_kernel(m)*weights(n)*weights(m)
            this%gamma_nm(:, n, m) = cmplx(real(gamma_sp(:), real32), aimag(gamma_sp(:)), kind=rp)
         end do
      end do
      deallocate(acos_x, sqrt_term, wscale, chebyshev_poly, cn, cm, g_kernel, weights, gamma_sp)
      call g_kpm_profile%stop('P_gamma_generation')
   end subroutine calculate_gamma_nm_fp32

   subroutine calculate_conductivity_tensor(this)
      implicit none
      ! Input
      class(conductivity), intent(inout) :: this
      ! Local variables
      integer :: i, m, n, l1, l2, ntype, loop_over, ne, m_cond, k_cond
      complex(rp), dimension(:, :, :), allocatable :: integrand
      complex(rp), dimension(:, :, :, :), allocatable :: integrand_at
      complex(rp), dimension(:, :), allocatable :: mu_diag, gamma_mu
      complex(rp), dimension(:, :, :), allocatable :: gpu_reconstruction
      complex(real32), dimension(:, :), allocatable :: gamma_sp, mu_diag_sp, gamma_mu_sp
      real(rp), dimension(:, :), allocatable :: integrand_l_im, integrand_l_real
      real(rp), dimension(:), allocatable :: integrand_tot_real, integrand_tot_im, fermi_f, wscale, real_part_l, im_part_l
      real(rp), dimension(:), allocatable :: integrated_tot_real, integrated_tot_im
      real(rp), dimension(:), allocatable :: raw_integrand_tot_real, raw_integrand_tot_im
      real(rp), dimension(:, :), allocatable :: integrated_l_real, integrated_l_im
      real(rp), dimension(:, :), allocatable :: integrated_type_real, integrated_type_im
      real(rp), dimension(:, :, :), allocatable :: integrated_type_l_real, integrated_type_l_im
      real(rp), dimension(:, :), allocatable :: integration_y, integration_result
      complex(rp), dimension(nb) :: temp
      real(rp) :: a, b, real_part, im_part, factor, volume, de
      real(rp) :: gpu_gamma_seconds, gpu_gamma_basis_seconds, gpu_gamma_fill_seconds
      real(rp) :: gpu_gemm_seconds, gpu_result_d2h_seconds
      integer(c_long_long) :: gpu_gamma_h2d_bytes, gpu_gamma_block_bytes, gpu_result_d2h_bytes
      integer :: gpu_energy_block
      integer(int64) :: gpu_complex_bytes
      logical :: use_gpu_reconstruction, no_write
      logical :: use_cpu_fp32_reconstruction
      character(len=48), allocatable :: output_debug(:), output_total(:), output_type(:,:)
      character(len=304), allocatable :: output_orb_real(:), output_orb_im(:)
      character(len=304), allocatable :: output_type_orb_real(:,:), output_type_orb_im(:,:)
      ! Printing variables
      character(len=*), parameter :: fname_cond_total = "cond_total.out"
      character(len=*), parameter :: fname_cond_orb_real = "cond_total_orb_real.out"
      character(len=*), parameter :: fname_cond_orb_im   = "cond_total_orb_im.out"
      character(len=sl) :: fname_r, fname_i, fname_orb_r, fname_orb_i, symbol

      allocate(integrand(nb, nb, this%en%channels_ldos + 10), real_part_l(nb), im_part_l(nb))
      allocate(integrand_tot_real(this%en%channels_ldos + 10), integrand_tot_im(this%en%channels_ldos + 10))
      allocate(wscale(this%en%channels_ldos + 10))
      allocate(integrand_l_real(nb, this%en%channels_ldos + 10), integrand_l_im(nb, this%en%channels_ldos + 10))
      allocate(integrand_at(nb, nb, this%en%channels_ldos + 10, this%lattice%ntype))

      ne = this%en%channels_ldos + 10
      m_cond = this%control%cond_ll
      k_cond = m_cond*m_cond
      allocate(integrated_tot_real(ne), integrated_tot_im(ne))
      allocate(raw_integrand_tot_real(ne), raw_integrand_tot_im(ne))
      allocate(integrated_l_real(nb, ne), integrated_l_im(nb, ne))
      allocate(integration_y(ne, 2*(nb + 1)), integration_result(ne, 2*(nb + 1)))
      allocate(output_debug(ne), output_total(ne), output_orb_real(ne), output_orb_im(ne))

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
      de = this%en%energy_max - this%en%energy_min 

      wscale(:) = (this%en%ene(:) - b)/a

      select case(this%control%cond_calctype)
      case('per_type')
         loop_over = this%lattice%ntype
         allocate(integrated_type_real(ne, loop_over), integrated_type_im(ne, loop_over))
         allocate(integrated_type_l_real(nb, ne, loop_over), integrated_type_l_im(nb, ne, loop_over))
         allocate(output_type(ne, loop_over), output_type_orb_real(ne, loop_over), output_type_orb_im(ne, loop_over))
      case('random_vec')
         loop_over = this%control%random_vec_num
      end select  
      no_write = g_kpm_profile%output_suppressed()

      use_gpu_reconstruction = .false.
      if (associated(this%recursion%gpu_backend)) then
         use_gpu_reconstruction = this%recursion%gpu_backend%resident_moments_available(loop_over)
      end if
      use_cpu_fp32_reconstruction = .not. use_gpu_reconstruction .and. &
         trim(this%control%cheb_backend) == 'fast' .and. &
         trim(this%control%cpu_reconstruction_precision) == 'fp32'
      gpu_energy_block = 0
      if (use_gpu_reconstruction) then
         gpu_complex_bytes = int(merge(8, 16, trim(this%control%gpu_precision) == 'fp32'), int64)
         allocate(gpu_reconstruction(ne, nb, loop_over))
         gpu_reconstruction = (0.0_rp, 0.0_rp)
         call g_kpm_profile%set_reconstruction_bytes(0_int64, int(k_cond*nb, int64)*gpu_complex_bytes)
      else
         ! Only the diagonal orbital blocks enter Gamma*mu.  These two
         ! workspaces are allocated once per conductivity call and reused for
         ! every trace.
         allocate(mu_diag(k_cond, nb), gamma_mu(ne, nb))
         call g_kpm_profile%set_reconstruction_bytes( &
            int(size(this%gamma_nm), int64)*int(storage_size(this%gamma_nm)/8, int64), &
            int(size(mu_diag), int64)*int(storage_size(mu_diag)/8, int64))
         if (use_cpu_fp32_reconstruction) then
            allocate(gamma_sp(ne, k_cond), mu_diag_sp(k_cond, nb), gamma_mu_sp(ne, nb))
            gamma_sp = cmplx(reshape(this%gamma_nm, [ne, k_cond]), kind=real32)
         end if
      end if

      volume = dot_product(this%lattice%a(:, 1), (cross_product(this%lattice%a(:, 2), this%lattice%a(:, 3))))
      !factor = 1 !(e_const**2) * (hbar_const / e_const) / (hbar_const * volume)
      !factor = (hbar_const / e_const) * ((4 * e_const * hbar_const) / (pi * volume))
      !write(*,*) (hbar_const / e_const), ((4 * e_const * hbar_const) / (pi * volume)) 
      factor = 16 / (pi * (de**2))
      write(*,*) factor, volume, de
      !write(*,*) (16 * hbar_const * (e_const**2)) / (pi * volume * ((de * ry2joule)**2))
      if (use_gpu_reconstruction) then
         call this%recursion%gpu_backend%reconstruct_conductivity( &
            this%en%ene, a, b, factor, m_cond, loop_over, 0, gpu_reconstruction, &
            gpu_gamma_seconds, gpu_gamma_basis_seconds, gpu_gamma_fill_seconds, &
            gpu_gemm_seconds, gpu_result_d2h_seconds, gpu_gamma_h2d_bytes, &
            gpu_gamma_block_bytes, gpu_result_d2h_bytes, gpu_energy_block)
         call g_kpm_profile%add_seconds('P_gamma_basis_setup', gpu_gamma_basis_seconds)
         call g_kpm_profile%add_seconds('P_gamma_generation', gpu_gamma_fill_seconds)
         call g_kpm_profile%add_seconds('D_gamma_basis', gpu_gamma_basis_seconds)
         call g_kpm_profile%add_seconds('D_gamma_fill', gpu_gamma_fill_seconds)
         call g_kpm_profile%add_seconds('P_reconstruction_total', gpu_gemm_seconds + gpu_result_d2h_seconds)
         call g_kpm_profile%add_seconds('D_reconstruction_BLAS', gpu_gemm_seconds)
         call g_kpm_profile%add_seconds('D_reconstruction_D2H', gpu_result_d2h_seconds)
         call g_kpm_profile%add_bytes('H2D', int(gpu_gamma_h2d_bytes, int64))
         call g_kpm_profile%add_bytes('RESULT_D2H', int(gpu_result_d2h_bytes, int64))
         call g_kpm_profile%set_gpu_energy_block(gpu_energy_block)
         call g_kpm_profile%set_reconstruction_bytes( &
            int(gpu_gamma_block_bytes, int64), int(k_cond*nb, int64)*gpu_complex_bytes)

         call g_kpm_profile%start('P_result_unpack')
         do ntype = 1, loop_over
            do l2 = 1, nb
               integrand(l2, l2, :) = integrand(l2, l2, :) + gpu_reconstruction(:, l2, ntype)
               if (this%control%cond_calctype == 'per_type') then
                  integrand_at(l2, l2, :, ntype) = gpu_reconstruction(:, l2, ntype)
               end if
            end do
         end do
         call g_kpm_profile%stop('P_result_unpack')
      else
      call g_kpm_profile%start('P_reconstruction_total')
      do ntype = 1, loop_over
         call g_kpm_profile%start('D_mu_pack')
         call pack_gamma_mu_diagonal(this%recursion%mu_nm_stochastic(:, :, :, :, ntype), mu_diag)
         call g_kpm_profile%stop('D_mu_pack')

         ! Do not wrap this call in an outer OpenMP region.  BLAS owns the
         ! parallelism for this dense reconstruction, avoiding nested
         ! OpenMP/BLAS oversubscription in the CPU benchmark.
         call g_kpm_profile%start('D_reconstruction_BLAS')
         if (use_cpu_fp32_reconstruction) then
            mu_diag_sp = cmplx(mu_diag, kind=real32)
            call gamma_mu_cblas(gamma_sp, mu_diag_sp, real(factor, real32), gamma_mu_sp)
            gamma_mu = cmplx(gamma_mu_sp, kind=rp)
         else
            call gamma_mu_blas(this%gamma_nm, mu_diag, factor, gamma_mu)
         end if
         call g_kpm_profile%stop('D_reconstruction_BLAS')

         do l2 = 1, nb
            integrand(l2, l2, :) = integrand(l2, l2, :) + gamma_mu(:, l2)
            if (this%control%cond_calctype == 'per_type') then
               integrand_at(l2, l2, :, ntype) = gamma_mu(:, l2)
            end if
         end do
      end do
      call g_kpm_profile%stop('P_reconstruction_total')
      end if

      call g_kpm_profile%start('P_tensor_postprocess')
      integrand_tot_real(:) = 0.0d0
      integrand_tot_im(:) = 0.0d0
      do l2 = 1, nb
         integrand_tot_real(:) = integrand_tot_real(:) + real(integrand(l2, l2, :))
         integrand_tot_im(:) = integrand_tot_im(:) + aimag(integrand(l2, l2, :))
         integrand_l_real(l2, :) = real(integrand(l2, l2, :))
         integrand_l_im(l2, :) = aimag(integrand(l2, l2, :))
      end do
      raw_integrand_tot_real = integrand_tot_real
      raw_integrand_tot_im = integrand_tot_im
      call g_kpm_profile%stop('P_tensor_postprocess')

      call g_kpm_profile%start('P_energy_integration')
      integration_y(:, 1) = integrand_tot_real
      integration_y(:, 2) = integrand_tot_im
      do l2 = 1, nb
         integration_y(:, 2 + l2) = integrand_l_real(l2, :)
         integration_y(:, 2 + nb + l2) = integrand_l_im(l2, :)
      end do
      call integrate_fermi_batch(wscale, integration_y, integration_result)
      integrated_tot_real = integration_result(:, 1)
      integrated_tot_im = integration_result(:, 2)
      do l2 = 1, nb
         integrated_l_real(l2, :) = integration_result(:, 2 + l2)
         integrated_l_im(l2, :) = integration_result(:, 2 + nb + l2)
      end do
      if (this%control%cond_calctype == 'per_type') then
         do ntype = 1, loop_over
            integrand_tot_real(:) = 0.0d0
            integrand_tot_im(:) = 0.0d0
            do l2 = 1, nb
               integrand_tot_real(:) = integrand_tot_real(:) + real(integrand_at(l2, l2, :, ntype))
               integrand_tot_im(:) = integrand_tot_im(:) + aimag(integrand_at(l2, l2, :, ntype))
               integrand_l_real(l2, :) = real(integrand_at(l2, l2, :, ntype))
               integrand_l_im(l2, :) = aimag(integrand_at(l2, l2, :, ntype))
            end do
            integration_y(:, 1) = integrand_tot_real
            integration_y(:, 2) = integrand_tot_im
            do l2 = 1, nb
               integration_y(:, 2 + l2) = integrand_l_real(l2, :)
               integration_y(:, 2 + nb + l2) = integrand_l_im(l2, :)
            end do
            call integrate_fermi_batch(wscale, integration_y, integration_result)
            integrated_type_real(:, ntype) = integration_result(:, 1)
            integrated_type_im(:, ntype) = integration_result(:, 2)
            do l2 = 1, nb
               integrated_type_l_real(l2, :, ntype) = integration_result(:, 2 + l2)
               integrated_type_l_im(l2, :, ntype) = integration_result(:, 2 + nb + l2)
            end do
         end do
      end if
      call g_kpm_profile%stop('P_energy_integration')

      ! Format all output records before opening any output file. This keeps
      ! formatting/buffer preparation measurable and lets the benchmark-only
      ! no-write mode execute the same preparation without filesystem writes.
      call g_kpm_profile%start('P_output_prepare')
      do i = 1, ne
         write(output_debug(i), '(3es16.6)') (a*wscale(i)+b) - this%en%fermi, &
            raw_integrand_tot_real(i), raw_integrand_tot_im(i)
         write(output_total(i), '(3es16.6)') (a*wscale(i)+b) - this%en%fermi, &
            integrated_tot_real(i) / real(loop_over), integrated_tot_im(i) / real(loop_over)
         write(output_orb_real(i), '(19es16.6)') (a*wscale(i)+b) - this%en%fermi, &
            integrated_l_real(:, i) / real(loop_over)
         write(output_orb_im(i), '(19es16.6)') (a*wscale(i)+b) - this%en%fermi, &
            integrated_l_im(:, i) / real(loop_over)
      end do
      if (this%control%cond_calctype == 'per_type') then
         do ntype = 1, loop_over
            do i = 1, ne
               write(output_type(i, ntype), '(3es16.6)') (a*wscale(i)+b) - this%en%fermi, &
                  integrated_type_real(i, ntype), integrated_type_im(i, ntype)
               write(output_type_orb_real(i, ntype), '(19es16.6)') (a*wscale(i)+b) - this%en%fermi, &
                  integrated_type_l_real(:, i, ntype)
               write(output_type_orb_im(i, ntype), '(19es16.6)') (a*wscale(i)+b) - this%en%fermi, &
                  integrated_type_l_im(:, i, ntype)
            end do
         end do
      end if
      call g_kpm_profile%stop('P_output_prepare')

      call g_kpm_profile%start('P_output_io')
      if (.not. no_write) then
         open(unit=3, file=fname_cond_total, status='replace', action='write')
         open(unit=32, file=fname_cond_orb_real, status='replace', action='write')
         open(unit=33, file=fname_cond_orb_im, status='replace', action='write')
         do i = 1, ne
            write(123, '(A)') output_debug(i)
            write(3, '(A)') output_total(i)
            write(32, '(A)') output_orb_real(i)
            write(33, '(A)') output_orb_im(i)
         end do
         close(3)
         close(32)
         close(33)

         if (this%control%cond_calctype == 'per_type') then
            do ntype = 1, loop_over
               fname_r = trim(this%lattice%symbolic_atoms(ntype)%element%symbol) // '_cond.out'
               fname_orb_r = trim(this%lattice%symbolic_atoms(ntype)%element%symbol) // '_cond_orb_real.out'
               fname_orb_i = trim(this%lattice%symbolic_atoms(ntype)%element%symbol) // '_cond_orb_im.out'
               open(unit=100+ntype, file=fname_r, status='replace', action='write')
               open(unit=300+ntype, file=fname_orb_r, status='replace', action='write')
               open(unit=400+ntype, file=fname_orb_i, status='replace', action='write')
               do i = 1, ne
                  write(100+ntype, '(A)') output_type(i, ntype)
                  write(300+ntype, '(A)') output_type_orb_real(i, ntype)
                  write(400+ntype, '(A)') output_type_orb_im(i, ntype)
               end do
               close(100+ntype)
               close(300+ntype)
               close(400+ntype)
            end do
         end if
      end if
      call g_kpm_profile%stop('P_output_io')

      deallocate(integrand, integrand_tot_real, integrand_tot_im, raw_integrand_tot_real, raw_integrand_tot_im, &
         wscale, real_part_l, im_part_l, integrand_l_real, integrand_l_im, integrand_at, &
         integrated_tot_real, integrated_tot_im, integrated_l_real, integrated_l_im, &
         integration_y, integration_result, &
         output_debug, output_total, output_orb_real, output_orb_im)
      if (allocated(integrated_type_real)) deallocate(integrated_type_real, integrated_type_im, &
         integrated_type_l_real, integrated_type_l_im, output_type, output_type_orb_real, output_type_orb_im)
      if (allocated(gpu_reconstruction)) deallocate(gpu_reconstruction)
      if (allocated(mu_diag)) deallocate(mu_diag, gamma_mu)
      if (allocated(gamma_sp)) deallocate(gamma_sp, mu_diag_sp, gamma_mu_sp)
   end subroutine calculate_conductivity_tensor

end module conductivity_mod
