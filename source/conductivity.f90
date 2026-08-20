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

      call g_kpm_profile%start('P_gamma')
      
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

      ! Initialize Gamma_nm
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
      call g_kpm_profile%stop('P_gamma')
   end subroutine 

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
      real(rp), dimension(:, :), allocatable :: integrand_l_im, integrand_l_real
      real(rp), dimension(:), allocatable :: integrand_tot_real, integrand_tot_im, fermi_f, wscale, real_part_l, im_part_l
      complex(rp), dimension(nb) :: temp
      real(rp) :: a, b, real_part, im_part, factor, volume, de
      real(rp) :: gpu_gamma_seconds, gpu_gamma_basis_seconds, gpu_gamma_fill_seconds
      real(rp) :: gpu_gemm_seconds, gpu_result_d2h_seconds
      integer(c_long_long) :: gpu_gamma_h2d_bytes, gpu_gamma_block_bytes, gpu_result_d2h_bytes
      integer :: gpu_energy_block
      integer(int64) :: gpu_complex_bytes
      logical :: use_gpu_reconstruction
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
      case('random_vec')
         loop_over = this%control%random_vec_num
      end select  

      use_gpu_reconstruction = .false.
      if (associated(this%recursion%gpu_backend)) then
         use_gpu_reconstruction = this%recursion%gpu_backend%resident_moments_available(loop_over)
      end if
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
         call g_kpm_profile%add_seconds('P_gamma', gpu_gamma_seconds)
         call g_kpm_profile%add_seconds('D_gamma_basis', gpu_gamma_basis_seconds)
         call g_kpm_profile%add_seconds('D_gamma_fill', gpu_gamma_fill_seconds)
         call g_kpm_profile%add_seconds('P_reconstruction_total', gpu_gemm_seconds + gpu_result_d2h_seconds)
         call g_kpm_profile%add_seconds('D_reconstruction_BLAS', gpu_gemm_seconds)
         call g_kpm_profile%add_seconds('D_reconstruction_D2H', gpu_result_d2h_seconds)
         call g_kpm_profile%add_bytes('H2D', int(gpu_gamma_h2d_bytes, int64))
         call g_kpm_profile%add_bytes('D2H', int(gpu_result_d2h_bytes, int64))
         call g_kpm_profile%set_gpu_energy_block(gpu_energy_block)
         call g_kpm_profile%set_reconstruction_bytes( &
            int(gpu_gamma_block_bytes, int64), int(k_cond*nb, int64)*gpu_complex_bytes)

         do ntype = 1, loop_over
            do l2 = 1, nb
               integrand(l2, l2, :) = integrand(l2, l2, :) + gpu_reconstruction(:, l2, ntype)
               if (this%control%cond_calctype == 'per_type') then
                  integrand_at(l2, l2, :, ntype) = gpu_reconstruction(:, l2, ntype)
               end if
            end do
         end do
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
         call gamma_mu_blas(this%gamma_nm, mu_diag, factor, gamma_mu)
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

      call g_kpm_profile%start('P_energy_integration')
      integrand_tot_real(:) = 0.0d0
      integrand_tot_im(:) = 0.0d0

      do l2 = 1, nb
         integrand_tot_real(:) = integrand_tot_real(:) + real(integrand(l2, l2, :))
         integrand_tot_im(:) = integrand_tot_im(:) + aimag(integrand(l2, l2, :))
         integrand_l_real(l2, :) = real(integrand(l2, l2, :))
         integrand_l_im(l2, :) = aimag(integrand(l2, l2, :))
      end do
      call g_kpm_profile%stop('P_energy_integration')

      ! Starting writing statements
      call g_kpm_profile%start('P_output_io')
      open(unit=3, file=fname_cond_total, status='replace', action='write')
      open(unit=32, file=fname_cond_orb_real, status='replace', action='write')
      open(unit=33, file=fname_cond_orb_im,   status='replace', action='write')

      do i = 1, this%en%channels_ldos + 10
         real_part = 0.0d0; im_part = 0.0d0; real_part_l(:) = 0.0d0; im_part_l(:) = 0.0d0
         write(123,'(3es16.6)') (a*wscale(i)+b) - this%en%fermi, integrand_tot_real(i), integrand_tot_im(i) 
         call simpson_f(real_part, wscale, wscale(i), this%en%nv1, integrand_tot_real(:), .true., .false., 0.0d0)
         call simpson_f(im_part, wscale, wscale(i), this%en%nv1, integrand_tot_im(:), .true., .false., 0.0d0)
         write(3, '(3es16.6)') (a*wscale(i)+b) - this%en%fermi, real_part / real(loop_over),  im_part / real(loop_over)
         do l2 = 1, nb
            call simpson_f(real_part_l(l2), wscale, wscale(i), this%en%nv1, integrand_l_real(l2, :), .true., .false., 0.0d0)
            call simpson_f(im_part_l(l2), wscale, wscale(i), this%en%nv1, integrand_l_im(l2, :), .true., .false., 0.0d0)
         end do
         write(32,'(19es16.6)') (a*wscale(i)+b) - this%en%fermi, real_part_l(1:nb) / real(loop_over)
         write(33,'(19es16.6)') (a*wscale(i)+b) - this%en%fermi, im_part_l(1:nb) / real(loop_over)
      end do
      close(3)
      close(32)
      close(33)
      call g_kpm_profile%stop('P_output_io')


      if (this%control%cond_calctype == 'per_type') then
         ! Loop over each atomic type
         do ntype = 1, loop_over

            call g_kpm_profile%start('P_energy_integration')
            integrand_tot_real(:) = 0.0d0
            integrand_tot_im(:)   = 0.0d0
            integrand_l_real(:, :) = 0.0d0
            integrand_l_im(:, :) = 0.0d0 

            do l2 = 1, nb
               integrand_tot_real(:) = integrand_tot_real(:) + real(integrand_at(l2, l2, :, ntype))
               integrand_tot_im(:)   = integrand_tot_im(:)   + aimag(integrand_at(l2, l2, :, ntype))
               integrand_l_real(l2, :) = real(integrand_at(l2, l2, :, ntype))
               integrand_l_im(l2, :) = aimag(integrand_at(l2, l2, :, ntype))
            end do
            call g_kpm_profile%stop('P_energy_integration')
         
            fname_r = trim(this%lattice%symbolic_atoms(ntype)%element%symbol) // "_cond.out"
            fname_orb_r = trim(this%lattice%symbolic_atoms(ntype)%element%symbol) // "_cond_orb_real.out"
            fname_orb_i = trim(this%lattice%symbolic_atoms(ntype)%element%symbol) // "_cond_orb_im.out"

            call g_kpm_profile%start('P_output_io')
            open(unit=100+ntype, file=fname_r, status='replace', action='write')
            open(unit=300+ntype, file=fname_orb_r, status='replace', action='write')
            open(unit=400+ntype, file=fname_orb_i, status='replace', action='write')
         
            do i = 1, this%en%channels_ldos + 10
               real_part = 0.0d0; im_part = 0.0d0; real_part_l(:) = 0.0d0; im_part_l(:) = 0.0d0
               ! Integrate over wscale for real and imaginary
               call simpson_f(real_part, wscale, wscale(i), this%en%nv1, integrand_tot_real(:), .true., .false., 0.0d0)
               call simpson_f(im_part,   wscale, wscale(i), this%en%nv1, integrand_tot_im(:), .true., .false., 0.0d0)
               write(100+ntype, '(3es16.6)') (a*wscale(i)+b) - this%en%fermi, real_part, im_part
               do l2 = 1, nb
                  call simpson_f(real_part_l(l2), wscale, wscale(i), this%en%nv1, integrand_l_real(l2, :), .true., .false., 0.0d0)
                  call simpson_f(im_part_l(l2), wscale, wscale(i), this%en%nv1, integrand_l_im(l2, :), .true., .false., 0.0d0)
               end do
               write(300+ntype,'(19es16.6)') (a*wscale(i)+b) - this%en%fermi, real_part_l(1:nb) 
               write(400+ntype,'(19es16.6)') (a*wscale(i)+b) - this%en%fermi, im_part_l(1:nb) 
               end do
            close(100+ntype)
            close(300+ntype)
            close(400+ntype)
            call g_kpm_profile%stop('P_output_io')
         end do  ! end do over ntype
      end if
      ! End writing statements

      deallocate(integrand, integrand_tot_real, integrand_tot_im, wscale, real_part_l, im_part_l, integrand_l_real, integrand_l_im, integrand_at)
      if (allocated(gpu_reconstruction)) deallocate(gpu_reconstruction)
      if (allocated(mu_diag)) deallocate(mu_diag, gamma_mu)
   end subroutine calculate_conductivity_tensor

end module conductivity_mod
