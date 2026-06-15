!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Green
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
!> Module to handle Greens functions calculation and related routines
!------------------------------------------------------------------------------

module green_mod

   use energy_mod
   use control_mod
   use lattice_mod
   use symbolic_atom_mod
   use recursion_mod
   use density_of_states_mod
   use chebyshev_fast_mod, only: cheb_green_fast
   use precision_mod, only: rp
   use math_mod
   use logger_mod, only: g_logger
   use timer_mod, only: g_timer
   use string_mod, only: fmt, int2str, real2str, log2str
   use recursion_gpu_mod, only: rsgpu
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Module´s main structure
   type, public :: green
      ! General variables
      !> Recursion
      class(recursion), pointer :: recursion
      !> Lattice
      class(lattice), pointer :: lattice
      !> Symbolic atom
      class(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Density of states
      class(dos), pointer :: dos
      !> Control
      class(control), pointer :: control
      !> Energy
      class(energy), pointer :: en
      !> Onsite Green Function
      complex(rp), dimension(:, :, :, :), allocatable :: g0
      !> Intersite Green Function
      complex(rp), dimension(:, :, :, :), allocatable :: gij, gji, ginmag, gjnmag, gix, giy, giz, gjx, gjy, gjz
      !> Intersite Green Function as a function of imaginary eta
      complex(rp), dimension(:, :, :, :), allocatable :: gij_eta, gji_eta, ginmag_eta, gjnmag_eta, gix_eta, giy_eta, giz_eta, gjx_eta, gjy_eta, gjz_eta
      !> Intersite Green Function with two indexes
      complex(rp), dimension(:, :, :, :), allocatable :: g00ij, g01ij, g00ji, g01ji, &
                                                       & gx1ij, gy1ij, gz1ij, gx0ij, gy0ij, gz0ij, &
                                                       & gx1ji, gy1ji, gz1ji, gx0ji, gy0ji, gz0ji

   contains
      procedure :: sgreen
      procedure :: bgreen
      procedure :: block_green
      procedure :: block_green_gpu
      procedure :: block_green_eta
      procedure :: block_green_ij
      procedure :: block_green_ij_gpu
      procedure :: block_green_ij_eta
      procedure :: calculate_intersite_gf
      procedure :: calculate_intersite_gf_twoindex
      procedure :: calculate_intersite_gf_eta
      procedure :: calculate_intersite_gf_eta_gpu
      procedure :: chebyshev_green
      procedure :: chebyshev_green_gpu
      procedure :: chebyshev_green_eta
      procedure :: chebyshev_green_ij
      procedure :: chebyshev_green_ij_gpu
      procedure :: chebyshev_green_ij_eta
      procedure :: chebyshev_dos_dispatch
      procedure :: restore_to_default
      procedure :: auxiliary_gij
      procedure :: transform_auxiliary_gij
      procedure :: gij_eta_to_gij
      final :: destructor
   end type

   interface green
      procedure :: constructor
   end interface green

contains

   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @param[in] fname Namelist file
   !> @return type(green)
   !---------------------------------------------------------------------------
   function constructor(dos_obj) result(obj)
      type(green) :: obj
      type(dos), target, intent(in) :: dos_obj

      obj%dos => dos_obj
      obj%recursion => dos_obj%recursion
      obj%en => dos_obj%en
      obj%symbolic_atom => dos_obj%recursion%hamiltonian%charge%lattice%symbolic_atoms
      obj%lattice => dos_obj%recursion%lattice
      obj%control => dos_obj%recursion%lattice%control

      call obj%restore_to_default()
   end function constructor

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Destructor
   !---------------------------------------------------------------------------
   subroutine destructor(this)
      type(green) :: this
#ifdef USE_SAFE_ALLOC
      if (allocated(this%g0)) call g_safe_alloc%deallocate('green.g0', this%g0)
      if (allocated(this%gij)) call g_safe_alloc%deallocate('green.gij', this%gij)
      if (allocated(this%gji)) call g_safe_alloc%deallocate('green.gji', this%gji)
      if (allocated(this%ginmag)) call g_safe_alloc%deallocate('green.ginmag', this%ginmag)
      if (allocated(this%gjnmag)) call g_safe_alloc%deallocate('green.gjnmag', this%gjnmag)
      if (allocated(this%gix)) call g_safe_alloc%deallocate('green.gix', this%gix)
      if (allocated(this%giy)) call g_safe_alloc%deallocate('green.giy', this%giy)
      if (allocated(this%giz)) call g_safe_alloc%deallocate('green.giz', this%giz)
      if (allocated(this%gjx)) call g_safe_alloc%deallocate('green.gjx', this%gjx)
      if (allocated(this%gjy)) call g_safe_alloc%deallocate('green.gjy', this%gjy)
      if (allocated(this%gjz)) call g_safe_alloc%deallocate('green.gjz', this%gjz)
      if (allocated(this%g00ij)) call g_safe_alloc%deallocate('green.g00ij', this%g00ij)
      if (allocated(this%g01ij)) call g_safe_alloc%deallocate('green.g01ij', this%g01ij)
      if (allocated(this%g00ji)) call g_safe_alloc%deallocate('green.g00ji', this%g00ji)
      if (allocated(this%g01ji)) call g_safe_alloc%deallocate('green.g01ji', this%g01ji)
      if (allocated(this%gx0ij)) call g_safe_alloc%deallocate('green.gx0ij', this%gx0ij)
      if (allocated(this%gy0ij)) call g_safe_alloc%deallocate('green.gy0ij', this%gy0ij)
      if (allocated(this%gz0ij)) call g_safe_alloc%deallocate('green.gz0ij', this%gz0ij)
      if (allocated(this%gx1ij)) call g_safe_alloc%deallocate('green.gx1ij', this%gx1ij)
      if (allocated(this%gy1ij)) call g_safe_alloc%deallocate('green.gy1ij', this%gy1ij)
      if (allocated(this%gz1ij)) call g_safe_alloc%deallocate('green.gz1ij', this%gz1ij)
      if (allocated(this%gx0ji)) call g_safe_alloc%deallocate('green.gx0ji', this%gx0ji)
      if (allocated(this%gy0ji)) call g_safe_alloc%deallocate('green.gy0ji', this%gy0ji)
      if (allocated(this%gz0ji)) call g_safe_alloc%deallocate('green.gz0ji', this%gz0ji)
      if (allocated(this%gx1ji)) call g_safe_alloc%deallocate('green.gx1ji', this%gx1ji)
      if (allocated(this%gy1ji)) call g_safe_alloc%deallocate('green.gy1ji', this%gy1ji)
      if (allocated(this%gz1ji)) call g_safe_alloc%deallocate('green.gz1ji', this%gz1ji)
#else
      if (allocated(this%g0)) deallocate (this%g0)
      if (allocated(this%gij)) deallocate (this%gij)
      if (allocated(this%gji)) deallocate (this%gji)
      if (allocated(this%ginmag)) deallocate (this%ginmag)
      if (allocated(this%gjnmag)) deallocate (this%gjnmag)
      if (allocated(this%gix)) deallocate (this%gix)
      if (allocated(this%giz)) deallocate (this%giz)
      if (allocated(this%gjx)) deallocate (this%gjx)
      if (allocated(this%gjy)) deallocate (this%gjy)
      if (allocated(this%gjz)) deallocate (this%gjz)
      if (allocated(this%g00ij)) deallocate (this%g00ij)
      if (allocated(this%g01ij)) deallocate (this%g01ij)
      if (allocated(this%g00ji)) deallocate (this%g00ji)
      if (allocated(this%g01ji)) deallocate (this%g01ji)
      if (allocated(this%gx0ij)) deallocate (this%gx0ij)
      if (allocated(this%gy0ij)) deallocate (this%gy0ij)
      if (allocated(this%gz0ij)) deallocate (this%gz0ij)
      if (allocated(this%gx1ij)) deallocate (this%gx1ij)
      if (allocated(this%gy1ij)) deallocate (this%gy1ij)
      if (allocated(this%gz1ij)) deallocate (this%gz1ij)
      if (allocated(this%gx0ji)) deallocate (this%gx0ji)
      if (allocated(this%gy0ji)) deallocate (this%gy0ji)
      if (allocated(this%gz0ji)) deallocate (this%gz0ji)
      if (allocated(this%gx1ji)) deallocate (this%gx1ji)
      if (allocated(this%gy1ji)) deallocate (this%gy1ji)
      if (allocated(this%gz1ji)) deallocate (this%gz1ji)
#endif
   end subroutine destructor

   ! Member functions
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Reset all members to default
   !---------------------------------------------------------------------------
   subroutine restore_to_default(this)
      use mpi_mod
      class(green) :: this

#ifdef USE_SAFE_ALLOC
      if (this%lattice%njij == 0) then
         call g_safe_alloc%allocate('green.g0', this%g0, (/nb, nb, this%en%channels_ldos + 10, atoms_per_process/))
      else
         call g_safe_alloc%allocate('green.g0', this%g0, (/nb, nb, this%en%channels_ldos + 10, 4/))
      end if
      call g_safe_alloc%allocate('green.gij', this%gij, (/nb, nb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gji', this%gji, (/nb, nb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.ginmag', this%ginmag, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjnmag', this%gjnmag, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gix', this%gix, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.giy', this%giy, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.giz', this%giz, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjx', this%gjx, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjy', this%gjy, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjz', this%gjz, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.g00ij', this%g00ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.g00ji', this%g00ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.g00ij', this%g01ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.g00ji', this%g01ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gx0ij', this%gx0ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gy0ij', this%gy0ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gz0ij', this%gz0ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gx1ij', this%gx1ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gy1ij', this%gy1ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gz1ij', this%gz1ij, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gx0ji', this%gx0ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gy0ji', this%gy0ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gz0ji', this%gz0ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gx1ji', this%gx1ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gy1ji', this%gy1ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gz1ji', this%gz1ji, (/norb, norb, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gij_eta', this%gij_eta, (/64, nb, nb, atoms_per_process/))
      call g_safe_alloc%allocate('green.gji_eta', this%gji_eta, (/64, nb, nb, atoms_per_process/))
      call g_safe_alloc%allocate('green.ginmag_eta', this%ginmag_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjnmag_eta', this%gjnmag_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.gix_eta', this%gix_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.giy_eta', this%giy_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.giz_eta', this%giz_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjx_eta', this%gjx_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjy_eta', this%gjy_eta, (/64, norb, norb, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjz_eta', this%gjz_eta, (/64, norb, norb, atoms_per_process/))
#else
      if (this%lattice%njij == 0) then
         allocate (this%g0(nb, nb, this%en%channels_ldos + 10, this%lattice%nrec))
      else
         allocate (this%g0(nb, nb, this%en%channels_ldos + 10, 4))
      end if
      allocate (this%gij(nb, nb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gji(nb, nb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%ginmag(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gjnmag(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gix(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%giy(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%giz(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gjx(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gjy(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gjz(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%g00ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%g00ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%g01ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%g01ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gx0ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gy0ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gz0ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gx1ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gy1ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gz1ij(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gx0ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gy0ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gz0ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gx1ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gy1ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gz1ji(norb, norb, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gij_eta(64, nb, nb, atoms_per_process))
      allocate (this%gji_eta(64, nb, nb, atoms_per_process))
      allocate (this%ginmag_eta(64, norb, norb, atoms_per_process))
      allocate (this%gjnmag_eta(64, norb, norb, atoms_per_process))
      allocate (this%gix_eta(64, norb, norb, atoms_per_process))
      allocate (this%giy_eta(64, norb, norb, atoms_per_process))
      allocate (this%giz_eta(64, norb, norb, atoms_per_process))
      allocate (this%gjx_eta(64, norb, norb, atoms_per_process))
      allocate (this%gjy_eta(64, norb, norb, atoms_per_process))
      allocate (this%gjz_eta(64, norb, norb, atoms_per_process))
#endif

      this%g0(:, :, :, :) = (0.0d0, 0.0d0)
      this%gij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gji(:, :, :, :) = (0.0d0, 0.0d0)
      this%ginmag(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjnmag(:, :, :, :) = (0.0d0, 0.0d0)
      this%gix(:, :, :, :) = (0.0d0, 0.0d0)
      this%giy(:, :, :, :) = (0.0d0, 0.0d0)
      this%giz(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjx(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjy(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjz(:, :, :, :) = (0.0d0, 0.0d0)
      this%g00ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%g00ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%g00ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%g00ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gx0ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gy0ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gz0ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gx1ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gy1ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gz1ij(:, :, :, :) = (0.0d0, 0.0d0)
      this%gx0ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gy0ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gz0ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gx1ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gy1ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gz1ji(:, :, :, :) = (0.0d0, 0.0d0)
      this%gij_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%gji_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%ginmag_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjnmag_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%gix_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%giy_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%giz_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjx_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjy_eta(:, :, :, :) = (0.0d0, 0.0d0)
      this%gjz_eta(:, :, :, :) = (0.0d0, 0.0d0)
   end subroutine restore_to_default

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the Green functions as a function of the imaginary part eta
   !  for the inter-site calculation
   !---------------------------------------------------------------------------
   subroutine block_green_ij_eta(this, istart, eta, fermi_point, g_ef)
      implicit none
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), intent(in) :: eta
      complex(rp), dimension(nb, nb, 4), intent(inout) :: g_ef
      integer, intent(in) :: fermi_point
      !
      integer :: ll, n, nw, ldim, na
      real(rp), dimension(4) :: a_inf0, b_inf0
      real(rp), dimension(nb, nb, 4) :: a_inf, b_inf
      complex(rp), dimension(nb, nb, this%en%channels_ldos + 10, 4) :: dum_g_ef
      !
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      ldim = nb
      na = 4
      nw = 10*ll

      call this%recursion%get_terminf(this%recursion%a_b(:, :, :, istart), this%recursion%b2_b(:, :, :, istart), na, ll, ldim, nw, &
                                      a_inf, b_inf, a_inf0, b_inf0)

      do n = 1, NA
         call this%bgreen(dum_g_ef(:, :, :, n), n + istart - 1, fermi_point, 1, a_inf(:, :, n), b_inf(:, :, n), eta)
         g_ef(:, :, n) = dum_g_ef(:, :, fermi_point, n)
      end do

   end subroutine block_green_ij_eta

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the Green functions for the inter-site calculation
   !---------------------------------------------------------------------------
   subroutine block_green_ij(this, istart)
      implicit none
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      !
      complex(rp) :: eta
      integer :: ll, n, nw, ldim, na
      real(rp), dimension(4) :: a_inf0, b_inf0
      real(rp), dimension(nb, nb, 4) :: a_inf, b_inf
      !
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      ldim = nb
      na = 4
      eta = (0.0d0, 0.0d0)
      nw = 10*ll

      call this%recursion%get_terminf(this%recursion%a_b(:, :, :, istart), this%recursion%b2_b(:, :, :, istart), na, ll, ldim, nw, &
                                      a_inf, b_inf, a_inf0, b_inf0)

      do n = 1, NA
         ! Initializing Shift and Scaling for rigid band shift
         !Dfac_mat=(1.0d0,0.0d0)
         !Cshi_mat=(0.0d0,0.0d0)
         call this%bgreen(this%g0(:, :, :, n), n + istart - 1, 1, this%en%channels_ldos + 10, a_inf(:, :, n), b_inf(:, :, n), eta)

      end do
   end subroutine block_green_ij

   !---------------------------------------------------------------------------
   !> GPU drop-in for block_green_ij: the 4 intersite-pair combos are
   !> reconstructed on the device over all nv energies in one batched call
   !> (rsrec_block_dos, natoms=4). Terminator stays on the CPU.
   !---------------------------------------------------------------------------
   subroutine block_green_ij_gpu(this, istart)
      implicit none
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp) :: eta
      integer :: ll, n, nw, ldim, na, nv, i
      real(rp), dimension(4) :: a_inf0, b_inf0
      real(rp), dimension(nb, nb, 4) :: a_inf, b_inf
      real(rp), dimension(nb, 4) :: a_inf_d, b_inf_d
      type(rsgpu), save :: gpu_backend
      logical, save :: first_call = .true.

      ll = this%control%lld
      ldim = nb
      na = 4
      eta = (0.0d0, 0.0d0)
      nw = 10*ll
      nv = this%en%channels_ldos + 10

      call this%recursion%get_terminf(this%recursion%a_b(:, :, :, istart), this%recursion%b2_b(:, :, :, istart), &
                                      na, ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)
      do n = 1, na
         do i = 1, nb
            a_inf_d(i, n) = a_inf(i, i, n)
            b_inf_d(i, n) = b_inf(i, i, n)
         end do
      end do

      if (first_call) then
         call gpu_backend%init(1, nb, 1, 1, 1)
         call gpu_backend%set_precision(0)  ! 0=fp32 (fast), 1=fp64 (validation)
         first_call = .false.
      end if

      call gpu_backend%block_dos(this%recursion%a_b(:, :, :, istart:istart + na - 1), &
                                 this%recursion%b2_b(:, :, :, istart:istart + na - 1), &
                                 a_inf_d, b_inf_d, this%en%ene(1:nv), eta, &
                                 this%control%sym_term, this%g0(:, :, 1:nv, 1:na))
   end subroutine block_green_ij_gpu

   subroutine calculate_intersite_gf_twoindex(this)
      use mpi_mod
   use basis_mod, only: nb, norb, spin_off
      implicit none
      class(green), intent(inout) :: this
      integer :: i, j, k, l1, l2, k0, j0, ia_glob, ia

      this%G00ij = 0.0d0; this%G01ij = 0.0d0; this%G00ji = 0.0d0; this%G01ji = 0.0d0; 
      this%Gx1ij = 0.0d0; this%Gy1ij = 0.0d0; this%Gz1ij = 0.0d0; this%Gx0ij = 0.0d0; this%Gy0ij = 0.0d0; this%Gz0ij = 0.0d0; 
      this%Gx1ji = 0.0d0; this%Gy1ji = 0.0d0; this%Gz1ji = 0.0d0; this%Gx0ji = 0.0d0; this%Gy0ji = 0.0d0; this%Gz0ji = 0.0d0

      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)
         do j = 1, norb
            do k = 1, norb
               l1 = int((k - 0.9)**0.5)
               l2 = int((j - 0.9)**0.5)
               k0 = l1*(l1 + 1) + 1
               j0 = l2*(l2 + 1) + 1
               this%G00ij(k, j, :, IA) = this%G00ij(k, j, :, IA) + 0.5d0*(this%Ginmag(k, j, :, IA) + (-1)**(k + j)*this%Gjnmag(2*j0 - j, 2*k0 - k, :, IA))
               this%G01ij(k, j, :, IA) = this%G01ij(k, j, :, IA) + 0.5d0*(this%Ginmag(k, j, :, IA) - (-1)**(k + j)*this%Gjnmag(2*j0 - j, 2*k0 - k, :, IA))
               this%G00ji(k, j, :, IA) = this%G00ji(k, j, :, IA) + 0.5d0*(this%Gjnmag(k, j, :, IA) + (-1)**(k + j)*this%Ginmag(2*j0 - j, 2*k0 - k, :, IA))
               this%G01ji(k, j, :, IA) = this%G01ji(k, j, :, IA) + 0.5d0*(this%Gjnmag(k, j, :, IA) - (-1)**(k + j)*this%Ginmag(2*j0 - j, 2*k0 - k, :, IA))
               this%Gx1ij(k, j, :, IA) = this%Gx1ij(k, j, :, IA) + 0.5d0*(this%Gix(k, j, :, IA) - (-1)**(k + j)*this%Gjx(2*j0 - j, 2*k0 - k, :, IA))
               this%Gy1ij(k, j, :, IA) = this%Gy1ij(k, j, :, IA) + 0.5d0*(this%Giy(k, j, :, IA) - (-1)**(k + j)*this%Gjy(2*j0 - j, 2*k0 - k, :, IA))
               this%Gz1ij(k, j, :, IA) = this%Gz1ij(k, j, :, IA) + 0.5d0*(this%Giz(k, j, :, IA) - (-1)**(k + j)*this%Gjz(2*j0 - j, 2*k0 - k, :, IA))
               this%Gx0ij(k, j, :, IA) = this%Gx0ij(k, j, :, IA) + 0.5d0*(this%Gix(k, j, :, IA) + (-1)**(k + j)*this%Gjx(2*j0 - j, 2*k0 - k, :, IA))
               this%Gy0ij(k, j, :, IA) = this%Gy0ij(k, j, :, IA) + 0.5d0*(this%Giy(k, j, :, IA) + (-1)**(k + j)*this%Gjy(2*j0 - j, 2*k0 - k, :, IA))
               this%Gz0ij(k, j, :, IA) = this%Gz0ij(k, j, :, IA) + 0.5d0*(this%Giz(k, j, :, IA) + (-1)**(k + j)*this%Gjz(2*j0 - j, 2*k0 - k, :, IA))
               this%Gx1ji(k, j, :, IA) = this%Gx1ji(k, j, :, IA) + 0.5d0*(this%Gjx(k, j, :, IA) - (-1)**(k + j)*this%Gix(2*j0 - j, 2*k0 - k, :, IA))
               this%Gy1ji(k, j, :, IA) = this%Gy1ji(k, j, :, IA) + 0.5d0*(this%Gjy(k, j, :, IA) - (-1)**(k + j)*this%Giy(2*j0 - j, 2*k0 - k, :, IA))
               this%Gz1ji(k, j, :, IA) = this%Gz1ji(k, j, :, IA) + 0.5d0*(this%Gjz(k, j, :, IA) - (-1)**(k + j)*this%Giz(2*j0 - j, 2*k0 - k, :, IA))
               this%Gx0ji(k, j, :, IA) = this%Gx0ji(k, j, :, IA) + 0.5d0*(this%Gjx(k, j, :, IA) + (-1)**(k + j)*this%Gix(2*j0 - j, 2*k0 - k, :, IA))
               this%Gy0ji(k, j, :, IA) = this%Gy0ji(k, j, :, IA) + 0.5d0*(this%Gjy(k, j, :, IA) + (-1)**(k + j)*this%Giy(2*j0 - j, 2*k0 - k, :, IA))
               this%Gz0ji(k, j, :, IA) = this%Gz0ji(k, j, :, IA) + 0.5d0*(this%Gjz(k, j, :, IA) + (-1)**(k + j)*this%Giz(2*j0 - j, 2*k0 - k, :, IA))
            end do
         end do
      end do
   end subroutine calculate_intersite_gf_twoindex

   subroutine calculate_intersite_gf(this)
      use mpi_mod
      implicit none
      class(green), intent(inout) :: this
      integer :: ia, ja_temp, j, i, ia_glob

      this%gij = 0.0d0; this%gji = 0.0d0; this%ginmag = 0.0d0; this%giz = 0.0d0; this%giy = 0.0d0; this%gix = 0.0d0
      this%gjnmag = 0.0d0; this%gjx = 0.0d0; this%gjy = 0.0d0; this%gjz = 0.0d0

      if (this%control%recur == 'block') call this%recursion%zsqr()

      !do ia=1,this%lattice%njij
      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)
         ja_temp = (ia - 1)*4 + 1
         select case (this%control%recur)
         case ('block')
            if (this%control%gpu_plugin) then
               call this%block_green_ij_gpu(ja_temp)
            else
               call this%block_green_ij(ja_temp)
            end if
         case ('chebyshev')
            if (this%control%gpu_plugin) then
               call this%chebyshev_green_ij_gpu(ja_temp)
            else
               call this%chebyshev_green_ij(ja_temp)
            end if
         end select
         if (this%lattice%ijpair(ia, 1) .eq. this%lattice%ijpair(ia, 2)) then
            this%gij(:, :, :, ia) = this%g0(:, :, :, 1)
            this%gji(:, :, :, ia) = this%g0(:, :, :, 1)
         else
            this%gij(:, :, :, ia) = this%g0(:, :, :, 1) - this%g0(:, :, :, 2) + (1.0d0/i_unit*this%g0(:, :, :, 3) - 1.0d0/i_unit*this%g0(:, :, :, 4))
            this%gji(:, :, :, ia) = this%g0(:, :, :, 1) - this%g0(:, :, :, 2) - (1.0d0/i_unit*this%g0(:, :, :, 3) - 1.0d0/i_unit*this%g0(:, :, :, 4))
            this%gij(:, :, :, ia) = this%gij(:, :, :, ia)*0.5d0
            this%gji(:, :, :, ia) = this%gji(:, :, :, ia)*0.5d0
         end if
         do i = 1, norb
            do j = 1, norb
               this%Ginmag(j, i, :, iA) = this%Ginmag(j, i, :, iA) + (this%gij(j, i, :, ia) + this%gij(j +spin_off, i +spin_off, :, ia))*0.5d0
               this%Giz(j, i, :, iA) = this%Giz(j, i, :, iA) + 0.5d0*(this%gij(j, i, :, ia) - this%gij(j +spin_off, i +spin_off, :, ia))        !+Ginmag(k,j,i,iIA)*0.5d0
               this%Giy(j, i, :, iA) = this%Giy(j, i, :, iA) + 0.5d0*(i_unit*this%gij(j, i +spin_off, :, ia) - i_unit*this%gij(j +spin_off, i, :, ia))  !+Ginmag(k,j,i,iIA)*0.5d0
               this%Gix(j, i, :, iA) = this%Gix(j, i, :, iA) + 0.5d0*(this%gij(j, i +spin_off, :, ia) + this%gij(j +spin_off, i, :, ia))        !+Ginmag(k,j,i,iIA)*0.5d0
               !
               this%Gjnmag(j, i, :, iA) = this%Gjnmag(j, i, :, iA) + (this%gji(j, i, :, ia) + this%gji(j +spin_off, i +spin_off, :, ia))*0.5d0
               this%Gjz(j, i, :, iA) = this%Gjz(j, i, :, iA) + 0.5d0*(this%gji(j, i, :, ia) - this%gji(j +spin_off, i +spin_off, :, ia))         !+Gjnmag(k,j,i,iIA)*0.5d0
               this%Gjy(j, i, :, iA) = this%Gjy(j, i, :, iA) + 0.5d0*(i_unit*this%gji(j, i +spin_off, :, ia) - i_unit*this%gji(j +spin_off, i, :, ia))   !+Gjnmag(k,j,i,iIA)*0.5d0
               this%Gjx(j, i, :, iA) = this%Gjx(j, i, :, iA) + 0.5d0*(this%gji(j, i +spin_off, :, ia) + this%gji(j +spin_off, i, :, ia))         !+Gjnmag(k,j,i,iIA)*0.5d0
            end do
         end do
      end do
   end subroutine calculate_intersite_gf

   subroutine calculate_intersite_gf_eta(this)
      use mpi_mod
      implicit none
      class(green), intent(inout) :: this
      integer :: ia, ja_temp, j, i, fermi_point
      complex(rp), dimension(nb, nb, 4) :: g0_ef
      complex(rp) :: eta
      complex(rp), dimension(64, nb, nb, 4) :: y
      real(rp) :: res, t
      real(rp), dimension(64) :: x, w

      integer :: ia_glob

      ! GPU path: the whole per-pair/per-eta GF loop (block continued fraction
      ! or Chebyshev moment contraction) is reconstructed on the device in one
      ! batched call over all local pairs.
      if (this%control%gpu_plugin .and. &
          (this%control%recur == 'block' .or. this%control%recur == 'chebyshev')) then
         call this%calculate_intersite_gf_eta_gpu()
         return
      end if

      ! Find the Gauss Legendre roots and weights
      call gauss_legendre(64, 0.0_rp, 1.0_rp, x, w)

      call this%en%e_mesh()

      do i = 1, this%en%channels_ldos + 10
         if ((this%en%ene(i) - this%en%fermi) .le. 0.000001d0) fermi_point = i
      end do

      write (*, *) this%en%fermi, fermi_point
      this%gij_eta = 0.0d0; this%gji_eta = 0.0d0; this%ginmag_eta = 0.0d0; this%giz_eta = 0.0d0; this%giy_eta = 0.0d0; this%gix_eta = 0.0d0
      this%gjnmag_eta = 0.0d0; this%gjx_eta = 0.0d0; this%gjy_eta = 0.0d0; this%gjz_eta = 0.0d0

      if (this%control%recur == 'block') call this%recursion%zsqr()

      !do ia=1,this%lattice%njij
      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)
         ja_temp = (ia - 1)*4 + 1
         y = (0.0_rp, 0.0_rp)
         do i = 1, 64
            eta = (0.0_rp, 0.0_rp)
            g0_ef = (0.0_rp, 0.0_rp)
            res = (1 - x(i))/x(i)
            eta = cmplx(0.0_rp, res)
            select case (this%control%recur)
            case ('block')
               call this%block_green_ij_eta(ja_temp, eta, fermi_point, g0_ef)
            case ('chebyshev')
               call this%chebyshev_green_ij_eta(ja_temp, eta, fermi_point, g0_ef)
            end select
            y(i, :, :, :) = g0_ef(:, :, :)
         end do
         this%gij_eta(:, :, :, ia) = y(:, :, :, 1) - y(:, :, :, 2) + (1.0d0/i_unit*y(:, :, :, 3) - 1.0d0/i_unit*y(:, :, :, 4))
         this%gji_eta(:, :, :, ia) = y(:, :, :, 1) - y(:, :, :, 2) - (1.0d0/i_unit*y(:, :, :, 3) - 1.0d0/i_unit*y(:, :, :, 4))

         this%gij_eta(:, :, :, ia) = this%gij_eta(:, :, :, ia)*0.5d0
         this%gji_eta(:, :, :, ia) = this%gji_eta(:, :, :, ia)*0.5d0
         do i = 1, norb
            do j = 1, norb
               this%Ginmag_eta(:, j, i, iA) = this%Ginmag_eta(:, j, i, iA) + (this%gij_eta(:, j, i, ia) + this%gij_eta(:, j +spin_off, i +spin_off, ia))*0.5d0
               this%Giz_eta(:, j, i, iA) = this%Giz_eta(:, j, i, iA) + 0.5d0*(this%gij_eta(:, j, i, ia) - this%gij_eta(:, j +spin_off, i +spin_off, ia))        !+Ginmag(k,j,i,iIA)*0.5d0
               this%Giy_eta(:, j, i, iA) = this%Giy_eta(:, j, i, iA) + 0.5d0*(i_unit*this%gij_eta(:, j, i +spin_off, ia) - i_unit*this%gij_eta(:, j +spin_off, i, ia))  !+Ginmag(k,j,i,iIA)*0.5d0
               this%Gix_eta(:, j, i, iA) = this%Gix_eta(:, j, i, iA) + 0.5d0*(this%gij_eta(:, j, i +spin_off, ia) + this%gij_eta(:, j +spin_off, i, ia))        !+Ginmag(k,j,i,iIA)*0.5d0
               !
               this%Gjnmag_eta(:, j, i, iA) = this%Gjnmag_eta(:, j, i, iA) + (this%gji_eta(:, j, i, ia) + this%gji_eta(:, j +spin_off, i +spin_off, ia))*0.5d0
               this%Gjz_eta(:, j, i, iA) = this%Gjz_eta(:, j, i, iA) + 0.5d0*(this%gji_eta(:, j, i, ia) - this%gji_eta(:, j +spin_off, i +spin_off, ia))         !+Gjnmag(k,j,i,iIA)*0.5d0
               this%Gjy_eta(:, j, i, iA) = this%Gjy_eta(:, j, i, iA) + 0.5d0*(i_unit*this%gji_eta(:, j, i +spin_off, ia) - i_unit*this%gji_eta(:, j +spin_off, i, ia))   !+Gjnmag(k,j,i,iIA)*0.5d0
               this%Gjx_eta(:, j, i, iA) = this%Gjx_eta(:, j, i, iA) + 0.5d0*(this%gji_eta(:, j, i +spin_off, ia) + this%gji_eta(:, j +spin_off, i, ia))         !+Gjnmag(k,j,i,iIA)*0.5d0
            end do
         end do
      end do
   end subroutine calculate_intersite_gf_eta

   !---------------------------------------------------------------------------
   !> GPU port of calculate_intersite_gf_eta: the per-pair x per-eta bgreen loop
   !> (nij x 64 x 4 continued fractions at the Fermi energy) is reconstructed on
   !> the device in ONE batched call (rsrec_block_gf_eta), over all local pairs.
   !> The terminator (get_terminf) stays on the CPU; only diagonals go down.
   !---------------------------------------------------------------------------
   subroutine calculate_intersite_gf_eta_gpu(this)
      use mpi_mod
      implicit none
      class(green), intent(inout) :: this
      integer :: ia, ja_temp, j, i, k, c, fermi_point, ia_glob, n_pairs, natoms, ll, ldim, nw
      complex(rp), dimension(64, nb, nb, 4) :: y
      real(rp), dimension(64) :: x, w
      complex(rp), dimension(64) :: eta_list
      real(rp), allocatable :: a_inf(:, :, :), b_inf(:, :, :), a_inf0(:), b_inf0(:)
      real(rp), allocatable :: a_inf_d(:, :), b_inf_d(:, :)
      complex(rp), allocatable :: g0_ef_all(:, :, :, :)
      ! Static GPU backend (one per program lifetime), as in chebyshev_green_gpu
      type(rsgpu), save :: gpu_backend
      logical, save :: first_call = .true.

      ll = this%control%lld
      ldim = nb
      nw = 10*ll

      call gauss_legendre(64, 0.0_rp, 1.0_rp, x, w)
      call this%en%e_mesh()
      fermi_point = 1
      do i = 1, this%en%channels_ldos + 10
         if ((this%en%ene(i) - this%en%fermi) .le. 0.000001d0) fermi_point = i
      end do

      this%gij_eta = 0.0d0; this%gji_eta = 0.0d0; this%ginmag_eta = 0.0d0; this%giz_eta = 0.0d0; this%giy_eta = 0.0d0; this%gix_eta = 0.0d0
      this%gjnmag_eta = 0.0d0; this%gjx_eta = 0.0d0; this%gjy_eta = 0.0d0; this%gjz_eta = 0.0d0

      if (this%control%recur == 'block') call this%recursion%zsqr()

      ! Gauss-Legendre eta contour (same map as the CPU path: eta = i*(1-x)/x)
      do k = 1, 64
         eta_list(k) = cmplx(0.0_rp, (1.0_rp - x(k))/x(k), kind=rp)
      end do

      ! All local pairs share the contiguous coefficient layout (4 combos/pair).
      n_pairs = end_atom - start_atom + 1
      natoms = n_pairs*4
      allocate (g0_ef_all(nb, nb, 64, natoms))

      if (first_call) then
         call gpu_backend%init(1, nb, 1, 1, 1)
         call gpu_backend%set_precision(0)  ! 0=fp32 (fast), 1=fp64 (validation)
         first_call = .false.
      end if

      ! One device call: Gij at the Fermi energy for all (pair x combo) atoms
      ! and all 64 contour etas -> g0_ef_all(nb,nb,64,natoms).
      select case (this%control%recur)
      case ('block')
         allocate (a_inf(nb, nb, natoms), b_inf(nb, nb, natoms), a_inf0(natoms), b_inf0(natoms))
         allocate (a_inf_d(nb, natoms), b_inf_d(nb, natoms))
         call this%recursion%get_terminf(this%recursion%a_b, this%recursion%b2_b, &
                                         natoms, ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)
         do k = 1, natoms
            do i = 1, nb
               a_inf_d(i, k) = a_inf(i, i, k)
               b_inf_d(i, k) = b_inf(i, i, k)
            end do
         end do
         call gpu_backend%block_gf_eta(this%recursion%a_b(:, :, :, 1:natoms), &
                                       this%recursion%b2_b(:, :, :, 1:natoms), &
                                       a_inf_d, b_inf_d, this%en%ene(fermi_point), &
                                       eta_list, this%control%sym_term, g0_ef_all)
         deallocate (a_inf, b_inf, a_inf0, b_inf0, a_inf_d, b_inf_d)
      case ('chebyshev')
         block
            real(rp) :: a, b, emin_win, emax_win
            call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
            a = (emax_win - emin_win)/(2 - 0.3_rp)
            b = (emax_win + emin_win)/2.0_rp
            call gpu_backend%chebyshev_gf_eta(this%recursion%mu_n(:, :, :, 1:natoms), &
                                              this%en%ene(fermi_point), eta_list, a, b, g0_ef_all)
         end block
      end select

      ! Recombine per pair (identical algebra to the CPU loop).
      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)
         ja_temp = (ia - 1)*4 + 1
         do c = 1, 4
            do k = 1, 64
               y(k, :, :, c) = g0_ef_all(:, :, k, ja_temp + c - 1)
            end do
         end do
         this%gij_eta(:, :, :, ia) = (y(:, :, :, 1) - y(:, :, :, 2) + (1.0d0/i_unit*y(:, :, :, 3) - 1.0d0/i_unit*y(:, :, :, 4)))*0.5d0
         this%gji_eta(:, :, :, ia) = (y(:, :, :, 1) - y(:, :, :, 2) - (1.0d0/i_unit*y(:, :, :, 3) - 1.0d0/i_unit*y(:, :, :, 4)))*0.5d0
         do i = 1, norb
            do j = 1, norb
               this%Ginmag_eta(:, j, i, iA) = this%Ginmag_eta(:, j, i, iA) + (this%gij_eta(:, j, i, ia) + this%gij_eta(:, j +spin_off, i +spin_off, ia))*0.5d0
               this%Giz_eta(:, j, i, iA) = this%Giz_eta(:, j, i, iA) + 0.5d0*(this%gij_eta(:, j, i, ia) - this%gij_eta(:, j +spin_off, i +spin_off, ia))
               this%Giy_eta(:, j, i, iA) = this%Giy_eta(:, j, i, iA) + 0.5d0*(i_unit*this%gij_eta(:, j, i +spin_off, ia) - i_unit*this%gij_eta(:, j +spin_off, i, ia))
               this%Gix_eta(:, j, i, iA) = this%Gix_eta(:, j, i, iA) + 0.5d0*(this%gij_eta(:, j, i +spin_off, ia) + this%gij_eta(:, j +spin_off, i, ia))
               !
               this%Gjnmag_eta(:, j, i, iA) = this%Gjnmag_eta(:, j, i, iA) + (this%gji_eta(:, j, i, ia) + this%gji_eta(:, j +spin_off, i +spin_off, ia))*0.5d0
               this%Gjz_eta(:, j, i, iA) = this%Gjz_eta(:, j, i, iA) + 0.5d0*(this%gji_eta(:, j, i, ia) - this%gji_eta(:, j +spin_off, i +spin_off, ia))
               this%Gjy_eta(:, j, i, iA) = this%Gjy_eta(:, j, i, iA) + 0.5d0*(i_unit*this%gji_eta(:, j, i +spin_off, ia) - i_unit*this%gji_eta(:, j +spin_off, i, ia))
               this%Gjx_eta(:, j, i, iA) = this%Gjx_eta(:, j, i, iA) + 0.5d0*(this%gji_eta(:, j, i +spin_off, ia) + this%gji_eta(:, j +spin_off, i, ia))
            end do
         end do
      end do

      deallocate (g0_ef_all)
   end subroutine calculate_intersite_gf_eta_gpu

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the Greens function as a function of the imaginary energy increment
   !  eta. Used for the Gauss-Legendre integration
   !---------------------------------------------------------------------------
   subroutine block_green_eta(this, eta, fermi_point, g_ef)
      use mpi_mod
      implicit none
      class(green), intent(inout) :: this
      complex(rp), intent(in) :: eta
      integer, intent(in) :: fermi_point
      complex(rp), dimension(nb, nb, atoms_per_process), intent(inout) :: g_ef
      !
      !
      integer :: nw
      integer :: ll
      integer :: ldim
      real(rp), dimension(this%lattice%nrec) :: a_inf0, b_inf0
      real(rp), dimension(nb, nb, this%lattice%nrec) :: a_inf, b_inf
      complex(rp), dimension(nb, nb, this%en%channels_ldos + 10, atoms_per_process) :: dum_g_ef
      integer :: n, n_glob
      !
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      ldim = nb

      ! Unchanged code
      nw = 10*ll

      ! Calculate terminator coefficients
      call this%recursion%get_terminf(this%recursion%a_b, this%recursion%b2_b, atoms_per_process, &
                                      ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)

      do n_glob = start_atom, end_atom
         n = g2l_map(n_glob)

         call this%bgreen(dum_g_ef(:, :, :, n), n, fermi_point, 1, a_inf(:, :, n), b_inf(:, :, n), eta)
         g_ef(:, :, n) = dum_g_ef(:, :, fermi_point, n)

      end do

   end subroutine block_green_eta

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the Greens function
   !---------------------------------------------------------------------------
   subroutine block_green(this)
      use mpi_mod
      implicit none
      class(green), intent(inout) :: this
      !
      integer :: nw
      integer :: ll
      integer :: ldim
      complex(rp) :: eta
      real(rp), dimension(this%lattice%nrec) :: a_inf0, b_inf0
      real(rp), dimension(nb, nb, this%lattice%nrec) :: a_inf, b_inf
      integer :: n, n_glob
      !
      ! GPU path: reconstruct the block continued fraction on the device.
      if (this%control%gpu_plugin) then
         call this%block_green_gpu()
         return
      end if
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      ldim = nb
      !na = this%lattice%nrec
      eta = (0.0d0, 0.0d0)

      ! Unchanged code
      nw = 10*ll

      ! Calculate terminator coefficients
      call this%recursion%get_terminf(this%recursion%a_b, this%recursion%b2_b, atoms_per_process, &
                                      ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)

      ! print *, "AB check:", shape(this%g0)
      do n_glob = start_atom, end_atom
         n = g2l_map(n_glob)

         call this%bgreen(this%g0(:, :, :, n), n, 1, this%en%channels_ldos + 10, a_inf(:, :, n), b_inf(:, :, n), eta)
      end do

   end subroutine block_green

   !---------------------------------------------------------------------------
   !> GPU port of block_green: the block (Haydock) continued-fraction Green's
   !> function is reconstructed on the device (rsrec_block_dos), replacing the
   !> per-atom/per-energy bgreen loop. The terminator (get_terminf) stays on the
   !> CPU; only its diagonals are passed down. Falls back transparently to the
   !> legacy path if the GPU plugin is unavailable (handled inside the backend).
   !---------------------------------------------------------------------------
   subroutine block_green_gpu(this)
      use mpi_mod
      implicit none
      class(green), intent(inout) :: this
      !
      integer :: nw, ll, ldim, nv, n_loc, n, n_glob, k, i
      complex(rp) :: eta
      real(rp), dimension(this%lattice%nrec) :: a_inf0, b_inf0
      real(rp), dimension(nb, nb, this%lattice%nrec) :: a_inf, b_inf
      complex(rp), allocatable :: ab_loc(:, :, :, :), b2_loc(:, :, :, :), g0_loc(:, :, :, :)
      real(rp), allocatable :: ainf_d(:, :), binf_d(:, :)
      ! Static GPU backend (one per program lifetime), as in chebyshev_green_gpu
      type(rsgpu), save :: gpu_backend
      logical, save :: first_call = .true.

      ll = this%control%lld
      ldim = nb
      eta = (0.0d0, 0.0d0)
      nw = 10*ll
      nv = this%en%channels_ldos + 10

      ! Terminator coefficients on the CPU (a_inf/b_inf), as in block_green.
      call this%recursion%get_terminf(this%recursion%a_b, this%recursion%b2_b, &
                                      atoms_per_process, ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)

      ! Gather this rank's local atoms into contiguous batch buffers (mirrors the
      ! n = g2l_map(n_glob) indexing of block_green's bgreen loop).
      n_loc = end_atom - start_atom + 1
      allocate (ab_loc(nb, nb, ll, n_loc), b2_loc(nb, nb, ll, n_loc), g0_loc(nb, nb, nv, n_loc))
      allocate (ainf_d(nb, n_loc), binf_d(nb, n_loc))
      do n_glob = start_atom, end_atom
         n = g2l_map(n_glob)
         k = n_glob - start_atom + 1
         ab_loc(:, :, :, k) = this%recursion%a_b(:, :, :, n)
         b2_loc(:, :, :, k) = this%recursion%b2_b(:, :, :, n)
         do i = 1, nb
            ainf_d(i, k) = a_inf(i, i, n)
            binf_d(i, k) = b_inf(i, i, n)
         end do
      end do

      if (first_call) then
         call gpu_backend%init(1, nb, 1, 1, 1)
         call gpu_backend%set_precision(0)  ! 0=fp32 (fast), 1=fp64 (validation)
         first_call = .false.
      end if

      call gpu_backend%block_dos(ab_loc, b2_loc, ainf_d, binf_d, &
                                 this%en%ene(1:nv), eta, this%control%sym_term, g0_loc)

      ! Scatter back to the local g0 slots.
      do n_glob = start_atom, end_atom
         n = g2l_map(n_glob)
         k = n_glob - start_atom + 1
         this%g0(:, :, 1:nv, n) = g0_loc(:, :, :, k)
      end do

      deallocate (ab_loc, b2_loc, g0_loc, ainf_d, binf_d)
   end subroutine block_green_gpu

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the Greens function
   !---------------------------------------------------------------------------
   subroutine sgreen(this)
      use mpi_mod
      !this, bdos, g0, e, nv, ldim, ll, na, nmdir, mom, ef)
      class(green) :: this
      ! Input
      !integer, intent(in) :: nv, ldim, ll, na, nmdir
      !real(rp), intent(in) :: ef
      ! Output
      !real(rp), dimension(nv), intent(inout) :: e
      !real(rp), dimension(nv, ldim, ldim, na), intent(out) :: bdos
      !complex(rp), dimension(nv, ldim, ldim, na), intent(out) :: g0
      !real(rp), dimension(na, 3), intent(inout) :: mom
      ! Local variables
      integer :: ia, ja, mdir, nw, ll_t, ie, j, i
      real(rp), dimension(this%en%channels_ldos + 10, this%lattice%ntype) :: dx, dy, dz
      real(rp), dimension(nb, this%en%channels_ldos + 10) :: doso
      real(rp), dimension(nb, this%en%channels_ldos + 10) :: dmag, dnmag
      complex(rp) :: dfac, sfac, impi
      complex(rp), dimension(4) :: gspinor
      complex(rp), dimension(3) :: gmask, lmask
      complex(rp), dimension(2, 3) :: gfac
      integer, dimension(4, 3) :: goff
      real(rp), dimension(this%control%lld, this%control%lld)  :: Sm
      real(rp), dimension(nb)  :: q_int
      real(rp) :: e_start, e_stop

      integer :: ia_glob

      impi = (pi, 0.0d0)
      gfac(1, 1) = 1.0d0; gfac(1, 2) = -i_unit; gfac(1, 3) = 1.0d0
      gfac(2, 1) = 1.0d0; gfac(2, 2) = i_unit; gfac(2, 3) = -1.0d0
      gmask(1) = 0.0d0; gmask(2) = 0.0d0; gmask(3) = 1.0d0
      !   if(lrot) then
      !      lmask(1)=0.0d0;lmask(2)=0.0d0;lmask(3)=1.0d0
      !   else
      lmask(1) = 1.0d0/3.0d0; lmask(2) = 1.0d0/3.0d0; lmask(3) = 1.0d0/3.0d0
      !   end if

      goff(1, 1) = 0; goff(2, 1) = 9; goff(1, 2) = 0; goff(2, 2) = 9; goff(1, 3) = 0; goff(2, 3) = 0
      goff(3, 1) = 9; goff(4, 1) = 0; goff(3, 2) = 9; goff(4, 2) = 0; goff(3, 3) = 9; goff(4, 3) = 9
      dfac = i_unit*impi/(2.0d0, 0.0d0)
      dx = 0; dy = 0; dz = 0
      nw = 10*this%lattice%ntype*this%control%lld
      ll_t = this%control%lld

      doso = 0.0d0
      this%g0 = 0.0d0
      !do ia = 1, this%lattice%nrec
      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)
         do mdir = 1, this%control%nmdir
            doso = 0.0d0
            call this%dos%density(doso, ia, mdir)
            if (this%control%nmdir == 1) then
               do ie = 1, this%en%channels_ldos + 10
                  do j = 1, nb
                     this%g0(j, j, ie, ia) = -i_unit*doso(j, ie)*impi
                  end do
                  write (300 + ia, *) this%en%ene(ie), sum(doso(1:norb, ie)), sum(doso(norb+1:nb, ie))
               end do
            else
               do ie = 1, this%en%channels_ldos + 10
                  do j = 1, norb
                     ! Charge, from main direction.. (not z-component)
                     this%g0(j, j, ie, ia) = this%g0(j, j, ie, ia) - (doso(j, ie) + doso(j +spin_off, ie))*dfac*lmask(mdir)!*1.0d0 /3.0d0 !mom(ja, mdir)**2
                     this%g0(j +spin_off, j +spin_off, ie, ia) = this%g0(j +spin_off, j +spin_off, ie, ia) - (doso(j, ie) + doso(j +spin_off, ie))*dfac*lmask(mdir)!*1.0d0 /3.0d0 !mom(ja, mdir)**2
                     ! Spin dependent part
                     this%g0(j + goff(1, mdir), j + goff(2, mdir), ie, ia) = &
                        this%g0(j + goff(1, mdir), j + goff(2, mdir), ie, ia) - (doso(j, ie) - doso(j +spin_off, ie))*gfac(1, mdir)*dfac

                     this%g0(j + goff(3, mdir), j + goff(4, mdir), ie, ia) = &
                        this%g0(j + goff(3, mdir), j + goff(4, mdir), ie, ia) - (doso(j, ie) - doso(j +spin_off, ie))*gfac(2, mdir)*dfac
                  end do
               end do
            end if
         end do
      end do
   end subroutine sgreen

   ! Routines that clones the g_eta to g. For use in exchange calculations
   subroutine gij_eta_to_gij(this)
      use mpi_mod
      !
      class(green) :: this
      integer :: echan, ijij

      if (allocated(this%Ginmag)) deallocate (this%Ginmag)
      if (allocated(this%Gjnmag)) deallocate (this%Gjnmag)
      if (allocated(this%Gix)) deallocate (this%Gix)
      if (allocated(this%Gjx)) deallocate (this%Gjx)
      if (allocated(this%Giy)) deallocate (this%Giy)
      if (allocated(this%Gjy)) deallocate (this%Gjy)
      if (allocated(this%Giz)) deallocate (this%Giz)
      if (allocated(this%Gjz)) deallocate (this%Gjz)

      allocate (this%Ginmag(norb, norb, 64, atoms_per_process))
      allocate (this%Gjnmag(norb, norb, 64, atoms_per_process))
      allocate (this%Gix(norb, norb, 64, atoms_per_process))
      allocate (this%Gjx(norb, norb, 64, atoms_per_process))
      allocate (this%Giy(norb, norb, 64, atoms_per_process))
      allocate (this%Gjy(norb, norb, 64, atoms_per_process))
      allocate (this%Giz(norb, norb, 64, atoms_per_process))
      allocate (this%Gjz(norb, norb, 64, atoms_per_process))

      do ijij = 1, atoms_per_process
         do echan = 1, 64
            this%Ginmag(:, :, echan, ijij) = this%Ginmag_eta(echan, :, :, ijij)
            this%Gjnmag(:, :, echan, ijij) = this%Gjnmag_eta(echan, :, :, ijij)
            this%Gix(:, :, echan, ijij) = this%Gix_eta(echan, :, :, ijij)
            this%Gjx(:, :, echan, ijij) = this%Gjx_eta(echan, :, :, ijij)
            this%Giy(:, :, echan, ijij) = this%Giy_eta(echan, :, :, ijij)
            this%Gjy(:, :, echan, ijij) = this%Gjy_eta(echan, :, :, ijij)
            this%Giz(:, :, echan, ijij) = this%Giz_eta(echan, :, :, ijij)
            this%Gjz(:, :, echan, ijij) = this%Gjz_eta(echan, :, :, ijij)
         end do
      end do

      return

   end subroutine gij_eta_to_gij
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Given an input Gij (physical site-resolved GF in the LMTO), the function
   !> returns the correspondent auxiliary GFs, considering here the
   !> orthogonal representation. Takes into account all njij pairs directly.
   !> It uses the fact that the phyisical GF´s are invariant with respect to
   !> screening constants (Turek´s book, page 72-73).
   !> Implemented by Ivan Miranda on 29.09.2023
   !---------------------------------------------------------------------------
   subroutine auxiliary_gij(this, green_ij, aux_gij, atom_i, atom_j)
      !
      class(green) :: this
      !
      integer, intent(in) :: atom_i, atom_j ! Input atoms for the Gij
      complex(rp), intent(in) :: green_ij(:, :, :) ! Input Gij (physical site-resolved GFs)
      complex(rp), intent(inout) :: aux_gij(:, :, :) ! Output auxiliary gij matrix (same size as Gij)
      !
      ! Local Variables
      integer :: nv, l, m, mls, s, lmaxi, lmaxj
      complex(rp), allocatable, dimension(:, :) :: cdelta_i, cdelta_j

      allocate (cdelta_i(size(green_ij, 1), size(green_ij, 2)))
      allocate (cdelta_j(size(green_ij, 1), size(green_ij, 2)))

      ! Set all matrices initially to zero
      aux_gij(:, :, :) = czero
      ! Calculate the lmax values of the input atoms
      lmaxi = this%symbolic_atom(this%lattice%iz(atom_i))%potential%lmax
      lmaxj = this%symbolic_atom(this%lattice%iz(atom_j))%potential%lmax

      do nv = 1, size(this%en%ene) ! Energy channel

         cdelta_i(:, :) = czero ! As zero in the beginning
         cdelta_j(:, :) = czero ! As zero in the beginning

         do s = 1, 2
            do l = 0, lmaxi ! Transform delta_i in complex
               do m = 1, 2*l + 1
                  mls = l*l + m + ((lmaxi + 1)**2)*(s - 1) ! composed index
                  ! Sqrt(Delta) matrix of atom i
                  cdelta_i(mls, mls) = cmplx(this%symbolic_atom(this%lattice%iz(atom_i))%potential%dele(l, s), 0.0_rp)
               end do
            end do
         end do
         do s = 1, 2
            do l = 0, lmaxj ! Transform delta_j in complex
               do m = 1, 2*l + 1
                  mls = l*l + m + ((lmaxj + 1)**2)*(s - 1) ! composed index
                  ! Sqrt(Delta) matrix of atom j
                  cdelta_j(mls, mls) = cmplx(this%symbolic_atom(this%lattice%iz(atom_j))%potential%dele(l, s), 0.0_rp)
               end do
            end do
         end do
         aux_gij(:, :, nv) = matmul(cdelta_i, matmul(green_ij(:, :, nv), cdelta_j))
      end do

      deallocate (cdelta_i)
      deallocate (cdelta_j)

   end subroutine auxiliary_gij

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Given an auxiliary GF´s (in the LMTO language), the function
   !> returns the correspondent auxiliary GF´s transformed from
   !> one representation (alpha - in) to another (beta - out).
   !> The different representations are here denoted by their
   !> set of screening constants.
   !> See Turek´s book, Eq. (3.57)
   !> Implemented by Ivan Miranda on 25.10.2023
   !---------------------------------------------------------------------------
   subroutine transform_auxiliary_gij(this, pmat_in_atom_i, pmat_out_atom_i, pmat_in_atom_j, pmat_out_atom_j, &
                                      aux_gij_in, aux_gij_out, screening_in, screening_out, atom_i, atom_j)
      !
      class(green) :: this
      !
      integer, intent(in) :: atom_i, atom_j ! Input atoms i and j
      complex(rp), intent(in) :: aux_gij_in(:, :, :) ! Auxiliary GF´s (in)
      complex(rp), intent(inout) :: aux_gij_out(:, :, :) ! Auxiliary GF´s (out), in another representation
      complex(rp), intent(in) :: pmat_in_atom_i(:, :, :) ! Potential P matrix (in) - in the same representation as GF¨s (in) for the atom i
      complex(rp), intent(in) :: pmat_out_atom_i(:, :, :) ! Potential P matrix (out) - in the same representation as GF´s (out) for the atom i
      complex(rp), intent(in) :: pmat_in_atom_j(:, :, :) ! Potential P matrix (in) - in the same representation as GF´s (in) for the atom j
      complex(rp), intent(in) :: pmat_out_atom_j(:, :, :) ! Potential P matrix (out) - in the same representation as GF´s (out) for the atom j
      real(rp), intent(in) :: screening_in(0:, :) ! Screening constants relative to the GF´s (in)
      real(rp), intent(in) :: screening_out(0:, :) ! Screening constants relative to the GF´s (out)
      !
      ! Local Variables
      integer :: l, s, m, mls, ie
      complex(rp) :: temp1, temp2
      complex(rp), allocatable, dimension(:, :, :) :: pmat_resc1, pmat_resc2, pmat_resc3 ! Re-scale P matrices

      allocate (pmat_resc1(size(pmat_in_atom_i, 1), size(pmat_in_atom_i, 2), size(pmat_in_atom_i, 3)))
      allocate (pmat_resc2(size(pmat_in_atom_i, 1), size(pmat_in_atom_i, 2), size(pmat_in_atom_i, 3)))
      allocate (pmat_resc3(size(pmat_in_atom_i, 1), size(pmat_in_atom_i, 2), size(pmat_in_atom_i, 3))) ! Only necessary if i = j

      ! First, calculate the pmat_resc1 = P^{alpha}/P^{beta} matrix (for atom i), the
      ! pmat_resc2 =  P^{alpha}/P^{beta} matrix (for atom j), and the
      ! pmat_resc3 = (beta - alpha)*P^{alpha}/P^{beta} matrix (added when i = j)

      pmat_resc1 = czero
      pmat_resc2 = czero
      pmat_resc3 = czero

      do s = 1, size(screening_in, 2) ! spin index
         do l = 0, size(screening_in, 1) - 1 ! lmax
            do m = 1, 2*l + 1
               mls = l*l + m + ((size(screening_in, 1))**2)*(s - 1) ! Composed diagonal index
               ! The screening constants can be of any atoms (i or j), because they will act only
               ! when i = j.
               temp1 = cmplx(screening_in(l, s), 0.0_rp) ! Screening constants (in)
               temp2 = cmplx(screening_out(l, s), 0.0_rp) ! Screening constants (out)
               do ie = 1, size(pmat_in_atom_i, 3) ! Energy channel
                  pmat_resc1(mls, mls, ie) = pmat_in_atom_i(mls, mls, ie)/pmat_out_atom_i(mls, mls, ie)
                  pmat_resc2(mls, mls, ie) = pmat_in_atom_j(mls, mls, ie)/pmat_out_atom_j(mls, mls, ie)
                  ! For the pmat_resc3 matrix, does not matter to consider the P potential function
                  ! for atom i or j, because it´s only acting when i = j
                  if (atom_i .eq. atom_j) then
                     pmat_resc3(mls, mls, ie) = (temp2 - temp1)*(pmat_in_atom_i(mls, mls, ie)/pmat_out_atom_i(mls, mls, ie))
                  end if
               end do
            end do
         end do
      end do

      ! Now do the rescaling of the auxiliary GF´s (in)

      do ie = 1, size(aux_gij_in, 3) ! Energy channel
         aux_gij_out(:, :, ie) = matmul(pmat_resc1(:, :, ie), matmul(aux_gij_in(:, :, ie), pmat_resc2(:, :, ie))) + &
                                 pmat_resc3(:, :, ie) ! If i != j, then the sum 0.
      end do

      deallocate (pmat_resc1)
      deallocate (pmat_resc2)
      deallocate (pmat_resc3)

   end subroutine transform_auxiliary_gij

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments form the Chebyshev recursions to calculate the onsite GF
   !---------------------------------------------------------------------------
   subroutine chebyshev_green_ij(this, istart)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      ! Local variables
      real(rp), dimension(:), allocatable :: kernel
      real(rp), dimension(:, :), allocatable :: polycheb
      real(rp), dimension(:), allocatable :: w, wscale
      real(rp) :: wstep, eps, wmin, wmax, a, b, emin_win, emax_win
      complex(rp) :: exp_factor
      integer :: ie, i, j, k, l, m, n

      this%g0 = 0.0d0

      allocate (kernel(this%control%lld*2 + 2), polycheb(this%en%channels_ldos + 10, 0:this%control%lld*2 + 2), w(this%en%channels_ldos + 10), &
                wscale(this%en%channels_ldos + 10))
      ! Defining rescaling coeficients
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp

      wscale(:) = (this%en%ene(:) - b)/a

      ! Calculating the Jackson Kernel
      call jackson_kernel((this%control%lld)*2 + 2, kernel)

      ! Calculating the Lorentz Kernel
!    call lorentz_kernel(this%control%lld, kernel, 4.0d0)

      do n = 1, 4 ! Loop on the number of on-site GFs to calculate the inter-site GFs
         ! Multiply the moments with the kernel
         do l = 1, nb
            do m = 1, nb
               this%recursion%mu_ng(l, m, :, n + istart - 1) = this%recursion%mu_n(l, m, :, n + istart - 1)*kernel(:)
            end do
         end do
         this%recursion%mu_ng(:, :, 2:size(kernel), n + istart - 1) = this%recursion%mu_ng(:, :, 2:size(kernel), n + istart - 1)*2.0_rp

         ! Calculate the Chebyshev polynomials
         call t_polynomial(size(w), size(kernel), wscale(:), polycheb)

         ! Calculate the density of states
         !$omp parallel do default(shared) private(ie, i, exp_factor, l,m)
         do ie = 1, this%en%channels_ldos + 10
            do i = 1, size(kernel)
               exp_factor = -i_unit*exp(-i_unit*(i - 1)*acos(wscale(ie)))
               do l = 1, nb
                  do m = 1, nb
                     this%g0(l, m, ie, n) = this%g0(l, m, ie, n) + this%recursion%mu_ng(l, m, i, n + istart - 1)*exp_factor
                  end do
               end do
            end do
            do l = 1, nb
               do m = 1, nb
                  this%g0(l, m, ie, n) = this%g0(l, m, ie, n)/((sqrt((a**2) - ((this%en%ene(ie) - b)**2))))
               end do
            end do
         end do
         !$omp end parallel do
      end do  ! End loop on n

      deallocate (kernel, polycheb, w, wscale)
   end subroutine chebyshev_green_ij

   !---------------------------------------------------------------------------
   !> GPU drop-in for chebyshev_green_ij: the 4 intersite-pair combos are
   !> reconstructed on the device over all nv energies in one batched moment
   !> contraction (rsrec_chebyshev_dos), reusing the on-site Chebyshev engine.
   !---------------------------------------------------------------------------
   subroutine chebyshev_green_ij_gpu(this, istart)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), allocatable :: mu_local(:, :, :, :), g0_local(:, :, :, :)
      real(rp) :: a, b, emin_win, emax_win
      integer :: nv, n_mom
      type(rsgpu), save :: gpu_backend
      logical, save :: first_call = .true.

      this%g0 = 0.0d0
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp
      nv = this%en%channels_ldos + 10
      n_mom = size(this%recursion%mu_n, 3)

      allocate (mu_local(nb, nb, n_mom, 4), g0_local(nb, nb, nv, 4))
      mu_local(:, :, :, 1:4) = this%recursion%mu_n(:, :, :, istart:istart + 3)

      if (first_call) then
         call gpu_backend%init(1, nb, 1, 1, 1)
         call gpu_backend%set_precision(0)  ! 0=fp32 (fast), 1=fp64 (validation)
         first_call = .false.
      end if

      call gpu_backend%chebyshev_dos(mu_local, this%en%ene(1:nv), a, b, g0_local)
      this%g0(:, :, 1:nv, 1:4) = g0_local(:, :, :, 1:4)
      deallocate (mu_local, g0_local)
   end subroutine chebyshev_green_ij_gpu

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments form the Chebyshev recursions to calculate the onsite GF
   !> using a small complex eta into the energy channel
   !---------------------------------------------------------------------------
   subroutine chebyshev_green_ij_eta(this, istart, eta, fermi_point, g_ef)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), dimension(nb, nb, 4), intent(inout) :: g_ef
      complex(rp), intent(in) :: eta
      integer, intent(in) :: fermi_point
      ! Local variables
      real(rp), dimension(:), allocatable :: kernel
      real(rp), dimension(:, :), allocatable :: polycheb
      real(rp), dimension(:), allocatable :: w, wscale
      real(rp) :: wstep, eps, wmin, wmax, a, b, emin_win, emax_win
      complex(rp) :: exp_factor
      integer :: ie, i, j, k, l, m, n

      g_ef = 0.0d0

      allocate (kernel(this%control%lld*2 + 2), polycheb(this%en%channels_ldos + 10, 0:this%control%lld*2 + 2), w(this%en%channels_ldos + 10), &
                wscale(this%en%channels_ldos + 10))
      ! Defining rescaling coeficients
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp

      wscale(:) = (this%en%ene(:) - b)/a

      ! Calculating the Jackson Kernel
      call jackson_kernel((this%control%lld)*2 + 2, kernel)

      ! Calculating the Lorentz Kernel
!    call lorentz_kernel(this%control%lld, kernel, 4.0d0)

      do n = 1, 4 ! Loop on the number of on-site GFs to calculate the inter-site GFs
         ! Multiply the moments with the kernel
         do l = 1, nb
            do m = 1, nb
               this%recursion%mu_ng(l, m, :, n + istart - 1) = this%recursion%mu_n(l, m, :, n + istart - 1)*kernel(:)
            end do
         end do
         this%recursion%mu_ng(:, :, 2:size(kernel), n + istart - 1) = this%recursion%mu_ng(:, :, 2:size(kernel), n + istart - 1)*2.0_rp

         ! Calculate the Chebyshev polynomials
         call t_polynomial(size(w), size(kernel), wscale(:), polycheb)

         ! Calculate the density of states
         !$omp parallel do default(shared) private(ie, i, exp_factor, l,m)
         do ie = fermi_point, fermi_point
            do i = 1, size(kernel)
               exp_factor = -i_unit*exp(-i_unit*(i - 1)*acos(((this%en%ene(ie) + eta) - b)/a))
               do l = 1, nb
                  do m = 1, nb
                     g_ef(l, m, n) = g_ef(l, m, n) + this%recursion%mu_ng(l, m, i, n + istart - 1)*exp_factor
                  end do
               end do
            end do
            do l = 1, nb
               do m = 1, nb
                  g_ef(l, m, n) = g_ef(l, m, n)/((sqrt((a**2) - (((this%en%ene(ie) + eta) - b)**2))))
               end do
            end do
         end do
         !$omp end parallel do
      end do  ! End loop on n

      deallocate (kernel, polycheb, w, wscale)
   end subroutine chebyshev_green_ij_eta

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments form the Chebyshev recursions to calculate the onsite GF
   !---------------------------------------------------------------------------
   subroutine chebyshev_green(this)
      use mpi_mod
      class(green), intent(inout) :: this
      ! Local variables
      complex(rp), allocatable :: mu_local(:, :, :, :)
      complex(rp), allocatable :: g0_local(:, :, :, :)
      real(rp) :: a, b, emin_win, emax_win
      integer :: n, n_glob, n_local, nv

      this%g0 = 0.0d0

      ! Defining rescaling coeficients
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp

      ! Number of DOS points
      nv = this%en%channels_ldos + 10

      n_local = end_atom - start_atom + 1
      if (n_local <= 0) return

      allocate (mu_local(nb, nb, size(this%recursion%mu_n, 3), n_local))
      allocate (g0_local(nb, nb, nv, n_local))

      do n_glob = start_atom, end_atom
         n = g2l_map(n_glob)
         mu_local(:, :, :, n) = this%recursion%mu_n(:, :, :, n)
      end do

      call cheb_green_fast(mu_local, nb, size(mu_local, 3), n_local, this%en%ene(1:nv), nv, a, b, g0_local)

      do n_glob = start_atom, end_atom
         n = g2l_map(n_glob)
         this%g0(:, :, 1:nv, n) = g0_local(:, :, :, n)
      end do

      deallocate (g0_local, mu_local)
!!! ! MPI moved to moment section
!!! #ifdef USE_MPI
!!!     call g_timer%start(´MPI DOS communication´)
!!!     call MPI_Allgather(this%g0(:,:,:,start_atom:end_atom), (end_atom-start_atom+1)*nv*nb*nb, MPI_DOUBLE_COMPLEX, &
!!!                        this%g0, (end_atom-start_atom+1)*nv*nb*nb, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
!!!     call g_timer%stop(´MPI DOS communication´)
!!! #endif
   end subroutine chebyshev_green

   !---------------------------------------------------------------------------
   !> Chebyshev Green function with GPU CUDA acceleration (cuBLAS)
   !> 20-100x faster than scalar loop for large nv
   !> Falls back to CPU BLAS if GPU unavailable
   !---------------------------------------------------------------------------
   subroutine chebyshev_green_gpu(this)
      use mpi_mod
      class(green), intent(inout) :: this

      ! Local variables
      real(rp), dimension(this%control%lld*2+2) :: kernel
      real(rp), dimension(this%en%channels_ldos + 10) :: wscale
      real(rp) :: a, b, emin_win, emax_win
      integer :: ie, i, j, k, l, m, n, nv, n_glob

      ! Static GPU backend (one per program lifetime)
      type(rsgpu), save :: gpu_backend
      logical, save :: first_call = .true.

      this%g0 = 0.0d0

      ! === Setup (identical to legacy path) ===
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp

      wscale(:) = (this%en%ene(:) - b)/a

      ! Number of DOS points
      nv = this%en%channels_ldos + 10

      ! Calculating the Jackson Kernel
      call jackson_kernel((this%control%lld)*2 + 2, kernel)

      ! === Kernel application loop (identical to legacy path) ===
      do n_glob = start_atom, end_atom
         n = g2l_map(n_glob)

         ! Multiply the moments with the kernel
         do l = 1, nb
            do m = 1, nb
               this%recursion%mu_ng(l, m, :, n) = this%recursion%mu_n(l, m, :, n)*kernel(:)
            end do
         end do

         this%recursion%mu_ng(:, :, 2:size(kernel), n) = &
            this%recursion%mu_ng(:, :, 2:size(kernel), n)*2.0_rp
      end do

      ! === Initialize GPU backend (first call only) ===
      if (first_call) then
         call gpu_backend%init(1, nb, 1, 1, 1)
         call gpu_backend%set_precision(0)  ! 0=fp32 (fast), 1=fp64 (validation)
         first_call = .false.
      end if

      ! === DOS reconstruction via GPU GEMM (one call replaces entire ie/i/l/m loop nest) ===
      ! Transfer matrix F(i,ie) absorbs Jackson kernel, ×2 coeff, phase, and prefactor
      ! Then G₀ = moments · F via cublasZgemmStridedBatched (GPU) or cgemm (CPU fallback)
      ! Automatic fallback to CPU if GPU unavailable
      call gpu_backend%chebyshev_dos(this%recursion%mu_n(:, :, :, start_atom:end_atom), &
                                      this%en%ene(1:nv), a, b,                         &
                                      this%g0(:, :, 1:nv, start_atom:end_atom))

   end subroutine chebyshev_green_gpu

   !---------------------------------------------------------------------------
   !> Dispatcher for Chebyshev DOS reconstruction
   !> Selects implementation based on control flags (gpu_plugin)
   !> - gpu_plugin=.true.  → chebyshev_green_gpu (GPU with CPU fallback)
   !> - otherwise          → chebyshev_green (legacy Fortran)
   !---------------------------------------------------------------------------
   subroutine chebyshev_dos_dispatch(this)
      class(green), intent(inout) :: this

      ! Try GPU path first (if compiled and enabled)
      if (this%control%gpu_plugin) then
         if (this%control%nsp > 4) then
            call g_logger%warning('gpu_plugin requested for Chebyshev DOS, '// &
               'but nsp > 4. Using legacy path.', __FILE__, __LINE__)
            call this%chebyshev_green()
            return
         end if
         call this%chebyshev_green_gpu()
         return
      end if

      ! Fall back to legacy path
      call this%chebyshev_green()

   end subroutine chebyshev_dos_dispatch

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments from the Chebyshev recursions to calculate the onsite GF
   !> using a small complex eta into the energy channel
   !---------------------------------------------------------------------------
   subroutine chebyshev_green_eta(this, eta, fermi_point, g_ef)
      use mpi_mod
      class(green), intent(inout) :: this
      !complex(rp), dimension(nb, nb,this%lattice%nrec), intent(inout) :: g_ef
      complex(rp), dimension(nb, nb, atoms_per_process), intent(inout) :: g_ef
      complex(rp), intent(in) :: eta
      integer, intent(in) :: fermi_point
      ! Local variables
      real(rp), dimension(this%control%lld*2 + 2) :: kernel
      real(rp), dimension(this%en%channels_ldos + 10, 0:this%control%lld*2 + 2) :: polycheb
      real(rp), dimension(this%en%channels_ldos + 10) :: w, wscale
      real(rp) :: wstep, eps, wmin, wmax, a, b, emin_win, emax_win
      complex(rp) :: exp_factor
      integer :: ie, i, j, k, l, m, n, nv
      integer :: n_glob

      g_ef = 0.0d0

      ! Defining rescaling coeficients
      call this%recursion%resolve_chebyshev_window(emin_win, emax_win)
      a = (emax_win - emin_win)/(2 - 0.3_rp)
      b = (emax_win + emin_win)/2.0_rp

      wscale(:) = (this%en%ene(:) - b)/a

      ! Number of DOS points
      nv = this%en%channels_ldos + 10

      ! Calculating the Jackson Kernel
      call jackson_kernel((this%control%lld)*2 + 2, kernel)

      ! Calculating the Lorentz Kernel
      ! call lorentz_kernel(this%control%lld, kernel, 4.0d0)

      do n_glob = start_atom, end_atom ! Loop on self-consistent atoms
         n = g2l_map(n_glob)
         ! Multiply the moments with the kernel
         do i = 1, nb
            this%recursion%mu_ng(i, i, :, n) = this%recursion%mu_n(i, i, :, n)*kernel(:)
         end do

         this%recursion%mu_ng(:, :, 2:size(kernel), n) = this%recursion%mu_ng(:, :, 2:size(kernel), n)*2.0_rp

         !do i=1, size(kernel)
         !  write(400+n, *) i, sum(this%recursion%mu_n(n, i, 1:18, 1:18))
         !end do

         ! Calculate the Chebyshev polynomials
         call t_polynomial(size(w), size(kernel), wscale(:), polycheb)

         ! Calculate the density of states
         !$omp parallel do default(shared) private(ie, i, exp_factor, l,m)
         do ie = fermi_point, fermi_point
            do i = 1, size(kernel)
               exp_factor = -i_unit*exp(-i_unit*(i - 1)*acos(((this%en%ene(ie) + eta) - b)/a))
               do l = 1, nb
                  do m = 1, nb
                     g_ef(l, m, n) = g_ef(l, m, n) + this%recursion%mu_ng(l, m, i, n)*exp_factor
                  end do
               end do
            end do
            do l = 1, nb
               do m = 1, nb
                  g_ef(l, m, n) = g_ef(l, m, n)/((sqrt((a**2) - (((this%en%ene(ie) + eta) - b)**2))))
               end do
            end do
         end do
         !$omp end parallel do
      end do  ! End loop on self-consistent atoms
   end subroutine chebyshev_green_eta

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Calculates the Greens function
   !---------------------------------------------------------------------------
   subroutine bgreen(this, g_out, i_site, ie_start, ie_len, a_inf, b_inf, eta)
      implicit none
      class(green), intent(inout) :: this
      integer, intent(in) :: i_site
      integer, intent(in) :: ie_start
      integer, intent(in) :: ie_len
      complex(rp), dimension(nb, nb, this%en%channels_ldos + 10), intent(inout) :: g_out
      real(rp), dimension(nb, nb), intent(in) :: a_inf
      real(rp), dimension(nb, nb), intent(in) :: b_inf
      complex(rp), intent(in) :: eta
      !
      integer :: nv, ldim, ll, ll_t, llinf
      real(rp) :: factor_z
      real(rp), dimension(this%en%channels_ldos + 10) :: e
      integer :: i, j, l, ei, info, ln, lwork
      integer :: ii, jj
      real(rp) :: recMaxA, recMaxB, valrec
      real(rp) :: maxQ, maxB2r
      real(rp) :: qMax, wMax, b2Max
      real(rp) :: qval, wval, b2val
      logical :: qHasNaN
      logical :: recNaN
      logical :: found_nan
      complex(rp) :: vv
      integer, dimension(nb) :: ipiv
      complex(rp) :: etop, ebot, ea, eb
      real(rp), dimension(this%lattice%nrec) :: a_inf0, b_inf0
      complex(rp), dimension(nb, nb) :: Q, Qp, Q2p, Z, one, W, B2z, Qt, P
      complex(rp) :: zoff, cone, czero, ze, zterm, det, im
      complex(rp) :: coff
      complex(rp), dimension(nb*nb) :: work
      real(rp), dimension(this%en%channels_ldos + 10) :: ene
      complex(rp), dimension(nb, nb) :: Dfac_mat, Cshi_mat
      real(rp) :: a_diag, b_diag
      !
      integer, dimension(nb) :: m_tab
      !
      integer :: n_glob
      logical, save :: dump_done = .false.
      character(len=200) :: dump_fname
      integer :: dump_unit, dump_iostat
      !
      g_out = (0.0d0, 0.0d0)
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      nv = this%en%channels_ldos
      ldim = nb
      e = this%en%ene
      factor_z = 1.0d0

      ! Lightweight runtime shapes/logging (once per green instance call)
      ! call g_logger%info('DEBUG:bgreen shapes nb='//int2str(nb)//' ldim='//int2str(ldim), __FILE__, __LINE__)
      ! call g_logger%info('DEBUG:bgreen a_b dims=('//int2str(size(this%recursion%a_b,1))//','//int2str(size(this%recursion%a_b,2))//','//int2str(size(this%recursion%a_b,3))//','//int2str(size(this%recursion%a_b,4))//')', __FILE__, __LINE__)

      ! Unchanged code
      ll_t = ll
      llinf = ll!+ll/2!0!3*ll!/2!+30 !300
      cone = (1.0d0, 0.0d0)
      czero = (0.0d0, 0.0d0)
      im = (0.0d0, 1.0d0)
      zoff = (0.0d0, 0.01d0)*0.0d0
      coff = 0.000d0*(0.0d0, 1.0d0)
      one = (0.0d0, 0.0d0)

      ! Initializing Shift and Scaling for rigid band shift
      Dfac_mat = (1.0d0, 0.0d0)
      Cshi_mat = (0.0d0, 0.0d0)

      do i = 1, ldim
         m_tab(i) = i
      end do

      do i = 1, ldim
         one(i, i) = (1.0d0, 0.0d0)
      end do
      if (ldim >= 10) then
         a_diag = (a_inf(1, 1) + a_inf(10, 10))*0.5d0
         b_diag = (b_inf(1, 1) + b_inf(10, 10))*0.5d0
      else
         a_diag = a_inf(1, 1)
         b_diag = b_inf(1, 1)
      end if
      ! Standard semi-circle construction of Greens functions
      !write(12334,´(18f8.4)´) a_inf
      !write(12335,´(18f8.4)´) b_inf
      !write(12336,´(18f8.4)´) a_inf
      !write(12337,´(18f8.4)´) this%recursion%a_b(:,:,2,i_site)
      !write(12337,´(18f8.4)´) this%recursion%b2_b(:,:,2,i_site)
    !!!$omp parallel do default(shared) private(ei, Z, Ze, Q, B2z, i, ea, eb, det, zoff, l, ln, j, P , ipiv, info, work, lwork, W)
      !$omp parallel do default(shared) &
      !$omp          private(ei, Z, Ze, Q, B2z, i, etop, ebot, ea, eb, det, zoff, l, ln, j, P, ipiv, info, work, lwork, W, &
      !$omp                  found_nan, qHasNaN, qMax, wMax, b2Max, qval, wval, b2val, ii, jj, &
      !$omp                  recMaxA, recMaxB, recNaN, valrec, vv, maxQ, maxB2r)
      do ei = ie_start, ie_start + ie_len - 1
         Z = e(ei)*one
         ze = e(ei)!-zoff
         Q = (0.0d0, 0.0d0)
         B2z = 0.0d0
         found_nan = .false.
         if (this%control%sym_term) then
            do i = 1, ldim !  Orbital-independent
               etop = a_diag + 2.0d0*b_diag
               ebot = a_diag - 2.0d0*b_diag
               ea = e(ei) - etop
               eb = e(ei) - ebot
               det = ea*eb
               zoff = sqrt(det)
               Q(i, i) = (ze + eta - a_diag - zoff)*0.5d0
            end do
         else
            do i = 1, ldim  ! Orbital-dependent
               if (i == 1 .or. (ldim >= 10 .and. i == 10)) then  !(now the s-bands are broadened earlier)
                  etop = a_inf(i, i) - Cshi_mat(i, i) + 2*b_inf(i, i)*1.025d0/Dfac_mat(i, i)
                  ebot = a_inf(i, i) - Cshi_mat(i, i) - 2*b_inf(i, i)*1.025d0/Dfac_mat(i, i)
               else
                  etop = a_inf(i, i) - Cshi_mat(i, i) + 2*b_inf(i, i)!*1.01d0/Dfac_mat(i,i)
                  ebot = a_inf(i, i) - Cshi_mat(i, i) - 2*b_inf(i, i)!*1.01d0/Dfac_mat(i,i)
               end if
               ea = e(ei)/Dfac_mat(i, i) - etop - 1.0d0*Cshi_mat(i, i)
               eb = e(ei)/Dfac_mat(i, i) - ebot - 1.0d0*Cshi_mat(i, i)
               det = ea*eb
               zoff = sqrt(det)!*0.5d0
               zoff = zoff*(1.0)**ei
           !!! if(nv<=10) zoff=zoff*1.0d-3   ! Hack for Efermi evaluation
               ! Below for orbital dependent terminator
               Q(i, i) = ((e(ei) + eta)/Dfac_mat(i, i) - 1.0d0*Cshi_mat(i, i) - a_inf(i, i) - real(zoff) - aimag(zoff)*factor_z*im)*0.5d0
            end do
         end if
         do l = llinf - 1, 1, -1
            ln = min(l, ll)
            do j = 1, ldim
               do i = 1, ldim
                  if (real(Z(i, j)) .ne. 0) then
                     P(i, j) = (Z(i, j) + eta)
                  else
                     P(i, j) = Z(i, j)
                  end if
                  if ((abs(real(Q(i, j))) .lt. 1.0d-12) .and. (abs(aimag(Q(i, j))) .lt. 1.0d-12)) then
                     Q(i, j) = (0.0d0, 0.0d0)
                  end if
                  Q(i, j) = P(i, j)/Dfac_mat(i, j) - Cshi_mat(i, j) - cone*this%recursion%a_b(i, j, ln, i_site) - Q(i, j)
                  B2z(i, j) = cone*this%recursion%b2_b(i, j, ln, i_site)
               end do
            end do
            ! Check Q for NaNs before LU factorization
            qMax = 0.0_rp
            qHasNaN = .false.
            do jj = 1, ldim
               do ii = 1, ldim
                  qval = abs(Q(ii, jj))
                  if (qval == qval) then
                     if (qval > qMax) qMax = qval
                  else
                     qHasNaN = .true.
                  end if
               end do
            end do
            if (qHasNaN) then
               call g_logger%error('Pre-LU Q contains NaN (site='//int2str(i_site)//' ie='//int2str(ei)//' Q_max='//real2str(qMax)//')', __FILE__, __LINE__)
               !$omp critical
               if (.not. dump_done) then
                  dump_done = .true.
                  dump_fname = 'debug_bgreen_preLU_site'//int2str(i_site)//'_ie'//int2str(ei)//'.txt'
                  dump_unit = 901
                  open(dump_unit, file=trim(dump_fname), status='replace', action='write', iostat=dump_iostat)
                  if (dump_iostat == 0) then
                     write(dump_unit, '(A)') 'PRE-LU DEBUG DUMP: site='//int2str(i_site)//' ie='//int2str(ei)
                        write(dump_unit, '(A)') 'SHAPES: nb='//int2str(nb)//' ldim='//int2str(ldim)//' ll='//int2str(ll)
                        write(dump_unit, '(A)') 'Small P (real imag) first 4x4:'
                        do jj = 1, min(4,ldim)
                           do ii = 1, min(4,ldim)
                              write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(P(ii,jj)), aimag(P(ii,jj))
                           end do
                        end do
                        write(dump_unit, '(A)') 'Small Dfac_mat (real imag) first 4x4:'
                        do jj = 1, min(4,ldim)
                           do ii = 1, min(4,ldim)
                              write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(Dfac_mat(ii,jj)), aimag(Dfac_mat(ii,jj))
                           end do
                        end do
                        write(dump_unit, '(A)') 'Small Cshi_mat (real imag) first 4x4:'
                        do jj = 1, min(4,ldim)
                           do ii = 1, min(4,ldim)
                              write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(Cshi_mat(ii,jj)), aimag(Cshi_mat(ii,jj))
                           end do
                        end do
                        write(dump_unit, '(A)') 'Sample recursion contributions (ln=1..min(4,ll)):'
                        do ln = 1, min(4,ll)
                           write(dump_unit, '(A,I0)') 'ln=', ln
                           do jj = 1, min(4,ldim)
                              do ii = 1, min(4,ldim)
                                 write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16,1X,ES24.16,1X,ES24.16)') ii, jj, real(this%recursion%a_b(ii,jj,ln,i_site)), aimag(this%recursion%a_b(ii,jj,ln,i_site)), real(this%recursion%b2_b(ii,jj,ln,i_site)), aimag(this%recursion%b2_b(ii,jj,ln,i_site))
                              end do
                           end do
                        end do
                        write(dump_unit, '(A)') '--- continuing with full Q/B2z dump below ---'
                     write(dump_unit, '(A)') 'Q matrix (real imag):'
                     do jj = 1, ldim
                        do ii = 1, ldim
                           write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(Q(ii,jj)), aimag(Q(ii,jj))
                        end do
                     end do
                     write(dump_unit, '(A)') 'B2z matrix (real imag):'
                     do jj = 1, ldim
                        do ii = 1, ldim
                           write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(B2z(ii,jj)), aimag(B2z(ii,jj))
                        end do
                     end do
                     close(dump_unit)
                  end if
               end if
               !$omp end critical
               found_nan = .true.
               exit
            end if

            ! LU factorization
            call zgetrf(ldim, ldim, Q, ldim, ipiv, info)
            if (info /= 0) then
               ! Collect diagnostics and attempt tiny regularization then retry once
               maxQ = 0.0_rp
               maxB2r = 0.0_rp
               do jj = 1, ldim
                  do ii = 1, ldim
                     if (abs(Q(ii, jj)) > maxQ) maxQ = abs(Q(ii, jj))
                     if (abs(B2z(ii, jj)) > maxB2r) maxB2r = abs(B2z(ii, jj))
                  end do
               end do
               recMaxA = 0.0_rp
               recMaxB = 0.0_rp
               recNaN = .false.
               do ln = 1, ll
                  do ii = 1, ldim
                     do jj = 1, ldim
                        valrec = abs(this%recursion%a_b(ii, jj, ln, i_site))
                        if (valrec > recMaxA) recMaxA = valrec
                        if (IsNaN(real(this%recursion%a_b(ii, jj, ln, i_site))) .or. IsNaN(aimag(this%recursion%a_b(ii, jj, ln, i_site)))) recNaN = .true.
                        valrec = abs(this%recursion%b2_b(ii, jj, ln, i_site))
                        if (valrec > recMaxB) recMaxB = valrec
                        if (IsNaN(real(this%recursion%b2_b(ii, jj, ln, i_site))) .or. IsNaN(aimag(this%recursion%b2_b(ii, jj, ln, i_site)))) recNaN = .true.
                     end do
                  end do
               end do
               call g_logger%warning('LU factorization failed (zgetrf info='//fmt('I0', info)//') site='//int2str(i_site)//' ie='//int2str(ei)//' maxQ='//real2str(maxQ)//' maxB2='//real2str(maxB2r)//' recMaxA='//real2str(recMaxA)//' recMaxB='//real2str(recMaxB), __FILE__, __LINE__)
               ! Try tiny regularization
               Q = Q + (1.0d-12, 0.0d0)*one
               call zgetrf(ldim, ldim, Q, ldim, ipiv, info)
               if (info /= 0) then
                  call g_logger%error('zgetrf retry failed (info='//fmt('I0', info)//') skipping energy slice site='//int2str(i_site)//' ie='//int2str(ei), __FILE__, __LINE__)
                  found_nan = .true.
                  exit
               end if
            end if
            ! Inverse of Q from LU factorization
            lwork = ldim*ldim
            call zgetri(ldim, Q, ldim, ipiv, work, lwork, info)
            ! Check Q for NaNs after inversion
            qMax = 0.0_rp
            qHasNaN = .false.
            do jj = 1, ldim
               do ii = 1, ldim
                  qval = abs(Q(ii, jj))
                  if (qval == qval) then
                     if (qval > qMax) qMax = qval
                  else
                     qHasNaN = .true.
                  end if
               end do
            end do
            if (qHasNaN) then
               call g_logger%error('Post-inverse Q contains NaN (site='//int2str(i_site)//' ie='//int2str(ei)//' Q_max='//real2str(qMax)//' zgetri_info='//fmt('I0', info)//')', __FILE__, __LINE__)
               !$omp critical
               if (.not. dump_done) then
                  dump_done = .true.
                  dump_fname = 'debug_bgreen_postInv_site'//int2str(i_site)//'_ie'//int2str(ei)//'.txt'
                  dump_unit = 902
                  open(dump_unit, file=trim(dump_fname), status='replace', action='write', iostat=dump_iostat)
                  if (dump_iostat == 0) then
                     write(dump_unit, '(A)') 'POST-INV DEBUG DUMP: site='//int2str(i_site)//' ie='//int2str(ei)
                     write(dump_unit, '(A)') 'Q matrix (real imag):'
                     do jj = 1, ldim
                        do ii = 1, ldim
                           write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(Q(ii,jj)), aimag(Q(ii,jj))
                        end do
                     end do
                     write(dump_unit, '(A)') 'B2z matrix (real imag):'
                     do jj = 1, ldim
                        do ii = 1, ldim
                           write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(B2z(ii,jj)), aimag(B2z(ii,jj))
                        end do
                     end do
                     close(dump_unit)
                  end if
               end if
               !$omp end critical
               found_nan = .true.
               exit
            end if
            if (info /= 0) then
               call g_logger%warning('zgetri (inverse) failed (info='//fmt('I0', info)//') site='//int2str(i_site)//' ie='//int2str(ei), __FILE__, __LINE__)
               ! Attempt small regularization and retry inverse
               Q = Q + (1.0d-12, 0.0d0)*one
               call zgetrf(ldim, ldim, Q, ldim, ipiv, info)
               if (info == 0) then
                  call zgetri(ldim, Q, ldim, ipiv, work, lwork, info)
               end if
               if (info /= 0) then
                  call g_logger%error('zgetri retry failed (info='//fmt('I0', info)//') skipping energy slice site='//int2str(i_site)//' ie='//int2str(ei), __FILE__, __LINE__)
                  found_nan = .true.
                  exit
               end if
            end if
            call zgemm('n', 'n', ldim, ldim, ldim, cone, Q, ldim, B2z, ldim, czero, W, ldim)
            call zgemm('c', 'n', ldim, ldim, ldim, cone, B2z, ldim, W, ldim, czero, Q, ldim)
         end do
         !
         if (.not. found_nan) then
            do j = 1, ldim
               do i = 1, ldim
                  !bdos(ei,m_tab(i),m_tab(j),n)=bdos(ei,m_tab(i),m_tab(j),n) + abs(-aimag(Q(i,j))/3.14159265359d0)/Dfac_mat(i,j)
                  if (aimag(eta) .eq. 0) then
                     g_out(m_tab(i), m_tab(j), ei) = g_out(m_tab(i), m_tab(j), ei) + &
                                                     real(Q(i, j)/Dfac_mat(i, j)) + aimag(Q(i, j)/Dfac_mat(i, j))*(1.0d0)**ei*(0.0d0, 1.0d0)
                  else
                     g_out(m_tab(i), m_tab(j), ei) = g_out(m_tab(i), m_tab(j), ei) + (Q(i, j)/Dfac_mat(i, j))
                  end if
               end do
            end do
         end if
         ! Quick NaN check for this energy slice (log first occurrence)
         found_nan = .false.
         do j = 1, ldim
            do i = 1, ldim
               vv = g_out(m_tab(i), m_tab(j), ei)
               if (real(vv) /= real(vv) .or. aimag(vv) /= aimag(vv)) then
                  call g_logger%warning('NaN detected in g_out at site '//fmt('I0', i_site)//' ie='//fmt('I0', ei), __FILE__, __LINE__)
                  found_nan = .true.
                        ! Collect recursion coefficient summaries for this site
                        recMaxA = 0.0_rp
                        recMaxB = 0.0_rp
                        recNaN = .false.
                        do ln = 1, ll
                           do ii = 1, ldim
                              do jj = 1, ldim
                                 valrec = abs(this%recursion%a_b(ii, jj, ln, i_site))
                                 if (valrec > recMaxA) recMaxA = valrec
                                 if (IsNaN(real(this%recursion%a_b(ii, jj, ln, i_site))) .or. IsNaN(aimag(this%recursion%a_b(ii, jj, ln, i_site)))) recNaN = .true.
                                 valrec = abs(this%recursion%b2_b(ii, jj, ln, i_site))
                                 if (valrec > recMaxB) recMaxB = valrec
                                 if (IsNaN(real(this%recursion%b2_b(ii, jj, ln, i_site))) .or. IsNaN(aimag(this%recursion%b2_b(ii, jj, ln, i_site)))) recNaN = .true.
                              end do
                           end do
                        end do
                        ! call g_logger%info('DEBUG:bgreen recursion site='//int2str(i_site)//' ie='//int2str(ei)//' maxA='//real2str(recMaxA)//' maxB='//real2str(recMaxB)//' recNaN='//log2str(recNaN), __FILE__, __LINE__)
                        ! ! Log diagonal of a_inf/b_inf and Dfac_mat for quick inspection
                        ! do ii = 1, ldim
                        !    call g_logger%info('DEBUG:bgreen a_inf('//int2str(ii)//')='//real2str(a_inf(ii,ii))//' b_inf='//real2str(b_inf(ii,ii))//' Dfac='//real2str(real(Dfac_mat(ii,ii))), __FILE__, __LINE__)
                        ! end do
                        ! Additional diagnostics: check Q, W and B2z for NaNs/large values
                        qMax = 0.0_rp
                        qHasNaN = .false.
                        wMax = 0.0_rp
                        b2Max = 0.0_rp
                        do jj = 1, ldim
                           do ii = 1, ldim
                              qval = abs(Q(ii,jj))
                              if (qval == qval) then
                                 if (qval > qMax) qMax = qval
                              else
                                 qHasNaN = .true.
                              end if
                              wval = abs(W(ii,jj))
                              if (wval == wval) then
                                 if (wval > wMax) wMax = wval
                              end if
                              b2val = abs(B2z(ii,jj))
                              if (b2val == b2val) then
                                 if (b2val > b2Max) b2Max = b2val
                              end if
                           end do
                        end do
                        ! call g_logger%info('DEBUG:bgreen Q_max='//real2str(qMax)//' Q_hasNaN='//log2str(qHasNaN)//' W_max='//real2str(wMax)//' B2z_max='//real2str(b2Max), __FILE__, __LINE__)
                       !$omp critical
                       if (.not. dump_done) then
                          dump_done = .true.
                          dump_fname = 'debug_bgreen_site'//int2str(i_site)//'_ie'//int2str(ei)//'.txt'
                          dump_unit = 900
                          open(dump_unit, file=trim(dump_fname), status='replace', action='write', iostat=dump_iostat)
                          if (dump_iostat == 0) then
                             write(dump_unit, '(A)') 'DEBUG DUMP: site='//int2str(i_site)//' ie='//int2str(ei)
                             write(dump_unit, '(A)') 'Q matrix (real imag):'
                             do jj = 1, ldim
                                do ii = 1, ldim
                                   write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(Q(ii,jj)), aimag(Q(ii,jj))
                                end do
                             end do
                             write(dump_unit, '(A)') 'W matrix (real imag):'
                             do jj = 1, ldim
                                do ii = 1, ldim
                                   write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(W(ii,jj)), aimag(W(ii,jj))
                                end do
                             end do
                             write(dump_unit, '(A)') 'B2z matrix (real imag):'
                             do jj = 1, ldim
                                do ii = 1, ldim
                                   write(dump_unit, '(I4,1X,I4,1X,ES24.16,1X,ES24.16)') ii, jj, real(B2z(ii,jj)), aimag(B2z(ii,jj))
                                end do
                             end do
                             write(dump_unit, '(A)') 'recursion a_b (ln, i, j, real imag):'
                             do ln = 1, ll
                                do jj = 1, ldim
                                   do ii = 1, ldim
                                      write(dump_unit, '(I4,1X,I4,1X,I4,1X,ES24.16,1X,ES24.16)') ln, ii, jj, real(this%recursion%a_b(ii,jj,ln,i_site)), aimag(this%recursion%a_b(ii,jj,ln,i_site))
                                   end do
                                end do
                             end do
                             write(dump_unit, '(A)') 'recursion b2_b (ln, i, j, real imag):'
                             do ln = 1, ll
                                do jj = 1, ldim
                                   do ii = 1, ldim
                                      write(dump_unit, '(I4,1X,I4,1X,I4,1X,ES24.16,1X,ES24.16)') ln, ii, jj, real(this%recursion%b2_b(ii,jj,ln,i_site)), aimag(this%recursion%b2_b(ii,jj,ln,i_site))
                                   end do
                                end do
                             end do
                             close(dump_unit)
                          end if
                       end if
                       !$omp end critical
                        exit
               end if
            end do
            if (found_nan) exit
         end do
         !write(12333,´(19f8.4)´) e(ei), (real(g_out(i,i,ei)),i=1,18)
      end do
      !$omp end parallel do

   end subroutine bgreen

end module
