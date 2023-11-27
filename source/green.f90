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
   use precision_mod, only: rp
   use math_mod
   use logger_mod, only: g_logger
   use timer_mod, only: g_timer
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Module's main structure
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
   contains
      procedure :: sgreen
      procedure :: bgreen
      procedure :: block_green
      procedure :: block_green_eta
      procedure :: block_green_ij
      procedure :: block_green_ij_eta
      procedure :: calculate_intersite_gf
      procedure :: calculate_intersite_gf_eta
      procedure :: chebyshev_green
      procedure :: chebyshev_green_eta
      procedure :: chebyshev_green_ij
      procedure :: chebyshev_green_ij_eta
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
         call g_safe_alloc%allocate('green.g0', this%g0, (/18, 18, this%en%channels_ldos + 10, atoms_per_process/))
      else
         call g_safe_alloc%allocate('green.g0', this%g0, (/18, 18, this%en%channels_ldos + 10, 4/))
      end if
      call g_safe_alloc%allocate('green.gij', this%gij, (/18, 18, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gji', this%gji, (/18, 18, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.ginmag', this%ginmag, (/9, 9, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjnmag', this%gjnmag, (/9, 9, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gix', this%gix, (/9, 9, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.giy', this%giy, (/9, 9, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.giz', this%giz, (/9, 9, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjx', this%gjx, (/9, 9, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjy', this%gjy, (/9, 9, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjz', this%gjz, (/9, 9, this%en%channels_ldos + 10, atoms_per_process/))
      call g_safe_alloc%allocate('green.gij_eta', this%gij_eta, (/64, 18, 18, atoms_per_process/))
      call g_safe_alloc%allocate('green.gji_eta', this%gji_eta, (/64, 18, 18, atoms_per_process/))
      call g_safe_alloc%allocate('green.ginmag_eta', this%ginmag_eta, (/64, 9, 9, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjnmag_eta', this%gjnmag_eta, (/64, 9, 9, atoms_per_process/))
      call g_safe_alloc%allocate('green.gix_eta', this%gix_eta, (/64, 9, 9, atoms_per_process/))
      call g_safe_alloc%allocate('green.giy_eta', this%giy_eta, (/64, 9, 9, atoms_per_process/))
      call g_safe_alloc%allocate('green.giz_eta', this%giz_eta, (/64, 9, 9, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjx_eta', this%gjx_eta, (/64, 9, 9, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjy_eta', this%gjy_eta, (/64, 9, 9, atoms_per_process/))
      call g_safe_alloc%allocate('green.gjz_eta', this%gjz_eta, (/64, 9, 9, atoms_per_process/))
#else
      if (this%lattice%njij == 0) then
         allocate (this%g0(18, 18, this%en%channels_ldos + 10, this%lattice%nrec))
      else
         allocate (this%g0(18, 18, this%en%channels_ldos + 10, 4))
      end if
      allocate (this%gij(18, 18, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gji(18, 18, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%ginmag(9, 9, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gjnmag(9, 9, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gix(9, 9, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%giy(9, 9, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%giz(9, 9, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gjx(9, 9, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gjy(9, 9, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gjz(9, 9, this%en%channels_ldos + 10, atoms_per_process))
      allocate (this%gij_eta(64, 18, 18, atoms_per_process))
      allocate (this%gji_eta(64, 18, 18, atoms_per_process))
      allocate (this%ginmag_eta(64, 9, 9, atoms_per_process))
      allocate (this%gjnmag_eta(64, 9, 9, atoms_per_process))
      allocate (this%gix_eta(64, 9, 9, atoms_per_process))
      allocate (this%giy_eta(64, 9, 9, atoms_per_process))
      allocate (this%giz_eta(64, 9, 9, atoms_per_process))
      allocate (this%gjx_eta(64, 9, 9, atoms_per_process))
      allocate (this%gjy_eta(64, 9, 9, atoms_per_process))
      allocate (this%gjz_eta(64, 9, 9, atoms_per_process))
#endif

      this%g0(:, :, :, :) = (0.0D0, 0.0D0)
      this%gij(:, :, :, :) = (0.0D0, 0.0D0)
      this%gji(:, :, :, :) = (0.0D0, 0.0D0)
      this%ginmag(:, :, :, :) = (0.0D0, 0.0D0)
      this%gjnmag(:, :, :, :) = (0.0D0, 0.0D0)
      this%gix(:, :, :, :) = (0.0D0, 0.0D0)
      this%giy(:, :, :, :) = (0.0D0, 0.0D0)
      this%giz(:, :, :, :) = (0.0D0, 0.0D0)
      this%gjx(:, :, :, :) = (0.0D0, 0.0D0)
      this%gjy(:, :, :, :) = (0.0D0, 0.0D0)
      this%gjz(:, :, :, :) = (0.0D0, 0.0D0)
      this%gij_eta(:, :, :, :) = (0.0D0, 0.0D0)
      this%gji_eta(:, :, :, :) = (0.0D0, 0.0D0)
      this%ginmag_eta(:, :, :, :) = (0.0D0, 0.0D0)
      this%gjnmag_eta(:, :, :, :) = (0.0D0, 0.0D0)
      this%gix_eta(:, :, :, :) = (0.0D0, 0.0D0)
      this%giy_eta(:, :, :, :) = (0.0D0, 0.0D0)
      this%giz_eta(:, :, :, :) = (0.0D0, 0.0D0)
      this%gjx_eta(:, :, :, :) = (0.0D0, 0.0D0)
      this%gjy_eta(:, :, :, :) = (0.0D0, 0.0D0)
      this%gjz_eta(:, :, :, :) = (0.0D0, 0.0D0)
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
      complex(rp), dimension(18, 18, 1, 4), intent(inout) :: g_ef
      integer, intent(in) :: fermi_point
      !
      integer :: ll, n, nw, ldim, na
      real(rp), dimension(4) :: a_inf0, b_inf0
      real(rp), dimension(18, 18, 4) :: a_inf, b_inf
      !
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      ldim = 18
      na = 4
      nw = 10 * ll

      call this%recursion%get_terminf(this%recursion%a_b(:, :, :, istart), this%recursion%b2_b(:, :, :, istart), na, ll, ldim, nw, &
                                      a_inf, b_inf, a_inf0, b_inf0)

      do n = 1, NA

         call this%bgreen(g_ef(:, :, :, n), n + istart - 1, 1, this%en%channels_ldos + 10, a_inf(:, :, n), b_inf(:, :, n), eta)

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
      real(rp), dimension(18, 18, 4) :: a_inf, b_inf
      !
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      ldim = 18
      na = 4
      eta = (0.0D0, 0.0D0)
      nw = 10 * ll

      call this%recursion%get_terminf(this%recursion%a_b(:, :, :, istart), this%recursion%b2_b(:, :, :, istart), na, ll, ldim, nw, &
                                      a_inf, b_inf, a_inf0, b_inf0)

      do n = 1, NA
         ! Initializing Shift and Scaling for rigid band shift
         !Dfac_mat=(1.0d0,0.0d0)
         !Cshi_mat=(0.0d0,0.0d0)
         call this%bgreen(this%g0(:, :, :, n), n + istart - 1, 1, this%en%channels_ldos + 10, a_inf(:, :, n), b_inf(:, :, n), eta)

      end do
   end subroutine block_green_ij

   subroutine calculate_intersite_gf(this)
      use mpi_mod
      implicit none
      class(green), intent(inout) :: this
      integer :: ia, ja_temp, j, i, ia_glob

      this%gij = 0.0D0; this%gji = 0.0D0; this%ginmag = 0.0D0; this%giz = 0.0D0; this%giy = 0.0D0; this%gix = 0.0D0
      this%gjnmag = 0.0D0; this%gjx = 0.0D0; this%gjy = 0.0D0; this%gjz = 0.0D0

      if (this%control%recur == 'block') call this%recursion%zsqr()

      !do ia=1,this%lattice%njij
      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)
         ja_temp = (ia - 1) * 4 + 1
         select case (this%control%recur)
         case ('block')
            call this%block_green_ij(ja_temp)
         case ('chebyshev')
            call this%chebyshev_green_ij(ja_temp)
         end select
         if (this%lattice%ijpair(ia, 1) .eq. this%lattice%ijpair(ia, 2)) then
            this%gij(:, :, :, ia) = this%g0(:, :, :, 1)
            this%gji(:, :, :, ia) = this%g0(:, :, :, 1)
         else
            this%gij(:, :, :, ia) = this%g0(:, :, :, 1) - this%g0(:, :, :, 2) + (1.0D0 / i_unit * this%g0(:, :, :, 3) - 1.0D0 / i_unit * this%g0(:, :, :, 4))
            this%gji(:, :, :, ia) = this%g0(:, :, :, 1) - this%g0(:, :, :, 2) - (1.0D0 / i_unit * this%g0(:, :, :, 3) - 1.0D0 / i_unit * this%g0(:, :, :, 4))
            this%gij(:, :, :, ia) = this%gij(:, :, :, ia) * 0.5D0
            this%gji(:, :, :, ia) = this%gji(:, :, :, ia) * 0.5D0
         end if
         do i = 1, 9
            do j = 1, 9
               this%Ginmag(j, i, :, iA) = this%Ginmag(j, i, :, iA) + (this%gij(j, i, :, ia) + this%gij(j + 9, i + 9, :, ia)) * 0.5D0
               this%Giz(j, i, :, iA) = this%Giz(j, i, :, iA) + 0.5D0 * (this%gij(j, i, :, ia) - this%gij(j + 9, i + 9, :, ia))        !+Ginmag(k,j,i,iIA)*0.5d0
               this%Giy(j, i, :, iA) = this%Giy(j, i, :, iA) + 0.5D0 * (i_unit * this%gij(j, i + 9, :, ia) - i_unit * this%gij(j + 9, i, :, ia))  !+Ginmag(k,j,i,iIA)*0.5d0
               this%Gix(j, i, :, iA) = this%Gix(j, i, :, iA) + 0.5D0 * (this%gij(j, i + 9, :, ia) + this%gij(j + 9, i, :, ia))        !+Ginmag(k,j,i,iIA)*0.5d0
               !
               this%Gjnmag(j, i, :, iA) = this%Gjnmag(j, i, :, iA) + (this%gji(j, i, :, ia) + this%gji(j + 9, i + 9, :, ia)) * 0.5D0
               this%Gjz(j, i, :, iA) = this%Gjz(j, i, :, iA) + 0.5D0 * (this%gji(j, i, :, ia) - this%gji(j + 9, i + 9, :, ia))         !+Gjnmag(k,j,i,iIA)*0.5d0
               this%Gjy(j, i, :, iA) = this%Gjy(j, i, :, iA) + 0.5D0 * (i_unit * this%gji(j, i + 9, :, ia) - i_unit * this%gji(j + 9, i, :, ia))   !+Gjnmag(k,j,i,iIA)*0.5d0
               this%Gjx(j, i, :, iA) = this%Gjx(j, i, :, iA) + 0.5D0 * (this%gji(j, i + 9, :, ia) + this%gji(j + 9, i, :, ia))         !+Gjnmag(k,j,i,iIA)*0.5d0
            end do
         end do
      end do
   end subroutine calculate_intersite_gf

   subroutine calculate_intersite_gf_eta(this)
      use mpi_mod
      implicit none
      class(green), intent(inout) :: this
      integer :: ia, ja_temp, j, i, fermi_point
      complex(rp), dimension(18, 18, 4) :: g0_ef
      complex(rp) :: eta
      complex(rp), dimension(64, 18, 18, 4) :: y
      real(rp) :: res
      real(rp), dimension(64) :: x, w

      integer :: ia_glob

      ! Find the Gauss Legendre roots and weights
      call gauss_legendre(64, 0.0_RP, 1.0_RP, x, w)

      call this%en%e_mesh()

      do i = 1, this%en%channels_ldos + 10
         if ((this%en%ene(i) - this%en%fermi) .le. 0.000001D0) fermi_point = i
      end do

      write (*, *) this%en%fermi, fermi_point
      this%gij_eta = 0.0D0; this%gji_eta = 0.0D0; this%ginmag_eta = 0.0D0; this%giz_eta = 0.0D0; this%giy_eta = 0.0D0; this%gix_eta = 0.0D0
      this%gjnmag_eta = 0.0D0; this%gjx_eta = 0.0D0; this%gjy_eta = 0.0D0; this%gjz_eta = 0.0D0

      if (this%control%recur == 'block') call this%recursion%zsqr()

      !do ia=1,this%lattice%njij
      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)
         ja_temp = (ia - 1) * 4 + 1
         y = (0.0_RP, 0.0_RP)
         do i = 1, 64
            eta = (0.0_RP, 0.0_RP)
            g0_ef = (0.0_RP, 0.0_RP)
            res = (1 - x(i)) / x(i)
            eta = cmplx(0.0_RP, res)
            select case (this%control%recur)
            case ('block')
               call this%block_green_ij_eta(ja_temp, eta, fermi_point, g0_ef)
            case ('chebyshev')
               call this%chebyshev_green_ij_eta(ja_temp, eta, fermi_point, g0_ef)
            end select
            y(i, :, :, :) = g0_ef(:, :, :)
         end do
         this%gij_eta(:, :, :, ia) = y(:, :, :, 1) - y(:, :, :, 2) + (1.0D0 / i_unit * y(:, :, :, 3) - 1.0D0 / i_unit * y(:, :, :, 4))
         this%gji_eta(:, :, :, ia) = y(:, :, :, 1) - y(:, :, :, 2) - (1.0D0 / i_unit * y(:, :, :, 3) - 1.0D0 / i_unit * y(:, :, :, 4))

         this%gij_eta(:, :, :, ia) = this%gij_eta(:, :, :, ia) * 0.5D0
         this%gji_eta(:, :, :, ia) = this%gji_eta(:, :, :, ia) * 0.5D0
         do i = 1, 9
            do j = 1, 9
               this%Ginmag_eta(:, j, i, iA) = this%Ginmag_eta(:, j, i, iA) + (this%gij_eta(:, j, i, ia) + this%gij_eta(:, j + 9, i + 9, ia)) * 0.5D0
               this%Giz_eta(:, j, i, iA) = this%Giz_eta(:, j, i, iA) + 0.5D0 * (this%gij_eta(:, j, i, ia) - this%gij_eta(:, j + 9, i + 9, ia))        !+Ginmag(k,j,i,iIA)*0.5d0
               this%Giy_eta(:, j, i, iA) = this%Giy_eta(:, j, i, iA) + 0.5D0 * (i_unit * this%gij_eta(:, j, i + 9, ia) - i_unit * this%gij_eta(:, j + 9, i, ia))  !+Ginmag(k,j,i,iIA)*0.5d0
               this%Gix_eta(:, j, i, iA) = this%Gix_eta(:, j, i, iA) + 0.5D0 * (this%gij_eta(:, j, i + 9, ia) + this%gij_eta(:, j + 9, i, ia))        !+Ginmag(k,j,i,iIA)*0.5d0
               !
               this%Gjnmag_eta(:, j, i, iA) = this%Gjnmag_eta(:, j, i, iA) + (this%gji_eta(:, j, i, ia) + this%gji_eta(:, j + 9, i + 9, ia)) * 0.5D0
               this%Gjz_eta(:, j, i, iA) = this%Gjz_eta(:, j, i, iA) + 0.5D0 * (this%gji_eta(:, j, i, ia) - this%gji_eta(:, j + 9, i + 9, ia))         !+Gjnmag(k,j,i,iIA)*0.5d0
               this%Gjy_eta(:, j, i, iA) = this%Gjy_eta(:, j, i, iA) + 0.5D0 * (i_unit * this%gji_eta(:, j, i + 9, ia) - i_unit * this%gji_eta(:, j + 9, i, ia))   !+Gjnmag(k,j,i,iIA)*0.5d0
               this%Gjx_eta(:, j, i, iA) = this%Gjx_eta(:, j, i, iA) + 0.5D0 * (this%gji_eta(:, j, i + 9, ia) + this%gji_eta(:, j + 9, i, ia))         !+Gjnmag(k,j,i,iIA)*0.5d0
            end do
         end do
      end do
   end subroutine calculate_intersite_gf_eta

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
      !complex(rp), dimension(18,18,this%lattice%nrec), intent(inout) :: g_ef
      complex(rp), dimension(18, 18, 1, atoms_per_process), intent(inout) :: g_ef
      !
      !
      integer :: nw
      integer :: ll
      integer :: ldim
      real(rp), dimension(this%lattice%nrec) :: a_inf0, b_inf0
      real(rp), dimension(18, 18, this%lattice%nrec) :: a_inf, b_inf
      integer :: n, n_glob
      !
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      ldim = 18

      ! Unchanged code
      nw = 10 * ll

      ! Calculate terminator coefficients
      call this%recursion%get_terminf(this%recursion%a_b, this%recursion%b2_b, atoms_per_process, &
                                      ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)

      do n_glob = start_atom, end_atom
         n = g2l_map(n_glob)

         call this%bgreen(g_ef(:, :, :, n), n, 1, this%en%channels_ldos + 10, a_inf(:, :, n), b_inf(:, :, n), eta)

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
      real(rp), dimension(18, 18, this%lattice%nrec) :: a_inf, b_inf
      integer :: n, n_glob
      !
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      ldim = 18
      !na = this%lattice%nrec
      eta = (0.0D0, 0.0D0)

      ! Unchanged code
      nw = 10 * ll

      ! Calculate terminator coefficients
      call this%recursion%get_terminf(this%recursion%a_b, this%recursion%b2_b, atoms_per_process, &
                                      ll, ldim, nw, a_inf, b_inf, a_inf0, b_inf0)

      do n_glob = start_atom, end_atom
         n = g2l_map(n_glob)

         call this%bgreen(this%g0(:, :, :, n), n, 1, this%en%channels_ldos + 10, a_inf(:, :, n), b_inf(:, :, n), eta)
      end do

   end subroutine block_green

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
      integer :: ia, mdir, nw, ll_t, ie, j
      real(rp), dimension(this%en%channels_ldos + 10, this%lattice%ntype) :: dx, dy, dz
      real(rp), dimension(18, this%en%channels_ldos + 10) :: doso

      complex(rp) :: dfac, impi

      complex(rp), dimension(3) :: gmask, lmask
      complex(rp), dimension(2, 3) :: gfac
      integer, dimension(4, 3) :: goff




      integer :: ia_glob

      impi = (pi, 0.0D0)
      gfac(1, 1) = 1.0D0; gfac(1, 2) = -i_unit; gfac(1, 3) = 1.0D0
      gfac(2, 1) = 1.0D0; gfac(2, 2) = i_unit; gfac(2, 3) = -1.0D0
      gmask(1) = 0.0D0; gmask(2) = 0.0D0; gmask(3) = 1.0D0
      !   if(lrot) then
      !      lmask(1)=0.0d0;lmask(2)=0.0d0;lmask(3)=1.0d0
      !   else
      lmask(1) = 1.0D0 / 3.0D0; lmask(2) = 1.0D0 / 3.0D0; lmask(3) = 1.0D0 / 3.0D0
      !   end if

      goff(1, 1) = 0; goff(2, 1) = 9; goff(1, 2) = 0; goff(2, 2) = 9; goff(1, 3) = 0; goff(2, 3) = 0
      goff(3, 1) = 9; goff(4, 1) = 0; goff(3, 2) = 9; goff(4, 2) = 0; goff(3, 3) = 9; goff(4, 3) = 9
      dfac = i_unit * impi / (2.0D0, 0.0D0)
      dx = 0; dy = 0; dz = 0
      nw = 10 * this%lattice%ntype * this%control%lld
      ll_t = this%control%lld

      doso = 0.0D0
      this%g0 = 0.0D0
      !do ia = 1, this%lattice%nrec
      do ia_glob = start_atom, end_atom
         ia = g2l_map(ia_glob)
         do mdir = 1, this%control%nmdir
            doso = 0.0D0
            call this%dos%density(doso, ia, mdir)
            if (this%control%nmdir == 1) then
               do ie = 1, this%en%channels_ldos + 10
                  do j = 1, 18
                     this%g0(j, j, ie, ia) = -i_unit * doso(j, ie) * impi
                  end do
                  write (300 + ia, *) this%en%ene(ie), sum(doso(1:9, ie)), sum(doso(10:18, ie))
               end do
            else
               do ie = 1, this%en%channels_ldos + 10
                  do j = 1, 9
                     ! Charge, from main direction.. (not z-component)
                     this%g0(j, j, ie, ia) = this%g0(j, j, ie, ia) - (doso(j, ie) + doso(j + 9, ie)) * dfac * lmask(mdir)!*1.0d0 /3.0d0 !mom(ja, mdir)**2
                     this%g0(j + 9, j + 9, ie, ia) = this%g0(j + 9, j + 9, ie, ia) - (doso(j, ie) + doso(j + 9, ie)) * dfac * lmask(mdir)!*1.0d0 /3.0d0 !mom(ja, mdir)**2
                     ! Spin dependent part
                     this%g0(j + goff(1, mdir), j + goff(2, mdir), ie, ia) = &
                        this%g0(j + goff(1, mdir), j + goff(2, mdir), ie, ia) - (doso(j, ie) - doso(j + 9, ie)) * gfac(1, mdir) * dfac

                     this%g0(j + goff(3, mdir), j + goff(4, mdir), ie, ia) = &
                        this%g0(j + goff(3, mdir), j + goff(4, mdir), ie, ia) - (doso(j, ie) - doso(j + 9, ie)) * gfac(2, mdir) * dfac
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

      allocate (this%Ginmag(9, 9, 64, atoms_per_process))
      allocate (this%Gjnmag(9, 9, 64, atoms_per_process))
      allocate (this%Gix(9, 9, 64, atoms_per_process))
      allocate (this%Gjx(9, 9, 64, atoms_per_process))
      allocate (this%Giy(9, 9, 64, atoms_per_process))
      allocate (this%Gjy(9, 9, 64, atoms_per_process))
      allocate (this%Giz(9, 9, 64, atoms_per_process))
      allocate (this%Gjz(9, 9, 64, atoms_per_process))

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
   !> It uses the fact that the phyisical GF's are invariant with respect to
   !> screening constants (Turek's book, page 72-73).
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
               do m = 1, 2 * l + 1
                  mls = l * l + m + ((lmaxi + 1)**2) * (s - 1) ! composed index
                  ! Sqrt(Delta) matrix of atom i
                  cdelta_i(mls, mls) = cmplx(this%symbolic_atom(this%lattice%iz(atom_i))%potential%dele(l, s), 0.0_RP)
               end do
            end do
         end do
         do s = 1, 2
            do l = 0, lmaxj ! Transform delta_j in complex
               do m = 1, 2 * l + 1
                  mls = l * l + m + ((lmaxj + 1)**2) * (s - 1) ! composed index
                  ! Sqrt(Delta) matrix of atom j
                  cdelta_j(mls, mls) = cmplx(this%symbolic_atom(this%lattice%iz(atom_j))%potential%dele(l, s), 0.0_RP)
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
   !> Given an auxiliary GF's (in the LMTO language), the function
   !> returns the correspondent auxiliary GF's transformed from
   !> one representation (alpha - in) to another (beta - out).
   !> The different representations are here denoted by their
   !> set of screening constants.
   !> See Turek's book, Eq. (3.57)
   !> Implemented by Ivan Miranda on 25.10.2023
   !---------------------------------------------------------------------------
   subroutine transform_auxiliary_gij(this, pmat_in_atom_i, pmat_out_atom_i, pmat_in_atom_j, pmat_out_atom_j, &
                                      aux_gij_in, aux_gij_out, screening_in, screening_out, atom_i, atom_j)
      !
      class(green) :: this
      !
      integer, intent(in) :: atom_i, atom_j ! Input atoms i and j
      complex(rp), intent(in) :: aux_gij_in(:, :, :) ! Auxiliary GF's (in)
      complex(rp), intent(inout) :: aux_gij_out(:, :, :) ! Auxiliary GF's (out), in another representation
      complex(rp), intent(in) :: pmat_in_atom_i(:, :, :) ! Potential P matrix (in) - in the same representation as GF's (in) for the atom i
      complex(rp), intent(in) :: pmat_out_atom_i(:, :, :) ! Potential P matrix (out) - in the same representation as GF's (out) for the atom i
      complex(rp), intent(in) :: pmat_in_atom_j(:, :, :) ! Potential P matrix (in) - in the same representation as GF's (in) for the atom j
      complex(rp), intent(in) :: pmat_out_atom_j(:, :, :) ! Potential P matrix (out) - in the same representation as GF's (out) for the atom j
      real(rp), intent(in) :: screening_in(0:, :) ! Screening constants relative to the GF's (in)
      real(rp), intent(in) :: screening_out(0:, :) ! Screening constants relative to the GF's (out)
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
            do m = 1, 2 * l + 1
               mls = l * l + m + ((size(screening_in, 1))**2) * (s - 1) ! Composed diagonal index
               ! The screening constants can be of any atoms (i or j), because they will act only
               ! when i = j.
               temp1 = cmplx(screening_in(l, s), 0.0_RP) ! Screening constants (in)
               temp2 = cmplx(screening_out(l, s), 0.0_RP) ! Screening constants (out)
               do ie = 1, size(pmat_in_atom_i, 3) ! Energy channel
                  pmat_resc1(mls, mls, ie) = pmat_in_atom_i(mls, mls, ie) / pmat_out_atom_i(mls, mls, ie)
                  pmat_resc2(mls, mls, ie) = pmat_in_atom_j(mls, mls, ie) / pmat_out_atom_j(mls, mls, ie)
                  ! For the pmat_resc3 matrix, does not matter to consider the P potential function
                  ! for atom i or j, because it's only acting when i = j
                  if (atom_i .eq. atom_j) then
                     pmat_resc3(mls, mls, ie) = (temp2 - temp1) * (pmat_in_atom_i(mls, mls, ie) / pmat_out_atom_i(mls, mls, ie))
                  end if
               end do
            end do
         end do
      end do

      ! Now do the rescaling of the auxiliary GF's (in)

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
      real(rp) :: a, b
      complex(rp) :: exp_factor
      integer :: ie, i, l, m, n

      this%g0 = 0.0D0

      allocate (kernel(this%control%lld * 2 + 2), polycheb(this%en%channels_ldos + 10, 0:this%control%lld * 2 + 2), w(this%en%channels_ldos + 10), &
                wscale(this%en%channels_ldos + 10))
      ! Defining rescaling coeficients
      a = (this%en%energy_max - this%en%energy_min) / (2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min) / 2

      wscale(:) = (this%en%ene(:) - b) / a

      ! Calculating the Jackson Kernel
      call jackson_kernel((this%control%lld) * 2 + 2, kernel)

      ! Calculating the Lorentz Kernel
!    call lorentz_kernel(this%control%lld, kernel, 4.0d0)

      do n = 1, 4 ! Loop on the number of on-site GFs to calculate the inter-site GFs
         ! Multiply the moments with the kernel
         do l = 1, 18
            do m = 1, 18
               this%recursion%mu_ng(l, m, :, n + istart - 1) = this%recursion%mu_n(l, m, :, n + istart - 1) * kernel(:)
            end do
         end do
         this%recursion%mu_ng(:, :, 2:size(kernel), n + istart - 1) = this%recursion%mu_ng(:, :, 2:size(kernel), n + istart - 1) * 2.0_RP

         ! Calculate the Chebyshev polynomials
         call t_polynomial(size(w), size(kernel), wscale(:), polycheb)

         ! Calculate the density of states
         !$omp parallel do default(shared) private(ie, i, exp_factor, l,m)
         do ie = 1, this%en%channels_ldos + 10
            do i = 1, size(kernel)
               exp_factor = -i_unit * exp(-i_unit * (i - 1) * acos(wscale(ie)))
               do l = 1, 18
                  do m = 1, 18
                     this%g0(l, m, ie, n) = this%g0(l, m, ie, n) + this%recursion%mu_ng(l, m, i, n + istart - 1) * exp_factor
                  end do
               end do
            end do
            do l = 1, 18
               do m = 1, 18
                  this%g0(l, m, ie, n) = this%g0(l, m, ie, n) / ((sqrt((a**2) - ((this%en%ene(ie) - b)**2))))
               end do
            end do
         end do
         !$omp end parallel do
      end do  ! End loop on n

      deallocate (kernel, polycheb, w, wscale)
   end subroutine chebyshev_green_ij

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments form the Chebyshev recursions to calculate the onsite GF
   !> using a small complex eta into the energy channel
   !---------------------------------------------------------------------------
   subroutine chebyshev_green_ij_eta(this, istart, eta, fermi_point, g_ef)
      class(green), intent(inout) :: this
      integer, intent(in) :: istart
      complex(rp), dimension(18, 18, 4), intent(inout) :: g_ef
      complex(rp), intent(in) :: eta
      integer, intent(in) :: fermi_point
      ! Local variables
      real(rp), dimension(:), allocatable :: kernel
      real(rp), dimension(:, :), allocatable :: polycheb
      real(rp), dimension(:), allocatable :: w, wscale
      real(rp) :: a, b
      complex(rp) :: exp_factor
      integer :: ie, i, l, m, n

      g_ef = 0.0D0

      allocate (kernel(this%control%lld * 2 + 2), polycheb(this%en%channels_ldos + 10, 0:this%control%lld * 2 + 2), w(this%en%channels_ldos + 10), &
                wscale(this%en%channels_ldos + 10))
      ! Defining rescaling coeficients
      a = (this%en%energy_max - this%en%energy_min) / (2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min) / 2

      wscale(:) = (this%en%ene(:) - b) / a

      ! Calculating the Jackson Kernel
      call jackson_kernel((this%control%lld) * 2 + 2, kernel)

      ! Calculating the Lorentz Kernel
!    call lorentz_kernel(this%control%lld, kernel, 4.0d0)

      do n = 1, 4 ! Loop on the number of on-site GFs to calculate the inter-site GFs
         ! Multiply the moments with the kernel
         do l = 1, 18
            do m = 1, 18
               this%recursion%mu_ng(l, m, :, n + istart - 1) = this%recursion%mu_n(l, m, :, n + istart - 1) * kernel(:)
            end do
         end do
         this%recursion%mu_ng(:, :, 2:size(kernel), n + istart - 1) = this%recursion%mu_ng(:, :, 2:size(kernel), n + istart - 1) * 2.0_RP

         ! Calculate the Chebyshev polynomials
         call t_polynomial(size(w), size(kernel), wscale(:), polycheb)

         ! Calculate the density of states
         !$omp parallel do default(shared) private(ie, i, exp_factor, l,m)
         do ie = fermi_point, fermi_point
            do i = 1, size(kernel)
               exp_factor = -i_unit * exp(-i_unit * (i - 1) * acos(((this%en%ene(ie) + eta) - b) / a))
               do l = 1, 18
                  do m = 1, 18
                     g_ef(l, m, n) = g_ef(l, m, n) + this%recursion%mu_ng(l, m, i, n + istart - 1) * exp_factor
                  end do
               end do
            end do
            do l = 1, 18
               do m = 1, 18
                  g_ef(l, m, n) = g_ef(l, m, n) / ((sqrt((a**2) - (((this%en%ene(ie) + eta) - b)**2))))
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
      real(rp), dimension(this%control%lld*2 + 2) :: kernel
      real(rp), dimension(this%en%channels_ldos + 10, 0:this%control%lld*2 + 2) :: polycheb
      real(rp), dimension(this%en%channels_ldos + 10) :: w, wscale
      real(rp) :: a, b
      complex(rp) :: exp_factor
      integer :: ie, i, l, m, n, nv

      integer :: n_glob

      this%g0 = 0.0D0

      ! Defining rescaling coeficients
      a = (this%en%energy_max - this%en%energy_min) / (2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min) / 2

      wscale(:) = (this%en%ene(:) - b) / a

      ! Number of DOS points
      nv = this%en%channels_ldos + 10

      ! Calculating the Jackson Kernel
      call jackson_kernel((this%control%lld) * 2 + 2, kernel)

      ! Calculating the Lorentz Kernel
      ! call lorentz_kernel(this%control%lld, kernel, 4.0d0)

      ! Determine how many atoms each process should handle
      ! Already done earlier
      !call get_mpi_variables(rank,this%lattice%nrec)

      do n_glob = start_atom, end_atom ! Loop on self-consistent atoms
         n = g2l_map(n_glob)
         ! Multiply the moments with the kernel
         do l = 1, 18
            do m = 1, 18
               this%recursion%mu_ng(l, m, :, n) = this%recursion%mu_n(l, m, :, n) * kernel(:)
            end do
         end do

         this%recursion%mu_ng(:, :, 2:size(kernel), n) = this%recursion%mu_ng(:, :, 2:size(kernel), n) * 2.0_RP

         !do i=1, size(kernel)
         !  write(400+n, *) i, sum(this%recursion%mu_n(n, i, 1:18, 1:18))
         !end do

         ! Calculate the Chebyshev polynomials
         call t_polynomial(size(w), size(kernel), wscale(:), polycheb)

         ! Calculate the density of states
         !$omp parallel do default(shared) private(ie, i, exp_factor, l,m)
         do ie = 1, this%en%channels_ldos + 10
            do i = 1, size(kernel)
               exp_factor = -i_unit * exp(-i_unit * (i - 1) * acos(wscale(ie)))
               do l = 1, 18
                  do m = 1, 18
                     this%g0(l, m, ie, n) = this%g0(l, m, ie, n) + this%recursion%mu_ng(l, m, i, n) * exp_factor
                  end do
               end do
            end do
            do l = 1, 18
               do m = 1, 18
                  this%g0(l, m, ie, n) = this%g0(l, m, ie, n) / ((sqrt((a**2) - ((this%en%ene(ie) - b)**2))))
               end do
            end do
         end do
         !$omp end parallel do
      end do  ! End loop on self-consistent atoms
!!! ! MPI moved to moment section
!!! #ifdef USE_MPI
!!!     call g_timer%start('MPI DOS communication')
!!!     call MPI_Allgather(this%g0(:,:,:,start_atom:end_atom), (end_atom-start_atom+1)*nv*18*18, MPI_DOUBLE_COMPLEX, &
!!!                        this%g0, (end_atom-start_atom+1)*nv*18*18, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
!!!     call g_timer%stop('MPI DOS communication')
!!! #endif
   end subroutine chebyshev_green

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Uses the moments from the Chebyshev recursions to calculate the onsite GF
   !> using a small complex eta into the energy channel
   !---------------------------------------------------------------------------
   subroutine chebyshev_green_eta(this, eta, fermi_point, g_ef)
      use mpi_mod
      class(green), intent(inout) :: this
      !complex(rp), dimension(18,18,this%lattice%nrec), intent(inout) :: g_ef
      complex(rp), dimension(18, 18, atoms_per_process), intent(inout) :: g_ef
      complex(rp), intent(in) :: eta
      integer, intent(in) :: fermi_point
      ! Local variables
      real(rp), dimension(this%control%lld*2 + 2) :: kernel
      real(rp), dimension(this%en%channels_ldos + 10, 0:this%control%lld*2 + 2) :: polycheb
      real(rp), dimension(this%en%channels_ldos + 10) :: w, wscale
      real(rp) :: a, b
      complex(rp) :: exp_factor
      integer :: ie, i, l, m, n, nv
      integer :: n_glob

      g_ef = 0.0D0

      ! Defining rescaling coeficients
      a = (this%en%energy_max - this%en%energy_min) / (2 - 0.3)
      b = (this%en%energy_max + this%en%energy_min) / 2

      wscale(:) = (this%en%ene(:) - b) / a

      ! Number of DOS points
      nv = this%en%channels_ldos + 10

      ! Calculating the Jackson Kernel
      call jackson_kernel((this%control%lld) * 2 + 2, kernel)

      ! Calculating the Lorentz Kernel
      ! call lorentz_kernel(this%control%lld, kernel, 4.0d0)

      do n_glob = start_atom, end_atom ! Loop on self-consistent atoms
         n = g2l_map(n_glob)
         ! Multiply the moments with the kernel
         do i = 1, 18
            this%recursion%mu_ng(i, i, :, n) = this%recursion%mu_n(i, i, :, n) * kernel(:)
         end do

         this%recursion%mu_ng(:, :, 2:size(kernel), n) = this%recursion%mu_ng(:, :, 2:size(kernel), n) * 2.0_RP

         !do i=1, size(kernel)
         !  write(400+n, *) i, sum(this%recursion%mu_n(n, i, 1:18, 1:18))
         !end do

         ! Calculate the Chebyshev polynomials
         call t_polynomial(size(w), size(kernel), wscale(:), polycheb)

         ! Calculate the density of states
         !$omp parallel do default(shared) private(ie, i, exp_factor, l,m)
         do ie = fermi_point, fermi_point
            do i = 1, size(kernel)
               exp_factor = -i_unit * exp(-i_unit * (i - 1) * acos(((this%en%ene(ie) + eta) - b) / a))
               do l = 1, 18
                  do m = 1, 18
                     g_ef(l, m, n) = g_ef(l, m, n) + this%recursion%mu_ng(l, m, i, n) * exp_factor
                  end do
               end do
            end do
            do l = 1, 18
               do m = 1, 18
                  g_ef(l, m, n) = g_ef(l, m, n) / ((sqrt((a**2) - (((this%en%ene(ie) + eta) - b)**2))))
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
      complex(rp), dimension(18, 18, ie_len), intent(inout) :: g_out
      real(rp), dimension(18, 18), intent(in) :: a_inf
      real(rp), dimension(18, 18), intent(in) :: b_inf
      complex(rp), intent(in) :: eta
      !
      integer :: nv, ldim, ll, ll_t, llinf
      real(rp) :: factor_z
      real(rp), dimension(this%en%channels_ldos + 10) :: e
      integer :: i, j, l, ei, info, ln, lwork
      integer, dimension(18) :: ipiv
      real(rp) :: etop, ebot, ea, eb

      complex(rp), dimension(18, 18) :: Q, Z, one, W, B2z, P
      complex(rp) :: zoff, cone, czero, ze, det, im
      complex(rp) :: coff
      complex(rp), dimension(18*18) :: work

      complex(rp), dimension(18, 18) :: Dfac_mat, Cshi_mat
      real(rp) :: a_diag, b_diag
      !
      integer, dimension(18) :: m_tab
      !

      !
      g_out = (0.0D0, 0.0D0)
      !
      ! Definitions so it is not necessary to change the code
      ll = this%control%lld
      nv = this%en%channels_ldos
      ldim = 18
      e = this%en%ene
      factor_z = 1.0D0

      ! Unchanged code
      ll_t = ll
      llinf = ll!+ll/2!0!3*ll!/2!+30 !300
      cone = (1.0D0, 0.0D0)
      czero = (0.0D0, 0.0D0)
      im = (0.0D0, 1.0D0)
      zoff = (0.0D0, 0.01D0) * 0.0D0
      coff = 0.000D0 * (0.0D0, 1.0D0)
      one = (0.0D0, 0.0D0)

      ! Initializing Shift and Scaling for rigid band shift
      Dfac_mat = (1.0D0, 0.0D0)
      Cshi_mat = (0.0D0, 0.0D0)

      do i = 1, ldim
         m_tab(i) = i
      end do

      do i = 1, ldim
         one(i, i) = (1.0D0, 0.0D0)
      end do
      do i = 1, ldim
         !a_diag= a_inf(i,i,n) !(a_inf(1,1,n) + a_inf(10,10,n))*0.5d0
         !b_diag= b_inf(i,i,n) !(b_inf(1,1,n) + b_inf(10,10,n))*0.5d0
         a_diag = (a_inf(1, 1) + a_inf(10, 10)) * 0.5D0
         b_diag = (b_inf(1, 1) + b_inf(10, 10)) * 0.5D0
      end do
      ! Standard semi-circle construction of Greens functions
      !write(12334,'(18f8.4)') a_inf
      !write(12335,'(18f8.4)') b_inf
      !write(12336,'(18f8.4)') a_inf
      !write(12337,'(18f8.4)') this%recursion%a_b(:,:,2,i_site)
      !write(12337,'(18f8.4)') this%recursion%b2_b(:,:,2,i_site)
    !!!$omp parallel do default(shared) private(ei, Z, Ze, Q, B2z, i, ea, eb, det, zoff, l, ln, j, P , ipiv, info, work, lwork, W)
      !$omp parallel do default(shared) private(ei, Z, Ze, Q, B2z, i, etop, ebot, ea, eb, det, zoff, l, ln, j, P , ipiv, info, work, lwork, W)
      do ei = ie_start, ie_start + ie_len - 1
         Z = e(ei) * one
         ze = e(ei)!-zoff
         Q = (0.0D0, 0.0D0)
         B2z = 0.0D0
         if (this%control%sym_term) then
            do i = 1, ldim !  Orbital-independent
               etop = a_diag + 2.0D0 * b_diag
               ebot = a_diag - 2.0D0 * b_diag
               ea = e(ei) - etop
               eb = e(ei) - ebot
               det = ea * eb
               zoff = sqrt(det)
               Q(i, i) = (ze + eta - a_diag - zoff) * 0.5D0
            end do
         else
            do i = 1, ldim  ! Orbital-dependent
               if (i == 1 .or. i == 10) then  !(now the s-bands are broadened earlier)
                  etop = a_inf(i, i) - Cshi_mat(i, i) + 2 * b_inf(i, i) * 1.025D0 / Dfac_mat(i, i)
                  ebot = a_inf(i, i) - Cshi_mat(i, i) - 2 * b_inf(i, i) * 1.025D0 / Dfac_mat(i, i)
               else
                  etop = a_inf(i, i) - Cshi_mat(i, i) + 2 * b_inf(i, i)!*1.01d0/Dfac_mat(i,i)
                  ebot = a_inf(i, i) - Cshi_mat(i, i) - 2 * b_inf(i, i)!*1.01d0/Dfac_mat(i,i)
               end if
               ea = e(ei) / Dfac_mat(i, i) - etop - 1.0D0 * Cshi_mat(i, i)
               eb = e(ei) / Dfac_mat(i, i) - ebot - 1.0D0 * Cshi_mat(i, i)
               det = ea * eb
               zoff = sqrt(det)!*0.5d0
               zoff = zoff * (1.0)**ei
           !!! if(nv<=10) zoff=zoff*1.0d-3   ! Hack for Efermi evaluation
               ! Below for orbital dependent terminator
               Q(i, i) = ((e(ei) + eta) / Dfac_mat(i, i) - 1.0D0 * Cshi_mat(i, i) - a_inf(i, i) - real(zoff) - aimag(zoff) * factor_z * im) * 0.5D0
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
                  if ((abs(real(Q(i, j))) .lt. 10**(-12)) .and. (abs(aimag(Q(i, j))) .lt. 10**(-12))) then
                     Q(i, j) = (0.0D0, 0.0D0)
                  end if
                  Q(i, j) = P(i, j) / Dfac_mat(i, j) - Cshi_mat(i, j) - cone * this%recursion%a_b(i, j, ln, i_site) - Q(i, j)
                  B2z(i, j) = cone * this%recursion%b2_b(i, j, ln, i_site)
               end do
            end do
            ! LU factorization
            call zgetrf(ldim, ldim, Q, ldim, ipiv, info)
            ! Inverse of Q from Chol. fac.
            lwork = ldim * ldim
            call zgetri(ldim, Q, ldim, ipiv, work, lwork, info)
            call zgemm('n', 'n', ldim, ldim, ldim, cone, Q, ldim, B2z, ldim, czero, W, ldim)
            call zgemm('c', 'n', ldim, ldim, ldim, cone, B2z, ldim, W, ldim, czero, Q, ldim)
         end do
         !
         do j = 1, ldim
            do i = 1, ldim
               !bdos(ei,m_tab(i),m_tab(j),n)=bdos(ei,m_tab(i),m_tab(j),n) + abs(-aimag(Q(i,j))/3.14159265359d0)/Dfac_mat(i,j)
               if (aimag(eta) .eq. 0) then
                  g_out(m_tab(i), m_tab(j), ei) = g_out(m_tab(i), m_tab(j), ei) + &
                                                  real(Q(i, j) / Dfac_mat(i, j)) + aimag(Q(i, j) / Dfac_mat(i, j)) * (1.0D0)**ei * (0.0D0, 1.0D0)
               else
                  g_out(m_tab(i), m_tab(j), ei) = g_out(m_tab(i), m_tab(j), ei) + (Q(i, j) / Dfac_mat(i, j))
               end if
            end do
         end do
         !write(12333,'(19f8.4)') e(ei), (real(g_out(i,i,ei)),i=1,18)
      end do
      !$omp end parallel do

   end subroutine bgreen

end module
