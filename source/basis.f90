!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Basis
!
!> @author
!> Anders Bergman
!
! DESCRIPTION:
!> Module that defines the orbital-basis dimensions used throughout the code.
!> All size parameters are set once at initialisation via basis_init and then
!> available as protected module-level variables.
!>
!> Supported bases (lmax value → number of orbitals per spin channel):
!>
!>   lmax=1  sp      norb= 4   nb= 8
!>   lmax=2  spd     norb= 9   nb=18   (default, fully supported)
!>   lmax=3  spdf    norb=16   nb=32
!>
!> For Bogoliubov-de Gennes (BdG / Nambu) mode the full block size is
!> nb = 4*norb instead of 2*norb.
!>
!> Usage
!> -----
!> Call basis_init() once, early in the program, before any arrays that depend
!> on nb or norb are allocated.  After that, simply USE this module and read
!> norb, nb, spin_off as needed.
!>
!>   use basis_mod, only: nb, norb, spin_off
!>   ...
!>   complex(rp) :: ham(nb, nb)   ! automatic array in a subroutine
!>   allocate(work(nb, nb, natom)) ! heap array
!------------------------------------------------------------------------------
module basis_mod

   use precision_mod, only: ip
   use logger_mod, only: g_logger
   use string_mod, only: int2str
   implicit none

   private

   !> Maximum angular momentum quantum number of the highest occupied shell.
   !> lmax=1 → sp, lmax=2 → spd, lmax=3 → spdf.
   integer, public, protected :: lmax_basis = 2

   !> Number of orbitals per spin channel: norb = (lmax_basis + 1)**2
   integer, public, protected :: norb = 9

   !> Full spinor block size: nb = 2*norb in normal mode, 4*norb in BdG mode.
   integer, public, protected :: nb = 18

   !> Spin-down block offset within the nb×nb spinor block.
   !> In normal mode spin_off = norb (spin-down block starts at index norb+1).
   integer, public, protected :: spin_off = 9

   !> BdG / Nambu-spinor mode flag.
   !> When .true., nb = 4*norb and the Hamiltonian lives in Nambu space.
   logical, public, protected :: bdg_mode = .false.

   public :: basis_init

contains

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief Initialise basis dimensions from the maximum angular momentum.
   !
   !> Call this once, early in the calculation driver, before any module that
   !> allocates arrays depending on nb or norb.  Typically called from
   !> control%build_from_file.
   !>
   !> @param[in] lmax   Maximum angular momentum (1=sp, 2=spd, 3=spdf)
   !> @param[in] bdg    (optional) Enable Nambu/BdG mode (default .false.)
   !---------------------------------------------------------------------------
   subroutine basis_init(lmax, bdg)
      integer, intent(in) :: lmax
      logical, intent(in), optional :: bdg

      logical :: bdg_

      bdg_ = .false.
      if (present(bdg)) bdg_ = bdg

      if (lmax < 1 .or. lmax > 3) then
         write (*, '(a,i0)') 'basis_init: unsupported lmax = ', lmax
         error stop 'basis_init: lmax must be 1 (sp), 2 (spd), or 3 (spdf)'
      end if

      lmax_basis = lmax
      norb       = (lmax + 1)**2
      spin_off   = norb

      if (bdg_) then
         nb       = 4*norb
         bdg_mode = .true.
      else
         nb       = 2*norb
         bdg_mode = .false.
      end if

      call g_logger%info('Basis initialised: lmax=' // int2str(lmax_basis) // &
         ' (' // basis_label(lmax_basis) // ')' // &
         '  norb=' // int2str(norb) // &
         '  nb=' // int2str(nb) // &
         '  spin_off=' // int2str(spin_off) // &
         '  bdg=' // merge('T', 'F', bdg_mode), __FILE__, __LINE__)

   end subroutine basis_init

   !---------------------------------------------------------------------------
   ! Returns a short orbital-basis label string for a given lmax.
   !---------------------------------------------------------------------------
   pure function basis_label(lmax) result(label)
      integer, intent(in) :: lmax
      character(len=4) :: label
      select case (lmax)
      case (1); label = 'sp  '
      case (2); label = 'spd '
      case (3); label = 'spdf'
      case default; label = '?   '
      end select
   end function basis_label

end module basis_mod
