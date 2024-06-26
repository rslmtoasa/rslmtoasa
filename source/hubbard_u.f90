!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Hubbard_u
!
!> @author
!> Viktor Frilén
!> Emil Beiersdorf
!
! DESCRIPTION:
!> Summerproject 2024 LDA+U, module for all implementations
!------------------------------------------------------------------------------

module hubbard_u_mod
   use energy_mod
   use control_mod
   use lattice_mod
   use symbolic_atom_mod
   use recursion_mod
   use density_of_states_mod
   use precision_mod, only : rp
   use math_mod
   use green_mod
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Modules main structure
   type, public :: hubbard_u
      !> General variables
      !> Test
      !> Orbital name
      character(len=12), allocatable :: o_type
      !> Hubbard U
      real(rp), allocatable :: u

   contains
      procedure :: restore_to_default
      final :: destructor
   end type

   interface hubbard_u
      procedure constructor
   end interface hubbard_u

contains

   function constructor() result(obj)
      type(hubbard_u) :: obj

      call obj%restore_to_default()
   end function constructor

   subroutine destructor(this)
      type(hubbard_u) :: this

      if(allocated(this%o_type)) deallocate(this%o_type)
      if(allocated(this%u)) deallocate(this%u)
   end subroutine destructor

   subroutine restore_to_default(this)
      class(hubbard_u) :: this

      allocate(this%o_type)
      allocate(this%u)

      this%o_type = 'orbital_test'
      this%u = 4.0d0

      print *, this%o_type
   
   end subroutine restore_to_default

end module hubbard_u_mod