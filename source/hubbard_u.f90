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
   use hamiltonian_mod
   use green_mod
#ifdef USE_SAFE_ALLOC
   use safe_alloc_mod, only: g_safe_alloc
#endif
   implicit none

   private

   !> Modules main structure
   type, public :: hubbard_u
      !> General variables
      !> Pointer to Greens object to have access to all information. Maybe change?
      class(green), pointer :: green
      !> Hamiltonian object
      class(hamiltonian), pointer :: hamiltonian
      !> Test
      !> Orbital name
      character(len=12) :: o_type
      !> Hubbard U
      real(rp) :: u

   contains
      procedure :: restore_to_default
      procedure :: calc_test
      final :: destructor
   end type

   interface hubbard_u
      procedure constructor
   end interface hubbard_u

contains

   function constructor(green_obj) result(obj)
      type(hubbard_u) :: obj
      type(green), target, intent(in) :: green_obj

      obj%green => green_obj
      obj%hamiltonian => green_obj%recursion%hamiltonian

      call obj%restore_to_default()
   end function constructor

   subroutine destructor(this)
      type(hubbard_u) :: this

   end subroutine destructor

   subroutine restore_to_default(this)
      class(hubbard_u) :: this

      this%o_type = 'orbital_test'
      this%u = 4.0d0

      print *, this%o_type
   
   end subroutine restore_to_default

   subroutine calc_test(this)
      class(hubbard_u) :: this
      real(rp) :: ak_test
      integer :: k, l, m1, m2, m3, m4

      do k = 0, 4
         do l = 0, 2
            do m1 = -l, l
               do m2 = -l, l
                  do m3 = -l, l
                     do m4 = -l, l
                        ak_test = a_k(k,l,m1,m2,m3,m4)
                        print *, 'k = ', k, ' l = ', l, ' m1 = ', m1, ' m2 = ', m2, ' m3 = ', m3, ' m4 = ', m4
                        print *, 'a_k = ', ak_test
                     end do
                  end do
               end do
            end do
         end do
      end do
                     
      stop
   end subroutine

   subroutine calc_density_matrix(this)
      class(hubbard_u) :: this
      real(rp), dimension(:,:,:,:,:), allocatable :: density_mat
      integer :: i, j, max_nmb_orb, nrec, m_max
      character(len=1) :: orb

      nrec = this%hamiltonian%lattice%nrec
      ! max_nmb_orb = size(this%hamiltonian%hubbard_u,2)
      m_max = 0
      orb = ''

      do i = 1, nrec
         
      end do
      ! Test
      ! do i = 1, nrec
      !    do j = 1, max_nmb_orb
      !       print *, this%hamiltonian%hubbard_orb(i)
      !       print *, this%hamiltonian%hubbard_orb(i)(j,j)
      !       orb = this%hamiltonian%hubbard_orb(i)(j,j)
      !       if (orb == '') then
      !          !nothing
      !          print *, 'No hubbard u!!'
      !       else if (orb == 's') then
      !          print *, 's'
      !       else if (orb == 'p') then
      !          print *, 'p'
      !          if (m_max == 0) then
      !             m_max = 1
      !          end if
      !       else if (orb == 'd') then
      !          print *, 'd'
      !          m_max = 2
      !       end if
      !    end do
      ! end do
      print *, 'Total m = ', 2*m_max + 1
      allocate (density_mat(size(this%hamiltonian%hubbard_u, 1), size(this%hamiltonian%hubbard_u, 2), 2, 2*m_max + 1, 2*m_max + 1))

   end subroutine

end module hubbard_u_mod