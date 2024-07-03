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
      procedure :: save_density_matrix_to_file
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
      integer :: k, l, m1, m2, m3, m4, count

      do l = 0, 10
         print *, '-----------------------------------------------------------------------------------------------------'
         print *, 'l = ', l
         print *, '-----------------------------------------------------------------------------------------------------'
         do k = 0, 50
            print *, 'k = ', k
            count = 0
            do m1 = -l, l
               do m2 = -l, l
                  do m3 = -l, l
                     do m4 = -l, l
                        ak_test = a_k(k,l,m1,m2,m3,m4)
                        if (ak_test .ne. 0.0d0) then
                           count = count + 1
                        end if
                        ! print *, 'k = ', k, ' l = ', l, ' m1 = ', m1, ' m2 = ', m2, ' m3 = ', m3, ' m4 = ', m4
                        ! print *, 'a_k = ', ak_test
                     end do
                  end do
               end do
            end do
            print *, 'Nmb of non-zero aks: ', count
            print *, '-----------------------------------------------------------------------------------------------------'
         end do
         
      end do

      stop
   end subroutine

   !> Implemented for single atom (i.e choose na) and hubbard_orb = spd
   subroutine save_density_matrix_to_file(this, matrix, channels_ldos)
      class(hubbard_u) :: this
      real(rp), intent(in) :: matrix(:,:,:,:)
      integer, intent(in) :: channels_ldos
      character(len=30) :: filename
      integer :: i, j, l
      integer :: n, m, lsize, ssize
      
      write(filename, '(A,I0)') 'ld_matrix_ch_ldos_', channels_ldos 
      
      lsize = size(matrix, 1)
      ssize = size(matrix, 2)
      n = size(matrix, 3)
      m = size(matrix, 4)
      
      ! Save the matrix to a file
      ! s-orbital
      open(unit=10, file=trim(filename) // '_orb_s_spin_u', status='replace')
      do i = 1, n
         write(10, *) (matrix(1, 1, i, j), j = 1, m)
      end do
      close(10)

      open(unit=10, file=trim(filename) // '_orb_s_spin_d', status='replace')
      do i = 1, n
         write(10, *) (matrix(1, 2, i, j), j = 1, m)
      end do
      close(10)
      !-------------------------------------------------------------------------
      ! p-orbital
      open(unit=10, file=trim(filename) // '_orb_p_spin_u', status='replace')
      do i = 1, n
         write(10, *) (matrix(2, 1, i, j), j = 1, m)
      end do
      close(10)

      open(unit=10, file=trim(filename) // '_orb_p_spin_d', status='replace')
      do i = 1, n
         write(10, *) (matrix(2, 2, i, j), j = 1, m)
      end do
      close(10)
      !-------------------------------------------------------------------------
      ! d-orbital
      open(unit=10, file=trim(filename) // '_orb_d_spin_u', status='replace')
      do i = 1, n
         write(10, *) (matrix(3, 1, i, j), j = 1, m)
      end do
      close(10)

      open(unit=10, file=trim(filename) // '_orb_d_spin_d', status='replace')
      do i = 1, n
         write(10, *) (matrix(3, 2, i, j), j = 1, m)
      end do
      close(10)
  
      print *, 'Matrix data saved to ', trim(filename)
  end subroutine save_density_matrix_to_file

end module hubbard_u_mod