!------------------------------------------------------------------------------
! RS-LMTO-ASA -- unit test
!------------------------------------------------------------------------------
!
!> @brief Verify the default state of the two-index Green-function arrays.
!> @details This uses only the public lifecycle contract with small dummy
!>          lattice and energy objects; no physical system is constructed.
!------------------------------------------------------------------------------
program test_green_lifecycle
   use basis_mod, only: basis_init
   use energy_mod, only: energy
   use green_mod, only: green
   use lattice_mod, only: lattice
   use mpi_mod, only: atoms_per_process
   use precision_mod, only: rp
   implicit none

   type(lattice), target :: lattice_obj
   type(energy), target :: energy_obj
   type(green) :: green_obj
   logical :: failed
   complex(rp), parameter :: zero = (0.0_rp, 0.0_rp)

   failed = .false.
   call basis_init(1)
   atoms_per_process = 1
   lattice_obj%njij = 0
   lattice_obj%nrec = 1
   energy_obj%channels_ldos = 0

   green_obj%lattice => lattice_obj
   green_obj%en => energy_obj
   call green_obj%restore_to_default()

   call require(allocated(green_obj%g00ij), 'g00ij is allocated')
   call require(allocated(green_obj%g00ji), 'g00ji is allocated')
   call require(allocated(green_obj%g01ij), 'g01ij is allocated')
   call require(allocated(green_obj%g01ji), 'g01ji is allocated')
   call require(all(green_obj%g00ij == zero), 'g00ij is zero')
   call require(all(green_obj%g00ji == zero), 'g00ji is zero')
   call require(all(green_obj%g01ij == zero), 'g01ij is zero')
   call require(all(green_obj%g01ji == zero), 'g01ji is zero')

   if (failed) then
      write (*, '(a)') 'RESULT: FAIL'
      error stop 1
   else
      write (*, '(a)') 'RESULT: PASS'
   end if

contains

   subroutine require(condition, description)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: description

      if (.not. condition) then
         write (*, '(a,a)') 'FAIL: ', trim(description)
         failed = .true.
      end if
   end subroutine require

end program test_green_lifecycle
