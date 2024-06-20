!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: MPI
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
!> Module to handle the MPI variables for parallelization
!------------------------------------------------------------------------------
module mpi_mod
#ifdef USE_MPI
   use mpi
#endif
   integer :: ierr, rank, numprocs
   integer :: remainder_atoms, atoms_per_process, start_atom, end_atom
   integer, dimension(:), allocatable :: l2g_map
   integer, dimension(:), allocatable :: g2l_map

contains

   subroutine get_mpi_variables(rank_in, number_of_atoms)
      implicit none
      integer, intent(in) :: number_of_atoms, rank_in

#ifdef USE_MPI
      atoms_per_process = number_of_atoms / numprocs
      remainder_atoms = mod(number_of_atoms, numprocs)

      if (rank_in < remainder_atoms) then
         atoms_per_process = atoms_per_process + 1
         start_atom = rank_in * atoms_per_process + 1
      else
         start_atom = rank_in * atoms_per_process + remainder_atoms + 1
      end if
      end_atom = start_atom + atoms_per_process - 1
#else
      remainder_atoms = 0
      atoms_per_process = number_of_atoms
      start_atom = 1
      end_atom = number_of_atoms
#endif

      call get_mpi_mapping(rank_in, number_of_atoms)

      return

   end subroutine get_mpi_variables

   subroutine get_mpi_mapping(rank_in, number_of_atoms)
      implicit none
      integer, intent(in) :: number_of_atoms, rank_in

      integer :: i_g, i_l

      if (allocated(l2g_map)) deallocate (l2g_map)
      if (allocated(g2l_map)) deallocate (g2l_map)

      allocate (l2g_map(atoms_per_process))
      l2g_map = 0
      allocate (g2l_map(number_of_atoms))
      g2l_map = 0

      do i_g = start_atom, end_atom
         i_l = i_g - start_atom + 1
         l2g_map(i_l) = i_g
         g2l_map(i_g) = i_l
      end do

      return

   end subroutine get_mpi_mapping

end module mpi_mod
