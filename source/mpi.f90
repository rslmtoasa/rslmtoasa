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
   use logger_mod, only: g_logger
#ifdef USE_MPI
   use mpi
#endif
   integer :: ierr, rank, numprocs
   integer :: remainder_atoms, atoms_per_process, start_atom, end_atom
   logical :: oversubscribe_warned = .false.
   integer, dimension(:), allocatable :: l2g_map
   integer, dimension(:), allocatable :: g2l_map

contains

   subroutine get_mpi_range(rank_in, n_items, start_idx, end_idx, items_per_process, l2g_map_out, g2l_map_out, region_tag)
      implicit none
      integer, intent(in) :: rank_in, n_items
      integer, intent(out) :: start_idx, end_idx, items_per_process
      integer, dimension(:), allocatable, intent(out), optional :: l2g_map_out
      integer, dimension(:), allocatable, intent(out), optional :: g2l_map_out
      character(len=*), intent(in), optional :: region_tag

      integer :: remainder_items, i_g, i_l
      character(len=160) :: warn_msg
      character(len=32) :: region_name

      region_name = 'work'
      if (present(region_tag)) region_name = trim(region_tag)

#ifdef USE_MPI
      if ((numprocs > n_items) .and. (rank_in == 0) .and. (.not. oversubscribe_warned)) then
         write (warn_msg, '(a,a,a,i0,a,i0,a)') 'MPI oversubscription in ', trim(region_name), &
            '-parallel region (ranks=', numprocs, ', items=', n_items, '). Some ranks will stay idle.'
         call g_logger%warning(trim(warn_msg), __FILE__, __LINE__)
         oversubscribe_warned = .true.
      end if

      items_per_process = n_items/numprocs
      remainder_items = mod(n_items, numprocs)

      if (rank_in < remainder_items) then
         items_per_process = items_per_process + 1
         start_idx = rank_in*items_per_process + 1
      else
         start_idx = rank_in*items_per_process + remainder_items + 1
      end if
      end_idx = start_idx + items_per_process - 1

      if (items_per_process == 0) then
         start_idx = 1
         end_idx = 0
      end if
#else
      items_per_process = n_items
      start_idx = 1
      end_idx = n_items
#endif

      if (present(l2g_map_out)) then
         if (allocated(l2g_map_out)) deallocate(l2g_map_out)
         allocate(l2g_map_out(items_per_process))
         l2g_map_out = 0
      end if
      if (present(g2l_map_out)) then
         if (allocated(g2l_map_out)) deallocate(g2l_map_out)
         allocate(g2l_map_out(n_items))
         g2l_map_out = 0
      end if

      do i_g = start_idx, end_idx
         i_l = i_g - start_idx + 1
         if (present(l2g_map_out)) l2g_map_out(i_l) = i_g
         if (present(g2l_map_out)) g2l_map_out(i_g) = i_l
      end do

   end subroutine get_mpi_range

   subroutine get_mpi_variables(rank_in, number_of_atoms)
      implicit none
      integer, intent(in) :: number_of_atoms, rank_in
      call get_mpi_range(rank_in, number_of_atoms, start_atom, end_atom, atoms_per_process, region_tag='atom')
      remainder_atoms = max(0, number_of_atoms - numprocs*(number_of_atoms/ max(1, numprocs)))
      call get_mpi_mapping(rank_in, number_of_atoms)

      return

   end subroutine get_mpi_variables

   subroutine get_mpi_mapping(rank_in, number_of_atoms)
      implicit none
      integer, intent(in) :: number_of_atoms, rank_in
      call get_mpi_range(rank_in, number_of_atoms, start_atom, end_atom, atoms_per_process, l2g_map, g2l_map, 'atom')

      return

   end subroutine get_mpi_mapping

end module mpi_mod
