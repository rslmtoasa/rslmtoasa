!> Global parameters
!> /note Keep the contents of this module to a minimum
module Parameters

   implicit none

   integer, parameter :: snglprec = selected_real_kind(6, 37)  !< define precision for single reals
   !integer, parameter :: dblprec = selected_real_kind(6, 37)  !< define precision for double reals
   integer, parameter :: dblprec = selected_real_kind(15, 307)  !< define precision for double reals
   integer, parameter :: qdprec = selected_real_kind(33, 4931)  !< define precision for quad reals
   integer, parameter :: long = selected_int_kind(range(1) * 2)   !< define long integer
   real(dblprec) :: dbl_tolerance = 1E-14 !!! Note this is used for both ne and eq as > or < comparisons

   integer, parameter :: ifileno = 55 !< File handle number for input files
   integer, parameter :: ofileno = 66 !< File handle number for output files
   integer, parameter :: mfileno = 77 !< File handle number for memory log file

   integer, parameter :: ofileno2 = 67 !< File handle number for output files
   integer, parameter :: ofileno3 = 68 !< File handle number for output files

   public

end module
