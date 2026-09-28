program test_path_operator_main
   use test_path_operator_contour_mod, only: run_test_path_operator_contour
   use test_path_operator_fixed_z_mod, only: run_test_path_operator_fixed_z
   implicit none

   character(len=128) :: case_name

   call get_command_argument(1, case_name)
   select case (trim(case_name))
   case ('path_operator_contour')
      call run_test_path_operator_contour()
   case ('path_operator_fixed_z')
      call run_test_path_operator_fixed_z()
   case default
      error stop 'unknown LMTO path-operator unit-test case'
   end select
end program test_path_operator_main
