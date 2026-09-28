program test_rotation_main
   use test_rotation_finite_q_torque_mod, only: run_test_rotation_finite_q_torque
   use test_rotation_production_adapter_mod, only: run_test_rotation_production_adapter
   use test_rotation_second_order_mod, only: run_test_rotation_second_order
   use test_rotation_static_contour_mod, only: run_test_rotation_static_contour
   use test_rotation_static_lkag_mod, only: run_test_rotation_static_lkag
   use test_rotation_torque_hessian_mod, only: run_test_rotation_torque_hessian
   implicit none

   character(len=128) :: case_name

   call get_command_argument(1, case_name)
   select case (trim(case_name))
   case ('rotation_finite_q_torque')
      call run_test_rotation_finite_q_torque()
   case ('rotation_production_adapter')
      call run_test_rotation_production_adapter()
   case ('rotation_second_order')
      call run_test_rotation_second_order()
   case ('rotation_static_contour')
      call run_test_rotation_static_contour()
   case ('rotation_static_lkag')
      call run_test_rotation_static_lkag()
   case ('rotation_torque_hessian')
      call run_test_rotation_torque_hessian()
   case default
      error stop 'unknown LR rotation unit-test case'
   end select
end program test_rotation_main
