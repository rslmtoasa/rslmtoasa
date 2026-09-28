program test_kernel_dyson_main
   use test_alsda_kernel_mod, only: run_test_alsda_kernel
   use test_compact_dyson_oracle_mod, only: run_test_compact_dyson_oracle
   use test_compact_gsr_action_mod, only: run_test_compact_gsr_action
   use test_compact_interaction_mod, only: run_test_compact_interaction
   use test_dyson_mod, only: run_test_dyson
   use test_lcmm_sumrule_mod, only: run_test_lcmm_sumrule
   use test_projected_interacting_mod, only: run_test_projected_interacting
   use test_projected_lcmm_mod, only: run_test_projected_lcmm
   implicit none

   character(len=128) :: case_name

   call get_command_argument(1, case_name)
   select case (trim(case_name))
   case ('alsda_kernel')
      call run_test_alsda_kernel('')
   case ('alsda_provenance')
      call run_test_alsda_kernel('provenance')
   case ('alsda_zero')
      call run_test_alsda_kernel('zero')
   case ('compact_dyson_oracle')
      call run_test_compact_dyson_oracle()
   case ('compact_gsr_action')
      call run_test_compact_gsr_action()
   case ('compact_interaction')
      call run_test_compact_interaction()
   case ('dyson')
      call run_test_dyson()
   case ('lcmm_sumrule')
      call run_test_lcmm_sumrule()
   case ('projected_interacting')
      call run_test_projected_interacting()
   case ('projected_lcmm')
      call run_test_projected_lcmm()
   case default
      error stop 'unknown LR kernel/Dyson unit-test case'
   end select
end program test_kernel_dyson_main
