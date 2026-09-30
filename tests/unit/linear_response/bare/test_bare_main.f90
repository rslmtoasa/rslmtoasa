program test_bare_main
   use test_eigenpair_gf_oracle_mod, only: run_test_eigenpair_gf_oracle
   use test_eigenpair_gf_quadrature_mod, only: run_test_eigenpair_gf_quadrature
   use test_ks_susceptibility_mod, only: run_test_ks_susceptibility
   use test_product_eigenpair_gf_oracle_mod, only: run_test_product_eigenpair_gf_oracle
   use test_product_ks_susceptibility_mod, only: run_test_product_ks_susceptibility
   use test_projected_chi0_mod, only: run_test_projected_chi0
   use test_juelich_d_reference_mod, only: run_test_juelich_d_reference
   use test_rs_gf_susceptibility_mod, only: run_test_rs_gf_susceptibility
   implicit none

   character(len=128) :: case_name

   call get_command_argument(1, case_name)
   select case (trim(case_name))
   case ('eigenpair_gf_oracle')
      call run_test_eigenpair_gf_oracle()
   case ('eigenpair_gf_quadrature')
      call run_test_eigenpair_gf_quadrature('')
   case ('eigenpair_gf_mixed')
      call run_test_eigenpair_gf_quadrature('mixed')
   case ('ks_susceptibility')
      call run_test_ks_susceptibility()
   case ('product_eigenpair_gf_oracle')
      call run_test_product_eigenpair_gf_oracle('')
   case ('product_eigenpair_gf_reject_integration_eta')
      call run_test_product_eigenpair_gf_oracle('reject_integration_eta')
   case ('product_ks_susceptibility')
      call run_test_product_ks_susceptibility()
   case ('projected_chi0')
      call run_test_projected_chi0()
   case ('juelich_d_reference')
      call run_test_juelich_d_reference()
   case ('rs_gf_susceptibility')
      call run_test_rs_gf_susceptibility()
   case default
      error stop 'unknown LR bare unit-test case'
   end select
end program test_bare_main
