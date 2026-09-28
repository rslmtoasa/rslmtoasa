program test_basis_main
   use test_endpoint_branch_action_mod, only: run_test_endpoint_branch_action
   use test_gf_endpoint_augmentation_mod, only: run_test_gf_endpoint_augmentation
   use test_pauli_transition_vertex_mod, only: run_test_pauli_transition_vertex
   use test_product_response_mod, only: run_test_product_response
   use test_product_response_basis_mod, only: run_test_product_response_basis
   use test_projected_site_spin_mod, only: run_test_projected_site_spin
   use test_response_basis_mod, only: run_test_response_basis
   use test_response_space_mod, only: run_test_response_space
   implicit none

   character(len=128) :: case_name

   call get_command_argument(1, case_name)
   select case (trim(case_name))
   case ('endpoint_branch_action')
      call run_test_endpoint_branch_action()
   case ('gf_endpoint_augmentation')
      call run_test_gf_endpoint_augmentation()
   case ('pauli_transition_vertex')
      call run_test_pauli_transition_vertex('')
   case ('pauli_soc')
      call run_test_pauli_transition_vertex('soc')
   case ('pauli_overlap')
      call run_test_pauli_transition_vertex('overlap')
   case ('pauli_noncollinear')
      call run_test_pauli_transition_vertex('noncollinear')
   case ('pauli_spdf')
      call run_test_pauli_transition_vertex('spdf')
   case ('pauli_additive')
      call run_test_pauli_transition_vertex('additive')
   case ('product_response')
      call run_test_product_response('')
   case ('product_response_strict')
      call run_test_product_response('strict')
   case ('product_response_basis')
      call run_test_product_response_basis()
   case ('projected_site_spin')
      call run_test_projected_site_spin('')
   case ('projected_site_spin_spdf')
      call run_test_projected_site_spin('spdf')
   case ('response_basis')
      call run_test_response_basis()
   case ('response_space')
      call run_test_response_space()
   case default
      error stop 'unknown LR basis unit-test case'
   end select
end program test_basis_main
