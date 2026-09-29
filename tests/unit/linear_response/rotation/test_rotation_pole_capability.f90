module test_rotation_pole_capability_mod
   use linear_response_mod, only: lmto_live_hamiltonian_fixture, lmto_fixture_init, &
      validate_rotation_pole_workflow_capability
   implicit none
   private
   public :: run_test_rotation_pole_capability, run_test_rotation_pole_multisite

contains

   subroutine run_test_rotation_pole_capability()
      type(lmto_live_hamiltonian_fixture) :: fixture

      call lmto_fixture_init(fixture,1,1,1,.true.)
      call validate_rotation_pole_workflow_capability(fixture)
      call fixture%clear()
      write(*,'(a)') 'one-site rotation pole capability: PASS'
   end subroutine run_test_rotation_pole_capability

   subroutine run_test_rotation_pole_multisite()
      type(lmto_live_hamiltonian_fixture) :: fixture

      call lmto_fixture_init(fixture,2,1,1,.true.)
      call validate_rotation_pole_workflow_capability(fixture)
      error stop 'multisite rotation pole capability unexpectedly accepted'
   end subroutine run_test_rotation_pole_multisite

end module test_rotation_pole_capability_mod
