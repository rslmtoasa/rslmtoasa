program test_lmto_path_operator_legacy_alpha_guard
   use lmto_path_operator_mod, only: native_screening_alpha
   use precision_mod, only: rp
   use symbolic_atom_mod, only: symbolic_atom
   implicit none

   type(symbolic_atom) :: atom
   real(rp) :: alpha(0:4)

   atom%potential%lmax = 4
   call native_screening_alpha(atom, alpha)
   error stop 'native_screening_alpha accepted legacy screening constants beyond l=3'
end program test_lmto_path_operator_legacy_alpha_guard
