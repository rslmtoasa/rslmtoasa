integer :: lld, llsp, nlim, npold, nsp, idos, mext, txc, partype, terminator, random_vec_num, cond_ll
integer :: lmax
real(rp) :: conca, concb, ruban
logical :: lrot, incorb, do_asd, svac, blockrec, do_cochg, asd_jij, do_comom, hyperfine, sym_term
logical :: export_hamiltonian, validate_backend_roundtrip
character(len=sl) :: calctype, recur, cond_type, cond_calctype
character(len=32) :: recur_backend, recur_precision
character(len=sl) :: recur_backend_plugin, recur_backend_library, export_hamiltonian_path
! Constraints are provided in a separate namelist file `constraints.f90`

namelist /control/ recur, lld, llsp, nlim, npold, nsp, &
   idos, lrot, incorb, do_asd, mext, &
   svac, calctype, txc, blockrec, partype, do_cochg, asd_jij, terminator, conca, &
        concb, ruban, do_comom, hyperfine, sym_term, random_vec_num, cond_ll, cond_type, cond_calctype, &
        lmax, recur_backend, recur_backend_plugin, recur_backend_library, recur_precision, export_hamiltonian, export_hamiltonian_path, &
        validate_backend_roundtrip
