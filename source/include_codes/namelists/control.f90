integer :: lld, llsp, nlim, npold, nsp, idos, mext, txc, partype, terminator, random_vec_num, cond_ll
integer :: lmax
real(rp) :: conca, concb, ruban, dipole_mix
logical :: lrot, incorb, do_asd, svac, blockrec, do_cochg, asd_jij, do_comom, hyperfine, sym_term
logical :: cpp_plugin
logical :: gpu_plugin
logical :: dipole_electrostatics
character(len=sl) :: calctype, recur, cond_type, cond_calctype, gpu_backend, cheb_backend
character(len=sl) :: linear_in, linear_out
! Constraints are provided in a separate namelist file `constraints.f90`

namelist /control/ recur, lld, llsp, nlim, npold, nsp, &
   idos, lrot, incorb, do_asd, mext, &
   svac, calctype, txc, blockrec, partype, do_cochg, asd_jij, terminator, conca, &
   concb, ruban, do_comom, hyperfine, sym_term, cpp_plugin, gpu_plugin, gpu_backend, cheb_backend, &
   random_vec_num, cond_ll, cond_type, cond_calctype, linear_in, linear_out, lmax, &
   dipole_electrostatics, dipole_mix
