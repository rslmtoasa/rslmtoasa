integer :: lld, llsp, nlim, npold, nsp, idos, mext, txc, partype, terminator, random_vec_num, cond_ll, gpu_stochastic_block
integer :: lmax
real(rp) :: conca, concb, ruban, dipole_mix
logical :: lrot, incorb, do_asd, svac, blockrec, do_cochg, asd_jij, do_comom, hyperfine, sym_term
logical :: cpp_plugin
logical :: gpu_plugin
logical :: dipole_electrostatics
character(len=sl) :: calctype, recur, cond_type, cond_calctype, gpu_backend, gpu_precision, cheb_backend
logical :: gpu_moment_download
character(len=sl) :: linear_in, linear_out
character(len=sl) :: density_policy
! Constraints are provided in a separate namelist file `constraints.f90`

namelist /control/ recur, lld, llsp, nlim, npold, nsp, &
   idos, lrot, incorb, do_asd, mext, &
   svac, calctype, txc, blockrec, partype, do_cochg, asd_jij, terminator, conca, &
   concb, ruban, do_comom, hyperfine, sym_term, cpp_plugin, gpu_plugin, gpu_backend, gpu_precision, gpu_moment_download, cheb_backend, &
   gpu_stochastic_block, &
   random_vec_num, cond_ll, cond_type, cond_calctype, linear_in, linear_out, lmax, &
   dipole_electrostatics, dipole_mix, density_policy
