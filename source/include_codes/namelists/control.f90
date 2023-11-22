integer :: lld, llsp, nlim, npold, nsp, idos, mext, txc, partype, terminator
real(rp) :: conca, concb, ruban
logical :: lrot, incorb, do_asd, svac, blockrec, do_cochg, asd_jij, do_comom, hyperfine, sym_term
character(len=sl) :: calctype, recur

namelist /control/ recur, lld, llsp, nlim, npold, nsp,  &
idos, lrot, incorb, do_asd, mext , &
svac, calctype, txc, blockrec, partype, do_cochg, asd_jij, terminator, conca, &
concb, ruban, do_comom, hyperfine, sym_term
