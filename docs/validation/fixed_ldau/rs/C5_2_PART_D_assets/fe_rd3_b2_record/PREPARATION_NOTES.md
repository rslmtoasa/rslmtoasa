# Fe RD3 endpoint capture preparation (private; not launched)

The accepted Fe U3.8/J1.1 bulk endpoint P_OUT trace and accepted local impurity endpoint P_OUT trace both preserve six q0 occupations (s/p/d by spin), but neither contains energy-resolved Green matrices or a d-LDOS curve. Existing `Fe_dos.out`/`Fe_orbital_dos.out` files are later-live outputs and are excluded from the endpoint comparison.

Bulk endpoint: iteration 92, state SHA `6e93e4141e8a84b1d4c3bac0d1f62997bf3c702e56e923f5c043e7becb86d16b`; trace SHA `a8297b63b4a84392e81dd17e6cc025204fad07481107b78954e085f5420fd69d`. P_OUT d occupations are 4.6519194389087062 / 1.8573696818368057; summed-l q0 is 5.30701807491765389 / 2.69298192508234499, giving `mtot=|sum(up)-sum(down)|=2.61403614983530890 μB`.

Impurity endpoint: iteration 20, state SHA `cb4b3b2f2eca12d48030e6e2447a619bbe53ed8f2ba61988818fde64dda165aa`; trace SHA `1646f4dcfd4a` is the suffix of the named full trace (full hash is in its preflight/package sums). P_OUT d occupations are 4.6519207879528413 / 1.8573618222099855; summed-l q0 is 5.30702431778329014 / 2.69297837348953865, giving `mtot=2.61404594429375149 μB`.

No matched Fe U=0 control is accepted, so no RD3 ratio is computed. Bulk calctype B has free EF; impurity calctype I defaults to fixed host EF. Both serialized values are −0.085152 Ry at six decimal precision; exact higher precision identity is not established.

The prepared capture pair uses the same frozen RD1 Release executable SHA `4bb76cc2bbd251d95cdc434eec0348b9082069af2deeb4f4859f6d8056779a4f`, same bcc/RS/block/rc80/lld21/window[−1,0.5], local mode and U/J values. Each deck changes only `self.nstep` to 1 and uses the exact accepted recursive state as restart. Captures remain PREPARED ONLY. The d-LDOS parser follows `bands.f90:492-496` (`−Im(Gii)/π`) and requires complete 18×18 basis matrices plus an identical energy grid; it does no interpolation and introduces no gate. A self-comparison smoke on the existing C3 test capture validates parser shape/zero-difference behavior only, not Fe physics.
