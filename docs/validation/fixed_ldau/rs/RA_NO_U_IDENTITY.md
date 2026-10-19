# RA no-U identity result

**Measured: pass.** No-U bcc-Fe RS and KS output snapshots from the RA worktree
are byte-identical to runs from a clean Debug/OpenMP executable built from the
pinned `8b8f7ef7bfc026fd34a56b37a3ad3281fd0a6832` source. The RS fixture ran
one iteration and the KS fixture ran its configured two iterations. The comparison covers
`Fe_out.nml`, `Fe_scf.nml`, `Fe_dos.out`, `Fe_orbital_dos.out`, `totaldos.out`,
and `magneticdos.out` where present. This checks the protocol’s no-U output
boundary for these tested input decks.

The final baseline executable SHA256 is
`8307387ecead3893f63b5348cb84a75ce1bd886becdec088e0bfa7b5e895e9fd`; the
final RA executable SHA256 is
`6ec0e5820c711115570577ddf3b2ffe9fe8c9e222e6a50061fadf773cffc0041`.
Input and output file SHA256 values and the raw output path are in
[`RA_NO_U_IDENTITY.json`](RA_NO_U_IDENTITY.json); raw files remain under
`/tmp/rslmto-ra-no-u-latest/`.

The exact source archive, configure, build, and paired-run command is recorded
in the JSON provenance. Both builds used GNU Fortran 16.2, Debug, OpenMP on,
MPI and libXC off. The measured scope is these RS/KS fixtures, not every
possible no-U input.
