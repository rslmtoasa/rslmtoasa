# MnTe tutorial starting state

These six input/checkpoint files were supplied by the instructor in
`mntenotbroken.tar`. The archive SHA-256 is
`c575273c8039b592dc70578b73362a15f39c3e49c787a03267f7e2d66da77da0`.
The crystallographic attribution in `input.nml` is retained. Trailing spaces
and the final extra blank line were removed; parsed namelist values are unchanged.

The four electronic types are Mn1, Mn2, Te1, Te2 (atomic numbers 25, 25,
52, 52), with `izp = 1,2,3,4`. Mn1 and Mn2 have opposite starting directions
along x. `nsp=4` includes noncollinearity and SOC. Both Te types remain in
the electronic Hamiltonian; their `mom` entries are orientation vectors,
not spin-moment magnitudes.

The supplied settings are Chebyshev `lld=200`, PBC `15x15x15`,
`alat=4.1501998901` Angstrom, and `hoh=.false.`. No report or convergence log
was included, so the files alone do not establish the original SCF residual.
The duplicate `_scf.nml` files in the upload were identical to `_out.nml`
and are omitted here. These `_out.nml` files are starting inputs for the
lesson, not newly generated tutorial results.

Copy this entire directory to `runs/alternative-mnte-scf`; never run in this
tracked directory. Follow the MnTe section of
[the alternative tutorial](../../docs/tutorial-alternative-workflows.md).
Its staging commands bound the SCF replay to one step, clear the supplied
unused `njij=162`, and explicitly disable Morton reordering. Geometry,
energy settings, spin treatment, and checkpoints are retained. The teaching
copy reduces PBC to 11x11x7 (3,388 sites); the original files here retain
15x15x15 (13,500 sites). This changes sampling and requires a new convergence
check; do not assume the supplied state stays self-consistent after resizing.
