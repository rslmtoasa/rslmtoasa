# GBT WP01 q=0 fixture

The production fixture is the committed fixed-potential bcc-Fe deck at
`tests/gbt_wp01/cases/bccFe`.  It uses the Fe potential carried by the earlier
WP5 fixed-potential deck and the `strux_lib` structure-constant backend for
both runs; `val18_gbt_wp01_q0.py` patches only the magnetic representation and
sets `q_ss = (0, 0, 0)`.  The runner compares the ordinary `periodic_nc` and
`gbt_single_q` outputs from the production `kspace_green` stack.

The direct operator gate is `UnitGbtWp01Q0`, which checks directed real-space
blocks, onsite blocks, reciprocal `H(k)`, global SU(2), a two-sublattice frame,
and a negative endpoint-frame control at relative tolerance `1e-12`.
