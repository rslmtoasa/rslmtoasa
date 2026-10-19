# RA1 exploratory mapping amendment

Recorded 2026-10-18 19:42 UTC. The original RA1 preregistration states that
`lattice%ib(ia)` supplies the energy-evaluator type. An exploratory B2FeCo
run falsified that assumption: for its one recursive atom, the log showed
`lattice%ib(1)=16`, while the three symbolic atom types are the two bulk types
and the impurity type at `nbulk+ia=3`. The `ib` value is the geometric
representative index; the Hubbard energy evaluator arrays are keyed by
symbolic-atom type.

The first independent comparison therefore failed: the diagnostic emitted
`E_int=-0.9080520743969145 Ry`, `E_DC=0`, and `E_Hub=0`; evaluation of the
serialized `Fe_out.nml` d-occupation matrix with the standalone production
energy evaluator gave `E_int=2.6637383372890966 Ry`,
`E_DC=2.5582858248167537 Ry`, and `E_Hub=0.10545251247234289 Ry`. This was
treated as a mapping defect, not a tolerance failure. The diagnostic now
indexes its energy arrays with `nbulk+ia`; the separate radial snapshot path
continues to use `ib(ia)` for its lattice-type arrays.

No acceptance threshold changed: the preregistered maximum difference remains
`1e-12 Ry` per energy component. The final scored rerun reported symbolic type
3 and passed with maximum absolute difference `8.891e-13 Ry`. Raw exploratory
and scored run products are in `/tmp/rslmto-ra1-final2/` and
`/tmp/rslmto-ra1-verified/`, respectively.
