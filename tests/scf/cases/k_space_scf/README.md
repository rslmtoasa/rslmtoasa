# Canonical lean k-space SCF fixtures

These cases are functional baselines for the reciprocal-space SCF branch, not
converged material-physics benchmarks.

`Example_k_space_scf_diamondSi_sp` reuses the canonical two-site diamond-Si
sp fixture from `bulk/diamondSi`: `lmax=1`, `nsp=1`, one SCF step, a full
2x2x2 Gamma-centered mesh, Gaussian DOS broadening `sigma=0.01 Ry`, and a
200-point DOS grid. The equal Si spin channels make this a physically
nonmagnetic case, while `nsp=1` remains the scalar-relativistic collinear
spin-polarized route.

`Example_k_space_scf_bccFe_spd_magnetic` uses the bcc-Fe spd fixture with
`lmax=2`, `nsp=2`, one SCF step, the same 2x2x2 mesh, Gaussian broadening
`sigma=0.01 Ry`, and a 200-point DOS grid. The canonical checks include the
nonzero site spin moment in addition to occupations, energies, site charge,
and selected DOS values.

Both cases select `strux_backend='strux_lib'` explicitly and register serial
and 2-rank MPI launches against one shared reference. Reference values were
generated from clean production runs with the recompiled `build-gc-serial` /
`build-gc-mpi` binaries on 2026-08-13. Wall-clock timings were 1.43 s serial
and 1.96 s on 2 ranks for Si, and 4.02 s serial and 4.62 s on 2 ranks for Fe.
These timings are provenance for the lean fixture, not performance targets.
