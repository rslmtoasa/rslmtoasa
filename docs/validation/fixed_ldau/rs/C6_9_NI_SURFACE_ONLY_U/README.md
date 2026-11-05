# Ni(001): U only on the surface Ni site

This is a two-step density-response probe of the surface-impurity model. The frozen bulk Ni host and all starting potentials were copied from the user's no-U surface run. The `Ni-S` recursion type (type 3, surface site) alone receives U=4/J=1 on its d channel; `Ni`, `ES`, and `Ni-S-1` have U=J=0. The surface geometry, `alat`, `wav`, energy window [-1.0, 1.5] Ry, fixed Fermi level -0.062420 Ry, and mixing parameters are unchanged from the source no-U deck. The probe sets `nstep=2`, `conv_thr=5e-8`, and `hubbard_radial_moments='local'`.

The U=0 control is the user's original one-step run, preserved in `seed_noU_control/`. It reports Diff=7.71e-8 and spin moments (ES / surface Ni / subsurface Ni) 0.007355/0.689005/0.640933 μB. Its fixed Fermi energy is -0.062420 Ry.

## U-only run

Launched with `.venv/bin/python tools/run_rslmto_checked.py ... --binary build/bin/rslmto.x --host-ef-trace ../host_ef_trace.csv`. The wrapper's `preflight.json` records PASS for the host-EF check (delta 0) and RUN_COMPLETED (return code 0). The solver reports `HUBBARD_RADIAL_MODE mode=local`; no provenance error or NaN/fatal error occurred. Its fixed Fermi energy remains -0.062420 Ry.

| SCF step | Diff | spin moments atoms 1 / 2 / 3 (μB) |
|---|---:|---:|
| 1 | 0.1149764639 | 0.037554 / 1.081094 / 0.489236 |
| 2 | 0.0389620957 | 0.010834 / 0.112752 / 0.129423 |

At step 2, the local d-shell occupation printed for `Ni-S` is n(up)=2.644274, n(down)=1.573365, total=4.217638. Types 1, 2, and 4 have zero Hubbard potential; type 3 has max |V_U|=0.0468 Ry. The moments and occupation are still changing substantially, and the requested Diff <= 5e-8 was not reached within two steps. Treat this as an initial response probe, not a converged surface result.

## Artifacts

- `surface_u4j1_only/input.nml`: deck used for the U run.
- `surface_u4j1_only/preflight.json`: preflight and launch record.
- `surface_u4j1_only/testrun.log`: solver output.
- `seed_noU_control/`: original no-U input, output, and potential files used to seed the run.
- `host_ef_trace.csv`: fixed host Fermi energy extracted from the original no-U output.

The output recursion-site mapping is atom 1 = ES, atom 2 = surface Ni (`Ni-S`), and atom 3 = subsurface Ni (`Ni-S-1`); the frozen bulk host is type 1 and is not one of these reported recursion atoms. See `C6_10_RESULTS.md` for the matched two-step U=0 control, per-step comparison, occupation results, and current-binary gates.
