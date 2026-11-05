# C6.4 — Stage RE smoke validation disposition

Branch `fix/ldau-surface`, starting tip `1ebee1ef7a3b320d50ded696ca6d4b4815191a0c`.
The user-approved C6.4 amendment is recorded verbatim in `C6_AMENDMENT.md`.

## T1 gates

All mandatory T1 gates passed on the working tree containing the C6.4 change. Python was
`.venv/bin/python3` (f90nml 1.4.4); the Release binary was
`/private/tmp/c511-release-build/bin/rslmto.x` (SHA256
`46524061a637a1c9ec7f25df74d26d1d0a74cec65222a2241a2842c990366082`).

- Single-thread no-U identity: RS 27/27 and KS 24/24; zero mismatches. Record:
  `/private/tmp/c64_t1_identity/C5_7_NO_U_IDENTITY.json` (SHA256
  `3fa90c3cf784a92075867a6aedcaf9b39f2ceaa6dd26c77981c00dc739d1a36a`).
- Guard fixtures: 13/13; record `/private/tmp/c64_t1_guards/RA2_provenance.json` (SHA256
  `8808e5262d32c05d238dce88ae850efe3e9c649430b6e9295f612c824a8f541d`).
- Impurity auto and surface auto fixture: `UnitAutoImpurityLocalMode` passed 1/1.
- Frozen-plus-recursive-site oracle, including the surface-like frozen-U / ES / recursive-U
  layout: `UnitLdaUProjectorMultisiteFrozenH` passed 1/1; direct-construction moment difference
  was within `1e-12` and the U-free ES site was skipped.

The identity and guard scripts ran single-threaded with C5.4 preflight before each solver
invocation. Together they launched 15 validation solver processes, outside the four-run surface
campaign budget. The focused CTest targets launched none.

## T2 Ni(001), RS host

The copied surface geometry creates layer positions by projecting dimensionless lattice
coordinates (`cr`) without modifying the geometry; Cartesian output is scaled by `alat`. Thus
the fixture positions are in alat units. The staged decks were audited against the authorized
changes: Ni `alat`/`wav`, Ni labels and accepted host/surface seeds, Ni-host `fermi` and window,
Ni d-channel U/J, `nsp=1`, `nstep=200`, and the accepted Ni bulk `&mix` parameters. ES remained
the fixture empty-sphere potential.

### T2(a), U=0 control

The run was launched from `/private/tmp/c64_t2a_ni_u0_run` without the required C5.4 preflight.
This was an execution-procedure error; it is counted as one of the four authorized SCF launches,
but is not valid preflighted baseline evidence. It stopped at iteration 8 with the B1 non-finite
guard at `atom=ES`, `field=PL`, `value=NaN` (log SHA256
`0fc6c9b445e6da22ee0b73ec500f5dde99f2f246042794019755266ff87a29f5`). No baseline layer
occupations or moments are reported.

### T2(b), U=3/J=0.7 auto

C5.4 preflight passed in `/private/tmp/c64_t2b_ni_u3j07_run`; deck and accepted RS-local Ni host
endpoint E_F both equal `-0.061141` Ry (delta 0). The runtime selected
`HUBBARD_RADIAL_MODE mode=local`, and the log contains no provenance error. The run did not
converge: its outer Diff remained `0.0211502370` at the last completed iteration before the
moment routine stopped with a negative rotating-frame density eigenvalue at site 2, l=0
(`lambda_min = -3.15867E-06`). It therefore fails the C5.13 smoke gate (final Diff `<=1e-7` and
20-step spread `<=1e-5`); no layer-resolved converged occupations/moments or U-dependent
surface-minus-bulk comparison can be reported. Log SHA256:
`c055e03991b46d7011cc33ffa0b0e5630e077ae3a393f8b4caca925f22253a12`.

T2 consumed two campaign launches out of four. Since the T2(b) smoke gate failed, T3(a/b) and
the T4 merge checks were not run. No T3 runs were launched. The surface workflow documentation
and merge are deferred until the local-moment physicality failure is resolved.

## Disposition

The C6.4 local site/type mapping code change and T1 gates pass. Stage RE smoke validation is
incomplete because the U=0 Ni baseline failed its non-finite guard and the U>0 Ni smoke failed
the C5.13 convergence/physicality gate. No converged surface occupations or moments are claimed.
The two attempted SCF launches consumed 2/4 of the authorized campaign budget. This is an OPEN
validation issue; see `tests/KNOWN_ISSUES.md`.
