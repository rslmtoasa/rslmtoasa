# C6.6 follow-up — current-tip gates and Ni(001) response probe

Date: 2026-10-25. Branch: `fix/ldau-surface`; source HEAD before this record:
`7681fc44fad63e3a6849aa4d06711f0a6e7ecb97`.
The user authorized exceeding the prior run budget for this investigation and selected the
surface response criterion `Diff <= 5e-8 within two steps`. All solver launches below were
single-threaded. The two Ni surface launches used `tools/run_rslmto_checked.py` and wrote matching
`preflight.json` records before launch.

## T1 gates, rerun on current tip

- Single-thread RS identity: **PASS, 27/27 artifacts**, compared with the versioned refreshed
  current-runtime oracle `C6_6_RS_ORACLE/artifacts`.
- Single-thread KS identity: **PASS, 24/24 artifacts**.
- Provenance guard fixtures: **PASS, 13/13** (`rs_embedding_ra2.py`).
- Impurity auto and surface auto unit fixture: **PASS, 1/1** (`UnitAutoImpurityLocalMode`).
- Frozen-plus-recursive-site oracle: **PASS, 1/1** (`UnitLdaUProjectorMultisiteFrozenH`).

The identity rerun used two solver processes; RA2 used thirteen. Both CTest cases were run
serially from `/private/tmp/c66-nolibxc-build`. Machine-readable records are archived as
[`C6_6_T1_IDENTITY_RERUN.json`](C6_6_T1_IDENTITY_RERUN.json) and
[`C6_6_T1_GUARD_RERUN.json`](C6_6_T1_GUARD_RERUN.json).

## Geometry and deck audit

Both Ni decks were copied from `tests/scf/cases/surface/fccCu001`. Source inspection and the
actual `surfclu.out` confirm that generated layer coordinates are in `alat` units; no geometry
lines were edited. The active Ni surface recursion site is at projected layer coordinate 0.0 and
the frozen Ni host representative at 0.5. The ES recursion type remains the fixture ES.

The exact line-by-line deck audit with reasons and archived preflight-bound input decks is in [`C6_6_surface_deck_audit/DECK_AUDIT.md`](C6_6_surface_deck_audit/DECK_AUDIT.md). The only `input.nml` edits were: Ni `alat=3.52340` and `wav=1.37692`; Ni host/surface labels;
`nstep=2` and `conv_thr=5e-8` for the approved response probe; `fermi=-0.061141` and
`energy_max=0.5` from the accepted Ni host input; `nsp=1`; and `beta=0.1` from the Ni host deck.
Both decks set `energy_min=-1.0` as in the fixture/host. The appended Hubbard block applies U/J
only to d channels of host type 1 and surface type 3: U=J=0 for the control, U=3/J=0.7 with
`hubbard_radial_moments='auto'` for the test. The ES potential file is unchanged. The surface
seed is a copy of the accepted same-material host potential with only its `symbol` rekeyed to
`Ni-imp-S`.

Accepted host sources:

- U=0: `C5_2_PART_D_assets/ni_u0_bulk_rc80_host/accepted/iter0024_state.nml.gz`
  (SHA256 `f6787f8afb39b3067ea01e9ab4a8657ca030cdf74c6f69db60951d7e51255e81`).
- RS-local U=3/J=0.7: `C5_2_PART_D_assets/ni_rd3_ni_record/bulk_hosts/local_u3/iter0074_state.nml`
  (SHA256 `bc3340a8ad20d2467eaeccd43929a659d4cfff2a9f28b748fb6554edc4ee3ac9`).

For the mandatory host-EF check, each checked run's one-row `host_ef_trace.csv` was derived from
the accepted host input's exact `fermi=-0.061141` with `fix_fermi=.true.`. Preflight reported
PASS with delta 0.0 Ry for both runs. The trace-file SHA256 is
`5cf0db91849e60289f67e8671628b187416eb5a781c8a1978eac0d6e497c51b3`.
Both records report `status=RUN_COMPLETED`, preflight PASS, and return code 0. Preflight JSON SHA256: U=0 `84ccd4e2857672ed4c57b1c269eca03c0967d39cd7ed2314d5d4c9042a0772e6`; U=3/J=0.7 `8cae5a4547cd7a01126cf8f31c0cebefb97e2fef73e51d2f1f8133a2a0f0375e` (each embeds the host-EF check and solver return status).

## T2 Ni(001), RS host

The active site is the Ni surface recursion type (`Ni-imp-S`, type 3); the frozen Ni host
(`Ni-imp`, type 1) provides the bulk comparison. Occupancies below are the final two-step output
`ql(1,l=2,spin)` in electrons. Spin moments are from the surface `report.out` and rounded there
to 1e-6 μB; the two spin-resolved d occupations are equal to the displayed precision.

| Run | Mode | Layer/type | d up | d down | d total | reported spin moment (μB) |
|---|---|---|---:|---:|---:|---:|
| U=0 | spectral (no U) | frozen Ni host | 4.701519526337 | 4.701519599652 | 9.403039125989 | — |
| U=0 | spectral (no U) | active Ni surface | 4.231290265016 | 4.231290328218 | 8.462580593233 | 0.000000 |
| U=3/J=0.7 | local | frozen Ni host | 4.738761219470 | 4.738761362201 | 9.477522581671 | — |
| U=3/J=0.7 | local | active Ni surface | 4.265437775585 | 4.265437884558 | 8.530875660143 | 0.000001 |

Surface minus bulk d occupation is `-0.940458532756 e` at U=0 and `-0.946646921527 e` at
U=3/J=0.7. The U-dependent change in that difference is `-0.006188388771 e` (per spin:
`-0.003094182563` and `-0.003094206208 e`). These are finite two-step response measurements,
not converged surface observables.

| Run | Diff step 1 | Diff step 2 | Gate |
|---|---:|---:|---|
| U=0 control | 0.0162421175 | 0.3840255893 | control only; not a gate pass |
| U=3/J=0.7 auto | 0.0163260002 | 0.3869175281 | **FAIL** (`mode=local`, no provenance error, but Diff exceeds 5e-8) |

The archived evidence bundles, including each matching `preflight.json`, are [`T2a_U0`](C6_6_surface_deck_audit/T2a_U0/) and [`T2b_U3J07_auto`](C6_6_surface_deck_audit/T2b_U3J07_auto/). Original run directories were `/private/tmp/c671-ni-surface-rs/T2a_U0` and `/private/tmp/c671-ni-surface-rs/T2b_U3J07_auto`. Run-log SHA256 values are respectively
`eba61571c7cbc5110ebd0f94797a015ec2fbc8c24b4089cd704a2edf79450120` and
`f6924085dfad976868e00139cf863b7aff619cf4902528dca339c726ac9db6bc`. The auto run selected
`mode=local source=default`; it emitted no provenance error. Surface Hubbard-energy routing was
not exercised and remains OFF.

## Disposition

The T1 gates pass, and the C6 local-moment correction reaches the U-carrying Ni surface site in
local mode. The approved two-step Diff gate fails decisively, so execution stopped here: Cu T3,
workflow documentation, merge, and adoption were not performed. The user-provided six-layer
no-U surface snapshot remains a candidate rather than an accepted RS reference: its manifest
records only one Diff value (6.68e-8), no matching preflight, and a host/deck E_F mismatch. The
C5.3 item 2 reference requirement therefore remains uncovered.
