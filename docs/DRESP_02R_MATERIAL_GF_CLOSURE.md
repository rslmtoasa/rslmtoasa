# DRESP-02R — Material GF/Lehmann Closure Remediation

Status: `BLOCKED — MATERIAL GF CLOSURE UNRESOLVED`; DRESP-02 algebraic
certification remains valid and DRESP-02P has removed the performance blocker.

The exact controlled material rerun and performance evidence are recorded in
[`DRESP_02P_PROJECTED_GF_PERFORMANCE.md`](DRESP_02P_PROJECTED_GF_PERFORMANCE.md).

This remediation stays inside DRESP-02. It does not open a new response rung,
change the DRESP-01 operator, change the product oracle, change a tolerance,
or enter KXC, Dyson, Goldstone, loss, or mode fitting.

## Locked facts

The following results are established and are not re-litigated here:

| check | result |
| --- | ---: |
| direct projected Lehmann versus projected product Lehmann | `3.18e-22` max residual |
| direct projected GF versus projected product GF | `9.91e-18` max residual |
| direct projected site response allocation | none |
| required two-site matrix and q/−q covariance unit tests | PASS |

The independent product GF remains an oracle only. The remediation adds
observability to the direct site GF path without routing the production GF
through that oracle.

## Material-state identity gate

The accepted material state is one frozen bcc-Fe `4x4x4` full Monkhorst-Pack
state: 64 k points, `ham_only`, second-order Hamiltonian, orthogonal,
collinear, no SOC, no extra operator, 300 K, `EF=-0.0612124078445383 Ry`,
and `response_lmax=4`. The GF ladders reuse the same left eigenpairs,
occupations, radial/LMTO snapshots, k weights, and exact folded endpoint;
the driver does not perform SCF, rebuild the Hamiltonian, or regenerate an
endpoint inside a ladder.

The output now writes deterministic state fingerprints and eigenvector
unitarity residuals. In the accepted Γ audit artifact they were:

```text
eigenvalue checksum       4.2910724817163382e3 Ry
eigenvector checksum      8.5882320444735899e2 - 1.2856885642972964e2 i
occupation checksum       2.3482472833262864e3
unitarity max residual    2.2204462702395403e-15
radial checksum           1.4274436288316949e4
radial L2 norm             2.6450730612151682e2
```

The Γ endpoint fingerprint equals the left-state fingerprint, as required for
the exact folded endpoint. Every request also validates the folded k+q
coordinates, dimensions, weights, EF, temperature, energy zero, occupation
semantics, reciprocal mode, Hamiltonian order, spin flags, and DRESP-01
projection dimensions before contraction. The state artifact records the
full eigenvalue/occupation and occupation-weighted projector data.

## Moment reconciliation

The earlier DRESP-01 material gate (`8x8x8` accepted state) reported
`d=2.096590` and `spd=2.003489 mu_B`, with accepted total
`2.003488859 mu_B`. DRESP-02 uses the separately accepted `4x4x4` response
state and reports `d=1.970011188` and `spd=1.950787628 mu_B`, with accepted
total `1.950787628 mu_B`. The mismatch is therefore a state/mesh provenance
difference, not evidence that the DRESP-01 selector or DRESP-02 trace was
silently rewritten. The two states have different accepted SCF moment
artifacts and must not be mixed in a closure comparison.

On the DRESP-02 state, the accepted band-moment `spd` sum reproduces the
accepted total within `3.3e-11 mu_B`. The direct DRESP operator moments are
reported separately (`d=2.020544770`, `spd=2.001437103`); they are a
projection/operator diagnostic and are not substituted for the accepted SCF
band-moment gate. Core channels remain excluded by policy.

## DRESP-02R observables

`projected_chi0_request%diagnostics=.true.` now records, without changing the
response equation:

- `eta_response`, `eta_int`, `Emin`, `Emax`, `NE`, `h`, and `h/eta_int`;
- both Kubo terms separately, with an internal term-sum residual;
- per-k projected site matrices and per-k closure residuals;
- finite-window zeroth, first, and Fermi-weighted spectral-moment residuals
  for both left and right endpoint states;
- the largest eight direct Lehmann transitions, including k/band indices,
  occupations, transition energy, matrix-element weight, and score.

The material driver writes three separate campaigns, all against the same
Gamma Lehmann result:

1. fixed `eta_int` with `NE=1001,2001,4001` to isolate h/mesh error;
2. an `eta_int` ladder `4,2,1` times the base value, with points selected to
   hold `h/eta_int` near `0.40`;
3. an energy-window ladder with fixed base mesh and `gf_energy_margin`
   values `0.30,0.60,1.00 Ry`.

The exact machine-readable rows are emitted as `gf_audit`,
`gf_spectral_moments`, `gf_kubo_terms`, `gf_k`, and `dominant_lehmann` lines
in the projected response artifact.

## Material evidence obtained

The frozen-state 4x4x4 Γ fixed-`eta_int=0.002 Ry` mesh campaign completed
the d projection through 4001 points. Its required fields were:

| selector | NE | h/eta_int | `||chi_L||` | `||chi_GF||` | abs difference | relative |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| d | 1001 | 1.4890 | 26.1428 | 25.9945 | 0.2990 | 0.01144 |
| d | 2001 | 0.7445 | 26.1428 | 25.9307 | 0.3773 | 0.01443 |
| d | 4001 | 0.3723 | 26.1428 | 26.0338 | 0.3180 | 0.01216 |

At the finest d sample, the independent Kubo terms summed internally to a
`2.1e-12` Frobenius residual. The finite-window spectral residuals were:

| moment | left | right |
| --- | ---: | ---: |
| zeroth | `1.47e-3` | `1.47e-3` |
| first | `2.36e-3` | `2.36e-3` |
| Fermi-weighted | `7.90e-2` | `7.90e-2` |

The finite-temperature audit used the immutable 300 K Fermi occupations in
both the Lehmann sum and the GF factor `f(E)`; no zero-temperature step
function was substituted. The already completed same-state DRESP-02
finite-frequency spot at Γ confirms that the unresolved scale is present at
`omega=0.015 Ry` as well:

| selector | omega (Ry) | chi11 Lehmann (Re, Im) | chi11 GF (Re, Im) | abs difference |
| --- | ---: | --- | --- | ---: |
| d | 0.015 | `(-28.6335,-1.88439)` | `(-74.9168,-6.81628)` | `46.5454` |
| spd | 0.015 | `(-35.4206,-8.77265)` | `(-161.622,-69.3080)` | `139.969` |

The existing completed material runs at the same accepted 4x4x4 state also
give the 4001-point `spd` fixed-window result `||chi_L-chi_GF||=1.1648`,
with `||chi_L||=30.785` and `||chi_GF||=29.978`. Earlier completed material
mesh samples were:

| NE | d difference | spd difference |
| ---: | ---: | ---: |
| 1001 | `2.9904e-1` | `2.1281` |
| 2001 | `3.7732e-1` | `1.3917` |
| 4001 | `3.1802e-1` | `1.1648` |

These values are far outside the unit-oracle residuals and are not a closure
pass. The fixed-eta d sequence is non-monotone; the spd sequence decreases
but remains macroscopic. The newly added controlled-eta and window campaigns
are implemented in the driver. The original material run stopped after the
4001-point d diagnostics because of the performance budget; DRESP-02P
subsequently completed the full two-selector campaign with the same frozen
state and controls. Its residual tables and timing profile are reported in
`DRESP_02P_PROJECTED_GF_PERFORMANCE.md`.

## Root-cause classification and exact return action

The evidence separates the layers as follows:

- state identity and endpoint construction pass machine checks;
- DRESP-01 projection algebra and both compact projection oracles pass;
- the GF Kubo terms close internally, so the reported mismatch is not a
  missing second term in the diagnostic accumulator;
- finite-window spectral moments, especially the Fermi-weighted moment, are
  not yet resolved at the material sample;
- the optimized material GF campaign is now executable at the controlled
  resolution; its remaining GF/Lehmann residual is therefore no longer
  classified as a performance failure.

Accordingly the performance blocker is cleared, but this gate is not `PASS`.
Use the completed h/eta, eta, window, per-k, term, moment, and transition
tables for the next DRESP-02 closure decision. The tolerance remains
unchanged. If the controlled ladders converge to the Lehmann matrix, return
`PASS — numerical closure`; if they converge to a nonzero offset after all
controls are resolved, return `BLOCKED — FORMULATION` with the offending
term/transition identified.

After either a genuine closure pass or a separately documented formulation
resolution, return directly to DRESP-02’s next authorized integration step.
Do not start a new DRESP rung from this remediation.

## Verification

The implementation was verified with:

```text
cmake --build build --target UnitLrProjectedReciprocalChi0 rslmto.x -j2
ctest --test-dir build -R 'UnitLrProjectedReciprocalChi0|UnitTddftProductionDriver|UnitLrProduct(KsSusceptibility|GfSusceptibility)$' --output-on-failure
```

The focused regression passed `8/8` tests. The material scratch campaign was
run from `/tmp/dresp02_fe_4k` with `backend=projected_chi0`, Γ, `omega=0`,
`eta=0.01 Ry`, base `eta_int=0.002 Ry`, and `gf_integration_points=4001`.
