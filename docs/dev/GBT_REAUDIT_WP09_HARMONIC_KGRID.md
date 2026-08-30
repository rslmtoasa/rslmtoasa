# GBT WP09 re-audit: harmonic cone and k-grid convergence

Date: 2026-08-27
Status: instrumentation and evidence workflow delivered; converged physical
WP09 result **not established**.

This report follows the acceptance prompt in
`docs/dev/gbt_luna_blueprint/LUNA_GBT_WP09_HARMONIC_CONE_AND_KGRID.md`.
The production campaign was run against baseline revision
`3a22de4a21f1dbb052a2f66f1d968964395396b0` with the WP09 source changes in
the working tree. The WP09 implementation commit is
`49a13abde91d9b8230cf78d602b4c0bd69e2e086`.

## Scope and bookkeeping

The production reciprocal frozen-magnon path now emits
`frozen_magnon_harmonic_diagnostics.dat` with schema
`gbt_wp09_harmonic_v1`. Each row retains:

- `E(q,theta)`, `E(qref,theta)`, raw `DeltaE`;
- bare-MFT same-q `E(q,0)`, `E(qref,0)`, gauge-subtracted `DeltaE`, and pure
  drift `E(q,0)-E(qref,0)`;
- `DeltaE/sin(theta)^2`, `Mtot`, and `omega`;
- Fermi level, canonical electron count, electron error, target count,
  occupation-weight sum, and the full k mesh;
- the q/2 mesh-mapping control.

For corrected constrained MFT, the schema marks the same-q gauge fields as
unavailable. Its raw occupied-eigenvalue difference is retained, while
constraint penalty and field-coupling terms remain diagnostics only.

`tools/analyze_gbt_wp09.py` reads this schema and the legacy WP08 files. It
reports raw/gauge/pure definitions separately, selects a contiguous harmonic
plateau, and labels k-grid comparisons `converged` only when the full mesh
spread is at most either `1e-6 Ry` absolute or `1%` relative. The default
plateau tolerance is `5%`, with at least two angle points.

## Reproducibility record

Compiler and tools:

- GNU Fortran 13.3.0 (`/usr/bin/gfortran`)
- CMake 3.28.3
- Python 3.12.3
- Fortran flags: `-ffree-line-length-0 -cpp -march=native`; Debug flags also
  include `-g -O0 -Wall -fbacktrace -fcheck=all,no-recursion`.

Build:

```text
cmake -S . -B build
cmake --build build -j2
```

Fixture and campaign:

- material: committed bcc Fe frozen-potential fixture;
- `alat = 2.86120`, q coordinates Cartesian in `2*pi/alat` units;
- harmonic q list: Gamma and `(0,0,0.05)`;
- finite-q control: Gamma, `(0,0,0.05)`, and `(0,0,0.5)`;
- cone angles: `2, 5, 10, 15, 20` degrees;
- meshes: `8^3`, `12^3`, `16^3`;
- reciprocal policy: `full_bz`; nonzero-q runs use the full mesh without
  symmetry or time-reversal reduction;
- q/2 mapping: uniform Cartesian fallback, with file-level metadata preferred
  for non-cubic lattices.

Bare-MFT production run:

```text
python3 tests/validation/val19_gbt_harmonic_kgrid.py \
  --binary "$PWD/build/bin/rslmto.x" \
  --scratch-root /tmp/gbt-wp09-full \
  --modes mft --require-plateau --timeout 600
```

This produced 39 rows and passed the requested angle-plateau gate. The
normalized sweep and analysis were retained at:

```text
/tmp/gbt-wp09-full/wp09_sweep.csv
/tmp/gbt-wp09-full/wp09_analysis_thresholded.json
```

Corrected constrained-MFT reference run:

```text
python3 tests/validation/val19_gbt_harmonic_kgrid.py \
  --binary "$PWD/build/bin/rslmto.x" \
  --scratch-root /tmp/gbt-wp09-cmft-reference \
  --modes mft_constrained --timeout 600
```

This uses the stable WP08 reference convention: `theta=90` degrees and
q = Gamma, `(0,0,+0.025)`, `(0,0,-0.025)`. It produced three converged
constraint rows. A direct small-angle corrected attempt at `theta=5` degrees
and q = Gamma, `(0,0,0.05)` was also attempted; it failed at q index 2 after
20 constraint iterations, with RMS angle about `1.56 rad` and the spin moment
collapsing to zero. That is classified as a solver-convergence boundary, not
as a valid corrected-MFT harmonic result.

## Results

### Harmonic angle response

At each fixed reciprocal mesh, the bare-MFT gauge-subtracted definition has a
contiguous admitted angle set of `2–20` degrees for q = `(0,0,0.05)` under the
5% plateau rule. Its extracted K values are:

| mesh | mean K (Ry) | relative angle spread |
|---|---:|---:|
| `8^3` | `8.071953e-4` | `2.3946%` |
| `12^3` | `2.257883e-4` | `1.5996%` |
| `16^3` | `4.349479e-4` | `2.3567%` |

The raw q-reference subtraction has no admitted q = `(0,0,0.05)` plateau.
This is the expected diagnostic signature of a q-only finite-mesh artifact,
not evidence that the raw and gauge definitions are interchangeable.

### k-grid convergence

For q = `(0,0,0.05)`, the mesh spread over `8^3`, `12^3`, and `16^3` is:

- raw `DeltaE`: `1.0204e-3` to `1.0889e-3 Ry`, **not converged**;
- pure-gauge drift: `1.0896e-3 Ry`, **not converged**;
- gauge-subtracted `DeltaE`: `7.0162e-7` to `6.9224e-5 Ry` across the angle
  set. The 2° comparison is within the `1e-6 Ry` absolute threshold, but the
  5–20° comparisons are **not converged**.

The finite-q pure drift itself is:

| mesh | `E(q,0)-E(Gamma,0)` (Ry) |
|---|---:|
| `8^3` | `-1.116810e-3` |
| `12^3` | `-2.717018e-5` |
| `16^3` | `-4.572961e-4` |

The non-monotonic drift confirms that mesh refinement has not yet reached the
requested stability threshold. The q = `(0,0,0.5)` control is classified as
mapped on all three cubic meshes, while q = `(0,0,0.05)` is not.

### Occupations and corrected-MFT reference

The bare reciprocal runs preserve the canonical electron count to a maximum
absolute error of `3.8488e-11` electrons. The reported Fermi levels span
`[-7.046812e-2, -6.703147e-2] Ry` across q, angle, and mesh; the variation is
retained as data rather than folded into an energy correction.

The stable corrected reference has `M_reference = 2.10370635`, all constraint
rows marked converged, maximum RMS-angle residual `4.5309e-9`, zero maximum
fixed-potential drift, and constraint penalty/field-coupling diagnostics at
about `1e-16 Ry`. It is a valid WP08 reference comparison, but it has only one
angle and one real-space mesh, so it does not establish a corrected-MFT WP09
harmonic plateau or reciprocal k-grid convergence.

## Acceptance classification

| requirement | classification | evidence |
|---|---|---|
| retain raw/gauge/pure energy definitions | established | production schema and analyzer |
| angles `2,5,10,15,20` | established for bare MFT | 39-row production sweep |
| explicit k-grid comparison | established as a diagnostic | `8^3/12^3/16^3` rows and status labels |
| bare-MFT harmonic plateau | nominally established per fixed mesh | gauge K spread below 5% |
| bare-MFT physical K convergence | not established | mesh means differ by factors of order unity |
| raw/pure finite-q convergence | not established | ~`1e-3 Ry` mesh spread and large q-only drift |
| corrected constrained-MFT characterization | reference established; WP09 cone not established | stable WP08 run plus recorded small-angle failure |
| q/2 mapping control | established for the cubic fixture | 21 mapped and 18 non-mapped rows |

The workflow therefore prevents a false WP09 pass: it exposes a useful
fixed-mesh harmonic window, but the production K and corrected small-angle
comparison remain open for the next work package or a repaired reciprocal
occupation/constraint campaign.

## Delivered files

- `source/calculation.f90`: machine-readable production WP09 diagnostic schema.
- `tools/analyze_gbt_wp09.py`: schema reader, plateau selector, mesh/drift/
  occupation/mapping analysis, JSON and text output.
- `tests/validation/val19_gbt_harmonic_kgrid.py`: reproducible bare and
  corrected-reference campaign runner.
- `tests/regression/wp9_validation/harmonic_kgrid/`: q fixtures and mapping
  documentation.
- `tests/unit/test_gbt_wp09_harmonic_kgrid.py`: synthetic analyzer and source
  contract checks.
- CTest registrations `GbtWp09HarmonicKgrid` and
  `Val19GbtHarmonicKgrid`.

Validation completed:

```text
cmake --build build -j2                         PASS
ctest --test-dir build -R '^GbtWp09HarmonicKgrid$' --output-on-failure  PASS
bare WP09 production campaign                  PASS (with non-convergence findings above)
corrected WP08 reference campaign               PASS
```
