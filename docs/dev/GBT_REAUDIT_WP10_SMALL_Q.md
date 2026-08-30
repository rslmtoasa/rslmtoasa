# GBT WP10 re-audit: small-q quadratic limit and q-reversal symmetry

Date: 2026-08-27

Status: q-reversal symmetry is established for the bare reciprocal MFT
fixture.  A converged small-q curvature is **not established**: the fitted
coefficient changes substantially between the tested k meshes, and the
fit-window reduction exposes strong higher-order/non-converged structure.

This is a measurement report, not a GBT-gauge change.  No SOC, empirical
energy shift, angle rescaling, or fitted correction was used.

## Workflow delivered

`tools/analyze_gbt_wp10.py` is the public WP10 entry point.  It reads the WP09
machine-readable sweep schema and performs the following without discarding
raw rows:

- pairs positive and negative q values and reports
  `E_even=(E(+q)+E(-q))/2` and `E_odd=(E(+q)-E(-q))/2`;
- fits the zero-intercept models `DeltaE = A q^2 + B q^4` and
  `DeltaE = A q^2`;
- reports coefficient uncertainties from the least-squares residual and
  fit-window stability as the largest q is reduced;
- reports raw and `sin(theta)^2`-normalized fits, plus mesh/cone sensitivity;
- keeps internal q coefficients and, when `alat` is known, converts them to
  inverse-Angstrom q coefficients.

The implementation is shared with `tools/analyze_gbt_wp09.py` through its
`analyze_wp10_rows` API.  The dedicated entry point makes the workflow
explicit while preserving the WP09 reader and output contract.

The q-sweep runner is
`tests/validation/val20_gbt_small_q.py`.  Its committed fixture starts at
Gamma and contains the matched pairs `+/-0.02`, `+/-0.05`, `+/-0.10`, and
`+/-0.20` along `[001]`.  It writes raw `wp10_sweep.csv` before the analysis
JSON.  CTest entries `GbtWp10SmallQ` and `Val20GbtSmallQ` are registered.

## Conventions

The code-internal q coordinate is Cartesian `2*pi/alat`, matching
`hamiltonian%q_ss`.  For the bcc-Fe fixture, `alat = 2.86120 Angstrom`, so

```text
q_physical [1/Angstrom] = 2*pi/alat [1/Angstrom] * q_internal
                         = 2.1959965424 * q_internal.
```

For bare reciprocal MFT the primary fit quantity is the harmonic response

```text
DeltaE(q,theta) = [E(q,theta) - E(q,theta=0)]
                 - [E(Gamma,theta) - E(Gamma,theta=0)].
```

This is the q-dependent version of the WP09 same-q gauge subtraction.  The
raw q-reference difference is also fitted and retained for comparison.  The
primary `sin2` quantity divides `DeltaE` by `sin(theta)^2`; this makes cone
angle dependence visible rather than silently folding it into A.

## Production campaign

The bare-MFT campaign used the scalar-relativistic reciprocal bcc-Fe fixture,
full Brillouin-zone integration, q direction `[001]`, cone angles 2, 5, 10,
and 15 degrees, and meshes `8^3`, `12^3`, and `16^3`:

```text
python tests/validation/val20_gbt_small_q.py \
  --binary build/bin/rslmto.x \
  --scratch-root /tmp/gbt-wp10-full \
  --angles 2 5 10 15 \
  --meshes 8 12 16 \
  --modes mft
```

The run produced 108 raw rows in
`/tmp/gbt-wp10-full/wp10_sweep.csv` and analysis JSON in
`/tmp/gbt-wp10-full/wp10_analysis.json`.  The paths are intentionally outside
the source tree; rerunning the command regenerates the complete raw table.

## q-reversal result

All 12 mesh/angle groups contained four complete positive/negative pairs.
The largest absolute odd component in each mesh was:

| mesh | max `|E_odd|` (Ry) | declared tolerance | result |
|---|---:|---:|---|
| `8^3` | `7.22e-15` | `1.00e-10 Ry + 1.00e-6 * max|E_even|` | within tolerance |
| `12^3` | `1.19e-14` | same | within tolerance |
| `16^3` | `1.97e-14` | same | within tolerance |

This supports q-sign, bond-reversal, and Fourier symmetry for this
SOC-free reciprocal fixture.  It does not establish curvature convergence.

## Fits

The table below gives the full-window `sin(theta)^2`-normalized
`A q^2 + B q^4` fit.  A is in Ry per internal-q-squared and B is in Ry per
internal-q-fourth.  The 5-degree uncertainties are shown because they are a
representative harmonic-angle slice; every angle and the raw fits remain in
the JSON file.

| mesh | theta | A +/- uncertainty (Ry) | B +/- uncertainty (Ry) |
|---|---:|---:|---:|
| `8^3` | 2° | 0.196313 +/- 0.035083 | -2.86877 +/- 0.903 |
| `8^3` | 5° | 0.196061 +/- 0.035217 | -2.86581 +/- 0.907 |
| `8^3` | 10° | 0.195924 +/- 0.035548 | -2.87208 +/- 0.916 |
| `8^3` | 15° | 0.197287 +/- 0.035810 | -2.91727 +/- 0.922 |
| `12^3` | 2° | 0.046544 +/- 0.009829 | 0.69478 +/- 0.253 |
| `12^3` | 5° | 0.046451 +/- 0.009867 | 0.70906 +/- 0.254 |
| `12^3` | 10° | 0.046570 +/- 0.009902 | 0.73473 +/- 0.257 |
| `12^3` | 15° | 0.047914 +/- 0.009713 | 0.70414 +/- 0.251 |
| `16^3` | 2° | 0.111401 +/- 0.017367 | -2.39653 +/- 0.378 |
| `16^3` | 5° | 0.111244 +/- 0.017437 | -2.36403 +/- 0.449 |
| `16^3` | 10° | 0.110829 +/- 0.017656 | -2.25354 +/- 0.456 |
| `16^3` | 15° | 0.110848 +/- 0.017883 | -2.12593 +/- 0.435 |

At 5 degrees the same three representative physical A coefficients are
`0.040656 +/- 0.007303`, `0.009632 +/- 0.002046`, and
`0.023068 +/- 0.003616 Ry Angstrom^2` for `8^3`, `12^3`, and `16^3`,
respectively.  The `2*pi/alat` conversion is therefore visible and does not
change the non-convergence conclusion.

Fit-window stability is also negative.  For example, at 5 degrees the
normalized quartic A changes from `0.38689` to `0.19606` as the `8^3` maximum
q changes from 0.10 to 0.20, from `0.10177` to `0.04645` on `12^3`, and from
`0.20733` to `0.11124` on `16^3`.  The utility reports these windows and marks
their A stability `unstable`; no single window is promoted to a stiffness.

The fixed-mesh cone-angle variation of the full-window normalized A is much
smaller than the mesh variation (about 0.7%, 3.1%, and 0.5% for `8^3`,
`12^3`, and `16^3`), but that does not rescue the mesh dependence.  The raw
fit is strongly cone-angle dependent, as expected before the harmonic
`sin(theta)^2` normalization.

## Mode separation and acceptance

This WP10 closure is bare MFT, the method intended for the WP11 LKAG
comparison.  The analyzer accepts corrected-MFT and SCF rows and keeps their
raw definition separate, but no corrected-MFT or fully self-consistent
small-q curvature is claimed here.  WP09’s corrected-MFT evidence remains a
large-angle reference run, not a small-angle harmonic closure.

| checklist item | result |
|---|---|
| symmetric +/- q points | established |
| even and odd components | established |
| odd component within declared tolerance | established for 12/12 groups |
| q units and physical conversion | established |
| q2+q4 and pure-q2 fits | established |
| fit-window stability | assessed; unstable |
| k-grid sensitivity | assessed; not converged |
| cone-angle sensitivity | assessed; normalized cone variation smaller than mesh variation |
| bare/corrected/SCF distinction | explicit; WP10 data are bare MFT only |
| stable long-wavelength coefficient | **not established; failure documented** |

The result is therefore a useful WP10 diagnostic and a valid q-reversal pass,
but not Gate G10 closure.  WP11 must not compare LKAG against one of these
mesh-dependent coefficients until the small-q/k-grid campaign is repaired or
extended.
