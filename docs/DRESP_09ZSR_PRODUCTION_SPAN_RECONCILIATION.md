# DRESP-09ZS-R production-span reconciliation

Verdict: **PASS-A**.

The accepted Fe audit establishes that the live production four-branch space
is a strict numerical subspace of the complete second-order six-branch space.
The required DRESP-10A migration is therefore **EXTEND**.

## Spaces and radial representation

The three spaces are distinct by construction:

```text
S4_prod = live lmto_product_response_basis phi/phidot product modes
          branches 00/10/01/11; 232 candidates

S4_2nd  = second-order branch_radial shadow
          branches 00/10/01/11; 232 candidates

S6_2nd  = second-order branch_radial shadow
          branches 00/10/01/11/20/02; 348 candidates
```

The production endpoints are

```text
e0 = phi - Enu*phidot
e1 = phidot
```

so its radial products are `e0_L e0_R`, `e1_L e0_R`, `e0_L e1_R`, and
`e1_L e1_R`, with the production origin interpolation. In contrast, the
second-order absolute-energy `00` coefficient includes the additional
`+1/2 phiddot_L*e0_R` and `+1/2 e0_L*phiddot_R` terms. The second-order `10`
and `01` coefficients likewise contain their corresponding `phiddot` terms;
the second-order `11` branch is the shared `phidot_L*phidot_R` term. Thus
`S4_prod == S4_2nd` is not assumed and is numerically false.

The primary projector is the actual initialized
`product_basis%blocks(1,K)%weighted_modes`. The secondary oracle independently
reconstructs the production candidate columns from the live `e0/e1` formulas,
performs the weighted direct SVD, and compares its retained modes and
`forward_transform` to the live basis.

## Rank and stability

All three spaces have the following candidate inventories and retained ranks;
each rank is stable at `tau1/tau10/tau100`:

| K | `S4_prod` | `S4_2nd` | `S6_2nd` |
|---:|---:|---:|---:|
| 0 | 12 / 12 / 12 | 12 / 12 / 12 | 18 / 18 / 18 |
| 1 | 16 / 16 / 16 | 16 / 16 / 16 | 24 / 24 / 24 |
| 2 | 16 / 16 / 16 | 16 / 16 / 16 | 24 / 24 / 24 |
| 3 | 8 / 8 / 8 | 8 / 8 / 8 | 12 / 12 / 12 |
| 4 | 4 / 4 / 4 | 4 / 4 / 4 | 6 / 6 / 6 |
| **total** | **232 / 232 / 232** | **232 / 232 / 232** | **348 / 348 / 348** |

The independent production oracle agrees with the live basis at every `K`:
the rank/subspace agreement is true, principal minimum singular values are
`1.0`, subspace residuals are at most `2.3e-15`, and the maximum
`forward_transform` difference is `2.3e-15`.

## Subspace measurements

The bidirectional measurements use the production weighted radial metric. The
Frobenius residuals are normalized by the square root of the number of basis
columns.

| comparison | minimum overlap singular value | max first outside second | max second outside first | Frobenius residuals (both directions) |
|---|---:|---:|---:|---:|
| `S4_prod` vs `S4_2nd` | `9.0639116468e-1` | `4.2079176312e-1` | `4.2208689955e-1` | `8.3481033180e-2`, `8.3481033180e-2` |
| `S4_prod` vs `S6_2nd` | `1.0000000000e0` | `1.9659532513e-8` | `9.9999991014e-1` | `2.3243009919e-9`, `5.7735026919e-1` |

Per-`K` principal minima and one-sided maximum residuals are:

| K | `S4_prod`/`S4_2nd` principal | prod outside `S4_2nd` | `S4_2nd` outside prod | `S4_prod`/`S6_2nd` principal | prod outside `S6_2nd` | `S6_2nd` outside prod |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | `.99526745` | `.09644621` | `.09532856` | `1.0` | `1.4968e-10` | `.99999991` |
| 1 | `.95345856` | `.25988048` | `.27293351` | `1.0` | `5.7207e-9` | `.99999903` |
| 2 | `.95856101` | `.28156052` | `.24857681` | `1.0` | `1.9660e-8` | `.99999797` |
| 3 | `.90639116` | `.42079176` | `.42208690` | `1.0` | `9.1727e-12` | `.99999917` |
| 4 | `.99849030` | `.05488555` | `.05484143` | `1.0` | `2.4834e-13` | `.99989537` |

The per-`K` maximum `S4_prod` residual outside `S6_2nd` is
`1.9659532513e-8`; the corresponding `S6_2nd` residual outside `S4_prod` is
`9.9999991014e-1`. The measured relation is therefore:

```text
                 S6_2nd
              ┌──────────┐
              │ S4_prod  │
              │          │
              └──────────┘

S4_2nd is a different, non-nested 232-dimensional space relative to S4_prod.
```

The complete second-order space contains the live production span to numerical
precision, but is materially larger than it. The six-branch physical content
outside production is significant, so this is a strict subspace extension,
not a claim that production is already complete.

## Corrected residuals against live production

The corrected `20/02` projections use `P_prod = U_prod U_prod^H`:

| diagnostic | maximum | RMS | median |
|---|---:|---:|---:|
| `R20` | `3.0599022205e-3` | `1.3777574477e-3` | `1.9709131534e-6` |
| `R02` | `9.6510267341e-3` | `4.3171570548e-3` | `2.0106157003e-6` |

The corrected full six-candidate Frobenius residual is
`1.1182586021e-3` globally. The per-`K` values are `6.0767630656e-5`,
`5.1774024009e-7`, `6.8008082059e-7`, `1.0786296504e-4`, and
`4.1800337726e-3` (the artifact contains the full
precision values). The old DRESP-09ZS values are retained in the artifact under
the explicit `S4_2nd` labels; for example, its global residual is
`4.5574142293e-4`.

The corrected arbitrary-`L` physical residuals are:

| L | physical SR maximum | physical SR RMS | delta-O maximum | delta-O RMS |
|---:|---:|---:|---:|---:|
| 0 | `1.1838878956e-6` | `1.1838878956e-6` | `1.9037710575e-2` | `1.9037710575e-2` |
| 1 | `9.2333615400e-8` | `7.4091208774e-8` | `1.2546334171e-2` | `1.0723162071e-2` |
| 2 | `2.6184031119e-8` | `1.6828056414e-8` | `1.9730900644e-2` | `1.4463253590e-2` |
| 3 | `2.1896530805e-5` | `1.6038398298e-5` | `6.4927361629e-2` | `6.2267579082e-2` |
| 4 | `7.8323815566e-1` | `2.6107949734e-1` | `6.0453937036e-1` | `2.0899094145e-1` |

Global maxima are `7.8323815566e-1` for the physical SR field and
`6.0453937036e-1` for delta-O. Corrected component maxima are:

```text
Pauli                  7.8323815566e-1
SR upper               7.8323815566e-1
SR lower-small         9.2005950131e-1
SR lower rank-0        9.3256089583e-1
SR lower rank-2        9.3256089583e-1
SR total               7.9219505541e-1
delta-O                6.0453937036e-1
```

The mixed-`L/M` physical fixture residual is `2.2619866630e-4`.

## Duality and frozen closures

The raw field/density pairings are preserved as a physical-vertex property.
For the four deterministic pairings, the raw pairing magnitudes are
`[2.7778637518e-2, 0, 0, 2.7600423411e-2]`; projected magnitudes agree to the
reported precision. Compact-coordinate duality residuals are at most
`6.25e-16`; the projected pairing loss reaches `1.4033e-12` for the active
pair. The intentional double-weighting control is `9.8853e-1`.
Projection loss is reported separately from the field/density duality error.

The frozen closures remain unchanged:

```text
DRESP-09X / DRESP-09Y / DRESP-09Z / DRESP-09ZR = preserved
ALSDA Ward = NOT RUN
BES/Halle = OFF
```

## Final classification

```text
DRESP-09ZS-R verdict = PASS-A
primary classification = PRODUCTION_SPAN_STRICT_SUBSPACE_EXTENSION_REQUIRED
Nesting = STRICT_SUBSPACE
Recommended DRESP-10A migration = EXTEND
Production modification = none
```

The implementation only initializes and reads the live production basis for
this audit. No production product dimension, Green-function moment, product
GF, product Lehmann path, angular/radial physics closure, or frozen source
file was changed.

## Reproduction

```text
cmake --build build --target rslmto.x -j2
ctest --test-dir build -R 'TddftDresp09ZSCompactSpan|Dresp09ZSFeArtifact' --output-on-failure
```

The live artifact is `/tmp/dresp09zs_fe_4k.dat`; generated Fe artifacts are not
committed.
