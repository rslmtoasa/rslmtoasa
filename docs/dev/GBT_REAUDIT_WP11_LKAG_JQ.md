# GBT WP11 re-audit: LKAG `J(q)` and frozen-spiral consistency

Date: 2026-08-27

Status: **G11 not qualified**.  The independent `J(q)` analysis path and
normalization derivation are implemented, and the bcc-Fe comparison is
reproducible.  The gate remains open because the required WP10 small-q
curvature is not k-mesh converged, and the available LKAG data contain only
one two-shell interaction range.  The measured curve and curvature therefore
do not constitute a production closure.

This is an evidence report, not a GBT gauge or physics-source change.

## Delivered workflow

The reusable analysis entry point is
[`tools/analyze_gbt_wp11.py`](../../tools/analyze_gbt_wp11.py).  It provides:

- production `jij.out` parsing with explicit mRy-to-Ry conversion;
- explicit shell multiplicities, so a representative bond is not silently
  treated as a complete Fourier sum;
- scalar one-sublattice `J(q)` construction, using the complex directed phase
  or the centrosymmetric cosine reduction;
- directed multi-sublattice `J_ab(q)` construction and a moment-weighted
  acoustic/optic matrix path;
- q-by-q curve residuals, zero-intercept quadratic/quartic curvature fits,
  q=0 Goldstone checks, and interaction-range comparison;
- machine-readable JSON suitable for future nested LKAG ranges and WP10
  campaigns.

The dependency-free regression checks are in
[`tests/unit/test_gbt_wp11_lkag_jq.py`](../../tests/unit/test_gbt_wp11_lkag_jq.py)
and are registered as `GbtWp11LkagJq` in `CMakeLists.txt`.

## Convention and normalization derivation

### What the code actually produces

The canonical exchange route in
[`source/exchange.f90`](../../source/exchange.f90) constructs the isotropic
kernel in `dGdG_Jnc`, takes `imtrace9`, integrates it with `simpson_f`, and
writes

```text
Jij_file [mRy] = Jij_internal [Ry] * 1000 / (4*pi).
```

The production `jij.out` writer emits one selected pair record as
`type_i type_j Rx Ry Rz Jij distance`; it does not emit shell multiplicities
or an executable spin-Hamiltonian consumer.  The repository's exchange
examples identify positive `J` with ferromagnetic coupling.  Consequently,
the pair-counting convention must be made explicit by the analysis rather
than inferred from the seven-column file alone.

### Convention used by WP11

For the scalar centrosymmetric comparison, WP11 adopts the standard full
directed relative-vector form

```text
H = -1/2 * sum_(i,R) J_R e_i dot e_(i+R),
J(q) = sum_R J_R exp(i*q.R).
```

The `1/2` removes the `+R/-R` double counting.  A `jij.out` row that
represents a symmetry shell therefore receives the *total* number of
equivalent directed vectors as its multiplicity; it must not receive another
factor of two merely because the cosine form is used.  For a centrosymmetric
table the implementation evaluates

```text
J(q) = sum_shell multiplicity * J_shell * cos(2*pi*q_internal.R_shell).
```

The sign is thus positive for ferromagnetic coupling, and all energies are
kept in Ry internally.  The code's `frozen_magnon` output is an energy-like
quantity in Ry; no gyromagnetic factor or conversion to a frequency is
introduced.

For a common cone,

```text
e_R = (sin(theta)*cos(q.R), sin(theta)*sin(q.R), cos(theta))
e_0 dot e_R = cos(theta)^2 + sin(theta)^2*cos(q.R).
```

With the Hamiltonian above, the energy change per chemical cell is

```text
DeltaE(q,theta) = 1/2 * sin(theta)^2 * [J(0) - J(q)].
```

The actual acoustic frozen-magnon implementation uses

```text
omega_GBT = 4*DeltaE / [M_reference*sin(theta)^2]
```

for bare MFT, where `M_reference` is the sum of the reference sublattice
moments.  Substitution gives the WP11 comparison relation

```text
omega_LKAG = 2*[J(0) - J(q)] / M_reference.
```

This factor of two is not an assumption about the Halilov formula: it follows
from the explicit full-directed-vector `1/2` Hamiltonian and the source's
explicit `4*DeltaE/(M*sin^2(theta))` output normalization.  An implementation
that instead supplies only one orientation from each pair must use half-shell
multiplicities; that representation gives the same physical result.  Using
the historical `4*[J(0)-J(q)]/M` expression together with total directed
shell multiplicities would double count the exchange.

The phase convention is the repository convention: `q` is Cartesian
`2*pi/alat` and the bond vectors printed by the bcc fixture are in `alat`
units.  Thus the phase is `2*pi*q_internal dot R_alat`, with no additional
`2*pi/alat` inserted into the dimensionless dot product.

### Multi-sublattice scope

For directed ordered records, the helper constructs `J_ab(q)`.  Under the
same Hamiltonian, the linearized moment-weighted matrix is formed from

```text
K_ab(q) = delta_ab * sum_c J_ac(0) - Re[J_ab(q)]
Omega(q) = 2*M^(-1/2) K(q) M^(-1/2).
```

The eigenvalue that vanishes at q=0 is the acoustic branch when the ordered
records are Hermitian and the reference state is a common ferromagnetic
rotation.  Reverse ordered records are required for a general
multi-sublattice table; the unit test covers the matrix and Goldstone/optic
branch algebra.  A complete B2 FeCo material comparison remains scoped until
matched, converged directed exchange data are available.  The current
single-acoustic frozen-magnon route does not establish an optic branch by
itself.

## Reproducible bcc-Fe inputs

The comparison used:

- `results/validation/VAL-18_bccFe/jij_fixed/jij.out`, the canonical two-row
  bcc-Fe LKAG output;
- shell multiplicities `8,6` for the first and second bcc shells;
- the bare reciprocal-MFT WP10 campaign in
  `/tmp/gbt-wp10-full/wp10_sweep.csv`;
- `/tmp/gbt-wp10-full/wp10_analysis.json` as the required WP10 dependency.

The command is:

```bash
python3 tools/analyze_gbt_wp11.py \
  --jij results/validation/VAL-18_bccFe/jij_fixed/jij.out \
  --multiplicities 8,6 \
  --sweep /tmp/gbt-wp10-full/wp10_sweep.csv \
  --wp10-analysis /tmp/gbt-wp10-full/wp10_analysis.json \
  --json-out /tmp/gbt-wp11-analysis.json
```

The temporary files are campaign artifacts outside the source tree; the
command regenerates the WP11 JSON when those WP10 inputs are present.

For this two-shell table,

```text
J(0) = 8*0.739154 + 6*0.467405 = 8.717662 mRy.
```

No range-convergence claim is made from this file: two rows are two shells in
one range, not two nested ranges.

## Results

The analyzer reports 108 GBT rows, spanning bare MFT meshes `8^3`, `12^3`,
and `16^3`, four cone angles, and q values `+/-0.02`, `+/-0.05`, `+/-0.10`,
and `+/-0.20` along `[001]`.

| criterion | measured result | status |
|---|---:|---|
| WP10 bare-MFT `sin^2(theta)` curvature spread | `1.277027` relative spread | **not converged** |
| LKAG interaction-range convergence | one two-shell set only | **not assessed** |
| q=0 Goldstone residual | `1.106e-13 Ry` (`1e-10 Ry` tolerance) | within tolerance |
| largest nonzero-q curve residual | `0.934822` relative | **not within tolerance** |
| full-window quartic curvature comparison | see table below | **not within tolerance** |
| overall G11 | dependency and comparisons incomplete | **not qualified** |

Representative 5-degree full-window `omega = A q^2 + B q^4` fits are shown
below.  Both coefficients use internal q units; only `A` is listed because
it is the curvature compared by the gate.

| k mesh | GBT `A` (Ry/q²) | LKAG `A` (Ry/q²) | relative difference | result |
|---|---:|---:|---:|---|
| `8^3` | 0.381286 | 0.082115 | 0.784637 | not within tolerance |
| `12^3` | 0.088921 | 0.080830 | 0.090982 | within 10% locally |
| `16^3` | 0.210452 | 0.079880 | 0.620435 | not within tolerance |

The apparent `12^3` agreement is not promoted: it is surrounded by the
`8^3`/`16^3` spread and is explicitly blocked by WP10's convergence result.
At q=`0.05` on the `12^3`, 5-degree slice, the GBT value is
`4.29928e-4 Ry`, while the two-shell LKAG prediction is `2.01060e-4 Ry`;
the relative curve residual is `0.532341`.  This also shows why one fitted
coefficient alone is insufficient for closure.

The q-reversal evidence carried from WP10 remains positive: the paired
positive/negative q odd components are numerical noise for all tested mesh
and angle groups.  It does not repair the curvature or interaction-range
failures.

## Bare MFT versus constrained/SCF response

The compared GBT rows are bare MFT: the ordinary potential is converged at
the reference point, the band-energy probe is evaluated at each q, and the
bare k-space gauge reference is subtracted.  This is the response convention
closest to standard collinear LKAG.

Corrected constrained MFT and fully self-consistent frozen magnons are not
substituted into this comparison.  They change the energy landscape and/or
include a q-dependent constraining response; a difference from bare LKAG
would not by itself identify a GBT defect.  Such modes remain separate
follow-up physics and require their own converged normalization and energy
bookkeeping.

## Acceptance checklist

| WP11 requirement | result |
|---|---|
| exact LKAG kernel/sign and output scale traced | established from `exchange.f90`; positive J is FM |
| executable spin-Hamiltonian pair convention | not encoded by the producer; explicit full-directed convention adopted and documented |
| frozen-magnon prefactor derived | established: source `4*DeltaE/(M*sin²θ)` plus the explicit `1/2` Hamiltonian gives `2[J(0)-J(q)]/M` |
| independent scalar `J(q)` | implemented and unit-tested |
| interaction-range convergence | **not established**; only one two-shell set is available |
| GBT k/q convergence | q reversal passes; WP10 curvature is **not converged** |
| bcc-Fe curve and curvature | measured; **not closed** |
| multi-sublattice matrix treatment | implemented/tested algebraically; material FeCo branch comparison scoped |
| bare/constrained/SCF distinction | documented; only bare MFT compared |
| quantitative tolerance | declared: 10% curve/curvature, 5% nested-range change, `1e-10 Ry` q=0 |
| Gate G11 | **NOT QUALIFIED** |

## Next admissible step

Repair or extend WP10 until the small-q coefficient is stable across nested
k meshes and fit windows.  Then produce at least two genuinely nested LKAG
interaction ranges (with their shell multiplicities), rerun this analyzer,
and only then add the fcc-Ni secondary or B2-FeCo matrix campaign.  No GBT
core gauge or SCF/radial change is justified by the current mismatch.
