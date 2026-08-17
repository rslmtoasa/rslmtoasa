# VAL-19: fcc-Ni transverse TDDFT response

## Result

**Failed validation gate; fcc Ni remains Experimental and is not promoted.**

This campaign followed TDDFT-07 with a one-atom fcc cell at \(a=6.650\) bohr,
no SOC, `ham_only` reciprocal response, site projection, 18 response bands,
300 K electronic smearing, \(\eta=0.001\) Ry, and a 0--0.020 Ry frequency grid
with 201 points. The sampled direct-coordinate q line was

```text
(0,0,0), (0.01,0,0), (0.02,0,0), (0.03,0,0), (0.04,0,0).
```

The raw products and machine-readable campaign record are in
[results/validation/VAL-19_fccNi](../../results/validation/VAL-19_fccNi), with
the summary in [evidence_summary.json](../../results/validation/VAL-19_fccNi/evidence_summary.json).
The deterministic checker reaches the intended scientific rejection:

```text
FAIL: q index 1: no accepted coherent-mode fit; do not use it for stiffness
```

No q² stiffness or intrinsic damping value is reported.

## Ground-state and reference status

The restarted fcc-Ni SCF state has a moment of 2.000 μB in the negative-z
orientation. A spin-reversed seed gives +2.000 μB. This is not the intended
lightweight fcc-Ni reference state (the existing example documentation records
approximately 0.595 μB), so these are implementation diagnostics rather than a
mature material benchmark. The orientation-reversed runs were retained because
they expose a separate response defect; they were not used to tune the result
to a reference curve.

## Convergence and sensitivity

| Axis | Measurements | Result |
| --- | --- | --- |
| k mesh | 4³, 8³, 12³ | Raw legacy Goldstone residual is 0.03237984 in all three runs within printed precision; the sampled response is likewise unchanged. This is numerical convergence of the present diagnostic, not physical validation. |
| Band window | bands 1:18 and 1:14 at 8³ | No change in the sampled response at the printed precision. |
| Smearing | 100, 300, 600 K | The apparent legacy q=0.01 feature moves from 0.0040287 to 0.0042178 to 0.0074031 Ry. It is therefore not a smearing-insensitive mode. |
| η | 0.0005, 0.001, 0.002 Ry | The observed q=0.01 legacy-trace feature widths are 0.0011, 0.0021, and 0.0041 Ry. The width follows numerical broadening; no intrinsic or Landau damping value is claimed. |
| Projection/basis | site projection; legacy and pair-potential Xi routes | Site projection is the current production-supported projection in this campaign. The two Xi routes are retained as a physically relevant implementation comparison; no alternate production basis was available. |

The response convention is recorded in each manifest: frequencies and η are in
Ry, q is in direct reciprocal-lattice coordinates, and the loss is the output
trace convention \(-\operatorname{Im}\chi/\pi\). These definitions are aligned
internally across the matrix.

## Stoner continuum, loss, and Xi/mode evidence

The bare \(\chi_0\) trace is already nonzero at the lower frequency edge for
finite η. No clean finite Stoner onset can therefore be resolved on this grid;
the continuum is treated as overlapping the low-energy window rather than being
assigned a sharp onset.

At q=0.01, the legacy route reports a formal
\(\operatorname{Re}\lambda_\Xi=1\) crossing near 0.0042178 Ry, but its imaginary
part is -0.0074219 Ry, its projected enhanced weight is negative
(-2.539×10³ in the recorded output units), and the mode fit is rejected because
there is no isolated local maximum. The mode file consequently classifies this
as a strongly damped/incoherent enhanced feature, not a magnon.

The pair-potential route has no accepted Xi crossing and classifies its feature
as noncollective Stoner evidence. The dynamic legacy and pair-potential loss
products both contain negative spectral weights. They are retained as a
causality/response-representation failure; they are not clipped, normalized, or
reinterpreted as positive physical loss.

The static orientation probe gives:

| Orientation | Moment | Pair Xi Goldstone residual | Closest pair Xi eigenvalue |
| --- | ---: | ---: | ---: |
| negative z | -2.000 μB | 2.000000 | -1.000000 |
| positive z | +2.000 μB | 1.33×10⁻¹⁵ | +1.000000 |

This identifies an orientation-dependent pair-Xi sign defect. The positive-z
dynamic products still have negative spectral weights and no coherent mode, so
fixing the static sign would not by itself establish a Ni magnon.

## Small-q decision and linewidth

No coherent quadratic branch can be identified over the tested q interval.
The rejection is not explained by an unresolved k-mesh trend: the 4³--12³
diagnostics are stable. The band-window sweep is also unchanged. Smearing moves
the apparent trace feature, while the overlapping Stoner response, negative
spectral weights, and absent positive-weight Xi mode prevent a collective-mode
interpretation. The static pair-Xi orientation defect is an additional
implementation blocker.

Because no q point has an accepted coherent fit, no q² fit is made. A quadratic
fit through the rejected trace peaks would be scientifically misleading.

The η sweep is an observed trace-width sweep, not a collective-mode linewidth
sweep. Since there is no accepted Xi mode and the width scales with η, the
campaign makes no intrinsic/Landau damping claim.

## Reference comparison and maturity

GBT, Jij, and trusted-literature comparisons are deferred. The internal units,
direct reciprocal coordinates, and loss definition are explicit, but there is
no valid collective branch and the current 2 μB reference state is not a mature
fcc-Ni benchmark. No code parameter was tuned to match literature.

Ni maturity remains **Experimental**. Before promotion, the ground-state
reference, orientation-covariant pair Xi, response spectral sign, Stoner/mode
separation, and an independently comparable coherent branch must be repaired
and revalidated.

## TDDFT-07 checklist

- [x] k-mesh convergence measured.
- [x] Band-window sensitivity measured.
- [x] Smearing sensitivity measured.
- [x] η sensitivity measured.
- [x] Stoner continuum mapped as overlapping at the lower finite-η grid edge.
- [x] Collective-feature absence explained.
- [x] Small-q coherence assessed.
- [x] q² fit used only where justified; no fit was made.
- [x] Linewidth interpretation supported by an η sweep; no intrinsic claim made.
- [x] Literature comparison kept to aligned definitions; comparison deferred because no valid branch exists.
- [x] Ni maturity updated conservatively; it remains Experimental.
