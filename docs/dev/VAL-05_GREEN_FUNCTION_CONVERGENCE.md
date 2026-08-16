# VAL-05 — Lehmann and Dyson Green-function convergence

Date: 2026-08-16
Binary: `build-rf-debug/bin/rslmto.x`
Configuration: GNU Debug, serial CPU, full unreduced reciprocal mesh

## Scope

This validation establishes a quantitative scope for the reciprocal Green-function producer on a fixed bcc-Fe potential. It uses the production `post_processing = 'kspace_green'` route and the selected pairs:

| pair | quantity |
|---|---|
| `(1,1)` | onsite G_ii(z) |
| `(1,335)` | first-neighbour intersite G_ij(z) |

The direct records contain selected complex spin-diagonal elements and full nb-by-nb block metrics. The trace observable is `rho_i(E) = -Im Tr G_ii(E+i eta)/pi`.

The campaign samples E = (-0.8, -0.4, 0, 0.4, 0.8) Ry, sweeps Nk = 4^3, 8^3, 12^3, 16^3, sweeps eta = 0.01, 0.02, 0.04 Ry at 12^3, and compares safe energy windows `[-1.2,1.2]` and `[-1.5,1.5]` Ry. It is implemented by [`tests/validation/val05_green_convergence.py`](../../tests/validation/val05_green_convergence.py) and registered as `Val05GreenConvergence` with the `validation` label.

## Kernel invariant

The existing tight tests were not changed and passed before and after the campaign:

- `UnitLehmannChain`
- `UnitGammaSupercell`
- `UnitDysonEquivalence`
- `UnitGreenLifecycle`

`UnitDysonEquivalence` retains the chain, intersite phase, multiband block, and self-energy sign pins. In the material campaign, the maximum full-array Dyson-versus-Lehmann difference for Sigma=0 was 2.0e-11 to 9.4e-11, depending on mesh and eta. This is solver/numerical noise, not a representation discrepancy.

## Direct-GF convergence evidence

The maximum difference in selected complex G_11 elements, using the 16^3, eta=0.02 run as the campaign reference, was:

| mesh | spin-up selected-G envelope | spin-down selected-G envelope |
|---|---:|---:|
| 4^3 vs 16^3 | 3.165 | 1.696 |
| 8^3 vs 16^3 | 1.030 | 0.545 |
| 12^3 vs 16^3 | 0.512 | 0.321 |

Thus 4^3 and 8^3 are not adequate for direct metallic G_ii/G_ij values in this energy window. Even 12^3 is an observable-level intermediate mesh; 16^3 is the minimum practical reference mesh in this campaign, and a finer mesh is required for tighter pointwise accuracy or smaller eta.

The eta sweep at 12^3 gave the following maximum changes in selected G_11 values:

| comparison | spin-up | spin-down |
|---|---:|---:|
| eta=0.01 vs 0.02 | 0.558 | 0.357 |
| eta=0.02 vs 0.04 | 0.637 | 0.229 |

These are physical resolution/broadening changes in a metal. Eta must therefore be reported with reciprocal-GF data; it is not a parameter to reduce silently.

The wider safe energy window changed selected spin-up values by 0.094 onsite and 0.046 intersite at the sampled energies. This is an energy-grid/window resolution effect. The Lehmann producer sums the complete band set; it has no band truncation parameter. The real-space Chebyshev producer does use the energy window for scaling: a deliberately too-narrow `[-0.8,0.8]` Ry window triggered a floating-point failure in the fast recursion kernel and is recorded as outside the supported recursion-window regime, not normalized away.

## RS comparison

The RS and reciprocal runs used the same fixed potential, pair list, sampled energy domain, and `lld=150`; the reciprocal route used eta=0.02 Ry, the closest explicit resolution match available because the on-mesh RS intersite GF uses its continued-fraction/Chebyshev representation with zero explicit eta.

At 12^3, eta=0.02, the report gave:

- onsite full-block max |G_RS - G_Lehmann| = 9.481;
- intersite full-block max = 2.166;
- DOS max difference = 7.203, DOS RMS difference = 1.776;
- integrated trace weight: RS = 18.00002, Lehmann = 17.78673.

The trace-weight agreement is the useful broad spectral check. The large pointwise differences are not a normalization license: they combine different broadening/resolution definitions and finite recursion truncation with finite-k sampling. They are not used as a machine-equality acceptance criterion.

## Error classification and practical scope

| observed difference | classification |
|---|---|
| Lehmann vs Dyson, Sigma=0: <=9.4e-11 | solver/numerical noise; invariant retained |
| 4^3/8^3/12^3 vs 16^3 direct G | finite-k error, amplified near metallic spectral features |
| eta sweep | broadening/resolution difference |
| RS vs reciprocal pointwise GF | broadening/resolution difference plus finite recursion truncation |
| safe energy-window change | energy-grid/scaling resolution; no reciprocal band cutoff |
| too-narrow recursion window FPE | recursion-window support limitation, retained honestly |

For this fixed bcc-Fe, orthonormal-basis, full-BZ scope, 16^3 with eta=0.02 Ry is a practical starting point for direct reciprocal-GF data. This is not a general metallic-material accuracy claim and does not extend the invariant to non-orthogonal/generalized-overlap Lehmann data.

## Downstream evidence retained

The existing bcc-Fe route triads remain unchanged and remain the downstream contracts:

- `triad_bccFe_jij` — J_ij;
- `triad_bccFe_sigma` — conductivity;
- `triad_bccFe_damping` — Gilbert damping.

They retain their documented route-specific envelopes; VAL-05 does not add new exchange formulations or empirical normalization.

## Checklist

- [x] Lehmann/Dyson invariant retained
- [x] onsite G compared
- [x] intersite G compared
- [x] k-mesh convergence measured
- [x] eta convergence measured
- [x] band-window sensitivity measured where applicable
- [x] RS comparison uses matched resolution as closely as the two producers allow
- [x] existing downstream triads retained
- [x] no empirical normalization added
- [x] reciprocal GF maturity updated by scope

## Files and campaign

Changed for VAL-05:

- `source/calculation_reciprocal.f90` — direct-G report records only;
- `tests/validation/val05_green_convergence.py` — validation campaign;
- `CMakeLists.txt` — `Val05GreenConvergence` registration;
- `docs/dev/VAL-05_GREEN_FUNCTION_CONVERGENCE.md` — this evidence record;
- `docs/dev/PHASE_II_VALIDATION.md` — scoped maturity update.

Campaign run:

```text
python3 tests/validation/val05_green_convergence.py \
  --binary build-rf-debug/bin/rslmto.x \
  --scratch-root /tmp/val05_green_convergence
VAL-05 PASS
```

Commit message: `Validate Lehmann and Dyson Green functions`
