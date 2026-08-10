# TDDFT-07: transverse LR-TDDFT validation campaign

## Current status

The validation machinery is in place and the small deterministic gate passes.
This is not a declaration that a material benchmark has passed: no bcc-Fe,
fcc-Ni, corrected-GBT, or `Jij` material output is committed with this change.
The transverse milestone therefore remains open until the real evidence below
is produced and reviewed.

- [x] Independent brute-force `chi_KS` oracle (`UnitTddftChiKS`).
- [x] Deterministic campaign fixture checks raw and corrected Goldstone records,
  KS/Stoner and enhanced spectra, quadratic fitting, η extrapolation, all
  convergence-axis names, and read-only GBT/`Jij` comparison inputs.
- [x] CI registration: `UnitTddftValidationCampaign` (`unit;tddft;validation;quick`).
- [ ] Raw no-SOC material `q=0, omega=0` residual and unity-mode/magnetization
  overlap are stable under response-basis, k-mesh, and band-window sweeps.
- [ ] Stable three-route small-q stiffness for a real material.
- [ ] bcc Fe credible magnon and KS/Stoner spectra, with damping evolution.
- [ ] fcc Ni credible magnon and KS/Stoner spectra, with itinerancy sensitivity.
- [ ] Optional Co campaign.

The implementation must not be marked transverse-complete while an unchecked
material gate remains unchecked.

## Campaign procedure

Start from a converged collinear no-SOC ground state and retain the original
response output files. Run `post_processing = 'susceptibility'` with
`goldstone_mode = 'diagnose'` for quality measurements. A separate
`goldstone_mode = 'correct'` run with `xi_backend = 'pair_potential'` (or
`'compare'`) may emit a controlled comparison spectrum. It derives real
static column scales from direct pair-potential `Xi`, rejects ill-conditioned,
complex-static, small-moment, and over-25%-change cases, and writes
`*_pair_dyson.dat` and `*_pair_corrected_dyson.dat` separately. The raw
`*_goldstone.dat` must be retained next to the corrected output: the raw
residual is the quality metric and correction is not evidence of raw
convergence. The old `goldstone_mode = 'sum_rule'` spelling migrates to
`correct` with a warning.

Make one `campaign.json` for every fixed physical setup. The checker at
[`tests/regression/tddft_validation/tddft_validation.py`](../../tests/regression/tddft_validation/tddft_validation.py)
reads the production text outputs and writes a compact, reproducible evidence
summary:

```bash
python3 tests/regression/tddft_validation/tddft_validation.py results/bccFe/campaign.json \
  --report results/bccFe/evidence.json
```

The accepted stiffness is the zero-intercept least-squares fit
`omega = D |q|^2`; the report includes the relative quadratic residual and
does not silently include failed/incoherent fits. `GBT` and `Jij` values may
be supplied in `independent_routes`, but they remain values from separately
run methods. TDDFT inputs are fixed before their comparison is inspected.

For a linewidth, record three or more values of the **observed** FWHM at fixed
physics and varied η. The reported intercept of a linear FWHM(η) regression is
the zero-η estimate only when the fit residual is acceptable. A finite width
at one η is numerical broadening, not Landau damping. A credible Landau-damping
claim additionally needs a collective Xi candidate, a resolvable enhanced
feature, and its evolution against the separately output KS/Stoner continuum.

## Required real-material matrix

Use the same response convention and fixed reference potential within each
row. High-accuracy artifacts belong in the non-CI benchmark store, not in the
small fixture directory.

| System | Required products | Minimum independent checks |
| --- | --- | --- |
| bcc Fe | `chi0`, Dyson/loss, mode file, Goldstone file, small-q fit | corrected GBT/frozen magnon and `Jij` stiffness, η sweep, all convergence axes |
| fcc Ni | same products | explicitly show sensitivity to band window, smearing, η, and Stoner onset before interpreting a mode |
| hcp Co (optional) | same products | add only after the Fe/Ni workflow is stable |

For every convergence point preserve: k mesh, selected band interval,
response projection, electronic temperature/smearing, η, frequency range and
count, q list, XC functional, lattice parameter, and whether controlled
pair-potential correction was applied. Do not compare a corrected TDDFT curve
with an uncorrected Goldstone quality metric.

## Literature frame and expected discrepancies

The qualitative reference for Fe/Co/Ni is Buczek, Ernst, and Sandratskii,
*Different dimensionality trends in the Landau damping of magnons in iron,
cobalt and nickel* ([arXiv:1109.6217](https://arxiv.org/abs/1109.6217)); it
provides the relevant LR-TDDFT distinction between a collective mode and the
Stoner continuum. The related local-susceptibility/Goldstone treatment is
described by Lounis et al. ([arXiv:1006.0963](https://arxiv.org/abs/1006.0963)).

Comparison is deliberately qualitative or semi-quantitative until like-for-like
conditions are established. Differences can arise from ASA versus full-potential
treatment, XC functional and lattice constant, site-only versus richer response
projection, k and unoccupied-state convergence, and η/frequency resolution.
None of these is a license to retune a TDDFT control after seeing GBT, `Jij`, or
literature values. Report them as uncertainty/provenance instead.

## Definition-of-done ledger

- [x] Toy oracle and deterministic CI fixture.
- [ ] Material Goldstone diagnostics pass with raw quality reported separately.
- [ ] Small-q LR-TDDFT stiffness is stable and independently compared to corrected GBT and `Jij`.
- [ ] At least one itinerant ferromagnet has a credible magnon plus Stoner spectrum.
- [x] The campaign checker prohibits interpreting a lone finite-η width as physical damping.
