# RA4 spectral impurity baseline preregistration

Registered before RA4 runs. Source baseline is `8b8f7ef7bfc026fd34a56b37a3ad3281fd0a6832`; the existing `build/bin/rslmto.x` is used without rebuilding while shared source edits are in progress. Raw inputs, runs, and outputs are kept under `/tmp/rslmto-ra4-*`; compact results and hashes will be recorded in `RA4_RESULTS.md`.

## Prediction and decision

In spectral mode, an RS impurity embedded in a same-material RS bulk host should reproduce the host d-shell trace per spin, sphere moment, and d-projected LDOS within the corresponding U=0 impurity-minus-bulk discrepancy plus measured repeat scatter. Test bcc Fe at U=3.8 eV/J=1.1 eV and hcp Co at U=3 eV/J=0.7 eV. Run U=0 controls at the same geometry and compare two impurity starts: a no-U potential seed and a +U potential seed. A disagreement exceeding the U=0 discrepancy plus repeat scatter is a measured failure of the spectral baseline. If the control discrepancy is smaller than repeat scatter, mark that observable's control null and use the independent printed occupation matrix as the occupation oracle. Missing matching-material geometry or input/seed assets block a baseline claim; no tolerance is to be widened after observing outputs.

## Error budget and measured quantities

For scalar observables, define the U=0 control discrepancy as `abs(I_0 - B_0)` using the same observable extraction for impurity and bulk. Define repeat scatter as the maximum absolute difference among repeated runs from the same input and potential seed. The allowed U>0 discrepancy is the corresponding U=0 discrepancy plus that measured scatter. This is tied directly to the decision whether adding U changes impurity/bulk equivalence beyond known embedding mismatch and numerical reproducibility. Use identical `channels_ldos` grids and energy bounds for paired bulk/impurity runs; compare the d-projected LDOS on the common energy grid and report the max absolute difference and its maximum U=0 control value plus repeat scatter. For occupation matrices, independently parse the emitted potential and evaluate the d trace per spin; no independent tolerances are added beyond the observed control and scatter.

## Run design

Use existing bcc Fe and hcp Co RS material structures and existing single-site impurity templates. Derive scratch inputs by changing only scalar spin setting where required, fixed-U values, and atom labels needed to represent the same element at the impurity site. Do not alter geometry or unsupported guards. Start bulk U=0 and U>0 runs from matching no-U material potentials; use the converged bulk potential files as frozen hosts for the paired impurity runs. For each U>0 impurity, run once with the no-U impurity potential seed and once with a +U impurity potential seed; repeat the same input/seed pair to measure repeat scatter. If any required matching host/impurity asset or run support is absent, stop that material's campaign and report the precise gap.

## Provenance

Each run record will include the source commit, executable SHA256, generated input and potential-seed SHA256 values, output path and SHA256, and exact command. This preregistration itself has SHA256 recorded with the results.
