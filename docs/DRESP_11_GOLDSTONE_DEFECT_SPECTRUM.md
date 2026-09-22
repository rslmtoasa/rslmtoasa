# DRESP-11 — Goldstone defect spectrum and restoration admissibility

Starting HEAD: `544a708a6df637f8a17662d1a653c0f64cf140ae`
Status: diagnostic-only; no correction is enabled in production dynamics.

## Scope

DRESP-11 assembles the exact static fixed-basis operator in the live 348-dimensional Pauli product space,

\[
A = R_P L_f(H) F_{SR} K_{xc} D_P,\qquad D=I-A.
\]

The authoritative action is the existing DRESP-10F chain. Dense assembly is performed by applying that
matrix-free action to every canonical compact basis vector. The completed dense columns are then traversed
in reverse order as a deterministic assembly-order check, while a separate seven-vector oracle reapplies the
matrix-free action. The raw source insertion, six-branch basis, static Fréchet response, P3 target, and ALSDA
kernel are unchanged.

## Required gates

The campaign consumes the completed frozen DRESP-10F and DRESP-09U/Y reports produced by their independent
CTest gates, then checks matrix/action agreement, P3 reconstruction, SVD/eigen diagnostics, a compact-block
unitary gauge transformation, and the explicit diagnostic corrections. The conservative isolation gates are

\[
|\hat m_G^\dagger v_0|\ge 0.99,\qquad \sigma_1/\sigma_0\ge 50,
\qquad \|Q_GDm_G\|/\|Dm_G\|\le 0.1.
\]

Failure of these gates is reported as a distributed finite-LMTO Ward defect; thresholds are not tuned.

## Machine-readable result

The accepted Fe run writes the report and spectrum table to:

```text
/tmp/dresp11_fe_4k.dat
/tmp/dresp11_fe_4k.dat.spectrum.csv
```

The CSV contains the lowest 20 singular values in ascending order and the lowest relevant eigenvalues
ordered by magnitude, with rigid-mode overlaps. The textual artifact contains the raw Ward residual,
defect decomposition, matrix character, gauge stability, correction norms, and corrected spectra.

## Accepted 4×4×4 Fe result

The optimized accepted-state run completed in 1126.6 s and reports:

- `PASS-B`, `DISTRIBUTED_FINITE_LMTO_WARD_DEFECT`.
- (|D m_G|/|m_G|=2.3630	imes10^{-1}), with physical weighted radial residual (2.3630	imes10^{-1}).
- (sigma_0=1.8406847	imes10^{-2}), (sigma_1=5.6599442	imes10^{-1}), and (sigma_1/sigma_0=30.7491).
- The smallest-right-singular-vector overlap with (m_G) is (0.959126), while the orthogonal Ward fraction is (0.906583).
- The compact gauge check is closed; matrix-free/assembled action and denominator residuals are (3.58	imes10^{-15}) and (4.02	imes10^{-16}).
- Rigid, SVD, and BES corrections are `NOT_RUN`; production correction and dynamics remain disabled.

The isolation gates therefore remain open, and the next milestone is finite-LMTO covariance rather than a Goldstone correction or dynamics run.

## Production boundary

`BES/Halle production = OFF`, `Correction status = DIAGNOSTIC_ONLY`, and `Dynamics = NOT RUN` are hard
report invariants. Juelich/GSR is labeled `NOT_COMPARABLE` because its site-local interaction is not
embedded into the full 348-dimensional spatial kernel.
