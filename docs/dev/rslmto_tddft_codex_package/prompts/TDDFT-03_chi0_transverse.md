# Codex Prompt — TDDFT-03: Bare Transverse Kohn-Sham Susceptibility

## Mission

Implement the first real physics output: `chi_KS^{+-}(q,omega)` and the corresponding Stoner spectrum.

Use TDDFT-01 eigenpairs and TDDFT-02 transition vertices.

Read first:

- `README.md`
- `00_MASTER_BLUEPRINT.md`
- `01_PHYSICS_CONVENTIONS.md`

## Reference formula

Implement the finite-temperature occupation-difference form:

\[
\chi_{\rm KS}^{AB}(\mathbf q,\omega)
=
\frac{1}{N_k}\sum_{\mathbf k n m}
\frac{f_{n\mathbf k}-f_{m,\mathbf k+\mathbf q}}
{\omega+\epsilon_{n\mathbf k}-\epsilon_{m,\mathbf k+\mathbf q}+i\eta}
M^A_{nm}M^{B*}_{nm},
\]

translated to the exact internal energy/frequency convention established by TDDFT-00.

## Initial scope

- collinear;
- no SOC;
- transverse circular channel;
- site-projected response.

The implementation engine should nevertheless consume generalized vertices rather than hard-coded spin-up/down band lists.

## Tasks

- [x] Implement Fermi occupations using existing Fermi level and electronic smearing/temperature conventions.
- [x] Use the `f_n-f_m` all-band reference expression.
- [x] Make eta explicit and user-configurable.
- [x] Batch over frequency because vertices/numerators are frequency independent.
- [x] Accumulate response matrices on the fly.
- [x] Add optional band/occupation pruning only after the exact reference path passes.
- [ ] Output at least:
  - [x] `Re chi_KS`;
  - [x] `Im chi_KS`;
  - [x] site-diagonal spectra;
  - [x] trace spectral quantity;
  - [x] a clearly labeled KS/Stoner spectral map.
- [x] Write eta, units, k-mesh, band count/window, smearing and Fermi level to output metadata.
- [x] Define/test the positive-frequency spectral sign convention.

## Independent oracle

Create a one-orbital spin-split tight-binding fixture with a deliberately independent brute-force pair sum.

- [x] compare several q;
- [x] compare several omega;
- [x] compare finite smearing;
- [x] compare multiple eta;
- [x] require machine-precision agreement for tiny cases.

## Convergence tests

- [x] k mesh;
- [x] band count/window;
- [x] electronic smearing;
- [x] eta.

## Non-goals

- no ALDA enhancement;
- no Goldstone correction;
- no linewidth claim from eta.

## Completion protocol

Tick checkboxes and commit once:

`TDDFT-03: implement bare transverse Kohn-Sham susceptibility`
