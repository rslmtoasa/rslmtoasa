# TDDFT-06 — implement K-space GF chi0 backend

## Agent mission

Implement `chi0` from products of validated K-space Green functions while using the exact same response vertices/projectors and kernel conventions as the eigenpair backend.

## Context

This backend is a numerical bridge: it shares eigenpairs internally with the Lehmann GF but exercises the GF convolution/integration formulation independently from the explicit transition-denominator code.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Implement the retarded/advanced transverse GF convolution in k space for the supported collinear SOC-free case.
2. Start with a direct near-real-axis energy integration mode using explicit eta for transparent validation.
3. Use k and k+q with exact endpoint provenance and q mapping.
4. Reuse the common source/measurement vertices; do not create a second convention.
5. Compare matrix-valued chi0 against the explicit eigenpair backend on analytic fixtures and small Fe/Ni calculations.
6. Add frequency/channel symmetry tests.
7. Profile and batch over k and energy nodes where straightforward, but keep the implementation readable.

## Acceptance checklist

- [x] K-GF chi0 passes analytic fixtures.
- [x] K-GF and eigenpair chi0 agree pointwise under matched eta/temperature.
- [x] Static chi0(0,0) agrees.
- [x] Ward residual computed from K-GF matches eigenpair behavior.
- [x] No backend-specific kernel convention was introduced.
- [x] Provenance records energy integration settings.

## Required evidence / deliverables

Provide pointwise comparison tables for selected complex chi0 matrix elements, not just plotted spectra.

Evidence: `docs/dev/rslmto_tddft_blueprint_pack/TDDFT-06_KSPACE_GF_CHI0_EVIDENCE.md` records the analytic/static/Ward fixtures and pointwise complex comparison tables. The executable is `UnitTddftGreenChiKS`; the source contract is additionally checked by `UnitKspaceGFValidationSource`.

## One-line commit message

`tddft: add k-space Green-function chi0 backend`
