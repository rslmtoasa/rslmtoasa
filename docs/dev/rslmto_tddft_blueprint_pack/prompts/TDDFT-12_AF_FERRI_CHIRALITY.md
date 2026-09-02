# TDDFT-12 — AF/ferrimagnetic channels and chirality validation

## Agent mission

Extend the validated collinear transverse implementation from simple ferromagnets to multi-sublattice AF/ferrimagnetic systems while keeping `+-` and `-+` channels physically distinct.

## Context

Halle work on FeRh and CrMnSb shows that opposite precession/chirality channels can have qualitatively different Landau damping. This is a direct test of whether the response implementation has hidden FM-only assumptions.

## Non-negotiable constraints

- Work on the current development branch supplied by the user. Do not assume historical TD-DFT fixes are present.
- Preserve the codebase’s modular/class-like Fortran style; prefer existing type/factory patterns over gratuitous new abstractions.
- Do not edit legacy atomic LMTO routines unless the task demonstrates that the response derivation requires it.
- Do not weaken tests or “fix” spectra with empirical shifts.
- Keep unsupported response combinations explicitly rejected rather than silently approximated.
- Tick checklist boxes only after concrete evidence exists.

## Work plan

1. Audit all FM-specific assumptions in projection, mode tracking, moment signs and Goldstone vector construction.
2. Generalize the static Ward identity to signed sublattice magnetizations in the chosen response basis.
3. Preserve and output separate circular channels.
4. Validate FeRh FM and AF phases: acoustic/optical character, branch count, mode chirality and damping trends.
5. Validate CrMnSb or an equivalent compensated half-metal fixture for chirality-dependent access to the Stoner continuum.
6. Compare adiabatic low-energy branches with independent exchange/LSWT references where meaningful.
7. Add multi-sublattice regression fixtures that would fail if site moments are replaced by absolute values.

## Acceptance checklist

- [x] No hidden absolute-moment/FM assumption remains; signed `P_site sigma_z` is separate from the legacy magnitude population.
- [x] `+-` and `-+` are separate in data/output, with explicit channel metadata and reverse-channel filenames.
- [ ] FeRh FM/AF branch structure is physically sensible; the named-material campaign is retained as an optional higher-cost case because no FeRh input/evidence set is present in this checkout.
- [x] CrMnSb-like chirality asymmetry is reproduced qualitatively under converged settings by the committed compensated two-sublattice fixture; the named-material campaign remains optional.
- [x] Multi-sublattice Goldstone test passes with a signed `[+2,-1]` vector.
- [x] Backend equivalence still holds for both ordered channels on the AF/ferri fixture.

## Required evidence / deliverables

Add a concise validation report and retain the literature cases as optional higher-cost regression examples.

Validation report: `docs/TDDFT_AF_FERRI_CHIRALITY.md`.

## One-line commit message

`tddft: validate ferri and antiferromagnetic chirality`
