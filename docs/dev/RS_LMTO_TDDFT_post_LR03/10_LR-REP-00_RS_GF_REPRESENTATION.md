# LR-REP-00 — Certify native RS Green-function representation transformations

## Timing

This task is **not required** for the first reciprocal k-space TD-DFT validation.

Run before implementing a native real-space GF susceptibility backend.

## Goal

Establish exactly how the native recursion/Chebyshev Green functions relate to:

- orthogonal Hamiltonian resolvent;
- physical LMTO Green function;
- screened/auxiliary/path-operator representations used by legacy routines.

No susceptibility implementation.

## Primary targets

Audit the live equivalents of:

- `auxiliary_gij`;
- `transform_auxiliary_gij`;
- any `sqrt(delta)` endpoint factors;
- screening-constant transformations;
- onsite additive terms.

## Derive

For every transformation identify:

- input representation;
- output representation;
- site/orbital ordering;
- onsite versus offsite formula;
- energy dependence;
- screening parameters;
- exact inverse where one exists.

Use published LMTO identities as the formal target.

## Required tests

### Round trip

For deterministic matrices:

\[
g^\alpha\rightarrow g^\beta\rightarrow g^\alpha.
\]

### Onsite negative control

Deliberately omit the additive onsite screening term in the test oracle and require
the comparison to fail.

This protects against an offsite-only formula accidentally passing.

### Reciprocal/real-space bridge

Where both representations are available for a small periodic fixture, Fourier
transform the certified coefficient-space reciprocal GF and compare with the
native RS representation after the derived transformations.

Do not compare matrices with different representation labels directly.

### High-energy behavior

Check representation-specific asymptotics.

## Terminator policy

The accepted block terminator is not under audit.

Do not reopen it.

## Deliverable

`docs/LR_RS_GF_REPRESENTATION_AUDIT.md`

## PASS criterion

Every representation used by a future RS response backend must have an explicit,
tested conversion path into the same physical/Pauli augmented endpoint object used
by the reciprocal response.

If not, declare the affected RS backend BLOCKED.

## Commit

`docs: certify real-space Green-function representations`

## Completion checklist

- [x] live `auxiliary_gij` and `transform_auxiliary_gij` paths audited;
- [x] endpoint, ordering, screening, onsite, inverse, and energy contracts documented;
- [x] round-trip and onsite-negative-control tests pass;
- [x] reciprocal/RS coefficient-space bridge passes after matching representation conversions;
- [x] high-energy behavior passes for raw, endpoint-scaled, and screened paths;
- [x] accepted terminator left unchanged;
- [x] affected native RS response backend explicitly declared `BLOCKED` pending radial/Pauli endpoint augmentation;
- [x] deliverable `docs/LR_RS_GF_REPRESENTATION_AUDIT.md` created;
- [x] commit prepared with the prescribed message.
