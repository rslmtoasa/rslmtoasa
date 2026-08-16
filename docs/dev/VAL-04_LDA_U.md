# VAL-04 — Validate onsite LDA+U physics and support scope

## Result

The ordinary fixed-U/J Liechtenstein route is validated for the scoped
collinear bcc-Fe reciprocal-SCF envelope below. This is a physical response
validation of the onsite correction, not a claim that every magnetic/SOC route
or the total-energy functional is mature.

The current Hubbard-V route remains **Development**. It uses a local,
diagonal occupation proxy; the intersite density matrix needed for full
DFT+U+V is not wired and is not promoted here.

## Current onsite formulation

`source/hamiltonian_hubbard.f90::calculate_hubbard_u_potential_general`
constructs a real, spin-resolved orbital potential from
`ldm(l,ispin,m,m')`, with Coulomb direct/exchange factors generated from the
input U/J. For the default `hubbard_u_potential_form='liechtenstein'`, the
diagonal part includes the implemented FLL-like shift

\[
  -U(n-1/2) + J(n_\sigma-1/2),
  \qquad n=n_\uparrow+n_\downarrow.
\]

This is the double-counting convention actually used by the Hamiltonian
route. U and J are entered in eV and converted to Ry during namelist load.
The correction is added to the onsite block in `build_bulkham` and
`build_locham`; reciprocal assembly Fourier-transforms the same onsite block.

The current Hamiltonian path does not accumulate a separate scalar
`E_U-E_dc` contribution. Consequently, `etot` and `sumev` below are finite
total/band-energy response observables under the active Hamiltonian
convention. They must not be interpreted as validation of a fully variational
DFT+U total-energy functional.

## Evidence and commands

The compact physical campaign is
[`tests/validation/val04_lda_u.py`](../../tests/validation/val04_lda_u.py),
using the existing reciprocal bcc-Fe fixture
[`tests/scf/cases/k_space_scf/bccFe`](../../tests/scf/cases/k_space_scf/bccFe).
It runs no-U, explicit U=0, U=2 eV, U=4 eV, and an eight-step U=2 eV
convergence envelope from fresh copies of the same fixture.

```text
python3 tests/validation/val04_lda_u.py \
  --binary build-rf-debug/bin/rslmto.x \
  --scratch-root /tmp/val04-campaign
```

Measured output from the passing campaign:

| U (eV) | d trace (up, down) | final onsite split (Ry) | d moment | d charge | etot | sumev |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | not retained by inactive-U output path | 0.000000 | 2.199707 | 6.499088 | -2541.981939717 | -1.997839785 |
| 2 | (4.530563, 2.031650) | -0.073217 | 2.472995 | 6.555128 | -2542.202118541 | -2.554184746 |
| 4 | (4.711245, 1.933002) | -0.163708 | 2.635038 | 6.597677 | -2542.379525867 | -3.029020096 |

The response is qualitative and monotonic in this finite-SCF envelope:
increasing U increases the magnitude of the correlated d-level spin
splitting, increases the d moment and charge, and changes both reported
energy observables. Explicit U=0 recovers the no-U run for the checked energy
and d-charge observables; the internal U is zero.

The eight-step U=2 eV residual envelope decreases from `2.300013e-02` to a
minimum of `4.657013e-03`. This establishes improving convergence behavior,
but not a converged material ground state. No fixed-point equality between
separate U calculations is claimed.

## Occupation matrices and common-state RS/KS check

- The algebra unit
  [`tests/unit/test_lda_u_hamiltonian.f90`](../../tests/unit/test_lda_u_hamiltonian.f90)
  checks Hermiticity of the p-shell occupation/potential block, correlated
  shell traces, and the U=J=0 limit.
- The physical campaign checks d-shell traces remain in `[0, 5]` per spin.
  Production Green-function occupation matrices report maximum
  anti-Hermitian residuals of `2.0e-8`–`5.9e-8`; the campaign bounds these at
  `1e-6`, consistent with the energy-integration/serialization tolerance.
- The same unit constructs one frozen onsite potential, inserts it into the
  shared Hamiltonian block, and checks that the reciprocal Fourier assembly at
  a controlled common state returns that onsite block. This is the meaningful
  RS/KS comparison here: the block is shared before solving. It does not
  assert equality between independent RS and KS SCF fixed points.

## Support classification

| Combination | VAL-04 scope | Classification and reason |
|---|---|---|
| Collinear onsite U/J | Fixed-U/J bcc-Fe d shell | **Validated, scoped**. Algebra, U=0 recovery, physical U response, magnetic/charge response, and finite-SCF convergence are covered. |
| RS | Shared onsite block and existing bcc-Fe Hubbard smoke path | **Experimental**. The common block path is checked, but this campaign does not promote a separate converged RS material response. |
| KS / reciprocal | bcc-Fe 12³ finite-mesh reciprocal-SCF campaign | **Validated, scoped** for the recorded finite-step envelope and common-state assembly. |
| Noncollinear onsite U/J | Not promoted | **Development / unvalidated**. The stored LDM and Hubbard potential are two spin-resolved blocks; the correction is added block-diagonally while the noncollinear Hamiltonian has spin-mixing terms. A spinor, rotation-covariant occupation/DC contract is not established. |
| SOC + onsite U/J | Not promoted | **Development / unvalidated**. No spinor occupation and double-counting validation exists for this combination; GBT with SOC is explicitly rejected in the Hamiltonian builder. |
| Hubbard-V | Separate from ordinary onsite U/J | **Development**. Current V is a local diagonal-occupation proxy; full intersite `n^(JI)` is not wired. |

The noncollinear and SOC classifications identify the root cause and scope
boundary; they are not silently converted into broad unsupported claims for
all other Hamiltonian functionality.

## Checklist

- [x] Current onsite formulation documented
- [x] Double-counting convention documented
- [x] Physical U response established
- [x] Occupation matrices checked
- [x] U=0 recovery checked
- [x] Magnetic and charge response checked
- [x] RS/KS common-state comparison made where meaningful
- [x] Noncollinear/SOC support classified with root cause
- [x] Hubbard-V kept separate
- [x] No proxy promoted as full `+V`
- [x] Maturity ledger updated
