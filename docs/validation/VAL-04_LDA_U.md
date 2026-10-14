# VAL-04 — fixed-parameter onsite LDA+U validation status

## Status

**PARTIAL / BLOCKED.** The repository contains a complex orbital occupation
matrix, fixed-parameter Liechtenstein/FLL potential and energy evaluator, plus
algebraic and short runtime checks. This evidence does not certify a converged
material calculation or a physically correct production `ETOT`. The production
energy accounting through spectral moments, `PL`, `NEWRHO`, `SUMEV`, `SUMEC`
and `ETOT` remains unresolved (Gate E); the active-U SCF stopping criterion
also does not yet require full complex occupation and potential residuals
(Gate C). Finite-iteration response data must not be described as a converged
ground state.

The supported fixed-U/J input scope remains deliberately narrow: collinear
scalar-relativistic `nsp=1`, an onsite d channel, exact `control%lmax=2`
(`spd` basis), and the validated orthogonal Hamiltonian path. Fixed-U/J with
`nsp=2` SOC, `nsp=3/4` noncollinear modes, `lmax=3` (`spdf`), or `hoh=.true.`
is rejected by the input guard. `spdf` is excluded because its 16-orbital
`hcpx` path currently uses the identity rather than a cubic-to-spherical
transformation, while the d-shell Coulomb tensor is in complex spherical
`m=-2..2` order. The fixed-U projector remains in basis order, so accepting
`spdf` would contract a cubic-harmonic occupation matrix with a spherical
tensor. This guard applies only when fixed Liechtenstein U/J is active; absent
U/J and explicit U=J=0 remain inactive. No-U calculations retain their
existing mode behavior; input tests confirm no-U `nsp=2` and `lmax=3` remain
available. A full
orbital-spinor correlated density and a validated radial/XC contract are
required before fixed-U/J SOC can be supported.

## Implementation contract

The fixed-parameter Liechtenstein route stores complex spin-resolved
`potential%ldm(l,spin,m_a,m_b)` matrices in the standard orientation
`P_ab = psi_a psi_b*`. The corresponding complex Hermitian Hubbard potential
is inserted into the onsite Hamiltonian block. The existing energy evaluator
computes interaction, FLL double-counting, Hubbard correction and `Tr(VP)`
from explicit density/potential snapshots; that evaluator is not evidence that
the production radial total energy includes the required contribution.

U and J inputs are in eV and converted once to Ry. Supplying U without J means
J=0. Explicit U=J=0 and absent U/J leave the correction inactive. These input
contracts are covered by
[`tests/validation/val04_lda_u_inputs.py`](../../tests/validation/val04_lda_u_inputs.py).

Hubbard-V, SC-U and ACBN0 are separate code paths and are not promoted by this
document. Do not infer their support or behavior from the fixed-parameter
Liechtenstein checks.

## Historical campaign

The original [`tests/validation/val04_lda_u.py`](../../tests/validation/val04_lda_u.py)
forces `nsp=2` and runs a small finite-iteration U sweep. That campaign
predates the current fixed-U mode guard. Its U>0 inputs are now expected to be
rejected as SOC+U, so its previously recorded numerical values are historical
observations from an unsupported route, not validation of the current code.
The script and original input fixtures are retained for provenance; do not use
its output as a current passing test or as an energy reference.

Earlier nsp=1 finite-step smokes, tensor references, and unit tests remain
useful evidence for their specific layers, but they do not establish
production-SCF convergence or close the energy/convergence blockers. The
current acceptance ledger is
[`fixed_ldau/REMEDIATION_GATES.md`](fixed_ldau/REMEDIATION_GATES.md); see also
[`fixed_ldau/PHASE3.md`](fixed_ldau/PHASE3.md) for the unresolved energy
identity and [`fixed_ldau/PROJECTOR_CONTRACT.md`](fixed_ldau/PROJECTOR_CONTRACT.md)
for the production projector comparison status.

## Required evidence before promotion

- The finite-H source-connected occupation-extraction probe and its measured
  Simpson/broadening error are documented in
  [`fixed_ldau/PROJECTOR_CONTRACT.md`](fixed_ldau/PROJECTOR_CONTRACT.md). Gate P
  still requires a periodic physical LMTO operator through the real RS and KS
  assembly paths, site-map and radial-Ql relation checks, and MPI/restart
  agreement (Gate P / R02, R14).
- An analytic and numerically verified production energy identity connecting
  the Hubbard functional to the actual radial `ETOT` in RS and KS (Gate E / R04).
- Active-U convergence criteria that include full complex occupation and
  potential residuals, plus restart reproducibility (Gate C / R05, R14).
- Converged scalar bcc-Fe and noncubic J≠0 material campaigns with input,
  executable, compiler, convergence and output provenance (Gate M1 / R06–R07).

Only label scalar fixed-U/J **validated** after the scalar gates required by
the remediation protocol pass. SOC+U remains unsupported until the complete
10×10 spinor functional and production path, radial/XC scope, and converged
SOC material gates pass.
