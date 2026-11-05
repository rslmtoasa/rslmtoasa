# GBT oracles

Developer-owned reference data for the k-space GBT tests. Read-only to AI sessions (CLAUDE.md, "Verification integrity").
Changing any file here is a developer decision, recorded in the commit message.

| File | Content |
|---|---|
| `make_gbt_references.py` | Generator (numpy only) with built-in self-checks; writes the two `.nml` files |
| `o1_one_band.nml` | O1: one-orbital spiral model, eigenvalues (read by `UnitGbtKspaceOneBand`) |
| `o2_force_theorem.nml` | O2: force theorem on the O1 model, three sets (read by `UnitGbtKspaceForceTheorem`) |
| `tolerances.nml` | Every tolerance the GBT tests use, with its basis |

Regenerate: `python3 tests/gbt/oracles/make_gbt_references.py` from the repository root. Every self-check must print `ok`.

## Model
One s-like orbital, simple cubic, a = 1, nearest-neighbour hopping t = -0.05 Ry, W_up = W_dn = 1,
local-frame onsite C = -Delta/2 (majority) and +Delta/2, Delta = 0.15 Ry. Conventions as in the GBT restart spec:
spin-up (majority) at k - q/2; D(alpha) = diag(e^{-i alpha/2}, e^{+i alpha/2}); U = D(phi) Ry(theta) D(-phi).
k and q are given in direct coordinates of each cell, with Cartesian values (2 pi / a) alongside.

## O1
Closed-form cone eigenvalues for one site; a doubled cell (sc doubled along x), parallel and antiparallel
(theta_2 = pi - theta, phi_2 = phi + pi, so m_2 = -m_1 at every theta).
Self-checks, independent of the GBT algebra: explicit lab-frame supercells (spin-independent hopping,
onsite -(Delta/2) sigma.m_i) for q = 1/2, 1/3, 1/4 along x against the closed form and the folded two-site sets;
the doubled cell at the same Cartesian q against the two folded one-site spectra. All at round-off (<= 3.3e-16 Ry).

## O2
Gamma-centred 24^3 mesh, k_B T = 0.002 Ry; mu solved once per set for (q, theta) = (0, 0).
Sets: 1, N = 0.6; 2, N = 0.9; 3, set 1 with both C values shifted by +0.5 Ry.
Written per set: Omega(q, theta) (Mermin grand potential) at fixed mu, N(q, theta), and J(0) - J(q).
J(q): model form of the J(q) expression accompanying Eq. (7) of Pajda et al., PRB 64, 174402 (2001),
with P_up - P_dn -> Delta. The paper places spin-up at k+q; the two placements are equal here because
eps(k) = eps(-k). The finite-T Fermi-function form is our extension; the paper uses a T = 0 contour.
Self-checks: Richardson (0.25, 0.5, 1 deg) DeltaOmega/sin^2(theta) equals J(0) - J(q) to <= 2.3e-8 of max|J(0)-J(q)|;
Omega(q, 0) = Omega(0, 0) to round-off; set 3 reproduces set 1's Omega, N and J(0) - J(q).
Omega(0, theta) = Omega(0, 0) holds by construction (global rotation) and is not counted.
Mutation (band energy at fixed mu without -mu N) violates the relation by 18, 0.66 and 100 times max|J(0)-J(q)|
in sets 1, 2 and 3.

