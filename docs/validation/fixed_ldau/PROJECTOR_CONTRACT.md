# Fixed-U scalar projector contract (Phase P)

## Scope and status

The supported fixed-parameter Hubbard scope currently remains scalar-relativistic, orthogonal, collinear `nsp=1`. Active fixed-U input in `nsp=2/3/4` remains rejected; no-U `nsp=2` remains available. The protocol's Phase-P material gate is **BLOCKED** until an identified periodic material calculation compares the same frozen physical operator through both production backends and demonstrates the projector/charge relationship with controlled quadrature and broadening refinement.

`tests/validation/fixed_ldau/projector_common_frozen_h.f90` is a finite frozen-operator, source-connected contract probe. It constructs one 18×18 Hermitian matrix (9 orbitals × two collinear spin blocks) and its exact eigenpairs. It supplies the eigenpairs to the public `reciprocal%calculate_ldm_from_projected_dos` wrapper (which calls the internal eigenvector extractor) and supplies the corresponding retarded Lehmann resolvent to `bands%calculate_moments`. It then checks against an independently constructed unitary projector and a finite-window Lorentzian integral. It contains deliberate conjugation and transpose negative controls. This is evidence class **B** for these production occupation-extraction routines on the supplied finite operator, paired with class **A** analytic projector/integration oracles. It bypasses the RS recursion and reciprocal Hamiltonian assembly; it is not class **C** material evidence and does not certify self-consistency, radial QL normalization, or an energy identity. The fixture has one site, so it does not test global site/type mapping or detect multi-site aliasing.

The focused serial Debug CTest passed from the dirty worktree at HEAD `0ef82220d882ab8f5feb5f1e87f7d33284fa8bdf` (GNU Fortran `/opt/homebrew/bin/gfortran`, MPI/OpenMP off). The exact commands were `cmake --build /tmp/rslmto-ldau-remediation-debug --target UnitLdaUProjectorCommonFrozenH rslmto.x -j3` (then rebuilt the test target with `-j2`) and `ctest --test-dir /tmp/rslmto-ldau-remediation-debug -R '^UnitLdaUProjectorCommonFrozenH$' --output-on-failure`; exit 0, 1/1 in 7.22 s. The tested executable SHA256 is `380e08f8a4d850f200594837d878ba55b7034242f34ae6aac21320c2449a899d`. The final command log is `/tmp/rslmto-remediation-phasep-final.log`; CTest detail is `/tmp/rslmto-ldau-remediation-debug/Testing/Temporary/LastTest.log`. The measured KS error against the exact projector was `0`; the actual `bands%calculate_moments` output differed from the independent 6000-channel production-Simpson reconstruction by `1.11e-16`. Maximum elementwise error against the analytic finite-window Lorentz oracle decreased with mesh refinement: 750 channels `3.593160e-2`, 1500 `7.141594e-3`, 3000 `4.743046e-4`, and 6000 `3.267817e-5`. The acceptance budget is `1e-4` absolute occupation per matrix element at 6000 channels (one part in ten thousand of a unit occupation); it is satisfied and more than two orders below the measured sharp-KS versus broadened/windowed-RS difference (`3.668847e-3`). Negative controls produced mismatches `0.380423` (conjugated KS eigenvectors) and `0.379280` (transposed RS resolvent), both above the predeclared `0.1` failure threshold. `max|H0+V_U-H|=2.78e-17`, `max|H-H†|=7.85e-17`, and the maximum imaginary d coherence / V_U entry were `0.190211` / `0.00336391`. These are finite-model routine-consistency results, not material observables.

## Matrix and basis convention

For a local orbital-spin basis index `aσ`, the fixed-U density is

\[
P_{a\sigma,b\tau}=\langle c^\dagger_{b\tau}c_{a\sigma}\rangle.
\]

The scalar collinear implementation stores the two diagonal spin blocks independently as `potential%ldm(l+1,spin,m_a,m_b)`, with `P_ab = ψ_a conjg(ψ_b)`. The Hamiltonian and projector use the same row/column order. Hermitian physical densities satisfy `P_ab = conjg(P_ba)` and real diagonal entries. The off-diagonal orbital entries, including their imaginary parts, are retained for fixed Liechtenstein U/J.

Within the spherical onsite basis, orbital channels are contiguous with offset `l*l`; d indices are `l*l+1 : (l+1)*(l+1)` (one-based indices 5:9 for the spd basis). The documented spd spherical ordering is s `(l=0,m=0)`, p `(1,-1),(1,0),(1,+1)`, d `(2,-2),(2,-1),(2,0),(2,+1),(2,+2)`, as recorded in `source/electrostatics_multipole.f90` and produced by `math.f90::hcpx`. Collinear spin blocks are up then down, separated by `spin_off=norb=9`, as declared in `source/basis.f90`.

This contract assumes the current orthogonal onsite LMTO representation and excludes HOH/generalized-overlap fixed-U calculations. It does not claim that the radial QL channel is an identical metric projector: the relation between that scalar radial charge and the complex onsite P must be demonstrated for a same-state production calculation before Gate P can pass.

## Production extraction and occupations

The KS producer `source/reciprocal_projection.f90::calculate_ldm_from_eigenvectors` accumulates, for each site, k point, spin block and occupied band,

\[
P^{KS}_{ab\sigma}=\sum_{k n}w_k f(\epsilon_{nk},\mu,T)
\psi_{a\sigma,nk}\psi^*_{b\sigma,nk}.
\]

For fixed Liechtenstein mode, this routine uses Fermi-Dirac occupations with `kT=max(T*kB,1e-10 Ry)` and preserves complex entries; the legacy path retains its historical real/transposed convention. K-point weights must sum to one over the full mesh, including distributed meshes.

The RS producer `source/bands.f90::calculate_moments` extracts fixed-U occupations from the retarded Green matrix. Its fixed-Liechtenstein branch forms the matrix spectral density `A_ab(E)=i(G_ab-G^*_{ba})/(2π)` and integrates it to the Fermi level by the production Simpson integrator. For the finite-H probe, the independent analytic comparison is the corresponding Lorentzian integral over the same lower energy bound, broadening `η`, and Fermi endpoint. RS and KS differ by finite energy-window/broadening integration; those errors must be measured and refined rather than folded into a loose tolerance.

The test's deliberately wrong variants must fail: conjugating all KS eigenvectors changes the sign of imaginary orbital coherences, while transposing the spectral matrix by constructing the Green function of `H^T` does the same for RS. These variants test orientation-sensitive complex handling; they do not substitute for real production Hamiltonian or material validation.

## Evidence needed to close Gate P

The finite-H probe is useful for a common-oracle check and regression against the original real/transposed failure mode. Gate P still requires a frozen physical periodic LMTO operator whose exact spectrum is accessible, measured RS energy-window/broadening and KS k-weight/Fermi convergence, site/orbital/spin trace comparison, nonzero imaginary d entries, and serial/MPI/restart agreement. The production LDM trace must then be compared to same-state radial QL with the radial projector definition and any metric factors derived from source equations, not inferred from agreement at six-decimal log precision or from differently mixed SCF iterations.
