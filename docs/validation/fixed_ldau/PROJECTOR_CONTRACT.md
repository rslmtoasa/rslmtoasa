# Fixed-U scalar projector contract (Phase P)

## Scope and status

The supported fixed-parameter Hubbard scope currently remains scalar-relativistic, orthogonal, collinear `nsp=1`. Active fixed-U input in `nsp=2/3/4` remains rejected; no-U `nsp=2` remains available. The protocol's Phase-P material gate is **PARTIALLY EVIDENCED / BLOCKED**: a one-site bcc-Fe production run compares the same frozen periodic operator through RS and KS, and a separate finite two-site probe exercises site-resolved RS/KS extraction. These do not establish a multisite crystal calculation, nontrivial overlap metric, radial-projector relation, converged material occupations, or same-state restart agreement.

`tests/validation/fixed_ldau/projector_common_frozen_h.f90` is a finite frozen-operator, source-connected contract probe. It constructs one 18×18 Hermitian matrix (9 orbitals × two collinear spin blocks) and its exact eigenpairs. It supplies the eigenpairs to the public `reciprocal%calculate_ldm_from_projected_dos` wrapper (which calls the internal eigenvector extractor) and supplies the corresponding retarded Lehmann resolvent to `bands%calculate_moments`. It then checks against an independently constructed unitary projector and a finite-window Lorentzian integral. It contains deliberate conjugation and transpose negative controls. This is evidence class **B** for these production occupation-extraction routines on the supplied finite operator, paired with class **A** analytic projector/integration oracles. It bypasses the RS recursion and reciprocal Hamiltonian assembly; it is not class **C** material evidence and does not certify self-consistency, radial QL normalization, or an energy identity. The fixture has one site, so it does not test global site/type mapping or detect multi-site aliasing.

The focused serial Debug CTest passed from the dirty worktree at HEAD `0ef82220d882ab8f5feb5f1e87f7d33284fa8bdf` (GNU Fortran `/opt/homebrew/bin/gfortran`, MPI/OpenMP off). The exact commands were `cmake --build /tmp/rslmto-ldau-remediation-debug --target UnitLdaUProjectorCommonFrozenH rslmto.x -j3` (then rebuilt the test target with `-j2`) and `ctest --test-dir /tmp/rslmto-ldau-remediation-debug -R '^UnitLdaUProjectorCommonFrozenH$' --output-on-failure`; exit 0, 1/1 in 7.22 s. The tested executable SHA256 is `380e08f8a4d850f200594837d878ba55b7034242f34ae6aac21320c2449a899d`. The final command log is `/tmp/rslmto-remediation-phasep-final.log`; CTest detail is `/tmp/rslmto-ldau-remediation-debug/Testing/Temporary/LastTest.log`. The measured KS error against the exact projector was `0`; the actual `bands%calculate_moments` output differed from the independent 6000-channel production-Simpson reconstruction by `1.11e-16`. Maximum elementwise error against the analytic finite-window Lorentz oracle decreased with mesh refinement: 750 channels `3.593160e-2`, 1500 `7.141594e-3`, 3000 `4.743046e-4`, and 6000 `3.267817e-5`. The acceptance budget is `1e-4` absolute occupation per matrix element at 6000 channels (one part in ten thousand of a unit occupation); it is satisfied and more than two orders below the measured sharp-KS versus broadened/windowed-RS difference (`3.668847e-3`). Negative controls produced mismatches `0.380423` (conjugated KS eigenvectors) and `0.379280` (transposed RS resolvent), both above the predeclared `0.1` failure threshold. `max|H0+V_U-H|=2.78e-17`, `max|H-H†|=7.85e-17`, and the maximum imaginary d coherence / V_U entry were `0.190211` / `0.00336391`. These are finite-model routine-consistency results, not material observables.

### Two-inequivalent-site finite operator probe

[`projector_multisite_frozen_h.f90`](../../../tests/validation/fixed_ldau/projector_multisite_frozen_h.f90) constructs one 36×36 Hermitian operator for two distinct 18-orbital sites, with inequivalent U/J, seed occupations, local d-level spectra, and complex noncubic orbital mixing. The same frozen operator is passed through the public KS projected-DOS wrapper and production `bands%calculate_moments`; each site/spin d block is checked against its independently constructed exact projector or Lorentzian finite-window oracle. The operator includes the two production Hubbard onsite blocks once; it checks `max|(H0+V_U)-H|=2.78e-17` and `max|H-H†|=4.16e-17`. The two local d traces differ, so site aliasing is detectable. Actual KS and RS site outputs are compared against the opposite site's independent oracle as a site-swap negative control.

The source SHA256 is `7425e2308145723052e204d4316968a7ced5bd0eaed9b4ff458458c84abfc98e`. On serial Debug, the focused one-site and two-site CTests passed 2/2; log `/private/tmp/rslmto-remediate-projector-debug-test-final.log`, SHA256 `114e0248ff2afcbb3a8c2e6a569df064ee3630525ac1cdce7371e3d4f1637f96`. On the MPI-linked Debug build, one-site and two-site serial/MPI tests plus selected convergence and parallel-context tests passed 9/9; log `/private/tmp/rslmto-remediate-mpi-targeted-current.log`, SHA256 `530db567e459e3398416ff4ebaa8ffe05458693d45c77aba0d3b364df760226b`.

For the two-site case, KS matches the exact site-resolved projectors to zero at printed precision. RS max elementwise errors against the finite-window Lorentz oracle are `2.259199e-5` (site 1) and `2.300895e-5` (site 2); production RS differs from the independent 6000-channel Simpson reconstruction by `1.24e-16`. Lorentz-oracle quadrature errors decrease monotonically from `5.704945e-2`, `6.516675e-3`, `4.317598e-4` to `2.300895e-5` for 750, 1500, 3000 and 6000 channels. The RS-to-sharp-KS broadening/window difference is `3.050803e-3`. Conjugated-KS and transposed-RS negative controls mismatch by `0.364816` and `0.363804`; the swapped-site oracle comparison mismatches by `0.184084`. Maximum imaginary d coherence is `0.182408`, and the maximum imaginary Hubbard-potential entry is `0.00263048`.

This is finite class-B evidence for production site-resolved occupation extraction and distinct-site mapping in an orthogonal 36×36 model. It is not a real-space crystal/periodic-operator comparison, an overlap-metric test (the metric is identity), or a material result. The production LDM/potential API is keyed by symbolic-atom type: it stores one occupation and constructs one potential per type. Independent inequivalent occupations on multiple sites sharing the same type are not represented; this probe uses separate types.

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

## Production periodic operator comparison

[`projector_periodic_operator.py`](../../../tests/validation/fixed_ldau/projector_periodic_operator.py)
reconstructs every captured production KS `H(k)` from the full real-space
neighbor EE blocks captured during an RS run, using the captured fractional
bond vectors and the production `exp(+i 2π k·R)` convention. It also requires
the RS and KS captures to agree elementwise on every real-space block and bond
map, and compares their local orbital/site basis maps. The one-site bcc-Fe
captures at `/tmp/rslmto-energy-trace/common_rs_09` and `common_ks_09` pass:
15 complete 18×18 blocks (4860 matrix elements), 15 unique bonds (45 rows,
three fractional-vector components per bond), 1728 KS k points, exact RS/KS
real-space agreement, and maximum reconstructed-H(k) error `3.330669e-16`
against the `2e-11` threshold. The result record is
[`periodic_operator_bccfe_common_09.json`](periodic_operator_bccfe_common_09.json).

The paired pre-run decks have the same Fe starting-state SHA256
`fd5d907628cc4ee41ac3b1147aeb1f9e8ec1fe8bb2f38069afe1cf63384ae233`; their
input files differ only in `use_kspace=.false.` versus `.true.` (RS SHA256
`7b14254bfa60ed4e4d42a9fbc5c2f01ae9dcf6ba3081f2e43c707dbb2b8c2f41`, KS
`78ccdc39a93c1b3b45a1a50cc5fc6ca188fe5f0af62196b4dbe184308a50056d`). Both
use `nsp=1`, `hoh=.false.`, input U=1.8 eV and J=0 (0.9× a 2.0 eV reference
value). The CSV `lambda=1` is a trace-run annotation; it is not that 0.9
physical U scale or a finite-difference slope. The executable SHA256 is
`ffa895abc32af3cce2123c06efc703ad01297bfb22c3ea474da99201cd9c6764`.

As a separate output-density check, the captured KS eigenvectors, occupations,
and actual k weights reproduce the captured KS d-shell `P_OUT` to
`5.56e-15`; the 1728 weights sum to one and the weighted occupation is 8.
The RS and KS d-shell `P_OUT` matrices are not equal: their maximum
elementwise difference is `3.419410e-3` and RMS difference is `8.137314e-4`.
Their d traces are RS `(4.506326, 2.049989)` and KS `(4.501409, 2.053512)`.
The captures use different integration occupations (RS trace manifest reports
zero temperature; KS uses `kT=0.00190010442 Ry`), and no broadening/energy-grid
refinement envelope was run for this pair. Therefore the P comparison is an
observed finite-run difference, not a proven error bound or a claim of
backend-equivalent occupations.

The checker passes a negative Fourier-phase control (changing `+R` to `−R`
produces a `1.942e-1` H(k) mismatch) and rejects a corrupted RS basis-map site
label. This is a real periodic-operator comparison, class B evidence; the
one-site cell cannot detect permutations or aliasing between inequivalent
sites. The raw CSVs are retained under `/tmp` and are deliberately not copied
into the repository (112 MB RS and 240 MB KS); the JSON record preserves their
hashes and the runnable command. Re-run with:

```text
python3 tests/validation/fixed_ldau/projector_periodic_operator.py \
  /tmp/rslmto-energy-trace/common_ks_09/records.csv \
  --run-id run --iteration 1 --lambda 1.0 \
  --rs-trace /tmp/rslmto-energy-trace/common_rs_09/records.csv \
  --rs-run-id run --rs-iteration 1 --rs-lambda 1.0
```

This closes only the same-frozen-periodic-operator RS/KS comparison subset of
Gate P. It is not a multi-site test, a controlled RS/KS occupation-envelope
test, an RS/KS restart/MPI test, a radial QL normalization proof, a converged
material result, or class C physical validation.

## Evidence needed to close Gate P

The finite-H probe is useful for a common-oracle check and regression against the original real/transposed failure mode. Gate P still requires a frozen physical periodic LMTO operator whose exact spectrum is accessible, measured RS energy-window/broadening and KS k-weight/Fermi convergence, site/orbital/spin trace comparison, nonzero imaginary d entries, and serial/MPI/restart agreement. The production LDM trace must then be compared to same-state radial QL with the radial projector definition and any metric factors derived from source equations, not inferred from agreement at six-decimal log precision or from differently mixed SCF iterations.
