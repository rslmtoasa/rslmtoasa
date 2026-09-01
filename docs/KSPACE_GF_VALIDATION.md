# K-space Lehmann Green-function validation

## Scope and status

TDDFT-05 validates the newer K-space Lehmann one-electron propagator in the
orthogonal `reciprocal_mode='ham_only'` representation.  The validation is
independent of TD-DFT response construction and uses deterministic finite
Hermitian matrices.  It does not claim generalized-overlap support or
material-level convergence.

The reusable pure kernel is
`lehmann_kernel_mod::lehmann_kspace_resolvent`:

```text
G(k,z) = sum_n |psi_n(k)><psi_n(k)| / (z - epsilon_n(k))
```

The existing response-side eigenpair Green provider now calls this same
primitive.  `lehmann_pair_block` retains the production inverse-Bloch phase
and batched-energy path for site/intersite blocks.

## Direct-inverse comparison

`tests/unit/test_kspace_gf_validation.f90` builds five deterministic 4-by-4
Hermitian `H(k)` matrices with complex off-diagonal phases, diagonalizes each
with LAPACK `zheev`, and compares the Lehmann result against the independent
`dyson_kspace_inverse` oracle with `Sigma=0`:

```text
G_direct(k,z) = [z I - H(k)]^-1
```

Five complex energies include both upper- and lower-half-plane values.  The
maximum elementwise difference over all 25 `(k,z)` pairs is:

```text
1.0549e-14
```

This is the one-electron resolvent test; it does not compare two TD-DFT
response implementations.

## Identities, spin limits, and DOS

The same fixture checks the following:

| Check | Result |
| --- | ---: |
| `G^A(E) = G^R(E)†` | `1.3878e-16` max error |
| `A=i(G^R-G^A)` Hermiticity | `2.4980e-16` max error |
| minimum diagonal of `A` | `5.4488e-01` |
| `zG -> I`, `z=i 1e10` | `8.3500e-11` max error |
| total DOS vs eigenvalue Lorentzian path | `2.2204e-15` max error |
| projected local DOS vs eigenvector-weighted path | `1.3323e-15` max error |
| integrated four-state spectral weight | `2.0373e-03` absolute error |

The DOS checks use the existing `bsf_spectral_trace` convention,
`A=-Im Tr G/pi`, and compare both total and spin-up projected traces with
the corresponding eigenvalue/eigenvector spectral weights.  The normalization
integration uses `[-100,100] Ry`, `eta=0.08 Ry`, and 10,001 trapezoid points;
the remaining error is the analytically expected finite Lorentzian tail.

Spin ordering is `(orbital-up..., orbital-down...)`.  A duplicate up/down
spectrum gives zero x/y/z Pauli components, while a split collinear spectrum
gives zero x/y and a nonzero z component.  Thus the tested resolvent uses the
same collinear spin-block convention as the current TD-DFT response vertices.

## Generalized-overlap boundary

The validated kernel assumes Euclidean-orthonormal eigenvectors and therefore
only represents `S=I`.  For a generalized eigenproblem,

```text
H C = S C epsilon,       C† S C = I,
G(k,z) = [z S(k) - H(k)]^-1
```

the metric-consistent spectral representation and observable projectors must
be derived together.  In the usual complete generalized-eigenvector
normalization the matrix inverse can be written with `C(z-epsilon)^-1 C†`,
but the response weights and local projections still carry the metric
contract.  The current Lehmann/Dyson Green path does not implement that
contract.

Consequently, both the low-level `reciprocal%fill_green` entry point and the
TD-DFT susceptibility driver reject every
`reciprocal_mode` other than `ham_only` with an explicit diagnostic.  The
diagnostic states that `G=(z*S-H)^-1` and metric-consistent spectral weights
are not implemented.  No generalized-overlap case is approximated by
`S=I`, and no empirical correction is applied.

## Batched-energy profile

The production `lehmann_pair_block` call was profiled on the same executable
with `nmat=12`, `nk=48`, and a caller-owned output buffer.  This is an
informational GNU Debug CPU profile; it has no CI timing threshold.

| complex energies | elapsed seconds | output bytes | automatic scratch bytes |
| ---: | ---: | ---: | ---: |
| 1 | `5.7900e-04` | 2,304 | 2,304 |
| 8 | `1.7200e-03` | 18,432 | 2,304 |
| 32 | `5.5780e-03` | 73,728 | 2,304 |
| 128 | `2.0841e-02` | 294,912 | 2,304 |

The kernel has no allocatable work inside the k/band/energy loops.  The
output is supplied by the caller and scales as
`16*nmat*nmat*nenergy` bytes; the automatic complex scratch block is
`16*nmat*nmat` bytes.  The source-contract test also checks this no-inner-loop
allocation property.

## Reproduction

From a configured build with standalone unit tests enabled:

```text
cmake --build build --target UnitKspaceGFValidation UnitTddftGreenChiKS -j2
ctest --test-dir build --output-on-failure \
  -R '^(UnitKspaceGFValidation|UnitKspaceGFValidationSource|UnitLehmannChain|UnitGammaSupercell|UnitTddftGreenChiKS|UnitTddftBackend|UnitTddftDispatch)$'
```

Observed result on the current development build: **7/7 tests passed**.
The focused regression set retains the pre-existing chain, Gamma-supercell,
response-bubble, backend-interface, dispatch, and new source-contract tests.

## Checklist

- [x] Lehmann and direct-inverse `G(k,z)` agree.
- [x] Retarded/advanced identities pass.
- [x] Spectral normalization and large-`|z|` asymptotics pass.
- [x] Spin convention matches the collinear TD-DFT convention.
- [x] Generalized-overlap cases are rejected with a precise explanation.
- [x] A basic batched-energy performance/allocation profile exists.
