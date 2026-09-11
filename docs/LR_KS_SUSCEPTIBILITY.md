# LR-06 collinear transverse KS susceptibility

## Status and claim level

The first bare-response backend is implemented in
[`source/lr_ks_susceptibility.f90`](../source/lr_ks_susceptibility.f90) and is
registered in the `rslmto` library.  It is the explicit spectral/Lehmann route
in the certified LR-03 Pauli/no-SOC response space.

This establishes the LR-06 claim:

> the bare spectral KS susceptibility is implemented and passes its algebraic
> and numerical oracles.

It does not establish an XC kernel, Dyson enhancement, Goldstone correctness,
magnon frequencies, converged Fe/Ni validation, or literature agreement.

The implementation consumes the authoritative contracts in
[`LR_TDDFT_CONVENTIONS.md`](LR_TDDFT_CONVENTIONS.md),
[`LR_RESPONSE_SPACE_ALGEBRA.md`](LR_RESPONSE_SPACE_ALGEBRA.md), and
[`LR_PAULI_TRANSITION_VERTEX.md`](LR_PAULI_TRANSITION_VERTEX.md).  No legacy
TD-DFT susceptibility code or convention is used.

## Supported baseline

The public evaluator fails closed outside the initial production tuple:

```text
reciprocal mode       = ham_only
Hamiltonian order     = second
overlap               = orthogonal
spin                  = collinear, two-channel, no SOC
extra response op.    = absent
angular basis         = sp or spd
response mesh         = complete direct LR-04 radial/angular space
```

Only one certified circular sector is evaluated per request: `chi_plus` or
`chi_minus`.  Allocating a nominal four-component response is not treated as
support for an unimplemented sector.

## Spectral equation and circular convention

For a canonical pointwise response coordinate (I=(a,L,M,i,mu)), the raw
retarded kernel is

\[
 \chi^{KS,\mu\nu}_{IJ}(\mathbf q,\omega)
 =\frac{1}{N_k}\sum_{\mathbf k,nm}
 \frac{f_{n\mathbf k}-f_{m,\mathbf k+\mathbf q}}
 {\omega+\epsilon_{n\mathbf k}-\epsilon_{m,\mathbf k+\mathbf q}+i\eta}
 T^\mu_{nm;I}
 (T^\nu_{nm;J})^* .
\]

The code uses the LR-03 specialization for the halved circular Pauli
operators.  For `chi_plus`, the transition vertex is built with

\[
 \sigma^+=(\sigma_x+i\sigma_y)/2,
\]

and for `chi_minus` with (sigma^-).  The retarded source is the adjoint
circular operator, so the same channel vertex occurs in the conjugated outer
product.  The published circular quantity therefore contributes the exact
factor 2:

\[
 \chi^{KS,+}_{IJ}=2\frac{1}{N_k}\sum_{\mathbf k,nm}
 \frac{f_n-f_m}{\omega+\epsilon_n-\epsilon_m+i\eta}
 T^+_{nm;I}(T^+_{nm;J})^*,
\]

with the analogous equation for `chi_minus`.  The measurement operator remains
first and the source operator second.  No sign, factor, endpoint order, or
frequency broadening is inferred from Fe behavior.

The explicit `eta` is required to be positive and is always inserted as
`+i*eta`.  Zero frequency is not a separate approximation: `omega=0` is passed
through the same denominator and accumulation path.

The k-point normalization is the LR-GF-01/reciprocal convention: raw
nonnegative k weights are divided by their total weight, which is the finite
mesh realization of (1/N_k).  Every band in the finite eigensystem is used;
there is no default energy or band window.

## Response-space representation

The accumulator first forms the sampled raw kernel and then calls the LR-04
conversion exactly once per frequency:

\[
 B_{IJ}=\chi^{KS,\mathrm{raw}}_{IJ}W_J.
\]

The result stores (B), the canonical right-weighted operator.  Consequently
`response_apply_operator` performs the physical source contraction with ordinary
matrix multiplication, and no consumer may add another radial weight.  The
origin column is zero through the LR-04 conversion because its physical metric
weight is exactly zero.

The radial/angular response is never projected to site-only data.  LR-05 emits
the complete pointwise transition vector, including all configured product
harmonics and radial points; LR-06 forms its full outer product.

## Electronic-state provenance and q handling

`lr_electronic_state` is an immutable snapshot contract for the evaluator.  Its
initializer copies:

- the complete finite eigenvalue/eigenvector arrays;
- folded fractional k-points and their raw nonnegative weights;
- explicit Fermi occupations for every band and k-point;
- EF, electronic temperature, and energy-zero metadata;
- reciprocal mode, Hamiltonian order, overlap, spin, SOC, and extra-operator
  provenance.

The evaluator accepts both states as `intent(in)` and never solves EF, derives a
new occupation array, rephases an eigenvector, or truncates bands.  The
`lr_snapshot_from_reciprocal` adapter copies an accepted canonical reciprocal
state.  It computes the explicit Fermi occupations from the already stored EF
and temperature only; it never recomputes EF.  The
`lr_q_endpoint_from_reciprocal` adapter calls the existing arbitrary-k
eigenpair service for the exact list `k+q`, stores its folded endpoints, and
passes the returned eigenvectors through unchanged.

Before accumulation, the evaluator verifies that each endpoint is the exact
half-open-BZ fold of `k+q`, that the state dimensions and EF/temperature/energy
zero agree, and that both snapshots carry the certified baseline provenance.
There is no nearest-mesh substitution or hand-written BZ-crossing phase.

## Public API

The preferred request/result form is:

```fortran
type(lr_ks_susceptibility_request) :: request
type(lr_ks_susceptibility_result) :: result

request%q = [qx, qy, qz]
request%frequencies = [omega1, omega2]
request%eta = eta
request%channel = lr_channel_plus       ! or lr_channel_minus
request%response_space => space
request%radial_bases => radial_bases
request%electronic_state => left_state
request%q_endpoint_state => k_plus_q_state
call evaluate_lr_ks_susceptibility(request, result)
```

`result%susceptibility(:, :, ifrequency)` is the LR-04 canonical complex
matrix.  It retains `q`, frequencies, eta, channel, and response-space
metadata, and has no plotting or mode-extraction dependency.  An explicit
argument overload is available for callers that do not want to store the
object references in the request.

The static identity diagnostic is deliberately separate:

```fortran
call evaluate_lr_static_residual(space, chi_static, field, magnetization, abs_residual, rel_residual)
```

It applies the supplied canonical (chi^{KS}(q=0,\omega=0)) to independently
stored field and magnetization vectors with the LR-04 metric and reports the
raw residual.  It does not rescale the susceptibility, shift an eigenvalue, or
enforce Goldstone behavior.

## Tests and evidence

`UnitLrKsSusceptibility` is a standalone Fortran test registered with CTest.
Its independent levels are:

- a two-level analytic finite-model matrix element with the LR-03 circular
  factor and LR-04 source-column metric;
- a deliberately straightforward nested (k,n,m,\omega,I,J) loop on a
  deterministic two-site production-style fixture, independent of both the
  production susceptibility and production transition evaluator;
- q=0 support and an arbitrary off-mesh q with a folded second endpoint;
- combined q/-q retarded/advanced frequency covariance with circular-sector exchange;
- global spin reversal, verifying `chi_plus`/`chi_minus` exchange;
- simultaneous unitary rotations of complete degenerate occupied and empty
  subspaces;
- a raw q=0 static residual using separately filled field and magnetization
  vectors.

The focused verification was:

```text
cmake --build build --target UnitLrKsSusceptibility -j2
ctest --test-dir build --output-on-failure -R '^UnitLrKsSusceptibility$'
```

The executable reported:

```text
explicit-loop production oracle error =   6.2946E-10
Analytic two-level finite-model error =   6.5362E-09
q=0 and off-mesh/folded k+q support: PASS
q/-q retarded/advanced covariance error =   3.5394E-10
retarded/advanced circular covariance error =   3.5394E-10
global spin-reversal circular exchange error =   0.0000E+00
degenerate-subspace invariance error =   2.0720E-25
Raw q=0 static residual (absolute, relative) =   7.94988083E-04   1.00000000E+00
Raw static residual: diagnostic only; susceptibility unchanged
UnitLrKsSusceptibility: PASS (Lehmann, explicit loop, covariance, degeneracy)
```

These are algebraic and independent numerical-oracle claims.  The nonzero raw
static residual is retained as evidence and is not repaired.

### Fe fixture diagnostic envelope

The real converged bcc-Fe radial fixture used by LR-01/LR-02N is
`tests/scf/cases/bulk/bccFe`.  The LR-06 diagnostic envelope recorded for a
future production sweep is:

| item | diagnostic setting |
|---|---|
| k mesh | full (8\times8\times8) Monkhorst--Pack mesh, inherited from LR-02N |
| eta | explicit `0.02 Ry` response request; not inherited from DOS/GF settings |
| radial representation | complete direct LR-04 logarithmic radial mesh, pointwise (T(r_i)), canonical (\chi_{raw}W) output |
| angular cutoff | `spd`, complete product cutoff (L_{max}=4) |

This is a provenance/convergence diagnostic envelope, not a material
susceptibility result or Fe validation.  A real Fe response sweep remains part
of TDVAL-01 and must report systematic k, eta, radial, and angular changes.

## Limitations and explicit exclusions

- The only electronic-state backend implemented here is spectral/Lehmann.
  Reciprocal-GF bubbles remain an independent LR-GF-02 task.
- The baseline is collinear, no SOC, orthogonal, second-order/HOH, `ham_only`,
  and `sp`/`spd`.
- The Pauli/no-SOC transition vertex inherits the controlled SR-to-Pauli
  approximation documented by LR-02N/LR-03.
- No Kxc, ALSDA interaction, Dyson enhancement, Goldstone correction, mode
  fitting, site projection, or hidden band truncation is present.
- Individual degenerate eigenvectors are not physical observables.  The
  evaluator does not gauge-fix or rotate them; invariance is required only for
  the complete summed susceptibility.
- The standalone unit fixture is not a converged Fe/Ni calculation and does
  not establish literature agreement.
