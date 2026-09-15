# TDDY-01 enhanced susceptibility and loss matrix

## Status and claim level

The LR-04 canonical Dyson and loss layer is implemented in
[`source/tddft_dyson.f90`](../source/tddft_dyson.f90) and registered as
`UnitTddftDyson`.  The focused evidence establishes the discrete algebra,
independent finite-matrix solves, metric placement, retarded loss sign,
near-pole diagnostics, covariance metadata, and interaction-route separation.
It does not claim a converged Fe/Ni response, a magnon branch, stiffness, or
literature agreement.

The layer consumes a completed LR-06 `chiKS` and a canonical interaction from
one explicitly selected KXC-01, GSR-01, or GCR-01 route.  It does not select or
reconstruct an interaction route, and it does not use the retired pair-`Xi`
enhancement path.

## TDVK-07 compact representation contract

TDVK-07 certifies the same Dyson service in the accepted orthonormal LMTO
product representation.  If `U` contains the retained weighted-orthonormal
product modes and `W` is the point radial metric, the compact coordinate maps
are

\[
 c=U^H W^{1/2}x,\qquad x=W^{-1/2}Uc,
 \]

and a local point operator is supplied as

\[
 K_c=U^H K_{point}U.
\]

For `request%compact_orthonormal=.true.`, the service therefore uses ordinary
compact multiplication with no point-space metric insertion:

\[
 D_c=I-\chi^{KS}_cK_c,\qquad D_c\chi_c=\chi^{KS}_c.
\]

The solve remains LAPACK `zgesv` with the complete `chiKS` matrix as the
right-hand side; no inverse, denominator shift, pseudoinverse, or additional
eta is used.  Compact loss is the ordinary orthonormal form
`L=-(chi-chi^H)/(2*i*pi)`.  The legacy point-space path and its LR-04 metric
adjoint remain unchanged.

The result now records the minimum-magnitude denominator eigenvalue and the
Frobenius, relative-Frobenius, and infinity-norm Dyson residuals.  TDVK-07's
independent point-space projection oracle uses a separate pivoted Gaussian
solve and separate explicit projection; it is not a second call to the Dyson
helper.

## Discrete LR-04 Dyson equation

LR-04 stores a nonlocal raw pointwise operator `A` as the canonical
right-weighted matrix `B=A*W`.  A local radial interaction is supplied in the
canonical local form defined by LR-04.  Both are operators on pointwise
response vectors, so canonical composition is ordinary matrix multiplication:

\[
 \Xi_{\rm can}=\chi_{KS,\rm can}K_{\rm can},
 \qquad
 D=I-\Xi_{\rm can}.
\]

For every frequency the production equation is

\[
 \boxed{D(\mathbf q,\omega)\,\chi(\mathbf q,\omega)
        =\chi_{KS}(\mathbf q,\omega).}
\]

The right-hand side contains all LR-04 source-column quadrature already
present in `chiKS`.  The interaction contains its own canonical source-column
measure when it is nonlocal.  No `W` is inserted in the Dyson module and no
raw/canonical conversion is performed implicitly.  The solve uses LAPACK
`zgesv` with all columns of `chiKS` as simultaneous right-hand sides; an
explicit inverse is never formed.

The request/result API is:

```fortran
type(tddft_dyson_request) :: request
type(tddft_dyson_result) :: result

request%response_space => space
request%q = q
request%frequencies = omega
request%eta = eta
request%channel = 'chi_plus'
request%ks_susceptibility = chiKS_canonical
request%canonical_interaction = kxc_canonical
request%interaction_route = lr_dyson_route_direct_alsda
request%interaction_provenance = 'KXC-01 direct ALSDA LR-03'
request%electronic_state_provenance = 'LR-06 ham_only complete eigenpair snapshot'
call evaluate_tddft_dyson(request, result)
```

The route must be one of:

| route | required provenance prefix | supplied object |
|---|---|---|
| `direct_alsda` | `KXC-01` | direct canonical ALSDA interaction |
| `goldstone_sumrule` | `GSR-01` | independent canonical sum-rule interaction |
| `direct_alsda_goldstone_corrected` | `GCR-01` | explicitly selected corrected canonical interaction |

This is metadata validation, not a claim that the interaction was generated
by this layer.  The caller must only pass a route after its own prerequisite
route has passed.  The selected route, interaction provenance, electronic-state
provenance, q, frequency grid, eta, response representation, and response-space
metadata are copied to every result object.

## Conditioning and pole policy

The denominator singular values are evaluated independently with SVD and the
result retains the minimum/maximum singular values and 2-norm condition number
for every frequency.  The result also retains the LAPACK solve `info` code and
per-frequency status.

The diagnostics separate three cases:

* a solved denominator with condition number at least `1e8` is reported as a
  near-singular collective-pole candidate;
* a solved denominator with condition number at least `1e12` is additionally
  reported as numerically ill-conditioned; and
* a precision-limited or failed solve is reported as numerically singular and
  `BLOCKED`.

The thresholds only classify the result.  No pole is shifted, clipped,
regularized, or replaced by a pseudoinverse.  A finite-eta near pole remains
available in the result so later validation can decide whether it is physical.

## Retarded loss matrix

The LR-03 loss convention is

\[
 L=-\frac{\chi-\chi^\dagger}{2i\pi}.
\]

For a canonical LR-04 matrix the dagger in this equation is represented by the
LR-04 metric adjoint:

\[
 L_{\rm can}=-\frac{\chi_{\rm can}-\chi_{\rm can}^{\ddagger_W}}{2i\pi}.
\]

This is the canonical form of the physical/raw anti-Hermitian part.  It is
metric-Hermitian, and for a scalar diagonal response it is
`-Im(chi)/pi`; positive-frequency absorption therefore has positive loss in a
channel where `Im(chi)<0`.  Off-diagonal entries use the paired adjoint entry,
not an elementwise imaginary part.  The circular channel is never negated to
make a plot positive, and `chi_plus` and `chi_minus` retain their LR-03
frequency conventions independently.

The exact null-measure origin is given a zero active-space loss extension.  It
does not change a physical contraction and avoids inventing an inverse origin
weight.  The ordinary `tddft_loss_matrix(chi)` overload is retained only for
isolated orthonormal algebra fixtures; production LR-04 callers use
`tddft_loss_matrix(space, chi)`.

## Test evidence

`UnitTddftDyson` covers:

1. the scalar closed-form Dyson expression;
2. a noncommuting 2x2 `chiKS`/`K` fixture with an independently evaluated
   solve;
3. a metric-sensitive continuum-discretized quadrature fixture, including a
   deliberate weight-free wrong result;
4. an exact static q=0 Goldstone-style near pole at controlled finite eta;
5. the LR-03 q/frequency/circular covariance and global spin-reversal
   covariance identities;
6. ordinary and metric loss Hermiticity/sign/origin checks; and
7. separate direct ALSDA, Goldstone sum-rule, and corrected-ALSDA route objects
   with independent metadata.

Focused verification:

```text
cmake --build build --target UnitTddftDyson -j2
ctest --test-dir build --output-on-failure -R '^UnitTddftDyson$'
```

Result: `UnitTddftDyson` passed (1/1 test).

`UnitTddftCompactDysonOracle` is the TDVK-07 compact gate.  Its two-site,
multi-sector complex/noncommuting fixture compares an independently solved
point-space response projected into compact space with the live compact Dyson
solve.  It reports `dF=3.7740745696804097e-15`, relative `dF=6.2590899944262202e-15`,
`dInf=2.9736722724136411e-16`, loss `dF=3.2737677930040462e-16`, and Dyson
residual relative `2.7566246189202793e-16`.

## Scope boundary

This task implements enhanced susceptibility and loss products only.  It does
not fit Lorentzians, track peaks, assign magnon branches, compute stiffness,
diagonalize a mode matrix, add a pair-`Xi` route, or claim material/literature
validation.  Those are downstream validation or post-processing tasks.
