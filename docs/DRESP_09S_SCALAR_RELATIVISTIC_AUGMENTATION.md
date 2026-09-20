# DRESP-09S scalar-relativistic L=0 augmentation

Status: **PASS-B**  
Primary classification: **`BASIS_ROTATION_TERM_REQUIRED`**

This milestone recovers and tests the scalar-relativistic transverse
augmentation contract for a global rigid rotation of accepted bcc Fe.  It is
restricted to the accepted 4x4x4, 64-k-point, collinear, orthogonal HOH state
and an L=0 field.  It does not certify a general L>0 SR vertex, ALSDA, a
Dyson solve, or BES/Halle.

## Legacy representation

The authoritative implementation is `RSEQSR`/`RSQSR1`/`RSQSR2`/`GINTSR` in
`source/self.f90`.  The first-order variables used by the legacy solver map
to the stored arrays as

\[
 G_l = G(:,1)=r g_l,\qquad
 \Phi_l = G(:,2)=\frac{r g'_l}{2Mc},\qquad TMC=2Mc.
\]

The code's `C=274.074` is the atomic-unit value of `2Mc`.  With

\[
 TMC(r)=C-\frac{V(r)-2Z/r-E}{C},
\]

the scalar-relativistic spinor is represented, up to the legacy phase
convention, by

\[
 \Psi_{\kappa m}(r,\Omega)=\frac1r
 \begin{pmatrix}
 G_l\,Y_{lm}\chi\\
 -i\sigma_r\left[\Phi_l-
 \frac{\boldsymbol\sigma\!\cdot\!\mathbf L}{TMC\,r}G_l\right]Y_{lm}\chi
 \end{pmatrix}.
\]

This is the Kölling--Harmon scalar-relativistic reduction, with the
variables mapped to the live `RSEQSR` storage rather than imported as a new
radial equation.  The corresponding Kölling--Harmon reference is:

* [Koelling and Harmon, J. Phys. C 10, 3107 (1977)](https://doi.org/10.1088/0022-3719/10/16/019)

For a same-spin scalar overlap, angular reduction gives

\[
 G_l^2\left[1+\frac{l(l+1)}{(TMC\,r)^2}\right]+\Phi_l^2.
\]

Thus `GFAC` is not an empirical multiplier: it is the squared angular
gradient part of the large numerator, while `G(:,2)` is the independent
minor/lower radial numerator.  `GINTSR` applies this metric with its legacy
log-mesh Simpson rule.

The accepted radial snapshot now retains `potential(:,spin)` and
`tmc(:,l,spin)` in addition to the existing large/small, derivative,
`gfac`, and energy arrays.  No radial equation is recomputed by DRESP-09S.

## Same-spin hard oracle

The unit fixture compares the recovered pointwise expression and its
quadrature against the live `GINTSR` wrapper for `phi-phi`, `phi-phidot`,
`phidot-phi`, and `phidot-phidot` over several l channels.  It also checks
that the accepted up=down limit reproduces the same derived angular metric.
The closure is at machine precision; the test prints the maximum residuals.

## Mixed-spin L=0 contract

For unlike endpoints, the two endpoint factors are retained separately:

\[
 \mathcal R^{L=0}_{pq}(u,d;r)=\frac1{r^2}\left[
 G_{u,p}G_{d,q}
 -\frac13\left(
 \Phi_{u,p}\Phi_{d,q}+
 \frac{l(l+1)G_{u,p}G_{d,q}}
 {TMC_uTMC_d r^2}\right)\right].
\]

The \(-1/3\) is the scalar angular average of the transverse Pauli
operator between the lower spinor pieces,
\(\langle\sigma_r\sigma_\pm\sigma_r\rangle_{L=0}=-\sigma_\pm/3\).
There is no GFAC averaging and no replacement of
`TMC_up*TMC_down` by a spin average.  The symmetric up=down limit is tested
explicitly.

For the four energy branches, define

\[
 X^{(0)}=\phi-\varepsilon_{\rm work}\dot\phi,\qquad
 X^{(1)}=\dot\phi.
\]

The source stores all `00`, `01`, `10`, and `11` endpoint products.  With
`O_pq` the corresponding radial matrix, the field operator is assembled as

\[
 \delta H=\sum_{p,q=0}^1 H^q O_{pq}H^p,
\]

which is the live endpoint orientation used by the native DRESP-08 product
algebra.  The density-side dual uses the effective density
\(H^p\delta\rho H^q\), so that the trace pairing is exact by construction
of the recovered radial observable, not by reusing the DRESP-09R Pauli
metric adjoint.

The source and density paths are implemented independently in
`lr_dresp09s_scalar_relativistic.f90`.

The controlled symmetric fixture also checks the plus/minus circular-channel
Hermitian relation, with the energy-branch indices exchanged under the
transpose.

## Field/density pairing

The controlled radial fixture uses different up/down potentials and
linearisation energies, nonzero l, nonzero large and lower components, and
nonzero `phidot`.  It checks all four branches, the plus/minus channel
construction, the symmetric limit, and

\[
 \operatorname{Tr}(\delta\rho\,\delta H)
 =\int d^3r\,\delta m(\mathbf r)\delta B(\mathbf r).
\]

The accepted Fe 64-k bridge reports a maximum relative field/density pairing
residual of **1.3236e-15**.

## Fe hierarchy

The independent native oracle is the DRESP-08 rigid tangent.  The fixed-basis
SR hierarchy is:

| route | weighted RMS versus native | max relative Frobenius |
|---|---:|---:|
| Pauli large component | 1.1911710463e-1 | 1.6369497261e-1 |
| + lower radial component | 1.1911710463e-1 | — |
| + full SR angular/GFAC term | 1.1914630717e-1 | 1.6372452129e-1 |
| native DRESP-08 tangent | 0 | 0 |

The full SR field is therefore a valid recovered fixed-basis physical
operator, but it does not reproduce the native LMTO tangent by itself.  The
residual is not fitted into the radial metric.  A separate
basis/potential-parameter representation contribution is required.

The finite-rotation oracle independently rotates the accepted native
collinear operators and forms
\([H(+\theta)-H(-\theta)]/(2\theta)\).  Its relative residuals against the
native commutator are:

| theta | residual |
|---:|---:|
| 1.0e-2 | 1.6666583334e-5 |
| 5.0e-3 | 4.1666614587e-6 |
| 2.5e-3 | 1.0416663415e-6 |

The factor-of-four convergence identifies the expected \(\theta^2\) central
difference behavior.  The native commutator itself closes at
**3.3475e-15**.  This is the independent evidence for the
`BASIS_ROTATION_TERM_REQUIRED` classification.

## Density comparison

The exact global coefficient-space rotation is contracted with the same
mixed-spin SR observable and compared with the accepted radial targets:

| density target | relative mismatch |
|---|---:|
| Pauli projected valence | 4.0715438245e-2 |
| full SR valence | 4.0820980180e-2 |
| full SR plus accepted core bookkeeping | 3.0242379630e-2 |

The spectral response contains the valence eigensystem, so the primary
comparison is valence-only.  The core is retained as a separate diagnostic;
the accepted core integral is -1.9267e-5 in the artifact.  The full SR
density remains open until the representation response is added, consistent
with the PASS-B classification.

## Capability boundary

| capability | status |
|---|---|
| Pauli large-component transverse vertex | CERTIFIED; insufficient for native H tangent |
| Same-spin legacy GINTSR metric | CERTIFIED |
| Mixed-spin SR radial/angular L=0 vertex | CERTIFIED internally and dual |
| Full SR L=0 rigid-rotation vertex including representation response | PASS-B; basis term required |
| Full SR arbitrary-L vertex | NOT YET CERTIFIED |
| Native LMTO tangent | CERTIFIED independently by DRESP-08 |
| BES/Halle | OFF |
| ALSDA | OFF |

The next bounded milestone is to derive the explicit LMTO basis/potential-
parameter response term, then rerun this same L=0 field and density gate.
Only after that closes should an L>0 spatial SR vertex be considered.

## Files and tests

Implementation:

* `source/lmto_radial_augmentation.f90` — accepted SR provenance snapshot;
* `source/lr_dresp09s_scalar_relativistic.f90` — radial/angular field and
  density contractions;
* `source/lr_dresp09s_bridge.f90` — independent 64-k material gate;
* `tests/unit/test_dresp09s_scalar_relativistic.f90` — radial, branch,
  symmetric-limit, angular, and duality fixture;
* `tests/integration/tddft_driver_smoke/input_dresp09s_fe.nml` — accepted Fe
  input;
* `tests/validation/dresp09s_fe_artifact.py` — independent artifact checks.

The production backend is `sr_l0_scalar_relativistic`; it requires one
Gamma point, static omega=0, `response_lmax=4`, and the accepted k-space SCF
handoff.  DRESP-08 and DRESP-09R remain separate oracle paths.  No generated
material artifact is committed.

Starting HEAD: `a533b444c9f7d8ad256bc377c57b36749289689c`  
Primary classification: **`BASIS_ROTATION_TERM_REQUIRED`**  
Verdict: **PASS-B**
