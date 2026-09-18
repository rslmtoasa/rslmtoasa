# DRESP-04 — Projected Mills/Stoner Site-Space RPA

Status: implementation complete on `fable_v4`; synthetic algebra and the
bounded bcc-Fe production campaign pass.  The material result is classified
`PASS-B`: the interaction and Dyson algebra are closed, while the Fe local
scalarization is explicitly retained as a projected approximation.

This deliverable is a controlled site-space response layer.  It consumes the
accepted reciprocal no-SOC collinear SCF state, the DRESP-01 projected spin
contract, and the DRESP-02 direct site-matrix `chi0`.  It does not reopen the
DRESP-03TG native exchange construction and does not add a Goldstone repair,
Jülich/LCMM sum rule, ALSDA kernel, or force-theorem exchange vertex.

## Closed equations and conventions

The projected field is defined by

\[
 H_{\uparrow}(k)-H_{\downarrow}(k)=2B_\sigma(k),\qquad
 H(k)=H_0(k)I+B_\sigma(k)\sigma_z .
\]

The DRESP-01 `Vz` operator is `sigma_z`, so the fitted splitting is
`Delta_i = B_sigma`, not the full up/down difference.  The local scalar Mills
fit is

\[
 B_\sigma(k) \simeq \sum_i U_i V_i(k),\qquad
 U_i=\frac{\Delta_i}{M_i},
\]

where `M_i` is the same DRESP-01 projected moment used by the site contract.
The fit is a weighted complex least-squares problem over all selected matrix
entries and accepted reciprocal samples.  The implementation reports the
fit condition, scalarization residual, off-site locality residual, and any
imaginary coefficient residual.  It classifies the result as:

| classification | meaning |
| --- | --- |
| `EXACT_SCALAR` | local scalar representation passes the configured tolerance |
| `PROJECTED_SCALAR_APPROXIMATION` | a finite projection residual is exposed |
| `UNSUPPORTED` | the mapping is rank-deficient, ill-conditioned, or invalid |

The site response is formed directly from DRESP-02:

\[
 D(\omega)=I-\chi^0(\omega)U,\qquad
 D(\omega)\chi(\omega)=\chi^0(\omega).
\]

`U` is diagonal in site space.  Production calls the existing certified
LAPACK solve in `source/tddft_dyson.f90`; it does not duplicate an explicit
matrix inversion.  Each frequency reports solve information, relative and
absolute Dyson residuals, singular-value conditioning, and the minimum
magnitude eigenvalue of `D`.

The loss convention is the ordinary site matrix

\[
 L(\omega)=-\frac{\chi(\omega)-\chi^\dagger(\omega)}{2i\pi},
 \qquad -\pi\,\operatorname{Im}\operatorname{Tr}\chi
 \]

with no Hermitian projection or pole repair.  The `chi_plus` production path
also audits the non-self-inverse relation between `(+q,+omega,+)` and
`(-q,-omega,-)`.

## Implementation seam

The new service is
`source/lr_projected_interacting_response.f90`, registered immediately after
the DRESP-02 projected chi0 service in `source/CMakeLists.txt`.  Its two
public operations are:

1. `evaluate_projected_mills_from_reciprocal`, which extracts the selected
   `d` or `spd` coefficient field from the accepted reciprocal Hamiltonian and
   fits the local interaction;
2. `evaluate_projected_dyson`, which consumes the DRESP-02 site matrix and
   returns the enhanced response, denominator diagnostics, and loss data.

The production branch is `backend = 'projected_mills'`.  The input selector
can be `d`, `spd`, or `both`; the bounded fixture uses both so the requested
primary `Fe spd` result and the `d` diagnostic are produced in one artifact.
The general `interaction_route` namelist field remains syntactically required
by the shared parser, but it is ignored by this backend and is recorded as
Mills/Stoner in the output metadata.

The production writer records provenance for the accepted state, DRESP-02
bare response, projected interaction, Dyson convention, broadening, channel,
classification, and every matrix element.  It records `RAW_GAMMA` diagnostics
at `q=0`; these are measurements of the uncorrected denominator, not a
Goldstone-enforcing operation.

## Synthetic acceptance oracles

`tests/unit/test_lr_projected_interacting_response.f90` contains independent
one-site and two-site fixtures.

The one-site fixture uses `U=-0.6`, `M=0.7`, and `B_sigma=-0.42`; the fit
recovers both `U` and `B_sigma` to machine precision.  The two-site fixture
uses a hopping Hamiltonian with `U=-0.8`, `M_i=0.5`, and
`B_sigma=-0.4`.  Its bare circular response is assembled independently from
finite Hamiltonian eigenpairs.  The test verifies:

- self-consistent uniform static Goldstone algebra,
- enhancement by the site Dyson solve,
- agreement with an independent `zgetrf/zgetri` inverse oracle,
- Hermiticity of the loss matrix,
- non-self-inverse q/channel covariance,
- and explicit `PROJECTED_SCALAR_APPROXIMATION` classification when an
  off-site field is injected.

No synthetic test changes the production result or repairs its denominator.

## Bounded Fe campaign

The reproducible input is
`tests/integration/tddft_driver_smoke/input_dresp04_fe.nml`.  It uses the
certified bcc-Fe database, a converged reciprocal `4x4x4` accepted state,
`q = Gamma, (0.03,0,0), (0.125,0,0)`, `omega = 0..0.20 Ry`, and
`eta = 0.04, 0.02, 0.01 Ry`.  The run completed with exit code zero and
produced both selectors.

Representative Mills diagnostics from that run are:

| selector | projected moment | `B_sigma` | `U_Mills` | scalarization residual |
| --- | ---: | ---: | ---: | ---: |
| `d` | 2.3104702388 | -0.0794941878 | -0.0344060644 | 0.1323955344 |
| `spd` | 2.2034798883 | -0.0630997637 | -0.0286364146 | 0.4168281989 |

Both fits are local (`locality_residual = 0`) for this one-site bcc primitive
cell, but neither is classified exact: the finite scalarization residual is
the material approximation being measured.  The requested `spd` primary
result is therefore not presented as a faithful scalar exchange kernel.

At Gamma, the raw denominator minimum singular values for `d` are about
`0.2152`, `0.1118`, and `0.0601` for the three broadenings; for `spd` they
are about `0.3018`, `0.2341`, and `0.2128`.  The enhanced response and loss
change substantially as `eta` is reduced, demonstrating broadening-sensitive
low-energy spectral enhancement.  The bounded run does not fit or claim a
converged collective pole: the raw denominator diagnostics remain the
reported evidence.

There is no same-q DRESP-03TG artifact consumed by this bounded campaign, so
no numerical DRESP-03TG equality is asserted here.  DRESP-03TG source files
remain frozen.  A future comparison must use the same accepted state, q,
frequency, channel, and broadening before any cross-backend conclusion is
drawn.

## Verification and next gate

The focused CTest target is `UnitLrProjectedInteractingResponse`; the
production binary also completed the Fe fixture above.  The committed
changes add the interaction service, its independent unit oracle, the
production input branch, the namelist selector, and this record.  The next
scientific gate is not another algebraic repair: it is a same-state
cross-backend comparison and, separately, a physical pole/residue analysis
with explicit convergence in q, omega, and eta.
