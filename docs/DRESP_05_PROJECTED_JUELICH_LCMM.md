# DRESP-05 — Projected Jülich/LCMM interaction and same-state response

Status: implemented on `fable_v4` as a separate projected site-space rung.
The synthetic algebra is certified by `UnitLrProjectedJuelichInteraction`.
The Fe result is a controlled projected-material result: its classification is
determined by rank, conditioning, real-constrained residual, and static-eta
stability, never by whether a spectrum looks attractive.

The literature contract is the local-response construction of Lounis, Costa,
Muniz, and Mills, *Phys. Rev. B* **83**, 035109 (2011), together with the
Katsnelson–Lichtenstein dynamical-response distinction,
*J. Phys.: Condens. Matter* **16**, 7439 (2004).  The implementation derives
the factors from the live RS-LMTO conventions below; it does not import the
radial GSR prefactors by analogy.

## 1. Site Ward identity

DRESP-01 uses

\[
 \sigma^+=(\sigma_x+i\sigma_y)/2,\qquad
 \sigma^-=(\sigma_x-i\sigma_y)/2,\qquad V^z=\sigma_z.
\]

DRESP-02 accumulates the certified circular site response, including its
existing factor of two.  For the exact self-consistent scalar Stoner fixture,

\[
 B_{\sigma,i}=U_iM_i,\qquad
 [I-\chi_0(0,0)U]M=0.
\]

Consequently the projected Ward equation is

\[
 M_i=\sum_j\chi^0_{ij}(0,0)U_jM_j,
 \qquad
 \Gamma_{ij}=\chi^0_{ij}(0,0)M_j,
 \qquad
 \boxed{\Gamma U=M}.
\]

There is no additional `2`, `1/2`, `4*pi`, `mu_B`, or `g` in this site
equation.  The DRESP-04 oracle is the normalization check: its field is
`B_sigma=(H_up-H_down)/2` and its Mills value is `B_sigma/M`.

The new service is
[`source/lr_projected_juelich_interaction.f90`](../source/lr_projected_juelich_interaction.f90).
It keeps `U_ij=U_i delta_ij`, while `chi0_ij` and `Gamma_ij` remain fully site
coupled.  It solves the finite equation independently with LAPACK SVD:

* complex `ZGELSS` solves `Gamma U_complex=M`;
* real `DGELSS` solves the stacked system
  `[Re Gamma; Im Gamma] U_real=[M;0]`.

The production Dyson route remains the DRESP-04
`evaluate_projected_dyson` service.  No second interacting-response solver was
introduced.

## 2. Radial GSR normalization audit

The existing radial service uses, in its own product/point convention,

\[
 \Gamma^{rad}_{ab}=4\pi\chi^0_{ab}(0)m_{00,b},
 \qquad
 \Gamma^{rad}U_{LCMM}=m_{00},
 \qquad K_{eff}=4\pi U_{LCMM}.
\]

The apparent `4*pi` difference is bookkeeping, not a missing site factor.
The radial unknown is a positive-radius function and `m_00(r)` is the radial
coefficient used by the LR response-space metric.  The site quantity `M_i` is
the DRESP-01 site integration functional applied to the Pauli operator; it is
not the pointwise value of `m_00(r)`.  `Y_00`, the radial volume element, and
the site functional are therefore consumed at different stages.

A separable bridge makes the mapping explicit.  Let

\[
 m_{00}(r)=M w(r),\qquad \int dr\,r^2 w(r)=1,
 \qquad U(r)=U_{site}
\]

and apply the same normalized radial test functional to both sides of the
radial response action.  The radial action becomes the site action
`M = chi_site U_site M` after the radial volume, `Y_00`, and the response
functional have been absorbed into `chi_site`.  If instead the full radial
kernel is retained, the object that is comparable to the site scalar is the
functional action of `4*pi U(r)`, not the pointwise radial `U(r)`.  Thus
`U_site=U_radial` is asserted only for the explicitly separable constant-
kernel fixture.  For a general radial interaction, equality of two numbers
called `U` is neither expected nor required.

No production code copies the radial `4*pi` into `Gamma U=M`.

## 3. Complex and real interactions

At zero broadening in a collinear no-SOC state the physical local interaction
is real.  The finite-width DRESP-02 static response may be complex, so the
service reports both constructions:

\[
 U_C=\Gamma^{-1}M,
 \qquad
 \frac{\|\operatorname{Im}U_C\|}{\max(\|\operatorname{Re}U_C\|,\epsilon)},
\]

and the constrained real solution from the stacked system.  Production uses
`interaction_U_real`; it does not silently take `real(U_complex)`.

The static eta ladder is the configured `eta_grid`.  The finest configured
point is selected, the preceding point is retained as a stability comparison,
and the selected real `U` is then frozen for every dynamic eta.  If the two
finest points differ by more than the bounded relative tolerance, the result
is classified `ETA_LIMITED_LOCAL_SUMRULE`.

The construction residual is evaluated at the exact static response used to
solve for `U`.  A holdout residual evaluates frozen `U` against another static
response.  They are separate fields.  The construction residual alone is not
called Goldstone validation.

## 4. Synthetic acceptance oracles

`tests/unit/test_lr_projected_juelich_interaction.f90` contains the required
independent fixtures:

* one-site known-`U` recovery and Mills equality;
* a genuinely coupled two-site matrix with unequal local interactions;
* a separate equal-`U` two-site recovery, so a uniform scalar cannot pass both
  cases accidentally;
* explicit LAPACK LU/inverse comparison for a nonsingular Gamma;
* machine-rank deficiency rejection;
* finite-broadening complex `U`, real-constrained residual, and convergence as
  eta decreases;
* a frozen-`U` holdout/model-mismatch residual.

The two-site fixture uses a full off-diagonal Gamma generated from a coupled
finite-Hamiltonian-style response, not two disconnected one-site copies.
The holdout negative control is deliberately evaluated outside the construction
response, which is the meaningful place to expose an approximation that a
square static equation can otherwise absorb into a diagonal fit.

## 5. Production backend and same-state ledger

The new route is `backend = 'projected_juelich'`, with the bounded input at
[`tests/integration/tddft_driver_smoke/input_dresp05_fe.nml`](../tests/integration/tddft_driver_smoke/input_dresp05_fe.nml).
It performs, in one accepted state:

1. DRESP-01 `d` and/or `spd` projected moments;
2. DRESP-02 static Gamma response over the eta ladder;
3. independent projected Jülich solves and diagnostics;
4. Mills projection from the same reciprocal Hamiltonian;
5. frozen-`U` dynamic DRESP-02 response;
6. shared site Dyson evaluation for bare, Mills, and Jülich routes.

The artifact reports the accepted-state provenance, selector, mesh, Fermi
level, temperature, projected moments, complex and real interactions, rank,
condition number, imaginary ratio, construction residual, holdout residual,
eta stability, raw Mills/Jülich Ward residuals, relative interaction
difference, and explicit dynamic route rows.  Dynamic rows retain the loss,
minimum singular value, condition number, minimum denominator eigenvalue, and
Dyson residual.  The loss field is named `minus_im_trace_over_pi` and means
`-Im Tr(chi)/pi`.

For the first nonzero q, the backend also writes `FREQUENCY_REFINEMENT` rows:
the coarse-grid route peak, a locally refined peak, the refined denominator
minimum, and the refined loss.  These are grid diagnostics only; they do not
fit linewidths or promote a peak to a collective-mode claim.

The primary Fe selector is `spd`; `d` is a controlled projection diagnostic.
Mills and Jülich values are not forced to agree.  The comparison is precisely
the distinction between compressing the Hamiltonian splitting and satisfying
the projected static Ward identity.

## 6. Negative controls and covariance

Rank deficiency is reported as `UNSUPPORTED_RANK_DEFICIENT`.  Invalid or
ill-conditioned cases are `UNSUPPORTED`.  A well-conditioned finite-residual
site-diagonal projection is `PROJECTED_LOCAL_SUMRULE`; an eta-unstable result
is `ETA_LIMITED_LOCAL_SUMRULE`.  `EXACT_LOCAL_SUMRULE` is reserved for a
full-rank real-constrained construction with a small residual and stable
static interaction.

The existing non-self-inverse q/channel test is retained.  For a non-self-
inverse q, the Jülich route checks

\[
 (q,+,\omega)\leftrightarrow(-q,-,-\omega),
\]

with the same real diagonal frozen interaction.  Since the interaction is
static and real, it cannot introduce a new symmetry violation into the bare
response.

The full selected site × site field is now carried by the DRESP-04 reciprocal
adapter.  Intersite field content is consequently visible in Mills
`locality_residual`; it is not silently discarded.  No nonlocal production
Jülich interaction is added in this rung.

## 7. Fe scope and limitations

The accepted first material route is the DRESP-04 bounded reciprocal `4^3`
state with `Gamma`, `(0.03,0,0)`, and `(0.125,0,0)`, with the same dynamic eta
sequence.  This is a same-state comparison, not a converged Fe dispersion.
The static eta ladder and resulting classification are written into the
artifact.  A denser `6^3` or `8^3` repeat is optional evidence only; it does
not silently promote the result to a converged material interaction.

The executed `4^3` artifact gives the following frozen-interaction ledger:

| selector | `M` | `U_Mills` | `U_Juelich_real` | construction residual | holdout residual | eta status | relative Jülich/Mills difference |
| --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| `d` | 2.06232549 | -0.03396582 | -0.03507710 | 0.0635986 | 0.126351 | stable | +0.0327176 |
| `spd` | 2.05455046 | -0.02559760 | -0.03148645 | 0.0969592 | 0.179640 | ETA-limited | +0.230055 |

The corresponding unconstrained complex interactions have imaginary ratios
`0.0637277` (`d`) and `0.0974182` (`spd`).  Mills scalarization residuals are
`0.121600` and `0.440639`.  These numbers are diagnostic comparisons, not
fitted shifts and not evidence that either route is a full ALSDA kernel.
The primary `spd` result therefore makes this material rung `PASS-B` while the
synthetic algebra remains PASS-A quality.

The route does not fit a Goldstone shift, use DRESP-03TG as a target, invoke
ALSDA/BES/GCR, extract intrinsic linewidths, or start longitudinal/charge
response.  DRESP-03TG source files remain frozen.  DRESP-03TG and radial GSR
comparisons remain external normalization/physics diagnostics, not tuning
conditions.

Frequency-grid refinement and correlated pole evidence remain the next
analysis step when a finite-q feature is selected.  A loss maximum alone is
not a collective-mode claim.

## 8. Next rung

After this rung, choose one bounded next milestone: full projected/radial
ALSDA comparison (DRESP-06A), or systematic pole/mode extraction comparing
Mills, Jülich, and DRESP-03TG (DRESP-06B).  Neither is started by DRESP-05.
