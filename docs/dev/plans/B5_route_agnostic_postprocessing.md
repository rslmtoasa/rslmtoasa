# B5 — Route-agnostic G(E) post-processing: J_ij, damping, transport

**Effort:** S–M. **Lead:** OPUS for the moment generator and the damping
audit, SONNET for dispatch/regression plumbing. **Depends on:** B2 (+B4 for
speed). **Design stance:** most of this milestone falls out of B2's
fill-the-canonical-arrays contract **by construction**; this plan covers
only what doesn't, plus the regression triads that lock everything down.

## 1. Scope (per BLUEPRINTS §B5 — no invented abstractions)

1. **Exchange & Gilbert damping: no code changes.** They consume
   `gij/gji/gij_eta` + torque components; B2 fills them. The only new code
   is the `gf_route = recursion | lehmann | dyson` dispatch key — already
   delivered by task B2.5. This plan's job is the *validation triads*, not
   plumbing.
2. **Conductivity (moment-native): exact moment generator.** The KPM
   pipeline consumes Chebyshev moments, not gij. In the eigenbasis,
   T_n(H̃(k)) = U T_n(ε̃) U†, so the moments μ_n are computed **exactly**
   from eigenpairs and written into the *same* moment arrays consumed by
   `finish_conductivity_moments` (`conductivity.f90`). The entire
   downstream Kubo machinery is reused unchanged; the Lehmann route becomes
   an exact moment producer. Respect the Chebyshev energy rescaling
   (`chebfermi`/bounds conventions in `recursion_chebyshev.f90` /
   `bounds.f90`) — moments live on the scaled spectrum; document the
   mapping in the routine header.
3. **Damping audit (report first, code second).** The Kamberský
   torque-correlation implementation runs on the eta ladder; after B2.3 it
   runs on Lehmann input as-is. Task B5.3 first *audits* whether the
   current implementation contains the SOC-derivative torque operators or
   only the exchange-torque flavor — the outcome is a short report; new
   torque operators are implemented only if the audit says they're missing
   and Anders approves the scope.

## 2. Validation — the triads (this is the deliverable)

Same-system triads on bcc Fe (and one noncollinear case for the torque
components): each observable computed via **recursion vs Lehmann vs
Dyson(Σ=0)**, pinned as regression cases with documented agreement
envelopes:

| Observable | Envelope basis |
|---|---|
| J_ij(R), first 5 shells | k-mesh convergence from B2.6 (G-B2-2 thresholds) |
| Gilbert α | η-extrapolation band; the η → 0 protocol is documented once, here, in `docs/dev/route_agnostic_estimators.md` |
| σ (Kubo, KPM) | exact-vs-recursion moment error bound — the direct KPM error measurement, itself a regression test |

## 3. Tasks
- **B5.1 [OPUS]** Exact moment generator + wiring into
  `finish_conductivity_moments`; exact-vs-recursion moment test (bcc Fe;
  report the KPM error bound in the test log). *Kit:* this file; B2 module
  interface; `conductivity.f90` moment-consumer region;
  `recursion_chebyshev.f90` moment array format + rescaling conventions.
- **B5.2 [SONNET]** The J_ij and σ triad regression cases + envelopes doc.
  *Kit:* B2.5/B2.6 outputs; B5.1; existing regression case as template.
- **B5.3 [OPUS]** Damping audit report (SOC-torque presence); α triad with
  the documented η protocol; escalate scope if operators are missing.
  *Kit:* the damping/exchange torque routines (locate via DEVELOPER_MAP —
  grep `gij_eta` consumers in `exchange.f90`); this file §1.3, §2.

## 4. Checklist
- [x] B5.1 exact moments + KPM error bound test **(2026-07-16)**
- [ ] B5.2 J_ij + σ triads pinned
- [ ] B5.3 damping audit + α triad (scope decision if needed)

### B5.1 notes (done)

- New dependency-free `moment_kernel.f90` (`moment_onsite_block`): exact
  eigenbasis double-Chebyshev moment `<site| T_m(H~) v_a T_n(H~) v_b |site>`
  via T_p(H~(k)) = U diag(T_p(eps~)) U†. New `reciprocal_moments` submodule
  (`reciprocal%fill_moments`) fills the SAME `mu_nm_stochastic(:,:,n,m,ntype)`
  the recursion route fills, so `calculate_conductivity_tensor` is untouched.
- The k-space velocity is `v_{a,b}(k) = FT[v_{a,b}(R)]` — the recursion route's
  own real-space velocity blocks (`build_realspace_velocity_operators`)
  Fourier-summed with the SAME neighbor map/phase as H(k)=FT[ee]. So both routes
  evaluate ONE operator in two representations; they agree in the joint
  large-cluster / dense-mesh limit.
- Rescaling: `a,b` passed from the recursion route's `resolve_chebyshev_window`
  (a=(emax-emin)/(2-0.3), b=(emax+emin)/2) — identical scaling to the recursion
  moments; only H is scaled, velocities are applied unscaled.
- Scope: cond_type='charge', cond_calctype='per_type', first-order H(k) (S=I /
  strict-Lehmann regime); dispatched by `gf_route` in
  `post_processing_conductivity` (default 'recursion' = bit-identical legacy).
  'dyson' here coincides with 'lehmann' (Σ=0); Σ≠0 transport needs vertex
  corrections (B8) and is out of scope.
- **KPM error bound.** Unit test `UnitMomentKernel` pins the eigenbasis identity
  against a direct matrix-Chebyshev reference to ~5e-16 (machine precision).
  End-to-end on bcc Fe (6³ mesh, cond_ll=30): diagonal σ_xx from
  recursion vs exact k-space agree to ~0.15% near E_F (3.240 vs 3.235), well
  inside the coarse-mesh envelope; off-diagonal σ_xy is a symmetry-forced zero
  (no SOC) and is noise on an unreduced coarse mesh — NOT a route mismatch. The
  pinned σ-triad regression case is B5.2.
