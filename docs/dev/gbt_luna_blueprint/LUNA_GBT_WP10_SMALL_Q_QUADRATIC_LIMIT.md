# Luna Prompt — GBT-WP10: Small-q Quadratic Limit and q-Reversal Symmetry

## Mission

Demonstrate the controlled long-wavelength limit of the SOC-free GBT energy and extract a stable quadratic curvature before comparing with LKAG.

## Dependency

WP09 must identify a converged k mesh and harmonic cone-angle regime.

## Required q sampling

Choose a lattice direction with simple symmetry first. Use symmetric points:

\[
\pm q_1,\pm q_2,\ldots,\pm q_n
\]

including sufficiently small q values to probe the quadratic regime and somewhat larger values to expose q^4 corrections.

Do not fit only positive q.

## Symmetric/antisymmetric decomposition

For every pair define:

\[
E_{even}(q)=\frac{E(q)+E(-q)}{2},
\]

\[
E_{odd}(q)=\frac{E(q)-E(-q)}{2}.
\]

In the scalar-relativistic SOC-free centrosymmetric benchmark, `E_odd` must be numerical noise.

A systematic odd component is a red flag for q sign, Fourier convention, incomplete bond reversal, or k-grid asymmetry.

## Fit model

Fit the even energy difference to

\[
\Delta E_{even}(q)=Aq^2+Bq^4
\]

and compare against a pure quadratic fit over progressively smaller windows.

Report:

- A and uncertainty;
- B and uncertainty;
- residuals;
- stability of A as the maximum fitted q is reduced;
- sensitivity to k mesh and cone angle.

Do not declare a stiffness from one arbitrary fit window.

## Dimensional conventions

Document q units carefully. Convert to physical inverse-length units where needed before reporting a stiffness.

Keep code-internal units and physical units both visible so no `2*pi/alat` factor can hide.

## Modes

Primary closure should use the method intended for comparison with LKAG in WP11, most likely bare MFT first because standard LKAG is a force-theorem response around the collinear state. Also show corrected-MFT/SCF curvature as separate physics, not as a silent replacement.

## Optional anisotropy

After one high-symmetry direction is established, repeat along at least one additional direction if useful to test cubic isotropy / stiffness tensor behavior.

## Deliverables

- q-sweep and fit utility;
- raw tabular data;
- `docs/dev/GBT_REAUDIT_WP10_SMALL_Q.md`.

## Completion checklist

- [ ] Symmetric ±q points used.
- [ ] Even and odd energy components reported.
- [ ] Odd component is within declared numerical tolerance.
- [ ] q units and conversions documented.
- [ ] q^2+q^4 fits performed.
- [ ] Fit-window stability assessed.
- [ ] k-grid sensitivity assessed.
- [ ] cone-angle sensitivity assessed.
- [ ] Bare/corrected/SCF curvatures clearly distinguished where available.
- [ ] Stable long-wavelength coefficient identified or failure documented.

## Commit

`test: establish quadratic small-q GBT limit`
