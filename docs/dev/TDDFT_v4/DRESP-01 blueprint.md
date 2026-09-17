# DRESP-01 — Projected Site-Spin Operator and Matching Moment Contract

## Status

**Type:** physics-contract certification + narrow implementation
**Branch:** `fable_v4`, live HEAD to be verified
**Depends on:** DRESP-00 architectural rebase
**Produces:** certified projected site-spin vertices and matching projected moments
**Does not yet produce:** projected dynamical susceptibility

---

# 1. Parent scientific question

Can RS-LMTO-ASA define and evaluate, without an uncontrolled approximation, a site-resolved spin operator

$$
\hat S_i^+
=
\int_{\Omega_i}d^3r\,
\hat\psi^\dagger(\mathbf r)\sigma^+\hat\psi(\mathbf r)
$$

and its longitudinal partner using the **same LMTO basis, radial normalization, orbital selection, metric and frozen-core policy** as will later be used to construct

$$
\bar\chi_{0,ij}^{+-}(q,\omega)?
$$

The output must be a certified operator/moment contract for the supported orbital selections:

$$
d,\qquad spd.
$$

`spdf` remains deferred until the one-electron capability gap identified by DRESP-00 is closed.

---

# 2. Why DRESP-01 is separate from projected chi0

Do not yet build or validate the final site × site susceptibility.

The logical hierarchy is

$$
\text{LMTO states}
\rightarrow
\text{site-spin operator}
\rightarrow
\text{projected transition amplitude}
\rightarrow
\bar\chi_0.
$$

If the operator itself is wrong, agreement between two susceptibility implementations can merely reproduce the same wrong projection.

DRESP-01 therefore certifies the operator and matching moment first.

DRESP-02 will subsequently use that certified contract in independent Lehmann and GF calculations.

---

# 3. Required operator definitions

For a site \(i\), define the projected spin operators from the actual field operator inside the site domain \(\Omega_i\):

$$
S_i^+
=
\int_{\Omega_i}
\psi^\dagger(\mathbf r)\sigma^+\psi(\mathbf r)\,d^3r,
$$

$$
S_i^-
=
(S_i^+)^\dagger,
$$

and the corresponding longitudinal/site magnetization operator in the repository's locked spin convention.

Do not insert factors of

* \(2\),
* \(1/2\),
* \(\mu_B\),
* \(g\),
* \(4\pi\),

from memory.

All such factors must come from the live convention ledger and the actual definitions of the existing radial spin density and response operators.

---

# 4. Orbital-selection meaning

`d` and `spd` are **one-electron Hilbert-space selections**, not response-harmonic truncations.

For example:

* `d`: retain the supported \(l=2\) one-electron components;
* `spd`: retain \(l=0,1,2\).

The response created by products of those states may contain higher response harmonics through angular products.

Do not implement `d` by simply keeping a response \(L=2\) channel.

Do not implement `spd` by truncating response harmonics at \(L=2\).

Use the exact Gaunt/product structure implied by the selected one-electron endpoint orbitals.

---

# 5. Matching projected moment

Define

$$
M_i^{(\mathcal P)}
$$

using exactly the same:

* site domain;
* one-electron orbital selector \(\mathcal P\);
* radial functions;
* LMTO representation;
* spin convention;
* normalization;
* valence/core policy

as the transverse projected operator.

The eventual Jülich sum rule will require the same \(M_i^{(\mathcal P)}\) as the projection used to construct \(\bar\chi_0\).

Therefore DRESP-01 must not use the total SCF moment merely because it is readily available.

If `d` is used for response, produce a `d`-projected moment.

If `spd` is used, produce an `spd`-projected moment.

A total/core-inclusive moment may be reported separately for comparison, but it is not automatically the projected moment.

---

# 6. Preferred implementation architecture

DRESP-00 reports that the exact LMTO product-response representation survives as certified infrastructure/oracle.

Use this where appropriate.

A particularly desirable architecture is:

$$
\text{one-electron transition}
\rightarrow
z_{\alpha}
\quad\text{(compact product coordinates)}
$$

followed by a certified site-integration functional

$$
p_{i,\alpha}
$$

so that

$$
T_i^{mn}(k,q)
=
\sum_\alpha
p_{i,\alpha}^{*}
z_{\alpha}^{mn}(k,q).
$$

Here \(T_i^{mn}\) is the matrix element of the projected site-spin operator needed by the future site susceptibility.

This is an **exact contraction of the existing transition representation**, not construction of a 232×232 susceptibility followed by projection.

If another direct formulation in the augmented LMTO basis is cleaner, it may also be implemented, but its equivalence to the response-basis contraction must be demonstrated.

---

# 7. Site-integration functional

Derive the functional that maps the accepted response representation to

$$
\int_{\Omega_i}d^3r\,\delta m^+(\mathbf r).
$$

The derivation must explicitly account for:

* radial measure;
* response metric;
* radial quadrature;
* spherical-harmonic normalization;
* the angular component surviving full-sphere integration;
* site indexing;
* any factors implied by complex versus real spherical harmonics.

Do not assume the surviving angular coefficient is numerically equal to the site integral.

In particular, derive any \(Y_{00}\), \(\sqrt{4\pi}\), or equivalent normalization from the live basis conventions.

---

# 8. Finite-q operator contract

The operator must support the same finite-q convention as the certified reciprocal transition machinery.

Do not invent a new phase convention.

Prefer:

$$
T_i(k,q)
=
p_i^\dagger z(k,q)
$$

using the already certified transition vector carrying the established \(k\rightarrow k+q\) phase/gauge information.

Document whether the site position \(\tau_i\) phase is already contained in the transition representation.

At minimum certify:

* \(q=0\);
* one nonzero mesh-compatible q on a controlled fixture.

Full \(\chi_0(q,\omega)\) testing remains DRESP-02.

---

# 9. Direct one-electron operator oracle

Where the existing augmented LMTO machinery permits, construct the corresponding operator matrix directly:

$$
V_i^\mu
$$

with matrix elements obtained from the selected augmented basis functions inside \(\Omega_i\).

This direct route should be independent of the compact product-coordinate contraction as far as practical.

Use it as an oracle for transition amplitudes:

$$
\langle n k|V_i^+|m,k+q\rangle.
$$

If an independent direct construction cannot be made without introducing new unproven physics, document this and retain the product-space construction as the certified route.

Do not manufacture a fake independent oracle.

---

# 10. Required algebraic identities

For every supported projected space test:

## 10.1 Adjoint

$$
V_i^-=(V_i^+)^\dagger.
$$

## 10.2 Longitudinal Hermiticity

The projected longitudinal spin/magnetization operator must be Hermitian.

## 10.3 Site support

A vertex assigned to site \(i\) must carry only the intended site-domain contribution.

## 10.4 Orbital selector

`d` and `spd` must correspond to their declared one-electron orbital content.

No material-dependent pruning.

## 10.5 Representation consistency

Where both direct and compact-product routes are available,

$$
T_{i,\mathrm{direct}}^{mn}
\approx
T_{i,\mathrm{product}}^{mn}
$$

to numerical precision appropriate to the already certified LMTO/product transformation.

---

# 11. Ground-state moment certification

For an accepted collinear material state calculate the projected moment in at least two genuinely equivalent but independently evaluated ways where possible:

### Route A — ground-state density

Project/integrate the accepted spin density with the same site/orbital definition.

### Route B — one-electron density matrix

Evaluate the expectation value of the certified longitudinal vertex:

$$
M_i^{(\mathcal P)}
=
\mathrm{Tr}
[
\rho\,V_i^z
]
$$

with the exact repository convention restored.

These two must agree within controlled numerical tolerance.

If factors such as \(\mu_B\) distinguish a spin-density moment from a physical magnetic moment, report both explicitly rather than absorbing the conversion.

---

# 12. Core policy

DRESP-01 must produce an explicit answer to:

> Is the projected response a valence-only object, and if so, is the matching \(M_i\) also valence-only?

If frozen-core spin density exists in the ground state, it must not silently enter the projected Jülich moment.

Report separately:

* selected-valence projected moment;
* other valence contribution outside the selector, if any;
* frozen-core contribution, if available;
* total SCF moment.

Do not combine them unless the future Ward identity is derived for that combined quantity.

---

# 13. Material validation

Use a real accepted collinear material state after exact fixture tests.

Primary material:

**bcc Fe**

For DRESP-01 we require only operator/moment validation—not a response spectrum.

For each supported selection:

```text
d
spd
```

report:

* accepted electronic-state provenance;
* site moment from direct density projection;
* site moment from operator/density-matrix evaluation;
* difference;
* total SCF moment for context;
* core contribution where available.

Do not tune the projection so that the projected moment equals the total SCF moment.

That equality is not the acceptance criterion.

---

# 14. Required tests

Create lean independent tests covering:

1. site-integration functional normalization;
2. angular normalization/full-sphere integral;
3. `d` selector;
4. `spd` selector;
5. \(V^-=(V^+)^\dagger\);
6. longitudinal Hermiticity;
7. direct-versus-product transition amplitudes where available;
8. q=0 transition projection;
9. one finite-q transition projection;
10. projected moment from density versus operator;
11. real bcc-Fe `d` projected moment;
12. real bcc-Fe `spd` projected moment.

Synthetic tests may certify algebra.

They do not replace the Fe execution.

---

# 15. Explicit non-goals

Do not implement:

* site × site Lehmann susceptibility;
* projected GF susceptibility;
* K/L dynamic exchange;
* Mills \(U=\Delta/M\);
* Jülich \(U\);
* GSR modifications;
* ALSDA modifications;
* BES/GCR;
* Dyson/loss changes;
* magnon dispersion;
* linewidth extraction;
* native-RSGF production integration;
* `spdf`;
* longitudinal response.

DRESP-01 stops at certified projected operators, transition amplitudes and matching moments.

---

# 16. Hard blockers

Stop and return evidence if:

1. site integration cannot be defined in the current response metric;
2. the orbital selector cannot be mapped unambiguously to augmented LMTO states;
3. direct density and operator definitions use irreconcilable normalizations;
4. core and valence cannot be separated sufficiently to define the projected moment;
5. finite-q projection requires inventing a new phase convention;
6. product-coordinate contraction loses information needed by the site integral;
7. the projected moment cannot be defined with the same Hilbert-space selector as the transverse operator.

Do not patch around these.

---

# 17. Required documentation

Create:

`docs/DRESP_01_PROJECTED_SITE_SPIN_CONTRACT.md`

It must document:

1. convention mapping;
2. operator definition;
3. orbital-selector semantics;
4. site-integration functional;
5. response-metric treatment;
6. finite-q convention;
7. direct operator construction if available;
8. matching projected-moment definition;
9. core policy;
10. algebraic test results;
11. Fe `d` and `spd` evidence;
12. blockers/limitations;
13. exact API exposed for DRESP-02.

---

# 18. Completion criteria

DRESP-01 passes only if, for `d` and `spd`:

$$
\boxed{
\text{same projector}
+
\text{same normalization}
+
\text{same site domain}
+
\text{same core policy}
}
$$

are used consistently for transverse transition vertices and projected moments.

We must leave DRESP-01 with a routine/API capable of providing

$$
T_i^{mn}(k,q)
=
\langle nk|S_i^+|m,k+q\rangle
$$

and

$$
M_i^{(\mathcal P)}
$$

without constructing a large response matrix.

---

# 19. Expected next milestone

If DRESP-01 passes:

**DRESP-02 — projected reciprocal bare susceptibility**

will accumulate

$$
\bar\chi_{0,ij}^{+-}(q,\omega)
$$

from the certified site transition amplitudes and compare independent Lehmann and reciprocal-GF evaluations.

Do not begin DRESP-02 automatically.
