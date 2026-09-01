# RS-LMTO-ASA linear-response TD-DFT master blueprint

## 1. Executive objective

Build a scientifically defensible and numerically efficient linear-response TD-DFT subsystem for RS-LMTO-ASA, initially for collinear magnets without spin-orbit coupling, with a single response/Dyson layer fed by **three independently testable Kohn-Sham response backends**:

- explicit eigenpair transitions;
- K-space Green functions `G(k,z)` from the newer Lehmann GF path;
- native real-space Green functions `G(R,z)` from the mature RS machinery.

The real-space GF path is a first-class production route, not a conversion step to K space. For periodic systems it should construct

`G(R,z) -> chi0(R,omega) -> Fourier transform -> chi0(q,omega)`

and for films/embedded/finite systems it should retain the natural mixed or fully real-space representation.

The project must solve the present scientific issue — especially Goldstone/Ward consistency and the noisy long-wavelength spectrum — before adding feature breadth.

---

## 2. Scientific sources and what they imply for this code

### 2.1 Buczek 2025 habilitation (user-supplied `TillAnders.pdf`)

The habilitation is the most useful current overview of the Halle line of LRTDDFT work. The design consequences for RS-LMTO-ASA are:

- The physical target is the full generalized charge/spin response, with transverse spin response as the clean first sector for collinear SOC-free magnets.
- The KKR/Halle method evaluates the Kohn-Sham response from products/convolutions of one-electron Green functions.
- Complex-energy evaluation moves the expensive electronic GF work away from real-axis singularities and can improve convergence very strongly, at the cost of careful continuation/integration strategy.
- Real-space GF formulations are intrinsically suitable for surfaces, interfaces, defects, and systems without full translation symmetry.
- Longitudinal response is a meaningful later target in its own right (order-parameter fluctuations / “magnetic second sound”), not merely an LLB fitting trick.
- Noncollinear and relativistic response are substantially harder and should not be obtained by merely feeding SOC eigenvectors into a scalar transverse formula.

### 2.2 Buczek, Ernst & Sandratskii, PRB 84, 174418 (2011)

This is the primary Halle implementation reference for bulk/film magnons and Landau damping. Important lessons:

- Use a full first-principles KS electronic structure; do not inherit the simplifying d-only/energy-independent basis ansatz of some earlier local-susceptibility schemes unless explicitly chosen as an approximation.
- The physical susceptibility follows a Dyson-like equation with a KS response and the Hartree/xc kernel.
- In collinear SOC-free systems the transverse sector decouples from charge/longitudinal response.
- The magnon is best identified through the **loss matrix/eigenmodes**, not solely through a scalar trace.
- The exact theory contains the Goldstone mode through the identity linking static `chi_KS`, the ground-state exchange-correlation field, and the ground-state magnetization.
- Their implementation corrects a small numerical residual in the Goldstone eigenvalue; in RS-LMTO-ASA this must be treated as a controlled numerical fallback/diagnostic, not as permission for arbitrary fitting.
- Their complex-plane integration and analytic continuation are powerful but continuation is numerically delicate.
- KS response evaluation is the dominant cost and is naturally highly parallel.

### 2.3 Lounis, Costa, Muniz & Mills, PRL 105, 187205 (2010) and PRB 83, 035109 (2011)

These papers are the main Jülich/local-KKR reference for a real-space implementation and Goldstone consistency.

Key consequences:

- The KS transverse susceptibility can be written directly as a convolution of retarded/advanced one-electron Green functions between sites.
- The energy integral can be decomposed into an analytic part that can be moved to a regular complex contour and a small nonanalytic near-real-axis part. This is attractive because it can reduce real-axis work **without requiring a global analytic continuation of the final susceptibility**.
- They derive a sum rule of the form `chi0(0) * B_eff = m` (with spatial/site integration and convention-dependent factors) and use it to construct an effective local interaction compatible with the Goldstone mode.
- Ad-hoc scalar adjustments become especially dangerous for multi-component systems.

For RS-LMTO-ASA, use their sum-rule structure as a **Ward identity and discretization test**. If a site-scalar kernel cannot satisfy the identity under systematic numerical convergence, that is evidence that the response basis/kernel representation is too coarse; do not simply tune it.

### 2.4 Savrasov, PRL 81, 2570 (1998)

Savrasov provides the LMTO-native alternative: a variational time-dependent Sternheimer formulation in a muffin-tin-orbital representation. It avoids explicit high-energy state summation and large response-matrix inversion by solving directly for first-order wavefunction/density/magnetization changes under a monochromatic perturbation.

This is strategically valuable for RS-LMTO-ASA because it is an LMTO-family response formulation and naturally treats the self-consistent first-order potential. It should, however, be a **later independent response engine**, after the current Dyson/GF route is made trustworthy.

### 2.5 Binci, Marzari & Timrov, npj Comput. Mater. 11, 100 (2025)

Modern TDDFPT reinforces the Savrasov direction and adds two lessons:

- Sternheimer/Liouville-Lanczos avoids explicit empty-state sums.
- A Liouville-Lanczos formulation can obtain a broad frequency spectrum from one iterative response problem rather than a self-consistent solve at every frequency.
- Their Goldstone behavior is tied to using the same ground-state and response ingredients consistently rather than empirical kernel rescaling.

This should inform the later design study, particularly if Hubbard response is eventually required.

### 2.6 Relativistic and longitudinal references

- dos Santos Dias et al., PRB 91, 075405 (2015): SOC breaks global spin-rotation invariance and introduces a different sum-rule structure; use this before any relativistic TD-DFT claim.
- Buczek et al., PRB 101, 214420 (2020): longitudinal charge-spin response, collective order-parameter fluctuations, and separation from plasmon energy scales.
- Wysocki et al., PRB 96, 184418 (2017): complete wave-vector/frequency spin-density fluctuation tensor and fluctuation-dissipation perspective.
- Sandratskii & Buczek, PRB 85, 020406(R) (2012), and Odashima et al., PRB 87, 174420 (2013): AF/ferri chirality and strongly channel-dependent Landau damping.

---

## 3. Physics contract

### 3.1 General response

Use the generalized density basis `alpha in {0,x,y,z}` with

`chi = chi_KS + chi_KS (v_C + K_xc) chi`.

For a collinear magnet along z and SOC=0, the response separates into:

- transverse block `(x,y)` or equivalently circular `(+,-)` channels;
- coupled charge/longitudinal block `(0,z)`.

The first production target is the transverse block.

### 3.2 Circular-channel conventions are a release-critical invariant

The implementation must document one exact convention for:

- Pauli matrices;
- `sigma_+/-` versus spin operators `S_+/-` and any factors of 1/2;
- magnetization sign and units;
- Bohr magneton factors;
- retarded Fourier transform sign;
- definition of the loss matrix;
- which of `chi^{+-}` or `chi^{-+}` corresponds to up-to-down versus down-to-up transitions in the code convention.

Do not infer these from paper notation. Derive them from the current RS-LMTO Hamiltonian and magnetization conventions and lock them with unit tests.

### 3.3 Keep the source/measurement vertex separate from the interaction kernel

This is crucial.

The KS susceptibility is the response of a chosen spin-density observable to an external magnetic perturbation. The interaction kernel is the functional derivative of the self-consistent Hartree/xc field with respect to density/magnetization. They are related by Ward identities but are not interchangeable.

Do **not** put `dH_xc/dm` on both ends of `chi0` and then also multiply by `Kxc` in the Dyson equation unless a derivation explicitly proves that this is the chosen response definition; doing so risks double-counting the xc kernel.

Instead define explicitly:

- source vertex `Gamma_ext = dH/dB_ext` in the LMTO basis;
- measurement projector/operator for local spin density or site moment;
- interaction kernel `K = dB_Hxc/dn` or the corresponding operator in the same response basis.

A “kernel-weighted vertex” may be a useful computational composite for forming `chi_KS K`, but it must not obscure this separation.

### 3.4 LMTO consistency requirement

The response objects must correspond to the **same Hamiltonian actually used in the ground state**. In particular:

- same energy zero and Fermi level;
- same spin convention;
- same orthogonalization/metric assumptions;
- same XC bundle and local exchange-correlation field;
- same site/orbital projection weights;
- same Hamiltonian branch (`ham_only` versus generalized-overlap/HOH etc.).

A powerful local test is a finite-difference rotation of the spin-dependent Hamiltonian: rotate a local exchange field by a tiny angle and verify that the analytic transverse perturbation operator reproduces the linear change of H.

### 3.5 Initial support boundary

Production-supported initially:

- collinear magnetic ground states;
- SOC off;
- orthogonal/`ham_only` electronic problem;
- standard DFT XC paths for which a consistent transverse kernel is available;
- site response projection, with optional site-orbital diagnostic basis if needed.

Explicitly reject or guard until separately derived:

- generalized-overlap/HOH response;
- GBT response states;
- CCOR if it changes the response operator and is not differentiated consistently;
- Hubbard-corrected TD-DFT;
- general noncollinear response;
- SOC/relativistic response.

Do not “support” these cases by running them through the collinear formula.

---

## 4. Common TD-DFT architecture

The code should have one physics layer and multiple `chi0` providers.

Conceptually:

```text
Ground-state electronic structure
        |
        +-------------------------------+
        |                               |
 response basis / vertices       one-electron propagators
        |                               |
        |         +---------------------+-------------------+
        |         |                     |                   |
        |   eigenpair chi0         G(k,z) chi0         G(R,z) chi0
        |         |                     |                   |
        +---------+---------------------+-------------------+
                                  |
                            common chi_KS
                                  |
                    Ward / sum-rule diagnostics
                                  |
                         common Kxc / Dyson
                                  |
                           physical chi(q,w)
                                  |
                    loss matrix / mode analysis
```

Recommended logical modules/types (adapt names to existing code style rather than imposing these literally):

- response configuration / validation guards;
- response basis and projector;
- source-vertex/kernel provider;
- abstract or factory-selected `chi0` backend;
- Dyson solver;
- loss-matrix/spectral analyzer;
- Fourier/mixed-representation helper;
- convergence/provenance report.

Follow the codebase’s current modular, class-like Fortran style. Avoid creating an inheritance hierarchy solely for aesthetic reasons; use the lightest interface that allows backend interchange and batched evaluation.

### 4.1 Backend API must permit batching

Do not lock the interface to a scalar `evaluate(q,omega)` call if that would prevent efficient implementations.

The contract should support at least:

- one `(q,omega)` point;
- a batch of frequencies at fixed q;
- a batch of q points at fixed frequency;
- for real-space GF, building `chi0(R,omega)` once and Fourier transforming to many q points.

The backend should return explicit metadata/provenance: broadening/contour settings, k mesh or R cutoff, temperature/smearing, convergence status, and whether any Goldstone correction was applied.

---

## 5. Backend A — explicit eigenpair transitions

Keep the existing eigenpair backend as the transparent reference implementation.

For the transverse channel, it should reduce to the usual occupied/unoccupied (or finite-temperature occupation-difference) transition expression connecting `(n,k)` to `(m,k+q)` with the correct site/spin vertices.

Requirements:

- endpoint provenance for both k and k+q states must never be lost;
- finite-temperature occupation differences and degeneracies handled safely;
- exact q mapping in the reciprocal basis used by the Hamiltonian;
- static limit evaluated as an actual static response, not inferred from a large artificial broadening;
- no hidden Goldstone energy shift;
- observables derived from the full complex response matrix.

This backend is the first reference for backend-equivalence tests.

---

## 6. Backend B — K-space GF via Lehmann

The newer K-space GF is useful but must first be validated as an electronic propagator.

For an orthogonal basis:

`G(k,z) = [z - H(k)]^{-1}`

and its Lehmann form must agree with direct inversion at complex energies to numerical precision.

For generalized overlap the correct object is instead based on `z S(k) - H(k)` and requires metric-consistent spectral weights. Until this is implemented and tested, TD-DFT must keep the generalized-overlap guard.

The KS response is obtained from products of `G(k,z)` and `G(k+q,z+omega)` with retarded/advanced structure and the same response vertices used by the other backends.

Why keep this backend even though it uses eigenpairs internally?

- It validates the GF convolution/integration machinery independently of the explicit transition denominator implementation.
- It provides a bridge between direct band-sum and native RS-GF formulations.
- It allows direct comparison of complex-energy algorithms at fixed k.

It is not as algorithmically independent as the native RS GF and should not be the sole “GF validation”.

---

## 7. Backend C — native real-space GF

This is the key new production path.

### 7.1 No mandatory G(R)->G(k) transform

For sites/basis atoms a,b and lattice vector R, construct the real-space KS response directly from products of one-electron GFs connecting the two sites, schematically

`chi0_ab(R,w) ~ integral dE [ G^down_ab(R,E+w) Im G^up_ba(-R,E) + Im G^down_ab(R,E) G^A_up_ba(-R,E-w) ]`

with the exact factors and vertex contractions fixed by the convention audit.

Then, for a periodic bulk crystal,

`chi0_ab(q,w) = sum_R exp(-i q.R) chi0_ab(R,w)`.

This is the q-space object passed to the common Dyson solver.

The electronic GF itself need never be Fourier transformed to k space.

### 7.2 Natural generalization

- bulk: full 3D FT of `chi0(R,w)`;
- film/slab: keep layer/site indices and FT only R_parallel -> q_parallel;
- impurity/cluster: retain `chi_ij(w)` directly without q;
- embedded/nonperiodic geometries: use the natural real-space response matrix.

### 7.3 Real-space convergence

In metals the range grows as the imaginary energy offset becomes small. Therefore production calculations must expose and test:

- R-shell/Rmax cutoff;
- real-space norm tail of `chi0(R,w)`;
- convergence of q-space spectra under Rmax;
- contour/broadening dependence.

Converge the response, not merely the one-electron GF.

### 7.4 Important amortization advantage

For a periodic q path, `chi0(R,w)` is independent of q. Once built at a given frequency, the same real-space response can be Fourier transformed to all requested q points cheaply. This can make the native RS route particularly attractive for dense q paths when the R cutoff is controlled.

---

## 8. Energy integration strategy

Implement the integration in stages so there is always a slow transparent reference.

### 8.1 Stage 1: direct near-real-axis reference

Use retarded/advanced GF values with explicit small eta and a sufficiently dense real-energy quadrature. This path is not expected to be fastest. Its purpose is:

- correctness;
- debugging signs/factors;
- validating contour paths;
- backend equivalence on small fixtures.

### 8.2 Stage 2: mixed contour decomposition — recommended first production method

Adapt the Lounis/Mills decomposition of the GF convolution:

- analytic same-half-plane contribution -> deform onto a regular complex contour;
- nonanalytic contribution -> integrate only over the small near-Fermi interval controlled by omega.

Do this using the full RS-LMTO basis; do **not** import the d-only/energy-independent-wavefunction approximation unless explicitly offered as a separate approximate mode.

Advantages:

- expensive GF evaluations can be moved away from poles;
- only a short energy window needs near-real-axis resolution;
- no global analytic continuation of the final susceptibility is required.

### 8.3 Stage 3: Halle complex-frequency + analytic continuation — optional advanced mode

The Halle implementation evaluates susceptibility away from the real axis using finite complex contours/Matsubara terms and continues it back to the real axis, commonly with rational/Padé-like approximations.

Implement this only after direct and mixed-contour paths are trusted because analytic continuation is intrinsically ill-conditioned.

Requirements if implemented:

- never be the only available production route;
- report the complex sampling grid and continuation method;
- estimate stability by varying sampling points/order and comparing with direct/mixed-contour results;
- reject continuations that introduce noncausal poles or violate response symmetries/sum rules;
- distinguish numerical imaginary offset from physical linewidth.

---

## 9. Kernel, Ward identity, and Goldstone policy

This is the highest-priority physics issue.

### 9.1 Exact structural identity

For the static transverse response of a collinear SOC-free magnet, the ground-state magnetization vector/function is the Goldstone eigenvector. In operator notation the core relation is equivalent to

`B_xc = chi_KS^{-1}(q=0,w=0) m_GS = K_xc m_GS`

and therefore

`D = I - chi_KS(0,0) K_xc`

must have a zero eigenvalue associated with `m_GS`.

The Lounis real-space sum rule is the same consistency statement in a projected representation: the static KS response acting on the effective exchange field must reproduce the magnetization.

### 9.2 RS-LMTO policy

1. Derive the response kernel/perturbation from the same ground-state Hamiltonian representation.
2. Check the uncorrected Ward residual and show that it converges systematically with k/R/energy integration and response-basis resolution.
3. If a site-only response basis cannot make the residual converge, test a richer local orbital response basis before applying a correction.
4. Provide controlled optional Ward-restoring schemes only as diagnostics/numerical cleanup:
   - Lounis-style sum-rule kernel reconstruction in the chosen response basis;
   - Halle-style projection of only the small spurious Goldstone eigenvalue.
5. Never apply an empirical global lambda or arbitrary energy shift.
6. Always print/report correction magnitude. A large correction is a failure of the underlying discretization, not a successful run.

Do not hard-code a universal acceptable percentage before evidence exists. Establish tolerances from convergence studies on Fe/Ni and require the residual/correction to shrink under systematic numerical refinement.

### 9.3 Spectral sum rule

The integrated transverse loss spectral weight is tied to the ground-state magnetization (sign/factor depends on the circular-channel convention). Implement this as a diagnostic. It is an excellent detector of missing factors of two, incorrect frequency sign, or inconsistent projection weights.

---

## 10. Dyson solver and response observables

The common Dyson layer should accept a complex KS response matrix and kernel in the same response basis.

Use numerically stable linear solves/factorizations rather than explicit matrix inversion wherever possible.

Track conditioning of

`D(q,w) = I - chi_KS(q,w) K`.

Near collective modes D becomes nearly singular by physics; the solver must distinguish this from numerical breakdown.

### 10.1 Loss matrix

Define and document the retarded loss matrix convention, e.g. the anti-Hermitian part with sign chosen so physical spectral weight is positive in the selected channel.

Diagonalize the loss matrix to obtain:

- mode spectral weights;
- mode eigenvectors/site character;
- acoustic/optical branches;
- chirality/channel character;
- peak position and width where a quasiparticle-like peak exists.

Do not force a Lorentzian fit when the mode has merged strongly with the Stoner continuum. Report it as overdamped/non-quasiparticle instead.

### 10.2 Stoner continuum

Keep both KS and enhanced loss spectra. The KS spectrum is the Stoner/electron-hole baseline; the enhanced spectrum shows the collective modes and their Landau damping.

Optional later observable: q-resolved “Landau maps” similar to the Halle analysis.

---

## 11. Validation hierarchy

### 11.1 Analytic/software fixtures

Before material calculations:

- one/two-level analytic GF and chi0 fixtures;
- direct inverse versus Lehmann G(k,z);
- Fourier consistency G(R,z) <-> G(k,z) for a periodic toy/electronic fixture (validation only; not a runtime requirement);
- direct energy integral versus mixed contour;
- causal/retarded symmetry identities;
- spectral sum rule;
- Ward/Goldstone static identity.

### 11.2 bcc Fe and fcc Ni — release gate

These are the primary production validation systems.

Required evidence:

- all three chi0 backends agree at selected q,w points under converged settings;
- SOC=0 Goldstone mode tends to zero without empirical tuning;
- small-q acoustic dispersion is quadratic for the FM;
- stiffness converges with k mesh/R cutoff/energy integration;
- TD-DFT low-q stiffness is consistent with the code’s independent LKAG and frozen-magnon/GBT references within the expected methodological/numerical uncertainty;
- Stoner continuum and linewidth trends are physically sensible;
- results are stable against eta/contour parameters.

For Ni, explicitly check the previously problematic near-Gamma behavior and rule out path-coordinate/backfolding mistakes.

### 11.3 FeRh — FM/AF transition test

Use the same material in FM and AF phases to test:

- AF transverse-channel structure;
- acoustic and optical modes;
- Landau damping versus Stoner continuum;
- induced-moment/multi-sublattice character where representable in the current model.

### 11.4 CrMnSb — chirality test

Use the compensated half-metallic ferrimagnet as a targeted test that `chi^{+-}` and `chi^{-+}` are kept distinct. The two chiral branches should exhibit qualitatively different access to the Stoner continuum and hence different damping.

### 11.5 Longitudinal Fe/Ni/FeSe

Only after the transverse release gate:

- coupled `(0,z)` response;
- Hartree screening in the charge channel;
- longitudinal collective fluctuations;
- energies/FWHM and their momentum dependence;
- no false Goldstone expectation in the longitudinal channel.

---

## 12. Performance design

### 12.1 Optimize only validated kernels

Do not begin with a GPU port. First establish backend equivalence and profile real workloads.

### 12.2 Parallel decomposition

Potential MPI levels:

- q points;
- frequency points;
- k points for eigenpair/k-GF paths;
- R blocks and/or contour-energy nodes for native RS GF.

Avoid duplicating large GF caches across MPI ranks without measuring the tradeoff.

The Halle experience shows KS-response evaluations at distinct q/frequency points are highly parallel. The native RS route adds a different opportunity: build `chi0(R,w)` once and transform to many q points, so q-level parallelism should not destroy that reuse.

### 12.3 BLAS/batching

- Contract orbital/site blocks with GEMM-like kernels where possible.
- Batch contour energies and q points.
- Avoid repeated small allocations in inner loops.
- Reuse projectors/vertices and static ground-state quantities.
- Cache one-electron GFs only when memory-bandwidth tradeoffs are favorable.

### 12.4 GPU policy

Only port kernels with sufficient arithmetic intensity/batch size. The recent broader RS-LMTO GPU experience shows that launch/setup/transfer overhead can dominate small matrix problems. A GPU path must therefore be justified by measured crossover data, not by feature parity.

---

## 13. Input/API policy

Preserve existing inputs. Add backend selection without silently changing historical semantics.

Suggested conceptual options (adapt names to current parser):

```text
chi0_backend = eigenpairs | gf_k_lehmann | gf_r_native

gf_integration = direct | mixed_contour | complex_frequency_ac
response_projection = site | site_orbital   # latter may begin as diagnostic

goldstone_policy = diagnose | sum_rule | projected
```

Rules:

- `diagnose` should be the default scientific mode during development.
- `sum_rule`/`projected` must report exactly what was changed and by how much.
- unsupported ground-state features fail early with a precise reason.
- provenance output must record backend, integration scheme, grid/cutoff, eta/contour, and correction policy.

---

## 14. Sternheimer/Liouville-Lanczos future track

After the three Dyson/GF backends are validated, perform a design study for a fourth independent engine.

### 14.1 Savrasov/Sternheimer

Solve resonant and antiresonant first-order equations for the response orbitals under a monochromatic magnetic perturbation, self-consistently updating the induced Hartree/xc fields. This directly produces induced charge/magnetization without explicit empty-state sums.

RS-LMTO-specific research questions:

- how to formulate the Sternheimer equation in the orthogonalized LMTO representation used by modern RS-LMTO;
- how to handle metallic occupation/Fermi-surface terms;
- preconditioning/linear solver appropriate to the LMTO Hamiltonian;
- reuse of recursion/GF technology as a linear solver/preconditioner;
- exact relation to existing noncollinear Hamiltonian variations;
- later inclusion of Hubbard response.

### 14.2 Liouville-Lanczos

Modern TDDFPT shows that a Liouvillian Lanczos approach can obtain the spectrum over many frequencies after a single iterative construction. This is potentially attractive once a stable first-order Hamiltonian action is available.

Do not implement this before the design prompt demonstrates that it can share the response vertices/kernel and validation suite with the existing TD-DFT subsystem.

---

## 15. Relativistic/noncollinear future track

SOC invalidates the simple transverse decoupling and ordinary Goldstone condition. Before enabling:

- derive the full spin/charge response structure in the noncollinear spinor basis;
- follow the relativistic sum-rule logic of dos Santos Dias et al.;
- include anisotropy/external-field contributions consistently;
- validate the SOC-induced magnon gap and linewidth against independent static anisotropy and known test systems;
- retain charge coupling where symmetry allows it.

Until then, keep strict guards.

---

## 16. Definition of done for the core campaign

The core TD-DFT campaign is complete only when:

- [ ] the response conventions are documented and unit-tested;
- [ ] the transverse kernel/vertices satisfy the static Ward identity under convergence without empirical tuning;
- [ ] the existing eigenpair backend is stable and reproducible;
- [x] `G(k,z)` Lehmann GF is independently validated;
- [x] a K-space GF chi0 backend exists;
- [ ] a native real-space GF chi0 backend exists without mandatory `G(k)` conversion;
- [ ] direct and contour integration paths agree on controlled fixtures;
- [ ] all three chi0 backends agree for Fe/Ni within converged numerical uncertainty;
- [ ] q=0/near-q=0 Goldstone and q^2 dispersion are correct for SOC=0 FM cases;
- [ ] stiffness is consistent with independent LKAG/frozen-magnon/GBT references;
- [ ] loss-matrix mode analysis yields stable energies and linewidths;
- [ ] all unsupported response combinations fail explicitly;
- [ ] performance evidence identifies the production backend regimes rather than assuming one backend is always fastest;
- [ ] documentation records the physics approximations and validation evidence.

Only then should the code describe the transverse TD-DFT path as production-validated.
