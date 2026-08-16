# VAL-17 — GBT harmonic cone-angle and Goldstone behavior

Date: 2026-08-16
Status: A/B physical gates pass; C is measured but the fine-mesh stiffness is not yet converged. GBT maturity is therefore not promoted.

## Scope and route

This campaign uses the current VAL-16 GBT/supercell foundation, scalar-relativistic nsp=3 decks, strux_backend='strux_lib', and the reciprocal ham_only force-theorem route. The independent FeCo comparison uses the same branch_mode='auto' deck through real-space recursion. No SOC, diagonal shift, empirical angle correction, or fitted rescaling is used.

The reproducible campaign is tests/validation/val17_gbt_harmonic_goldstone.py, registered as Val17GbtHarmonicGoldstone in CMakeLists.txt.

## Part A — bcc-Fe cone angle

The production sweep uses q=(0,0,0.05), nk=12^3, and theta_ss=5,10,15,20 degrees. The diagnostic file records both the historical fixed-Gamma subtraction and the same-q collinear reference:

| theta | E(q,theta) | E(0,theta) | DeltaE raw | sin2(theta) | M | E(q,0) | DeltaE gauge | DeltaE gauge/sin2 | omega |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | -1.998855050127 | -1.998829585938 | -2.5464189e-5 | 7.5961235e-3 | 2.089530404 | -1.998856756118 | 1.7059910e-6 | 2.25e-4 | 4.29928e-4 |
| 10 | -1.998849958802 | -1.998829586033 | -2.0372769e-5 | 3.01536896e-2 | 2.089530404 | -1.998856756213 | 6.7974111e-6 | 2.25e-4 | 4.31533e-4 |
| 15 | -1.998841573890 | -1.998829586033 | -1.1987857e-5 | 6.6987298e-2 | 2.089530404 | -1.998856756213 | 1.5182323e-5 | 2.27e-4 | 4.33867e-4 |
| 20 | -1.998830091404 | -1.998829586055 | -5.0534915e-7 | 1.1697778e-1 | 2.089530404 | -1.998856756235 | 2.6664831e-5 | 2.28e-4 | 4.36362e-4 |

The raw DeltaE/sin2 control is strongly angle-dependent, reproducing the historical failure. The physical same-q subtraction has a 1.49% omega spread across the sweep and a stable DeltaE/sin2 slope. The Gamma row is zero to 2.3e-16 Ry in the diagnostic energy.

Occupation and integration diagnostics were recorded independently. The target electron count is 8.0 for every theta; canonical occupation residuals are below 2.7e-11, raw k-point weight sum is 1.0, and the moment is constant to the printed precision. The q=0 EF is -0.0670314745 Ry; finite-q EF values drift from about -0.0679538495 to -0.0678141765 Ry across the cone sweep while the electron count remains exact. Thus the Fermi search is not losing charge, but its finite-mesh state is not identical between q and Gamma.

The requested deck state was q_symmetry_policy='full_bz', use_time_reversal=true, and nk=12^3. The nonzero-q route rebuilt the full 1728-point chemical BZ; occupation logs report normalized weights. The q=0 GBT common-frame check is max|G-I|=0. These observations rule out an under-reduced BZ and a frame mismatch as the cause.

### Root cause and code change

At theta=0 the finite-q GBT phase is a pure gauge. In the continuum, E(q,0)=E(0,theta), but the finite k mesh does not preserve that identity when q is not a mesh translation. The old E(q,theta)-E(0,theta) subtraction therefore retained a q-only integration offset and amplified it by 1/sin2(theta). The reciprocal MFT path now evaluates E(q,0) on the same fixed potential and same k-point set, then uses E(q,theta)-E(q,0). The raw and gauge-subtracted values are both written to frozen_magnon_diagnostics.dat for auditability. The real-space subtraction and output contract are unchanged.

This is a physical gauge-reference correction, not a theta-dependent empirical correction. No potential is reused across different q values without the existing reciprocal mesh policy; within each q, the force-theorem potential and moment normalization remain fixed.

## Part B — FeCo Goldstone behavior

The converged reciprocal FeCo run uses nstep=20, nk=8^3, N=17, and the q list 0, 0.02, 0.05. Its active moments are Fe 2.78814385 and Co 1.75514900 uB. Canonical occupation logs keep N=17, |dN|<8.2e-11, and normalized weights.

At Gamma:

- reciprocal acoustic branch: omega_ac=6.48e-15 Ry;
- reciprocal acoustic eigenvector real components: -0.78338,-0.62154, hence in phase;
- reciprocal optical branch: 3.8963e-2 Ry;
- q=0 GBT common-frame check: 0 to machine precision.

The independent real-space run also has a zero acoustic candidate, -4.16e-17 Ry, with same-phase Fe/Co components. Branch numbering differs between routes, so the acoustic branch is identified by the smallest absolute Gamma eigenvalue, not by a fixed branch index. The real-space finite-q values are not used as a stiffness gate because the finite-cluster route has its own long-wavelength limitations.

No diagonal shift or sum-rule tuning was applied. The reciprocal auto-branch Hessian uses the same-q collinear energy as its q-dependent force-theorem reference; at q=0 this is the ordinary collinear reference and the Goldstone result follows from the energy surface.

## Part C — small q and mesh convergence

The reciprocal FeCo acoustic values after the same-q subtraction are:

| mesh | omega/q2 at q=.02 | omega/q2 at q=.05 | interpretation |
|---:|---:|---:|---|
| 8^3 | 9.06e-3 | 1.01e-1 | visibly mesh-sensitive at the smallest q |
| 12^3 | 5.25e-2 | 6.20e-2 | internally close to quadratic |
| 16^3 | 1.09e-1 | 9.88e-2 | internally close to quadratic |

The 12-to-16 change in the inferred stiffness is about 55%, so the material stiffness is not promoted as converged. The finer meshes show the expected small-q quadratic trend, while the 8^3 result demonstrates why a single low-resolution point is insufficient. This is an open convergence axis, not a reason to tune the Gamma value.

For the independent bcc-Fe mesh probe at q=.5, the eband values remain -1.98668388, -1.98714329, and -1.98785909 Ry for nk=8,12,16. The refinement steps are 4.59e-4 and 7.16e-4 Ry, respectively: the mild non-monotonicity remains visible. At q=1.0 the corresponding steps shrink, 5.72e-4 to 5.25e-4 Ry. The campaign reports this behavior rather than hiding it behind a monotonicity gate.

## Checklist

- [x] theta sweep reproduced
- [x] DeltaE~sin2(theta) tested directly, including the raw fixed-Gamma control
- [x] occupations, EF, moments, k-point set/weights, and symmetry state recorded
- [x] frame conventions, q=0 GBT identity, fixed potential, and q-only offset checked
- [x] root cause of the angle dependence identified
- [x] no empirical theta rescaling added
- [x] k-space FeCo Gamma remeasured
- [x] acoustic eigenvector checked
- [x] independent RS Gamma comparison recorded
- [x] small-q q2 behavior assessed
- [x] k-mesh convergence and non-monotonicity assessed
- [ ] GBT maturity promoted — held because the 12-to-16 small-q stiffness is not converged
