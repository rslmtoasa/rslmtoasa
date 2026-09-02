# TDDFT-07: native real-space Green-function chi0

TDDFT-07 adds the native real-space route

```text
G_ab(R,z), G_ba(-R,z)  ->  chi0_ab(R,omega)  ->  chi0_ab(q,omega)
```

The response layer consumes the mature RS-LMTO `green%gij` and `green%gji`
blocks directly. It does not call `reciprocal%fill_green`, construct a
`G(k,z)` intermediate, or require k+q eigenpairs for the bare response.
The concrete provider is `tddft_native_realspace_gf_provider` in
`source/tddft_chi0_realspace.f90`; its `initialize_from_green` bridge copies
the native blocks and lattice-pair geometry into a response-owned source.

## Exact dynamic contraction

All energies and frequencies are in Rydberg. For response channels A and B,
the implementation uses the same unhalved site/circular operators as the
eigenpair and K-space GF backends. For each pair `(a,b,R)` it evaluates

```text
chi0_AB(R,w) = -1/(2*pi*i) integral dE f(E) {
    Tr[A G^R_ab(R,E+w) B (G^R_ba(-R,E) - G^A_ba(-R,E))]
  + Tr[A (G^R_ab(R,E) - G^A_ab(R,E)) B G^A_ba(-R,E-w)]
}
```

`G^A_ba(-R,E) = transpose(conjg(G^R_ab(R,E)))` and
`G^A_ab(R,E) = transpose(conjg(G^R_ba(-R,E)))`. Keeping these two advanced
blocks distinct is important for non-Hermitian directional pair data. The
native implementation linearly interpolates the supplied retarded blocks and
integrates them with the trapezoidal rule on their native energy grid. The
energy interval for a frequency is clipped so both `E+w` and `E-w` remain on
the supplied grid. `green_eta=0` is reported as `eta/2`; it does not regenerate
or alter the supplied native Green data.

The circular vertex direction is retained: a left `PLUS` channel and a right
`MINUS` channel use the `G_ab` / `G_ba` ordering above. No empirical sign,
factor, or spectral shift is applied. The output spectral convention remains
`-Im chi0 / pi`.

## Exact static contraction

The native provider has a separate `evaluate_static_realspace` operation and
the common backend exposes it through `evaluate_static_grid`. It does not call
the dynamic routine at `omega = 0` with finite `eta`. For a directed pair its
static definition is

```text
T^R(E) = Tr[A G^R_ab(R,E) B G^R_ba(-R,E)]
T^A(E) = Tr[A G^A_ab(R,E) B G^A_ba(-R,E)]
chi0_AB(R,0) = -1/(2*pi*i) integral dE f(E) [T^R(E) - T^A(E)]
```

The spectral representation of this retarded/advanced contour identity is the
divided difference `(f_n-f_m)/(e_n-e_m)`, including its equal-energy Fermi
derivative. The RGF implementation evaluates this identity directly on the
native real-axis energy grid, keeps the `G_ab`/`G_ba` orientation, then
Fourier-transforms `chi0(R,0)`. Static metadata records zero response
broadening and the static contour implementation explicitly.

## Susceptibility Fourier transform

For periodic bulk response the transform is

```text
chi0_AB(q,w) = sum_R exp(-i 2*pi q.R) chi0_AB(R,w),
```

where both `q` and `R` use the project fractional-coordinate convention.
The code transforms `chi0(R,w)`, not the one-electron Green function. A
separate `fourier_transform_realspace_green` routine exists only as a
validation utility for checking phase conventions against a trusted K-space
Green function.

The representation selector is explicit:

* `bulk`: transform over all configured axes;
* `film`: transform only over `realspace_fourier_axes` (normally the two
  in-plane axes), retaining the out-of-plane real-space coordinate;
* `finite` or `local`: use the direct real-space response, with unit phase;
  q is accepted only as a common batch envelope and has no physical meaning.

The generic provider API accepts arbitrary pair tables, so the same response
contraction can be used for finite, impurity, embedded, and future film
sources. The validated collinear, SOC-free transverse bulk path supports both
the direct real-axis reference and `gf_integration='mixed_contour'` when a
native complex-energy source is attached. Native RGF also supplies the exact
static bare response needed by the transverse Ward diagnostic. Unsupported
longitudinal/full-response combinations and Xi/Dyson enhancement remain
capability-gated independently of the contour implementation.

## R cutoff and tail diagnostics

`&tddft` exposes:

```fortran
realspace_rmax          = huge(1.0)
realspace_tail_tolerance = 1.0e-3
realspace_representation = 'bulk'
realspace_fourier_axes   = 1, 2, 3
```

The cutoff is measured with the supplied lattice metric. The provider reports
the number of input, retained, and omitted pair points; effective cutoff;
distinct shell count; retained outer-shell norm; omitted-tail norm; tail
ratio; tolerance; and whether convergence was assessed. The tail norm is
computed from `chi0(R,w)` pair contributions, not from the magnitude of
`G(R,z)`. When no pair is omitted, the status is “all supplied real-space
pairs retained” and convergence is reported as not assessed, which avoids
turning a finite source list into an unjustified convergence claim.

These fields are propagated into `tddft_chi0_metadata` and written by
`write_chi_ks_text`, making every q/frequency product self-describing.

## Reuse and validation

The provider builds `chi0(R,omega)` once per frequency batch. It then performs
one susceptibility transform over all requested q points. The unit fixture
`UnitTddftRealspaceGF` verifies:

* native response evaluation without a K-space source;
* separate forward/reverse pair geometry and advanced-block conventions;
* bulk phase sign and finite/local direct representation;
* Rmax omission and chi0 tail metrics;
* equality of a converged one-site periodic native bubble and the K-space GF
  bubble on the same real-axis grid;
* exact static native-RGF agreement with the eigenpair divided-difference
  reference at q=0 and finite q, including the uncorrected raw Ward residual;
* one real-space build reused for a three-q batch.

The K-space Green transform test in that unit is deliberately isolated as a
phase-validation utility. It is not part of the native production path.
