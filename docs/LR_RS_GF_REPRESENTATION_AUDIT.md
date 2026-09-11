# LR-REP-00 — Native RS Green-function representation audit

## Disposition

The representation algebra used by the live RS Green-function consumers is
certified for the orthogonal, collinear, no-SOC coefficient-space baseline:

```text
orthogonal coefficient GF
    -> sqrt(Delta) endpoint auxiliary GF
    -> screened auxiliary/path-operator GF
```

The reciprocal/RS coefficient-space bridge and both representation conversions
are covered by `UnitLrRsGfRepresentation`. The accepted block terminator and
Chebyshev reconstruction are not changed or re-opened by this audit.

The native RS *response* backend remains **BLOCKED**. The reason is precise:
`green%auxiliary_gij` is an LMTO endpoint scaling operation, not the radial/Pauli
augmentation used by the reciprocal TD-DFT response. A future RS response route
must add and test that endpoint adapter before it can claim the same physical
response object as LR-06/LR-GF-02.

This is a representation audit, not a susceptibility implementation.

## Audit identity

| item | value |
|---|---|
| branch | `fable_v4` |
| starting HEAD | `68b7fc9736a7c1460b2b3b377d3a28ce9ab82559` |
| audit date | 2026-09-11 |
| compiler/build | GNU Fortran 13.3.0, existing `build`, `RUN_UNIT_TESTS=ON` |
| production scope | orthogonal `ham_only`, second-order/HOH, `sp`, collinear, no SOC/additive operator |
| prescribed commit | `docs: certify real-space Green-function representations` |

The worktree already contained the untracked post-LR03 prompt pack at preflight;
those prompt files are not part of this implementation except for the completion
checklist added to task 10.

## Live representation graph

### Reciprocal coefficient-space GF

The reciprocal Lehmann backend fills the same `green%gij/gji` arrays used by the
RS consumers. Its block is

\[
 G^c_{ij}(z)=\frac{1}{N_k}\sum_{k,n}
 e^{+i k\cdot(R_i-R_j)}
 \frac{c_{i,nk}c_{j,nk}^{\dagger}}{z-\epsilon_{nk}}.
\]

The live path is:

```text
reciprocal%fill_green
  -> fill_green_lehmann
     -> lehmann_pair_block
        -> green%gij / green%gji
```

`source/reciprocal_green.f90` documents this as an orthogonal LMTO
coefficient-space block. It is not a radial spatial Green function and is not
silently a screened/path-operator representation.

### Native RS coefficient-space GF

The native path is:

```text
recursion%recur_b_ij or recursion%chebyshev_recur_ij
  -> green%calculate_intersite_gf
     -> four phase combinations
        -> green%gij / green%gji
```

For a pair `(i,j)`, the four seeds are the normalized block combinations

\[
 (a,b)=\frac{1}{\sqrt 2}(1,1),\quad
 \frac{1}{\sqrt 2}(1,-1),\quad
 \frac{1}{\sqrt 2}(1,+i),\quad
 \frac{1}{\sqrt 2}(1,-i).
\]

The live reconstruction in `source/green_block.f90` is

\[
 G^c_{ij}=\frac12\left[g_{++}-g_{+-}
       +\frac{1}{i}(g_{+i}-g_{-i})\right].
\]

The RS result is therefore the same coefficient-space object as the reciprocal
block when the same resolvent is used. Native block recursion and Chebyshev add
their documented finite-chain/finite-polynomial approximations; that numerical
convergence is outside this representation audit.

## Basis and ordering contract

The live coefficient block is site-local and has shape `(nb,nb)`. Within a site,
the ordering is:

```text
spin 1: l=0,m=1; l=1,m=1..3; ...
spin 2: the same orbital sequence
```

The flattened orbital index is `l*l + m`, and the spin offset is
`(lmax+1)**2`. For the certified `sp`/`spd` normal basis this is the site-major,
spin-blocked ordering already used by the reciprocal Hamiltonian and radial
augmentation code. `gij` has endpoint `i` on the left and endpoint `j` on the
right; `gji` reverses both the site order and the block orientation.

The inverse-Bloch phase is also fixed by the live code:

\[
 e^{+i k\cdot(R_i-R_j)},
\]

with fractional `k` and bond coordinates and the `2*pi` factor applied by the
implementation.

## Transformation 1 — coefficient GF to auxiliary endpoint GF

### Live operation

`green%auxiliary_gij` constructs diagonal endpoint matrices

\[
 D_a=\operatorname{diag}_{lms}\sqrt{\Delta_{a l\sigma}}
\]

from `potential%dele`, whose source-level contract is `sqrt(Delta)`, and applies

\[
 g^{\mathrm{aux}}_{ij}(z)=D_i\,G^{\mathrm{in}}_{ij}(z)\,D_j.
\]

`dele` and the endpoint ordering are energy independent. The operation is
performed independently for every stored energy channel. It is an endpoint
scaling identity; it does not calculate radial functions, spherical harmonics,
small components, GFAC factors, or Pauli response channels.

### Representation label

The routine’s input comments call `green_ij` a physical site-resolved GF, but the
live reciprocal and RS fillers populate `green%gij/gji` with the orthogonal
coefficient resolvent documented above. Therefore the safe live label is:

```text
auxiliary_gij(Gc) = coefficient-space endpoint-scaled auxiliary candidate
```

It may be called a physical LMTO auxiliary GF only when its input has first been
certified as the physical LMTO GF. This audit does not promote the raw `gij`
array across that boundary.

## Transformation 2 — screened auxiliary/path-operator representation

### Potential-function identity

For one diagonal channel, with `alpha` equal to `screening_in` and `beta` equal
to `screening_out`, the live `transform_pmatrix` operation is

\[
 P^\beta(z)=
 \frac{P^\alpha(z)}{1+[\alpha-\beta]P^\alpha(z)}.
\]

The live `p_matrix` contract is

\[
 P^\alpha_l(z)=\frac{z-C_l}{\Delta_l},
\]

stored in the diagonal site/orbital/spin matrix with one value per energy
channel. Screening constants are real, channel dependent, and energy
independent.

### Auxiliary-GF identity

Define the diagonal endpoint ratios

\[
 A_i^{\alpha\to\beta}(z)=P_i^\alpha(z)[P_i^\beta(z)]^{-1}.
\]

The live `transform_auxiliary_gij` formula is

\[
 g^\beta_{ij}(z)=
 A_i^{\alpha\to\beta}(z)g^\alpha_{ij}(z)A_j^{\alpha\to\beta}(z)
 +\delta_{ij}[\beta_i-\alpha_i]A_i^{\alpha\to\beta}(z).
\]

The additive term is a full diagonal matrix in the onsite block and is absent
for `i /= j`. The source currently assumes diagonal/spherical `P` matrices and
does not validate the denominator or matrix shape; the certified production
scope supplies those conditions.

### Exact inverse

The inverse is the same operation with input and output exchanged:

\[
 g^\alpha_{ij}(z)=
 A_i^{\beta\to\alpha}(z)g^\beta_{ij}(z)A_j^{\beta\to\alpha}(z)
 +\delta_{ij}[\alpha_i-\beta_i]A_i^{\beta\to\alpha}(z).
\]

The two endpoint factors and the onsite additive terms cancel exactly. This is
not an offsite-only identity; omitting the onsite term fails the inverse.

## Physical/radial/Pauli boundary

The certified reciprocal response uses the separate radial endpoint operator

\[
 \Psi_{n,a}(r)\simeq
 [\Phi_a(r)+\dot\Phi_a(r)h_{\gamma,a}]c_{a,n},
\]

with the scalar-relativistic/no-SOC Pauli projection and angular Gaunt map. The
live `lr_gf_susceptibility` route carries this endpoint dependence through its
energy-moment vertex tensor (`p,q=0,1`); it does not call
`green%auxiliary_gij`.

The current RS `gij/gji` fillers have no equivalent call that maps both GF
endpoints through `lmto_radial_basis`, the Pauli operator, and the LR-04 radial
response coordinate space. In particular:

```text
RS coefficient GF -> dele endpoint scaling       certified below
RS coefficient GF -> radial/Pauli response GF    not implemented
```

This is the blocking condition for RSGF-01. It prevents a future RS response
implementation from comparing a screened coefficient/path matrix directly with
the reciprocal physical/Pauli response matrix.

## Required tests and results

Test source: `tests/unit/test_lr_rs_gf_representation.f90`.

| oracle | result | claim kept separate |
|---|---|---|
| `sqrt(Delta)` endpoint scaling against explicit diagonal matrices | PASS, `0.0` max error | algebraic endpoint scaling only |
| screened offsite `alpha -> beta -> alpha` | PASS, `4.46e-16` | independent numerical inverse check |
| screened onsite `alpha -> beta -> alpha` | PASS, `6.72e-16` | includes additive onsite term |
| onsite negative control with additive term deliberately omitted | PASS, wrong oracle differs by `2.94e-1` | protects onsite formula from offsite-only regression |
| reciprocal Lehmann block versus four-phase RS projector reconstruction | PASS, `1.67e-15` | coefficient-space bridge, not terminator validation |
| transformed reciprocal/RS bridge after endpoint and screening conversions | PASS, `6.81e-16` | same representation labels are compared |
| raw coefficient onsite high-energy normalization | PASS, `1.90e-4` at `|z|~1000`, improved from `2.36e-3` at `|z|~80` | resolvent asymptotic |
| `sqrt(Delta)` onsite high-energy normalization | PASS, `1.11e-4` at `|z|~1000`, improved from `1.38e-3` | endpoint-scaled asymptotic |
| `sqrt(Delta)` offsite high-energy decay | PASS, `1.12e-4` for `|z|G_{ij}^{aux}` | offsite resolvent asymptotic |
| screened high-energy onsite identity | PASS, `4.41e-14` | additive-term cancellation |

The bridge test uses a two-site periodic Hermitian fixture. The reciprocal block
is produced by the independent Lehmann kernel; the native-RS side uses direct
complex resolvent blocks and the same four seed combinations as
`green_block.f90`. This certifies phase, orientation, and representation
conversion algebra without treating a finite recursion terminator as an exact
inverse.

## Precision correction made during the audit

The live representation routines used two-argument `cmplx(real,0)` calls. With
GNU Fortran this selects the default single-precision complex kind before
assignment to `complex(rp)`, producing an endpoint-scaling error of about
`2.9e-8` in the deterministic oracle. The representation paths now use the
explicit `cmplx(real,0,rp)` form in:

- `source/green.f90` (`auxiliary_gij` and `transform_auxiliary_gij`);
- `source/symbolic_atom.f90` (`p_matrix` and `transform_pmatrix`).

The correction changes precision only; it does not change the LMTO identities,
screening convention, terminator, or response implementation.

## Terminator policy

The accepted block terminator remains outside this audit. No terminator branch,
tail formula, recursion depth, or Chebyshev damping convention was changed.
The bridge test intentionally uses an exact finite fixture to isolate
representation and phase algebra from terminator convergence.

## Claims gate

| claim | status |
|---|---|
| orthogonal reciprocal/RS coefficient-space labels agree | **certified for the stated fixture and existing LR-GF-01 scope** |
| `sqrt(Delta)` endpoint conversion is algebraically implemented and tested | **certified** |
| screened/path-operator conversion, onsite term, and inverse | **certified for diagonal ASA P matrices** |
| native RS block/Chebyshev numerical convergence | **existing approximation contract; not re-audited here** |
| native RS GF already equals the reciprocal radial/Pauli response object | **not claimed** |
| native RS TD-DFT response backend | **BLOCKED pending the radial/Pauli GF endpoint adapter** |

## Verification

```text
cmake -S . -B build -DRUN_UNIT_TESTS=ON -DENABLE_SPGLIB=OFF
cmake --build build --target UnitLrRsGfRepresentation -j2
ctest --test-dir build --output-on-failure -R '^UnitLrRsGfRepresentation$'
```

Result:

```text
100% tests passed, 0 tests failed out of 1
```

## Completion checklist

- [x] live `auxiliary_gij` path audited;
- [x] live `transform_auxiliary_gij` path audited;
- [x] `sqrt(Delta)` endpoint factors and `dele` ordering documented;
- [x] screening-constant and potential-function identities derived;
- [x] onsite versus offsite formulas documented;
- [x] exact inverse documented and round-trip tested;
- [x] onsite negative control requires the additive term;
- [x] reciprocal/RS coefficient-space bridge tested;
- [x] transformed bridge compares matching representation labels;
- [x] high-energy behavior tested for raw, endpoint-scaled, and screened paths;
- [x] accepted terminator left unchanged;
- [x] precision loss in the audited conversion paths removed;
- [x] affected native RS response backend explicitly declared `BLOCKED`;
- [x] `docs/LR_RS_GF_REPRESENTATION_AUDIT.md` created;
- [x] focused unit test registered with CMake.

The audit therefore certifies the representation algebra needed for a future
RS backend, while preserving the separate blocker on radial/Pauli response
endpoint augmentation.

## Formal reference

The formal target is the screened-LMTO/KKR representation framework described in
[Andersen et al., *Third-Generation TB-LMTO*](https://arxiv.org/abs/cond-mat/9804166),
together with the Eq. (3.56)/(3.57) convention cited by the live source comments
to Turek et al., *Electronic Structure of Disordered Alloys, Surfaces and
Interfaces*.

## Subsequent RSGF-CLOSE-01 closure

The `BLOCKED` native-response statement above was correct at this audit’s
starting point (`68b7fc9736a7c1460b2b3b377d3a28ce9ab82559`): the physical
radial/Pauli endpoint adapter and native response service did not yet exist.
Subsequent tasks added and tested those seams:

- RSGF-00 supplies the four-branch two-endpoint Pauli augmentation, including
  onsite/offsite resolvent contact terms;
- RSGF-01R calls that adapter for directed native GF blocks, assembles the full
  real-space bubble and q phase, and returns the LR-04 canonical response;
- the focused finite/provider tests compare the complete response against the
  independent spectral and reciprocal-GF references.

Therefore the historical blocker is **superseded for R0–R2** within the
documented finite/provider baseline. Production-driver registration remains R3
and is intentionally pending TDRUN-02; Fe/Ni material validation remains R4
and is intentionally pending TDVAL-01R. The authoritative current ledger is
[`RSGF_CAPABILITY_CLOSURE.md`](RSGF_CAPABILITY_CLOSURE.md).
