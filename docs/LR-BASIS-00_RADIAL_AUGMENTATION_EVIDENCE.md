# LR-BASIS-00: radial augmentation evidence

Status: baseline certified for the scope below. This report separates code-derived facts, test-confirmed facts, formal LMTO interpretation, and deferred extensions.

## Audit identity and scope

The source audit was started on branch `fable_v4` at commit `497ff422d38c93c6bb82f5c077cc868de54f046a6` (`Certify one-electron Green-function contracts for linear response`). The working-tree additions are the focused tests, test-only seams in `source/self.f90`, and `lmto_radial_augmentation_mod`. The legacy radial routine bodies remain in `source/self.f90`.

| condition | status |
|---|---|
| reciprocal `ham_only` | supported |
| orthogonal second-order / HOH | supported |
| collinear two-channel radial basis | supported |
| SOC | guarded and deferred |
| Hubbard or other additive operator | guarded and deferred |
| generalized overlap | guarded and deferred |
| `sp` / `spd` | supported |
| `spdf` complex-harmonic interpretation | guarded and deferred |

No TD-DFT, `chi0`, Dyson response, Ward correction, terminator, or radial solver rewrite was added.

## Call graph: eigenvector to radial SCF density

```text
self%run
  -> run_recursion
       -> symbolic_atom%build_pot
       -> hamiltonian%build_bulkham
            -> build_obarm / build_enim / eeo = ee * obarm
  -> run_dos
       -> reciprocal%build_kspace_hamiltonian [use_kspace path]
            -> reciprocal assembler -> h(k), eeo(k), hoh(k)=eeo(k)h(k)
            -> eigensolver
       -> reciprocal%accumulate_spin_density_kspace
            -> spin blocks and M0,M1,M2
            -> fill_band_moments_from_spin_density
            -> bands/scf QL and gravity-center parameters
  -> run_scf
       -> atomsc -> NEWRHO
            -> RSEQSR
            -> PHDFSR -> RSEQSR at shifted boundary slopes
            -> GINTSR metric -> radial density polynomial
```

The ordinary real-space DOS branch has the analogous `bands%accumulate_spin_density_rs -> bands%calculate_moments` route. The radial endpoint is shared: the next atomic update calls `NEWRHO` with `PL`, `QL`, `EV`, potential, and mesh.

Relevant live locations are `source/self.f90` (`run_recursion`, `run_dos`, `run_scf`, `NEWRHO`, `RSEQSR`, `PHDFSR`, `GINTSR`), `source/reciprocal_fourier.f90`, `source/reciprocal_spin_density.f90`, `source/bands.f90`, and `source/spin_density.f90`.

## Radial solver contract

### Mesh, units, components — derived from code

`atomsc` constructs the mesh by a logarithmic recurrence of the form

```text
r_i + B = B exp(A(i-1)),   r_1 = 0.
```

The test fixture uses the same form with `A=0.03`, `B=0.10`, `NR=101`, `Z=1`, and zero external radial potential. Energies and Hamiltonian entries use the code's atomic/Ry convention; no Fermi-energy-zero assumption is made.

`RSEQSR` solves the scalar-relativistic equation. Its live comment states that the large-component boundary data are for `U(R)` with `psi=(U/R)Y_lm` and that the output is normalized to one. Its two columns are large and small radial components; angular harmonics are not in `G`.

The normalization metric is, at each non-origin point,

```text
(r+B) [ GFAC*G_large**2 + G_small**2 ]
GFAC = 1 + l(l+1)/(TMC*r)**2
TMC  = C - (V - 2Z/r - E)/C
C    = 274.074
```

`GINTSR`, `RSEQSR`, and `NEWRHO` use the same logarithmic-coordinate Simpson weighting. The origin is skipped because of `1/r`; the new API records diagnostic `GFAC=1` there but does not integrate that point. `ISP=1,2` are local scalar radial channels selected per site/species through the symbolic-atom potential.

The radial equation uses the internal atomic/Ry convention (the Coulomb term
is `-2Z/r` and `C=274.074`); the lattice `alat` field, documented elsewhere
in Å, is not the radial mesh unit. `RSEQSR` performs volume normalization with
the metric above. The later `VAL/SQRT(SUM)` operation rescales boundary data
passed to `PHDFSR`; it is not a replacement for that volume normalization.

### Provenance of `phi`, `phidot`, `phiddot` — derived from code

`NEWRHO` calls `RSEQSR` at `EVAL=EV(IVAL)`, rescales its boundary value and slope using the returned norm, then calls `PHDFSR`. `PHDFSR` is numerical differentiation: it sets `DELE=0.003`, perturbs the boundary slope using `DDDE=-RMAX/G(NR)**2`, calls `RSEQSR` twice, and applies three-point Lagrange first- and second-derivative weights to all large/small entries. `GP` and `GPP` are therefore `dG/dE` and `d2G/dE2` at the `E` passed to `PHDFSR`, with finite-difference noise.

`source/self.f90::legacy_radial_fixture` is only a seam around unchanged `RSEQSR`/`PHDFSR`; `legacy_newrho_fixture` calls unchanged `NEWRHO` for the live comparison.

## Normalization identities

For the real scalar-relativistic pair and the `GINTSR` metric,

```text
<phi|phi> = 1
Re <phi|phidot> = 0
<phidot|phidot> + Re <phi|phiddot> = 0.
```

`tests/unit/test_lr_basis_radial.f90` uses a real solution from the live solver and the same Simpson metric. It reports:

```text
<phi|phi>-1                                      -1.11022302E-16
Re <phi|phidot>                                  -1.42096174E-07
<phidot|phidot> + Re <phi|phiddot>                2.07827473E-07
```

Tolerances are `2e-10` for normalization and `5e-5` for derivative identities. The latter reflects `PHDFSR` finite-difference noise.

## SCF density identity

### Live `NEWRHO` expression

For each `(l,ISP)`, `NEWRHO` reads `Q0=QL(1)`, `Q1=QL(2)`, `Q2=QL(3)` and evaluates

```text
rho += Q0 [ GFAC*G**2 + Gs**2 ]
     + 2*Q1 [ GFAC*G*GP + Gs*GPs ]
     + Q2 [ GFAC*(GP**2 + G*GPP) + GPs**2 + Gs*GPPs ].
```

This is the exact live expression in `source/self.f90::NEWRHO`.

### Independent derivation

With `delta=epsilon-E_nu_work`,

```text
phi(E_nu_work+delta) = phi + delta*phidot + (delta**2/2)*phiddot + ...
```

the density through second order is

```text
cross0 = |phi|**2
cross1 = Re(phi^dagger*phidot)
cross2 = |phidot|**2 + Re(phi^dagger*phiddot)
rho = Q0*cross0 + 2*Q1*cross1 + Q2*cross2,
```

with `GFAC` applied to the large-large products and unit weight to small-small products. Thus the density has `phiddot` even though the certified LMTO state amplitude is the first-order `phi+delta*phidot`; an explicit second-order `phiddot/2` term is not inferred for every state.

### `M` to `Q`

`reciprocal_spin_density.f90` accumulates, per site/channel/spin,

```text
M0 = sum(w_k f_nk)                 c c^dagger
M1 = sum(w_k f_nk epsilon_nk)      c c^dagger
M2 = sum(w_k f_nk epsilon_nk**2)   c c^dagger.
```

It forms the spin block from `u*conjg(u)`, `d*conjg(d)`, `u*conjg(d)`, and `d*conjg(u)`, then multiplies by `wk*occ*epsilon**p`. Relative to a working energy,

```text
Q0 = M0
Q1 = M1 - E_nu_work*M0
Q2 = M2 - 2*E_nu_work*M1 + E_nu_work**2*M0.
```

`spin_density%radial_band_moments` reports occupation, first-energy centre, and central second moment. `bands%calculate_moments` writes occupation to `QL(1)`, writes the measured centre minus `VMAD` to `gravity_center`, writes the central variance to `QL(3)`, and sets `QL(2)=0` in that SCF parameter interface because its selected centre is represented by the scalar `EV`/gravity-centre parameter. The augmentation API keeps the general `Q1` mapping explicit.

## Working linearization energy

`symbolic_atom%predls` reads raw `potential%enu` and constructs

```text
center_band  = (C-ENU)*X + ENU + VMAD
shifted_band = (C-ENU)*X
obar         = Y.
```

These become `potential%cx`, `potential%cex`, and `potential%obx`. `build_enim` then forms `eu=cx-cex` and `ed=cx-cex` for the two spin channels. Therefore the onsite quantity carried by `enim` is exactly

```text
E_nu_work = center_band - shifted_band = ENU + VMAD.
```

This is distinct from raw `potential%enu`, the measured band gravity centre, and `shifted_band`. `bands` writes the measured gravity centre as the band centre minus `VMAD`; it is used to update potential parameters and `PL`, not silently substituted for `enim`.

## HOH and the baseline wavefunction

### Hamiltonian identity — code-derived and test-confirmed

`hamiltonian_build` links each directed block as `eeo_ij=ee_ij*obarm_j`. The reciprocal assembler forms

```text
h(k)   = sum_R ee(R)  exp(i k.R)
eeo(k) = sum_R eeo(R) exp(i k.R)
hoh(k) = eeo(k)*h(k)
H(k)   = h(k)-hoh(k)+enim+lsham(+ccor).
```

With `lsham=0` and no CCOR this gives `hgamma=h-h*Obar*h` and `H=E_nu_work+hgamma`. `eeoee` is only the same-bond historical diagnostic, not the global HOH contraction.

`tests/unit/test_lr_basis_augmentation.f90` obtains `h(k)` and `eeo(k)` through the production Fourier routine, obtains `H(k)` through `reciprocal%build_hamiltonian_at_kpoint`, and compares them. The `spd` result is `1.73472348E-18`. All 18 production eigenvectors give maximum residual `1.58401437E-15` for

```text
(h-h*Obar*h)c_n = (epsilon_n I-E_nu_work)c_n.
```

### Formal LMTO interpretation

For the certified baseline,

```text
hgamma*c_n = (epsilon_n I-E_nu_work)*c_n
Psi_n(r) ~= Phi(r)*c_n + Phidot(r)*hgamma*c_n
          = [Phi(r)+Phidot(r)*(epsilon_n-E_nu_work)]*c_n.
```

At site `a`, this is

```text
Psi_n(r_a) ~= sum_{L,sigma} c_{aLsigma,n} Y_L(rhat_a)
              [phi_{al sigma}(r_a)
               +(epsilon_n-E_{nu,al sigma}^work)*phidot_{al sigma}(r_a)] |sigma>.
```

`reconstruct_state` produces the large/small radial factors and deliberately does not multiply by `Y_L`; angular, site, and Bloch assembly remain caller responsibilities.

## State-wise density versus live `NEWRHO`

The `spd` test uses production eigenpairs, deterministic positive weights, live radial solutions/derivatives, and both spin blocks. It independently:

1. accumulates the per-eigenstate second-order density polynomial over orbitals, energies, weights, and spins;
2. constructs `M0,M1,M2` from the production eigenvectors and evaluates `density_from_moments`;
3. passes the same `QL` and `EV` data to unchanged `NEWRHO` through `legacy_newrho_fixture` and compares summed `l` channels.

Observed maximum differences:

```text
state-wise API density - independent state density     8.88178420E-16
state-wise density - moment polynomial                 8.88178420E-16
state-wise density - live NEWRHO                       8.88178420E-16
```

This is a state/eigenvector-to-radial-density check, not only a moments-versus-moments check.

## Angular basis and transforms

The live spherical ordering is

```text
sp:  (00), (1,-1), (1,0), (1,1)
spd: (00), (1,-1), (1,0), (1,1), (2,-2), (2,-1), (2,0), (2,1), (2,2).
```

`basis.f90` gives `(norb,nb)=(4,8)`, `(9,18)`, and `(16,32)` for `sp`, `spd`, and `spdf`. The real/cubic ordering documented by the code is `(s,px,py,pz,dxy,dyz,dzx,x2-y2,3z2-r2)`. `hcpx` performs `H_sph=V^dagger H_cart V` and the inverse with explicit `1/sqrt(2)` matrices; `hamiltonian_build` applies it to all spin blocks of `obarm` and `enim`.

`UnitLrBasisAngular` round-trips deterministic Hermitian matrices with errors `1.66533454E-16` (`sp`) and `4.44089210E-16` (`spd`).

The `hcpx` `n=16` branch is an explicit identity fallback; its source comment says that no Cartesian-to-complex f transform is implemented and records a real f ordering. Coefficients 10–16 are therefore not certified as complex `Y_3m`. The test verifies the identity fallback and the augmentation API rejects `lmax=3`.

## Spinor and local-frame semantics

Normal site-local storage is site-major and spin-blocked:

```text
site_start+1       : site_start+norb       -> spin component 1
site_start+norb+1  : site_start+2*norb     -> spin component 2.
```

`reciprocal_spin_density` uses `u` and `d` from those blocks and stores `rho_uu=|u|^2`, `rho_dd=|d|^2`, `rho_ud=u*conjg(d)`, `rho_du=d*conjg(u)`. `spin_density` defines `rho=(n I+m.sigma)/2` and projects after accumulation with `up=(n+m.axis)/2`, `down=(n-m.axis)/2`. The stored density is rotating-frame; the axis is explicit. Collinear z-axis calculations reduce to the two diagonal radial channels certified here.

The live noncollinear `obarm`/`enim` builders create off-diagonal spin blocks from the local moment and transform them angularly. That contract is documented but excluded from the baseline API.

## SOC and additive operators

The assembler adds `lsham` and may add CCOR/other blocks, so the general relation is

```text
H = E_nu_work+hgamma+H_SO+H_extra
hgamma*c = (epsilon I-E_nu_work-H_SO-H_extra)c.
```

The simple `(epsilon-E_nu_work)c` relation is invalid there. The new API requires `has_soc=.false.` and `has_extra_operator=.false.`; `UnitLrBasisAugmentation` checks SOC, noncollinear, generalized-overlap, additive, and spdf rejection. Full SOC radial-response/augmentation is deferred.

## Minimal API and physical-GF seam

`source/lmto_radial_augmentation.f90` adds the narrow `lmto_radial_basis` type. It captures live large/small arrays and derivatives, stores radial and working energies, reconstructs the first-order radial state, evaluates the second-order density polynomial state-wise or from `M0/M1/M2`, exposes the legacy normalization metric, and fails closed outside the certified capability tuple.

The capture is an audit snapshot because legacy radial arrays are local to the atomic update; it does not replace or migrate the atomic solver. A future production seam can pass existing arrays directly.

The later physical Lehmann seam is explicit but not implemented:

```text
G_phys(k;r,r';z) = sum_n Psi_nk(r) Psi_nk(r')^dagger/(z-epsilon_nk)
                  = A G_c A^dagger,
A = Phi+Phidot*hgamma.
```

No response kernel, transverse operator, or Dyson implementation belongs here.

## Executable evidence and checklist

```text
cmake -S . -B build -DRUN_UNIT_TESTS=ON -DENABLE_SPGLIB=OFF
cmake --build build --target UnitLrBasisRadial UnitLrBasisAngular UnitLrBasisAugmentation -j2
ctest --test-dir build -R '^UnitLrBasis(Radial|Angular|Augmentation)$' --output-on-failure
```

All three tests passed: `UnitLrBasisRadial`, `UnitLrBasisAngular`, and `UnitLrBasisAugmentation`.

- [x] Current post-purge code traced.
- [x] `phi`, `phidot`, `phiddot` provenance and radial metric documented.
- [x] `M0,M1,M2`, `Q0,Q1,Q2`, and `NEWRHO` polynomial established.
- [x] `E_nu_work=ENU+VMAD` and the HOH identity established.
- [x] Baseline wavefunction and eigenstate identity tested against production matrices.
- [x] State-wise density agrees with both the moment route and live `NEWRHO`.
- [x] `sp`/`spd` transforms tested; `spdf` explicitly guarded.
- [x] Spinor/local-frame semantics documented; noncollinear excluded.
- [x] SOC/additive limitations derived, guarded, and deferred.
- [x] No legacy atomic numerical machinery or TD-DFT response code was rewritten.

Certified stopping statement:

```text
Psi_nk(r) = [Phi(r)+Phidot(r)*(h-h*Obar*h)] c_nk
          = [Phi(r)+Phidot(r)*(epsilon_nk-E_nu_work)] c_nk
```

for reciprocal orthogonal second-order `ham_only`, collinear, no-SOC, no-additive-operator `sp`/`spd` calculations, with direct numerical agreement against the production SCF radial-density machinery.
