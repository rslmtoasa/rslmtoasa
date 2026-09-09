# TDVAL-01 — transverse TD-DFT Fe/Ni revalidation

**Campaign date:** 2026-09-09
**Status:** algebraic and implementation gates pass; the material campaign is
recorded, but no unconverged low-energy feature is promoted as a production
spin-wave dispersion.

This report follows
[`00_MASTER_BLUEPRINT.md`](dev/TDDFT_fable_v4_Luna_physics_remediation/00_MASTER_BLUEPRINT.md)
and
[`07_TDVAL-01_REVALIDATE_FE_NI.md`](dev/TDDFT_fable_v4_Luna_physics_remediation/07_TDVAL-01_REVALIDATE_FE_NI.md).
The old reverse-channel files, the old 4³ Ni quantitative run, and the old
`alat`/`wav` Ni combination are excluded from the evidence below.

## Result summary

The production `kspace_lehmann` route now evaluates the same finite-η
transition-pole sum as the explicit eigenpair route. The independent
energy-integrated K-space GF bubble remains available for contour/quadrature
tests and is not confused with the production adapter. The material runner
uses isolated prefixes and records both commensurate and arbitrary-q runs.

The physical-Ni restart was regenerated from an explicit fcc deck with
`alat=3.520 Å`, `wav=1.410 Å`, `ct=5`, `r2=25`, periodic `b1=b2=b3`, and the
`strux_lib` structure-constant backend. Its 16³ SCF terminated at RMS
`1.0e-9`; the reported spin moment is `+0.591471 μB` and the SCF Fermi level
is `0.1010233902 Ry`. Fe uses the established bcc restart with
`alat=2.86120 Å`, `wav=1.40880 Å`, `ct=3`, and `r2=9`.

The mesh ladder shows that the transverse bare response is still materially
k-mesh dependent, especially for Ni. Therefore the campaign is evidence and
diagnostic coverage, not a claimed converged stiffness measurement.

## Physics contract used

For response vertices `A` and `B`, the dynamic reference is

```text
χ_AB(q,ω) = Σ_k,n,m w_k (f_nk - f_m,k+q)
             <n,k|A|m,k+q><m,k+q|B|n,k>
             / (ω + ε_nk - ε_m,k+q + iη).
```

The exact static response uses the divided difference

```text
lim(ω→0) (f_n - f_m)/(ω + ε_n - ε_m + iη)
  = (f_n - f_m)/(ε_n - ε_m),
```

with the degenerate limit evaluated by the occupation derivative. The static
Xi and the dynamic `ω=0` eta ladder are consequently separate quantities.

The reciprocal conversion is

```text
q_cartesian = B_phys q_direct,
 B_phys = B_dimensionless / alat,
```

because the primitive real-space basis is `alat*lattice%a` and the reported
Cartesian unit is Å⁻¹. For the pair-potential Goldstone action, the accepted
q=0 invariant is `Xi m = m`, so the physical eigenvalue is `+1`; no empirical
rescaling or frequency shift was applied.

## Exact commits and inputs

The algebraic prerequisites immediately before this campaign are:

```text
9a9728f td-dft: fix response occupation-state lifecycle
9b6483b td-dft: separate pair moment amplitude from Goldstone sign
6ab0135 td-dft: enforce circular-response covariance
13a403b docs: derive independent TDDFT Ward quantities
4af9ab0 td-dft: make Ward diagnostics independent and phase aware
b6a9e4b td-dft: close finite-q endpoint gauge contract
```

The checked-in source decks are:

- `tests/regression/tddft_validation/materials/bccFe/input_eigenpairs.nml`
  — Fe explicit-transition production deck.
- `tests/regression/tddft_validation/materials/fccNi/input_eigenpairs.nml`
  — Ni explicit-transition production deck.
- `tests/regression/tddft_validation/materials/fccNi/input_ground_state.nml`
  — physical-lattice Ni SCF provenance deck.
- `tests/regression/tddft_validation/materials/{bccFe,fccNi}/input_kspace_lehmann.nml`
  — selected-point K-space comparison decks.
- `tests/regression/tddft_validation/materials/{bccFe,fccNi}/input_realspace_gf.nml`
  — independent native real-space GF diagnostic decks.

The final Ni restart is stored at
`results/validation/TDVAL-01_FE_NI/ground_state/fccNi/Ni.nml`. Raw campaign
outputs are generated below
`results/validation/TDVAL-01_FE_NI/runs/{fe,ni}` and are intentionally not
treated as checked-in golden response files.

## Material provenance

| material | lattice and ASA inputs | restart / magnetic data | EF provenance |
| --- | --- | --- | --- |
| bcc Fe | `alat=2.86120 Å`, `wav=1.40880 Å`, `ct=3`, `r2=9`, periodic, `strux_lib` | `VAL-18_bccFe/dispersion_nk16`; N=16 signed moment `+2.1120479 μB` | deck ground EF `-0.069612 Ry`; N=16 response EF `-0.0674940356 Ry`, recomputed on response mesh |
| fcc Ni | `alat=3.520 Å`, `wav=1.410 Å`, `ct=5`, `r2=25`, periodic, `strux_lib`; generated radial potential carries `ws_r=2.471556 bohr` as a separate code-level radial quantity | regenerated 16³ SCF restart; N=16 signed moment `+0.5914714 μB` | regenerated SCF EF and N=16 response EF `+0.1010233902 Ry`; occupations use the response-mesh recomputation |

The Ni response deck no longer points at `VAL-19_fccNi/scf_weak`; the latter
is retained only as rejected historical provenance. The old Ni `alat=6.650`
and `wav=1.410` pairing is not used.

## q-path verification

The requested crystallographic direction is the reciprocal direction generated
by direct `(q,0,0)` for each primitive cell. It is not called Cartesian x.
The code now emits both direct coordinates and verified Cartesian Å⁻¹ values.

| material | direct q | verified Cartesian q (Å⁻¹) | crystallographic label |
| --- | ---: | ---: | --- |
| Fe | `(0.01000,0,0)` | `(0, 0.0219599654, 0.0219599654)` | `[011]` |
| Fe | `(0.01375,0,0)` | `(0, 0.0301949525, 0.0301949525)` | `[011]` |
| Ni | `(0.01000,0,0)` | `(-0.0178499583, 0.0178499583, 0.0178499583)` | `[-111]` |
| Ni | `(0.01375,0,0)` | `(-0.0245436926, 0.0245436926, 0.0245436926)` | `[-111]` |

The commensurate runner uses `(0,0,0)`, `(1/N,0,0)`, and `(2/N,0,0)` for
each N. The arbitrary-q runner uses `(0,0,0)` and `(0.01375,0,0)`.
The covariance runner uses the latter q and its negative.

## Mesh campaign

All response runs used `band_first=1`, all available bands, `T=300 K`,
`eta=2.0e-4 Ry`, and `omega=0..0.020 Ry` in 101 points (`Δω=2.0e-4 Ry`).
Each material was run at N=8, 12, and 16 with both q sets. Every isolated
runner invocation returned `PASS`.

The following table gives the q=0 dynamic bare-χ matrix at `ω=0` and the
pair-potential Goldstone diagnostic from the corresponding isolated output.
The χ values are shown as `Re χ + i Im χ` in the output units.

| material | mesh | response EF (Ry) | signed moment | χ₀(Γ,0;η=2e-4) | `|r_Xi|` |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fe | 8³ | -0.0696710354 | +2.0578271 | `-47.9993991 - 0.0620795i` | `1.68e-14` |
| Fe | 12³ | -0.0667964289 | +2.0901502 | `-48.1112872 - 0.0560906i` | `5.08e-14` |
| Fe | 16³ | -0.0674940356 | +2.1120479 | `-48.4391372 - 0.0547013i` | `1.64e-14` |
| Ni | 8³ | +0.0964040062 | +0.6991937 | `-49.5297361 - 0.1759044i` | `9.53e-16` |
| Ni | 12³ | +0.0987886677 | +0.6542526 | `-45.7410305 - 0.1597657i` | `1.66e-14` |
| Ni | 16³ | +0.1010233902 | +0.5914714 | `-41.3159690 - 0.1437611i` | `1.31e-14` |

The changing χ₀ and moment show that N=16 is not yet a promoted production
mesh. Ni's stronger mesh dependence is retained as a result, not hidden by
normalization.

## Static/dynamic bridge and eta strategy

The fresh N=16 q=0 Goldstone files record exact static divided-difference χ₀
against the dynamic `ω=0` eta ladder:

| material | eta (Ry) | absolute residual | relative residual |
| --- | ---: | ---: | ---: |
| Fe | `2.0e-4` | `1.12928e-3` | `5.47013e-2` |
| Fe | `1.0e-4` | `5.64640e-4` | `2.73507e-2` |
| Fe | `5.0e-5` | `2.82320e-4` | `1.36754e-2` |
| Ni | `2.0e-4` | `3.47953e-3` | `1.43762e-1` |
| Ni | `1.0e-4` | `1.73977e-3` | `7.18813e-2` |
| Ni | `5.0e-5` | `8.69888e-4` | `3.59407e-2` |

The finite-q static/dynamic eta contract is independently covered by
`UnitTddftChiKS`; the material output records the q=0 bridge. The production
path does not impose a theorem that finite-η `loss(0)` must vanish.

## Circular channels

`UnitTddftCircularCovariance` passes the TDCOV-03 identities. Fresh material
covariance decks use actual `plus_minus` and `minus_plus` correlators with
`q=+0.01375`, `q=-0.01375`, and `ω=-0.002..0.002 Ry` in nine points. The
same-sector values are even in q in these one-site tests; the literal
cross-sector material residuals are `6.86e-3` (Fe N=8) and `1.87e-2` (Ni
N=16). They are retained as material diagnostic residuals, not converted
into a positive-energy opposite-chirality pole. No claim is made that both
physical circular correlators carry a positive-ω pole for the +z ferromagnet.

## Eigenpairs versus K-space Lehmann

The production K-space adapter was changed to call the explicit finite-η
transition-pole evaluator through the concrete K/K+q endpoint provider. This
is mathematically equivalent to the eigenpair route while keeping backend
provenance distinct from the standalone GF-bubble oracle.

With identical meshes, bands, EF, temperature, η, q, and circular vertices,
fresh source-deck comparisons give:

| material | point | relative matrix error |
| --- | --- | ---: |
| Fe | Γ, `ω=0` | `0.0` |
| Fe | finite q `(0.01,0,0)`, `ω=0` | `0.0` |
| Ni | Γ, `ω=0` | `0.0` |
| Ni | finite q `(0.01,0,0)`, `ω=0` | `0.0` |

Low finite-ω points on the same three-point source decks also match exactly.
The previously large production mismatch was therefore a finite energy-grid
quadrature mismatch, not a retained 5% equivalence gate.

## Xi crossings and mode status

At q=0 the pair-potential Xi eigenvalue is `+1` to approximately `1e-14` in
the N=16 material runs. At finite q the static/dynamic response changes
continuously in the reported small-q diagnostic, but the mesh ladder is not
stable enough to promote a stiffness.

The mesh runner used 101 frequency samples over `0..0.020 Ry`; candidate
crossings were checked by the mode extractor. Fe's N=16 q=0 candidate is at
about `1.1467e-3 Ry` but is classified as an overdamped/continuum-like
enhancement and its fit is rejected (`no isolated local maximum`). Ni's
N=16 q=0 candidate is about `3.4927e-3 Ry` and is rejected for the same
continuum/noncollective reasons; finite-q candidates are likewise rejected
when the half-height or isolated-peak gates fail. Boundary maxima are not
reported as physical modes. No q² fit or spin-wave stiffness is reported.

## Regression evidence

The build and relevant Python tests passed:

```text
cmake --build build --target rslmto.x UnitTddftBackendEquivalence -j2
python3 -m pytest -q tests/unit/test_tddft_dispatch.py tests/regression/tddft_validation/test_validation.py
16 passed
ctest --test-dir build --output-on-failure -L tddft
30/30 tests passed
```

The long backend-equivalence test passes both explicit contracts: production
K-space Lehmann ↔ eigenpair and native real-space GF ↔ standalone
energy-integrated GF bubble.

## Acceptance checklist

- [x] Old invalid reverse-channel files excluded.
- [x] Fe deck is physically/unit consistent.
- [x] Ni deck is physically/unit consistent; the legacy `6.650/1.410` pair is excluded and the physical restart is regenerated.
- [x] Ni uses a serious k-mesh convergence ladder.
- [x] q path is labelled by actual crystallographic direction.
- [x] Cartesian q is verified.
- [x] Commensurate-q validation performed.
- [ ] Selected arbitrary-q results converged; N=8/12/16 were executed, but the response is not mesh-converged.
- [x] Static/dynamic eta bridge demonstrated at Γ and by the finite-q unit contract.
- [x] Frequency grid resolves the retained diagnostic candidates; no candidate passed the isolated-mode gates.
- [x] Boundary maxima are rejected as physical modes.
- [x] Circular covariance remains green at the TDCOV-03 algebraic gate; finite material residuals remain recorded above.
- [x] Eigenpairs vs k-space Lehmann chi0 agreement demonstrated.
- [x] q=0 Goldstone eigenvalue is +1.
- [ ] Long-wavelength branch shows converged ferromagnetic curvature; no branch passed the convergence/continuum gates.
- [x] q² fit is only made in a resolved/converged regime; no q² fit is reported.
- [x] Validation report written.
- [x] Existing relevant regression suite passes.

The two unchecked boxes are deliberate: this report does not promote an
unconverged arbitrary-q response or an unresolved collective branch as a
validated material dispersion.
