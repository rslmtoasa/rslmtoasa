# VAL-02 — Reciprocal SCF convergence and magnetic scope

Date: 2026-08-16
Binary: `build-rf-debug/bin/rslmto.x`
Configuration: GNU Debug, serial CPU, `strux_backend='strux_lib'`

## Scope and result

This campaign establishes the present reciprocal-SCF envelope for the two
canonical fixtures:

- diamond Si, `sp`, `nsp=1`, nonmagnetic;
- bcc Fe, `spd`, `nsp=2`, collinear magnetic.

Within that envelope, reciprocal SCF is **Validated** for finite-temperature
(`T=300 K`) Gaussian/tetrahedron DOS evaluation, automatic canonical EF, and
the stated finite k-mesh ranges. Electron counting and site charges are stable.
The Fe moment is converged with respect to SCF iterations, but remains
metallic-mesh sensitive; the 12³→16³ result is the useful compact magnetic
envelope, not an infinite-mesh materials benchmark.

`nsp=3` and `nsp=4` both execute the reciprocal spinor path. A compact
second-order (`hoh=.true.`) check with an x-axis Fe seed passed Hermiticity,
electron count, spin-density orientation, and moment magnitude checks. They
remain **Experimental/conditional**: the ordinary first-order path explicitly
omits `enim+lsham`, so `nsp=4, hoh=.false.` is not established as a fully
relativistic reciprocal-SCF combination; the extended checks also retain a DOS
state-count residual and do not exercise a genuinely transverse multi-site
texture.

## Campaign protocol

The reproducible, non-CTest driver is
[`val02_reciprocal_scf.py`](../../tests/validation/val02_reciprocal_scf.py).
It uses isolated scratch directories and is not registered in `quick`.

The base campaign used unreduced Gamma-centered meshes, `ham_only`, automatic
EF, Gaussian `sigma=0.01/0.02 Ry`, and a 200-point DOS grid. Tetrahedron runs
used 2001 points. The SCF axis used meshes 4³ and `nstep=1,3,6,12,24`; mesh
sweeps used `nstep=12`. Additional Fe mesh checks used 8³, 10³, 12³, and 16³.
Site charge is the reported occupation/charge-transfer pair; the Fe moment is
the reported Cartesian spin-moment magnitude/vector.

## k-mesh convergence

The DOS state count is the trapezoidal integral of the production
`totaldos.out` file; the expected counts are 16 for the two-site Si/sp basis
and 18 for the one-site Fe/spd basis.

| fixture | mesh | EF (Ry) | N | EBAND (Ry) | site charge | Fe moment (μB) | report Etot (Ry) | DOS states |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Si/sp | 2³ | -0.012519 | 8.000000 | -3.769470 | 4.000000 | — | -1157.366484 | 16.000239 |
| Si/sp | 4³ | 0.002134 | 8.000000 | -3.887129 | 3.999998 | — | -1157.495434 | 16.000168 |
| Si/sp | 6³ | 0.003934 | 8.000000 | -3.888459 | 4.000001 | — | -1157.499928 | 15.999748 |
| Si/sp | 8³ | 0.005573 | 8.000000 | -3.888241 | 4.000000 | — | -1157.500342 | 16.000054 |
| Si/sp | 10³ | 0.001849 | 8.000000 | -3.888164 | 4.000012 | — | -1157.500408 | 16.000032 |
| Fe/spd | 2³ | -0.005630 | 8.000000 | -1.679822 | 8.000000 | 3.499997 | -2541.847414 | 18.012832 |
| Fe/spd | 4³ | -0.061169 | 8.000000 | -1.994591 | 8.000000 | 1.950775 | -2541.981342 | 17.992897 |
| Fe/spd | 6³ | -0.076879 | 8.000000 | -1.981084 | 8.000000 | 2.222220 | -2541.983934 | 17.989567 |
| Fe/spd | 8³ | -0.066983 | 8.000000 | -1.986686 | 8.000000 | 2.003490 | -2541.979291 | 18.010067 |
| Fe/spd | 10³ | -0.061495 | 8.000000 | -1.994350 | 8.000000 | 2.129566 | -2541.981902 | 18.003035 |
| Fe/spd | 12³ | -0.067574 | 8.000000 | -2.000063 | 8.000000 | 2.087162 | -2541.981909 | 17.997579 |
| Fe/spd | 16³ | -0.068065 | 8.000000 | -1.999080 | 8.000000 | 2.124505 | -2541.981251 | 18.005218 |

The Si charge, EBAND, total energy, and DOS state count settle strongly by
4³–10³. Its EF is not a unique mesh-invariant observable because Si is
gapped and the finite-temperature chemical potential can move inside the
gap; at fixed 4³ it is stable to 2.7×10⁻⁸ Ry between 12 and 24 SCF steps.

Fe is not converged by the canonical 2³–4³ fixture. The denser 12³→16³
change is 4.91×10⁻⁴ Ry in EF, 9.84×10⁻⁴ Ry in EBAND, 6.59×10⁻⁴ Ry in total
energy, and 0.0373 μB in moment. That finite-mesh spread is retained as a
scope limitation.

## SCF iteration convergence

At fixed 4³ Gaussian `sigma=0.01 Ry`:

| fixture | nstep | final reported RMS diff | EF (Ry) | EBAND (Ry) | site charge | Fe moment (μB) |
|---|---:|---:|---:|---:|---:|---:|
| Si/sp | 1 | 5.19×10⁻² | 0.001456 | -4.104472 | 3.550000 | — |
| Si/sp | 3 | 2.44×10⁻³ | 0.007609 | -3.827620 | 3.971877 | — |
| Si/sp | 6 | 1.41×10⁻³ | 0.001820 | -3.878272 | 4.000221 | — |
| Si/sp | 12 | 2.91×10⁻⁷ | 0.002134 | -3.887129 | 3.999998 | — |
| Si/sp | 24 | 2.46×10⁻⁸ | 0.002134 | -3.887133 | 4.000000 | — |
| Fe/spd | 1 | 1.74×10⁻² | -0.063611 | -1.996867 | 8.000000 | 1.986713 |
| Fe/spd | 3 | 2.64×10⁻³ | -0.059857 | -1.973923 | 8.000000 | 1.945045 |
| Fe/spd | 6 | 5.83×10⁻⁴ | -0.061198 | -1.995589 | 8.000000 | 1.952531 |
| Fe/spd | 12 | 6.13×10⁻⁶ | -0.061169 | -1.994591 | 8.000000 | 1.950775 |
| Fe/spd | 24 | 2.48×10⁻⁸ | -0.061166 | -1.994566 | 8.000000 | 1.950773 |

N is 8.000000 with canonical residual below 3×10⁻¹¹ in every reported
reciprocal solve; the raw k-weight sum is 1.000000. At fixed mesh, EF and
EBAND are unchanged at the shown precision once the SCF residual is below the
10⁻⁶–10⁻⁷ level. The Fe site charge is exactly 8.000000 throughout the
iteration sweep.

## DOS method and normalization

At 4³, the spectral observables from Gaussian `sigma=0.01`, Gaussian
`sigma=0.02`, and tetrahedron have the same canonical EF, N, and EBAND. The
DOS state-count evidence is:

| fixture | DOS route | DOS states from file | site charge | Fe moment (μB) |
|---|---|---:|---:|---:|
| Si/sp | Gaussian, σ=.01 | 16.000168 | 3.999998 | — |
| Si/sp | Gaussian, σ=.02 | 16.000010 | 3.999998 | — |
| Si/sp | tetrahedron | 15.999650 | 3.999998 | — |
| Fe/spd | Gaussian, σ=.01 | 17.992897 | 8.000000 | 1.950775 |
| Fe/spd | Gaussian, σ=.02 | 18.000129 | 8.000000 | 1.950775 |
| Fe/spd | tetrahedron, 4³ | 18.109108 | 8.000000 | 1.950775 |
| Fe/spd | tetrahedron, 8³ | 17.995531 | 8.000000 | 2.003490 |

The 4³ Fe tetrahedron residual is a coarse-metal mesh effect; the 8³ check
restores the state count to within 0.025%. Gaussian normalization is within
0.07% over the Fe mesh sweep. The nsp=3/4 second-order spin checks below are
not used to promote DOS normalization: their finite-window Gaussian integrals
are 19.038237 and 18.990386 states, respectively, against the 18-state
spanning basis.

Meaningful energy components were also retained from the production `par`
records. Values below are per site; the Si report total is the sum over its
two equivalent sites.

| route | mesh | sumec | sumev | etot | utot | ekin | rhoeps |
|---|---:|---:|---:|---:|---:|---:|---:|
| Si/sp | 4³ | -320.238310 | -1.943564 | -578.747717 | -1117.412207 | 578.910848 | -40.246358 |
| Si/sp | 10³ | -320.222277 | -1.944060 | -578.750204 | -1117.429971 | 578.927075 | -40.247309 |
| Fe/spd | 4³ | -1486.027409 | -1.994559 | -2541.981342 | -5008.472295 | 2574.696870 | -108.205917 |
| Fe/spd | 16³ | -1485.935424 | -1.999080 | -2541.981251 | -5008.538906 | 2574.789800 | -108.232145 |

## RS comparison

The existing production comparison was rerun with
[`test_si_dos_equivalence.py`](../../tests/scf/test_si_dos_equivalence.py):
Chebyshev order 200 versus KS 8³ with fixed Gaussian σ=0.020 Ry. It passed
with:

- integrated states: RS 16.00000568, KS 16.00000232;
- central-window states: RS 15.99988193, KS 15.99749758;
- central DOS relative RMS: 0.139495;
- maximum absolute DOS difference: 11.465596;
- principal peak positions: 0.549142 Ry (RS), 0.544106 Ry (KS).

This is a spectral comparison only. No Block/Lanczos DOS equality was
demanded. Reciprocal and real-space charges/moments were not asserted equal:
the available comparison driver does not establish identical converged
Hamiltonian/potential conventions for those two SCF routes.

## Extended spin scope: nsp=3 and nsp=4

The control documentation defines nsp=3 as noncollinear scalar-relativistic
and nsp=4 as noncollinear fully relativistic. The reciprocal code stores a
two-component spinor, computes Cartesian spin-density components, disables
time-reversal reduction for nsp≥3, and enforces Hermiticity before
eigensolution. The compact functional check used `hoh=.true.`, a 2³ mesh, 12
SCF steps, an x-axis Fe seed, and the wider `[-4,4] Ry` DOS window:

| class | EF (Ry) | N | EBAND (Ry) | site charge | moment vector (μB) | moment magnitude (μB) | transverse density | Hermiticity |
|---|---:|---:|---:|---:|---|---:|---|---|
| nsp=3 | -0.068425 | 8.000000 | -1.932166 | 8.000000 | (3.500158, -0.000000, 0.000000) | 3.500158 | 0.000000 | pass |
| nsp=4 | -0.068434 | 8.000000 | -1.932754 | 8.000000 | (3.499573, -0.000000, -0.000000) | 3.499573 | 0.000000 | pass |

The x-axis orientation is preserved and the moment magnitude is finite. The
zero transverse component is expected for this one-site uniformly rotated
fixture; it is not evidence for a textured noncollinear material calculation.
The DOS state-count residual and the first-order SOC omission limit both
classes to the conditional scope recorded above.

## Validation checklist

- [x] Si k-mesh convergence measured.
- [x] Fe k-mesh convergence measured, including 12³/16³ magnetic envelope.
- [x] Electron counting stable (`N=8`, `weight_sum=1`, canonical residual <3×10⁻¹¹).
- [x] EF stability characterized: fixed-mesh SCF stability established; Si gap and coarse-metal mesh dependence recorded.
- [x] Charges/moments stability characterized; the remaining Fe dense-mesh spread is recorded.
- [x] DOS normalization retained for the established Si/nsp=1 and Fe/nsp=2 envelope; coarse tetrahedron and nsp=3/4 exceptions are explicit.
- [x] RS comparison uses Chebyshev with controlled KS Gaussian broadening.
- [x] Intended nsp=3 scope established narrowly as conditional second-order spinor SCF.
- [x] Intended nsp=4 scope established narrowly as conditional second-order/SOC spinor SCF; default first-order fully relativistic support is not claimed.
- [x] Validation evidence recorded.
- [x] Maturity ledger updated narrowly.

## Commands and files

Campaigns run:

```text
ctest -L kspace --output-on-failure                         # 10/10 passed
python3 tests/validation/val02_reciprocal_scf.py ...       # base sweeps
python3 tests/validation/val02_reciprocal_scf.py --high-mesh-only ...
python3 tests/validation/val02_reciprocal_scf.py --dense-fe-only ...
python3 tests/validation/val02_reciprocal_scf.py --extended-only ...
python3 tests/validation/val02_reciprocal_scf.py --tetra-high-only ...
python3 tests/scf/test_si_dos_equivalence.py ...            # pass
```

Changed files are this report, the narrow Phase-II ledger entry, and the
non-CTest campaign driver. No production source or `quick` case was changed.
