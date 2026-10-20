# Scientific audit: time-odd longitudinal Dx in P4/mmm Mn2F2

## Result and scope

The earliest **measured** violation is nonreciprocity of the supplied screened
structure constants, before Hamiltonian construction, GF decomposition or CC/SC
contraction. For the first Mn1->Mn2 bond,
`max|Sbar_ij - Sbar_ji^t| = 1.1189163346969622e-6`. These matrices produce a
non-Hermitian bare hopping block. An independent reconstruction from the binary
structure constants and saved potential widths agrees with the code's captured
bare hopping to `2.57e-16` Ry or better elementwise. HoH propagates the discrepancy
into onsite blocks as well.

The GF leakage has the specific parity that explains the observed odd Dx:
forbidden **same-mirror-parity** G01 is predominantly bare-T **even**, and forbidden
same-parity Gx1 is predominantly bare-T **odd**. Their products with the large allowed
Gx0 (odd) and G00 (even), respectively, are odd. This is visible in individual
orbital matrices and traces, before integration. It is not a finding that Gx0
should vanish.

This identifies a concrete upstream numerical reciprocity failure and its path
into the reported contractions. It does **not** establish whether the original
screening approximation, its convergence parameters, or an implementation issue
created the nonreciprocal structure constants. No production Hamiltonian, SCF
parameter, exchange formula, or standard output format was changed. No inference
of a physical higher-order interaction is justified by the odd output alone.

Baseline: current `main` at `e604f4029e59d35b6e53551dd2d95919a9baf8e0`.
Diagnostic branch: `codex/mn2f2-gf-symmetry-audit`; no merge into main.
The phase is the supplied P4/mmm monolayer, not P4/nmm FeSe-type Mn2F2.

## Derivation and actual basis

Let `a=(l,m)`, with source index `k=l(l+1)+m+1`. The numerical `hcpx`
transformation V has `cart=V*sph`. Direct tests of its entries and the source L
operators give

$$C=V^\dagger V^*,\quad C_{(l,m),(l',m')}=(-1)^m\delta_{ll'}\delta_{m,-m'},
\qquad R(A)=CA^tC.$$

In [green.f90](source/green.f90#L386), `k0=l(l+1)+1`, so
`(-1)^(k+j) A(2*j0-j,2*k0-k)` is exactly R, including **both** reversed m indices.
C is real, symmetric, involutory; R reverses multiplication order and is complex
linear, not conjugate linear. The numerical Lz eigenvalues in this convention
are `(0,1,0,-1,2,1,0,-1,-2)`, opposite the signs of the printed m ordering in the
`hcpx` comment. The mirror and time-reversal matrices were therefore checked
against numerical operators rather than inferred from that comment.

For `m_z`, reflection is about the Mn plane, not absolute z=0:

$$M_{aa}=(-1)^{l+m},\quad p=(1,1,-1,1,1,-1,1,-1,1),\quad [C,M]=0.$$

In the spin-block basis `(up:9,down:9)`, define

$$U_T=i\sigma_y\otimes C,\quad U_m=\sigma_z\otimes M,\quad
U_h=U_m U_T\simeq\sigma_x\otimes MC.$$

An overall unitary phase is irrelevant. `U_T U_T^*=-1` and `U_h U_h^*=+1`.
The source L operators, transformed with V, obey
`M Lx,y M=-Lx,y`, `M Lz M=Lz`, `C Lmu^t C=-Lmu`.
Consequently the implemented scalar-prefactor SOC `L.s` respects T and `m_z*T`
when orbital polarization is off, as in these inputs. This checks SOC conventions,
not just a formal parity array.

### Why a retarded reciprocity relation uses transpose

Write z in the upper half-plane. Antiunitary invariance gives

$$H=U_h H^*U_h^\dagger,\qquad
G(z)=U_h[G(z^*)]^*U_h^\dagger.$$

The stronger same-z transpose relation additionally requires **Hermiticity**:
`H^*=H^t`, so

$$G_{ij}(z)=U_h[G_{ji}(z)]^tU_h^\dagger.$$

Transpose is over the complete spin/orbital matrix, and exchanges site indices.
It is not elementwise conjugation of the retarded GF. For a non-Hermitian
Hamiltonian, the first relation does not imply the second. A synthetic
non-Hermitian counterexample included in the tests demonstrates this explicitly.
This missing assumption matters quantitatively in the supplied calculation. For an exact resolvent, put
`H'=U_h H^t U_h^dagger` and `G'=(z-H')^-1`. Then

$$G-G'=G(H-H')G'.$$

If magnetic antiunitary covariance holds, `H-H'` is controlled by the
non-Hermitian part. The high-energy expansion already contains this discrepancy
at order `H/z^2`. The production accelerated Chebyshev moment identities also
assume Hermitian H, so a non-Hermitian input can affect subsequent moments in
addition to the exact-resolvent effect. A separate synthetic model with no SOC
shows that numerical nonreciprocity alone can produce odd CCx/SCx; that test is
not an SOC-off calculation for Mn2F2.

With Hermiticity and an individually fixed in-plane directed Mn pair, Pauli
components satisfy

$$g^\mu_{ij}=s_\mu M R(g^\mu_{ji})M,\qquad
(s_0,s_x,s_y,s_z)=(1,1,1,-1).$$

For the code's `G^{mu,0/1}=(g^mu_ij +/- R(g^mu_ji))/2`, R's involution and
commutation with mirror conjugation imply

$$M G^{\mu,0}M=s_\mu G^{\mu,0},\qquad
M G^{\mu,1}M=-s_\mu G^{\mu,1}.$$

| Sector | Mirror sign | Allowed orbital blocks | Bare-T sign at matched ordered bond |
| --- | --- | --- | --- |
| G00 | + | same parity | + |
| G01 | - | opposite parity | - |
| Gx0, Gy0 | + | same parity | - |
| Gx1, Gy1 | - | opposite parity | + |
| Gz0 | - | opposite parity | - |
| Gz1 | + | same parity | + |

Neither Gx0 nor Gx1 is required to vanish. Synthetic matrices make both nonzero
while their forbidden traces vanish exactly.

The source `d_matrix()` is diagonal and constant within each l shell:

$$\Delta_l(E)=\frac{(c_{l\downarrow}+v_M)w_{l\uparrow}^2
 -(c_{l\uparrow}+v_M)w_{l\downarrow}^2
 +(w_{l\downarrow}^2-w_{l\uparrow}^2)E}
 {w_{l\uparrow}w_{l\downarrow}}.$$

Here w denotes the `dele` value used by that routine, not the separate hopping
width parameter. Thus `[Delta,M]=0`. For any two factors of opposite mirror sign,
cyclic invariance of trace gives `Tr(Delta_i A Delta_j B)=-Tr(...)`, hence

$$SC_x(E)=\Re\mathrm{Tr}[\Delta_iG00_{ij}\Delta_jGx1_{ji}]=0,$$
$$CC_x(E)=\Re\mathrm{Tr}[\Delta_iG01_{ij}\Delta_jGx0_{ji}]=0.$$

The same holds for y, **per bond and per energy**. For z, the respective signs
are `+*+` and `-*-`, so these traces can be nonzero.

### Bare T is a separate comparison

For independently saved reversed domains, with preserved local channel labels,
`H_- = U_T H_+^* U_T^dagger`. If both are Hermitian,

$$G^-_{ij}(z)=U_T[G^+_{ji}(z)]^t U_T^\dagger.$$

The scalar Pauli component has raw sign +, all vector components sign -. Splitting
with R gives the last table column. Delta is nominally even when local up/down
channels are preserved. Consequently CCx and SCx are each nominally even. Calling
a current operator time odd does not change the parity of the two-factor product.
The actual leakage violates these nominal sector parities as documented below.

## Input inventory and original-output reproduction

All four supplied jobs contain the same 24 ordered cluster-index pairs, identical
relative displacements, 1600 cluster atoms, four atomic types, a 2000-channel grid
(2010 allocated evaluation entries), `lld=250`, Chebyshev recursion, `hoh=T`,
`post_processing='exchange'`, and default `local_axis=F`, `orb_pol=F`.
There are 12 Mn1->Mn2 bonds among the 24. Actual execution logs confirm reading
`Mn1_out.nml`, `Mn2_out.nml`, `F1_out.nml`, `F2_out.nml`; these take precedence over
the plain atom files. No SCF was executed.

The binary `sbar.` is byte-identical in all four cases, SHA256
`951f587460c1ade64d6676e1a6ab591f3bd7de831fe5ff54c716cfa25f818c4a`.
Detailed inventories, input hashes, raw matrix dumps and results remain outside Git.

Original printed sums over the 12 intersublattice bonds, in mRy:

| Domain | Dx | Dy | Dz | CCx | SCx |
| --- | ---: | ---: | ---: | ---: | ---: |
| [100] | 1.82060157e-6 | 2.07e-15 | 0 | 6.06640940e-7 | 1.213960558e-6 |
| [-100] | -1.82060135e-6 | -1.13e-15 | 0 | -6.06640800e-7 | -1.213960550e-6 |
| [110] | 1.28719687e-6 | 1.28719688e-6 | 4.3708678e-3 | 4.29725510e-7 | 8.57471204e-7 |
| [-1-10] | -1.28719678e-6 | -1.28719678e-6 | 4.3706562e-3 | -4.29725460e-7 | -8.57471318e-7 |

`dijparts.out` prints CCxyz followed by SCxyz, already with the factor of two
used in second-order D. At printed precision `CC+SC=D`, with maximum rounding
difference `5e-10` mRy over the saved records.

Across all 24 matched directed bonds, [100]/[-100] have
`max|Dx_-+Dx_+|=2.76e-11` mRy and `max|Dz_- - Dz_+|=2e-9` mRy.
The maximum even-Dx discrepancy is `7.773483e-7` mRy. Thus the small Dx is a
reproducible odd signal, not a bond-order mismatch.

Inputs are approximately, not exactly, symmetry-related. For [100], top/bottom
F c-parameter differences reach `3.1444e-10` Ry. Under domain reversal, Mn1 c
parameters differ by up to `1.37348e-8` Ry and F c parameters by `2.22909e-8` Ry.
Mn moments reverse exactly at saved precision; the [110] F moment-direction
oddness residual reaches `2.20819e-4` (these are unit-direction parameters).
The nominal geometric mirror differs by about `1e-10` in lattice-coordinate
units. None of these differences alone establishes the cause of Dx; the separate
binary reciprocity and hopping reconstruction do.

## Earliest intermediate: screened constants and Hamiltonian action

Matrices below are read from unformatted binary records, not the four-decimal
matrix preview. For `r/alat=(-1/2,+1/2,0)`, compare source type 1 with type 2 at
opposite displacement:

| Quantity | Maximum element residual | Frobenius residual |
| --- | ---: | ---: |
| Sbar_ij - Sbar_ji^t | 1.1189163347e-6 | 2.2423591469e-6 |
| Bare h_ij - h_ji^dagger (Ry) | 2.8921200022e-8 | 9.9330128271e-8 |
| HoH_ij - HoH_ji^dagger (Ry) | 5.5173235879e-8 | 2.6884553685e-7 |
| Effective H_ij - H_ji^dagger (Ry) | 5.5026866503e-8 | 3.0223039565e-7 |
| Effective onsite H_11 - H_11^dagger (Ry) | 4.2645863618e-8 | 1.4071679717e-7 |

The source's bare block is reproduced as

$$h_{ij}=W_i(1_2\otimes\bar S_{ij}^{\,t})W_j,
\quad W_i=w_i^0+\boldsymbol\sigma\cdot\mathbf n_i\,w_i^1,$$

with the `hcpx` transform applied to orbital blocks. Here w0,w1 are the shell
expanded half-sum/difference of `width_band` channels. Onsite radial terms are
added separately and are Hermitian. Reconstructed offsite h agrees with the
capture at maximum residuals `2.57e-16` and `1.60e-16` Ry for ij and ji.

Already in bare h, using the same R split as the GF, the forbidden blocks are:

| Quantity | Forbidden-block Frobenius norm | Largest leaking element |
| --- | ---: | --- |
| scalar 01 | 3.4885027421e-8 Ry | orbital (2,1), both mirror +, magnitude 1.4460600011e-8 Ry |
| x1 | 4.0427910711e-9 Ry | orbital (2,7), both mirror +, magnitude 1.3876980081e-9 Ry |

Bare onsite h is Hermitian to about `8.4e-17` Ry. Captured energy-centre and SOC
terms have zero Hermiticity residual at output precision; HoH introduces the
observed onsite discrepancy. The same actual effective H columns nevertheless
obey antiunitary magnetic-mirror covariance to maximum `1.20358e-10` Ry. There is
therefore no logical contradiction: antiunitary covariance alone does not supply
the Hermiticity needed to replace advanced conjugation by retarded transpose.

The largest screened-constant reciprocity residual anywhere in the saved 86
blocks is about `2.14741e-5`, in F-F blocks. Onsite Sbar is symmetric to `1e-14`.
The original main-pair bare-H defect is reproduced in all four domains with the
same structure-constant input.

`shldch` solves central columns of a finite screening problem for each independent
centre and retains a restricted vector set. This provides a concrete convergence
question. The measurements do not distinguish that approximation from a defect
in its construction/mapping. No screening radius was changed in this audit.

## Matrix-level GF and trace evidence

The following initial numerical measurements use the first original directed
bond, the unchanged 2000-channel calculation, and 64 sampled full matrices.
Absolute GF residuals are in the code's GF units. These are `||MGM-sign*G||`,
not merely norms of G:

| ij sector | Maximum element | Maximum Frobenius norm | Maximum relative Frobenius residual |
| --- | ---: | ---: | ---: |
| G00 | 2.53183e-8 | 6.46713e-8 | 3.67028e-7 |
| G01 | 1.72725e-6 | 3.76882e-6 | 2.28222e-1 |
| Gx0 | 2.48771e-8 | 6.38450e-8 | 8.92454e-6 |
| Gx1 | 1.54271e-6 | 3.67176e-6 | 1.76376e-2 |
| Gy0 | 1.33168e-8 | 2.87994e-8 | 8.10561e-4 |
| Gy1 | 1.66757e-8 | 4.63387e-8 | 1.61579e-4 |
| Gz0 | 2.21762e-8 | 5.54928e-8 | 1.76671e-3 |
| Gz1 | 1.54942e-8 | 3.07252e-8 | 1.03767e-4 |

In these numerical tables, “forbidden” means forbidden under the derived
**Hermitian-H retarded-transpose** selection rules. Leakage is not, by itself,
proof that antiunitary magnetic mirror covariance of a non-Hermitian H is lost.
Maxima in different columns need not occur at the same energy. The raw spinor
mirror-transpose comparison already fails: maximum element `1.73792e-6`, maximum
Frobenius `6.63247e-6`. Raw Pauli decomposition and the implemented R-sector
projection agree with independent Python calculations **exactly** at dump
precision. Delta mirror commutators are exactly zero. Independent CC/SC
contractions agree to maximum absolute residual `2.78e-17`.

At E=-0.3475763306 Ry (sample 577), SCx is `6.2716205008e-8` before integration:

| First / second factor block | SCx contribution |
| --- | ---: |
| allowed G00 / allowed Gx1 | 0 |
| allowed G00 / forbidden Gx1 | 6.2717817267e-8 |
| forbidden G00 / allowed Gx1 | -1.6122584997e-12 |
| forbidden / forbidden | 0 |

The largest same-parity Gx1_ji element is `(9,9)`, d-d, magnitude
`2.62868988e-7`, with complex value approximately
`2.62778756e-7 + i*6.88694965e-9`. The d(+)->d(+) trace contribution is
`5.9933560054e-8`.

At E=-0.4756234050 Ry (sample 513), CCx is `-9.0748963122e-9`:
`-9.0635349122e-9` comes from forbidden G01 with allowed Gx0;
`-1.1361399977e-11` comes from allowed G01 with forbidden Gx0.
The leading angular contributions are s(+)->d(+) and d(+)->s(+), approximately
`-2.82904e-9` and `-2.81501e-9`.

At sample 577, compare matched ji matrices between [100] and [-100]:

| Orbital block | Plus-domain norm | Residual if even | Residual if odd |
| --- | ---: | ---: | ---: |
| G01, same parity (forbidden) | 6.57201e-7 | 2.07854e-11 | 1.31441e-6 |
| G01, opposite parity (allowed) | 5.19073e-2 | 1.03815e-1 | 6.55980e-9 |
| Gx1, same parity (forbidden) | 5.32008e-7 | 1.06402e-6 | 1.64671e-11 |
| Gx1, opposite parity (allowed) | 4.33685e-2 | 7.11113e-9 | 8.67371e-2 |
| Gx0, same parity (allowed) | 4.29157 | 8.58315 | 3.29975e-7 |

Thus the dominant allowed matrices follow nominal T parity while the leaking
blocks have the opposite parity. Their products explain both observed odd CCx
and odd SCx before any energy integration. No spin-axis output rotation is involved. Delta's inter-domain difference is
small but nonzero: maximum Frobenius differences are `3.31e-8` and `1.38e-6`
for sites i/j across samples. Its mirror commutator is still zero, and a diagonal
Delta difference cannot create a trace between exact opposite-parity blocks.

## Other operations and execution status

The analyzer uses orbital C2x derived from Cartesian polynomial parities and V,
`C2z=I*M`, orbital inversion `(-1)^l`, and the corresponding spinor transformations.
C2z*T uses the full transpose and the reverse raw GF; inversion and C2x are
unitary. Directed bond maps use types and displacements, with the actual cluster
indices retained. For the first four original pairs `(1,2)`, `(1,592)`, `(1,611)`,
`(1,612)`:

| Operation | Mapped target indices in the same order |
| --- | --- |
| C2x | 592, 2, 612, 611 |
| C2z*T | 611, 612, 2, 592 |
| inversion | 611, 612, 2, 592 |

The first index is 1 in each mapped pair. These are geometric maps, not assumed
list-order relations.
Independent reconstruction of **all 86 bare hopping/radial blocks**, including
F1<->F2 source/destination maps, gives direct-H maximum covariance residuals:

| Operation | Maximum element residual (Ry) | Missing mapped blocks |
| --- | ---: | ---: |
| C2x | 2.33670e-10 | 0 |
| C2z*T (direct antiunitary H test) | 1.15186e-15 | 0 |
| inversion | 2.33670e-10 | 0 |
| m_z*T (direct antiunitary H test) | 2.33670e-10 | 0 |

For the 18 Mn-Mn bare blocks, antiunitary **transpose** comparisons reach
`7.17862e-8` Ry, while direct-H covariance remains at `5.06e-11` Ry or better.
This independently separates spatial covariance from reciprocal Hermiticity.
The full effective-H action is checked through its captured columns rather than
asserted identical to this bare-block reconstruction. The detailed analyzer separately reports missing mapped
pairs, absolute/relative raw covariance residuals at every sampled energy, and
matrix-level bare-T residuals.

Execution log and full-run numerical summary are completed below after the
24-pair calculations. Initial validation: 15 synthetic/source/analyzer tests pass;
unmodified baseline and instrumented serial first-pair [100]/[-100] runs succeed;
`dijso.out` and `dijparts.out` are identical at printed precision for each domain.
GNU Fortran 13.3, OpenBLAS 0.3.26 and OpenMPI 4.1.6 are available through a local
toolchain. An MPI two-rank/two-pair smoke calculation succeeds, with separate rank
and global-pair diagnostic files. Baseline serial smoke runtimes were about 93 s
per domain; instrumented runtimes about 91-95 s under that load. Timing is not
an overhead benchmark.

## Decision tree and remaining uncertainty

1. **Do raw screened constants obey reciprocal transpose?** No. This is the first
   observed failure. Binary input and independent reconstruction locate it before
   exchange formulas or GF projections.
2. **Does it reach the actual Hamiltonian?** Yes. Captured h/HoH/effective-H defects
   and forbidden scalar-01/x1 blocks match the reconstructed h quantitatively.
3. **Can magnetic mirror invariance alone justify the proposed retarded test?**
   No. The actual H approximately preserves antiunitary mirror covariance while
   failing Hermiticity at a larger scale. The transpose step's assumption fails.
4. **Where do the nonzero contractions originate?** Same-parity G01 and Gx1
   leakage, predominantly through the block products listed above. Projection,
   Delta's mirror commutation, and matrix multiplication do not create it.
5. **Why is Dx odd?** The measured leakage has reversed nominal bare-T parity;
   multiplying it by the dominant allowed partner gives odd CCx and SCx. Their
   original energy integration then reports odd Dx.
6. **Is this a physical global Cartesian pairwise DMI or a higher-order term?**
   This calculation does not establish either interpretation. First establish a
   reciprocal/Hermitian Hamiltonian and the symmetry of the exchange Hessian.
7. **What further work is needed?** A separate controlled screening-convergence
   study and full-H Hermiticity/covariance checks can distinguish approximation,
   finite-neighborhood effects and implementation defects. Do not silently
   symmetrize or increase cutoffs as part of this diagnostic PR.

No additional original SCF input is required to reproduce the supplied jobs.
To establish screening convergence, request the original structure-constant
construction provenance (source revision, screening/cutoff settings, generation
logs) and authorize new screening-only comparisons. The code revision that first
produced the supplied outputs is not recorded in the archive. Complete matrix
JSONs and raw dumps are private artifacts and must accompany any independent
numerical review through an appropriate private channel, not this public PR.

## Source-code references

| Stage | Source |
| --- | --- |
| Orbital transformation V | [math.f90:1508](source/math.f90#L1508) |
| Numerical L operators | [math.f90:133](source/math.f90#L133) |
| R transpose and GF sectors | [green.f90:386](source/green.f90#L386) |
| Raw spinor GF reconstruction and Pauli components | [green.f90:425](source/green.f90#L425) |
| Chebyshev GF kernel | [green.f90:890](source/green.f90#L890) |
| Diagonal Delta | [symbolic_atom.f90:241](source/symbolic_atom.f90#L241) |
| Hopping widths and channel expansion | [symbolic_atom.f90:163](source/symbolic_atom.f90#L163) |
| SOC construction | [hamiltonian.f90:1707](source/hamiltonian.f90#L1707) |
| Spin-block Hamiltonian assembly | [hamiltonian.f90:2017](source/hamiltonian.f90#L2017) |
| Bare spinor hopping | [hamiltonian.f90:2689](source/hamiltonian.f90#L2689) |
| Structure-constant extraction/transpose | [hamiltonian.f90:2835](source/hamiltonian.f90#L2835) |
| Finite screening and retained blocks | [lattice.f90:2425](source/lattice.f90#L2425) |
| Actual first HoH/SOC action | [recursion.f90:2317](source/recursion.f90#L2317) |
| Four intersite polarization seeds | [recursion.f90:2448](source/recursion.f90#L2448) |
| Accelerated moment construction (uses Hermitian-H identities) | [recursion.f90:2678](source/recursion.f90#L2678) |
| CC/SC multiplication | [exchange.f90:917](source/exchange.f90#L917) |
| CC/SC DMI traces and factor of two | [exchange.f90:1176](source/exchange.f90#L1176) |
| Unchanged Simpson integration and mRy conversion | [exchange.f90:1295](source/exchange.f90#L1295) |

Reproduction commands and dump schema: [diagnostic usage](diagnostics/gf_symmetry/README.md).
