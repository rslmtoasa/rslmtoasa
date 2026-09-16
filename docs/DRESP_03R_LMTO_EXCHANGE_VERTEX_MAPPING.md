# DRESP-03R — Orthogonal-Hamiltonian ↔ LMTO Exchange-Vertex Mapping

Status: **partial PASS for the offsite canonical/auxiliary contraction map;
BLOCKED for the physical path-operator → orthogonal-resolvent map and the
finite-Hamiltonian local-torque bridge**.

This report is a representation remediation of DRESP-03.  The earlier finite
K/L contraction is not reopened and its scalar audit is not reinterpreted.
The new fixture only proves live LMTO algebra and the local derivative of the
live magnetic Hamiltonian builder.  No later TDDFT physics is added.

## 1. Exact live `H_orth`

The repository does not assemble a textbook single-spin
`C + sqrt(Delta) S (1-gamma S)^(-1) sqrt(Delta)` expression directly.  Its
live Hamiltonian is assembled in two layers.

First, `lmto_bond_value` in `source/lmto_magnetic_tangent.f90:18-53`
constructs the four Pauli components of every directed block.  With
`hhh=hhh_ij`, `w0_i=potential%wx0`, `w1_i=potential%wx1`, and local unit
moments `m_i,m_j`, the live equations are

\[
 h^0_{ij}=w^0_i hhh_{ij}w^0_j+w^1_i hhh_{ij}w^1_j(m_i\cdot m_j),
\]
\[
 h^a_{ij}=w^1_i hhh_{ij}w^0_jm_{i,a}
          +w^0_i hhh_{ij}w^1_jm_{j,a}
          +i w^1_i hhh_{ij}w^1_j(m_i\times m_j)_a.
\]
For an onsite block the live builder adds diagonal `cx0` and `cx1` terms.
The four components are embedded as `H0+Hz`, `H0-Hz`, `Hx-iHy`, and
`Hx+iHy` by `source/hamiltonian_build.f90:1254-1265`.

Define the reciprocal sums actually used by the assembler,

\[
 B(k)=\sum_R ee(R)e^{i2\pi k\cdot R},\qquad
 Q(k)=\sum_R eeo(R)e^{i2\pi k\cdot R}.
\]

For first order, the live reciprocal path is

\[
 H_{\rm live}^{(1)}(k)=B(k).
\]

For second-order `hoh`, `source/reciprocal_fourier.f90:149-182` assembles

\[
 H_{\rm live}^{(2)}(k)=B(k)-Q(k)B(k)
 +\operatorname{diag}_{\rm site}[enim+lsham]+H_{cc}(k),
\]

where `eeo(R)=ee(R) obarm_j` is built at
`source/hamiltonian_build.f90:1289-1305`.  `H_cc` is present only when the
optional combined correction is enabled.  `enim` carries the onsite
`cex0/cex1` spin dependence (`source/hamiltonian_build.f90:1162-1203`).
The normal `ham_only` reciprocal Green path solves this assembled matrix with
`S=I`; it does not solve a generalized eigenproblem.

Therefore the exact live `H_orth` for this audit is the assembled
`H_live^(1)` or `H_live^(2)` above, not a guessed textbook expression.

## 2. Exact live `P`, `S`, and `g`

The live potential-function routine is
`source/symbolic_atom.f90:477-506`.  In the live orbital ordering
`l*l+m`, with spin-up first and spin-down second, it constructs

\[
 P_l^\sigma(E)=\frac{E-[c_l^\sigma+vmad]}{dele_l^{\sigma 2}}.
\]

The live equivalents are:

| textbook object | live repository object | precise role |
|---|---|---|
| `C^sigma` | `potential%c(l,s)+potential%vmad` | center used by `p_matrix` and `d_matrix` |
| `Delta^sigma` | `potential%dele(l,s)^2` | denominator of the live `P` function |
| `sqrt(Delta^sigma)` | `potential%dele(l,s)` | endpoint factor in `green%auxiliary_gij` |
| transformed `C`, width, `O` | `center_band`, `width_band`, `obar` and their `cx/wx/obx` copies | tight-binding/HOH Hamiltonian construction |
| `gamma^sigma` in the auxiliary screened route | `potential%qpar(l,s)` | explicitly identified as orthogonal screening at `exchange.f90:487-490` |
| canonical screening | `screening_can=0` in the three-index auxiliary route | canonical endpoint of the live transform |
| optional strux screening | `potential%screening_alpha`, `potential%screening_sigma`; lattice copies have the same names | structure-constant backend controls |
| `S` | `lattice%sbar` screened structure-constant blocks, with `sbarvec` geometry | source structure data; reciprocal H consumes the already-built `ee` blocks |

The source contains no live array literally named `gamma`.  The
`qpar`/gamma identification is a source comment and an actual input to the
auxiliary screening route, not a claim about `obar`.

The mathematical path operator required by the question would be

\[
 g^\sigma(E)=[P^\sigma(E)-S]^{-1}.
\]

The current production tree does **not** expose a routine that constructs
this object from `p_matrix` and `lattice%sbar` and then feeds it into the
native reciprocal/recursive Green path.  `green%gij` is generated from the
assembled H by recursion or direct/Lehmann resolvents.  This is the first
blocking boundary.

## 3. Physical-GF transformation

The live coefficient-space endpoint routine is
`source/green.f90:346-396`:

\[
 g^{aux}_{ij}(E)=D_iG_{ij}(E)D_j,
 \qquad D_i=\operatorname{diag}[dele_i(l,\sigma)].
\]

The routine comments call its input a physical site-resolved GF and its
output an auxiliary GF.  The implementation is nevertheless only the
explicit endpoint multiplication above; it does not construct the local
`lambda(E)` term or a live `mu(E)` array.  The separate screened transform at
`source/green.f90:409-473` applies

\[
 g^{out}_{ij}=R_i g^{in}_{ij}R_j+\delta_{ij}A_i,
 \quad R=P^{in}/P^{out},
\]

with the additive onsite term `A_i` implemented as `pmat_resc3`.  For an
offsite pair the additive term is zero.

Thus the repository-equivalent offsite identity is the tested endpoint map
`g_aux = D_i G_H D_j`.  The stronger statement

\[
 G_H=\lambda+\mu[P-S]^{-1}\mu
\]

cannot be established from live code because the corresponding `g=(P-S)^-1`,
`lambda`, and `mu` construction is absent from this response path.

## 4. Numerical `G` proof and its limit

The available numerical proof is the coefficient-space identity, not a proof
of the missing physical/path-operator identity.

`UnitKspaceGFValidation` compares the orthogonal Lehmann sum

\[
 G_H(k,z)=\sum_n\frac{|\psi_{nk}\rangle\langle\psi_{nk}|}{z-\epsilon_{nk}}
\]

with `[zI-H(k)]^-1` on five k points and eight complex energies:

| check | observed maximum error |
|---|---:|
| Lehmann versus direct inverse | `7.6404e-15` |
| resolvent residual | `1.1411e-15` |
| retarded/advanced conjugacy | `2.2204e-16` |

`UnitLrRsGfRepresentation` separately verifies the live endpoint and screened
auxiliary algebra on deterministic blocks, including an offsite phase bridge
error of `1.13e-15` and endpoint/screening bridge error of `5.27e-16`.
Those checks do not create a native `P-S` inverse and therefore do not upgrade
the result to a physical-GF proof.

## 5. Native `d_matrix(E)`

The canonical native LKAG call is
`source/exchange.f90:1204-1210`, and its source definition is
`source/symbolic_atom.f90:313-341`:

\[
 d_l(E)=
 \frac{c_l^\downarrow dele_l^{\uparrow2}
       -c_l^\uparrow dele_l^{\downarrow2}
       +(dele_l^{\downarrow2}-dele_l^{\uparrow2})E}
      {dele_l^\uparrow dele_l^\downarrow},
\]

where both `c` values include `vmad`.  Because the live `P` function divides
by `dele^2`, the exact algebraic relation is

\[
 d_l(E)=dele_l^\uparrow dele_l^\downarrow
        [P_l^\uparrow(E)-P_l^\downarrow(E)].
\]

It is therefore not raw `DeltaP`.  It is an energy-dependent, endpoint-scaled
`DeltaP` with the source spin ordering up-minus-down.  It is diagonal in the
live orbital ordering and carries the code energy unit (the `E`, `c`, and
`vmad` arrays are in Ry; no extra conversion occurs in `d_matrix`).

The controlled test measured:

| comparison | result |
|---|---:|
| `d_matrix - dele_up*dele_down*DeltaP` | `1.4814e-7` maximum |
| `d_matrix - raw DeltaP` | `8.5571e-1` maximum |

The small first residual is from the source’s default-kind `cmplx` calls in
`d_matrix` (`c`/`dele` are converted without `kind=rp`), while `p_matrix`
uses explicit `kind=rp`.  No precision behavior was silently corrected in
this audit.

`d_matrix` has no explicit combined-correction, Hubbard, SOC, or screening
argument.  Any such physics can enter only through the Green function and
Hamiltonian state consumed by the surrounding exchange path; it is not part
of the vertex formula itself.

## 6. Local rotation in both representations

For an independent local rotation of site `i`, write

\[
 \delta m_i=\delta\theta_i\times m_i,
 \qquad \delta m_j=0\quad(i\ne j).
\]

Differentiating the live bond formula gives, with
`a=w0_i hhh w0_j`, `b=w1_i hhh w1_j`,
`c=w1_i hhh w0_j`, `d=w0_i hhh w1_j`, and `e=w1_i hhh w1_j`,

\[
 \delta h^0_{ij}=b[\delta m_i\cdot m_j+m_i\cdot\delta m_j],
\]
\[
 \delta h^a_{ij}=c\,\delta m_{i,a}+d\,\delta m_{j,a}
 +ie[(\delta m_i\times m_j)+(m_i\times\delta m_j)]_a.
\]

For an onsite block both occurrences of the site moment rotate, and the
explicit onsite term adds

\[
 \delta h^a_{ii}\supset c1_i\,\delta m_{i,a}.
\]

The full orthogonal torque is the spinor embedding of these derivatives in
every block touching site `i`.  In second order it also differentiates the
assembled `Q(k)B(k)` product, and any enabled `enim`/correction dependency.
This is the invariant finite-H object:

\[
 T_i^H=\left.\frac{\partial H_{live}(\{m\})}{\partial\theta_i}\right|_0,
\]

not a site block of `H_up-H_down`.

The native P-side route currently exposes diagonal `P(E)` and `DeltaP(E)`;
it has no production local-rotation derivative or P-space force-theorem
response to compare against.  Consequently the same-rotation P/H response
cannot yet be certified.

## 7. Transformed vertex and complete-contraction result

For an offsite collinear pair, the live endpoint routine gives, channel by
channel,

\[
 g^{aux,\sigma}_{ij}=D_i^\sigma G^\sigma_{H,ij}D_j^\sigma.
\]

Substitution into the auxiliary pair product and cyclicity of the trace give

\[
 \operatorname{Tr}[\Delta P_i g^{aux,\uparrow}_{ij}
                    \Delta P_j g^{aux,\downarrow}_{ji}]
 =
 \operatorname{Tr}[d_iG^\uparrow_{H,ij}
                    d_jG^\downarrow_{H,ji}],
\]

where `d_i=D_i^up*D_i^down*DeltaP_i` and similarly for `j`.  The reverse
spin ordering supplies the partner term.  This proves the allowed outcome
“contraction equivalent, vertices not identical” for the live offsite
canonical/auxiliary pair algebra, provided both routes consume the same
`green%gij` H-resolvent blocks.

The general candidate

\[
 T_i^H(E)=(\mu_i^\downarrow)^{-1}\Delta P_i(E)(\mu_i^\uparrow)^{-1}
\]

and its opposite-spin partner follows only if the missing
`G=lambda+mu*g*mu` contract is first established.  It is not installed as a
production vertex here.  The endpoint proof above is narrower and does not
prove the full path-operator representation or the local-rotation derivative
map.

## 8. Pair-integrand comparison

The native complete J integrand is the one in
`source/exchange.f90:640-666`:

\[
 \operatorname{Tr}\{d_iG^0_{ij}d_jG^0_{ji}
 -d_iG^x_{ij}d_jG^x_{ji}
 -d_iG^y_{ij}d_jG^y_{ji}
 -d_iG^z_{ij}d_jG^z_{ji}\}.
\]

The auxiliary route forms `DeltaP`, endpoint-scales `green%gij`, and applies
the screened transform at `source/exchange.f90:283-307`.  For offsite,
collinear blocks the algebra in Section 7 proves equality of the complete
contraction, not equality of isolated vertices.  The onsite auxiliary branch
has an additive local term (`source/exchange.f90:309-316`) and is not covered
by the offsite identity.

No native-P versus transformed-H energy-resolved pair table is claimed.  Such
a table would require the missing live `g=(P-S)^-1`/physical-GF construction
and a P-side local rotation response.  The bcc-Fe energy-resolved comparison
is therefore blocked rather than hidden behind an integrated number.

## 9. Contour-integrated comparison

The mature native route integrates `imtrace9` with `simpson_f` to the Fermi
energy and writes `10^3/(4*pi)` mRy output
(`source/exchange.f90:1226-1229`).  The finite DRESP-03 service evaluates a
different explicit operator contraction with its locked spectral prefactor.

The offsite canonical/auxiliary contraction identity can be integrated with
the same contour once a common H-resolvent state is selected.  That is an
algebraic consequence of the pointwise identity; it is not yet a comparison
of the mature native `d_matrix` route to the finite-H local torque.

No pairwise Fe shell table is produced in this rung because the missing
physical/P-side local-rotation map prevents a same-state, same-vertex
comparison.  Existing native Fe references remain provenance-only evidence;
they are not mixed with the DRESP-02/03 finite snapshots.

## 10. Representation order and corrections

The rigorously matched capability subset is:

* orthogonal `ham_only`, with the same first- or second-order selection on
  both sides;
* the same H-resolvent Green blocks and energy mesh;
* no generalized-overlap solve (`S=I` in the reciprocal Green backend);
* no SOC, Hubbard, or combined correction for the clean exchange-vertex
  comparison; and
* the same `hoh` status, including the global reciprocal `Q(k)B(k)` term when
  second order is enabled.

The native canonical vertex itself is only `d_matrix(E)`; it does not carry
`qpar`, `screening_alpha`, `screening_sigma`, `ccor`, or an SOC correction.
The auxiliary route can transform from `qpar`/orthogonal screening to
canonical screening zero.  A route with different screening, `hoh`, CCOR,
overlap, or correction order cannot be demanded to match the clean finite
bridge exactly.

## 11. Impact on DRESP-03

The earlier implementation in `source/lr_kl_static_bridge.f90` uses an
explicit supplied operator.  Its helper
`build_ham_only_exchange_vertex`:

* takes supplied `up_blocks` and `down_blocks`;
* forms their difference;
* masks selected orbitals; and
* embeds that difference only in the local plus-spin-flip block.

It is therefore a **projected onsite H-up-minus-H-down operator**, not the
full `H_up-H_down`, not `d_matrix(E)`, and not the local derivative
`T_i^H=partial H/partial theta_i`.  It does not include the offsite derivative
blocks proven nonzero by the new fixture, nor the second-order derivative of
`Q(k)B(k)`.

The earlier DRESP-03 PASS remains valid exactly as stated: the supplied
operator-valued finite-H contraction and its independent spectral oracle
agree, and scalar compression is audited rather than fitted.  Its physical
LMTO/LKAG interpretation now has a more precise status:

**DRESP-03 finite-H vertex requires revision before it can be called the
LMTO local-rotation torque.**

No repair is made in this rung because the P-side local-rotation response and
the full physical-GF map are not yet available.  The mature
`source/exchange.f90` path is unchanged.

## 12. Final mapping verdict

| question | verdict |
|---|---|
| bare `d_matrix(E)` equals raw `DeltaP(E)`? | **NO**; it is endpoint-scaled `dele_up*dele_down*DeltaP` |
| offsite canonical and auxiliary full contractions? | **PASS**, under the live `auxiliary_gij` endpoint map and common H-resolvent blocks |
| native `g=(P-S)^(-1)` to H-resolvent physical-GF map? | **BLOCKED**; no live `g`, `lambda`, or `mu` construction is exposed |
| true finite-H local rotation vertex? | **PROVEN LOCALLY** from `lmto_bond_value`; it contains nonzero offsite blocks |
| DRESP-03 current finite vertex equals that torque? | **NO**; revision required |
| energy-resolved and contour-integrated Fe comparison? | **BLOCKED** pending the missing common representation/rotation response |

The invariant quantity identified by this audit is the complete force-theorem
pair contraction built from the same local rotation: in the offsite collinear
canonical/auxiliary subset it is unchanged by moving the endpoint `dele`
factors between `d_matrix` and the Green blocks.  The corresponding invariant
for the finite-H route is the contraction with the full derivative
`T_i^H=partial H_live/partial theta_i`, but equality of that torque response
to native P/LKAG has not been established.  The permitted final status is
therefore **BLOCKED — PHYSICAL-GF TRANSFORM UNRESOLVED**, with the narrower
offsite contraction equivalence recorded as a PASS.

### Fixture and verification

The controlled two-site, four-orbital fixture is
`tests/unit/test_dresp03r_lmto_mapping.f90`.  It uses unequal orbital widths,
spin-dependent endpoint factors, a non-diagonal complex structure block, and
non-collinear moments.  Observed results are:

| check | result |
|---|---:|
| endpoint-scaled `d_matrix` relation | `1.4814e-7` |
| raw `DeltaP` mismatch | `8.5571e-1` |
| full local-rotation derivative central error | `5.1267e-12` |
| nonzero offsite derivative | `3.9809e-2` |
| width/structure noncommutator | `6.0736e-2` |

Verification run:

```text
UnitDresp03rLmtoMapping: PASS
UnitLrKlStaticBridge: PASS
UnitLrRsGfRepresentation: PASS
UnitKspaceGFValidation: PASS
```

No Mills/Jülich U, GSR, ALSDA, Goldstone, Dyson, or projected-susceptibility
work is included.  The next rung remains forbidden until the blocked mapping
and DRESP-03 torque vertex are closed.
