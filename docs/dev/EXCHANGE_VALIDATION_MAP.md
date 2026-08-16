# VAL-06 — Exchange formulations and validation oracles

Status: audit completed 2026-08-16 against the current working tree.

This document is an inventory, not a formula change.  It records what the
public exchange methods currently do, which callers can reach them, and which
comparisons are scientifically legitimate.  A method being public does not
mean that it is wired into a production driver or that its output is already
validated.

## Executive conclusions

The current production exchange workflow has three distinct consumers:

1. `calculate_exchange` evaluates the combined scalar Heisenberg term, DMI
   vector, and a full (3\times3) anisotropic term from the canonical Pauli
   decomposition of the intersite Green functions.  It is called for the
   ordinary `post_processing='exchange'` route.
2. `calculate_exchange_twoindex` evaluates a two-index perturbative
   decomposition (first/second order and charge/spin subparts) from the
   `g00/g01/g{x,y,z}{0,1}` arrays.  It is called immediately after
   `calculate_exchange`, but its first/second-order outputs are not the same
   observable as the combined full result.
3. `calculate_exchange_gauss_legendre` evaluates the same low-level matrix
   kernels on the 64-point Fermi-eta ladder.  It is called by
   `post_processing='exchange_p2rs'` after PAOFLOW Hamiltonian import.

`calculate_jij`, `calculate_jij_auxgreen`, `calculate_jijk`, and
`calculate_exchange_rs2pao` are public methods with no current production
callers.  `calculate_dij` and `calculate_aij` are empty public stubs.  The
Gilbert damping method is now reachable through `calculation%do_damping`; the
moment-of-inertia method remains public but has no caller.

The strongest already-established route oracle is:

```text
Lehmann GF == Dyson GF at Sigma = 0       tight, solver precision
real-space GF ~= Lehmann GF                only a documented convergence/broadening envelope
```

The second relationship is deliberately not an equality at fixed finite
broadening.  The existing bcc-Fe triad measures the exchange observable from
different GF producers, not a universal normalization constant.

## Source and data contract

The public type and all exchange kernels are in
[`source/exchange.f90`](../../source/exchange.f90), with damping and inertia
in [`source/exchange_dynamics.f90`](../../source/exchange_dynamics.f90).
The exchange object is constructed from a `bands` object and holds pointers to
the common `green`, `lattice`, `energy`, `symbolic_atom`, `hamiltonian`, and
`control` objects.

For the canonical exchange kernels, let

* (Δ_i(E)) be the (9\times9) `d_matrix` returned for atom (i);
* (G^0_{ij}) denote `Ginmag` and (G^a_{ij}), (a=x,y,z), denote
  `Gix/Giy/Giz`; and
* the corresponding (ji) quantities be `Gjnmag` and `Gjx/Gjy/Gjz`.

The active Simpson path forms, at each real energy, the matrix products

\[
 K_J = \operatorname{ImTr}\left[\Delta_iG^0_{ij}\Delta_jG^0_{ji}
       -\sum_a\Delta_iG^a_{ij}\Delta_jG^a_{ji}\right],
\]

\[
 K_D^a = \operatorname{ReTr}\left[\Delta_iG^0_{ij}\Delta_jG^a_{ji}
       -\Delta_jG^0_{ji}\Delta_iG^a_{ij}\right],
\]

and

\[
 K_A^{ab}=\operatorname{ImTr}\left[
 \tfrac12(\Delta_iG^a_{ij}\Delta_jG^b_{ji}
 +\Delta_jG^a_{ji}\Delta_iG^b_{ij})\right].
\]

These are the source-level conventions in `dGdG_Jnc`, `dGdG_Dnc`, and
`dGdG_Anc`; the energy integral is Simpson integration to `en%fermi`,
followed by the source output scale `10^3/(4*pi)`.  The literature convention
and the source convention must be kept separate when comparing signs and
units.

The real-space intersite producer fills these arrays from the four block
continued-fraction Green functions in
[`source/green_block.f90`](../../source/green_block.f90).  The reciprocal
producer fills the same canonical arrays from the Lehmann or Dyson k-space
engine in [`source/reciprocal_green.f90`](../../source/reciprocal_green.f90).
`calculate_intersite_gf_twoindex` derives the two-index families from the
canonical Pauli families; it is an algebraic change of representation, not a
new GF solver.

## Public formulation inventory

### 1. `calculate_jij`

* **Quantity:** scalar intersite Heisenberg exchange (J_{ij}).
* **Formulation:** the `dGdG_Jnc` matrix kernel above, followed by
  `imtrace9` and Simpson integration.  It is the scalar J-only version of the
  J part of `calculate_exchange`.
* **GF/input objects:** canonical `Ginmag/Gjnmag` and `Gix/Giy/Giz` plus the
  corresponding `j` arrays; `d_matrix` from each symbolic atom; `en%ene`,
  `en%fermi`, and `en%nv1`.
* **Spin/SOC assumptions:** uses the global-frame Pauli decomposition and
  spin-split Δ matrices.  It can consume the canonical arrays produced by a
  noncollinear/SOC-capable producer, but the method itself does not add an
  SOC operator or rotate frames.
* **Outputs:** only a diagnostic line on stdout, scaled as
  `jij*1e3/(4*pi)`; no file or per-pair array is retained beyond the scalar
  member.
* **Callers/tests:** no production caller.  No test invokes this method by
  name.  The ordinary exchange fixtures exercise the same kernel through
  `calculate_exchange`, not this entry point.
* **Documented origin:** the module-level documentation identifies the LKAG
  formalism; the auxiliary-GF method, rather than this J-only wrapper, cites
  Pajda et al., PRB **64**, 174402 (2001).
* **Legitimate oracle:** exact algebraic equality to the J column produced by
  `calculate_exchange` is valid when the same canonical GF arrays, Δ matrices,
  energy mesh, and integration settings are supplied.  It is a useful compact
  kernel oracle.  It is not an independent physical calculation.

### 2. `calculate_dij`

* **Quantity/formulation:** none in the current source; the method body is an
  empty stub.
* **Inputs/assumptions/outputs:** no arrays are read and no output is written.
* **Callers/tests:** no callers and no tests.
* **Origin/oracle:** no implemented formulation and therefore no validation
  oracle.  DMI validation must target the DMI branch of
  `calculate_exchange` or the corresponding two-index/eta route, not this
  name.

### 3. `calculate_aij`

* **Quantity/formulation:** none in the current source; the method body is an
  empty stub.
* **Inputs/assumptions/outputs:** no arrays are read and no output is written.
* **Callers/tests:** no callers and no tests.
* **Origin/oracle:** no implemented formulation and therefore no validation
  oracle.  The implemented anisotropic quantity is the `aij` branch of
  `calculate_exchange` and its variants.

### 4. `calculate_exchange`

* **Quantity:** the combined pair interaction: scalar (J_{ij}), DMI vector
  (mathbf D_{ij}), and a full (A_{ij}^{ab}) tensor.
* **Formulation:** uses `dGdG_Jnc`, `dGdG_Dnc`, and `dGdG_Anc` with the
  canonical arrays.  J uses `ImTr`, D uses `ReTr`, and A uses `ImTr` as shown
  above.  All three are Simpson-integrated to the Fermi energy.
* **GF/input objects:** `green%Ginmag/Gjnmag/Gi{x,y,z}/Gj{x,y,z}`;
  `symbolic_atom% d_matrix`; lattice pair list `ijpair`; real-energy mesh;
  optional constraint state.  The preceding caller first fills these arrays
  from recursion, Lehmann, or Dyson.
* **Spin/SOC assumptions:** the exchange kernel is spinor/Pauli decomposed and
  stores intersite blocks in the global spin frame.  SOC is present only if it
  was included by the Hamiltonian/GF producer; DMI and anisotropic terms have
  no nonzero guarantee without the relevant SOC and symmetry breaking.
* **Outputs:**
  `jij.out` (species, bond vector, J, distance), `dij.out` (bond vector and
  D), `aij.out` (bond vector and nine A entries), stdout diagnostics, and the
  assembled `J_tensor` printed to stdout.  `jtens.out` is opened and closed,
  but the current routine does not write a record to it.  The combined tensor
  assembled in memory is

  \[
  T_{ij}=J_{ij}I + D_{ij}^{\rm skew}+A_{ij},
  \quad D^{\rm skew}=\begin{pmatrix}0&D_z&-D_y\\-D_z&0&D_x\\D_y&-D_x&0\end{pmatrix}.
  \]

* **Callers/tests:** `calculation::post_processing_exchange` calls it when
  `njij>0` and `njijk==0`, after the selected GF route and before
  `calculate_exchange_twoindex`.  `Example_exchange_bccFe` and its HOH sibling
  exercise it; the bcc-Fe route triad exercises its J output for recursion,
  Lehmann, and Dyson.  Existing checks assert J and D columns, but do not
  assert `aij.out` or the complete assembled tensor.
* **Documented origin:** module documentation identifies LKAG; the DMI
  reference in the source documentation is Sandratskii, PRB **96**, 024450
  (2017).  The source itself does not attach a separate equation number to
  the combined kernel.
* **Legitimate oracle:** exact equality to `calculate_jij` for scalar J under
  identical inputs.  The D vector and `D^skew` assembly are exact source-level
  identities.  Do not infer D solely from the antisymmetric part of the
  printed full tensor unless (A) has first been shown symmetric: the source
  stores a full, not explicitly symmetrized, `aij`.

### 5. `calculate_exchange_twoindex`

* **Quantity:** first-order and second-order exchange decompositions, with
  charge/spin density/current parts for J, D, and A.
* **Formulation:** starts from the two-index families generated by
  `green%calculate_intersite_gf_twoindex`:
  `g00`, `g01`, `gx0/gy0/gz0`, and `gx1/gy1/gz1` for both bond directions.
  The source combines their traces as:

  ```text
  J_SO = CD - SD + CC - SC
  J_FO = CD + SD - CC - SC
  D_SO = 2*(SC + CC)       D_FO = 2*(SC - CC)
  A_SO = SD + SC           A_FO = -SD + SC
  ```

  Here the labels are the source comments/variables (`jcd`, `jsd`, `jcc`,
  `jsc`, `dcc`, `dsc`, `isd`, `isc`).  This is a perturbative decomposition,
  not a second implementation of the combined full-(J/D/A) observable.
* **GF/input objects:** the two-index arrays above, Δ matrices from
  `d_matrix`, real-energy mesh, pair list, and MPI pair partitioning.
* **Spin/SOC assumptions:** relies on the two-index orbital-reversal and spin
  decomposition in `green_block.f90`; the first/second-order interpretation is
  tied to the source's SOC/order convention.  It is not valid to compare an
  SO or FO file to full `jij.out` without stating the perturbative order and
  convention.
* **Outputs:** `jijso.out`, `jijfo.out`, `dijso.out`, `dijfo.out`,
  `aijso.out`, `aijfo.out`, and the parts files `jijparts.out`,
  `dijparts.out`, `aijparts.out`.  `jtensso.out` and `jtensfo.out` are opened
  and closed but are not populated by current write statements.
* **Callers/tests:** called immediately after `calculate_exchange` in
  `post_processing_exchange`.  The ordinary exchange fixtures therefore run
  it, but their reference checks do not inspect any of these files.  No unit
  test asserts the component recombinations.
* **Documented origin:** the source gives the order labels and the module
  documentation describes the LKAG/noncollinear formalism; no more specific
  literature equation is cited for this two-index split.
* **Legitimate oracle:** the two J, two D, and two A recombination equations
  above are exact internal algebraic identities.  A full result versus an SO
  or FO result is intentionally not an equality.  A future independent check
  should compare the two-index quantities against a separately constructed
  perturbative reference, not merely sum arbitrary output files and call the
  result full exchange.

### 6. `calculate_jij_auxgreen`

* **Quantity:** auxiliary-GF LKAG pair tensor `jij_aux(1:9)` and on-site
  `jij00_aux` for (i=j).  The source prints a derived `Dij_zz_aux` from
  `0.5*(jij_aux(xy)-jij_aux(yx))`.
* **Formulation:** builds spin-channel potential matrices (P_i^σ(E)),
  ΔP = (P^↑-P^↓), converts the physical GF with
  `green%auxiliary_gij`, and evaluates nine angular projections of products of
  auxiliary up/up and down/down blocks.  The on-site branch uses its separate
  `J00` expression.
* **GF/input objects:** `green%gij/gji`, `symbolic_atom%p_matrix`, atom
  `lmax`, `auxiliary_gij`, and the energy mesh.  The source comment says this
  works better with `hoh` because it supposes an orthogonal representation.
* **Spin/SOC assumptions:** explicitly channel-resolved up/down P matrices;
  therefore the natural scope is a collinear spin-channel calculation in a
  matched orthogonal representation.  It is not automatically the same
  object as the full spinor/SOC tensor from `calculate_exchange`.
* **Outputs:** updates `jij_aux` or `jij00_aux` and prints to stdout.  Contrary
  to the older Sphinx reference, the current routine does not write an
  `exchange.out` file.
* **Callers/tests:** no production caller and no test.
* **Documented origin:** source comment: Eq. 3 of Pajda, Kudrnovský, Turek,
  Drchal, and Bruno, PRB **64**, 174402 (2001).
* **Legitimate oracle:** in a collinear, SOC-free, orthogonal-representation
  limit, with identical GF, P/Δ convention, energy mesh, and normalization,
  the appropriate scalar auxiliary projection may be compared with the
  canonical LKAG J.  This is conditional, not a general equality for the
  nine auxiliary components or for arbitrary SOC/noncollinear data.

### 7. `calculate_jijk`

* **Quantity:** a nine-component spin-lattice coupling tensor for a trio
  ((i,j,k)), with a specified displacement direction for atom (k).  The
  source reports units as mRy/Å after its `wav`/unit conversion.
* **Formulation:** constructs orthogonal and canonical P matrices, obtains
  six auxiliary GFs for (ij,ik,jk) and reverse directions, constructs the
  displacement matrix `umatrix_k`, and evaluates the eight-term product used
  in `int_all`.  It Simpson-integrates each component to the Fermi energy.
* **GF/input objects:** three pair slots per trio in `green%gij/gji`, P matrices,
  screening constants, `udisp_matrix`, `lattice%ijktrio(:,1:6)`, `wav`, and
  the energy mesh.
* **Spin/SOC assumptions:** channel-resolved up/down P matrices and the
  orthogonal-to-canonical transformations; displacement is normalized before
  use.  No active route currently establishes a general SOC/noncollinear
  interpretation for this method.
* **Outputs:** updates `jijk(1:9)` and prints the tensor and displacement
  direction to stdout; no file is written by the method.
* **Callers/tests:** no caller.  The lattice lifecycle can expand a trio into
  three pair records, but `post_processing_exchange` explicitly runs the
  pair consumers only when `njijk==0`, and no driver calls `calculate_jijk`.
  No test exercises it.
* **Documented origin:** source comments call it spin-lattice interaction;
  the Sphinx reference describes it as a derivative of pair exchange but the
  source has no specific citation or equation number.
* **Legitimate oracle:** a controlled finite-difference derivative of a
  separately rebuilt pair-(J) calculation under (+u_k) and (-u_k), with
  the same representation and unit convention, is the appropriate independent
  oracle.  `jijk` is not expected to equal a pair J value.

### 8. `calculate_exchange_gauss_legendre`

* **Quantity:** combined J, D, and A outputs, using a 64-node quadrature on
  the Fermi-eta ladder.
* **Formulation:** calls `green%gij_eta_to_gij`, then reuses
  `dGdG_Jnc/Dnc/Anc` with 64 eta-indexed samples.  The Green producer maps
  (x\in(0,1)) to η = ((1-x)/x); this routine applies the source weight
  `w/x**2` and sums rather than calling `simpson_f`.  It also uses
  `rtrace9` for J and A and `imtrace9` for D, unlike the real-energy Simpson
  path's `imtrace9`, `rtrace9`, and `imtrace9` respectively.  That difference
  must be resolved by the contour derivation before treating the two methods
  as numerically equivalent.
* **GF/input objects:** `green%gij_eta/gji_eta` and their Pauli families,
  64-point Gauss-Legendre nodes/weights, Hamiltonian `ee` spin-up/down blocks
  for Δ, pair geometry, and MPI pair partitioning.
* **Spin/SOC assumptions:** same global-frame Pauli family convention as the
  canonical kernels, but sampled at a Fermi-centered imaginary-eta contour.
  It is not the same GF data as the full real-energy Simpson path.
* **Outputs:** `jij.out`, `dij.out`, and `aij.out`; stdout diagnostics.  The
  current source does not write `jtens.out` records.
* **Callers/tests:** `post_processing_exchange_p2rs` runs intersite moments,
  `calculate_intersite_gf_eta`, then this method.  The PAOFLOW exchange fixture
  `Example_paoflow_exchange_p2rs_bccFe` checks J and D only.  No test checks A
  or performs a same-stack Simpson-versus-Gauss comparison.
* **Documented origin:** the source uses the shared `gauss_legendre` utility;
  no separate literature equation is cited for this exchange contour.
* **Legitimate oracle:** a numerical-integration cross-check against the
  Simpson path is valid only after holding the Hamiltonian, Δ matrices, pair
  geometry, broadening/contour, and output convention fixed, and after deriving
  the (x\mapsto\eta) Jacobian and the real/imaginary trace conversion.  The
  existing PAOFLOW fixture is not such an oracle because its Hamiltonian,
  Fermi level, and energy window are not matched to the ordinary fixture.

### 9. `calculate_exchange_rs2pao`

* **Quantity:** combined J, D, and A, using the same three canonical matrix
  kernels as `calculate_exchange`.
* **Formulation:** the main difference is the source of Δ: it takes the real
  difference of the Hamiltonian `ee` up/down blocks for each species rather
  than calling `symbolic_atom%d_matrix`.  It then Simpson-integrates the same
  canonical GF products and applies the same `10^3/(4*pi)` scale.
* **GF/input objects:** canonical `green%Ginmag/Gi*` families, Hamiltonian
  `ee`, lattice pairs, and real-energy mesh.  The name indicates a real-space
  Hamiltonian imported/exported through the PAOFLOW format, but the current
  calculation driver does not call this method.
* **Spin/SOC assumptions:** inherits the canonical global-frame assumptions;
  equivalence to the ordinary method additionally requires the imported
  `ee` spin splitting to be identical to the symbolic-atom Δ convention.
* **Outputs:** writes `jij.out`, `dij.out`, and attempts to write `aij.out`;
  stdout also reports J/D/A.  The source explicitly opens units 20 and 30 but
  not unit 40 before the A write, so the A-file path is a validation blocker
  until its I/O contract is independently checked.  No production change is
  made here.
* **Callers/tests:** no callers.  The PAOFLOW post-processing driver uses
  `calculate_exchange_gauss_legendre`, not this method.  The existing PAOFLOW
  reference therefore does not exercise `calculate_exchange_rs2pao`.
* **Documented origin:** no separate literature citation; it is an import
  route around the same canonical exchange kernels.
* **Legitimate oracle:** exact equality to `calculate_exchange` is valid only
  for a controlled native-Hamiltonian export/import round trip with matched
  Δ, GF, energy mesh, and units.  A generic PAOFLOW input is not expected to
  match a native calculation.

## Damping and inertia cross-reference

These are public methods on the same `exchange` object, but they are not
exchange-formula oracles for J/D/A.

| Method | Physical quantity and formulation | Inputs/assumptions and outputs | Callers/tests and valid oracle |
|---|---|---|---|
| `calculate_gilbert_damping` | Kamberský torque-correlation tensor.  Forms (A_{ij}=G_{ij}-G_{ji}^\dagger), contracts it with `hamiltonian%tmat`, and evaluates the tensor at the energy mesh point nearest (E_F). | Requires `nsp=2` and SOC torque operators.  The source torque builder is the SOC commutator (T^a=[\sigma^a,H_{SOC}]).  Writes `damping-energy.out` and `alldampings.out`; the latter includes `0.5*(xx+yy)`. | Called from `post_processing_exchange` only when `do_damping=T`, using the same `green%gij/gji` route selected for exchange.  `Triad_triad_bccFe_damping` checks recursion/Lehmann/Dyson alpha: Lehmann and Dyson are tight; recursion is a broadening envelope.  It is not an oracle for J/D/A.  Literature documented in source docs: Miranda et al., PRM **2**, 013801 (2018), and Ebert et al., PRL **107**, 066603 (2011). |
| `calculate_moment_of_inertia` | Experimental torque-correlation/inertia quantity using anti-Hermitian and Hermitian GF parts plus a numerical second energy derivative of the Hermitian part. | Requires `nsp=2`; uses `tmat`, `gij/gji`, the energy step, and pair list.  Writes `example-real.out` and `example-imag.out`. | No caller and no test.  A future oracle would be a converged finite-difference/derivative campaign against the same production GF, not J/D/A or damping equality.  Source comment cites Sci. Rep. **7**, 931 (2017). |

## Call-site map

```text
calculation%post_processing='exchange'
  -> prepare_post_processing_stack(... use_exchange_pairs=.true.)
  -> gf_route='recursion': run_intersite_moments -> green%calculate_intersite_gf
  -> gf_route='lehmann'/'dyson': reciprocal%fill_green
  -> green%calculate_intersite_gf_twoindex
  -> exchange%calculate_exchange
  -> exchange%calculate_exchange_twoindex
  -> [do_damping] exchange%calculate_gilbert_damping

calculation%post_processing='exchange_p2rs'
  -> PAOFLOW Hamiltonian import
  -> run_intersite_moments -> green%calculate_intersite_gf_eta
  -> exchange%calculate_exchange_gauss_legendre

public but currently uncalled:
  calculate_jij, calculate_dij, calculate_aij, calculate_jij_auxgreen,
  calculate_jijk, calculate_exchange_rs2pao, calculate_moment_of_inertia
```

There is no current calculation-driver dispatch to `calculate_jijk`; the
`njijk` condition prevents the ordinary pair consumers from running in that
case.  There is also no driver dispatch to `calculate_exchange_rs2pao`; the
name `exchange_p2rs` selects the Gauss-Legendre method instead.

## Existing tests and what they actually establish

| Test/campaign | Formulation reached | Current assertion | What it does not establish |
|---|---|---|---|
| `Example_exchange_bccFe` | `calculate_exchange` + `calculate_exchange_twoindex` through recursion/block, no HOH | `jij.out` J columns for two pairs and `dij.out` D columns; D is zero for this bcc-Fe case | No A-file/tensor assertion; no two-index SO/FO/parts assertion; not an independent formula oracle |
| `Example_exchange_bccFe_hoh` | Same consumers with `hoh=.true.` | J and D reference columns for two pairs | Same A/two-index omissions; HOH/native agreement is not asserted |
| `Example_paoflow_exchange_p2rs_bccFe` | `calculate_exchange_gauss_legendre` after PAOFLOW import | J and D reference columns for two pairs | No A check and no Simpson-versus-Gauss equivalence; does not call `calculate_exchange_rs2pao` |
| `Triad_triad_bccFe_jij` | `calculate_exchange` for recursion, Lehmann, Dyson | two J shells; Lehmann==Dyson tight and recursion/Lehmann ratio envelope | No D/A/two-index validation; route comparison is not a fixed-broadening equality |
| `Triad_triad_bccFe_damping` | `calculate_gilbert_damping` on the same route dispatch | on-site alpha, with SOC and HOH | No inertia or J/D/A tensor assertion |
| `Val05GreenConvergence` | upstream RS/Lehmann/Dyson GF production and direct block diagnostics | GF k/eta/window convergence, DOS, and Sigma=0 Dyson==Lehmann | It does not run exchange consumers |
| `UnitLehmannChain`, `UnitDysonEquivalence`, `UnitGammaSupercell` | upstream reciprocal GF/kernel contracts | chain, backend equivalence, and intersite normalization identities | No end-to-end J/D/A output |

The focused baseline run on the current `build-rf-debug` tree passed all three
configured exchange fixtures: ordinary exchange (26.16 s), HOH exchange
(47.16 s), and PAOFLOW exchange import (20.45 s).  The CTest label `exchange`
currently selects only `Example_exchange_bccFe`; the HOH and PAOFLOW cases are
configured post-processing tests but are not under that label.  No source
files were changed for this audit.

## Oracle classification

| Relationship | Classification | Legitimate scope |
|---|---|---|
| `calculate_jij` ↔ J part of `calculate_exchange` | Exact algebraic equivalence | Same canonical GF arrays, Δ matrices, mesh, Simpson settings, and scaling.  This is a duplicate-kernel check, not an independent physical oracle. |
| `calculate_exchange` J/D/A ↔ printed `J_tensor` | Exact assembly identity | `T=JI+D^skew+A` as constructed.  Extracting D from `antisym(T)` additionally requires `A` symmetric; that is not enforced by the code. |
| Full tensor scalar trace ↔ scalar J | Generally inequivalent | `tr(T)/3 = J + tr(A)/3`; it equals J only when the anisotropic term has zero trace.  `tr(A)/3` is not a J oracle. |
| DMI vector ↔ antisymmetric tensor | Conditional algebraic equivalence | Exact for the explicitly constructed `D^skew`; not automatically recoverable from the total tensor if `A` has an antisymmetric part. |
| Symmetric tensor ↔ anisotropic exchange | Conditional physical projection | Use `sym(A)=(A+A^T)/2` as the symmetric anisotropic part only after documenting whether the model assigns antisymmetric pieces to D.  Do not silently equate the raw full `aij` with a symmetric tensor. |
| Two-index SO/FO ↔ combined full exchange | Intentionally inequivalent | Compare only the source component identities or a separately defined perturbative reference.  SO and FO are order-resolved quantities, not alternative full J/D/A values. |
| Gauss-Legendre eta method ↔ Simpson real-energy method | Conditional numerical-integration oracle | Requires matched GF/Hamiltonian data, contour/Jacobian derivation, and trace-convention reconciliation.  Current source uses different real/imaginary trace selectors between paths, so equality must not be presumed. |
| RS GF ↔ Lehmann GF | Convergence/envelope relationship | At matched contour and converged representation, compare physical observables with k/eta convergence.  Existing J triad permits a broad ratio band, not equality at fixed finite eta. |
| Lehmann GF ↔ Dyson GF | Exact solver-route equivalence at Σ=0 | Full unreduced BZ, same H(k), standard problem, and Σ=0.  This is the strongest current cross-route oracle and is already tested upstream and in the J/damping triads. |
| Native RS GF/H ↔ PAOFLOW import route | Conditional round-trip equivalence | Only after native export/import preserves H, Δ, geometry, Fermi level, contour/mesh, and units.  Existing PAOFLOW fixture is a regression reference, not a matched equivalence test. |
| Auxiliary GF ↔ canonical LKAG | Conditional representation equivalence | Collinear, matched orthogonal representation and P/Δ convention; not arbitrary SOC/noncollinear data. |
| `Jijk` ↔ pair J | No direct equality | Use a finite-difference displacement derivative of independently rerun pair J as the oracle. |
| Damping/inertia ↔ exchange | No direct equivalence | They share GF producers but use SOC torque operators and different energy/derivative constructions. |

## Ranked validation priorities

1. **P0 — output coverage and source-contract checks.** Add read-only checks
   for `aij.out`, the two-index files, and the actual output columns/units;
   explicitly retain the fact that `jtens.out` is currently not populated.
   This closes the present gap where a passing exchange fixture validates only
   J and D.
2. **P0 — preserve the reciprocal GF route oracle.** Keep the tight
   Lehmann==Dyson Σ=0 check and the existing RS-vs-Lehmann envelope.  Do not
   tighten the latter into equality without a joint (N_k,η), and energy-grid
   convergence study.
3. **P1 — exact canonical-kernel cross-check.** Exercise the J-only wrapper
   against the J column of `calculate_exchange` on the same in-memory arrays,
   then check the explicit D-skew/tensor assembly identity.  This is compact
   and protects the shared kernel contract.
4. **P1 — Simpson versus Gauss-Legendre.** Construct a matched native-H route
   and verify the contour change of variables, signs, (w/x^2) Jacobian, and
   `Re/Im` trace selectors before accepting the PAOFLOW path as an independent
   numerical-integration oracle.
5. **P1 — tensor physics.** Add SOC/inversion-controlled cases that inspect
   all D components and both symmetric/antisymmetric A projections.  Assert
   only symmetry-forced limits (for example, a DMI-zero limit when the full
   Hamiltonian and geometry restore the required symmetry), not generic D=0 or
   A symmetry.
6. **P2 — auxiliary-GF representation check.** Expose a controlled collinear
   `hoh` case and compare the correctly selected auxiliary scalar/tensor
   projection against the canonical result after matching representation and
   normalization.
7. **P2 — spin-lattice coupling.** If `calculate_jijk` is made reachable,
   validate its nine components against symmetric finite differences in the
   specified displacement direction, including the `wav` unit conversion.
8. **P3 — dormant import and inertia routes.** First repair/verify their
   call/output contracts, then validate native/PAOFLOW round trips and the
   inertia derivative campaign.  Current test absence is a scope limitation,
   not evidence of correctness.

## VAL-06 checklist

- [x] All public exchange formulations inventoried, including stubs and dormant methods.
- [x] Call sites mapped.
- [x] Physical outputs mapped.
- [x] Spin, SOC, GF-representation, contour, and integration assumptions documented.
- [x] Existing tests and their actual assertions mapped.
- [x] Valid cross-checks and independent-oracle candidates identified.
- [x] Intentionally inequivalent formulations distinguished.
- [x] Validation priorities ranked.
- [x] No production exchange formulas changed.
