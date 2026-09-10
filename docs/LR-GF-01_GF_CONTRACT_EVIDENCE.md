# LR-GF-01: one-electron Green-function contract evidence

## Disposition

LR-GF-01 certifies the live post-purge one-electron infrastructure. It does not
restore or depend on the removed TD-DFT response implementation.

The trusted clean-room primitive is the orthogonal, reciprocal,
coefficient-space Lehmann Green function

\[
G_c(k,z)=C_k(zI-\epsilon_k)^{-1}C_k^\dagger
       =[zI-H(k)]^{-1},
\]

for `reciprocal_mode='ham_only'`, a complete finite eigensystem, and
`z=E+i eta` for the retarded route. The raw matrix is an orthogonal LMTO
coefficient-space resolvent. It is not, by itself, the augmented physical
spatial/radial Green function.

## Audit identity

- Branch: `fable_v4`
- Baseline commit audited: `2ba7087473d72c4fd4ab60daa1ad36c1fc5b31df`
  (`td-dft: establish clean-room linear-response baseline`)
- Audit date: 2026-09-10
- Certification changes: this document, the focused tests below, the small
  complex Chebyshev reconstruction kernel, and terminology corrections.

The post-purge tree contains no restored TD-DFT implementation. LR-GF-01 does
not implement `chi0`, Dyson TD-DFT, Ward/Goldstone repair, mode extraction,
radial augmentation, or response gauge handling.

## Live call graph and data flow

### Reciprocal mesh eigenpairs

The normal mesh path is:

```text
calculation / reciprocal consumers
  -> reciprocal%build_kspace_hamiltonian
     -> reciprocal_fourier::execute_normal_mesh_tiles
        -> reciprocal_assembler::assemble_batch
        -> reciprocal_execution_backend::execute_batch
           -> standard Hermitian eigensolve (ZHEEV) for ham_only
  -> reciprocal%diagonalize_hamiltonian
     -> cache/generation check; fused builds already own the eigensystem
```

`build_kspace_hamiltonian` assembles the same first- or second-order operator
used by the public reciprocal path. The second-order operator includes the
active onsite, hopping, HOH, optional CCOR, and spin-orbit terms through the
reciprocal assembler. `reciprocal_bands::diagonalize_hamiltonian` checks
Hermiticity and the operator-generation fingerprint and exposes the complete
`eigenvalues(nband,nk)` and `eigenvectors(nbasis,nband,nk)` arrays.

For `ham_only`, `ZHEEV` returns

\[
H(k)C_k=C_k\epsilon_k,\qquad C_k^\dagger C_k=I.
\]

The live implementation does not truncate the band sum: the Lehmann kernel
loops over all `n=1:nmat`.

### Arbitrary-​`k` eigenpairs

```text
reciprocal%calculate_eigenpairs_at_kpoints
  -> fold_kpoint (canonical fractional point in [-1/2,1/2))
  -> exact folded-point deduplication
  -> reciprocal_assembler::assemble_batch or host H(k) handoff
  -> execution_backend::execute_batch
  -> caller-owned eigenvalues/eigenvectors
```

`build_hamiltonian_at_kpoint` uses the same single-point assembler and the same
first/second-order selection as the normal mesh. The service does not write
`k_points`, `hk_bulk`, bands, DOS, or normal-mesh eigensystem caches. It returns
the complete basis eigensystem, not only states near the Fermi level.

The arbitrary-​`k` test also checks a coincident mesh point, an off-mesh point,
folding by an integer reciprocal vector, exact duplicate folded-point reuse,
first-order and second-order/SOC assembly, tiling, and operator-generation
refresh. In the simple one-site fixture the folded eigenvectors agree exactly;
this is a representative gauge check, not the full multi-site basis-gauge
claim. The full site-dependent gauge

\[
c_a(k+G)=e^{\pm iG\cdot\tau_a}c_a(k)
\]

is explicitly deferred to LR-KQ-00.

### Reciprocal coefficient-space Green function

```text
calculation%post_processing_exchange or post_processing_kspace_green
  -> reciprocal%fill_green
     -> ham_only guard
     -> build_green_contour: z(E)=E+i*green_eta
     -> fill_green_lehmann
        -> existing eigenpairs, if available
        -> CPU lehmann_pair_block
           or backend contract / CUDA Lehmann contraction
        -> green%gij, green%gji and Pauli spin blocks
     -> fill_green_dyson (selected backend only)
        -> dyson_kspace_inverse per (k,z), Sigma provider
```

The strict Lehmann block is

\[
G_{ij}(z)=\frac1{N_k}\sum_{k,n}e^{i k\cdot(R_i-R_j)}
 \frac{c_{i,nk}c_{j,nk}^\dagger}{z-\epsilon_{nk}}.
\]

The direct Dyson kernel is an independent dense inverse for the orthogonal
case:

\[
G_D(k,z)=[zI-H(k)-\Sigma(z)]^{-1}.
\]

With `Sigma=0`, it is the oracle used for the Lehmann equivalence tests. The
reciprocal contour uses energies in Ry and `green_eta > 0` for retarded values;
the advanced value is used only by the mathematical validation test.

### Native RS block-recursion Green function

```text
calculation_reciprocal::run_intersite_moments
  -> recursion%recur_b_ij
     -> four phase seeds for each pair
  -> green%calculate_intersite_gf
     -> green_block::calculate_intersite_gf_core
        -> block_green_ij
           -> recursion%get_terminf
           -> green_lanczos::bgreen / bgreen_complex
```

The accepted block terminator is used as existing infrastructure. LR-GF-01
does not redesign it, change its branch handling, or treat its convergence as a
new acceptance criterion.

The recursion matvec applies the effective orthogonal real-space Hamiltonian.
For HOH/second-order mode the source-level action is

\[
H_{\rm eff}=E_\nu+h-h\bar O h+H_{\rm SO}+H_{\rm enabled\ corrections},
\]

implemented in the order `psi2 - hohpsi + enupsi + socpsi`, with optional CCOR
terms in the enabled paths. The recursion block metric is ordinary Euclidean
(the four seeds are constructed from ordinary diagonal identity blocks).

For a pair `(i,j)`, the four-seed reconstruction in
`green_block.f90::calculate_intersite_gf_core` is the projector orientation

\[
\texttt{gij}=P_i G P_j^\dagger,\qquad
\texttt{gji}=P_j G P_i^\dagger.
\]

The source confirms the orientation through the seed combinations and the
separate assignments at `green_block.f90:267-273`. This is a coefficient/block
statement; it does not imply physical radial augmentation.

### Native RS Chebyshev Green function

```text
calculation_reciprocal::run_intersite_moments
  -> recursion%chebyshev_recur_ij
     -> cheb_moments_cpu / legacy or selected fast backend
        -> ham_vec_matmul or ham_hoh_vec_matmul
  -> green%calculate_intersite_gf
     -> green_chebyshev::chebyshev_green_ij
        -> cheb_green_fast for real-energy reconstruction
     -> green_chebyshev::chebyshev_green_ij_eta
        -> cheb_green_complex for complex-energy reconstruction
```

Moments are generated as

\[
\mu_n=\langle\psi_0|T_n[(H_{\rm eff}-b)/a]|\psi_0\rangle,
\]

with the same effective operator contract. The complex reconstruction uses
the native Jackson convention and retarded transfer factor

\[
G(z)\simeq\sum_n g_n c_n
\frac{-i\exp[-in\arccos((z-b)/a)]}
{\sqrt{a^2-(z-b)^2}}.
\]

This is a controlled truncated approximation, not an exact finite-order
resolvent. LR-GF-01 adds `UnitChebyshevGFOracle`: moments are generated by an
independent dense matrix recurrence and the production reconstruction is
compared with a separate dense matrix inverse at increasing orders.

## Supported modes and formal generalized-overlap statement

| Mode or route | LR-GF-01 status | Contract |
|---|---|---|
| Reciprocal `ham_only` Lehmann | **proven** | Complete orthogonal eigensystem and `G_c(k,z)=[zI-H(k)]^{-1}`. |
| Reciprocal direct Dyson, `ham_only` | **proven** | Independent `zI-H-Sigma` inverse; `Sigma=0` agrees with Lehmann to solver tolerance. |
| Arbitrary-​`k` service | **proven** | Same H(k) contract, complete eigensystem, folding, and no normal-mesh state mutation. |
| RS block recursion | **supported with approximation** | Effective orthogonal H and accepted finite-chain/terminator Green function; terminator redesign is out of scope. |
| RS Chebyshev | **supported with approximation** | Same effective H, Jackson-truncated polynomial reconstruction; finite-order convergence is tested. |
| `generalized_overlap_proxy` in band/arbitrary-​`k` eigensolver | **supported with approximation** | `ZHEGV` and O-normalized eigenvectors are available for the eigenproblem itself. |
| Generalized-overlap reciprocal response GF | **unsupported** | `reciprocal%fill_green` fails closed unless `reciprocal_mode='ham_only'`. The guard is preserved. |
| `generalized_overlap_kanpur` | **unsupported** | Not implemented as a production generalized reciprocal response route. |

For an O-normalized generalized eigensystem

\[
HC=OC\epsilon,\qquad C^\dagger OC=I,
\]

the formal spectral identity is

\[
(zO-H)^{-1}=C(z-\epsilon)^{-1}C^\dagger.
\]

The unresolved LR issue is not this spectral representation. It is metric and
operator placement in later response contractions. Those semantics belong to
LR-METRIC-00, so LR-GF-01 does not broaden the production response guard.

## Basis and energy conventions established by the live code

- `H(k)` and eigenvectors are site-major, with each site occupying one `nb`
  coefficient block. Arbitrary-​`k` returns `(nb*nrec, nband, nk)` eigenvectors.
- Pair offsets are zero-based in the pure Lehmann kernel and are converted from
  `lattice%iz` by `pair_geometry`; site rows are
  `((site-1)*nb+1):(site*nb)`.
- Normal-mode spin blocks use the `up` block followed by the `down` block;
  `spin_off=norb`. Pauli decomposition is the production charge/x/y/z
  decomposition and is independently pinned by `UnitLehmannChain` and
  `UnitKspaceGFValidation`.
- Reciprocal k points are fractional/direct coordinates. Fourier phases use
  `2*pi*k_frac dot R_frac`; reciprocal GF pair phases use
  `exp(+i*k dot (R_i-R_j))`.
- Reciprocal retarded energies are `z=E+i*green_eta`, with E and eta in Ry;
  `green_eta` is not the Chebyshev scaled Fermi variable.
- Chebyshev uses the affine map `(H_eff-b)/a`, where the live code obtains
  `a=(emax-emin)/(2-0.3)` and `b=(emax+emin)/2`, and applies Jackson weights.

## Coefficient-space versus physical-space Green functions

The certified reciprocal object is `G_c(k,z)` in the orthogonal LMTO
coefficient basis. It is not a claim that the raw matrix is the full physical

\[
G(r,r';z)
\]

with radial functions, augmentation, or all physical-space normalization
factors. Physical augmentation is deferred to LR-BASIS-00.

Likewise, raw reciprocal coefficient blocks are not declared equal to a
screened or auxiliary LMTO representation. The live code has explicit
`green%auxiliary_gij` and `green%transform_auxiliary_gij` operations; those
representation transformations are deferred to LR-REP-00.

## Existing tests and independent evidence

| Test | What it proves | Independence / limitation |
|---|---|---|
| `tests/unit/test_kspace_gf_validation.f90` (`UnitKspaceGFValidation`) | Dense inverse equivalence for five Hermitian matrices and eight complex energies; direct residual; retarded/advanced identity; causal sign; high-energy normalization; degenerate-subspace invariance; spin algebra and DOS checks. Observed inverse error `1.05e-14`, residual `1.14e-15`, advanced error `1.39e-16`, and `zG-I` error `8.35e-11`. | Lehmann uses `lehmann_kspace_resolvent`; the inverse oracle uses independent LAPACK `zgetrf/zgetri` through `dyson_kspace_inverse`. The spectral and asymptotic checks are direct matrix properties. |
| `tests/unit/test_lehmann_chain.f90` (`UnitLehmannChain`) | One-band chain closed forms for onsite and `m=2` intersite phase, plus normalization and Pauli signs. | Analytic chain oracle; no LMTO object or response code. |
| `tests/unit/test_dyson_equivalence.f90` (`UnitDysonEquivalence`) | Direct reciprocal inverse/Dyson blocks and the production Lehmann pair-block accumulation agree, including intersite orientation and Sigma sign behavior. | Strong route equivalence, but both routes use the same assembled fixture H(k); it is not a second H(k) assembly implementation. |
| `tests/unit/test_arbitrary_k_eigenpairs.f90` (`UnitArbitraryKEigenpairs`) | Coincident mesh/off-mesh eigenpairs, folding, q+G spectra/eigenvectors in the fixture, complete eigenvectors, tiling, cache generation, normal-mesh compatibility, second-order/SOC, and generalized eigenproblem diagnostics. | The service and normal mesh share the reciprocal assembler by design; direct `ZHEEV` checks validate residuals/gauges, while the test is primarily a service/state contract. |
| `tests/unit/test_chebyshev_gf_oracle.f90` (`UnitChebyshevGFOracle`) | Complex-energy Chebyshev GF against independent dense inverse for orders 16, 32, 64, 128. Observed max errors `2.11`, `1.01`, `0.354`, `0.102`: convergent truncated approximation. | Dense moments use a direct matrix recurrence and the oracle uses an independent dense inverse. The tested reconstruction is the production Jackson/complex transfer kernel. |
| `tests/unit/test_acc04_arbitrary_k_source.py` and related reciprocal source-contract tests | Preserve backend seam, host assembly, eigenvector ownership, and generalized-overlap guards. | Static source contracts; they do not replace numerical tests. |
| `tests/validation/val05_green_convergence.py` and the k-space report path | Existing functional route comparison/convergence evidence for real material workflows. | Material, mesh, eta, and recursion approximation claims remain scoped to their recorded cases. |

The direct Dyson inverse is the independent oracle for the reciprocal Lehmann
identity. A wrapper around the same eigensum would not count as an oracle; the
tests above distinguish these cases.

## Required dispositions and deferrals

| Area | Disposition |
|---|---|
| Orthogonal reciprocal Lehmann resolvent | **proven** by dense inverse, residual, causal sign, asymptotic, and gauge tests. |
| Reciprocal arbitrary-​`k` service | **proven** for the live service contract and tested first/second-order fixtures. |
| RS block recursion | **supported with approximation**; the effective operator, Euclidean seeds, and `gij/gji` orientation are certified. The native terminator is accepted and unchanged. |
| RS Chebyshev | **supported with approximation**; independent finite-system convergence oracle added. |
| Generalized-overlap spectral identity | **proven as a formal mathematical identity**, but not enabled for reciprocal response contractions. |
| Physical radial/spatial augmentation | **deferred** to LR-BASIS-00. |
| Generalized-overlap response metric/operator placement | **deferred** to LR-METRIC-00. |
| Full multi-site reciprocal basis gauge under `k -> k+G` | **deferred** to LR-KQ-00. |
| Physical/screened/auxiliary representation equivalence | **deferred** to LR-REP-00. |
| Generalized-overlap reciprocal response GF | **unsupported** in LR-GF-01; production guard remains active. |
| Complete finite-q response or TD-DFT | **out of scope / unsupported** here. |

## Acceptance checklist

- [x] Current `fable_v4` post-purge tree traced.
- [x] Reciprocal `ham_only` eigensystem contract documented.
- [x] Complete-spectrum Lehmann construction verified.
- [x] Lehmann versus dense inverse tested.
- [x] Resolvent residual tested.
- [x] Retarded/advanced relation tested.
- [x] Spectral-sign convention tested.
- [x] High-energy normalization tested.
- [x] Degenerate-subspace invariance established.
- [x] Arbitrary-​`k` contract checked.
- [x] Full `k+G` site gauge explicitly deferred to LR-KQ-00.
- [x] Generalized-overlap formalism documented; production guard preserved.
- [x] RS block GF semantics documented without changing the accepted terminator.
- [x] RS `gij/gji` orientation confirmed.
- [x] RS Chebyshev operator contract confirmed.
- [x] Independent Chebyshev GF oracle added.
- [x] Coefficient-space versus physical-spatial GF terminology cleaned up.
- [x] `docs/LR-GF-01_GF_CONTRACT_EVIDENCE.md` created.
- [x] Relevant unit tests pass.
- [x] No purged TD-DFT implementation restored.

## Verification command and result

With `RUN_UNIT_TESTS=ON`, the focused suite passed:

```text
ctest --test-dir build --output-on-failure -R \
  'Unit(KspaceGFValidation|LehmannChain|DysonEquivalence|ArbitraryKEigenpairs|ChebyshevGFOracle)'
100% tests passed, 0 tests failed out of 5
```

The accepted conclusion is therefore:

\[
\boxed{G_c(k,z)=[zI-H(k)]^{-1}}
\]

for orthogonal `ham_only`, with the four explicitly deferred representation,
metric, augmentation, and reciprocal-gauge boundaries above.
