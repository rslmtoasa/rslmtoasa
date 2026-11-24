# Hamiltonian export acceptance report

The feature branch is based on `4916e1618f92fd764adf9e5f4e97d319b8d5125b`, including compact conductivity PR #22 and k-space symmetry PR #23. The historical reference was `2fccbffcde6f`. PR #20 restored the older first-order writer and lost its active two-hop assembly. This implementation recovers the directed `h_ik o_k h_kj` sum and keeps the native operator intact during export.

## Implemented operator

| Route | Exported operator |
|---|---|
| First order | The built `ee` or `hall`, including the native fixed-U insertion, plus onsite `lsham` and additive `eecc` or `hallcc` when active. |
| HOH | The auxiliary native hopping minus every two-hop product with the intermediate site's `obarm`, plus built onsite `enim` and `lsham` and additive CCOR. Fixed onsite U already in `enim` enters once. |
| Ordinary periodic | Sparse destination and integer image keys. The second hop adds its image to the first. Repeated paths to the same key are summed. |
| Finite cluster | Actual site IDs and wrapped images. Local and bulk blocks are selected independently on both legs. A removed edge contributes no path. |
| Reciprocal | Full native H(k) on the requested mesh, plus the native S(k) for generalized overlap. The same assembler supplies ordinary, orthogonal and GBT modes. |

Import extends the graph to all effective hopping shells, including separate links to distinct images of a wrapped endpoint. It uses the original spherical spin-orbital basis and ordinary non-HOH actions. It constructs no additional SOC, Hubbard or CCOR operator and fixes E_F from the metadata. The only transport plumbing supplies the imported image displacements to the unchanged velocity formula.

## Material results

The table reports the largest relative H and V action errors over HOH-off/on and all three velocity directions with CCOR off. Native matrix/action comparisons used complex double precision. Matrix comparisons use 64 k points per ordinary periodic case. Matrix and action thresholds were 1e-10. The band eigenvalue threshold was 1e-8 eV.

| Frozen system | Relative H action | Relative V action | Maximum H(k) error, Ry | Hermiticity certification |
|---|---:|---:|---:|---|
| fccCu | 5.95e-16 | 3.38e-16 | 4.45e-16 | **FAIL**. Native H(k) residual 1.21e-8 Ry. |
| fccPt | 8.25e-16 | 3.66e-16 | 1.12e-15 | **FAIL**. Native H(k) residual 1.57e-8 Ry. |
| bccFe | 4.85e-16 | 1.56e-16 | 4.45e-16 | **PASS**. |
| B2 FeCo | 5.13e-16 | 3.13e-16 | 2.50e-16 | **PASS**. |
| FeO AFM-II, U=4 eV | 1.48e-15 | 1.36e-15 | 4.45e-15 | **PASS**. |
| B2 FeCo impurity | 7.09e-16 | 5.44e-16 | Cluster | Periodic criterion does not apply. |

The Cu and Pt periodic exports reproduced the native matrices and their existing reciprocity defects. The calculation used the 6³ representative replica for reciprocal comparisons and the 4³ closed replica for real-space transport. The 4³ reciprocal case exceeded the native reciprocal neighbor table dimensions. A bounds-checked build located the overflow before an exporter guard converted the failure to a diagnostic. Ritz values from ordinary matrices matched to 1.88e-13 eV or better, using the same Hermitian triangle convention even for the two matrix sets with failed Hermiticity certification.

The original FeO archive SHA-256 matched `d6da31bb61056a67583258364dc340411daf79f406912d157dc9e8a2c76a7ea7`. Both HOH settings loaded the archived `*_out.nml` potentials and referenced sidecars with local radial provenance. Native assembly gave 28 valence electrons and a maximum nonzero U block of 0.1469971714457357 Ry for both builders. Separate frozen first- and second-order builds used the same state with SOC and CCOR off. The comparison used actual current native H(k) and Hψ. It does not certify convergence of the archived SCF state. The installer rebuilds native screened constants and preserves the archive unchanged.

The native and imported real-space DOS files matched byte for byte for mixed-boundary Cu, bccFe and the FeCo impurity, including the local orbital DOS files. The first bccFe DOS window produced nonconvergent moments in both runs. A shared [-3.5, 3.5] Ry window resolved that test. Ordinary reciprocal Gaussian DOS comparisons use the eigenvalues of the exported Fourier matrices and the native matrices on the same mesh.

## Special reciprocal modes

The matrix reference was reconstructed from the pinned baseline's eigensystem. In the generalized case S was also recovered from generalized-normalized eigenvectors and checked independently.

| Frozen test | H maximum, Ry | S maximum | Eigenvalue maximum, Ry |
|---|---:|---:|---:|
| bccFe exact orthogonal, strux_lib | 4.45e-15 | — | 3.56e-15 |
| bccFe generalized overlap, HOH on | 7.40e-13 | 2.67e-15 | 7.56e-13 |
| B2 FeCo GBT single-q, strux_lib | 5.64e-13 | — | 1.07e-13 |

The legacy screening attempt at exact orthogonal assembly was rejected by the native strux_lib prerequisite. Reciprocal-only output has no effective real-space import route.

## Conductivity and qualifications

All six fccPt charge, spin and orbital comparisons used fully closed PBC and the same deterministic per-type trace, frozen state, energy grid, Fermi level and FP64 legacy moments with factorized Gamma. Both HOH settings passed the total and projected tests in both real and imaginary columns. The largest relative difference in a canonical response file was 5.43e-15 at order 28. The legacy seven-column route also passed, including the Cartesian orbital conversion. Native response files from the candidate with export enabled or disabled were byte-identical to the pinned baseline. The native total imaginary charge signal in the small 4³ cell was finite and was reproduced by the import. Accordingly the expected near-zero imaginary charge check **FAILS** on this fixture.

On the legal bccFe CCOR case, with screened derivatives and the fixture that contains `vrmax`, H action agreed to 1.50e-16 relative. The maximum velocity differences along x, y and z were 0.20846, 0.19214 and 0.19831 in native Ry Å displacement units, about 10–10.5% relative. Ordinary tight-binding velocity includes the offsite CCOR commutator through the full H. Native real-space velocity omits it. The synthetic nonuniform-weight case checks that weighting the two original HOH legs separately is generally different from weighting their completed hopping. Full CCOR transport equality was outside the acceptance condition.

Order-50 fccPt HOH export-disabled CPU process measurements over three runs were median 5.09 s (5.06–5.48 s) and 854,632 KiB on the baseline, versus 5.42 s (5.27–5.45 s) and 855,252 KiB on the candidate. These process timings include startup, geometry and shared-host variability. The paired native response files matched exactly. The fccPt 50-order HOH export had 19 original bond records, 361 two-hop paths and 27 nonzero effective blocks. The canonical scalar file was about 335 KiB. The raw run manifest records process time, peak RSS, input and executable hashes, and output file sizes. It does not isolate exporter peak RSS from the process high-water mark.

## Acceptance gates

| Gate | Status | Evidence |
|---|---|---|
| 1. Shared assembly and periodic algebra | **PASS**, with separate Cu/Pt Hermiticity failures above | Dense mixed local/bulk oracle 5.56e-17 Ry maximum, native action 1.58e-16 Ry maximum, skew integer mapping and repeated-image test. All five periodic Fourier matrix comparisons passed the operator reproduction tolerances. |
| 2. Finite graph and writers | **PASS** on tested paths | Cu mixed boundary and FeCo impurity with host/local paths. Existing export smoke tests pass for both legacy names. Standalone real-material surface/interface and explicit-texture replay **NOT RUN**. |
| 3. Effective import and conductivity | **PASS** for round-trip equality | Six material H and V actions, canonical and legacy conductivity, malformed metadata/nonfinite/duplicate/truncated rejections and DOS comparisons. Native zero imaginary charge certification **FAIL** on the small Pt test. |
| 4. Native reciprocal export | **PASS** on tested legal modes | Pinned-baseline eigensystem checks for GBT, exact orthogonal and generalized overlap. |
| 5. Regression and documentation | **PASS** on focused CPU checks | Original export, external PAOFLOW and keyword documentation checks. The full CPU quick suite passed 71/85 candidate tests versus 70/84 baseline tests. The same 14 reciprocal SCF and DOS checks failed in both builds. |

CPU builds used GNU Fortran 14.2.0, CMake 3.31.6, system BLAS/LAPACK, OpenMP enabled, MPI and spglib disabled, and one rank with one thread for numerical comparisons. MPI configuration was attempted but the installed Open MPI wrapper advertised a missing Fortran module directory. MPI execution, Intel, CUDA, a 250-order timing run, full surface/interface transport, explicit-texture import and independent reciprocal band paths beyond the sampled mesh were **NOT RUN**. No additional material study, velocity export, k-space importer or SCF-from-import route was added.

Implementation review entry points

- `source/hamiltonian_export.f90` for native sparse block assembly, global spin frame and serialization.
- `source/hamiltonian_effective_io.f90` for validation and graph replacement.
- `source/recursion_transport.f90` for the imported displacement override in `setup_kubo_bonds`.
- `source/reciprocal_export.f90` for bounded native H(k)/S(k) export.
- `docs/HAMILTONIAN_EXPORT.md` and `tests/hamiltonian_export/README.md` for input examples and reproducible checks.

## Review follow-up

The mixed-boundary Cu CLI regression confirms that automatic scope selects
the finite cluster, native and effective conductivity agree, invalid orbital
order is rejected, and a payload changed without changing its record count
fails the v3 checksum. The effective reader also accepts older files with separate
response and mesh Fermi references. The regression supplies distinct values to
exercise that compatibility path.

The same frozen Cu Hamiltonian payload is byte-identical with
`recur='chebyshev'`, `'block'` and `'lanczos'`. The new native block DOS
case exercises the production block Green-function path and checks a finite
nonzero site DOS. Effective DOS continues to use Chebyshev. In a separate
fully closed 6³ Pt check, the v3 canonical and metadata-matched seven-column
imports agreed within 4.3e-14 S/cm on the total and orbital-resolved responses.
The standalone external `conductivity_p2rs` regression also passed.

## Cluster seven-column output and spglib check

Both real-space scopes now write the seven-column Cartesian view alongside
the canonical payload and matching metadata. The cluster view uses actual
site-orbital indices and records translations of the replica along closed
directions only. The mixed-boundary Cu test reimports `rs2paoham.dat` directly
and reproduces native total and orbital-resolved conductivity and DOS. The
fully open Cu replica has zero translation columns. The open HOH B2 FeCo
impurity round trip through this file agrees with native actions at
7.62e-16 relative error for H and at most 6.17e-16 for the three velocities.

The feature and baseline build both failed to link against upstream
spglib 1.6.3 (commit `cfd19b2`), which lacks
`spgms_get_symmetry_with_site_tensors`. No runtime tests ran with that
library. With spglib 2.6.0, GNU Fortran 14.2.0, system BLAS/LAPACK,
OpenMP enabled and MPI disabled, the CPU quick suite passed 88/88 on
the feature and 87/87 on main `4916e16`. All 14 reciprocal tests
that failed in the spglib-disabled run passed with spglib enabled.

## Final option audit

| Finding | Assessment and follow-up |
| --- | --- |
| `kspace_ham_order='second'` with HOH off can label a first-order matrix as second order. | Confirmed on frozen bcc Fe. The first- and second-order requests produced identical 2,592 H records on a 2×2×2 mesh, while the header reported different requested orders. The writer should reject this combination or report the resolved order separately. The manual requires a built HOH operator for this example. |
| Fully open bulk can resolve to compact periodic output under `auto`. | Confirmed on bcc Fe with `source_boundary F F F`. This retains the existing compact-bulk behavior. Use explicit `cluster` for finite-graph fidelity. Reporting the resolved scope at the console would make this easier to detect. |
| Conductivity previously had an optional writer call. | Native conductivity now runs the original driver from `main`. The CLI test supplies an enabled export flag to native conductivity and checks that it writes no Hamiltonian files. A separate frozen export supplies the replay data. |
| `export='python'` and `'rs2pao'` share assembly and Cartesian serialization. | The two values preserve historical basenames. Both also write the spherical canonical view, which makes the basis and graph directly auditable. Either enabled value triggers the same reciprocal writer. |
| `conductivity_p2rs` and `conductivity_effective` share the response driver. | Keep both entry points. The former reads standalone external bulk PAOFLOW files, and the latter validates a matched RS-LMTO-ASA export and rebuilds all image-resolved shells. The two DOS entry points share the recursion and Green-function driver. |

The four examples in the [manual](docs/HAMILTONIAN_EXPORT.md) were exercised
with spglib 2.6.0. For the mixed-boundary Cu example, the largest absolute
difference across the native and seven-column replay response tables was
7.13e-14 S/cm and the DOS tables matched numerically. The focused CTest
selection passed 7/7, including keyword documentation coverage, both legacy
writer names, external PAOFLOW conductivity, the effective algebra unit and
both effective CLI tests. These figures are for the pre-separation review.

## Native conductivity and export separation

The native conductivity driver now matches `main` byte for byte in source.
The original DOS reconstruction, Fermi determination, Kubo moments and
response reconstruction stay in the native response workflow. A separate
`hamiltonian_export` job writes the frozen operator and its input Fermi
reference. The importer retains support for the optional mesh reference in
earlier files.

On the mixed-boundary Cu fixture with a common fixed Fermi reference, native
conductivity wrote no Hamiltonian files even with `export='rs2pao'` selected.
The native response and Cu DOS tables were byte-identical to those from
`main` at `4916e16` with the same compiler and spglib 2.6.0. The independent
seven-column export and replay agreed with native conductivity within
7.13e-14 S/cm, and the native and replayed Cu DOS tables matched numerically.
