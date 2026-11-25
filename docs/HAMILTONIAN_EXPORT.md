# Effective Hamiltonian export and import

Use this workflow to save the Hamiltonian built from a fixed set of RS-LMTO-ASA potentials, inspect it outside the code, or replay real-space DOS, conductivity and pair exchange from its hopping matrix. For a frozen run, copy an existing working deck into a separate directory. Keep the geometry, atom labels and potential database, the intended atom potentials, and any Hubbard density-matrix sidecars. Check the run log to confirm that the intended saved potentials and density matrices were loaded. Edit the existing namelists in place. The examples below list changes to those namelists.

## 1. Choose the workflow

`&calculation` chooses what the executable does. For a standalone export or replay, set `pre_processing='none'` and `processing='none'`.

| `post_processing` | Operator used | Result |
| --- | --- | --- |
| `hamiltonian_export` | Builds the frozen native operator. | Writes the selected representation. Requires `export='rs2pao'` or `'python'`. Runs no SCF or recursion iterations. |
| `conductivity` | Builds the frozen native operator. `gf_route='recursion'` selects native real-space transport. | Computes the native response. Use a separate `hamiltonian_export` job for matrix files. |
| `dos_frozen` | Builds the frozen native real-space operator. | Computes site DOS with `recur='chebyshev'` or `'block'`. Use a separate frozen export for a hopping file. |
| `conductivity_import` | Loads a real-space export and its metadata. | Runs `linear_in` and `linear_out` through native Kubo reconstruction. Requires `recur='chebyshev'`. |
| `dos_import` | Loads the same type of real-space export. | Computes site DOS with `recur='chebyshev'` or `'block'`. |
| `exchange_import` | Loads the same export and uses the frozen atomic potentials for the LMTO exchange vertices. | Computes native pair J, DMI and exchange tensors with `njij` and `ijpair`. Uses Chebyshev or block recursion and the native spd exchange kernel. |
| `density_of_states` | Builds the native reciprocal operator. | Computes k-space DOS. |
| `conductivity_p2rs` | Loads a standalone external PAOFLOW `paoham.dat` through the existing bulk reader. | Computes real-space conductivity. Shares the response driver with `conductivity_import`. |

An SCF run can also export its last built operator. Set `export` and `export_scope` in the normal SCF deck and leave `post_processing='none'` if that is the only output you want. To export from the potentials saved at the end of SCF, use a second, frozen job. In a combined SCF and conductivity run, any Hamiltonian export belongs to the SCF stage.

`gf_route` controls the native conductivity Green-function route. For a native-to-effective transport comparison, use `gf_route='recursion'` in the native conductivity job and a real-space scope in the separate export job. Use a separate `export_scope='kspace'` job to inspect the reciprocal matrix.

## 2. Choose the export representation and filename

`&hamiltonian export_scope` selects which matrix is built for output. `&hamiltonian export` enables output and selects the seven-column basename for a real-space export.

| `export_scope` | Matrix | When to choose it |
| --- | --- | --- |
| `periodic` | Compact $H_{ab}(R)$ between sites in a repeating bulk basis. $R$ counts translations of the exported bulk cell. | Repeating bulk or sheets, with one cell along each open direction. Requires one recursion class per basis site and a consistent representative neighbor map. |
| `cluster` | $H_{ij}(R)$ for individual sites in the actual finite calculation. $R$ counts translations of the whole replica along its closed directions. | Open boundaries, mixed boundaries, surfaces, interfaces, impurities and explicit textures. Choose this for the exact finite operator even when all three directions are closed. |
| `kspace` | Native $H(k)$ and, for the generalized overlap proxy, $S(k)$ on the full configured reciprocal mesh. | Bulk reciprocal work, including supported exact orthogonal and generalized Bloch theorem (GBT) configurations. |
| `auto` (default) | Resolves to one of the scopes above. | GBT single-$q$ selects `kspace`. Surfaces, interfaces, impurities, explicit textures and replicas built with `pbc=.true.` and an open direction select `cluster`. A cut bulk cluster with `pbc=.false.` can select `periodic` to reconstruct a bulk model from its central representatives. Use explicit `cluster` to reproduce those open boundaries. |

For the cluster graph, effective periodicity in each direction is the logical AND of `pbc` and that direction’s flag (`b1`, `b2`, or `b3`). The presence of a cell vector never closes an open direction. With HOH, intermediate sites on a two-hop path must exist in the native graph. Cropping an already assembled bulk $H(R)$ does not generally reproduce an open cluster.

| `export` | Real-space files | Reciprocal file |
| --- | --- | --- |
| `none` (default) | No files. | No file. |
| `rs2pao` | `rs2paoham.dat`, `rs_tb.tb`, `rs_tb.meta` | `rs_tb.kspace` |
| `python` | `rs_tb_paoham.dat`, `rs_tb.tb`, `rs_tb.meta` | `rs_tb.kspace` |

The two enabled values use the same matrix assembler and seven-column layout. Both are readable from Python. Use `rs2pao` for the standard `rs2paoham.dat` name. For reciprocal export, either value simply enables the single k-space writer. The seven-column layout is the RS/PAOFLOW layout. A Wannier90 `_hr.dat` writer is outside this workflow.

## 3. Choose the input format for real-space replay

Set `post_processing` to `'conductivity_import'`, `'dos_import'`, or `'exchange_import'`, and use `gf_route='recursion'`. The earlier `conductivity_effective` and `dos_effective` names are aliases. Use `export='none'` in the import deck and choose the input format with `&hamiltonian effective_import`.

| `effective_import` | Files to copy |
| --- | --- |
| `auto` (default) | The matching `rs_tb.meta` and either payload. Prefers `rs_tb.tb` when present. |
| `canonical` | `rs_tb.tb` and the matching `rs_tb.meta` |
| `legacy` | A matching seven-column file and `rs_tb.meta`. Lookup order is `paoham.dat`, `rs2paoham.dat`, `rs_tb_paoham.dat`. |
| `external` | A seven-column file for the explicit cluster, using the same filename order. Read geometry and Fermi energy from the deck and generate `rs_tb.meta` if absent. |

Use a clean run directory so a stray `paoham.dat` cannot take precedence. `export_scope` has no role in import. For paired imports, the metadata determines the representation and boundary flags. Match the basis size, `nsp`, geometry, atom labels and site ordering. Keep the atom files and any geometry inputs needed to construct the same calculation lattice. For compact periodic import, build a matching bulk or sheet replica with `pbc=.true.` and the same boundary flags as the export. Each open direction must contain one cell. For example, explicitly select `export_scope='periodic'` for a three-site MoS₂ cell with `b1=b2=.true.`, `b3=.false.`, and `n3=1`. The compact file retains 54 spin-orbitals for an `spd` basis and can be imported at a different in-plane replica size. For cluster import, match the exported cluster including its boundary flags and replica dimensions. Every nonzero effective shell is recovered from the file, including separate images of the same endpoint.

Effective replay uses the ordinary non-HOH tight-binding action and the existing displacement-times-hopping velocity. HOH, onsite SOC, Hubbard and CCOR contributions come from the export. New exports store the onsite SOC decomposition in `rs_tb.meta` for the spin and SOC torque operators. The effective drivers select generic FP64 recursion with `cheb_backend='legacy'` and `transport_backend='legacy'` and disable the GPU plugin. They fix the Fermi energy from metadata. A frozen export records the input `fermi`. For a paired validation run, use the same `fermi` and `fix_fermi=.true.` in the native conductivity job. Match energy limits, grid size, recursion order, stochastic settings, velocity directions, normalization and response operators.

For ordinary unweighted or uniformly weighted currents with CCOR off, the two-hop displacement is the sum of the leg displacements, so native HOH and effective velocities agree. Nonuniform per-type velocity weights act on each leg in the native HOH current and cannot in general be reconstructed from the final bond. With offsite CCOR, the imported current includes the position commutator of those hoppings, which is absent from the native real-space velocity. Matching Hamiltonians in either case does not establish conductivity equivalence.

The export evaluates each HOH path with the actual intermediate site's overlap block. Some native DOS and exchange recursion kernels use precomputed products from representative neighbor slots. In finite inhomogeneous systems those slots can have a different destination type. The imported DOS, exchange and conductivity routes all use the same exported matrix.

`conductivity_import` uses the existing `linear_in` and `linear_out` selections, polarizations, and velocity directions. The local accumulations and torque operators use the same Kubo machinery. For spin torque or SOC torque on a spin-orbit-coupled Hamiltonian, use a new export with onsite SOC decomposition.

`exchange_import` uses the native LMTO pair-exchange driver. With the default `&calculation exchange_delta='lmto'`, keep the same frozen atomic potential files, including any referenced sidecars, so the energy-dependent exchange vertices refer to the state from which H was exported. For a collinear `spd` TB model, choose `exchange_delta='onsite'` to use the full real orbital matrix `Re(H_upup - H_downdown)` from each site's zero-image onsite block, following the `exchange_p2rs` convention. This assumes the spin quantization axis is the supplied up/down axis. With SOC, the formula uses the full onsite block, including any spin-dependent diagonal SOC contributions. Set `njij` and each `ijpair(i,:)` in `&lattice` to the actual cluster indices in that replay, and use `njijk=0`. Pair selection can be changed after export. The standard `jij.out`, `dij.out`, `aij.out` and tensor/decomposition outputs are written. DOS and conductivity use the site layout independently of the pair list.

Use the existing native reciprocal drivers for k-space observables. The standalone `paoflow2rs`, `conductivity_p2rs` and `exchange_p2rs` workflows remain available for external files.

### Start from an external TB Hamiltonian

Use a fresh run directory and the standard lattice and atom-type inputs.
For custom geometry, set `&lattice crystal_sym='file'` in `input.nml` and
supply the geometry in `lattice.nml`. The external matrix must describe the
explicit cluster produced by that deck. In a closed direction, image
translations are counted in **whole cluster-cell vectors**. A primitive-cell
hopping file can be used directly when the generated cluster is that cell.
Provide both directions of the Hermitian hopping and use an orthonormal,
atom-centered basis with the same complete orbital shells per site.

For example, these settings run a collinear pair-exchange calculation from
`paoham.dat`. Keep the rest of the matching geometry and atom-type deck.
Choose cluster site numbers for the pair and set the Fermi energy relative
to the supplied matrix's energy zero.

```fortran
&calculation
  post_processing = 'exchange_import'
  exchange_delta = 'onsite'
/
&hamiltonian
  effective_import = 'external'
/
&control
  nsp = 1
  lmax = 2
  recur = 'chebyshev'
/
&energy
  fermi = 0.0  ! Ry
/
```

Set `njij=1` and, for example, `ijpair(1,:)=1,2` in `&lattice`. For DOS
or conductivity, select `dos_import` or `conductivity_import` and omit
`exchange_delta`. The external mode derives the site map from the cluster
and writes `rs_tb.meta` when absent. In subsequent runs,
`effective_import='legacy'` reads the generated sidecar and verifies its
hopping checksum. External mode uses the current deck and leaves an existing
sidecar in place. Use a fresh directory when changing the lattice or payload.

The default flat indices follow spin, site, then the legacy Cartesian
orbital order. For a different index permutation of those same orbitals,
supply `rs_tb_orbital_map.dat` with one line per index in increasing
external-index order. Each line has two one-based integers,
`external_index rs_flat_index`. The map must cover all orbitals exactly
once. The generated sidecar stores it for later replay.

Conductivity constructs the velocity from the mapped atomic positions and
the hopping images. Charge and compatible spin/orbital operators use the
existing response code. A seven-column file alone carries no separate
onsite SOC decomposition for the SOC torque operators.

## 4. File conventions

Both real-space payloads are sparse scalar records in eV with one-based indices. The physical matrix is unscaled and its energy zero is unchanged. The metadata records lattice vectors and site positions in Å, spin and orbital ordering, boundary flags, provenance and the source Fermi level in Ry. Keep payloads and metadata from the same run. V3 metadata fingerprints the exact completed payloads, including whitespace. The reader also accepts older v2 metadata without checksums.

```text
rs2paoham.dat or rs_tb_paoham.dat
R1 R2 R3 flat_row flat_col real_eV imag_eV

rs_tb.tb
R1 R2 R3 row_site column_site record_index row_orb col_orb flat_row flat_col real_eV imag_eV
```

The column site is in the image $R$ relative to the row site. In cluster data, every open direction has $R=0$. The displacement from row to column is $\mathbf r_j-\mathbf r_i+A R$, where $A$ uses the cell vectors in the file's metadata. A nonzero image of the same site is a hopping. A Fourier transform in the native position gauge uses the positive phase of this displacement. The canonical file uses spherical orbitals ordered by increasing $l,m$ within each spin block. The seven-column file uses the legacy Cartesian orbital order. Flat indices enumerate all sites in one spin block before the next. `record_index` numbers the blocks in each row with zero-image onsite first, then by column site and image. Exact-zero scalar entries are omitted.

`rs_tb.kspace` is a text file with fractional reciprocal coordinates and weights on each `k` line, followed by one-based `H row col real imag` records in Ry and optional dimensionless `S row col real imag` records. Its matrix order is site, spin, then $l,m$. The header identifies the spin frame and, for GBT, the spiral vector in Cartesian units of $2\pi/a$ and the cone angle in radians. The footer reports the largest native Hermiticity residual and the number of k points. The writer preserves the native matrix as assembled, including any numerical asymmetry.

In native real space, `hoh` selects the built first- or second-order Hamiltonian. First order uses the native `ee` or `hall` with the applicable Hubbard term already assembled there, plus the built SOC and CCOR terms. HOH uses the native $h_{ij}-\sum_k h_{ik}o_k h_{kj}$ plus the built onsite `enim` and the applicable SOC and CCOR terms. The exporter uses saved global-frame blocks for real-space local-axis runs.

In reciprocal scope, `&reciprocal kspace_ham_order` and `reciprocal_mode` configure the native assembler. With order `auto`, the ordinary assembler follows the built HOH state. For explicit `kspace_ham_order='second'`, also set `hoh=.true.`. With HOH off, the assembler falls back to first order and the file header still records the requested order. Exact orthogonal output uses its native `strux_lib` route. The overlap proxy and exact orthogonal mode are mutually exclusive.

## 5. Worked examples

### A. Frozen bulk bcc Fe to `rs2paoham.dat`

Copy [the bcc Fe deck](../example/bulk/bccFe). Keep the example structure and atom files. Change these settings.

```fortran
&calculation
  pre_processing = 'none'
  processing = 'none'
  post_processing = 'hamiltonian_export'
/
&hamiltonian
  hoh = .true.
  export = 'rs2pao'
  export_scope = 'periodic'
/
```

Run `rslmto.x input.nml > run.log 2>&1` inside that directory. The output includes the full effective HOH hopping matrix and the metadata. The original deck's recursion choice does not change the exported matrix. For a compact periodic conductivity replay, use a separate fully closed matching replica and the settings in step 3.

### B. Mixed-boundary Cu with native and replayed conductivity

Copy [the fcc Cu conductivity deck](../example/conductivity/fccCu). For this small verification run, retain the example geometry and atom potential and apply these overrides.

```fortran
&calculation
  pre_processing = 'none'
  processing = 'none'
  post_processing = 'hamiltonian_export'
/
&lattice
  pbc = .true.
  b1 = .true.
  b2 = .true.
  b3 = .false.
  n1 = 4
  n2 = 4
  n3 = 4
/
&control
  nsp = 1
  recur = 'chebyshev'
  lld = 16
  cond_ll = 16
  cond_calctype = 'per_type'
  linear_in = 'charge'
  linear_out = 'charge'
  cheb_backend = 'legacy'
  transport_backend = 'legacy'
/
&hamiltonian
  hoh = .true.
  ccor_2c = .false.
  export = 'rs2pao'
  export_scope = 'cluster'
/
&energy
  fermi = -0.3
  energy_min = -2.0
  energy_max = 2.5
  channels_ldos = 400
  fix_fermi = .true.
/
```

Run the export deck first. In a separate copy, set `post_processing='conductivity'`, `gf_route='recursion'` and `export='none'` to compute the native reference. In a third directory, copy the input and atom files and place `rs2paoham.dat` and `rs_tb.meta` from the export job there. For replay, change the workflow and file selection as follows.

```fortran
&calculation
  pre_processing = 'none'
  processing = 'none'
  post_processing = 'conductivity_import'
/
&hamiltonian
  export = 'none'
  effective_import = 'legacy'
/
```

Keep the original `v_alpha`, `v_beta` and other response settings in the namelist when making these edits. Compare `response_total.out`, `response_orbitals_real.out`, `response_orbitals_imag.out` and `Cu_dos.out`. This low-order deck checks replay. Converge the replica and recursion orders separately for production.

### C. DOS only on the same Cu geometry

For native block DOS, use a fresh copy of the native Cu deck above. Set `post_processing='dos_frozen'`, `recur='block'`, and `export='none'`. The corresponding native Chebyshev run uses `recur='chebyshev'`.

For imported DOS, use `post_processing='dos_import'`, `recur='chebyshev'` or `'block'`, `export='none'` and `effective_import='canonical'`. Copy `rs_tb.meta` and `rs_tb.tb` from example B and preserve the cluster geometry and energy mesh settings. The result is `Cu_dos.out` in each job. Use the same recursion method and Fermi reference when comparing native and imported spectra.

### D. Exact orthogonal bcc Fe in reciprocal space

Start from the bcc Fe deck in example A. Keep the frozen `&calculation` block. Change these settings to write an eight-point diagnostic mesh.

```fortran
&control
  calctype = 'B'
  nsp = 1
/
&lattice
  strux_backend = 'strux_lib'
/
&hamiltonian
  hoh = .false.
  ccor_2c = .false.
  export = 'rs2pao'
  export_scope = 'kspace'
/
&reciprocal
  kspace_ham_order = 'orthogonal'
  reciprocal_mode = 'ham_only'
  nk1 = 2
  nk2 = 2
  nk3 = 2
/
```

Use a scalar-relativistic frozen potential and leave Hubbard corrections off for this example. Inspect `rs_tb.kspace` for `matrix_size`, the basis map and `end_kpoints 8`. To write the generalized overlap proxy instead, use a compatible fixed-potential input with `reciprocal_mode='generalized_overlap_proxy'` and a non-orthogonal order such as `kspace_ham_order='auto'`. The file then carries both $H(k)$ and $S(k)$. For an existing legal `magnetic_representation='gbt_single_q'` deck, retain the native $q$, spin and structure settings and choose `export_scope='kspace'`.

## 6. Check a run

| Symptom | Check |
| --- | --- |
| No real-space file | Check that the selected driver writes exports, `export` is enabled, and the resolved scope is real space. `dos_frozen` only writes DOS. |
| Geometry or boundary mismatch | Compare the source and destination basis, labels, site order, cell vectors and boundary flags to `rs_tb.meta`. For an exact open calculation, re-export using `cluster`. |
| Payload checksum mismatch | Copy the payload and metadata from the same run. With `effective_import='legacy'`, check which of the three filenames is found first. |
| Ambiguous wrapped neighbor or inconsistent reciprocal slot | Enlarge the source replica so the native neighbor images are resolved. |
| Native and replayed responses differ | First match Fermi references, energy mesh, recursion settings and current definitions. Check CCOR, nonuniform weights and native Hermiticity separately. |
| A reciprocal reference test fails with symmetry disabled | Use the same symmetry-enabled environment as the references. Spglib 2.6.0 passed the CPU quick suite. Spglib 1.6.3 lacks a magnetic-symmetry symbol required by this revision. |

See the [validation report](../HAMILTONIAN_EXPORT_REPORT.md) for measured errors and the [validation instructions](../tests/hamiltonian_export/README.md) for automated and independent comparisons.
