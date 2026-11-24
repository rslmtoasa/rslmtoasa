# Effective Hamiltonian export and import

Use a frozen `*_out.nml` potential or an initial `*.nml` deck with the same native post-processing semantics. Keep the density-matrix sidecars referenced by the potentials. The frozen export route builds the operator and does no SCF iterations. The matrix is assembled from native hopping blocks and does not depend on the recursion algorithm. A Chebyshev build may also supply spectral bounds in the metadata.

```fortran
&calculation
  pre_processing = 'none'
  post_processing = 'hamiltonian_export'
/
&hamiltonian
  hoh = .true.
  export = 'rs2pao'
  export_scope = 'cluster'
/
```

The export options also work at the end of the standard SCF preprocessing routes and during native real-space conductivity. An SCF export describes the last built operator. Use the frozen route when the reference is the set of potentials just saved to disk.

| `export_scope` | Output | Interpretation |
|---|---|---|
| `periodic` | `rs_tb.meta`, `rs_tb.tb`, seven-column file | Bulk cell with one recursion class per basis site. Exact integer column-image translations. |
| `cluster` | `rs_tb.meta`, `rs_tb.tb`, seven-column file | Actual finite bulk, surface, interface or impurity graph, including wrapped images. Use this for open or mixed boundaries and explicit textures. |
| `kspace` | `rs_tb.kspace` | Native complex matrices on the full reciprocal mesh. Includes S(k) in `generalized_overlap_proxy`. Honors `kspace_ham_order` and native restrictions on GBT and orthogonal assembly. |
| `auto` | Varies | Compact periodic for ordinary bulk, cluster for mixed-boundary bulk, surface, interface, impurity and explicit texture, reciprocal for GBT single-q. Select `kspace` explicitly for orthogonal and overlap modes. |

For the finite graph a path contributes only if both edges exist in the native graph. Cropping a finished compact periodic Hamiltonian generally cannot reproduce open boundaries. An ordinary reciprocal export requires a replica large enough to supply consistent representative bond vectors. The reciprocal writer checks this for wrapped replicas and reports any matrix Hermiticity defect it measures.

The first-order operator includes the built `ee` or `hall`, onsite SOC and additive CCOR. For HOH the two original directed blocks contribute `-h_ik o_k h_kj` at the sum of their unwrapped image indices. The exporter uses the stored `enim`, which contains fixed onsite Hubbard U when the native HOH builder routes it there. The first-order builder has already added U to its onsite hopping. Local-axis real-space runs export the saved global-frame blocks.

## Real-space files

Canonical scalar records use eV and one-based indices. The metadata records positions and cell vectors in Angstrom, boundary flags, orbital order, spin and Hubbard settings. V3 metadata includes checksums for the exact canonical and seven-column file bytes. Keep both files from the same run. The reader also accepts earlier v2 metadata without a checksum. The numeric fields in `rs_tb.tb` have this order.

```text
R1 R2 R3 row_site column_site record_index row_orb col_orb flat_row flat_col real_eV imag_eV
```

The record index starts with the zero-image onsite block, followed by lexicographically sorted column and image keys. Exact-zero scalar elements are omitted, so some indices may be absent. The canonical basis uses native spherical harmonics in increasing l and m with a contiguous block per spin. Flattened indices list all sites in the first spin block before the second.

In either real-space scope, `export='rs2pao'` writes `rs2paoham.dat` and `export='python'` writes `rs_tb_paoham.dat`. The format is `R1 R2 R3 i j real_eV imag_eV` in the existing Cartesian orbital basis with spin-major global numbering. With `periodic`, the indices refer to compact-cell orbitals and R counts compact-cell translations. With `cluster`, the indices refer to orbitals on individual cluster sites and R counts whole replica-cell translations along closed directions. The open-direction components of R are zero. Keep the matching `rs_tb.meta` to preserve geometry, boundary flags and file provenance. External readers must interpret those boundaries explicitly for cluster data. The canonical file uses native spherical orbitals. Export uses one MPI writer. Use a separate working directory for each run.

Reciprocal export stores matrices indexed by k in `rs_tb.kspace`, with the overlap metric where required. The seven-column layout describes real-space hoppings and has no fields for a k mesh, spin-spiral gauge or overlap metric.

## Import as ordinary tight binding

Use the same basis, `nsp`, lattice, atom labels and site ordering. For a compact periodic file use a fully periodic replica with `pbc=.true.` and all three boundary flags enabled. A finite export must match the source cluster and its boundary flags.

```fortran
&calculation
  pre_processing = 'none'
  post_processing = 'conductivity_effective'
/
&control
  recur = 'chebyshev'
  cheb_backend = 'legacy'
  transport_backend = 'legacy'
/
&hamiltonian
  effective_import = 'canonical'
/
```

Select `dos_effective` for real-space DOS. The matching native fixed-potential DOS route is `dos_frozen`, which supports `recur='chebyshev'` and `recur='block'`. The imported DOS route currently uses Chebyshev. For a seven-column round trip, set `effective_import='legacy'` and keep the Cartesian file with its matching `rs_tb.meta`. The reader checks for `paoham.dat`, then `rs2paoham.dat`, then `rs_tb_paoham.dat`. The existing `paoflow2rs`, `conductivity_p2rs` and `exchange_p2rs` routes keep their external-file workflow.

Effective import reconstructs every nonzero hopping shell and keeps separate wrapped images. The generic non-HOH kernels consume those blocks with the file's Fermi energy. Export-only records the input Fermi level. When native conductivity exports, it records the final Fermi level used in the response and the initial Fermi reference used to build the DOS mesh. The final value can equal the input if the native determination leaves it unchanged. Keep the energy limits and number of mesh points the same for a paired run. Charge, spin and orbital responses use the existing reconstruction. The import route rejects torque and reciprocal-only input.

For an unweighted or uniformly weighted current with CCOR off, the ordinary effective velocity agrees with the native HOH action because the two-hop displacement adds along each path. Nonuniform per-type weights enter each native HOH leg separately. With CCOR on, the ordinary effective velocity includes the offsite CCOR commutator omitted by the current native real-space velocity.
