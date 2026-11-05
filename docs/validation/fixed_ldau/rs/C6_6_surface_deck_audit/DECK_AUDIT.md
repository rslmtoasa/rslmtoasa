# C6.6 Ni surface deck audit

The two archived `input.nml` files were copied from the `fccCu001` fixture and reconstructed
from its exact bytes. Their hashes match the deck hashes captured by the checked preflight before
the solver rewrote runtime namelists. `diff -u` against the fixture showed only the lines below;
geometry controls (`nlay`, `surftype`, `ntype`, cutoffs, and all generated-coordinate inputs) were
unchanged.

| Fixture line | Staged line | Reason |
|---|---|---|
| `alat = 3.614` | `alat = 3.52340` | Accepted fcc Ni bulk lattice constant |
| `wav = 1.41237` | `wav = 1.37692` | Accepted fcc Ni bulk Wigner-Seitz radius |
| `label(1) = 'Cu'` | `label(1) = 'Ni-imp'` | Match the accepted Ni frozen host potential |
| `label(3) = 'Cu-S'` | `label(3) = 'Ni-imp-S'` | Match the Ni surface recursion seed |
| `nstep = 1` | `nstep = 2` | User-approved maximum two-step response probe |
| (none) | `conv_thr = 5.0d-8` | User-approved response-probe tolerance |
| `fermi = -0.088376` | `fermi = -0.061141000000` | Accepted Ni host fixed-Fermi endpoint |
| `energy_max = 1.2` | `energy_max = 0.5` | Accepted Ni host energy window |
| `nsp = 2` | `nsp = 1` | Required surface local-moment spin setting |
| `beta = 0.0020` | `beta = 0.1` | Accepted same-material Ni bulk mixing parameter |

The unchanged `energy_min = -1.0`, `mix` type, geometry, and ES definition already match the
host/fixture requirements. The appended Hubbard lines, each with its reason, are:

| Appended line(s) | Deck | Reason |
|---|---|---|
| `&hamiltonian` / `/` | both | Add only the Hubbard settings required for this surface test. |
| `hubbard_u_general(1,3) = 0.0`; `hubbard_j_general(1,3) = 0.0` | U=0 control | Explicit U=0 on frozen Ni host type 1 d channel. |
| `hubbard_u_general(3,3) = 0.0`; `hubbard_j_general(3,3) = 0.0` | U=0 control | Explicit U=0 on active Ni surface type 3 d channel. |
| `hubbard_u_general(1,3) = 3.0`; `hubbard_j_general(1,3) = 0.7` | U>0 test | Apply requested U/J to frozen Ni host type 1 d channel. |
| `hubbard_u_general(3,3) = 3.0`; `hubbard_j_general(3,3) = 0.7` | U>0 test | Apply requested U/J to active Ni surface type 3 d channel. |
| `hubbard_radial_moments = 'auto'` | U>0 test | Exercise surface auto selection using the local U-carrying host. |

The first type is frozen Ni; the third type is active Ni surface recursion. ES remains the copied
fixture `ES.nml` (SHA256 `3dc19c818e6da0d23ad9bccb2c8ecc38e7f96d1fcb837b96a6c84d6cab60c524`) and has no U/J. The accepted host potential was copied byte-for-byte to `Ni-imp.nml`; the surface seed was copied from the same host potential with only the namelist `symbol` changed to `Ni-imp-S`.

| Deck | SHA256 (preflight-bound) | Host restart SHA256 | Rekeyed surface seed SHA256 |
|---|---|---|---|
| [`T2a_U0/input.nml`](T2a_U0/input.nml) | `21dfb570b87698f96d244f7adfa8716324cd506067e9af44e15cc3ba6fdecd4a` | `0213bf83b07bd82446f293bf19be8ba37b0a2bc816011d29f2c7e8d5c483738a` | `bf0c5ff385ea2b766729b145203e0e61a06246436afa972a7156e9c560d9b901` |
| [`T2b_U3J07_auto/input.nml`](T2b_U3J07_auto/input.nml) | `d69fe0ac32fdec6f247a6ba0fdf1dfb664d126e69f969b7fea0758f6e56d3a30` | `bc3340a8ad20d2467eaeccd43929a659d4cfff2a9f28b748fb6554edc4ee3ac9` | `2e444b62ba9f8303f033b0cb3a77e630025284d182ae03235c367af4c448c7b9` |


The accepted U=3/J=0.7 host deck used to derive the second host-EF check is archived as
`T2b_U3J07_auto/accepted_host_input.nml` (SHA256
`c0722a059758b3d52387f918be640e34ed10f11fd0e9cab06cf5e10a6e0b79f2`). It has
`fermi=-0.061141` and `fix_fermi=.true.`. The U=0 accepted host input is already versioned at
`C5_2_PART_D_assets/ni_u0_bulk_rc80_host/inputs/input.nml`. Each T2 artifact directory includes
the checked `preflight.json`, its host-EF CSV, solver log, report, output potentials, and layer map;
`SHA256SUMS` covers the archived files.
