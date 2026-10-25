# k-space GBT spin spirals

Generalized-Bloch-theorem spin spirals on the k-space path (`use_kspace`, `magnetic_representation = 'gbt_single_q'`), 32^3 Gamma-centred mesh, 300 K, `kspace_ham_order = 'second'`, no SOC (`nsp = 3`). q must have q/2 on the mesh: with 32^3 every direct component is a multiple of 1/16. All decks start from the converged ground-state potentials of `tests/spin_response/{fe_bcc,ni_fcc}`. Small-q results need much denser meshes.

| directory | method | run |
| --- | --- | --- |
| `bccFe/mft` | force theorem along Gamma-H-N-Gamma-P-H (33 q) | `rslmto.x input.nml` |
| `fccNi/mft` | force theorem along Gamma-X-W-K-Gamma-L (29 q) | `rslmto.x input.nml` |
| `fccNi/scf_unconstrained` | unconstrained cone SCF (theta = 30 deg, no constraining field), 9 q along Gamma-X | `./run.sh /path/to/rslmto.x` |
| `fccNi/scf_sweep` | the same unconstrained SCF as `scf_unconstrained`, all 9 q in one run (`post_processing = 'frozen_magnon'`, `mode = 'scf'`) | `rslmto.x input.nml` |
| `fccNi/scf_constrained` | constrained cone SCF (theta = 30 deg, `constraints_mom_ref = 0,0,1`), 9 q along Gamma-X | `./run.sh /path/to/rslmto.x` |

**mft.** `post_processing = 'frozen_magnon'`, `mode = 'mft'`, `nstep = 0`, `theta_ss = 5`. The q path is `q_path.dat` (`q_coordinates = 'direct'`, q = 0 first; q/2 off the mesh stops the run and prints the nearest allowed q). Output `frozen_magnon_kspace.dat`, one row per q: q direct (3), q Cartesian in 2pi/alat (3), theta [deg], Omega(q,theta) [Ry], Omega(q,0) [Ry], N(q,theta), omega [Ry], omega [meV], N(q,0). omega = 4 (Omega(q,theta) - Omega(0,0)) / (M sin^2 theta) uses the summed site moment M and is meaningful for one magnetic site per cell. The potentials are the collinear ground states of tests/spin_response (32^3 mesh, the deck's own Hamiltonian settings); they were not re-converged with `kspace_ham_order = 'second'`.

**scf.** `input.template` holds `q_ss = 0, 0, @QZ@` (Cartesian, 2pi/alat); both variants use the same theta, so the energy difference between the two directories at each q is the effect of the constraining field. `run.sh` writes one directory per q (`QZS="0.0 0.5" ./run.sh ...` runs a subset) and prints the total energy of each. The energy difference to q = 0 is the spiral energy. With `post_processing = 'frozen_magnon'` and `mode = 'scf'` (no `&constraints`) the same SCF runs over a `q_path.dat` in one run; the constrained variant has no such sweep on the k-space path.

Measured at 32^3 on 8 cores (Release): bcc Fe mft 94 s, fcc Ni mft 87 s, fcc Ni SCF 90-100 s per q (about 15 min per directory for the 9 q), `scf_sweep` 763 s for the 9 q.
