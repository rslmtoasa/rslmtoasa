# k-space GBT spin spirals

Generalized-Bloch-theorem spin spirals on the k-space path (`use_kspace`, `magnetic_representation = 'gbt_single_q'`), 16^3 Gamma-centred mesh, 300 K, `kspace_ham_order = 'second'`, no SOC (`nsp = 3`). q must have q/2 on the mesh: with 16^3 every direct component is a multiple of 1/8. All decks start from the converged ground-state potentials of `tests/spin_response/{fe_bcc,ni_fcc}`. 16^3 is for a quick look; small-q results need much denser meshes.

| directory | method | run |
| --- | --- | --- |
| `bccFe/mft` | force theorem along Gamma-H-N-Gamma-P-H (17 q) | `rslmto.x input.nml` |
| `fccNi/mft` | force theorem along Gamma-X-W-K-Gamma-L (15 q) | `rslmto.x input.nml` |
| `fccNi/scf_unconstrained` | flat spiral SCF (theta = 90 deg), one run per q along Gamma-X | `./run.sh /path/to/rslmto.x` |
| `fccNi/scf_constrained` | constrained cone SCF (theta = 30 deg, `constraints_mom_ref = 0,0,1`), one run per q along Gamma-X | `./run.sh /path/to/rslmto.x` |

**mft.** `post_processing = 'frozen_magnon'`, `mode = 'mft'`, `nstep = 0`, `theta_ss = 5`. The q path is `q_path.dat` (`q_coordinates = 'direct'`, q = 0 first; q/2 off the mesh stops the run and prints the nearest allowed q). Output `frozen_magnon_kspace.dat`, one row per q: q direct (3), q Cartesian in 2pi/alat (3), theta [deg], Omega(q,theta) [Ry], Omega(q,0) [Ry], N(q,theta), omega [Ry], omega [meV], N(q,0). omega = 4 (Omega(q,theta) - Omega(0,0)) / (M sin^2 theta) uses the summed site moment M and is meaningful for one magnetic site per cell. The potentials were converged on a denser mesh, so the force theorem is applied to a potential that is not self-consistent on 16^3.

**scf.** `input.template` holds `q_ss = 0, 0, @QZ@` (Cartesian, 2pi/alat); `run.sh` writes one directory per q (`QZS="0.0 0.5" ./run.sh ...` runs a subset) and prints the total energy of each. The energy difference to q = 0 is the spiral energy. There is no k-space `mode = 'scf'` loop over q.

Measured at 16^3 on 8 cores (Release): bcc Fe mft 12 s, fcc Ni mft 13 s, fcc Ni SCF about 25 s (unconstrained) and 17 s (constrained) per q.
