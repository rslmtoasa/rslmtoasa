# Supplied hcp Co Chebyshev RS reference validation

This record documents two bounded validations of the user-supplied converged hcp Co reference. It does not resume C5 or change the stopped no-U RS byte-identity criterion; the Co sizing fix remains unmerged.

| Run | Executable | Result |
|---|---|---|
| Frozen Release crosscheck | `710f91b3…` (source commit `2686e4e`) | Exit 0; converged at step 16, diff 4.8e-9 Ry. |
| Producer Debug replay | `e4295134…` (`build/bin/rslmto.x`) | Exit 0; converged at step 25, diff 4.8e-9 Ry. |

Both runs used the unchanged B/HCP `strux_lib` spectral RS deck: Chebyshev, U=J=0, nsp=1, lmax=2, rc=60, lld=150, EF=-0.032603 Ry, energy window [-1, 1.3] Ry, Broyden beta=0.15, nstep=25, conv_thr=5e-9 Ry, and one OMP/BLAS thread. The supplied `Co1.nml` and `Co2.nml` were copied byte-for-byte to the required `_out.nml` restart aliases in isolated scratch directories; original inputs were unchanged. The historical supplied `out.scf` reports restart-loaded potentials and convergence at step 8 (diff 4.7e-9 Ry).

The four original input hashes, executable and build provenance, commands, process outcomes, endpoint comparisons, finite-output checks, and per-file raw-output hashes are in [`C5_CO_REFERENCE_PROVENANCE.json`](C5_CO_REFERENCE_PROVENANCE.json) and [`C5_CO_REFERENCE_assets`](C5_CO_REFERENCE_assets/). Logs and final potentials are retained there. A producer-copy mode error caused one pre-exec exit 126; no solver started in that attempt. Only the frozen scratch copy was chmodded 0555, retaining its SHA; the single authorized producer replay then completed with exit 0.

**Interpretation (4 lines):**

Both independent-build reference validations converged at the unchanged deck threshold and produced finite checked potentials and reports.
The producer endpoint's per-spin d-occupation and reported energy/moment comparisons are descriptive records, not new tolerances.
These reference checks do not overturn C5 STOP2.5(a), establish C5 acceptance, or adopt any numerical change.
No source, production binary, supplied input, or C5 gate criterion was changed; no further run was launched.
