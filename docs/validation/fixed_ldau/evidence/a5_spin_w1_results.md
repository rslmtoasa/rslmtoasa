# A5.3 spin-selective W1 response results

The preregistered 14 cases completed from the common RS-origin restart. Every case returned 0 and reported literal convergence, the native reciprocal-KS branch, and finite entropy/TS. The captured density is the actual `potential%ldm` at P_OUT; both spin q values are retained, and inactive-spin applied Tr(WP) is zero. The complete per-case matrices, thermodynamics, input/seed/output hashes, and run metadata are in [`a5_spin_w1_results.json`](a5_spin_w1_results.json).

At primary `h=0.0025 Ry`, `B=[+0.00178134,-0.01104575] Ry/Ry`, with registered component uncertainties `[0.00166546,0.00186671]`; the sum is `-0.00926441 ± 0.00353217`. This sum differs from the historical A5.1 local aggregate residual by `-4.95e-5`; that comparison is diagnostic only because A5.1 used an eigen-reconstructed density.

The response matrix, with rows indexed by output spin occupation and columns by perturbed W1 spin, is `R=[[-2.42485,+1.88811],[+1.86883,-2.51948]] electron/Ry`. Its measured componentwise uncertainty is `[[0.43435,0.48420],[0.44496,0.49976]] electron/Ry`, dominated by the registered `h=0.005` versus `0.0025` step sensitivity. The nominal solve `R^T g=B` gives `g=[0.0062595,0.0090751] Ry/electron`, determinant `2.58080`.

Enumerating the preregistered component uncertainty box for B and R gives determinant corner values from `-1.46879` to `+6.63362`; the allowed set therefore contains singular response matrices. The gradient and chain-rule allowance are unresolved under the registered uncertainty model. The nominal local `N_ASA` from `dq/dU=[+0.08327349,-0.08073297] electron/eV` is `0.00125391 Ry/eV`, while the first-order propagated W1 uncertainty alone is `0.00246649 Ry/eV`; these nominal values are descriptive and are not a decision allowance.

The external-field free energy was reconstructed from the same final state as `F=ETOT+lambda*Tr(WP)-TS`. The runtime external-energy accumulator omitted the lambda term in these logs; the logger/reconstruction discrepancy and raw trace remain recorded in the JSON. This is an extraction correction only; SCF inputs and physical endpoints are unchanged.

Provenance: executable SHA-256 `52eb94c6e6f6ae12fca629a4eea1b30a6800e29736c9fdfdfd836538077f2225`; common input `ca65d9695b049f10beb3ee9461c05572b48f3973885de99ea7f8fee3dfc5699c`; common restart `25c1c8260d8996360812d978aec09986c4ec1e853c3cf17635e6c62b4c62fec5`; summary JSON SHA-256 `3d61a0c4e114d51a5acbf2039ff82ac711de08f6719856ad50a40b15307b61dd`.
