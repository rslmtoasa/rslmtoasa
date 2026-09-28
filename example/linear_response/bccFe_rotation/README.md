# bcc Fe rotation response
Run this deck from the directory containing `input.nml`.
It uses `post_processing='linear_response'`.
The formulation is `rotation`.
The reciprocal demonstration mesh is 4x4x4; it is intentionally not a
converged Fe magnon calculation.
The q path contains Gamma and one finite direct-coordinate point.
The accepted second-order (`HOH`) k-space SCF state is consumed without a
rebuild. The deck is scalar-relativistic, collinear, and SOC-free.
`Fe.nml` is the canonical bcc-Fe fixture reused from
`example/exchange_q/bccFe/Fe.nml`.
The spectral finite-H response is selected explicitly.
The main output file is `rotation_dynamics.dat`.
Optional frequency-grid output is controlled by the rotation keys.
