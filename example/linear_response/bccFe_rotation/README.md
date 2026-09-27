# bcc Fe rotation response
Run this deck from the directory containing `input.nml`.
It uses `post_processing='linear_response'`.
The formulation is `rotation`.
The reciprocal production mesh is 12x12x12.
The q path contains Gamma and one finite direct-coordinate point.
The accepted k-space SCF state is consumed without a rebuild.
The spectral finite-H response is selected explicitly.
The main output file is `rotation_dynamics.dat`.
Optional frequency-grid output is controlled by the rotation keys.
