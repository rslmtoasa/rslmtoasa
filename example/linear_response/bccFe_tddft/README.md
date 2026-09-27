# bcc Fe linear-response TD-DFT
Run this deck from the directory containing `input.nml`.
It uses `post_processing='linear_response'`.
The formulation is `tddft` with the compact product representation.
The bare response is the Lehmann service.
The interaction is direct ALSDA.
The accepted reciprocal SCF state is reused by the response run.
The q and frequency grids are defined in `input.nml`.
The main response file is `tddft_Fe_smallq_dispersion.dat`.
State and provenance sidecars use the same output stem.
