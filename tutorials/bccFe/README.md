# bcc Fe with an explicit lattice

This starting deck uses `crystal_sym = 'file'` and `lattice.nml` for the
one-site primitive bcc cell. The vectors reproduce the built-in bcc preset;
`alat = 2.86120` Angstrom remains in `input.nml`.

- `input.nml`: bulk SCF starting input. Apply the lesson's iteration/mixing edits in the editor.
- `input.exchange.nml`: complete two-pair exchange starting input; copy to `input.nml` in a separate working directory and match the final SCF settings.
- `lattice.nml`: shared explicit geometry and site/type mapping.
- `Fe.nml`: original bulk Fe atom input.

Copy the entire directory under `runs/` before running. Transfer the SCF
checkpoint and its matching final input settings to the exchange job. Use
the same `lattice.nml` for SCF, exchange, and the Fe UppASD export.
The export CLI's lattice-only path retains a common lattice-parameter unit;
it does not apply `alat` automatically. Review the target's length scale.

The selected exchange pairs are representative bonds, not a complete
exchange-shell table or a production spin model.
