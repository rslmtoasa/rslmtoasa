# C4 merged verification

The requested no-fast-forward merge is `48787e5c71b3d337e6f364ddac6ffe9d3df8f450` (parents `18972d6` and `6b06b8f`). A clean isolated Debug build completed successfully. The seven targeted LDA+U/local-moments unit tests passed.

The merged executable (`df7980df…`) and frozen clean-8b executable (`8307387e…`) were run from byte-identical per-mode snapshots with `OMP_NUM_THREADS=1`, starting from fresh atomic `Fe.nml` potentials. The 26 RS and 23 KS physical output files matched byte-for-byte; only the recorded runtime-only files were excluded.

An earlier attempted comparison restarted from generated `Fe_out.nml` and is explicitly excluded. The corrected paired run used shared snapshots and confirmed fresh-atomic startup for all four runs. Exact commands, input and binary hashes, raw logs, test output, and file-level comparisons are in [C4_MERGE_PROVENANCE.json](C4_MERGE_PROVENANCE.json).
