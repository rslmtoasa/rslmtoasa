# C5.6 T2a — KS Fe U=2/J=0 default-mode host regeneration

**Verdict: FAIL before SCF; KS acceptance unevaluated.** The single authorized T2a invocation used the current production binary and exact accepted legacy Fe state. The run resolved omitted `hubbard_radial_moments` to `local` for native KS, then stopped when the unchanged U-carrying restart was interpreted as legacy `spectral`. No SCF iteration ran, no qL/Diff window exists, and no regenerated host was written.

The runtime log records `HUBBARD_RADIAL_MODE mode=local source=default use_kspace=.true. local_energy_adopted=.true.`, followed by the warning that the legacy potential has no mode field and is treated as spectral, then a fatal mode mismatch at `source/self.f90:2392`. The staged `Fe_scf.nml` and `Fe_out.nml` stayed byte-identical to the protected accepted checkpoint (`bb647254…ead6b`). The original accepted files were not modified.

Preflight passed for the bcc Fe native-KS recipe, exact accepted restart, parsed namelists, and single-thread caps. Full-trace watcher and disk guard were running before solver launch. The watcher captured zero iteration snapshots because failure occurred before the first SCF loop; the 1 GiB disk floor was never approached. The production binary SHA256 is `4554ea2b…fdc22e`; input, source, process, disk, and raw artifact hashes are in [C5_6_T2A_PROVENANCE.json](C5_6_T2A_PROVENANCE.json). The compressed raw packet is [t2a_ks_u2_auto_failed_raw.tar.gz](C5_6_T2A_assets/t2a_ks_u2_auto_failed_raw.tar.gz) (SHA256 `49c1562e…17cdec4`).

C5.6 campaign use is **1/2**, including this failed start. T2b was not launched because T2a did not produce an accepted serialized host. No retry, metadata edit, or source change was made; work stops pending the sole batched audit.
