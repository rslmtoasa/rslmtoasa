S0 (user-approved; record verbatim in docs/validation/fixed_ldau/S0_AMENDMENT.md).
OBJECTIVE: janitor pass before SOC. No physics change.
Base: fable_v4b_u at 90f57c3. Work branch: fix/ldau-s0-janitor.
Token rules: Sol plans in at most 1/2 page; ONE Luna worker sequentially; ONE batched audit at the
end. Python: .venv. Launch everything through tools/run_rslmto_checked.py.
Budget: at most 12 SHORT runs (4 SCF steps each), no SCF campaigns.

IDENTITY GATE (used by every task; all runs single-threaded, 4 steps). Before/after binaries must
give byte-identical scientific output files for:
  - no-U: RS bulk/bccFe and KS k_space_scf/bccFe;
  - +U: bulk/bccFe_hubbard_nsp1 (RS, auto -> local) and the KS bccFe deck with U=2;
  - +U impurity: impurity/B2FeCo_nsp1.
Exceptions: logs, and the T3 file-format change, which is checked by value as specified there.
Also: all LDA+U unit tests, the 13 guard fixtures, and the impurity and surface auto fixtures
pass.

T1 Diagnostics extraction
   - Move the fixed-U trace, replay and audit code out of self.f90 (and any copies in bands.f90 /
     hamiltonian_build.f90) into source/ldau_diagnostics.f90, keeping the runtime switches.
   - Delete the hooks of closed campaigns: RSLMTO_FIXED_U_A2_SCAN, RSLMTO_C2_HUBBARD_AUDIT and
     RSLMTO_C2_HUBBARD_AUDIT_FILE, plus any code reachable only through them.
   - Keep the TRACE family and hubbard_test_operator: RD4 still needs the W path.
   - Report the line counts of self.f90 and of the new module, before and after.
   GATE: the identity gate.

T2 Screen output
   - Inventory every script that parses log lines (tests/, tools/, docs/validation/**/*.py, the
     preflight wrapper). List the tokens each one reads.
   - Introduce hubbard_log_level = 'summary' (default) | 'detailed' | 'debug'.
   - summary:
       * one start-of-run block: form, U/J per type in eV, radial mode with source and reason,
         energy routing, provenance verdict;
       * one compact line per iteration per correlated site: n_up, n_dn, moment, E_Hub, P/V
         residuals.
   - detailed: the current per-site LDM / physicality / energy lines, printed only for channels
     that carry U/J.
   - debug: everything printed today.
   - Every token a parser reads stays available at summary level, or the parser is updated in the
     same commit.
   - Also remove repeated or redundant non-Hubbard lines that the LDA+U work added (for example
     identical messages printed every iteration).
   GATE: the identity gate; all log parsers pass.

T3 Density-matrix file
   - For atoms carrying U/J, write <label>_dm.nml next to <label>_out.nml.
   - Header: label, radial mode, backend, U/J in eV per channel, run id, iteration.
   - Data: square blocks dm_re / dm_im indexed by (l, spin_row, spin_col), one matrix row per line,
     spinor-ready. Collinear runs write only the spin-diagonal blocks.
   - The potential file gains dm_file and dm_sha256 and stops writing ldm_flatten /
     ldm_flatten_imag.
   - Read order: the referenced dm file (fatal if missing or the hash mismatches) -> legacy inline
     ldm_flatten(_imag) -> zero.
   - MPI: rank 0 reads and writes, then broadcasts as today.
   - Tests:
       * a round trip reproduces the LDM bitwise;
       * a legacy inline file is read correctly;
       * a missing dm file is fatal; a hash mismatch is fatal;
       * the provenance guard behaves unchanged;
       * the impurity workflow works on a host with a dm file: the impurity run reads the host
         occupation matrix from the dm file and reproduces the pre-change impurity outputs.
   GATE: the identity gate, with +U potential files compared by value: every non-LDM field byte-
   identical, and the dm file values equal the old inline values bitwise.

T4 Small corrections from the audit
   (a) fixed_hubbard_energy_pout: read P_out explicitly (or rename it if potential%ldm is provably
       P_out at the call site). The converged bulk +U energies must be unchanged.
   (b) Require hubbard_experimental = .true. to select 'acbn0' or any +V input; otherwise stop
       with an actionable message.
   (c) lda_u_workflow.rst:
       - update the validated-scope section: KS and RS bulk with energy; RS impurity density
         only; surfaces experimental;
       - correct the Dudarev wording (validated limit: J=0 with U);
       - document the dm file and hubbard_log_level.
   (d) KNOWN_ISSUES:
       - mark surfaces experimental, noting the Ni(001) window/NaN observation;
       - add an OPEN entry: Ni gating host was non-magnetic (recheck in S0b).
   GATE: the identity gate.

T5 Local --no-ff merge into fable_v4b_u; rerun the identity gate on the merged tree. Report the
   SHA and the line-count summary. Do not push.

Uncovered outcomes: record, leave production unchanged, continue (C5.3 item 2). Tooling failures
are fixed, not reported as gate results (C6.3).
