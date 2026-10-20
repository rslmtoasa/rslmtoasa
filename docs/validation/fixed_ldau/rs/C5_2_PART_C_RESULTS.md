# C5.2 Part C — preliminary offline audit; comparison gates not run

| Item | Status | Record |
|---|---|---|
| C1 fixed-U Chebyshev code audit | **PRELIMINARY_OFFLINE_AUDIT_ONLY** | Source traversal indicates deterministic full-matrix g0/P and local-mode M1; no stochastic trace is used on this SCF path. |
| C2 block/Chebyshev comparisons | **NOT_RUN** | U=0 and U=2 eV comparison pairs were not launched. |
| C3 carrier check | **NOT_RUN** | No local-versus-spectral replay was launched. |
| Co U>0 admission | **NOT_ADMITTED** | No C2/C3 compatibility gates passed. |
| Part C decision | **STOPPED** | Part A A4(iv) STOP precedes these gates. |

**Interpretation (5 lines):**

- C1 is a preliminary source audit only, not runtime validation or a Part C gate result.
- The Chebyshev recursion seeds a deterministic identity orbital block and reconstructs full matrix Green functions.
- Fixed Liechtenstein P uses full orbital-pair spectral integrals; local-mode M1 uses full pairwise energy-weighted integrals.
- The separate stochastic Chebyshev routine is for transport; its random-vector count does not apply to fixed-U P/M1.
- C2/C3 were not run, no compatibility gate passed, and Co U>0 is not admitted.

Passed: no Part C gate was evaluated.
Recorded: preliminary C1 source audit only.
Provisional: C1 findings remain unverified at runtime.
Stopped: Part A A4(iv); C2/C3 NOT_RUN.
Next: user review of the no-reproduction branch.
