# TDDFT-R2-06 — Adversarial closeout and production boundary

Work on `fable_v4` after TDDFT-R2-01 through R2-05.

This is **TDDFT-R2-06: scientific closeout of the transverse TD-DFT baseline**.

Do not implement new response physics.

## Goal

Decide exactly what RS-LMTO-ASA can now claim to support.

## Tasks

1. Re-audit all three backends:

\[
\text{eigenpairs},\qquad
G(k,z),\qquad
G(R,z).
\]

2. Review every R2 finding and verify it against tests and stored material evidence.

3. Run the complete TD-DFT test suite.

4. Run the bcc Fe and fcc Ni production validation fixtures.

5. Check documentation against actual code behavior.

6. Produce a compact support matrix, e.g.

| Physics | Supported | Validated | Notes |
|---|---|---|---|
| collinear FM, SOC=0 transverse | | | |
| eigenpair backend | | | |
| K-GF backend | | | |
| R-GF backend | | | |
| static Ward | | | |
| MPI R-GF | | | |
| AF/ferri | | | |
| longitudinal | | | |
| SOC | | | |
| noncollinear | | | |

7. Distinguish clearly:
   - **implemented**;
   - **unit-tested**;
   - **material-validated**;
   - **production-supported**.

8. Keep unsupported configurations guarded rather than silently falling through.

9. Identify remaining items for a future phase:
   - AF/ferri material validation;
   - longitudinal FeSe/second-sound validation;
   - Sternheimer/Savrasov;
   - SOC/noncollinear response.

Do not start them here.

## Acceptance checklist

- [x] all R2 tests green.
- [x] Fe/Ni material evidence reviewed.
- [x] MPI status explicitly stated.
- [x] RGF contour status explicitly stated.
- [x] raw Goldstone status explicitly stated.
- [x] documentation contains no unsupported production claims.
- [x] unsupported physics still guarded.
- [x] production boundary is explicit.
- [x] remaining scientific risks documented.

## Deliverable

Create/update:

`docs/TDDFT_STATUS.md`

with an executive summary and evidence links.

## Completion evidence

The complete serial TDDFT suite passed 27/27 tests. The native real-space MPI
comparison passed at 1, 2, and 4 ranks. The deterministic campaign validator
passed, and all six bcc Fe/fcc Ni material decks were rerun in isolated
temporary copies. Their response, Goldstone, and manifest outputs match the
stored evidence apart from expected profiling timing metadata.

The material release decision remains bcc Fe **PARTIAL** and fcc Ni **FAIL**.
The three bare-response material matrices fail the required equivalence gate,
raw Ward convergence is not established, and no Goldstone correction was used
to hide those results. The production boundary and remaining risks are
recorded in [`docs/TDDFT_STATUS.md`](../../TDDFT_STATUS.md).

## Commit

`Close out baseline TDDFT validation`
