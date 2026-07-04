---

## Review notes after T1–T6 (independent verification, 2026-07-03)

Status: T1–T6 verified by an independent clean build (gfortran 13.3, OpenBLAS,
no MKL/CUDA) and by running the new regression matrix — **8/8 runnable cases
PASS** against the pre-refactor references (`bccFe_chebyshev_{fast, batched,
legacy, fast_hoh, legacy_hoh, fast_ccor_2c}`, `bccFe_block_fast_{sp, dp}`).
The two MKL cases are correctly skipped when built without `ENABLE_MKL_*`.

Three follow-ups before starting T7/T8, in priority order:

### R1. Use the abstract interface for the moment-driver callback (one-line fix)
`chebyshev_fast.f90 :: cheb_three_term_moments` declares its callback as

```fortran
external :: apply_h
```

but the module already defines `abstract interface :: ham_apply_t` for exactly
this purpose. Change to

```fortran
procedure(ham_apply_t) :: apply_h
```

Rationale: `external` gives the dummy an implicit interface, discarding all
compile-time argument checking — the entire point of defining `ham_apply_t` —
and is fragile when the actual arguments are internal procedures. No behavior
change; regression suite must stay green. Fold into the next commit that
touches the file.

### R2. Executable-stack trampolines — accept and document (no code change now)
The linker warns:

```
chebyshev_fast.f90.o: requires executable stack
```

Cause: the `apply_step*` callbacks are internal procedures capturing host
state (`ld`, `nb`, operator pointers), so gfortran emits stack trampolines for
the indirect calls. This is functionally correct, but on hardened systems that
enforce non-executable stacks (SELinux `allow_execstack=off`, some cluster
images) the binary can SIGSEGV at the first callback invocation.

Decision (KISS): keep the internal-procedure design for now — it is the
simplest correct structure. Add a note to the build docs / README that
`chebyshev_fast` requires an executable stack with gfortran. **Only if** a
target machine actually rejects the binary, escalate to the heavier fix:
promote the `apply_step*` implementations to module-level procedures taking a
small context derived type instead of host association. Do not do this
preemptively.

### R3. Remove the seven `ensure_*` delegation shims (KISS)
T5 left thin module-level wrappers (`ensure_work_buffers`,
`ensure_hoh_buffer`, `ensure_block_products`, `ensure_mkl_batch_ptr_arrays`,
`ensure_scaled_hamiltonian_sp`, `ensure_scaled_ortho_sp`,
`ensure_scaled_bsr_sp`, `ensure_scaled_velocity_sp`) that only delegate to the
type-bound `cheb_cache%ensure_*` procedures. All call sites live inside
`chebyshev_fast.f90`, so the shims add a parallel naming layer with no
external consumers. Delete the wrappers and call
`cheb_cache%ensure_...` directly at the ~10 call sites. Pure rename-level
change; regression suite must stay green. Fold into the next commit that
touches the file (natural companion to R1).

### Reviewer confirmations (no action needed)
- Moment recurrence in `cheb_three_term_moments` is verbatim-equivalent to the
  pre-refactor v1 loop (same `cherk`/`cgemm` sequence and doubling formulas).
- All four CPU backends route through the shared driver; stochastic and
  orbital kernels reuse the shared `hsweep_sp`/`happly_sp`.
- OpenMP `schedule` clauses preserved exactly across all apply variants.
- `bpopt`/`emami` retained per the terminator-estimator call chain
  (`get_terminf → get_cinf → bpopt → emami`).
- T6 `associate`-based operator selection is a valid (and clean) alternative
  to the pointer approach originally suggested.
