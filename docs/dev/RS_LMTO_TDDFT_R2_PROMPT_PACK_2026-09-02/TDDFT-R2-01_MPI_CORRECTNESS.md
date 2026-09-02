# TDDFT-R2-01 — Native real-space GF MPI correctness

Work on the current pushed `fable_v4` branch of RS-LMTO-ASA.

This is **TDDFT-R2-01: prove and, if necessary, repair MPI correctness of the native real-space Green-function TD-DFT backend**.

Do not add new TD-DFT physics. Do not change the mathematical definition of the response unless an actual correctness defect is demonstrated.

## Background

The native real-space TD-DFT path is intended to implement

\[
G_{ab}(R,z),G_{ba}(-R,z)
\rightarrow
\chi^0_{ab}(R,\omega)
\rightarrow
\chi^0_{ab}(q,\omega),
\]

without constructing \(G(k,z)\).

The previous audit found that `tddft_chi0_realspace.f90` appears to consume MPI-local real-space GF data and mentions a reduction of rank-local contributions, but no obvious reduction was found inside that module. This must be resolved before multi-rank calculations can be trusted.

## Tasks

1. Trace the complete MPI ownership/reduction path for:
   - `green%gij`;
   - `green%gji`;
   - real-space pair ownership;
   - `chi_r`;
   - Fourier-transformed `chi_q`.

2. Determine conclusively whether every MPI rank currently returns:
   - the complete susceptibility;
   - only its local partial contribution;
   - or replicated complete data through an outer reduction.

3. If a reduction is missing, implement it at the mathematically correct level.

   Prefer reducing the smallest sensible object. Do not blindly replicate very large intermediates if a smaller reduction is possible.

4. Ensure that normalization is performed exactly once. Be especially alert for:
   - double division by number of ranks;
   - missing lattice multiplicities;
   - reduction before versus after Fourier phase factors;
   - rank-local omissions of \(R\) pairs.

5. Add deterministic MPI regression tests comparing the native-RGF result from

\[
N_{\rm MPI}=1,\;2,\;\text{and if practical }4.
\]

Compare the **complex \(\chi^0(q,\omega)\)** itself, not just a derived peak energy.

6. Test at:
   - \(q=0\);
   - one finite generic \(q\);
   - at least two frequencies.

7. Use a small deterministic fixture for CI, but also perform one check using the existing `examples/susceptibility` bcc Fe or fcc Ni input if affordable.

## MPI ownership and reduction path (verified)

- `get_mpi_variables(rank,njij)` partitions `lattice%ijpair` into contiguous global pair ranges `[start_atom,end_atom]`. The recursion and Green-function builders use `g2l_map`, so each rank fills only its MPI-local `green%gij/gji` pair slots.
- `initialize_from_green` maps those local slots back to their global pair geometry and passes only the owned pairs to the native real-space provider. The provider accumulates local `chi_r(R,omega)` and applies the Fourier phase while forming a local partial `chi_q(q,omega)`.
- `reduce_realspace_chi0_batch` performs one `MPI_ALLREDUCE(SUM)` on the complex local `chi_q` contributions, plus the additive derived spectra and diagnostics. Cutoff/shell summary fields use conservative `MAX` reductions. The production loop then writes only the q-owner files.
- The phase is applied before communication, each real-space pair is owned by exactly one rank, and no MPI-size division is applied. There is no second reduction of the native χ⁰ result.

## Non-negotiables

- Do not disable MPI for the real-space backend.
- Do not make every rank independently recompute the entire response merely to make the test pass.
- Do not loosen tolerances until an actual floating-point reduction-order difference has been quantified.
- Serial and MPI must represent the same mathematical susceptibility.

## Acceptance checklist

Tick these in the task notes as they are completed:

- [x] MPI ownership and reduction path documented.
- [x] Missing/double reduction ruled out or repaired.
- [x] Serial and 2-rank complex \(\chi^0\) agree.
- [x] 4-rank comparison performed if supported by the test environment.
- [x] \(q=0\) tested.
- [x] finite \(q\) tested.
- [x] multiple frequencies tested.
- [x] no normalization-by-rank bug.
- [x] regression test added.
- [x] existing serial TD-DFT tests remain green.

## Evidence

The deterministic CI fixture uses eight disjoint real-space pairs, `q=[0,0,0]` and
`q=[0.173,-0.219,0.071]`, and `omega=[0.07,0.23]` Ry. Its nonzero reference
response scale is `5.3294165139979253E+00`. The maximum absolute and relative
differences against the rank-one reference were:

| MPI ranks | max absolute difference | max relative difference |
|---:|---:|---:|
| 1 | `0.0000000000000000E+00` | `0.0000000000000000E+00` |
| 2 | `0.0000000000000000E+00` | `0.0000000000000000E+00` |
| 4 | `0.0000000000000000E+00` | `0.0000000000000000E+00` |

The pre-reduction local partial-response maximum was
`5.3294165139979253E+00` for both the 2- and 4-rank runs, confirming that the
MPI test uses disjoint local work rather than recomputing the complete response
on every rank. The test tolerances are `2E-11` absolute and `2E-12` relative.

The available bcc-Fe native-RGF material deck was also attempted in isolated
temporary directories for 1 and 2 MPI ranks. Both runs stopped during input
setup because this checkout refers to the unavailable external element database
`/home/lucas/projetos-novos/pesquisa/rslmto/rslmto/database/elements/Fe.nml`;
therefore that optional material check could not reach TD-DFT and was not used
as acceptance evidence.

## Commit

`Fix MPI reduction for real-space TDDFT response`
