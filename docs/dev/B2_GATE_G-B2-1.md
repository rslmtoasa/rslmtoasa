# Gate G-B2-1 — energy-contour / grid convention for `reciprocal_green`

**Status:** OPEN — awaiting maintainer (Anders) sign-off before task B2.2.
**Task:** B2.1 (module skeleton + `sigma_provider` + contour adoption).
**Branch:** `fable_v2`.

## What you are signing

That the energy grid + retarded contour the two k-space backends will use
matches what the recursion route produces and what `exchange.f90` (and the
damping route) integrate over — so the canonical `green` arrays can be filled
by the k-space engine with **zero consumer changes**.

## The adopted convention (pinned in code)

Implemented in `source/reciprocal_green.f90::build_green_contour` and
documented at the top of that file.

| Item | Pin | Evidence in RS route |
|---|---|---|
| **Grid** | `en%ene(1:size(en%ene))`, real-axis mesh from `energy%e_mesh` (`energy_min .. above fermi`, `en%channels_ldos + 10` points) | `bgreen` sets `e = this%en%ene`; `exchange` integrates `simpson_f(..., this%en%ene, ...)` (`exchange.f90:167,315`) |
| **Contour** | retarded, `z(ie) = en%ene(ie) + i*green_eta` | `bgreen` forms `ze = e(ei)`, `Q = (ze + eta - ...)` with `eta = i*eta_imag` (`green.f90:1555–1587`) |
| **Fermi** | physical `en%fermi` (integration upper limit); **never** the `chebfermi`-scaled variable | `simpson_f(..., this%en%fermi, this%en%nv1, ...)` (`exchange.f90:167`) |
| **Representation** | screened/auxiliary LMTO blocks, **pre-`auxiliary_gij`** | consumers call `green%auxiliary_gij` themselves after the arrays are filled |
| **Broadening default** | `green_eta = 0.01_rp` Ry (matches the default DOS smearing scale) | set in `reciprocal_lifecycle.f90::restore_to_default`; namelist wiring is B2.5 |

## Open questions for the maintainer

1. **Default `green_eta` value.** 0.01 Ry is a placeholder tied to the DOS
   smearing scale. The RS `gij` route's fixed broadening is passed into
   `bgreen` as `eta` by `block_green_ij`; if that value differs, the C1
   elementwise match in B2.2 will want the same number. Confirm the target
   default (or defer to the B2.5 `&reciprocal_green` namelist as user-set).
2. **Grid length.** The contour spans the full `en%ene` (`channels_ldos + 10`)
   to match the `green` array third dimension. The Simpson integral uses
   `en%nv1` points up to `en%fermi`; both backends fill the whole grid and let
   the consumer choose the sub-range. OK?
3. **On-site representation subtlety (flag for B2.2, not this gate).** The RS
   route stores the on-site `g0` as `-i*pi*rho` (anti-Hermitian / spectral
   only; `sgreen`, `green.f90:892`), while the intersite `gij` from `bgreen`
   is the full complex resolvent. The C1 test ("on-site block vs RS at large
   broadening") must compare against the correct one of these. Noted here so
   B2.2 picks the right reference; not part of the contour pin.

## Scope delivered in B2.1

- `source/sigma_provider.f90` — abstract `sigma_provider` + `sigma_zero`.
- `source/reciprocal_green.f90` — `build_green_contour` (concrete, above) and
  `fill_green` dispatcher (backend cores are staged: E→B2.2, D→B2.4, currently
  raise not-implemented).
- `reciprocal` type gains `green_eta`, `green_backend`; defaults set.
- CMake registers both new files.
- Regression matrix **10/10, bit-identical** (`fill_green` is not yet wired
  into any production path — dispatch key is B2.5).

## Note (out of scope, pre-existing)

Configuring with `RUN_EXAMPLE_TESTS=ON` currently fails: three frozen-magnon
cases in `tests/scf/cases.json` (indices 18–20, `Example_frozen_magnon_bccFe*`)
lack a `namelists.control` member that `register_example_tests` requires
(`CMakeLists.txt:389`). This is B1 test-data debt, independent of B2. Built
here with `-DRUN_EXAMPLE_TESTS=OFF`; flag restored afterward.
