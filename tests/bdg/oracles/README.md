# bdg oracles

Developer-owned reference data for the Stage 1 BdG kernels (`bdg_assemble`,
`bdg_accumulate`, `bdg_gap`; spec `cleanup/bdg_stage01_spec.md`, section 2).
Read-only to AI sessions after review (CLAUDE.md, "Verification integrity").

As in `tests/spin_response/oracles/`, the generated `.nml` files are not
committed. Run `python3 make_references.py OUTDIR` (about 1 s, numpy only); all
self-checks must print `ok` and the script exits non-zero otherwise. Commit the
generator, this README and `tolerances.nml`.

| File | Content |
|---|---|
| `make_references.py` | Generator for O1, O1b, O2, O3, O4, with built-in self-checks |
| `o1_single_band.nml` | O1: single band, Delta at kT = 0, 0.5 Tc, 0.9 Tc; Tc; normal state at 1.1 Tc |
| `o1b_uniform_k.nml` | O1b: three degenerate orbitals, K = (g/3) ones(3,3) |
| `o2_zeeman.nml` | O2: Zeeman field h = Delta0/2 at kT = 0 and 0.5 Tc |
| `o3_kodd.nml` | O3: k-odd, time-reversal-symmetric spin splitting; spectra at the Fermi points |
| `o4_two_band.nml` | O4: two decoupled bands, symmetric and asymmetric K |
| `tolerances.nml` | Every tolerance the BdG tests use (draft) |

## Where the numbers come from

None of the reference values use the spec's doubled 4-component Nambu
construction or any code under test:

- **Delta, chi, Tc** come from the scalar BCS gap equation (multi-channel for O1b and O4),
  chi = sum_k w Delta/(2E) * [1 - f(E+h) - f(E-h)], solved by Newton/bisection on
  the same discrete mesh.
- **Electron counts, n_up, n_dn** come from the reduced 2x2 Nambu basis
  (c_{k up}, c^dag_{-k dn}), diagonalized per k; this also re-derives chi.
- **O3 spectra** are closed form.

The doubled construction of spec section 1 appears only in the `spec_*`
self-checks, which confirm that the spec's formulas for H_BdG, chi and N
reproduce these references (agreement 1e-16). A `spec_*` failure means the spec
is wrong.

The same session drafted the spec and these references. The independence is
between formalisms (scalar and reduced-basis versus doubled Nambu), not between
authors, so please review the conventions in the generator's docstring:
chi = <c_up c_dn>, D = Delta (x) i sigma_y, and Delta = K chi with K > 0 attractive.

## Conventions the tests must follow

- **Units and mesh.** Energies are in Ry. The mesh is k_n = (n - 1/2)/nk, n = 1..nk, with w = 1/nk and nk = 2000. It is closed under k -> -k: H(-k_n) is the value at index nk + 1 - n.
- **Hamiltonian input.** H(k) is diagonal in the code's per-site index L + (sigma - 1)*norb. Build it from the band arrays (or the `band_rule` strings) and pass H(k) and H(-k) explicitly to `bdg_assemble`. The chemical potential is `mu`; the arrays do not include it.
- **Temperature.** kT = 0 means a step function. No mesh point sits exactly at E = 0.
- **Seeds and signs.** Seeds must be real and positive. Tests compare signed Delta and chi; comparing magnitudes would let a sign error pass.
- **Code path.** Tests drive the production kernels in the production iteration (assemble -> eigensolve -> accumulate -> gap -> mix). There is no special solver for tests.

## What each model anchors, and the bug it catches

| ID | Exact or reference statement | Catches |
|---|---|---|
| O1 | Delta0 = 5.018e-3 Ry, kTc = 2.850e-3 Ry, 2 Delta0/kTc = 3.521; Delta(T) at 0.5 and 0.9 Tc; N and N_normal; Delta -> 0 at 1.1 Tc | Wrong D sign or i sigma_y structure, missing 1 - f, wrong k weights, wrong N formula |
| O1b | Every Delta_L equals O1's kT = 0 value (exact) | K applied as a diagonal g/3 * I: then Delta = 0 |
| O2a | kT = 0, h = Delta0/2: Delta = Delta0, N = N0, n_up - n_dn = 0 (exact, the lowest quasiparticle energy is Delta0 - h > 0) | Field acting as a pair breaker at kT = 0, for example pairing placed inside one spin block |
| O2b | kT = 0.5 Tc, h = Delta0/2: Delta = 4.158e-3, n_up - n_dn = 5.74e-3 (single nonzero root on (0, 2 Delta0]) | Field dropped from the pair energies (gives the h = 0 value 4.801e-3) or occupations combined wrongly |
| O3 | Delta = 4.682e-3; spectra at the two Fermi points and one generic point | Hole block built from -H(k)* instead of -H(-k)*: the spectra differ by 4.3e-3 Ry and Delta collapses to 0 (`wrong_*` keys, informational) |
| O4a | Symmetric K: Delta_b = 2.25e-3, induced through g_ab | Cross-channel coupling dropped |
| O4b | Asymmetric K: Delta_b = 1.08e-3 | Transposed K (the transpose changes Delta by 1.9e-3) |

## Iteration budgets

The linear fixed-point rate near the solution (unmixed) is 0.765 at kT = 0,
0.804 at 0.5 Tc and 0.954 at 0.9 Tc; above Tc it is 0.978 at 1.1 Tc. With linear
mixing a the rate becomes 1 - a(1 - r). For the 0.9 Tc and 1.1 Tc cases a
kernel test can run with a = 1 to stay within the quick tier, or use
`delta_near_tc_rel` with fewer iterations.
