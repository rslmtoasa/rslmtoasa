#!/usr/bin/env python3
"""LdaULocalNsp4: hubbard_radial_moments = 'local' at nsp = 4 (noncollinear + SOC; real-space bulk, Chebyshev, U = 2 eV on Fe d, 2 steps, single thread).
The spinor correction of H0 = H - V_U acts on the radial C and Q2 parameters only; SOC enters H0 through xi and is not changed by it. The cases, each
against an independent route or a symmetry:
  soc0        two-atom Fe pair with moments along z and x, SOC held at zero (soc_scale = 0, fix_soc): nsp = 4 against nsp = 3 (E_F, |m1|, |m2|).
  nsp2        ferromagnet along z with SOC: nsp = 4 against nsp = 2 (the spinor route of the bulk).
  cubic       ferromagnet with SOC, moment along z against x (cubic equivalence): E_F and |m| agree.
  engaged     the (z, x) pair at nsp = 4 with SOC and radial moments 'spectral' differs from the local run by more than 1e-3 muB in |m| at some step
              (a skipped correction would make the other cases pass trivially).
Bug caught: a correction missing or applied in the wrong spin frame at nsp = 4, SOC lost from the moments.
Not covered: an independent reference for the corrected moments themselves, k-space, impurities and surfaces at nsp = 4, SOC-on noncollinear pairs.

usage: test_local_nsp4_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guard_helpers as gh  # noqa: E402
import test_hoh_u as hu  # noqa: E402
import test_nsp34_noncollinear as nc  # noqa: E402
import test_nsp34_u as t3  # noqa: E402

TOL = nc.TOL
ENGAGED = 1.0e-3  # muB
SOC_OFF = "soc_scale = 0.0\nfix_soc = .true.\n"


def run(binary: Path, work: Path):
    status, log = nc.run(binary, work)
    if status != 0:
        return f"{work.name}: run failed\n{log[-1500:]}"
    return t3.af_numbers(log)[0], nc.moments(log, 1), nc.moments(log, 2)


def fm(binary: Path, scratch: Path, name: str, nsp: int, mom=(0.0, 0.0, 1.0), radial: str = "local"):
    work = nc.fm_deck(scratch, name, mom, nsp=nsp, fixture=hu.SOC_FIXTURE)
    path = work / "input.nml"
    path.write_text(path.read_text().replace("&hamiltonian", f"&hamiltonian\nhubbard_radial_moments = '{radial}'\n", 1))
    return run(binary, work)


def worst(x, y) -> float:
    return max(abs(p - q) for a, b in zip(x, y) for p, q in zip(a, b))


def compare(runs: dict, a: str, b: str, tol: float, text: str):
    for v in (runs[a], runs[b]):
        if isinstance(v, str):
            return False, v
    d = worst(runs[a], runs[b])
    return d < tol, f"{text}: max |difference| of E_F (Ry) and |m| (muB) over 2 steps {d:.1e} (need < {tol:.0e})"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    pair = dict(radial="local", mom1=(0.0, 0.0, 1.0), mom2=(1.0, 0.0, 0.0))
    runs = {"p3": run(binary, t3.af_deck(base, "pair_nsp3", nsp=3, **pair)),
            "p4": run(binary, t3.af_deck(base, "pair_nsp4_soc0", nsp=4, self_extra=SOC_OFF, **pair)),
            "f2": fm(binary, base, "fm_nsp2", 2), "f4": fm(binary, base, "fm_nsp4_z", 4),
            "f4x": fm(binary, base, "fm_nsp4_x", 4, mom=(1.0, 0.0, 0.0)),
            "q4": run(binary, t3.af_deck(base, "pair_nsp4_soc", nsp=4, **pair)),
            "q4s": run(binary, t3.af_deck(base, "pair_nsp4_soc_spec", nsp=4, **{**pair, "radial": "spectral"}))}
    failed = False
    results = [("soc0", compare(runs, "p4", "p3", TOL, "nsp=4 SOC off against nsp=3, (z, x) pair")),
               ("nsp2", compare(runs, "f2", "f4", 5 * TOL, "nsp=4 against nsp=2, FM along z with SOC")),
               ("cubic", compare(runs, "f4", "f4x", 5 * TOL, "FM with SOC, z against x"))]
    if isinstance(runs["q4"], str) or isinstance(runs["q4s"], str):
        results.append(("engaged", (False, runs["q4"] if isinstance(runs["q4"], str) else runs["q4s"])))
    else:
        d = worst(runs["q4"][1:2], runs["q4s"][1:2])
        results.append(("engaged", (d > ENGAGED, f"local against spectral, max ||m| difference| over 2 steps {d:.1e} muB (need > {ENGAGED})")))
    for name, (ok, info) in results:
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
