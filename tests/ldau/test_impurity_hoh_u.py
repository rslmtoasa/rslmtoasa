#!/usr/bin/env python3
"""LdaUImpurityHoh: hubbard_radial_moments = 'local' with hoh on an impurity cluster (calctype I; B2FeCo impurity deck, block recursion, 4 steps,
single thread) embedded in a k-space Fe host written in the same mode (hoh, U = 2 eV, J = 0). With hoh V_U is an additive onsite term in enim,
and the moments of H0 = H - V_U need the V_in that was actually inserted; the cases:
  nsp2_vs_nsp1  nsp = 2 with SOC held at zero (spinor correction; V_in read by type) against nsp = 1 (scalar correction; V_in captured at the
                enim insertion): the per-iteration Diff, Fermi energy, spin moment and spin-diagonal d density-matrix entries must agree
                (tolerance of the first-order impurity soc0 rows, 1e-7). Bug caught: a V_in captured at the wrong site, type or spin, or not
                captured (the run then stops: 'RD1 exact V-in snapshot missing/stale').
  engaged       the same impurity in a spectral-mode host (hoh) gives spin moments that differ from the local run by more than 1e-2 muB at
                some step, so the correction is acting (a skipped correction would make the first case pass trivially).
Not covered: surface (calctype S) with hoh, nsp = 3 and 4 on impurities, an independent reference for the corrected moments themselves
(an Fe impurity in an Fe host does not reproduce the bulk, even at U = 0 and without hoh: KNOWN_ISSUES.md, embedding offset).

usage: test_impurity_hoh_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference_rows as rr  # noqa: E402

TOL = rr.IMP_SOC0_TOL
HOH = "hoh = .true.\n"
ENGAGED = 1.0e-2  # muB


def host_and_impurity(binary: Path, base: Path, tag: str, nsp: int, mode: str):
    """Fe host (k-space, hoh, U = 2 J = 0, radial moments `mode`), then the B2FeCo impurity (hoh, same mode); SOC held at zero at nsp = 2."""
    keys = f"hubbard_radial_moments = '{mode}'\n" + HOH
    if nsp == 2:
        host_edit = lambda t: rr.soc_off(rr.add_u(rr.ks_nsp2_deck(t), 2.0, 0.0, keys))
    else:
        host_edit = lambda t: rr.add_u(rr.use_kspace(rr.set_nstep(t)), 2.0, 0.0, keys)
    host = rr.prepare(rr.CASES / "k_space_scf/bccFe", base / f"{tag}_host", host_edit)
    if rr.run(binary, host) != 0:
        return f"{tag}: host run failed\n{rr.log_text(host)[-1200:]}"

    def imp_edit(t):
        t = rr.in_hamiltonian(rr.set_nstep(t).replace("nsp = 1", f"nsp = {nsp}").replace("hoh = false", "hoh = .true."),
                              f"hubbard_radial_moments = '{mode}'\n")
        return rr.soc_off(t) if nsp == 2 else t

    d = rr.prepare(rr.CASES / "impurity/B2FeCo_nsp1", base / tag, imp_edit)
    shutil.copy2(host / "Fe_out.nml", d / "Fe.nml")
    shutil.copy2(host / "Fe_dm.nml", d / "Fe_dm.nml")
    if rr.run(binary, d) != 0:
        return f"{tag}: impurity run failed\n{rr.log_text(d)[-1200:]}"
    return d


def nsp2_vs_nsp1(binary: Path, base: Path) -> tuple[bool, str]:
    a = host_and_impurity(binary, base, "imp_nsp1", 1, "local")
    b = host_and_impurity(binary, base, "imp_nsp2", 2, "local")
    for x in (a, b):
        if isinstance(x, str):
            return False, x
    return rr.compare_signatures(rr.signature(b, True), rr.signature(a, True), TOL)


def engaged(binary: Path, base: Path) -> tuple[bool, str]:
    loc = host_and_impurity(binary, base, "imp_nsp1", 1, "local")
    spec = host_and_impurity(binary, base, "imp_spec", 1, "spectral")
    for x in (loc, spec):
        if isinstance(x, str):
            return False, x
    m = lambda d: [float(v) for v in re.findall(r"Spin moment of atom\s+1 is\s+(\S+)", rr.log_text(d))]
    ml, ms = m(loc), m(spec)
    if not ml or len(ml) != len(ms):
        return False, "no spin moments found, or different numbers of steps"
    d = max(abs(a - b) for a, b in zip(ml, ms))
    return d > ENGAGED, f"local against spectral: max |dm| over {len(ml)} steps {d:.1e} muB (need > {ENGAGED})"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    failed = False
    for name, result in (("nsp2_vs_nsp1", nsp2_vs_nsp1(binary, base / "parity")), ("engaged", engaged(binary, base / "engaged"))):
        ok, info = result
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
