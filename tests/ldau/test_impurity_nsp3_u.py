#!/usr/bin/env python3
"""LdaUImpurityNsp3: hubbard_radial_moments = 'local' at nsp = 3 on an impurity cluster (calctype I; B2FeCo impurity deck, block recursion, 4 steps,
single thread, first order and hoh) embedded in a k-space Fe host written at the same nsp and hoh setting (U = 2 eV, J = 0), moments along +z. For a
collinear moment the nsp = 3 spinor correction (along the site axis, V by type) must reproduce the nsp = 1 scalar correction (with hoh: V_in captured at
the enim insertion): the per-iteration Diff, Fermi energy and spin moment agree (tolerance of the first-order impurity soc0 rows). Bug caught: the
correction missing, or taken in the wrong spin frame or with the wrong V, at nsp = 3 on a cluster.
  nsp3_vs_nsp1, nsp3_vs_nsp1_hoh   the parity, without and with hoh
  engaged, engaged_hoh             the same impurity in a spectral-mode host differs from the local run by more than 1e-2 muB in the spin moment at some step.
Not covered: surfaces, nsp = 4, noncollinear axes, an independent reference for the corrected moments themselves.

usage: test_impurity_nsp3_u.py --binary rslmto.x --scratch-root DIR
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
ENGAGED = 1.0e-2  # muB


def host_and_impurity(binary: Path, base: Path, tag: str, nsp: int, mode: str, hoh: bool = False):
    keys = f"hubbard_radial_moments = '{mode}'\n" + ("hoh = .true.\n" if hoh else "")
    host_edit = lambda t: rr.add_u(rr.use_kspace(rr.set_nstep(t)), 2.0, 0.0, keys).replace("nsp = 1", f"nsp = {nsp}")
    host = rr.prepare(rr.CASES / "k_space_scf/bccFe", base / f"{tag}_host", host_edit)
    if rr.run(binary, host) != 0:
        return f"{tag}: host run failed\n{rr.log_text(host)[-1200:]}"
    imp_edit = lambda t: rr.in_hamiltonian(rr.set_nstep(t).replace("nsp = 1", f"nsp = {nsp}").replace("hoh = false", "hoh = .true." if hoh else "hoh = false"),
                                           f"hubbard_radial_moments = '{mode}'\n")
    d = rr.prepare(rr.CASES / "impurity/B2FeCo_nsp1", base / tag, imp_edit)
    shutil.copy2(host / "Fe_out.nml", d / "Fe.nml")
    shutil.copy2(host / "Fe_dm.nml", d / "Fe_dm.nml")
    if rr.run(binary, d) != 0:
        return f"{tag}: impurity run failed\n{rr.log_text(d)[-1200:]}"
    return d


def signature(d: Path) -> dict:
    sig = rr.signature(d, True)
    return {k: v for k, v in sig.items() if not k.startswith("dm_") and k != "energy"}


def nsp3_vs_nsp1(binary: Path, base: Path, hoh: bool):
    a = host_and_impurity(binary, base, "imp_nsp1", 1, "local", hoh)
    b = host_and_impurity(binary, base, "imp_nsp3", 3, "local", hoh)
    for x in (a, b):
        if isinstance(x, str):
            return False, x
    return rr.compare_signatures(signature(b), signature(a), TOL)


def engaged(binary: Path, base: Path, hoh: bool):
    loc = host_and_impurity(binary, base, "imp_nsp3", 3, "local", hoh)
    spec = host_and_impurity(binary, base, "imp_spec", 3, "spectral", hoh)
    for x in (loc, spec):
        if isinstance(x, str):
            return False, x
    m = lambda d: [float(v) for v in re.findall(r"Spin moment of atom\s+1 is\s+(\S+)", rr.log_text(d))]
    ml, ms = m(loc), m(spec)
    if not ml or len(ml) != len(ms):
        return False, "no spin moments found, or different numbers of steps"
    d = max(abs(a - b) for a, b in zip(ml, ms))
    return d > ENGAGED, f"nsp = 3 local against spectral: max |dm| over {len(ml)} steps {d:.1e} muB (need > {ENGAGED})"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    failed = False
    for name, result in (("nsp3_vs_nsp1", nsp3_vs_nsp1(binary, base / "parity", False)), ("engaged", engaged(binary, base / "engaged", False)),
                         ("nsp3_vs_nsp1_hoh", nsp3_vs_nsp1(binary, base / "parity_hoh", True)), ("engaged_hoh", engaged(binary, base / "engaged_hoh", True))):
        ok, info = result
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
