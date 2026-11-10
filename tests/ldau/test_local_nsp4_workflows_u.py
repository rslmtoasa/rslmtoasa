#!/usr/bin/env python3
"""LdaULocalNsp4Workflows: hubbard_radial_moments = 'local' at nsp = 4 in k-space and on an impurity cluster (U = 2 eV, J = 0, single thread). At nsp = 4 the
spinor correction is the nsp = 3 one with SOC added to the Hamiltonian, so against nsp = 3 with SOC held at zero (soc_scale = 0, fix_soc) and against
nsp = 2 for a ferromagnet along z the runs must agree:
  ks_af_soc0      k-space two-atom Fe antiferromagnet (mom = +z, -z), 3 steps: nsp = 4 SOC off against nsp = 3 (E_F, both site moments)
  ks_af_soc0_hoh  the same with hoh (V_U in enim)
  ks_fm_nsp2      k-space bcc Fe ferromagnet with SOC, 3 steps: nsp = 4 against nsp = 2 (Diff, E_F, spin moment)
  imp_soc0        B2FeCo impurity (calctype I, 4 steps, block recursion) in an nsp = 4 k-space host, SOC off: nsp = 4 against nsp = 3
  imp_soc0_hoh    the same with hoh
  engaged         k-space AF at nsp = 4, local against spectral: E_F differs by more than 1e-3 Ry (a skipped correction would make the parities hold trivially)
Bug caught: the nsp = 4 k-space or cluster route missing the correction or applying it in the wrong spin frame, or SOC lost from the moments. Not covered:
surfaces (the Ni(001) deck differs between nsp = 3 and nsp = 1 from step 2 on, see KNOWN_ISSUES.md), SOC-on impurity clusters, noncollinear axes at nsp = 4
in k-space, an independent reference for the corrected moments.

usage: test_local_nsp4_workflows_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guard_helpers as gh  # noqa: E402
import reference_rows as rr  # noqa: E402
import test_impurity_nsp3_u as imp  # noqa: E402
import test_nsp34_u as t3  # noqa: E402

TOL = 2.0e-6  # Ry and muB for the k-space AF (six printed decimals)
SOC_OFF = "soc_scale = 0.0\nfix_soc = .true.\n"
NSTEP = 3


def ks_af(binary: Path, base: Path, name: str, nsp: int, hoh: bool, radial: str = "local"):
    work = t3.af_deck(base, name, nsp=nsp, mom2=(0.0, 0.0, -1.0), nstep=NSTEP, kspace=True, radial=radial, hoh=hoh,
                      self_extra=SOC_OFF if nsp == 4 else "")
    status, log = gh.run_binary(binary, work)
    log = re.sub(r"\x1b\[[0-9;]*m", "", log)
    if status != 0:
        return f"{name}: run failed\n{log[-1200:]}"
    ef = [float(x) for x in re.findall(r"N\(E_F=\s*(-?\d+\.\d+)", log)]
    m = [[float(x) for x in re.findall(rf"Spin moment of atom\s+{a} is\s+(\S+)", log)] for a in (1, 2)]
    return ef + m[0] + m[1] if len(ef) == len(m[0]) == len(m[1]) == NSTEP else f"{name}: expected {NSTEP} steps, got {len(ef)}, {len(m[0])}, {len(m[1])}"


def ks_af_case(binary: Path, base: Path, tag: str, hoh: bool):
    a, b = ks_af(binary, base, f"{tag}_nsp3", 3, hoh), ks_af(binary, base, f"{tag}_nsp4", 4, hoh)
    for v in (a, b):
        if isinstance(v, str):
            return False, v
    d = max(abs(x - y) for x, y in zip(a, b))
    return d < TOL, f"nsp = 4 SOC off against nsp = 3, over {NSTEP} steps: max |difference| of E_F (Ry) and moments (muB) {d:.1e} (need < {TOL:.0e})"


def engaged(binary: Path, base: Path):
    a, b = ks_af(binary, base, "eng_local", 4, False), ks_af(binary, base, "eng_spectral", 4, False, "spectral")
    for v in (a, b):
        if isinstance(v, str):
            return False, v
    d = max(abs(x - y) for x, y in zip(a[:NSTEP], b[:NSTEP]))
    return d > 1.0e-3, f"k-space AF at nsp = 4 (SOC off), local against spectral: max |dE_F| {d:.1e} Ry (need > 1e-3)"


def ks_fm(binary: Path, base: Path, tag: str, nsp: int):
    edit = lambda t: rr.add_u(rr.ks_nsp2_deck(t), 2.0, 0.0, "hubbard_radial_moments = 'local'\n").replace("nsp = 2", f"nsp = {nsp}", 1)
    d = rr.prepare(rr.CASES / "k_space_scf/bccFe", base / tag, edit)
    return d if rr.run(binary, d) == 0 else f"{tag}: run failed\n{rr.log_text(d)[-1200:]}"


def ks_fm_case(binary: Path, base: Path):
    a, b = ks_fm(binary, base, "ksfm_nsp2", 2), ks_fm(binary, base, "ksfm_nsp4", 4)
    for v in (a, b):
        if isinstance(v, str):
            return False, v
    return rr.compare_signatures(imp.signature(b), imp.signature(a), 5 * rr.IMP_SOC0_TOL)


def imp_case(binary: Path, base: Path, tag: str, hoh: bool):
    a = imp.host_and_impurity(binary, base, f"{tag}_nsp3", 3, "local", hoh)
    b = imp.host_and_impurity(binary, base, f"{tag}_nsp4", 4, "local", hoh, soc0=True)
    for v in (a, b):
        if isinstance(v, str):
            return False, v
    return rr.compare_signatures(imp.signature(b), imp.signature(a), rr.IMP_SOC0_TOL)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    failed = False
    for name, result in (("ks_af_soc0", ks_af_case(binary, base, "ksaf", False)), ("ks_af_soc0_hoh", ks_af_case(binary, base, "ksafh", True)),
                         ("ks_fm_nsp2", ks_fm_case(binary, base)), ("engaged", engaged(binary, base)),
                         ("imp_soc0", imp_case(binary, base / "imp", "imp", False)), ("imp_soc0_hoh", imp_case(binary, base / "imph", "imph", True))):
        ok, info = result
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
