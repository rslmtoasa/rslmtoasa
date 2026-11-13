#!/usr/bin/env python3
"""LdaUKspaceEnergy: the Hubbard energy of hubbard_radial_moments = 'local' in k-space (bcc Fe, first order, U = 2 eV, J = 0, single thread). The routing
adds the energy of the spinor P_out to physical_total_energy; for collinear setups the spin representations describe the same state, so E_Hub_P_OUT
(Ry, 16 digits) and physical_total_energy (Ry) must agree at every step (the decks of test_kspace_local_nsp3.py and test_local_nsp4_workflows_u.py):
  af_nsp3_nsp1    two-atom antiferromagnet, 3 steps: nsp = 3 (same file, mom +z, -z) against nsp = 1 (second site's file with swapped channels)
  af_nsp4_nsp3    the same antiferromagnet, nsp = 4 with SOC held at zero against nsp = 3
  fm_nsp4_nsp2    ferromagnet along z with SOC, 4 steps: nsp = 4 against nsp = 2
  engaged         |E_Hub_P_OUT| of the nsp = 1 antiferromagnet and of the nsp = 2 ferromagnet is above 1e-2 Ry at every step (a skipped routing prints 0)
Bug caught: the k-space energy not routed at nsp >= 2, or built from the wrong spin frame or site. Not covered: hoh (test_hoh_energy_u.py), impurities and surfaces
(test_cluster_energy_u.py), an independent reference for the energy itself (the nsp = 1 energy is the reference), a finite-difference dE/dU check.

usage: test_kspace_energy_u.py --binary rslmto.x --scratch-root DIR
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
import test_local_energy_u as le  # noqa: E402
import test_local_nsp4_workflows_u as wf  # noqa: E402
import test_nsp34_u as t3  # noqa: E402


def energies(log: str, steps: int):
    log = re.sub(r"\x1b\[[0-9;]*m", "", log)
    hub = [float(x) for x in re.findall(r"E_Hub_P_OUT=\s*(\S+)", log)]
    tot = [float(x) for x in re.findall(r"(?<!base_)physical_total_energy=\s*(\S+)", log)]
    return (hub, tot) if len(hub) == len(tot) == steps else f"expected {steps} energies, got {len(hub)}, {len(tot)}"


def af(binary: Path, base: Path, name: str, nsp: int):
    work = t3.af_deck(base, name, nsp=nsp, mom2=(0.0, 0.0, -1.0), swap2=(nsp == 1), nstep=wf.NSTEP, kspace=True, radial="local", hoh=False,
                      self_extra=wf.SOC_OFF if nsp == 4 else "")
    status, log = gh.run_binary(binary, work)
    return f"{name}: run failed\n{log[-1200:]}" if status else energies(log, wf.NSTEP)


def fm(binary: Path, base: Path, name: str, nsp: int):
    d = wf.ks_fm(binary, base, name, nsp)
    return d if isinstance(d, str) else energies(rr.log_text(d), rr.NSTEP)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    r = {"a1": af(binary, base, "af_nsp1", 1), "a3": af(binary, base, "af_nsp3", 3), "a4": af(binary, base, "af_nsp4", 4),
         "f2": fm(binary, base, "fm_nsp2", 2), "f4": fm(binary, base, "fm_nsp4", 4)}
    cases = [("af_nsp3_nsp1", le.compare(r["a3"], r["a1"], "k-space AF nsp=3 against nsp=1")),
             ("af_nsp4_nsp3", le.compare(r["a4"], r["a3"], "k-space AF nsp=4 SOC off against nsp=3")),
             ("fm_nsp4_nsp2", le.compare(r["f4"], r["f2"], "k-space FM with SOC nsp=4 against nsp=2"))]
    bad = next((v for v in (r["a1"], r["f2"]) if isinstance(v, str)), None)
    if bad:
        cases.append(("engaged", (False, bad)))
    else:
        low = min(abs(x) for k in ("a1", "f2") for x in r[k][0])
        cases.append(("engaged", (low > le.ENGAGED, f"smallest |E_Hub_P_OUT| of the nsp=1 AF and nsp=2 FM runs {low:.1e} Ry (need > {le.ENGAGED:.0e})")))
    failed = False
    for name, (ok, info) in cases:
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
