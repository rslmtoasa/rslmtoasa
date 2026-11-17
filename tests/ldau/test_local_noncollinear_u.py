#!/usr/bin/env python3
"""LdaULocalNoncollinear: hubbard_radial_moments = 'local' at nsp = 3 (bulk; real-space Chebyshev, and k-space for the rotation check) with site axes that are not collinear.
Two-atom Fe cell of test_nsp34_u.af_deck, U = 2 eV on both sites, 2 steps, single thread. The spinor correction of H0 = H - V_U projects
along each site axis, so a global spin rotation of the whole state, or the exchange of the two sites, must leave E_F and the site |m| unchanged
(no SOC: a global spin rotation is a symmetry). The expected values are those symmetries and the nsp = 1 ferromagnet; no run of the code is a reference.
  rotation   moments (z, x) against (x, -z) (rotation about y by 90 deg), (x, z) (sites exchanged) and a general rotation
             (60 deg polar, 45 deg azimuth) of the pair: E_F and |m1|, |m2| over 2 steps agree within TOL.
             Bug caught: a correction projected on the global z axis or on the wrong site's axis.
  kspace_rotation  the same rotation check in k-space (3 steps, |m1| and |m2| only: the log has no per-step E_F line there).
  parallel   two parallel moments at nsp = 3 against the nsp = 1 ferromagnet, both local: agree within TOL.
  engaged    the (z, x) pair with radial moments 'spectral' differs from the local run by more than 1e-3 muB in |m| at some step, so the
             correction acts (a skipped correction would make the first two cases pass trivially).
Not covered: SOC (nsp = 4), impurities and surfaces at nsp = 3, an independent reference for the corrected moments themselves.

usage: test_local_noncollinear_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import math
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_nsp34_noncollinear as nc  # noqa: E402
import test_nsp34_u as t3  # noqa: E402

TOL = nc.TOL
ENGAGED = 1.0e-3  # muB
TH, PH = math.radians(60.0), math.radians(45.0)
# (z, x) rotated by Rz(PH) Ry(TH)
N1 = (math.sin(TH) * math.cos(PH), math.sin(TH) * math.sin(PH), math.cos(TH))
N2 = (math.cos(TH) * math.cos(PH), math.cos(TH) * math.sin(PH), -math.sin(TH))


def numbers(binary: Path, scratch: Path, name: str, radial: str = "local", **kw) -> tuple[list[float], list[float], list[float]] | str:
    status, log = nc.run(binary, t3.af_deck(scratch, name, radial=radial, **kw))
    if status != 0:
        return f"{name}: run failed\n{log[-1500:]}"
    return t3.af_numbers(log)[0], nc.moments(log, 1), nc.moments(log, 2)


def worst(x, y) -> float:
    return max(abs(p - q) for a, b in zip(x, y) for p, q in zip(a, b))


def rotation(binary: Path, scratch: Path, kspace: bool = False) -> tuple[bool, str]:
    setups = {"A": dict(mom1=(0.0, 0.0, 1.0), mom2=(1.0, 0.0, 0.0)), "B": dict(mom1=(1.0, 0.0, 0.0), mom2=(0.0, 0.0, -1.0)),
              "C": dict(mom1=(1.0, 0.0, 0.0), mom2=(0.0, 0.0, 1.0)), "F": dict(mom1=N1, mom2=N2)}
    out = {k: numbers(binary, scratch, f"rot_{k}", nsp=3, kspace=kspace, nstep=3 if kspace else 2, **kw) for k, kw in setups.items()}
    for v in out.values():
        if isinstance(v, str):
            return False, v
    d = {k: worst(out["A"], out[k]) for k in "BCF"}
    return max(d.values()) < TOL, (f"{'|m1|, |m2| over 3 k-space steps' if kspace else 'E_F, |m1|, |m2| over 2 steps'} against the (z, x) pair: rotated 90 deg {d['B']:.1e}, sites exchanged "
                                   f"{d['C']:.1e}, general rotation {d['F']:.1e} (need < {TOL}); |m1| = {out['A'][1]}")


def parallel(binary: Path, scratch: Path) -> tuple[bool, str]:
    a = numbers(binary, scratch, "par_nsp3", nsp=3, mom1=(0.0, 0.0, 1.0), mom2=(0.0, 0.0, 1.0))
    b = numbers(binary, scratch, "par_nsp1", nsp=1, swap2=False)
    for v in (a, b):
        if isinstance(v, str):
            return False, v
    d = worst(a, b)
    return d < TOL, f"parallel nsp=3 against nsp=1 FM, E_F and |m| over 2 steps: {d:.1e} (need < {TOL})"


def engaged(binary: Path, scratch: Path) -> tuple[bool, str]:
    kw = dict(nsp=3, mom1=(0.0, 0.0, 1.0), mom2=(1.0, 0.0, 0.0))
    a = numbers(binary, scratch, "eng_loc", **kw)
    b = numbers(binary, scratch, "eng_spec", radial="spectral", **kw)
    for v in (a, b):
        if isinstance(v, str):
            return False, v
    d = worst(a[1:2], b[1:2])
    return d > ENGAGED, f"local against spectral, max ||m1| difference| over 2 steps {d:.1e} muB (need > {ENGAGED})"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    failed = False
    for name, check in (("rotation", rotation), ("kspace_rotation", lambda b, s: rotation(b, s, True)), ("parallel", parallel), ("engaged", engaged)):
        (base / name).mkdir()
        ok, info = check(binary, base / name)
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
