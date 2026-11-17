#!/usr/bin/env python3
"""LdaUImpurityEqualsBulk: an Fe impurity in an Fe host with identical potentials, cluster radius and settings reproduces the bulk atom
(bcc Fe, nsp = 1, block recursion, hoh, hubbard_radial_moments = 'local', U = 2 eV, J = 0, rc = 20, one step from the same converged state).
The impurity run (calctype I) keeps the Fermi energy fixed at its input value, so the deck must carry the bulk run's own Fermi energy
(the 'Free Fermi energy' of the bulk step); the spin moment and the potential parameters (Fe_out.nml against Fe-imp_out.nml) then agree.
  equals_bulk   |dm| < 1e-3 muB and max |parameter difference| < 1e-3 (measured 2e-4 and 2e-4). Bug caught: a wrong V_U / enim / moment
                correction in the impurity builder (build_locham) that the bulk builder does not have.
  engaged       the same impurity with the Fermi energy of another state (-0.042265 Ry, the value the unconverged deck carries) misses the
                bulk moment by more than 1e-2 muB (measured 9e-2), so the first case is sensitive to the embedding and does not pass trivially.
Not covered: spectral mode (3.5e-3 muB residual, KNOWN_ISSUES.md), U = 0 / first order (4e-3 residual), nsp = 2.

usage: test_impurity_equals_bulk.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference_rows as rr  # noqa: E402

RC = 20
STALE_FERMI = -0.042265
MOM_TOL, PAR_TOL, ENGAGED = 1.0e-3, 1.0e-3, 1.0e-2

IMPURITY_DECK = """&calculation
pre_processing = 'newclubulk'
verbose = T
/
&lattice
rc = {rc}
alat = 2.86120
crystal_sym = 'bcc'
wav = 1.40880
ntype = 2
nclu = 1
ct(:) = 2*3.0d0
r2 = 9.00d0
inclu(1, :) = 0.0, 0.0, 0.0
/
&atoms
database = './'
label(1) = 'Fe'
label(2) = 'Fe-imp'
/
&self
nstep = 1
/
&energy
fermi = {fermi}
energy_min = -2.0
energy_max = 0.8
channels_ldos = 2500
/
&control
calctype = 'I'
nsp = 1
lld = 21
recur = 'block'
/
&mix
beta = 0.1
mixtype = 'broyden'
/
&hamiltonian
hoh = .true.
hubbard_u_general(1,3) = 2.0
hubbard_j_general(1,3) = 0.0
hubbard_u_general(2,3) = 2.0
hubbard_j_general(2,3) = 0.0
hubbard_radial_moments = 'local'
/
"""


def moment(d: Path) -> float:
    m = re.findall(r"Spin moment of atom\s+1 is\s+(\S+)", rr.log_text(d))
    return float(m[-1])


def params(path: Path) -> list[float]:
    t = path.read_text()
    out = []
    for k in ("gravity_center", "center_band", "enu", "pl", "c"):
        for sp in (1, 2):
            m = re.search(rf"(?m)^ {k}\(:, {sp}\) = (.*)$", t)
            out += [float(x) for x in m.group(1).replace(",", " ").split()]
    return out


def impurity(binary: Path, base: Path, a: Path, name: str, fermi: float):
    c = base / name
    c.mkdir()
    shutil.copy2(a / "Fe_out.nml", c / "Fe.nml")
    shutil.copy2(a / "Fe_dm.nml", c / "Fe_dm.nml")
    t = (a / "Fe_out.nml").read_text()
    (c / "Fe-imp.nml").write_text(t.replace("symbol = 'Fe'", "symbol = 'Fe-imp'").replace("Fe_dm.nml", "Fe-imp_dm.nml"))
    (c / "Fe-imp_dm.nml").write_text((a / "Fe_dm.nml").read_text().replace("label = 'Fe'", "label = 'Fe-imp'"))
    (c / "input.nml").write_text(IMPURITY_DECK.format(rc=RC, fermi=fermi))
    return c if rr.run(binary, c) == 0 else f"{name}: impurity run failed\n{rr.log_text(c)[-1200:]}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    src = rr.CASES / "bulk/bccFe_hubbard_nsp1"

    def edit(t):
        t = re.sub(r"(?m)^hoh = .*", "hoh = .true.", t)
        t = re.sub(r"(?m)^rc = .*", f"rc = {RC}", t).replace("beta = 0.01", "beta = 0.1")
        return t.replace("hubbard_j_general(1,3) = 0.0", "hubbard_j_general(1,3) = 0.0\nhubbard_radial_moments = 'local'")

    a = rr.prepare(src, base / "A", lambda t: rr.set_nstep(edit(t), 2))
    if rr.run(binary, a) != 0:
        print(f"FAIL setup bulk A failed\n{rr.log_text(a)[-1200:]}")
        return 1
    b = rr.prepare(src, base / "B", edit)
    shutil.copy2(a / "Fe_out.nml", b / "Fe.nml")
    shutil.copy2(a / "Fe_dm.nml", b / "Fe_dm.nml")
    if rr.run(binary, b) != 0:
        print(f"FAIL setup bulk B failed\n{rr.log_text(b)[-1200:]}")
        return 1
    ef = re.findall(r"Free Fermi energy:\s+(\S+)", rr.log_text(b))
    if not ef:
        print("FAIL setup no 'Free Fermi energy' in the bulk log")
        return 1
    good = impurity(binary, base, a, "C_bulk_fermi", float(ef[-1]))
    bad = impurity(binary, base, a, "C_stale_fermi", STALE_FERMI)
    for x in (good, bad):
        if isinstance(x, str):
            print(f"FAIL setup {x}")
            return 1
    mb, mg, ms = moment(b), moment(good), moment(bad)
    pg = max(abs(x - y) for x, y in zip(params(b / "Fe_out.nml"), params(good / "Fe-imp_out.nml")))
    ok1 = abs(mg - mb) < MOM_TOL and pg < PAR_TOL
    ok2 = abs(ms - mb) > ENGAGED
    print(f"{'PASS' if ok1 else 'FAIL'} equals_bulk bulk {mb:.6f} impurity {mg:.6f} muB (E_F {ef[-1]} Ry), "
          f"|dm| {abs(mg - mb):.1e} (< {MOM_TOL}), max |parameter diff| {pg:.1e} (< {PAR_TOL})")
    print(f"{'PASS' if ok2 else 'FAIL'} engaged stale E_F {STALE_FERMI}: impurity {ms:.6f}, |dm| {abs(ms - mb):.1e} (need > {ENGAGED})")
    return 0 if ok1 and ok2 else 1


if __name__ == "__main__":
    sys.exit(main())
