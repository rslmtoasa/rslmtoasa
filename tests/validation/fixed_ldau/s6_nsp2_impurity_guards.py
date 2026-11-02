#!/usr/bin/env python3
"""nsp=2 impurity/surface guards: each scenario must stop with its message (fatal)."""

from __future__ import annotations

import argparse
import re
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rs_embedding_ra2 as ra2  # noqa: E402


def set_line(path: Path, pattern: str, replacement: str) -> None:
    text = path.read_text()
    new, n = re.subn(pattern, replacement, text, count=1, flags=re.M)
    if n != 1:
        raise RuntimeError(f"{path}: pattern {pattern!r} not found")
    path.write_text(new)


def tilt_mom(path: Path) -> None:
    set_line(path, r"^\s*mom\s*=.*$", " mom = 1.0, 0.0, 0.0")


def expect_fatal(name: str, status: int, log: str, message: str) -> None:
    ok = status != 0 and message in log.lower()
    print(f"{'PASS' if ok else 'FAIL'} {name}: status={status}, message {'found' if ok else 'MISSING'}: {message!r}")
    if not ok:
        print(log[-1500:])
        raise SystemExit(1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True, type=Path)
    ap.add_argument("--scratch-root", type=Path)
    args = ap.parse_args()
    scratch = args.scratch_root or Path(tempfile.mkdtemp(prefix="s6_guards_"))
    scratch.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve()

    # Hosts: KS bulk Fe, U=2/J=0.5, one step, at nsp=1 and at nsp=2 (only the nsp line differs).
    hosts = {}
    for nsp in (1, 2):
        w = ra2.fresh_work(scratch, f"host_nsp{nsp}")
        ra2.configure(w, backend="ks", mode="local")
        set_line(w / "input.nml", r"^nsp\s*=\s*\d+", f"nsp = {nsp}")
        status, log = ra2.run(binary, w)
        if status != 0:
            raise RuntimeError(f"nsp={nsp} host run failed:\n{log[-3000:]}")
        hosts[nsp] = w

    def impurity(name: str, host: Path, *, imp_u: bool = False) -> Path:
        w = ra2.fresh_impurity(scratch, name)
        ra2.copy_potential(host, "Fe.nml", w)
        ra2.configure_impurity(w, "local")
        set_line(w / "input.nml", r"^nsp\s*=\s*\d+.*$", "nsp = 2")
        if imp_u:
            text = (w / "input.nml").read_text()
            (w / "input.nml").write_text(re.sub(r"(?im)^\s*&hamiltonian\s*$",
                "&hamiltonian\nhubbard_u_impurity(1,3) = 2.0\nhubbard_j_impurity(1,3) = 0.5\n", text, count=1))
        return w

    # (g1) nsp=2 impurity, frozen host written at nsp=1
    w = impurity("g1_host_nsp1", hosts[1])
    status, log = ra2.run(binary, w)
    expect_fatal("g1 nsp=2 impurity with an nsp=1 host file", status, log, "regenerate the host at nsp=2")

    # (g2) nsp=2 surface with U
    w = ra2.fresh_surface(scratch, "g2_surface")
    text = (w / "input.nml").read_text() + "\n&hamiltonian\nhubbard_u_general(1,3) = 2.0\nhubbard_j_general(1,3) = 0.0\n/\n"
    (w / "input.nml").write_text(text)
    status, log = ra2.run(binary, w, restart_name="Cu.nml")
    expect_fatal("g2 nsp=2 surface with U", status, log, "soc+u surface runs")

    # (g3) frozen host written along z (hubbard_axis recorded), run axis along x (the U-carrying impurity atom)
    w = impurity("g3_host_axis_mismatch", hosts[2], imp_u=True)
    tilt_mom(w / "Fe-imp.nml")
    status, log = ra2.run(binary, w)
    expect_fatal("g3 frozen host axis differs from the run axis", status, log, "regenerate the host along the run axis")
    print("s6_nsp2_impurity_guards: PASS")


if __name__ == "__main__":
    main()
