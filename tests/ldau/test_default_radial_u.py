#!/usr/bin/env python3
"""LdaUDefaultRadial: hubbard_radial_moments = 'auto' (the default) with fixed U, and moving between 'spectral' and 'local' (bcc Fe nsp = 1 bulk, 2 steps,
U = 2 eV, J = 0; B2FeCo impurity in a converged k-space Fe host, 4 steps; single thread).
  default_is_local   a U deck without the key logs radial_mode=local radial_source=default, and its E_Hub_P_OUT and spin moments equal the explicit 'local' run
                     (bit for bit); the explicit 'spectral' run differs by more than 1e-4 muB in the moment (so the equality is not trivial)
  no_u_unchanged     a deck without U logs no Hubbard line, and its spin moments equal those of the explicit 'spectral' run
  switch_restart     a bulk restart from the potential and density-matrix files of the other mode, spectral to local and local to spectral, runs to the end, logs
                     'density matrix is reset', and a restart in the same mode does not
  switch_frozen_host the impurity run in 'local' with a frozen host written in 'spectral' runs to the end and logs the reset of the host's density matrix;
                     with a host in the same mode it does not
Bug caught: 'auto' resolving to the wrong mode or to local without U, a mode mismatch stopping the run (it did for frozen hosts), a mismatch accepted without
the reset and the warning, a reset applied when the modes agree. Not covered: the fallback of 'auto' to 'spectral' when local is not implemented for the input
(no runnable input reaches it: the inputs local rejects stop earlier), nsp >= 2, k-space, the quality of a state continued after the reset.

usage: test_default_radial_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference_rows as rr  # noqa: E402
import test_impurity_nsp3_u as imp  # noqa: E402

DIFFERS = 1.0e-4  # muB
RESET = "density matrix is reset"


def bulk(binary: Path, base: Path, name: str, mode: str | None, case: str = "bulk/bccFe_hubbard_nsp1", seed: Path | None = None):
    def edit(t):
        t = rr.set_nstep(t, 2)
        return t if mode is None else rr.in_hamiltonian(t, f"hubbard_radial_moments = '{mode}'\n")
    d = rr.prepare(rr.CASES / case, base / name, edit)
    if seed is not None:
        shutil.copy2(seed / "Fe_out.nml", d / "Fe.nml")
        shutil.copy2(seed / "Fe_dm.nml", d / "Fe_dm.nml")
    return d if rr.run(binary, d) == 0 else f"{name}: run failed\n{rr.log_text(d)[-1000:]}"


def numbers(d: Path):
    log = rr.log_text(d)
    return ([float(x) for x in re.findall(r"E_Hub_P_OUT=\s*(\S+)", log)], [float(x) for x in re.findall(r"Spin moment of atom\s+1 is\s+(\S+)", log)])


def bad(*runs):
    return next((r for r in runs if isinstance(r, str)), None)


def default_is_local(binary: Path, base: Path):
    d, loc, spec = bulk(binary, base, "u_default", None), bulk(binary, base, "u_local", "local"), bulk(binary, base, "u_spectral", "spectral")
    if bad(d, loc, spec):
        return False, bad(d, loc, spec)
    m = re.search(r"HUBBARD_RUN form=\S+ radial_mode=(\S+) radial_source=(\S+)", rr.log_text(d))
    (hd, md), (hl, ml), (_, ms) = numbers(d), numbers(loc), numbers(spec)
    dm = max(abs(a - b) for a, b in zip(md, ms))
    ok = bool(m) and m.groups() == ("local", "default") and (hd, md) == (hl, ml) and dm > DIFFERS
    return ok, f"log {m.groups() if m else 'no HUBBARD_RUN line'}; default equals explicit local: {(hd, md) == (hl, ml)}; default against spectral max |dm| {dm:.1e} muB (need > {DIFFERS})"


def no_u_unchanged(binary: Path, base: Path):
    d, spec = bulk(binary, base, "nou_default", None, "bulk/bccFe"), bulk(binary, base, "nou_spectral", "spectral", "bulk/bccFe")
    if bad(d, spec):
        return False, bad(d, spec)
    has_line = "HUBBARD_RUN" in rr.log_text(d)
    same = numbers(d)[1] == numbers(spec)[1]
    return not has_line and same, f"Hubbard line in the log: {has_line}; spin moments equal the explicit spectral run: {same}"


def switch_restart(binary: Path, base: Path):
    out = []
    for first, second in (("spectral", "local"), ("local", "spectral"), ("local", "local")):
        a = bulk(binary, base, f"seed_{first}_for_{second}", first)
        b = bulk(binary, base, f"{first}_to_{second}", second, seed=a) if not isinstance(a, str) else a
        if bad(a, b):
            return False, bad(a, b)
        out.append((f"{first} to {second}", RESET in rr.log_text(b)))
    ok = out[0][1] and out[1][1] and not out[2][1]
    return ok, "; ".join(f"{k}: reset logged {v}" for k, v in out)


def switch_frozen_host(binary: Path, base: Path):
    mismatch = imp.host_and_impurity(binary, base / "frozen", "mismatch", 1, "local", host_mode="spectral")
    same = imp.host_and_impurity(binary, base / "frozen", "same", 1, "local")
    if bad(mismatch, same):
        return False, bad(mismatch, same)
    a, b = RESET in rr.log_text(mismatch) and "Frozen U-carrying host" in rr.log_text(mismatch), RESET in rr.log_text(same)
    return a and not b, f"spectral host, local impurity: reset logged {a}; local host: reset logged {b}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    failed = False
    for name, case in (("default_is_local", default_is_local), ("no_u_unchanged", no_u_unchanged), ("switch_restart", switch_restart),
                       ("switch_frozen_host", switch_frozen_host)):
        ok, info = case(binary, base)
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
