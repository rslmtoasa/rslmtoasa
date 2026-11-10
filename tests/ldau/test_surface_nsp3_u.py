#!/usr/bin/env python3
"""LdaUSurfaceNsp3: hubbard_radial_moments = 'local' at nsp = 3 on the reduced Ni(001) surface (calctype S, deck tests/ldau/decks/ni001_soc_u3, U = 3,
J = 0.7, Chebyshev recursion, 4 steps, single thread, moments along +z, density_policy = 'constrained_spiral'). The nsp = 3 host is the nsp = 2 host
with hubbard_nsp relabelled. The spinor correction (V by type) must reproduce the nsp = 1 scalar correction: per-iteration Diff, Fermi energy, spin moment
and spin-diagonal d density-matrix entries agree to SURF_SOC0_TOL.
  nsp3_vs_nsp1      first order
  nsp3_vs_nsp1_hoh  hoh (V_U in enim; the nsp = 1 side captures V_in at the insertion)
  engaged           nsp = 3 local against spectral: spin moments differ by more than 1e-3 muB at some step
With the default density_policy ('relaxed_reference') the empty sphere's small moment turns negative after step 1 and its axis flips to -z at nsp = 3,
and the runs then differ from nsp = 1 by 2.6e-2 (also in spectral mode; identical at U = 0, where that moment stays positive): see KNOWN_ISSUES.md.
Bug caught: the correction missing or in the wrong spin frame on a surface cluster at nsp = 3. Not covered: nsp = 4, SOC, noncollinear axes, an
independent reference for the corrected moments.

usage: test_surface_nsp3_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference_rows as rr  # noqa: E402
import test_surface_hoh_u as sh  # noqa: E402


def surface(binary: Path, base: Path, tag: str, nsp: int, mode: str, hoh: bool):
    d = rr.surface_deck(base / tag, "host_soc0", nsp, False)
    for f in d.glob("*.nml"):  # the seed potentials carry the radial-moment mode they were written in
        f.write_text(f.read_text().replace("hubbard_radial_moments = 'local'", f"hubbard_radial_moments = '{mode}'"))
    if nsp == 3:
        pot = d / "Ni.nml"
        pot.write_text(pot.read_text().replace("hubbard_nsp = 2", "hubbard_nsp = 3"))
    deck = d / "input.nml"
    text = deck.read_text()
    if nsp == 3:
        text = text.replace("&control\n", "&control\ndensity_policy = 'constrained_spiral'\n", 1)
    if hoh:  # the hoh Chebyshev recursion needs the wider window
        text = re.sub(r"(?m)^energy_min\s*=.*$", "energy_min = -3.0", text)
        text = re.sub(r"(?m)^energy_max\s*=.*$", "energy_max = 3.0", text)
        text = rr.in_hamiltonian(text, "hoh = .true.\n")
    deck.write_text(text)
    if rr.run(binary, d) != 0:
        return f"{tag}: surface run failed\n{rr.log_text(d)[-1200:]}"
    return d


def parity(binary: Path, base: Path, tag: str, hoh: bool):
    a, b = surface(binary, base, f"{tag}_nsp1", 1, "local", hoh), surface(binary, base, f"{tag}_nsp3", 3, "local", hoh)
    for x in (a, b):
        if isinstance(x, str):
            return False, x
    return rr.compare_signatures(rr.signature(b, True), rr.signature(a, True), rr.SURF_SOC0_TOL)


def engaged(binary: Path, base: Path):
    loc, spec = surface(binary, base, "eng_local", 3, "local", False), surface(binary, base, "eng_spectral", 3, "spectral", False)
    for x in (loc, spec):
        if isinstance(x, str):
            return False, x
    return sh.engaged(loc, spec)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    failed = False
    for name, result in (("nsp3_vs_nsp1", parity(binary, base, "p", False)), ("nsp3_vs_nsp1_hoh", parity(binary, base, "ph", True)),
                         ("engaged", engaged(binary, base))):
        ok, info = result
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
