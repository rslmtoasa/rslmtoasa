#!/usr/bin/env python3
"""LdaUSurfaceHoh: hubbard_radial_moments = 'local' with hoh on the reduced Ni(001) surface (calctype S, deck tests/ldau/decks/ni001_soc_u3, U = 3, J = 0.7,
Chebyshev recursion, 4 steps, single thread). With hoh V_U is an additive onsite term in enim and the moments of H0 = H - V_U need the V_in that
was actually inserted. The cases:
  nsp2_vs_nsp1  nsp = 2 with SOC held at zero (spinor correction, V_in read by type) against nsp = 1 (scalar correction, V_in captured at the enim
                insertion): per-iteration Diff, Fermi energy, spin moment and spin-diagonal d density-matrix entries agree to SURF_SOC0_TOL.
                Bug caught: a V_in snapshot with its spin channels exchanged (the nsp = 1 run then differs by 2e-3). A run without the capture
                is bit-identical on this deck (V_in equals the stored V at the correction), so a missing capture is not caught here.
  engaged       the same surface with radial moments 'spectral' (hoh) gives spin moments that differ from the local run by more than 1e-3 muB at
                some step, so the correction acts (a skipped correction would make the first case pass trivially).
Not covered: an independent reference for the corrected moments themselves; nsp = 3 and 4 on surfaces.

usage: test_surface_hoh_u.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference_rows as rr  # noqa: E402

ENGAGED = 1.0e-3  # muB


def hoh_surface(binary: Path, base: Path, tag: str, nsp: int, mode: str):
    d = rr.surface_deck(base / tag, "host_soc0", nsp, nsp == 2)
    deck = d / "input.nml"
    for f in d.glob("*.nml"):  # the seed potentials carry the radial-moment mode they were written in
        f.write_text(f.read_text().replace("hubbard_radial_moments = 'local'", f"hubbard_radial_moments = '{mode}'"))
    text = deck.read_text()
    # the hoh Chebyshev recursion diverges in the deck's window [-1, 1.5] ('Chebyshev moments did not converge')
    text = re.sub(r"(?m)^energy_min\s*=.*$", "energy_min = -3.0", text)
    text = re.sub(r"(?m)^energy_max\s*=.*$", "energy_max = 3.0", text)
    deck.write_text(rr.in_hamiltonian(text, "hoh = .true.\n"))
    if rr.run(binary, d) != 0:
        return f"{tag}: surface run failed\n{rr.log_text(d)[-1200:]}"
    return d


def nsp2_vs_nsp1(a: Path, b: Path) -> tuple[bool, str]:
    return rr.compare_signatures(rr.signature(b, True), rr.signature(a, True), rr.SURF_SOC0_TOL)


def engaged(loc: Path, spec: Path) -> tuple[bool, str]:
    m = lambda d: [float(v) for v in re.findall(r"Spin moment of atom\s+\d+ is\s+(\S+)", rr.log_text(d))]
    ml, ms = m(loc), m(spec)
    if not ml or len(ml) != len(ms):
        return False, "no spin moments found, or different numbers of steps"
    dm = max(abs(a - b) for a, b in zip(ml, ms))
    return dm > ENGAGED, f"local against spectral: max |dm| over {len(ml)} values {dm:.1e} muB (need > {ENGAGED})"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    args = parser.parse_args()
    shutil.rmtree(args.scratch_root, ignore_errors=True)
    args.scratch_root.mkdir(parents=True)
    binary, base = args.binary.resolve(), args.scratch_root.resolve()
    runs = {tag: hoh_surface(binary, base, tag, nsp, mode)
            for tag, nsp, mode in (("nsp1_local", 1, "local"), ("nsp2_local", 2, "local"), ("nsp1_spectral", 1, "spectral"))}
    for tag, r in runs.items():
        if isinstance(r, str):
            print(f"FAIL {tag} {r}")
            return 1
    failed = False
    for name, (ok, info) in (("nsp2_vs_nsp1", nsp2_vs_nsp1(runs["nsp1_local"], runs["nsp2_local"])),
                             ("engaged", engaged(runs["nsp1_local"], runs["nsp1_spectral"]))):
        print(f"{'PASS' if ok else 'FAIL'} {name} {info}")
        failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
