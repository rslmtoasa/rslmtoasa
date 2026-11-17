#!/usr/bin/env python3
"""Make RS-LMTO-ASA *_out.nml restart files non-magnetic (spin-collinear averaged).

Reads a potential namelist file (e.g. Fe_out.nml) and rewrites every
spin-resolved quantity so that the up and down channels are identical
(their average), and zeroes out the magnetic-moment vectors. Everything
else (element block, scalars, hubbard_u/j, lmax, ...) is left untouched.

Usage:
    python3 demagnetize_nml.py Fe_out.nml Pt_out.nml
    python3 demagnetize_nml.py --in-place Fe_out.nml
    python3 demagnetize_nml.py -o Fe_out_nm.nml Fe_out.nml
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

# Arrays shaped (..., spin) with pattern "name(:, <spin>) = v1, v2, ..."
TWO_D_LAST_SPIN = {
    "center_band", "width_band", "gravity_center", "shifted_band", "obar",
    "c", "enu", "ppar", "qpar", "srdel", "vl", "pl",
}

# Scalars/vectors written as a single line "name = v_spin1, v_spin2"
ONE_LINE_SPIN_PAIR = {"vrmax", "xi_p", "xi_d", "rac", "hyper_field"}

# Magnetic-moment vectors/scalars zeroed out for a non-magnetic state
ZERO_OUT = {"mom", "lmom", "mx", "my", "mz", "mtot", "lx", "ly", "lz", "ltot"}

LINE_RE = re.compile(r"^(?P<indent>\s*)(?P<name>[A-Za-z_][A-Za-z0-9_]*)"
                      r"(?:\((?P<index>[^)]*)\))?\s*=\s*(?P<values>.*)$")


def parse_floats(values_str: str) -> list[float]:
    parts = [p.strip() for p in values_str.split(",")]
    return [float(p.replace("D", "E").replace("d", "e")) for p in parts]


def fmt_floats(values: list[float]) -> str:
    return ", ".join(f"{v:.16E}" for v in values)


def process(lines: list[str]) -> list[str]:
    parsed = []
    for line in lines:
        m = LINE_RE.match(line)
        parsed.append(m)

    out = list(lines)

    # --- name(:, spin) = v1, v2, ... arrays -------------------------------
    groups: dict[str, dict[str, int]] = {}
    for i, m in enumerate(parsed):
        if not m or m.group("name") not in TWO_D_LAST_SPIN or not m.group("index"):
            continue
        idx_parts = [p.strip() for p in m.group("index").split(",")]
        if len(idx_parts) != 2 or idx_parts[1] not in ("1", "2"):
            continue
        groups.setdefault(m.group("name"), {})[idx_parts[1]] = i

    for name, spins in groups.items():
        if "1" not in spins or "2" not in spins:
            continue
        i1, i2 = spins["1"], spins["2"]
        v1 = parse_floats(parsed[i1].group("values"))
        v2 = parse_floats(parsed[i2].group("values"))
        avg = [(a + b) / 2.0 for a, b in zip(v1, v2)]
        for i, m in ((i1, parsed[i1]), (i2, parsed[i2])):
            out[i] = f"{m.group('indent')}{m.group('name')}({m.group('index')}) = {fmt_floats(avg)}\n"

    # --- ql(orb, :, spin) --------------------------------------------------
    ql_groups: dict[str, dict[str, int]] = {}
    for i, m in enumerate(parsed):
        if not m or m.group("name") != "ql" or not m.group("index"):
            continue
        idx_parts = [p.strip() for p in m.group("index").split(",")]
        if len(idx_parts) != 3:
            continue
        orb, mid, spin = idx_parts
        if spin not in ("1", "2"):
            continue
        ql_groups.setdefault(orb, {})[spin] = i

    for orb, spins in ql_groups.items():
        if "1" not in spins or "2" not in spins:
            continue
        i1, i2 = spins["1"], spins["2"]
        v1 = parse_floats(parsed[i1].group("values"))
        v2 = parse_floats(parsed[i2].group("values"))
        avg = [(a + b) / 2.0 for a, b in zip(v1, v2)]
        for i in (i1, i2):
            m = parsed[i]
            out[i] = f"{m.group('indent')}{m.group('name')}({m.group('index')}) = {fmt_floats(avg)}\n"

    # --- ldm_flatten(l, spin, :) --------------------------------------------
    ldm_groups: dict[str, dict[str, int]] = {}
    for i, m in enumerate(parsed):
        if not m or m.group("name") != "ldm_flatten" or not m.group("index"):
            continue
        idx_parts = [p.strip() for p in m.group("index").split(",")]
        if len(idx_parts) != 3:
            continue
        l_idx, spin, tail = idx_parts
        if spin not in ("1", "2"):
            continue
        ldm_groups.setdefault(l_idx, {})[spin] = i

    for l_idx, spins in ldm_groups.items():
        if "1" not in spins or "2" not in spins:
            continue
        i1, i2 = spins["1"], spins["2"]
        v1 = parse_floats(parsed[i1].group("values"))
        v2 = parse_floats(parsed[i2].group("values"))
        avg = [(a + b) / 2.0 for a, b in zip(v1, v2)]
        for i in (i1, i2):
            m = parsed[i]
            out[i] = f"{m.group('indent')}{m.group('name')}({m.group('index')}) = {fmt_floats(avg)}\n"

    # --- name = v_spin1, v_spin2 (single line, no index) --------------------
    for i, m in enumerate(parsed):
        if not m or m.group("name") not in ONE_LINE_SPIN_PAIR or m.group("index"):
            continue
        vals = parse_floats(m.group("values"))
        if len(vals) != 2:
            continue
        avg = sum(vals) / 2.0
        out[i] = f"{m.group('indent')}{m.group('name')} = {fmt_floats([avg, avg])}\n"

    # --- moment vectors/scalars -> zero --------------------------------------
    for i, m in enumerate(parsed):
        if not m or m.group("name") not in ZERO_OUT or m.group("index"):
            continue
        vals = parse_floats(m.group("values"))
        out[i] = f"{m.group('indent')}{m.group('name')} = {fmt_floats([0.0] * len(vals))}\n"

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="Input *_out.nml file(s)")
    ap.add_argument("-o", "--output", help="Output path (only valid with a single input file)")
    ap.add_argument("--in-place", action="store_true", help="Overwrite the input file(s)")
    ap.add_argument("--suffix", default="_nm", help="Suffix for output files (default: _nm)")
    args = ap.parse_args()

    if args.output and len(args.files) > 1:
        ap.error("-o/--output can only be used with a single input file")

    for fname in args.files:
        path = Path(fname)
        lines = path.read_text().splitlines(keepends=True)
        new_lines = process(lines)

        if args.in_place:
            out_path = path
        elif args.output:
            out_path = Path(args.output)
        else:
            out_path = path.with_name(path.stem + args.suffix + path.suffix)

        out_path.write_text("".join(new_lines))
        print(f"{path} -> {out_path}")


if __name__ == "__main__":
    main()
