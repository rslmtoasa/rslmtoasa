#!/usr/bin/env python3
"""Nid in the B1 convention from the deck's mirrored Nid.

Mirrored encoding (example/bulk/fccNiO/AF): Nid has mom = +z and its spin-1 channel holds the *minority* values.
B1 encoding: mom = -z and spin-1 holds the *majority* values, as for a site with mom(3) < 0 in the k-space GBT glue.
Every spin-indexed field of the &par namelist is swapped (index 1 <-> 2) and vrmax(2) = V_up - V_down changes sign.
Lines are rewritten textually, numbers are carried as the original strings; a spin-indexed key not in SPIN_AXIS and
present in the file stops the script.
usage: make_b1_nid.py Nid.nml Nid_B1.nml"""
import re
import sys

SPIN_AXIS = {  # key -> position of the spin index in the subscript list (0-based)
    "center_band": 1, "width_band": 1, "gravity_center": 1, "shifted_band": 1, "obar": 1, "c": 1, "enu": 1,
    "ppar": 1, "qpar": 1, "srdel": 1, "vl": 1, "pl": 1, "pnu": 1, "qi": 1, "dele": 1, "ql": 2, "ldm_flatten": 1,
}
PAIR_KEYS = ("xi_p", "xi_d", "rac", "hyper_field")          # two numbers (spin 1, spin 2) on one line
NO_SPIN = {"symbol", "atomic_number", "core", "valence", "sumec", "sumev", "etot", "utot", "ekin", "rhoeps", "ws_r",
           "vmad", "q10", "mom", "lmom", "hubbard_u", "hubbard_j", "f_core", "lmax", "num_quant_s", "num_quant_p",
           "num_quant_d", "num_quant_f"}
LINE = re.compile(r"^(\s*)(\w+)(\(([^)]*)\))?\s*=\s*(.*)$")


def main(src, dst):
    out, pending = [], {}
    for line in open(src).read().splitlines():
        m = LINE.match(line)
        if not m or m.group(2) in ("element", "par"):
            out.append(line)
            continue
        indent, key, sub, values = m.group(1), m.group(2), m.group(4), m.group(5)
        if key in SPIN_AXIS:
            axes = [s.strip() for s in sub.split(",")]
            axis = SPIN_AXIS[key]
            if axes[axis] in ("1", "2"):                     # one spin per line: collect and swap the pair
                pending.setdefault((key, tuple(a for i, a in enumerate(axes) if i != axis)), {})[axes[axis]] = (axes, values)
                out.append(("hold", key, tuple(a for i, a in enumerate(axes) if i != axis), axes[axis]))
                continue
            sys.exit(f"{key}({sub}): spin index is not a literal 1 or 2")
        elif key in PAIR_KEYS:
            v = [x.strip() for x in values.split(",")]
            out.append(f"{indent}{key} = {v[1]}, {v[0]}")
        elif key == "vrmax":
            avg, diff = [x.strip() for x in values.split(",")]
            diff = diff[1:] if diff.startswith("-") else "-" + diff
            out.append(f"{indent}vrmax = {avg}, {diff}")
        elif key == "mom":
            out.append(f"{indent}mom = 0.0000000000000000, 0.0000000000000000, -1.0000000000000000")
        elif key == "symbol":
            out.append(re.sub(r"'[^']*'", "'Nid'", line))
        elif key in NO_SPIN:
            out.append(line)
        else:
            sys.exit(f"unclassified key {key}: decide whether it is spin-indexed")
    final = []
    for item in out:
        if isinstance(item, tuple):
            _, key, rest, spin = item
            axes, _ = pending[(key, rest)][spin]
            other = "2" if spin == "1" else "1"
            axes2 = list(axes)
            axes2[SPIN_AXIS[key]] = spin
            values = pending[(key, rest)][other][1]            # the values of the other spin go on this line
            final.append(f" {key}({', '.join(axes2)}) = {values}")
        else:
            final.append(item)
    open(dst, "w").write("\n".join(final) + "\n")


main(*sys.argv[1:3])
