#!/usr/bin/env python3
"""Offline exact-file and descriptive numeric comparison for C5.2 A2 outputs.

Usage: python3 compare_runs.py BASELINE_DIR CANDIDATE_DIR [--output PATH]

This is analysis tooling only. It does not establish an acceptance decision.
Relative differences are |candidate-baseline|/|baseline| for nonzero baselines;
zero baselines are counted and never divided by a floor.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from decimal import Decimal, InvalidOperation, getcontext
from pathlib import Path

getcontext().prec = 60

FILES = [
    "Fe_1_spinene.out",
    "Fe_dos.out",
    "Fe_orbene.out",
    "Fe_orbital_dos.out",
    "Fe_orbquadene.out",
    "Fe_out.nml",
    "Fe_scf.nml",
    "magneticdos.out",
    "radial_ground_state_Fe_1.dat",
]
NUM_RE = re.compile(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][+-]?\d+)?")
ASSIGN_RE = re.compile(r"^\s*([A-Za-z_]\w*(?:\([^)]*\))?)\s*=\s*(.*?)(?:\s*!.*)?$")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dec(token: str) -> Decimal:
    return Decimal(token.replace("D", "E").replace("d", "e"))


def numeric_tokens(text: str) -> list[str]:
    return NUM_RE.findall(text)


def summary(a: list[Decimal], b: list[Decimal]) -> dict:
    if len(a) != len(b):
        raise ValueError(f"numeric shape mismatch: {len(a)} != {len(b)}")
    ds = [abs(y - x) for x, y in zip(a, b)]
    rel = [d / abs(x) for x, d in zip(a, ds) if x != 0]
    zero_deltas = [d for x, d in zip(a, ds) if x == 0]
    idx = max(range(len(ds)), key=ds.__getitem__) if ds else None
    return {
        "count": len(ds),
        "max_abs": str(max(ds)) if ds else "0",
        "max_abs_index": idx,
        "baseline_at_max_abs": str(a[idx]) if idx is not None else None,
        "candidate_at_max_abs": str(b[idx]) if idx is not None else None,
        "max_rel_nonzero_baseline": str(max(rel)) if rel else None,
        "relative_defined_count": len(rel),
        "baseline_zero_count": len(zero_deltas),
        "baseline_zero_nonzero_candidate_count": sum(d != 0 for d in zero_deltas),
        "max_abs_at_baseline_zero": str(max(zero_deltas)) if zero_deltas else "0",
    }


def signed_zero_tokens(a: list[str], b: list[str]) -> dict:
    """Count spelling changes among corresponding numeric zero tokens."""
    if len(a) != len(b):
        raise ValueError(f"numeric token shape mismatch: {len(a)} != {len(b)}")
    changes = []
    for i, (x, y) in enumerate(zip(a, b)):
        dx, dy = dec(x), dec(y)
        if dx == 0 and dy == 0 and x != y:
            changes.append({"token_index": i, "baseline": x, "candidate": y})
    return {"numeric_token_count": len(a), "signed_zero_spelling_change_count": len(changes), "examples": changes[:20]}


def numeric_table(path: Path, names: list[str] | None = None) -> tuple[list[str], list[list[Decimal]], list[list[str]]]:
    rows: list[list[Decimal]] = []
    spellings: list[list[str]] = []
    for line in path.read_text().splitlines():
        t = line.strip()
        if not t or t.startswith("#"):
            continue
        tokens = t.split()
        try:
            row = [dec(x) for x in tokens]
        except (InvalidOperation, ValueError):
            continue
        rows.append(row)
        spellings.append(tokens)
    if not rows:
        raise ValueError(f"no numeric rows in {path}")
    lengths = {len(row) for row in rows}
    if len(lengths) != 1:
        raise ValueError(f"irregular numeric rows in {path}: {lengths}")
    ncol = len(rows[0])
    return names or [f"column_{i+1}" for i in range(ncol)], rows, spellings


def nml_fields(path: Path) -> dict[str, list[Decimal]]:
    out: dict[str, list[Decimal]] = {}
    for line in path.read_text().splitlines():
        m = ASSIGN_RE.match(line)
        if not m:
            continue
        vals = []
        for token in numeric_tokens(m.group(2)):
            try:
                vals.append(dec(token))
            except InvalidOperation:
                pass
        if vals:
            out[m.group(1).replace(" ", "")] = vals
    return out


def header_fields(path: Path) -> dict[str, list[Decimal]]:
    out = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#") or line.startswith("# columns:") or "=" not in line:
            continue
        k, v = line[1:].strip().split("=", 1)
        vals = []
        for token in numeric_tokens(v):
            try:
                vals.append(dec(token))
            except InvalidOperation:
                pass
        if vals:
            out[k.strip()] = vals
    return out


def read_report_moments(path: Path) -> dict:
    if not path.exists():
        return {"present": False}
    text = path.read_text()
    found = {}
    for line in text.splitlines():
        if "Total spin moment:" in line or re.search(r"Spin moment (?:of|projections of) atom", line):
            found[line.strip().split(":", 1)[0]] = [str(dec(x)) for x in numeric_tokens(line.split(":", 1)[1])]
    return {"present": True, "sha256": sha(path), "printed_precision_records": found}


def compare_pair(base: Path, cand: Path) -> dict:
    missing = [f for root in (base, cand) for f in FILES if not (root / f).is_file()]
    if missing:
        raise FileNotFoundError(f"required A3 files missing: {missing}")
    result = {
        "analysis_kind": "descriptive offline comparison; no acceptance verdict",
        "baseline_dir": str(base.resolve()),
        "candidate_dir": str(cand.resolve()),
        "method": "Decimal parsed fields; relative=|candidate-baseline|/|baseline| only when baseline != 0; baseline-zero cases reported with no denominator floor",
        "files": {},
        "quantity_metrics": {},
        "mapping": {
            "occupations": "Fe_out.nml/Fe_scf.nml ql(1,:,spin), entries s,p,d in array order",
            "spin_axis": "potential%mom is normalized/reference spin direction; not sphere-spin magnitude",
            "orbital_moment": "potential%lmom is non-normalized orbital angular-momentum vector",
            "sphere_spin": "radial_ground_state_Fe_1.dat header reported_moment_xyz_muB; source self.f90 copies potential%mom0 to magmom, and bands.f90 sets mom0=(mx,my,mz)",
            "energies": "Fe_out.nml/Fe_scf.nml: sumec,sumev,etot,utot,ekin,rhoeps",
            "fermi": "input.nml fermi assignment when present",
        },
    }
    for f in FILES:
        bp, cp = base / f, cand / f
        be = {"baseline_sha256": sha(bp), "candidate_sha256": sha(cp), "raw_bytes_identical": bp.read_bytes() == cp.read_bytes()}
        raw_tokens_a = numeric_tokens(bp.read_text())
        raw_tokens_b = numeric_tokens(cp.read_text())
        if len(raw_tokens_a) == len(raw_tokens_b):
            be["all_text_numeric_signed_zero_spelling"] = signed_zero_tokens(raw_tokens_a, raw_tokens_b)
        else:
            be["all_text_numeric_signed_zero_spelling"] = {
                "comparable": False,
                "baseline_numeric_token_count": len(raw_tokens_a),
                "candidate_numeric_token_count": len(raw_tokens_b),
                "reason": "whole-file numeric token counts differ; structured field/table comparisons remain separate",
            }
        if f.endswith(".nml"):
            x, y = nml_fields(bp), nml_fields(cp)
            common = sorted(set(x) & set(y))
            mismatch = {k: [len(x[k]), len(y[k])] for k in common if len(x[k]) != len(y[k])}
            missing_fields = {"baseline_only": sorted(set(x) - set(y)), "candidate_only": sorted(set(y) - set(x))}
            if mismatch:
                raise ValueError(f"numeric field shape mismatch in {f}: {mismatch}")
            be.update({"numeric_fields": {k: summary(x[k], y[k]) for k in common}, "field_names_only_one_side": missing_fields})
            for spin in (1, 2):
                k = f"ql(1,:,{spin})"
                if k in x and k in y:
                    result["quantity_metrics"][f"{f}:occupations_spin{spin}_s_p_d"] = summary(x[k], y[k])
            for field, label in (("mom", "spin_reference_axis"), ("lmom", "orbital_moment")):
                if field in x and field in y:
                    result["quantity_metrics"][f"{f}:{label}_components"] = summary(x[field], y[field])
            for field in ("sumec", "sumev", "etot", "utot", "ekin", "rhoeps"):
                if field in x and field in y:
                    result["quantity_metrics"][f"{f}:energy_{field}"] = summary(x[field], y[field])
        elif f == "radial_ground_state_Fe_1.dat":
            cols = next(l.split(":", 1)[1].strip().split() for l in bp.read_text().splitlines() if l.startswith("# columns:"))
            n1, r1, t1 = numeric_table(bp, cols)
            n2, r2, t2 = numeric_table(cp, cols)
            if n1 != n2 or len(r1) != len(r2) or any(len(x) != len(y) for x, y in zip(r1, r2)):
                raise ValueError("radial table row/column shape mismatch")
            be["rows"] = len(r1)
            be["columns"] = {n: summary([r[j] for r in r1], [r[j] for r in r2]) for j, n in enumerate(cols)}
            be["table_signed_zero_spelling"] = signed_zero_tokens([x for row in t1 for x in row], [x for row in t2 for x in row])
            h1, h2 = header_fields(bp), header_fields(cp)
            be["header_fields"] = {k: summary(h1[k], h2[k]) for k in sorted(set(h1) & set(h2)) if len(h1[k]) == len(h2[k])}
            for k in ("reported_moment_xyz_muB", "integrated_moment_muB"):
                if k in h1 and k in h2:
                    result["quantity_metrics"][k] = summary(h1[k], h2[k])
            be["header_fields_only_one_side"] = {"baseline_only": sorted(set(h1)-set(h2)), "candidate_only": sorted(set(h2)-set(h1))}
        else:
            names = None
            if f == "Fe_dos.out": names = ["energy", "site_dos_total", "site_dos_up", "site_dos_dw"]
            if f == "magneticdos.out": names = ["energy", "dos_up", "dos_dw", "dos_mx", "dos_my", "dos_mz", "dos_nmag"]
            n1, r1, t1 = numeric_table(bp, names)
            n2, r2, t2 = numeric_table(cp, names)
            if n1 != n2 or len(r1) != len(r2) or any(len(x) != len(y) for x, y in zip(r1, r2)):
                raise ValueError(f"numeric table row/column shape mismatch in {f}")
            be["rows"] = len(r1)
            be["columns"] = {n: summary([r[j] for r in r1], [r[j] for r in r2]) for j, n in enumerate(n1)}
            be["table_signed_zero_spelling"] = signed_zero_tokens([x for row in t1 for x in row], [x for row in t2 for x in row])
        result["files"][f] = be

    for run_name, root in (("baseline", base), ("candidate", cand)):
        p = root / "input.nml"
        if p.exists():
            m = re.search(r"(?im)^\s*fermi\s*=\s*([+-]?[\d.]+(?:[EeDd][+-]?\d+)?)", p.read_text())
            result.setdefault("deck_reference", {})[run_name] = {"input_nml_sha256": sha(p), "fermi_Ry": m.group(1) if m else None}
    if "deck_reference" in result:
        x = result["deck_reference"].get("baseline", {}).get("fermi_Ry")
        y = result["deck_reference"].get("candidate", {}).get("fermi_Ry")
        if x is not None and y is not None:
            dx, dy = dec(x), dec(y)
            result["quantity_metrics"]["E_F_Ry"] = {"baseline": str(dx), "candidate": str(dy), "absolute_difference": str(abs(dy-dx)), "relative_difference": str(abs(dy-dx)/abs(dx)) if dx else None, "relative_defined": dx != 0}
    result["report_out_spin_crosscheck"] = {"baseline": read_report_moments(base / "report.out"), "candidate": read_report_moments(cand / "report.out")}
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("baseline", type=Path)
    ap.add_argument("candidate", type=Path)
    ap.add_argument("--output", type=Path, help="write JSON here; defaults to stdout")
    args = ap.parse_args()
    try:
        data = compare_pair(args.baseline, args.candidate)
    except Exception as exc:
        print(f"compare_runs: ERROR: {exc}", file=sys.stderr)
        return 2
    raw = json.dumps(data, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(raw)
        print(args.output)
    else:
        sys.stdout.write(raw)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
