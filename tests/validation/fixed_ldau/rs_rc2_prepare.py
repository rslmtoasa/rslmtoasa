#!/usr/bin/env python3
"""Prepare (but do not launch) C3 RC2 RS endpoint decks.

RC2 is conditional on the final C3 RC1 gate.  This helper makes the preregistered
U/J central-difference deck set from one frozen RS U=0 Fe seed; it deliberately
has no solver-launch command.  Once RC1 passes, the frozen manifest can be used
by the campaign owner to launch the cases with the required C1 trace and ladder.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE = Path("/private/tmp/rslmto-ra4-frozen/c1_trace_qlhistory_smoke_compact")
DEFAULT_SEED = Path("/private/tmp/rslmto-ra4-frozen/c1_trace_qlhistory_smoke_compact/Fe_out.nml")
DEFAULT_OUT = Path("/private/tmp/c3-rc2-rs-hf/staged")
ASSETS = ("clust", "map", "sbar", "view.sbar", "str.out", "ves.out", "mad.mat")
RY_PER_EV = 1.0 / 13.605703976


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def edit_group(text: str, group: str, assignments: dict[str, str | None]) -> str:
    lines = text.splitlines()
    start = next((i for i, line in enumerate(lines)
                  if line.strip().lower() == f"&{group}"), None)
    if start is None:
        lines.extend([f"&{group}", "/"])
        start = len(lines) - 2
    end = next((i for i in range(start + 1, len(lines)) if lines[i].strip() == "/"), None)
    if end is None:
        raise ValueError(f"unterminated &{group} namelist")
    pending = dict(assignments)
    kept: list[str] = []
    for line in lines[start + 1:end]:
        stripped = line.lstrip()
        key = stripped.split("=", 1)[0].strip().lower() if "=" in stripped and not stripped.startswith("!") else ""
        if key not in pending:
            kept.append(line)
            continue
        value = pending.pop(key)
        if value is not None:
            kept.append(f"{key} = {value}")
    kept.extend(f"{key} = {value}" for key, value in pending.items() if value is not None)
    lines[start + 1:end] = kept
    return "\n".join(lines) + "\n"


def cases() -> list[dict]:
    out = []
    # U derivative: U=2 eV at h=.05 and .10 eV, J=0.
    for u in (1.90, 1.95, 2.00, 2.05, 2.10):
        for mode in ("local", "spectral"):
            out.append({"family": "U", "mode": mode, "u_eV": u, "j_eV": 0.0})
    # J derivatives: both prescribed centers and both h=.05/.10 stencils.
    for center in (0.5, 0.9):
        for j in sorted({round(center + delta, 8) for delta in (-0.10, -0.05, 0.0, 0.05, 0.10)}):
            for mode in ("local", "spectral"):
                out.append({"family": f"J{center:g}", "mode": mode, "u_eV": 2.0, "j_eV": j})
    return out


def stage(source: Path, seed: Path, outdir: Path) -> dict:
    if outdir.exists():
        raise FileExistsError(f"refusing to overwrite staged campaign: {outdir}")
    for asset in ASSETS:
        if not (source / asset).is_file():
            raise FileNotFoundError(source / asset)
    if not (source / "input.nml").is_file() or not seed.is_file():
        raise FileNotFoundError("source input or accepted seed is missing")
    outdir.mkdir(parents=True)
    base = (source / "input.nml").read_text()
    records = []
    for c in cases():
        name = f"{c['family'].lower()}_{c['mode']}_u{c['u_eV']:.2f}_j{c['j_eV']:.2f}".replace(".", "p")
        d = outdir / name
        d.mkdir()
        deck = edit_group(base, "self", {"nstep": "200", "conv_thr": "1.0d-30"})
        deck = edit_group(deck, "hamiltonian", {
            "hubbard_u_general(1,3)": f"{c['u_eV']:.8f}",
            "hubbard_j_general(1,3)": f"{c['j_eV']:.8f}",
            "hubbard_u_potential_form": "'liechtenstein'",
            "hubbard_radial_moments": f"'{c['mode']}'",
            # A zero-amplitude W1 probe causes the existing diagnostic writer
            # to emit the exact full d-shell P_OUT matrix.  lambda=0 inserts
            # no potential and leaves the physical Hamiltonian unchanged.
            "hubbard_test_operator": "'w1'",
            "hubbard_test_lambda": "0.0",
            "hubbard_test_spin": "0",
        })
        deck = edit_group(deck, "mix", {"beta": "0.01", "mixtype": "'broyden'"})
        (d / "input.nml").write_text(deck)
        shutil.copyfile(seed, d / "Fe_out.nml")
        for asset in ASSETS:
            shutil.copyfile(source / asset, d / asset)
        records.append({
            **c, "case": name, "directory": str(d),
            "input_sha256": sha(d / "input.nml"),
            "seed_sha256": sha(d / "Fe_out.nml"),
            "asset_sha256": {a: sha(d / a) for a in ASSETS},
            "status": "STAGED_NOT_LAUNCHED_RC1_GATE_PENDING",
        })
    manifest = {
        "schema": "c3-rc2-rs-prepared-decks-v1",
        "status": "STAGED_ONLY; DO NOT LAUNCH BEFORE RC1 GATE PASS",
        "source_input": str(source / "input.nml"),
        "source_input_sha256": sha(source / "input.nml"),
        "seed_source": str(seed), "seed_source_sha256": sha(seed),
        "source_asset_sha256": {a: sha(source / a) for a in ASSETS},
        "stencils": {
            "u_center_eV": 2.0, "u_h_eV": [0.05, 0.10],
            "j_centers_eV": [0.5, 0.9], "j_h_eV": [0.05, 0.10],
            "units": "inputs in eV; all finite-difference slopes in Ry/eV",
        },
        "execution": {
            "nstep_cap": 200, "measurement_conv_thr": 1e-30,
            "initial_mixing": "Broyden beta=.01; C1 fallback ladder remains at most beta=.005 then Broyden beta=.01, each up to 200 steps",
            "capture": "RSLMTO_FIXED_U_TRACE_QL_ONLY=1; preserve complete P_OUT trace and endpoint files",
            "threads": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"},
            "seed_policy": "Each case starts from an independent copy of the same accepted Fe U=0 RS potential. Confirm physical recipe and the RC1 gate before launch.",
        },
        "case_count": len(records), "cases": records,
    }
    path = outdir / "RC2_STAGED_DECKS.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--seed", type=Path, default=DEFAULT_SEED)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    manifest = stage(args.source, args.seed, args.outdir)
    print(json.dumps({"manifest": str(args.outdir / "RC2_STAGED_DECKS.json"),
                      "case_count": manifest["case_count"],
                      "status": manifest["status"]}, indent=2))


if __name__ == "__main__":
    main()
