#!/usr/bin/env python3
"""VAL-16: current-kernel GBT/commensurate-supercell comparison.

The historical WP9 decks contain potentials produced by an older GBT kernel.
This campaign first converges each explicit magnetic supercell with the
current executable, then uses those output potentials for both the frozen
force-theorem and independently converged SCF comparisons.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WP9 = ROOT / "tests" / "regression" / "wp9_validation" / "commensurate_supercell"

REPORT_RE = re.compile(r"Band energy of system:\s*([-+0-9.eE]+)")
FERMI_RE = re.compile(r"Free Fermi energy:\s*([-+0-9.eE]+)")
Q0_RE = re.compile(r"ql\(1,\s*:,\s*(\d+)\)\s*=\s*([^\n]+)")


def run(binary: Path, cwd: Path, timeout: int) -> str:
    result = subprocess.run([str(binary)], cwd=cwd, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=timeout)
    (cwd / "run.log").write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"run failed in {cwd}:\n{result.stdout[-4000:]}")
    return result.stdout


def parse_run(workdir: Path, nsite: int) -> dict:
    report = (workdir / "report.out").read_text()
    log = (workdir / "run.log").read_text()
    moments = []
    charges = []
    for path in sorted(workdir.glob("Fe*_out.nml")):
        text = path.read_text()
        q = {int(m.group(1)): [float(v.replace("D", "E")) for v in m.group(2).split(",")]
             for m in Q0_RE.finditer(text)}
        up = q[1][2]
        down = q[2][2]
        moments.append(up - down)
        charges.append(sum(q[1]) + sum(q[2]))
    if len(moments) != nsite:
        raise RuntimeError(f"expected {nsite} output potentials in {workdir}, got {len(moments)}")
    ef = [float(m.group(1)) for m in FERMI_RE.finditer(log)]
    return {
        "nsite": nsite,
        "charge": charges,
        "moment": moments,
        "eband_per_atom": float(REPORT_RE.search(report).group(1)) / nsite,
        "ef": ef[-1] if ef else None,
    }


def copy_reference_potentials(reference: Path, target: Path, nsite: int) -> None:
    for i in range(1, nsite + 1):
        shutil.copy2(reference / f"Fe{i}_out.nml", target / f"Fe{i}.nml")


def make_gbt_potential(reference: Path, target: Path) -> None:
    text = (reference / "Fe1_out.nml").read_text()
    text = re.sub(r"symbol\s*=\s*'Fe1'", "symbol = 'Fe'", text, count=1)
    text = re.sub(r"mom\s*=\s*[-+0-9.eEdD,\s]+(?=\n)", "mom = 0.0, 0.0, 1.0", text, count=1)
    (target / "Fe.nml").write_text(text)


def make_gbt_input(template: Path) -> str:
    return template.read_text()


def compare(gbt: dict, sup: dict) -> dict:
    return {
        "d_charge": max(abs(gbt["charge"][0] - x) for x in sup["charge"]),
        "d_moment": max(abs(gbt["moment"][0] - x) for x in sup["moment"]),
        "d_eband": abs(gbt["eband_per_atom"] - sup["eband_per_atom"]),
        "d_ef": None if gbt["ef"] is None or sup["ef"] is None else abs(gbt["ef"] - sup["ef"]),
    }


def one_case(binary: Path, scratch: Path, qname: str, nsite: int, qvalue: str) -> None:
    source = WP9 / "gbt_supercell" / qname
    shutil.rmtree(scratch, ignore_errors=True)
    scratch.mkdir(parents=True)
    reference = scratch / "reference_supercell"
    shutil.copytree(source / "super_scf", reference)
    # The historical input files are only starting guesses.  The state used
    # below is always the output of this current executable's SCF run.
    reference_input = reference / "input.nml"
    reference_input.write_text(re.sub(r"nstep\s*=\s*\d+", "nstep = 100", reference_input.read_text(), count=1))
    reference_log = run(binary, reference, 7200)
    if "Converged!" not in reference_log:
        raise RuntimeError(f"current-kernel reference did not converge in {reference}")
    ref = parse_run(reference, nsite)

    rows = []
    for leg in ("mft", "scf"):
        super_run = scratch / f"super_{leg}"
        gbt_run = scratch / f"gbt_{leg}"
        shutil.copytree(source / ("super" if leg == "mft" else "super_scf"), super_run)
        shutil.copytree(source / ("gbt" if leg == "mft" else "gbt_scf"), gbt_run)
        copy_reference_potentials(reference, super_run, nsite)
        make_gbt_potential(reference, gbt_run)
        (gbt_run / "input.nml").write_text(make_gbt_input(source / ("gbt" if leg == "mft" else "gbt_scf") / "input.nml"))
        run(binary, super_run, 7200)
        run(binary, gbt_run, 7200)
        super_state = parse_run(super_run, nsite)
        gbt_state = parse_run(gbt_run, 1)
        rows.append((leg, gbt_state, super_state))

    print(f"=== VAL-16 {qname}: q=(0,0,{qvalue}), {nsite}-atom explicit supercell ===")
    print("reference supercell: charge=%s moment=%s EF=%s Eband/atom=%.12f" %
          (" ".join(f"{x:.9f}" for x in ref["charge"]),
           " ".join(f"{x:.9f}" for x in ref["moment"]), ref["ef"], ref["eband_per_atom"]))
    for leg, gbt, sup in rows:
        d = compare(gbt, sup)
        print(f"{leg}: GBT Eband/atom={gbt['eband_per_atom']:.12f} super={sup['eband_per_atom']:.12f} dE={d['d_eband']:.3e}")
        print(f"{leg}: GBT charge={gbt['charge'][0]:.9f} super spread={max(sup['charge'])-min(sup['charge']):.3e} dQ={d['d_charge']:.3e}")
        print(f"{leg}: GBT moment={gbt['moment'][0]:.9f} super spread={max(sup['moment'])-min(sup['moment']):.3e} dM={d['d_moment']:.3e}")
        print(f"{leg}: EF GBT={gbt['ef']} super={sup['ef']} dEF={d['d_ef']}")
    print("reference run log: %s" % (reference / "run.log"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scratch-root", type=Path, required=True)
    args = parser.parse_args()
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    one_case(args.binary, args.scratch_root / "q050", "q050", 4, "0.5")
    one_case(args.binary, args.scratch_root / "q033", "q033", 6, "0.3333333333333333")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
