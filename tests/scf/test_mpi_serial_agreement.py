#!/usr/bin/env python3
"""MPI/serial agreement of the k-space SCF with a distributed k mesh.

Runs each listed case serially and with 2 and 3 MPI ranks through run_test.py (no reference
comparison) and requires byte-identical *_out.nml and final-DOS files and identical per-step Fermi level and
band energy log lines.  The k sums are reduced in global k order, so the result does not depend on
the rank count.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

CASES = ["Example_k_space_scf_tetra_diamondSi_sp_nsp1", "Example_k_space_scf_tetra_bccFe_u_nsp4_symon",
         "Example_k_space_scf_tetra_bccFe_nsp1_symon"]
RANKS = [2, 3]
DOS_FILES = ["totaldos.out", "magneticdos.out", "band_moments.dat", "projected_dos.dat"]
EF_LINE = re.compile(r"Canonical k-space occupations: EF=.*")


def run(binary: str, cases_json: str, case: str, scratch_root: Path, ranks: int) -> Path:
    name = f"{case}_np{ranks}"
    mode = ["--force-serial"] if ranks == 1 else ["--mpi-enabled", "--mpi-procs", str(ranks)]
    cmd = [sys.executable, str(Path(__file__).parent.parent / "run_test.py"), "--binary", binary,
           "--cases-json", cases_json, "--case-name", case, "--scratch-root", str(scratch_root),
           "--scratch-name", name] + mode
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise SystemExit(f"{name}: run failed\n{result.stdout[-2000:]}\n{result.stderr[-2000:]}")
    return scratch_root / name


def ef_lines(workdir: Path) -> list[str]:
    return [m.group(0) for line in (workdir / "testrun.log").read_text(errors="replace").splitlines()
            if (m := EF_LINE.search(line))]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True)
    parser.add_argument("--cases-json", required=True)
    parser.add_argument("--scratch-root", required=True)
    args = parser.parse_args()
    root = Path(args.scratch_root)

    failures = 0
    for case in CASES:
        serial = run(args.binary, args.cases_json, case, root, 1)
        serial_ef = ef_lines(serial)
        serial_files = sorted(serial.glob("*_out.nml")) + [serial / f for f in DOS_FILES if (serial / f).is_file()]
        if not serial_ef or not serial_files:
            raise SystemExit(f"{case}: serial run produced no EF log lines or *_out.nml")
        for ranks in RANKS:
            other = run(args.binary, args.cases_json, case, root, ranks)
            bad = [f.name for f in serial_files if not (other / f.name).is_file()
                   or f.read_bytes() != (other / f.name).read_bytes()]
            if ef_lines(other) != serial_ef:
                bad.append("EF/EBAND log lines")
            print(f"{case} np={ranks}: {len(serial_ef)} EF steps, {len(serial_files)} files, "
                  f"{'identical' if not bad else 'DIFFER ' + ', '.join(bad)}")
            failures += bool(bad)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
