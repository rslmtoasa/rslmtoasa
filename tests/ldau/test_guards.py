#!/usr/bin/env python3
"""LDA+U guard scenarios in one test: each must stop the run with its message.

Scenarios (each prints PASS/FAIL by name; the test fails if any fails):
  surface_nsp1_host              nsp = 2 surface with a frozen U host written at nsp = 1
  surface_u_less_host            U requested on a type whose host was converged without U
  tilted_axis_block              tilted spin axis with block recursion
  mixed_axes                     different spin axes on different sites at nsp = 2
  impurity_host_nsp1             nsp = 2 impurity with a host written at nsp = 1
  impurity_host_axis_mismatch    frozen host written along z, run axis along x

Bugs caught: a guard that stops firing (the run would continue with an unsupported U path), or one that
fires with a different message.

usage: test_guards.py --binary rslmto.x --scratch-root DIR
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guard_helpers as gh  # noqa: E402

CASES_JSON = gh.ROOT / "tests/scf/cases.json"
RUN_TEST = gh.ROOT / "tests/run_test.py"
# scenario name -> case name in tests/scf/cases.json (the case carries the expected fatal message)
CASE_GUARDS = {
    "surface_nsp1_host": "Example_surface_soc_u_nsp1_host_guard",
    "surface_u_less_host": "Example_surface_fccNi001_uhost_guard",
    "tilted_axis_block": "Example_bulk_bccFe_soc_u_tilted_block_guard",
    "mixed_axes": "Example_bulk_Pt2MnGa_soc_u_mixed_axes_guard",
}


def run_case_guard(case: str, binary: Path, scratch: Path) -> tuple[bool, str]:
    proc = subprocess.run(
        [sys.executable, str(RUN_TEST), "--binary", str(binary), "--cases-json", str(CASES_JSON),
         "--case-name", case, "--scratch-root", str(scratch), "--force-serial"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return proc.returncode == 0, proc.stdout


def expect_fatal(status: int, log: str, message: str) -> tuple[bool, str]:
    ok = status != 0 and message in log.lower()
    return ok, "" if ok else f"status={status}, message {message!r} not found\n{log[-1500:]}"


def impurity_guards(binary: Path, scratch: Path) -> dict[str, tuple[bool, str]]:
    results: dict[str, tuple[bool, str]] = {}
    hosts = {}
    for nsp in (1, 2):
        work = gh.fresh(gh.HOST_FIXTURE, scratch, f"host_nsp{nsp}")
        gh.configure_host(work)
        gh.set_line(work / "input.nml", r"^nsp\s*=\s*\d+", f"nsp = {nsp}")
        status, log = gh.run_binary(binary, work)
        if status != 0:
            fail = (False, f"host run at nsp={nsp} failed\n{log[-1500:]}")
            return {"impurity_host_nsp1": fail, "impurity_host_axis_mismatch": fail}
        hosts[nsp] = work

    work = gh.fresh(gh.IMPURITY_FIXTURE, scratch, "g1_host_nsp1")
    gh.copy_host_potential(hosts[1], work)
    gh.configure_impurity(work, impurity_u=False)
    results["impurity_host_nsp1"] = expect_fatal(*gh.run_binary(binary, work), "regenerate the host at nsp=2")

    work = gh.fresh(gh.IMPURITY_FIXTURE, scratch, "g3_host_axis_mismatch")
    gh.copy_host_potential(hosts[2], work)
    gh.configure_impurity(work, impurity_u=True)
    gh.set_line(work / "Fe-imp.nml", r"^\s*mom\s*=.*$", " mom = 1.0, 0.0, 0.0")
    results["impurity_host_axis_mismatch"] = expect_fatal(*gh.run_binary(binary, work), "regenerate the host along the run axis")
    return results


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True, type=Path)
    ap.add_argument("--scratch-root", type=Path)
    args = ap.parse_args()
    scratch = args.scratch_root or Path(tempfile.mkdtemp(prefix="ldau_guards_"))
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)
    binary = args.binary.resolve()

    results: dict[str, tuple[bool, str]] = {}
    for scenario, case in CASE_GUARDS.items():
        results[scenario] = run_case_guard(case, binary, scratch / "cases")
    results.update(impurity_guards(binary, scratch / "impurity"))

    for name, (ok, detail) in results.items():
        print(f"{'PASS' if ok else 'FAIL'} {name}")
        if not ok:
            print("   | " + detail.replace("\n", "\n   | "))
    failed = [n for n, (ok, _) in results.items() if not ok]
    print(f"guard scenarios: {len(results) - len(failed)} passed, {len(failed)} failed" + (f" ({', '.join(failed)})" if failed else ""))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
