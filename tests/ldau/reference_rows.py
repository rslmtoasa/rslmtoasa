#!/usr/bin/env python3
"""LDA+U stored-reference rows.

Each row runs a short (4 SCF steps, single thread) deck with the candidate binary and compares scalar
quantities by value against the JSON in tests/ldau/references/<row>.json, within the tolerances below:
per-iteration Diff, Fermi energy, spin moment and routed energy from the log, the dm sidecar entries, and
(for the no-U rows) every number of the *_out.nml potential files.

Two kinds of row:
  abs   the candidate signature against the stored signature of the same deck.
  soc0  the nsp=2 deck with SOC held at zero (soc_scale=0, fix_soc) against the nsp=1 deck with the same U,
        both run now: iteration 1 within SOC0_TOL, all iterations within the row tolerance. The stored JSON
        of a soc0 row is an archive of the nsp=1 run (it is not read by the test).

Rows registered with ctest are listed in tests/register_tests.cmake; the others stay here as archived data.
References are rewritten only with an explicit reason:
    reference_rows.py ROW --binary BIN --update "one-line reason" [--variant]
which runs this binary and stores its numbers (the reason, binary hash and date go into the JSON);
--variant writes <row>.<platform>.json, which this platform then uses instead of <row>.json.

usage: reference_rows.py ROW --binary rslmto.x [--scratch-root DIR] [--update REASON] | --list
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CASES = ROOT / "tests/scf/cases"
SURFACE = ROOT / "tests/ldau/decks/ni001_soc_u3"
REFERENCES = ROOT / "tests/ldau/references"
NSTEP = 4
LOCAL = "hubbard_radial_moments = 'local'\n"
ABS_TOL = 1e-10          # abs rows: stored by value (the original byte/occupation equality is not portable)
SOC0_TOL = 1e-9          # soc0 rows, iteration 1 (the logged Diff has 10 decimals)
IMP_SOC0_TOL = 1e-7      # soc0 rows, all 4 iterations: impurity (local-mode transient amplification, KNOWN_ISSUES)
BULK_SOC0_TOL = 1e-9
SURF_SOC0_TOL = 1e-7     # surface, 4 iterations
OFFDIAG_TOL = 1e-12      # max |P_up,dn| at SOC = 0
SINGLE_THREAD = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "RSLMTO_OMP_THREADS_SERIAL": "1"}


def platform_name() -> str:
    machine = platform.machine().lower()
    arch = "x86_64" if machine in {"x86_64", "amd64"} else "arm64" if machine in {"arm64", "aarch64"} else machine
    return f"{'macos' if sys.platform == 'darwin' else 'linux' if sys.platform.startswith('linux') else sys.platform}-{arch}"


def reference_path(row: str, variant: bool = False) -> Path:
    """<row>.<platform>.json when it exists (last-digit differences of BLAS and compiler builds), else <row>.json."""
    native = REFERENCES / f"{row}.{platform_name()}.json"
    return native if variant or native.exists() else REFERENCES / f"{row}.json"


def run(binary: Path, work: Path, restart: str | None = None) -> int:
    """One solver process in `work`; stdout and stderr go to testrun.log. `restart` is unused (kept for the callers)."""
    with open(work / "testrun.log", "w") as log:
        return subprocess.run([str(binary), "input.nml"], cwd=work, env={**os.environ, **SINGLE_THREAD},
                              stdout=log, stderr=subprocess.STDOUT).returncode


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fresh(src: Path, dest: Path) -> Path:
    shutil.copytree(src, dest)
    for p in [dest, *dest.rglob("*")]:
        p.chmod(p.stat().st_mode | 0o200)
    return dest


def set_nstep(text: str, n: int = NSTEP) -> str:
    return re.sub(r"(?m)^(\s*)nstep\s*=\s*\d+", rf"\g<1>nstep = {n}", text)


def use_kspace(text: str) -> str:
    return re.sub(r"(?m)^&self\s*$", "&self\nuse_kspace = .true.", text, count=1)


def add_u(text: str, u: float, j: float, extra: str = "") -> str:
    return text + f"\n&hamiltonian\nhubbard_u_general(1,3) = {u}\nhubbard_j_general(1,3) = {j}\n{extra}/\n"


def in_hamiltonian(text: str, keys: str) -> str:
    """Insert keys into the existing &hamiltonian group (a second group would be ignored)."""
    return text.replace("&hamiltonian", "&hamiltonian\n" + keys, 1)


def soc_off(text: str) -> str:
    """SOC held at zero: soc_scale alone is undone by the per-iteration xi update, fix_soc keeps it."""
    return text.replace("&self", "&self\nsoc_scale = 0.0\nfix_soc = .true.", 1)


def prepare(src: Path, dest: Path, edit=None) -> Path:
    fresh(src, dest)
    if edit:
        deck = dest / "input.nml"
        deck.write_text(edit(deck.read_text()))
    return dest


def log_text(work: Path) -> str:
    out = ""
    for name in ("testrun.log", "checked_launcher.log"):
        if (work / name).exists():
            out += (work / name).read_text(errors="replace")
    return re.sub(r"\x1b\[[0-9;]*m", "", out)


def numbers(text: str) -> list[float]:
    return [float(x) for x in text.replace("D", "E").replace("d", "e").replace(",", " ").split()]


def signature(work: Path, diag_only: bool) -> dict[str, list[float]]:
    """Per-iteration Diff, Fermi energy, spin moment, routed energy from the log; dm sidecar entries (spin-diagonal only if asked)."""
    text = log_text(work)
    sig = {"diff": numbers(" ".join(re.findall(r"Diff=\s*([0-9.Ee+-]+)", text))),
           "fermi": numbers(" ".join(re.findall(r"(?:Free Fermi energy:|Fixed Fermi energy:|Using canonical/fixed Fermi level =)\s*([0-9.Ee+-]+)", text))),
           "spin_moment": numbers(" ".join(re.findall(r"Spin moment of atom\s+1 is\s+([0-9.Ee+-]+)", text))),
           "energy": numbers(" ".join(re.findall(r"(?<!base_)physical_total_energy=\s*([0-9.Ee+-]+)", text)))}
    for side in sorted(work.glob("*_dm.nml")):
        for line in side.read_text().splitlines():
            m = re.match(r"\s*(dm_(?:re|im)\((\d+),(\d+),(\d+),(\d+),:\))\s*=\s*(.*)", line)
            if m and (not diag_only or m.group(3) == m.group(4)):
                sig[m.group(1).replace("dm_", "dm_" + side.stem + "_")] = numbers(m.group(6))
    return sig


def compare_signatures(a: dict, b: dict, tol: float) -> tuple[bool, str]:
    if set(a) != set(b) or any(len(a[k]) != len(b[k]) for k in a):
        return False, "signature keys or lengths differ"
    worst = max((abs(x - y) for k in a for x, y in zip(a[k], b[k])), default=float("inf"))
    if not a["diff"] or not a["fermi"]:
        return False, "no SCF log values found"
    return worst <= tol, f"{len(a['diff'])} iterations, {len(a)} quantities, max |delta| {worst:.2e} (tol {tol:.0e})"


def ks_host_and_impurity(base: Path, tag: str, cand: Path, nsp: int, soc0: bool, rs_host: bool = False) -> Path | None:
    """Fe host (KS U=2 J=0 local, or the RS bulk +U fixture) at the given nsp, then the B2FeCo impurity deck embedded in it at the same nsp."""
    if rs_host:
        case, restart = ("bulk/bccFe_hubbard" if nsp == 2 else "bulk/bccFe_hubbard_nsp1"), "Fe.nml"
        host_edit = lambda t: in_hamiltonian(set_nstep(t), LOCAL)
    else:
        case, restart = "k_space_scf/bccFe", None
        host_edit = lambda t: add_u(use_kspace(set_nstep(t)), 2.0, 0.0, LOCAL)
        if nsp == 2:
            host_edit = lambda t: add_u(ks_nsp2_deck(t), 2.0, 0.0, LOCAL)
    if nsp == 2 and soc0:
        host_edit = lambda t, f=host_edit: soc_off(f(t))
    host = prepare(CASES / case, base / f"{tag}_host", host_edit)
    if run(cand, host, restart) != 0:
        return None

    def imp_edit(t):
        t = in_hamiltonian(set_nstep(t).replace("nsp = 1", f"nsp = {nsp}"), LOCAL)
        return soc_off(t) if soc0 else t

    d = prepare(CASES / "impurity/B2FeCo_nsp1", base / tag, imp_edit)
    shutil.copy2(host / "Fe_out.nml", d / "Fe.nml")
    shutil.copy2(host / "Fe_dm.nml", d / "Fe_dm.nml")
    return d if run(cand, d, "Fe.nml") == 0 else None


def surface_deck(dest: Path, host: str, nsp: int, soc0: bool) -> Path:
    """Reduced Ni(001) U=3/J=0.7 surface from the committed seeds and the frozen host in <host>/.
    nsp=1 derives its host from the nsp=2 one: spinor blocks dropped, nsp=1, new dm hash (same numbers)."""
    fresh(SURFACE, dest)
    for f in (dest / host).iterdir():
        shutil.copy2(f, dest / f.name)
    deck = dest / "input.nml"
    text = set_nstep(deck.read_text())
    text = re.sub(r"(?m)^nsp\s*=.*$", f"nsp = {nsp}", text)
    text = re.sub(r"(?m)^lld\s*=.*$", "lld = 50", text)
    text = re.sub(r"(?m)^channels_ldos\s*=.*$", "channels_ldos = 500", text)
    text = re.sub(r"(?m)^rc\s*=.*$", "rc = 100\nstrux_backend = 'legacy'", text)
    text = re.sub(r"(?m)^fermi\s*=.*$", "fermi = " + ("-0.076743" if host == "host_soc0" else "-0.076843"), text)
    deck.write_text(soc_off(text) if soc0 else text)
    if nsp == 1:
        dm = dest / "Ni_dm.nml"
        dm.write_text(dm.read_text().replace(" spinor_d = 1\n", ""))
        pot = dest / "Ni.nml"
        pot.write_text(pot.read_text().replace("hubbard_nsp = 2", "hubbard_nsp = 1"))
    return dest


def ks_nsp2_deck(t: str) -> str:
    return use_kspace(set_nstep(t)).replace("nsp = 1", "nsp = 2").replace("&control", "&control\nlmax = 2", 1)


def set_mom(path: Path, axis: tuple[float, float, float]) -> None:
    path.write_text(re.sub(r"(?m)^(\s*mom\s*=).*$", r"\g<1> %.16f, %.16f, %.16f" % axis, path.read_text(), count=1))


def no_symmetry(text: str) -> str:
    """Full k mesh: the spglib reduction uses the crystal group, not the magnetic group of a tilted axis."""
    return text.replace("use_symmetry_reduction = .true.", "use_symmetry_reduction = .false.")


# ------------------------------------------------------------------ rows
def potential_numbers(work: Path) -> dict[str, list[float]]:
    out = {}
    for pot in sorted(work.glob("*_out.nml")):
        if pot.name.endswith("_dm.nml"):
            continue
        text = "\n".join(l.split("!")[0] for l in pot.read_text().splitlines() if "=" in l and "dm_file" not in l)
        vals = re.findall(r"[-+]?\d+\.\d*(?:[eEdD][-+]?\d+)?|[-+]?\d+[eEdD][-+]?\d+", text.split("=", 1)[-1] if False else text)
        out["pot:" + pot.name] = [float(v.replace("d", "e").replace("D", "e")) for v in vals]
    return out


def n(t): return set_nstep(t)
def ks(t): return use_kspace(set_nstep(t))
def ks_u(t): return add_u(ks(t), 2.0, 0.5, LOCAL)
def nsp1(t): return set_nstep(t).replace("nsp = 2", "nsp = 1")
def ks_nsp2(t): return ks_nsp2_deck(t)
def rs_local(t): return in_hamiltonian(set_nstep(t), LOCAL)
def rs_soc0(t): return soc_off(in_hamiltonian(set_nstep(t), LOCAL))
def with_beta(f): return lambda t: f(t).replace("&mix", "&mix\nldm_beta = 0.3", 1)


def deck_signature(binary, base, tag, case, edit, restart=None, diag_only=False, pot=False, mom=None):
    d = prepare(CASES / case, base / tag, edit)
    if mom:
        set_mom(d / "Fe.nml", mom)
    if run(binary, d, restart) != 0:
        raise RuntimeError(f"{tag}: run failed\n{(d / 'testrun.log').read_text(errors='replace')[-1500:]}")
    sig = signature(d, diag_only)
    return {**sig, **potential_numbers(d)} if pot else sig


def surface_signature(binary, base, tag, host, nsp, soc0):
    d = surface_deck(base / tag, host, nsp, soc0)
    if run(binary, d) != 0:
        raise RuntimeError(f"{tag}: run failed\n{log_text(d)[-1500:]}")
    return signature(d, nsp == 1 or soc0)


def impurity_signature(binary, base, tag, nsp, soc0, rs_host):
    d = ks_host_and_impurity(base, tag, binary, nsp, soc0, rs_host=rs_host)
    if d is None:
        raise RuntimeError(f"{tag}: host or impurity run failed")
    return signature(d, nsp == 1 or soc0)


# row -> dict(kind, run(binary, base) -> candidate signature, ref(binary, base) -> reference signature for --update, tolerance, first)
def abs_row(case, edit, restart=None, pot=False, mom=None, diag_only=False, description="", drop=()):
    """`drop`: log quantities left out of the stored reference (the routed energy is not in the old build's comparison)."""
    def f(b, base, tag="abs"):
        sig = deck_signature(b, base, tag, case, edit, restart, diag_only, pot, mom)
        return {k: v for k, v in sig.items() if k not in drop}
    return dict(kind="abs", cand=f, ref=f, tol=ABS_TOL, description=description)


def soc0_row(case1, edit1, case2, edit2, restart=None, tol=BULK_SOC0_TOL, description=""):
    return dict(kind="soc0",
                cand=lambda b, base: deck_signature(b, base, "nsp2_soc0", case2, edit2, restart, True),
                ref=lambda b, base: deck_signature(b, base, "nsp1", case1, edit1, restart, True),
                tol=tol, description=description)


ROWS = {
    "noU_rs": abs_row("bulk/bccFe", n, pot=True, description="no-U RS bcc Fe nsp=2, key outputs", drop=("energy",)),
    "noU_rs_nsp1": abs_row("bulk/bccFe", nsp1, pot=True, description="no-U RS bcc Fe nsp=1", drop=("energy",)),
    "noU_ks": abs_row("k_space_scf/bccFe", ks, pot=True, description="no-U KS bcc Fe nsp=1", drop=("energy",)),
    "noU_ks_nsp2": abs_row("k_space_scf/bccFe", ks_nsp2, pot=True, description="no-U KS bcc Fe nsp=2", drop=("energy",)),
    "U_rs": abs_row("bulk/bccFe_hubbard_nsp1", rs_local, restart="Fe.nml", description="+U RS bulk bcc Fe nsp=1", drop=("energy",)),
    "U_ks": abs_row("k_space_scf/bccFe", ks_u, description="+U KS bcc Fe nsp=1, U=2 J=0.5", drop=("energy",)),
    "U_imp": abs_row("impurity/B2FeCo_nsp1", n, description="+U impurity B2FeCo nsp=1", drop=("energy",)),
    "rs_soc_u2": abs_row("bulk/bccFe_hubbard", lambda t: in_hamiltonian(set_nstep(t), LOCAL), restart="Fe.nml",
                         description="+U RS bulk nsp=2, SOC on, local"),
    "ks_soc_u2": abs_row("k_space_scf/bccFe", lambda t: add_u(ks_nsp2(t), 2.0, 0.0, LOCAL),
                         description="+U KS bulk nsp=2, SOC on, z axis"),
    "ks_soc_u2_x": abs_row("k_space_scf/bccFe", lambda t: no_symmetry(add_u(ks_nsp2(t), 2.0, 0.0, LOCAL)),
                           mom=(1.0, 0.0, 0.0), description="+U KS bulk nsp=2, SOC on, tilted (x) axis"),
    "U_rs_soc0": soc0_row("bulk/bccFe_hubbard_nsp1", rs_local, "bulk/bccFe_hubbard", rs_soc0, "Fe.nml",
                          description="+U RS nsp=2 SOC held at zero against nsp=1"),
    "U_ks_soc0": soc0_row("k_space_scf/bccFe", lambda t: add_u(ks(t), 2.0, 0.0, LOCAL), "k_space_scf/bccFe",
                          lambda t: soc_off(add_u(ks_nsp2(t), 2.0, 0.0, LOCAL)),
                          description="+U KS nsp=2 SOC held at zero against nsp=1"),
    "U_rs_soc0_ldmb": soc0_row("bulk/bccFe_hubbard_nsp1", with_beta(rs_local), "bulk/bccFe_hubbard", with_beta(rs_soc0), "Fe.nml",
                               description="+U RS nsp=2 SOC held at zero against nsp=1, ldm_beta=0.3"),
    "imp_ks_soc0": dict(kind="soc0", tol=IMP_SOC0_TOL, description="+U impurity nsp=2 SOC held at zero against nsp=1, KS host",
                        cand=lambda b, base: impurity_signature(b, base, "imp_nsp2_soc0", 2, True, False),
                        ref=lambda b, base: impurity_signature(b, base, "imp_nsp1", 1, False, False)),
    "imp_rs_soc0": dict(kind="soc0", tol=IMP_SOC0_TOL, description="+U impurity nsp=2 SOC held at zero against nsp=1, RS host",
                        cand=lambda b, base: impurity_signature(b, base, "imp_rs_nsp2_soc0", 2, True, True),
                        ref=lambda b, base: impurity_signature(b, base, "imp_rs_nsp1", 1, False, True)),
    "imp_soc_u2": dict(kind="abs", tol=ABS_TOL, description="+U impurity nsp=2 SOC on, KS host",
                       cand=lambda b, base: impurity_signature(b, base, "imp_soc_u2", 2, False, False),
                       ref=lambda b, base: impurity_signature(b, base, "imp_soc_u2", 2, False, False)),
    "surf_soc0": dict(kind="soc0", tol=SURF_SOC0_TOL, description="+U Ni(001) surface nsp=2 SOC held at zero against nsp=1",
                      cand=lambda b, base: surface_signature(b, base, "surf_nsp2_soc0", "host_soc0", 2, True),
                      ref=lambda b, base: surface_signature(b, base, "surf_nsp1", "host_soc0", 1, False)),
    "surf_soc_u3": dict(kind="abs", tol=ABS_TOL, description="+U Ni(001) surface nsp=2 SOC on",
                        cand=lambda b, base: surface_signature(b, base, "surf_soc_u3", "host_soc", 2, False),
                        ref=lambda b, base: surface_signature(b, base, "surf_soc_u3", "host_soc", 2, False)),
}


def compare(row: str, ref: dict, cand: dict) -> tuple[bool, str]:
    spec = ROWS[row]
    # A reference from the old build 90f57c3 predates the dm sidecars: it holds the log quantities only.
    if not set(ref) <= set(cand) or any(len(ref[k]) != len(cand[k]) for k in ref):
        return False, f"quantities differ: reference-only {sorted(set(ref) - set(cand))[:3]}"
    cand = {k: cand[k] for k in ref}
    if not ref.get("diff") or not ref.get("fermi"):
        return False, "no SCF log values found"
    worst = max((abs(a - b) for k in ref for a, b in zip(ref[k], cand[k])), default=float("inf"))
    ok = worst <= spec["tol"]
    detail = f"{len(ref['diff'])} iterations, {len(ref)} quantities, max |delta| {worst:.2e} (tol {spec['tol']:.0e})"
    if spec["kind"] == "soc0":
        first = max(abs(ref[k][0] - cand[k][0]) for k in ("diff", "fermi", "spin_moment", "energy") if ref.get(k))
        offdiag = max((abs(x) for k, v in cand.items() if (m := re.search(r"dm_\S*\((\d+),(\d),(\d),", k)) and m[2] != m[3] for x in v), default=0.0)
        ok = ok and first <= SOC0_TOL and offdiag <= OFFDIAG_TOL
        detail += f"; iteration-1 max |delta| {first:.2e} (tol {SOC0_TOL:.0e}); max |P_up,dn| {offdiag:.1e} (tol {OFFDIAG_TOL:.0e})"
    return ok, detail


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("row", nargs="?")
    ap.add_argument("--binary", type=Path)
    ap.add_argument("--scratch-root", type=Path)
    ap.add_argument("--update", metavar="REASON", help="rewrite the stored reference from this binary; REASON is recorded")
    ap.add_argument("--variant", action="store_true", help="with --update: write the reference of this platform (<row>.<platform>.json)")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for name, spec in ROWS.items():
            print(f"{name:16s} {spec['kind']:5s} {reference_path(name).exists()!s:5s} {spec['description']}")
        return 0
    if args.row not in ROWS or args.binary is None:
        ap.error("ROW (one of --list) and --binary are required")
    binary = args.binary.resolve()
    base = (args.scratch_root or Path(tempfile.mkdtemp(prefix="ldau_row_"))).resolve()
    if base.exists() and args.scratch_root:
        shutil.rmtree(base)
    base.mkdir(parents=True, exist_ok=True)
    spec = ROWS[args.row]
    path = reference_path(args.row, variant=args.update is not None and args.variant)
    if args.update is not None and not args.variant:
        path = REFERENCES / f"{args.row}.json"
    if args.update is not None:
        sig = spec["ref"](binary, base)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"row": args.row, "description": spec["description"], "kind": spec["kind"],
                                    "tolerance": spec["tol"], "reason": args.update,
                                    "binary_sha256": sha256(binary), "date": datetime.date.today().isoformat(),
                                    "reference": sig}, indent=1) + "\n")
        print(f"{args.row}: reference written to {path.relative_to(ROOT)} ({len(sig)} quantities)")
        return 0
    try:
        if spec["kind"] == "soc0":
            # A relation between two runs of this build: the nsp=1 run is made now, so that build-to-build
            # differences do not enter. The JSON of such a row is the archived nsp=1 run of the reference build.
            stored = spec["ref"](binary, base)
        else:
            stored = json.loads(path.read_text())["reference"]
        cand = spec["cand"](binary, base)
    except OSError:
        print(f"FAIL {args.row}: no reference {path.relative_to(ROOT)}")
        return 1
    except RuntimeError as err:
        print(f"FAIL {args.row}: {err}")
        return 1
    ok, detail = compare(args.row, stored, cand)
    print(f"{'PASS' if ok else 'FAIL'} {args.row}: {detail}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
