#!/usr/bin/env python3
"""One-command LDA+U regression: reference binary versus candidate binary.

Run with the repository .venv interpreter:
    .venv/bin/python tools/ldau_regression.py --ref REF/rslmto.x --cand build/bin/rslmto.x

Every deck is copied fresh from tests/scf/cases into a new run directory and
launched single-threaded through tools/run_rslmto_checked.py (4 SCF steps).

The guard-fixture script (RA2) runs against the candidate binary only: the
reference binary predates the guards it checks, so a reference row cannot pass.

The SOC-on nsp=2 +U rows (z axis, and an x-axis KS row since S6b) compare the candidate with a committed snapshot
(tests/validation/fixed_ldau/data/soc_u_snapshot.txt) at 1e-10. Rewrite the
snapshot only with --update-snapshot "one-line reason".
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / "tests/scf/cases"
WRAPPER = ROOT / "tools/run_rslmto_checked.py"
RA2 = ROOT / "tests/validation/fixed_ldau/rs_embedding_ra2.py"
S6_GUARDS = ROOT / "tests/validation/fixed_ldau/s6_nsp2_impurity_guards.py"
SNAPSHOT = ROOT / "tests/validation/fixed_ldau/data/soc_u_snapshot.txt"
SNAPSHOT_TOL = 1e-10
SOC0_TOL = 1e-9  # the logged Diff has 10 decimals
AXIS_TOL = 1e-8  # S6b: x-axis against z-axis at SOC held at zero, after spin rotation
IMP_KS_SOC0_TOL = 1e-7  # the two impurity soc0 rows only, reasons next to their use
SURFACE = ROOT / "tests/validation/fixed_ldau/decks/ni001_soc_u3"
SURF_SOC0_TOL = 1e-7  # surface soc0 row, 4 steps (iteration 1 uses SOC0_TOL)
PY = str(ROOT / ".venv/bin/python")
NSTEP = 4
NOT_SCIENTIFIC = {"preflight.json", "checked_launcher.log", "testrun.log", "run_monitor.json"}
DM_LINE = re.compile(r"^\s*(ldm_flatten|ldm_flatten_imag|dm_file|dm_sha256|hubbard_nsp)\b", re.M)

rows: list[dict] = []


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(name: str, ok: bool, detail: str, gate: bool = True) -> None:
    rows.append({"check": name, "result": "PASS" if ok else "FAIL", "detail": detail, "gate": gate})


# ---------------------------------------------------------------- running
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


LOCAL = "hubbard_radial_moments = 'local'\n"


def in_hamiltonian(text: str, keys: str) -> str:
    """Insert keys into the existing &hamiltonian group (a second group would be ignored)."""
    return text.replace("&hamiltonian", "&hamiltonian\n" + keys, 1)


def soc_off(text: str) -> str:
    """SOC held at zero: soc_scale alone is undone by the per-iteration xi update, fix_soc keeps it."""
    return text.replace("&self", "&self\nsoc_scale = 0.0\nfix_soc = .true.", 1)


def run(binary: Path, work: Path, restart: str | None = None) -> int:
    cmd = [PY, str(WRAPPER), str(work), "--binary", str(binary)]
    if restart:
        cmd += ["--restart", str(work / restart)]
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True).returncode


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


# ------------------------------------------------------------- comparison
def outputs(work: Path) -> dict[str, Path]:
    return {str(p.relative_to(work)): p for p in sorted(work.rglob("*"))
            if p.is_file() and p.name not in NOT_SCIENTIFIC and not p.name.startswith(".")}


def normalized(p: Path) -> bytes:
    data = p.read_bytes()
    if p.name == "str.out":
        data = re.sub(rb"(?m)(cpu_s=\s*)[0-9.]+", rb"\1<RUNTIME>", data)
    return data


def strip_dm(path: Path) -> bytes:
    return "".join(l for l in path.read_text().splitlines(keepends=True) if not DM_LINE.match(l)).encode()


def numbers(text: str) -> list[float]:
    return [float(x) for x in text.replace("D", "E").replace("d", "e").replace(",", " ").split()]


def legacy_blocks(path: Path) -> dict:
    """ldm_flatten(l,s,:) rows of a legacy inline potential: {(part,l,s): [values]}."""
    blocks = {}
    for line in path.read_text().splitlines():
        m = re.match(r"\s*(ldm_flatten(?:_imag)?)\((\d+),\s*(\d+),\s*:\)\s*=\s*(.*)", line)
        if m:
            blocks[("im" if m.group(1).endswith("imag") else "re", int(m.group(2)), int(m.group(3)))] = numbers(m.group(4))
    return blocks


def sidecar_blocks(path: Path) -> dict:
    """dm_re/dm_im(l,sr,sc,m1,:) of a dm sidecar, flattened row-major per (part,l,spin) for sr==sc."""
    rowsd: dict = {}
    for line in path.read_text().splitlines():
        m = re.match(r"\s*dm_(re|im)\((\d+),(\d+),(\d+),(\d+),:\)\s*=\s*(.*)", line)
        if m:
            l, sr, sc, m1 = map(int, m.group(2, 3, 4, 5))
            if sr == sc:
                rowsd.setdefault((m.group(1), l, sr), {})[m1] = numbers(m.group(6))
    return {k: [v for m1 in sorted(r) for v in r[m1]] for k, r in rowsd.items()}


def compare_by_value(ref: Path, cand: Path) -> tuple[bool, str]:
    """Reference and candidate runs: all non-occupation bytes identical, occupations equal by value."""
    ro, co = outputs(ref), outputs(cand)
    problems = []
    extra = sorted(set(co) - set(ro))
    missing = sorted(set(ro) - set(co))
    sidecars = [e for e in extra if e.endswith("_dm.nml")]
    if missing or [e for e in extra if e not in sidecars]:
        problems.append(f"file sets differ: missing={missing} extra={[e for e in extra if e not in sidecars]}")
    n_nml = n_vals = 0
    for name in sorted(set(ro) & set(co)):
        r, c = ro[name], co[name]
        if name.endswith(".nml") and (DM_LINE.search(r.read_text()) or DM_LINE.search(c.read_text())):
            n_nml += 1
            if strip_dm(r) != strip_dm(c):
                problems.append(f"{name}: non-occupation lines differ")
            m = re.search(r"(?m)^\s*dm_file\s*=\s*'([^']*)'", c.read_text())
            if m:
                side = cand / m.group(1)
                digest = re.search(r"(?m)^\s*dm_sha256\s*=\s*'([^']*)'", c.read_text())
                if not side.is_file() or digest is None or digest.group(1) != sha256(side):
                    problems.append(f"{name}: dm_sha256 does not match {m.group(1)}")
                    continue
                leg, new = legacy_blocks(r), sidecar_blocks(side)
                for key, vals in leg.items():
                    side_vals = new.get(key, [])
                    # legacy rows are zero-padded to 25 entries
                    if key[2] <= 2 and (not side_vals or vals[:len(side_vals)] != side_vals or any(vals[len(side_vals):])):
                        problems.append(f"{name}: occupation block {key} differs")
                    n_vals += len(vals)
                if not leg:
                    problems.append(f"{name}: reference has no inline occupations")
        elif normalized(r) != normalized(c):
            problems.append(f"{name}: differs")
    if not n_nml:
        problems.append("no potential file with occupations was compared")
    return (not problems, "; ".join(problems[:4]) if problems else
            f"{len(set(ro) & set(co))} files, {n_nml} potential files, {n_vals} occupation values equal")


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


def soc0_vs_nsp1(name, tag, case1, edit1, case2, edit2, base, cand, restart=None):
    """Candidate nsp=2 with SOC held at zero against the candidate nsp=1 path, same U."""
    dirs = []
    for sub, case, edit in (("nsp1", case1, edit1), ("nsp2_soc0", case2, edit2)):
        d = prepare(CASES / case, base / f"{tag}_{sub}", edit)
        rc = run(cand, d, restart)
        if rc != 0:
            record(name, False, f"{sub} run failed rc={rc}")
            return
        dirs.append(d)
    ok, detail = compare_signatures(signature(dirs[0], True), signature(dirs[1], True), SOC0_TOL)
    offdiag = max((abs(x) for k, v in signature(dirs[1], False).items() if re.match(r"dm_\w+\(3,([12]),([12]),", k)
                   and k.split("(3,")[1][0] != k.split("(3,")[1][2] for x in v), default=0.0)
    record(name, ok and offdiag <= 1e-12, f"{detail}; max |P_up,dn| {offdiag:.1e} (tol 1e-12)")


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
    text = re.sub(r"(?m)^fermi\s*=.*$", "fermi = " + ("-0.076743" if host == "host_soc0" else "-0.076576"), text)
    deck.write_text(soc_off(text) if soc0 else text)
    if nsp == 1:
        dm = dest / "Ni_dm.nml"
        dm.write_text(dm.read_text().replace(" spinor_d = 1\n", ""))
        pot = dest / "Ni.nml"
        pot.write_text(re.sub(r"dm_sha256 = '[0-9a-f]+'", f"dm_sha256 = '{sha256(dm)}'",
                              pot.read_text().replace("hubbard_nsp = 2", "hubbard_nsp = 1")))
    return dest


def surface_rows(base: Path, cand: Path) -> dict[str, Path]:
    """nsp=2 surface: SOC held at zero against nsp=1 (same host numbers); returns the SOC-on run for the snapshot."""
    name = "+U surface nsp=2 soc_scale=0 fix_soc versus nsp=1, local"
    runs = {}
    for tag, host, nsp, soc0 in (("surf_nsp1", "host_soc0", 1, False), ("surf_nsp2_soc0", "host_soc0", 2, True),
                                 ("surf_soc_u3", "host_soc", 2, False)):
        d = surface_deck(base / tag, host, nsp, soc0)
        if run(cand, d) != 0:
            record(name if tag != "surf_soc_u3" else "+U surface nsp=2 SOC on versus snapshot", False, f"{tag} run failed")
            return {}
        runs[tag] = d
    s1, s0 = signature(runs["surf_nsp1"], True), signature(runs["surf_nsp2_soc0"], True)
    ok, detail = compare_signatures(s1, s0, SURF_SOC0_TOL)
    first = max(abs(s1[k][0] - s0[k][0]) for k in ("diff", "fermi", "spin_moment", "energy"))
    offdiag = max((abs(x) for k, v in signature(runs["surf_nsp2_soc0"], False).items()
                   if (m := re.match(r"dm_Ni-S(?:-1)?_dm_(?:re|im)\(3,(\d),(\d),", k)) and m[1] != m[2] for x in v), default=0.0)
    record(name, ok and first <= SOC0_TOL and offdiag <= 1e-12,
           f"{detail}; iteration-1 max |delta| {first:.2e} (tol {SOC0_TOL:.0e}); max |P_up,dn| {offdiag:.1e} (tol 1e-12)")
    return {"surf_soc_u3": runs["surf_soc_u3"]}


def ks_nsp2_deck(t: str) -> str:
    return use_kspace(set_nstep(t)).replace("nsp = 1", "nsp = 2").replace("&control", "&control\nlmax = 2", 1)


def snapshot_rows(runs: dict[str, Path]) -> dict[str, dict]:
    return {tag: signature(d, False) for tag, d in runs.items()}


def read_snapshot() -> dict[str, dict]:
    snap: dict[str, dict] = {}
    for line in SNAPSHOT.read_text().splitlines():
        if line and not line.startswith("#"):
            tag, key, *vals = line.split()
            snap.setdefault(tag, {})[key] = [float(v) for v in vals]
    return snap


def write_snapshot(rows_: dict[str, dict], reason: str, cand: Path) -> None:
    lines = ["# SOC-on nsp=2 +U regression snapshot (bcc Fe, U=2 J=0, local, 4 steps; imp_soc_u2: B2FeCo impurity in the KS Fe host; surf_soc_u3: reduced Ni(001) U=3 J=0.7 surface); columns: case key values",
             f"# reason: {reason}", f"# candidate sha256: {sha256(cand)}"]
    for tag, sig in rows_.items():
        lines += [f"{tag} {key} " + " ".join(f"{v:.17e}" for v in vals) for key, vals in sig.items()]
    SNAPSHOT.write_text("\n".join(lines) + "\n")


def compare_identical(ref: Path, cand: Path, skip: tuple[str, ...] = ()) -> tuple[bool, str]:
    ro = {k: v for k, v in outputs(ref).items() if k not in skip}
    co = {k: v for k, v in outputs(cand).items() if k not in skip}
    bad = sorted(set(ro) ^ set(co)) + [n for n in sorted(set(ro) & set(co)) if normalized(ro[n]) != normalized(co[n])]
    return (not bad, f"{len(set(ro) & set(co))} files byte-identical" if not bad else "differ: " + ", ".join(bad[:5]))


def set_mom(path: Path, axis: tuple[float, float, float]) -> None:
    path.write_text(re.sub(r"(?m)^(\s*mom\s*=).*$", r"\g<1> %.16f, %.16f, %.16f" % axis, path.read_text(), count=1))


def no_symmetry(text: str) -> str:
    """Full k mesh: the spglib reduction uses the crystal group, not the magnetic group of a tilted axis."""
    return text.replace("use_symmetry_reduction = .true.", "use_symmetry_reduction = .false.")


def spinor_p(work: Path) -> np.ndarray:
    """Full 10x10 d spinor P of the Fe dm sidecar, index (s-1)*5+m."""
    p = np.zeros((10, 10), complex)
    for line in (work / "Fe_dm.nml").read_text().splitlines():
        m = re.match(r"\s*dm_(re|im)\(3,(\d),(\d),(\d),:\)\s*=\s*(.*)", line)
        if m:
            s, t, m1 = int(m.group(2)), int(m.group(3)), int(m.group(4))
            p[(s - 1) * 5 + m1 - 1, (t - 1) * 5:t * 5] += np.array(numbers(m.group(5))) * (1 if m.group(1) == "re" else 1j)
    return p


def to_z_frame(p: np.ndarray, n: tuple[float, float, float]) -> np.ndarray:
    """U P U^dagger with U = rows <n+|, <n-| (|n+> = (cos t/2, e^{i f} sin t/2)) x 1_5."""
    th, ph = np.arccos(n[2]), np.arctan2(n[1], n[0])
    r2 = np.array([[np.cos(th / 2), np.exp(-1j * ph) * np.sin(th / 2)],
                   [-np.exp(1j * ph) * np.sin(th / 2), np.cos(th / 2)]])
    u = np.kron(r2, np.eye(5))
    return u @ p @ u.conj().T


def axis_rotation_row(base: Path, cand: Path, edit) -> None:
    """KS bcc Fe nsp=2, U=2 J=0 local, SOC held at zero: magnetisation along x against z, by value after rotation."""
    name = "+U KS nsp=2 soc_scale=0 fix_soc, axis x versus z (rotated)"
    dirs = {}
    for tag, axis in (("z", (0.0, 0.0, 1.0)), ("x", (1.0, 0.0, 0.0))):
        d = prepare(CASES / "k_space_scf/bccFe", base / f"axis_soc0_{tag}", edit)
        set_mom(d / "Fe.nml", axis)
        if run(cand, d) != 0:
            record(name, False, f"{tag} run failed")
            return
        dirs[tag] = (d, axis)
    sz, sx = ({k: v for k, v in signature(d, False).items() if not k.startswith("dm_")} for d, _ in dirs.values())
    ok, detail = compare_signatures(sz, sx, AXIS_TOL)
    dp = np.abs(to_z_frame(spinor_p(dirs["x"][0]), dirs["x"][1]) - spinor_p(dirs["z"][0])).max()
    record(name, ok and dp <= AXIS_TOL, f"{detail}; max |U P_x U^+ - P_z| {dp:.2e} (tol {AXIS_TOL:.0e})")


# --------------------------------------------------------------- the rows
def paired(name, tag, case, edit, comparator, base, ref, cand, restart=None):
    dirs = {}
    for label, binary in (("ref", ref), ("cand", cand)):
        d = prepare(CASES / case, base / f"{tag}_{label}", edit)
        rc = run(binary, d, restart)
        dirs[label] = d
        if rc != 0:
            record(name, False, f"{label} run failed rc={rc}")
            return None
    ok, detail = comparator(dirs["ref"], dirs["cand"])
    record(name, ok, detail)
    return dirs


def dm_runtime(base: Path, ref_dir: Path, cand_dir: Path, cand: Path) -> None:
    deck = CASES / "bulk/bccFe_hubbard_nsp1"
    legacy, backed = ref_dir / "Fe_out.nml", cand_dir / "Fe_out.nml"
    sidecar = cand_dir / "Fe_dm.nml"

    def restart_run(tag, potential, sidecar_file=None, edit=None):
        d = prepare(deck, base / tag, lambda t: set_nstep(t, 1))
        shutil.copy2(potential, d / "Fe.nml")
        if sidecar_file:
            shutil.copy2(sidecar_file, d / "Fe_dm.nml")
        if edit:
            edit(d)
        return d, run(cand, d, "Fe.nml")

    d_leg, rc1 = restart_run("dm_legacy", legacy)
    d_dm, rc2 = restart_run("dm_backed", backed, sidecar)
    if rc1 or rc2:
        record("dm: legacy inline vs dm-backed restart", False, f"run failed rc={rc1},{rc2}")
    else:
        ok, detail = compare_identical(d_leg, d_dm, skip=("Fe.nml",))
        record("dm: legacy inline vs dm-backed restart, all outputs", ok, detail)

    def perturb(d):
        text = (d / "Fe.nml").read_text()
        text = re.sub(r"(?m)^( ldm_flatten\(3, 1, :\) = )([0-9.eE+-]+)",
                      lambda m: m.group(1) + repr(float(m.group(2)) + 1e-3), text, count=1)
        (d / "Fe.nml").write_text(text)

    d_pert, rc3 = restart_run("dm_legacy_perturbed", legacy, edit=perturb)
    if rc3:
        record("dm: control, perturbed legacy occupation changes output", False, f"run failed rc={rc3}")
    else:
        same, _ = compare_identical(d_leg, d_pert, skip=("Fe.nml",))
        record("dm: control, perturbed legacy occupation changes output", not same,
               "outputs differ as required" if not same else "outputs identical: restart occupations have no effect")

    d_miss, rc4 = restart_run("dm_missing", backed)
    record("dm: missing sidecar is fatal", rc4 != 0 and "dm" in log_text(d_miss).lower(),
           f"rc={rc4}; " + next((l.strip()[:90] for l in log_text(d_miss).splitlines() if "dm" in l.lower()), "no dm message"))

    def corrupt(d):
        t = (d / "Fe_dm.nml").read_text()
        (d / "Fe_dm.nml").write_text(re.sub(r"(dm_re\(3,1,1,1,:\)\s*=\s*)(\S+)",
                                            lambda m: m.group(1) + "%.17E" % (float(m.group(2)) + 1e-6), t, count=1))

    d_hash, rc5 = restart_run("dm_hash_mismatch", backed, sidecar, edit=corrupt)
    record("dm: hash mismatch is fatal", rc5 != 0 and "sha" in log_text(d_hash).lower(),
           f"rc={rc5}; " + next((l.strip()[:90] for l in log_text(d_hash).splitlines() if "sha" in l.lower()), "no hash message"))


def ctest(name: str, build: Path, selector: list[str]) -> None:
    shutil.rmtree(build / "Testing/validation", ignore_errors=True)  # stale scratch from earlier ctest runs
    p = subprocess.run(["ctest", "--test-dir", str(build), *selector], capture_output=True, text=True)
    m = re.search(r"(\d+)% tests passed, (\d+) tests failed out of (\d+)", p.stdout)
    record(name, p.returncode == 0 and bool(m), f"{m.group(3)} tests, {m.group(2)} failed" if m else "ctest failed")


def guard_fixtures(binary: Path, base: Path, label: str) -> None:
    p = subprocess.run([PY, str(RA2), "--binary", str(binary), "--scratch-root", str(base / f"ra2_{label}")],
                       cwd=ROOT, capture_output=True, text=True)
    detail = "13 scenarios" if p.returncode == 0 else (p.stderr.strip().splitlines() or ["failed"])[-1][:160]
    record(f"guard fixtures (RA2), {label} binary", p.returncode == 0, detail)
    p = subprocess.run([PY, str(S6_GUARDS), "--binary", str(binary), "--scratch-root", str(base / f"s6_{label}")],
                       cwd=ROOT, capture_output=True, text=True)
    detail = f"{p.stdout.count('PASS ')} scenarios" if p.returncode == 0 else (p.stdout.strip().splitlines() or ["failed"])[-1][:160]
    record(f"guard fixtures (S6 nsp=2 impurity), {label} binary", p.returncode == 0, detail)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref", required=True, type=Path)
    ap.add_argument("--cand", required=True, type=Path)
    ap.add_argument("--cand-build", type=Path, default=ROOT / "build", help="build directory for ctest")
    ap.add_argument("--work", type=Path, help="scratch root (default: new temporary directory)")
    ap.add_argument("--json", type=Path, help="write a JSON sidecar here")
    ap.add_argument("--update-snapshot", metavar="REASON", help="rewrite the SOC-on +U snapshot; REASON is recorded in its header")
    args = ap.parse_args()
    ref, cand = args.ref.resolve(), args.cand.resolve()
    base = (args.work or Path(tempfile.mkdtemp(prefix="ldau_regression_"))).resolve()
    base.mkdir(parents=True, exist_ok=True)

    n = lambda t: set_nstep(t)
    ks = lambda t: use_kspace(set_nstep(t))
    ks_u = lambda t: add_u(ks(t), 2.0, 0.5)
    imp = lambda t: set_nstep(t)

    nsp1 = lambda t: set_nstep(t).replace("nsp = 2", "nsp = 1")
    ks_nsp2 = lambda t: use_kspace(set_nstep(t)).replace("nsp = 1", "nsp = 2").replace("&control", "&control\nlmax = 2", 1)
    paired("no-U RS bulk/bccFe (nsp=2), byte-identical", "noU_rs", "bulk/bccFe", n, compare_identical, base, ref, cand)
    paired("no-U RS bulk/bccFe, nsp=1, byte-identical", "noU_rs_nsp1", "bulk/bccFe", nsp1, compare_identical, base, ref, cand)
    paired("no-U KS k_space_scf/bccFe, nsp=2, byte-identical", "noU_ks_nsp2", "k_space_scf/bccFe", ks_nsp2,
           compare_identical, base, ref, cand)
    paired("no-U KS k_space_scf/bccFe, byte-identical", "noU_ks", "k_space_scf/bccFe", ks, compare_identical, base, ref, cand)
    rs_u = paired("+U RS bulk/bccFe_hubbard_nsp1, by value", "U_rs", "bulk/bccFe_hubbard_nsp1", n,
                  compare_by_value, base, ref, cand, restart="Fe.nml")
    paired("+U KS bccFe U=2 J=0.5, by value", "U_ks", "k_space_scf/bccFe", ks_u, compare_by_value, base, ref, cand)
    paired("+U impurity/B2FeCo_nsp1, by value", "U_imp", "impurity/B2FeCo_nsp1", imp, compare_by_value, base, ref, cand)

    # nsp=2 +U (flag on, local): SOC held at zero against nsp=1, then SOC on against the committed snapshot
    rs_soc0 = lambda t: soc_off(in_hamiltonian(set_nstep(t), LOCAL))
    rs_soc = lambda t: in_hamiltonian(set_nstep(t), LOCAL)
    soc0_vs_nsp1("+U RS nsp=2 soc_scale=0 fix_soc versus nsp=1, local", "U_rs_soc0", "bulk/bccFe_hubbard_nsp1",
                 lambda t: in_hamiltonian(set_nstep(t), LOCAL), "bulk/bccFe_hubbard", rs_soc0, base, cand, restart="Fe.nml")
    soc0_vs_nsp1("+U KS nsp=2 soc_scale=0 fix_soc versus nsp=1, local", "U_ks_soc0", "k_space_scf/bccFe",
                 lambda t: add_u(ks(t), 2.0, 0.0, LOCAL), "k_space_scf/bccFe",
                 lambda t: soc_off(add_u(ks_nsp2(t), 2.0, 0.0, LOCAL)), base, cand)
    with_beta = lambda f: (lambda t: f(t).replace("&mix", "&mix\nldm_beta = 0.3", 1))
    soc0_vs_nsp1("+U RS nsp=2 soc_scale=0 fix_soc versus nsp=1, local, ldm_beta=0.3", "U_rs_soc0_ldmb", "bulk/bccFe_hubbard_nsp1",
                 with_beta(lambda t: in_hamiltonian(set_nstep(t), LOCAL)), "bulk/bccFe_hubbard", with_beta(rs_soc0),
                 base, cand, restart="Fe.nml")
    # S6b: magnetisation axis x against z with SOC held at zero (KS), and an x-axis SOC-on snapshot row
    axis_rotation_row(base, cand, lambda t: no_symmetry(soc_off(add_u(ks_nsp2(t), 2.0, 0.0, LOCAL))))
    soc_runs = {}
    for tag, case, edit, restart in (("rs_soc_u2", "bulk/bccFe_hubbard", rs_soc, "Fe.nml"),
                                     ("ks_soc_u2", "k_space_scf/bccFe",
                                      lambda t: add_u(ks_nsp2(t), 2.0, 0.0, LOCAL), None),
                                     ("ks_soc_u2_x", "k_space_scf/bccFe",
                                      lambda t: no_symmetry(add_u(ks_nsp2(t), 2.0, 0.0, LOCAL)), None)):
        d = prepare(CASES / case, base / tag, edit)
        if tag == "ks_soc_u2_x":
            set_mom(d / "Fe.nml", (1.0, 0.0, 0.0))
        if run(cand, d, restart) == 0:
            soc_runs[tag] = d
        else:
            record(f"+U nsp=2 SOC on {tag.upper()} versus snapshot", False, "run failed")
    # nsp=2 impurity (KS Fe host, spinor dm): SOC held at zero against nsp=1, and SOC on into the snapshot
    imp1 = ks_host_and_impurity(base, "imp_nsp1", cand, 1, False)
    imp0 = ks_host_and_impurity(base, "imp_nsp2_soc0", cand, 2, True)
    name0 = "+U impurity nsp=2 soc_scale=0 fix_soc versus nsp=1, local, KS host"
    if imp1 and imp0:
        s1, s0 = signature(imp1, True), signature(imp0, True)
        # KS nsp=2 vs nsp=1 hosts differ at ~1e-10 (equivalent code paths); local-mode +U transient amplification
        # over 4 steps (documented S4 / KNOWN_ISSUES 2026-10-26) gives ~4e-8. Defects this row targets (missing
        # coherences, site map, unmixed occupation) produce 1e-2 to 1e-1 differences, so 1e-7 still discriminates.
        ok, detail = compare_signatures(s1, s0, IMP_KS_SOC0_TOL)
        first = max(abs(s1[k][0] - s0[k][0]) for k in ("diff", "fermi", "spin_moment", "energy"))
        record(name0, ok and first <= SOC0_TOL, f"{detail}; iteration-1 max |delta| {first:.2e} (tol {SOC0_TOL:.0e})")
    else:
        record(name0, False, "run failed")
    nameR = "+U impurity nsp=2 soc_scale=0 fix_soc versus nsp=1, local, RS host"
    imp1r = ks_host_and_impurity(base, "imp_rs_nsp1", cand, 1, False, rs_host=True)
    imp0r = ks_host_and_impurity(base, "imp_rs_nsp2_soc0", cand, 2, True, rs_host=True)
    if imp1r and imp0r:
        # Impurity local-mode +U transient amplifies roundoff-level differences between the nsp=2 (10x10) and
        # nsp=1 (5x5) paths over 4 steps (measured 7.15e-9 RS, ~4e-8 KS). Targeted defects produce 1e-2 to 1e-1.
        r1, r0 = signature(imp1r, True), signature(imp0r, True)
        ok, detail = compare_signatures(r1, r0, IMP_KS_SOC0_TOL)
        first = max(abs(r1[k][0] - r0[k][0]) for k in ("diff", "fermi", "spin_moment", "energy"))
        record(nameR, ok and first <= SOC0_TOL, f"{detail}; iteration-1 max |delta| {first:.2e} (tol {SOC0_TOL:.0e})")
    else:
        record(nameR, False, "run failed")
    imp_soc = ks_host_and_impurity(base, "imp_soc_u2", cand, 2, False)
    if imp_soc:
        soc_runs["imp_soc_u2"] = imp_soc
    else:
        record("+U impurity nsp=2 SOC on versus snapshot", False, "run failed")
    soc_runs.update(surface_rows(base, cand))
    if args.update_snapshot is not None and len(soc_runs) == 5:
        write_snapshot(snapshot_rows(soc_runs), args.update_snapshot, cand)
    for tag, d in list(soc_runs.items()):
        name = {"imp_soc_u2": "+U impurity nsp=2 SOC on versus snapshot",
                "surf_soc_u3": "+U surface nsp=2 SOC on versus snapshot"}.get(tag, f"+U nsp=2 SOC on {tag.upper()} versus snapshot")
        if not SNAPSHOT.is_file():
            record(name, False, "no snapshot; create it with --update-snapshot")
            continue
        snap = read_snapshot().get(tag, {})
        ok, detail = compare_signatures(snap, signature(d, False), SNAPSHOT_TOL) if snap else (False, "case not in snapshot")
        record(name, ok, detail)

    guard_fixtures(cand, base, "candidate")
    ctest("CTest label lda_u", args.cand_build, ["-L", "lda_u"])
    ctest("impurity/surface auto fixtures (UnitAutoImpurityLocalMode)", args.cand_build, ["-R", "UnitAutoImpurityLocalMode"])
    if rs_u:
        dm_runtime(base, rs_u["ref"], rs_u["cand"], cand)
    else:
        record("dm runtime cases", False, "needs the +U RS rows")

    width = max(len(r["check"]) for r in rows)
    print(f"{'check':<{width}}  result  detail")
    for r in rows:
        print(f"{r['check']:<{width}}  {r['result']:<6}  {r['detail']}" + ("" if r["gate"] else "  [info only]"))
    gate_ok = all(r["result"] == "PASS" for r in rows if r["gate"])
    print("OVERALL", "PASS" if gate_ok else "FAIL", f"(work: {base})")
    if args.json:
        args.json.write_text(json.dumps({
            "reference": {"binary": str(ref), "sha256": sha256(ref)},
            "candidate": {"binary": str(cand), "sha256": sha256(cand)},
            "work": str(base), "rows": rows, "overall": "PASS" if gate_ok else "FAIL"}, indent=2) + "\n")
    return 0 if gate_ok else 1


if __name__ == "__main__":
    sys.exit(main())
