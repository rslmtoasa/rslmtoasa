#!/usr/bin/env python3
"""Exercise fixed-U potential provenance, mismatch and legacy behavior."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import signal
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / "tests/scf/cases/k_space_scf/bccFe"
IMPURITY_FIXTURE = ROOT / "tests/scf/cases/impurity/B2FeCo_nsp1"
SURFACE_FIXTURE = ROOT / "tests/scf/cases/surface/fccCu001"
LAUNCHER = ROOT / "tests/validation/fixed_ldau/c5_7_launch.sh"
PREFLIGHT = ROOT / "docs/validation/fixed_ldau/rs/C5_4_PREFLIGHT.py"
THREAD_CAPS = {
    "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE", "RSLMTO_OMP_THREADS_SERIAL": "1",
    "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
    "BLIS_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
}


def configure(work: Path, *, backend: str, mode: str | None, uj: bool = True) -> None:
    path = work / "input.nml"
    text = path.read_text()
    text = re.sub(r"(?m)^nstep\s*=\s*\d+", "nstep = 1", text)
    if re.search(r"(?m)^\s*use_kspace\s*=", text):
        text = re.sub(r"(?m)^\s*use_kspace\s*=.*$", f"use_kspace = {'.true.' if backend == 'ks' else '.false.'}", text)
    else:
        text = re.sub(r"(?m)^&self\s*$", "&self\nuse_kspace = " + (".true." if backend == "ks" else ".false."), text, count=1)
    radial = f"hubbard_radial_moments = '{mode}'\n" if mode else ""
    hamiltonian = "&hamiltonian\n" + radial
    if uj:
        hamiltonian += "hubbard_u_general(1,3) = 2.0\nhubbard_j_general(1,3) = 0.5\n"
    hamiltonian += "/\n"
    path.write_text(text + "\n" + hamiltonian)


def run(binary: Path, work: Path, *, restart_name: str = "Fe.nml") -> tuple[int, str]:
    env = os.environ.copy()
    env.update(THREAD_CAPS)
    pre = subprocess.run(
        [sys.executable, str(PREFLIGHT), str(work), "--deck", str(work / "input.nml"),
         "--restart", str(work / restart_name), "--guard", str(LAUNCHER)],
        cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    if pre.returncode:
        raise RuntimeError(f"existing preflight failed for {work}:\n{pre.stdout[-4000:]}")
    (work / "preflight.json").write_text(pre.stdout)
    df_pre = subprocess.check_output(["df", "-k", "/private/tmp"], text=True)
    disk_samples = [shutil.disk_usage("/private/tmp").free]
    if disk_samples[-1] < 1_073_741_824:
        raise RuntimeError(f"1 GiB /private/tmp free-space guard tripped before {work}")
    proc = subprocess.Popen(
        ["/bin/bash", str(LAUNCHER), str(binary)], cwd=work, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, start_new_session=True,
    )
    while proc.poll() is None:
        time.sleep(1)
        free = shutil.disk_usage("/private/tmp").free
        disk_samples.append(free)
        if free < 1_073_741_824:
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            raise RuntimeError(f"1 GiB /private/tmp live guard stopped {work}")
    stdout = proc.stdout.read() if proc.stdout else ""
    log = work / "testrun.log"
    output = log.read_text() if log.exists() else stdout
    output = re.sub(r"\x1b\[[0-9;]*m", "", output)
    (work / "run_monitor.json").write_text(json.dumps({
        "preflight": str(work / "preflight.json"), "thread_caps": THREAD_CAPS,
        "disk_free_bytes_samples": disk_samples,
        "df_pre": df_pre, "df_post": subprocess.check_output(["df", "-k", "/private/tmp"], text=True),
        "disk_free_bytes_before": disk_samples[0], "disk_free_bytes_after": shutil.disk_usage("/private/tmp").free,
        "log_bytes": log.stat().st_size if log.exists() else 0,
        "solver_wrapper_pid": proc.pid, "returncode": proc.returncode,
    }, indent=2) + "\n")
    return proc.returncode, output


def copy_potential(source: Path, name: str, destination: Path) -> None:
    """Copy a potential and, when it names one, its dm sidecar."""
    shutil.copy2(source / "Fe_out.nml", destination / name)
    match = re.search(r"(?m)^\s*dm_file\s*=\s*'([^']*)'", (source / "Fe_out.nml").read_text())
    if match:
        shutil.copy2(source / match.group(1), destination / match.group(1))


def legacy_potential(source: Path) -> str:
    """Potential text with the dm sidecar inlined as ldm_flatten and no provenance fields."""
    text = (source / "Fe_out.nml").read_text()
    sidecar = (source / re.search(r"(?m)^\s*dm_file\s*=\s*'([^']*)'", text).group(1)).read_text()
    rows: dict = {}
    for m in re.finditer(r"dm_(re|im)\((\d+),(\d+),(\d+),(\d+),:\)\s*=\s*(.*)", sidecar):
        part, l, spin_row, spin_col, row = m.group(1), *(int(m.group(i)) for i in (2, 3, 4, 5))
        if spin_row == spin_col:
            rows.setdefault((part, l, spin_row), {})[row] = m.group(6).split()
    inline = ""
    for (part, l, spin), block in sorted(rows.items()):
        values = [v for row in sorted(block) for v in block[row]]
        values += ["0.0"] * (25 - len(values))
        name = "ldm_flatten" if part == "re" else "ldm_flatten_imag"
        inline += f" {name}({l}, {spin}, :) = " + ", ".join(values) + "\n"
    text = re.sub(r"(?m)^\s*dm_file\s*=.*\n", inline, text)
    return re.sub(r"(?im)^\s*(dm_sha256|hubbard_radial_moments|hubbard_backend)\s*=.*\n", "", text)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fresh_work(root: Path, name: str) -> Path:
    work = root / name
    shutil.copytree(FIXTURE, work)
    return work


def configure_impurity(work: Path, mode: str | None) -> None:
    path = work / "input.nml"
    text = path.read_text()
    text = re.sub(r"(?ims)(^\s*&hamiltonian\b.*?)(^\s*/\s*$)",
                  lambda m: re.sub(r"(?im)^\s*hubbard_radial_moments\s*=.*\n", "", m.group(1)) +
                  (f"hubbard_radial_moments = '{mode}'\n" if mode else "") + m.group(2), text, count=1)
    if not re.search(r"(?im)^\s*hubbard_u_general\(1,3\)", text):
        text = re.sub(r"(?im)^\s*&hamiltonian\s*$",
                      "&hamiltonian\nhubbard_u_general(1,3) = 2.0\nhubbard_j_general(1,3) = 0.5\n", text, count=1)
    path.write_text(text)


def fresh_impurity(root: Path, name: str) -> Path:
    work = root / name
    shutil.copytree(IMPURITY_FIXTURE, work)
    return work


def fresh_surface(root: Path, name: str) -> Path:
    work = root / name
    shutil.copytree(SURFACE_FIXTURE, work)
    return work


def configure_surface(work: Path) -> None:
    path = work / "input.nml"
    text = path.read_text().replace("nsp = 2", "nsp = 1")
    text += "\n&hamiltonian\nhubbard_radial_moments = 'spectral'\n"
    text += "hubbard_u_general(1,3) = 2.0\nhubbard_j_general(1,3) = 0.0\n"
    text += "hubbard_u_general(3,3) = 2.0\nhubbard_j_general(3,3) = 0.0\n/\n"
    path.write_text(text)


def stamp_potential(path: Path, mode: str) -> None:
    text = path.read_text()
    match = re.search(r"(?ims)(^\s*&par\b)(.*?)(^\s*/\s*$)", text)
    if not match:
        raise RuntimeError(f"potential {path} has no &par group")
    body = re.sub(r"(?im)^\s*hubbard_(?:u|j|radial_moments|backend)\s*=.*\n?", "", match.group(2))
    body += ("\n hubbard_u = 0.0, 0.0, 2.0\n hubbard_j = 0.0, 0.0, 0.0\n" +
             f" hubbard_radial_moments = '{mode}'\n hubbard_backend = 'ks'\n")
    text = text[:match.start()] + match.group(1) + body + match.group(3) + text[match.end():]
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--scratch-root", type=Path)
    args = parser.parse_args()
    scratch = args.scratch_root or Path(tempfile.mkdtemp(prefix="rs_embedding_ra2_"))
    scratch.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve()

    # A fresh U run writes its resolved mode and backend. This fixture is a
    # reciprocal KS bulk run using an ordinary no-U atomic seed.
    seed = fresh_work(scratch, "ks_local_seed")
    configure(seed, backend="ks", mode="local")
    status, log = run(binary, seed)
    if status != 0:
        raise RuntimeError(f"KS local seed failed:\n{log[-4000:]}")
    potential = (seed / "Fe_out.nml").read_text().lower()
    if "hubbard_radial_moments" not in potential or "'local'" not in potential:
        raise RuntimeError("U-carrying potential output omitted its resolved local mode")
    if "hubbard_backend" not in potential or "'ks'" not in potential:
        raise RuntimeError("U-carrying potential output omitted its KS backend")

    # Bulk auto restart converts both a provenance-free legacy seed and a
    # spectral-stamped seed to the selected local mode, then stamps this run.
    bulk_legacy = fresh_work(scratch, "bulk_legacy_auto_local")
    (bulk_legacy / "Fe_out.nml").write_text(legacy_potential(seed))
    configure(bulk_legacy, backend="ks", mode=None)
    status, log = run(binary, bulk_legacy)
    legacy_output = (bulk_legacy / "Fe_out.nml").read_text().lower()
    if (status != 0 or "seed in spectral, run uses local; seed will be rewritten" not in log.lower() or
            "hubbard_radial_moments" not in legacy_output or "'local'" not in legacy_output or "'ks'" not in legacy_output):
        raise RuntimeError(f"bulk legacy U restart did not write provenance:\n{log[-4000:]}")

    bulk_spectral = fresh_work(scratch, "bulk_spectral_seed")
    configure(bulk_spectral, backend="rs", mode="spectral")
    status, log = run(binary, bulk_spectral)
    if status != 0:
        raise RuntimeError(f"bulk spectral seed run failed:\n{log[-4000:]}")
    bulk_auto = fresh_work(scratch, "bulk_spectral_auto_local")
    copy_potential(bulk_spectral, "Fe_out.nml", bulk_auto)
    configure(bulk_auto, backend="ks", mode=None)
    status, log = run(binary, bulk_auto)
    rewritten = (bulk_auto / "Fe_out.nml").read_text().lower()
    if status != 0 or "seed in spectral, run uses local; seed will be rewritten" not in log.lower():
        raise RuntimeError(f"bulk spectral seed was not accepted and rewritten locally:\n{log[-4000:]}")
    if "'local'" not in rewritten or "'ks'" not in rewritten:
        raise RuntimeError("auto-local bulk restart did not write local/KS provenance")

    # A matching local-mode RS restart is accepted; the backend is informational.
    matching = fresh_work(scratch, "rs_matching_local")
    copy_potential(seed, "Fe_out.nml", matching)
    configure(matching, backend="rs", mode="local")
    status, log = run(binary, matching)
    if status != 0:
        raise RuntimeError(f"matching KS-host/RS local-mode run failed:\n{log[-4000:]}")

    # Backend provenance is informational: an RS consumer must accept the
    # same local-mode host even when the backend field is absent.
    missing_backend = fresh_work(scratch, "rs_missing_backend_local")
    copy_potential(seed, "Fe_out.nml", missing_backend)
    backendless_text = re.sub(r"(?im)^\s*hubbard_backend\s*=.*\n", "", (seed / "Fe_out.nml").read_text())
    (missing_backend / "Fe_out.nml").write_text(backendless_text)
    configure(missing_backend, backend="rs", mode="local")
    status, log = run(binary, missing_backend)
    if status != 0:
        raise RuntimeError(f"local-mode potential without backend provenance was rejected:\n{log[-4000:]}")

    # An updated bulk atom may convert a stored local seed to the selected
    # spectral mode, even if auto selected local for the KS seed run.
    mismatch = fresh_work(scratch, "rs_mode_mismatch")
    copy_potential(seed, "Fe_out.nml", mismatch)
    configure(mismatch, backend="rs", mode="spectral")
    status, log = run(binary, mismatch)
    rewritten = (mismatch / "Fe_out.nml").read_text().lower()
    if status != 0 or "seed in local, run uses spectral; seed will be rewritten" not in log.lower():
        raise RuntimeError(f"updated bulk seed did not convert to spectral:\n{log[-4000:]}")
    if "'spectral'" not in rewritten or "'rs'" not in rewritten:
        raise RuntimeError("updated bulk spectral run did not write spectral/RS provenance")

    # A legacy U-carrying file without the new fields is classified as spectral.
    legacy = fresh_work(scratch, "rs_legacy_spectral")
    (legacy / "Fe_out.nml").write_text(legacy_potential(seed))
    configure(legacy, backend="rs", mode="spectral")
    status, log = run(binary, legacy)
    legacy_output = (legacy / "Fe_out.nml").read_text().lower()
    if (status != 0 or "mode is unknown" in log.lower() or
            "seed in spectral, run uses spectral; seed will be rewritten" not in log.lower() or
            "'spectral'" not in legacy_output or "'rs'" not in legacy_output):
        raise RuntimeError(f"legacy updated atom seed was not reported and rewritten as spectral:\n{log[-4000:]}")

    # In an impurity calculation, only frozen host ownership is strict. A
    # local frozen host conflicts fatally with forced spectral, while a legacy
    # frozen host is warned about and classified spectral.
    impurity_mismatch = fresh_impurity(scratch, "impurity_frozen_local_forced_spectral")
    copy_potential(seed, "Fe.nml", impurity_mismatch)
    configure_impurity(impurity_mismatch, "spectral")
    status, log = run(binary, impurity_mismatch)
    if status == 0 or "fixed-u potential moment-mode mismatch" not in log.lower():
        raise RuntimeError(f"frozen local host did not reject forced spectral impurity:\n{log[-4000:]}")

    impurity_legacy = fresh_impurity(scratch, "impurity_legacy_frozen_host")
    (impurity_legacy / "Fe.nml").write_text(legacy_potential(seed))
    configure_impurity(impurity_legacy, None)
    status, log = run(binary, impurity_legacy)
    if (status != 0 or "mode is unknown, treating it as spectral" not in log.lower() or
            not re.search(r"radial_mode[ =]+(mode=)?spectral", log.lower())):
        raise RuntimeError(f"legacy frozen host did not warn and proceed in spectral mode:\n{log[-4000:]}")

    # Surface ownership uses its own recursion-site map. The frozen type 1
    # host remains strict; a recursion type may rewrite its stored mode.
    surface_host = fresh_surface(scratch, "surface_frozen_local_forced_spectral")
    configure_surface(surface_host)
    stamp_potential(surface_host / "Cu.nml", "local")
    status, log = run(binary, surface_host, restart_name="Cu.nml")
    if status == 0 or "fixed-u potential moment-mode mismatch" not in log.lower():
        raise RuntimeError(f"surface frozen local host did not reject spectral mode:\n{log[-4000:]}")

    surface_rec = fresh_surface(scratch, "surface_recursion_local_to_spectral")
    configure_surface(surface_rec)
    stamp_potential(surface_rec / "Cu.nml", "spectral")
    stamp_potential(surface_rec / "Cu-S.nml", "local")
    status, log = run(binary, surface_rec, restart_name="Cu.nml")
    rec_output = (surface_rec / "Cu-S_out.nml").read_text().lower()
    if status != 0 or "seed in local, run uses spectral; seed will be rewritten" not in log.lower():
        raise RuntimeError(f"surface recursion seed did not convert under SCF ownership:\n{log[-4000:]}")
    if "'spectral'" not in rec_output or "'rs'" not in rec_output:
        raise RuntimeError("surface recursion output omitted spectral/RS provenance")

    # No-U potential output must not acquire provenance fields.
    no_u = fresh_work(scratch, "no_u")
    configure(no_u, backend="ks", mode=None, uj=False)
    status, log = run(binary, no_u)
    if status != 0:
        raise RuntimeError(f"no-U run failed:\n{log[-4000:]}")
    no_u_text = (no_u / "Fe_out.nml").read_text().lower()
    if "hubbard_radial_moments" in no_u_text or "hubbard_backend" in no_u_text:
        raise RuntimeError("no-U output gained fixed-U provenance fields")

    scenarios = {}
    for name in ("ks_local_seed", "bulk_legacy_auto_local", "bulk_spectral_seed", "bulk_spectral_auto_local",
                 "rs_matching_local", "rs_missing_backend_local", "rs_mode_mismatch", "rs_legacy_spectral",
                 "impurity_frozen_local_forced_spectral", "impurity_legacy_frozen_host",
                 "surface_frozen_local_forced_spectral", "surface_recursion_local_to_spectral", "no_u"):
        work = scratch / name
        record = {"input": str(work / "input.nml"), "input_sha256": sha256(work / "input.nml"),
                  "log": str(work / "testrun.log"), "log_sha256": sha256(work / "testrun.log")}
        potential_path = work / "Fe_out.nml"
        if potential_path.exists():
            record["potential_output"] = str(potential_path)
            record["potential_output_sha256"] = sha256(potential_path)
        record["artifact_sha256"] = {
            str(path.relative_to(work)): sha256(path)
            for path in sorted(work.rglob("*"))
            if path.is_file() and path.name not in {"RA2_provenance.json"}
        }
        scenarios[name] = record
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unavailable"
    manifest = {
        "stage": "RA2",
        "source_commit": commit,
        "source_sha256": sha256(ROOT / "source/self.f90"),
        "fixture_runner_sha256": sha256(Path(__file__).resolve()),
        "executable": str(binary),
        "executable_sha256": sha256(binary),
        "command": ["python3", str(Path(__file__).resolve()), "--binary", str(binary), "--scratch-root", str(scratch)],
        "raw_output_location": str(scratch),
        "thread_caps": THREAD_CAPS,
        "preflight_path": str(PREFLIGHT),
        "scenarios": scenarios,
    }
    (scratch / "RA2_provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("RA2 PASS: bulk legacy/spectral seed rewrites, frozen-host strict mismatch, legacy host warning, "
          "cross-backend restart, no-U omission")


if __name__ == "__main__":
    main()
