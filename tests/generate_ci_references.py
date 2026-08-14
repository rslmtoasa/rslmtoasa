#!/usr/bin/env python3
"""Generate references from an existing CI-equivalent build.

This is intentionally a separate entry point from the generic reference
generator.  It makes the reference provenance explicit and keeps the stored
baselines tied to the supplied build directory rather than to an arbitrary
developer build.  Configure and build the binary separately; this script only
validates the build and runs the reference cases.

Usage:
  tests/generate_ci_references.py --build-dir BUILD [--case NAME ...]
"""

from __future__ import annotations

import argparse
import json
import os
import platform
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
CI_SERIAL_OMP_THREADS = 2


def run(command: list[str], *, cwd: Path | None = None, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=cwd, env=env, check=True)


def command_text(command: list[str]) -> str:
    try:
        return subprocess.check_output(command, text=True, stderr=subprocess.STDOUT).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def openmpi_environment() -> dict[str, str]:
    """Load the repository's supported Open MPI environment in a subprocess."""
    script = ROOT / "env" / "openmpi.sh"
    if not script.exists():
        raise RuntimeError(f"missing Open MPI environment script: {script}")
    dump = subprocess.check_output(
        ["bash", "-c", f"source {script} >/dev/null; env -0"],
        text=False,
    )
    environment: dict[str, str] = {}
    for item in dump.split(b"\0"):
        if b"=" in item:
            key, value = item.split(b"=", 1)
            environment[key.decode()] = value.decode()
    return environment


def runner_native_environment() -> dict[str, str]:
    """Return the environment used by the platform CI runner itself.

    This deliberately does not source env/openmpi.sh.  The runner workflow
    installs its dependencies in a clean environment, and CMake should resolve
    the compiler, MPI, and BLAS/LAPACK installation available there.
    """
    return os.environ.copy()


def resolve_mpi_launcher(environment: dict[str, str]) -> str:
    """Find an MPI launcher in the active runner environment."""
    requested = environment.get("RSLMTO_MPI_LAUNCHER")
    candidates = [requested] if requested else []
    candidates.extend(["mpirun.openmpi", "mpirun", "mpiexec"])
    for candidate in candidates:
        if candidate:
            resolved = shutil.which(candidate, path=environment.get("PATH"))
            if resolved:
                return resolved
    raise RuntimeError("no MPI launcher found in the runner environment")


def cache_values(cache_path: Path, keys: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    if not cache_path.exists():
        return values
    wanted = set(keys)
    for line in cache_path.read_text().splitlines():
        if "=" not in line or ":" not in line or line.startswith("//"):
            continue
        key, value = line.split("=", 1)
        key = key.split(":", 1)[0]
        if key in wanted:
            values[key] = value
    return values


def cmake_option_enabled(cache_path: Path, option: str) -> bool:
    return cache_values(cache_path, [option]).get(option, "OFF").upper() in {
        "ON", "TRUE", "1"
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build-dir", required=True,
        help=("Existing CMake build directory containing bin/rslmto.x; "
              "this script does not configure or build it"),
    )
    parser.add_argument(
        "--cases-json", default="tests/scf/cases.json",
        help="Case manifest to generate (default: tests/scf/cases.json)",
    )
    parser.add_argument(
        "--references-dir", default="tests/scf/references",
        help="Reference directory (default: tests/scf/references)",
    )
    parser.add_argument(
        "--scratch-root", default=None,
        help="Scratch root passed to the reference generator",
    )
    parser.add_argument(
        "--profile", choices=("linux-ci-equivalent", "runner-native"),
        default="linux-ci-equivalent",
        help=("Execution environment for reference runs. 'linux-ci-equivalent' "
              "uses the repository's local Open MPI environment; "
              "'runner-native' leaves the current runner environment untouched."),
    )
    parser.add_argument(
        "--case", nargs="+", default=[], metavar="NAME",
        help="Generate only the named cases (default: all cases)",
    )
    args = parser.parse_args()

    build_dir = Path(args.build_dir)
    if not build_dir.is_absolute():
        build_dir = ROOT / build_dir
    cases_json = Path(args.cases_json)
    if not cases_json.is_absolute():
        cases_json = ROOT / cases_json
    references_dir = Path(args.references_dir)
    if not references_dir.is_absolute():
        references_dir = ROOT / references_dir

    cache_path = build_dir / "CMakeCache.txt"
    binary = build_dir / "bin" / "rslmto.x"
    if not cache_path.is_file():
        raise RuntimeError(
            f"existing CMake build directory required: {cache_path} was not found; "
            "configure the build before running this script"
        )
    if not binary.is_file():
        raise RuntimeError(
            f"existing reference binary required: {binary} was not found; "
            "build the target before running this script"
        )

    cache = cache_values(
        cache_path,
        [
            "CMAKE_BUILD_TYPE", "CMAKE_Fortran_COMPILER", "ENABLE_MPI",
            "ENABLE_MARCH_NATIVE", "ENABLE_MKL_KERNELS", "BLA_VENDOR",
            "BLAS_LIBRARIES", "LAPACK_LIBRARIES",
        ],
    )
    mpi_enabled = cmake_option_enabled(cache_path, "ENABLE_MPI")

    # The committed Linux profile uses the repository-supported environment to
    # isolate local MPI installations.  The runner-native profile intentionally
    # leaves the runner environment untouched and lets the existing build's
    # runtime dependencies resolve in that environment.
    ci_env = (openmpi_environment() if args.profile == "linux-ci-equivalent"
              else runner_native_environment())
    if args.profile == "linux-ci-equivalent":
        for variable in ("CMAKE_PREFIX_PATH", "CPATH", "C_INCLUDE_PATH", "CPLUS_INCLUDE_PATH"):
            ci_env.pop(variable, None)
    if mpi_enabled:
        ci_env["RSLMTO_MPI_LAUNCHER"] = resolve_mpi_launcher(ci_env)
    ci_env["RSLMTO_OMP_THREADS_SERIAL"] = str(CI_SERIAL_OMP_THREADS)
    effective_compiler = cache.get("CMAKE_Fortran_COMPILER", "unknown")
    manifest = json.loads(cases_json.read_text())
    available_cases = {
        case["name"]: case for case in manifest.get("cases", [])
    }
    requested_cases = args.case or list(available_cases)
    selected_cases: list[str] = []
    for name in requested_cases:
        if name not in available_cases:
            raise RuntimeError(f"case {name!r} not found in {cases_json}")
        required = available_cases[name].get("requires_cmake_option")
        if required and not cmake_option_enabled(build_dir / "CMakeCache.txt", required):
            if args.case:
                raise RuntimeError(
                    f"case {name!r} requires {required}=ON, but the canonical "
                    "build has it disabled"
                )
            print(f"SKIP  [{name}]: {required}=OFF in canonical build")
            continue
        selected_cases.append(name)

    provenance = {
        "profile": args.profile,
        "platform": platform.platform(),
        "git_commit": command_text(["git", "-C", str(ROOT), "rev-parse", "HEAD"]),
        "compiler": effective_compiler,
        "compiler_version": (
            command_text([effective_compiler, "--version"])
            if effective_compiler != "unknown" else "unavailable"
        ),
        "underlying_compiler": (
            command_text([effective_compiler, "--showme:command"])
            if effective_compiler != "unknown" else "unavailable"
        ),
        "cmake_version": command_text(["cmake", "--version"]).splitlines()[0],
        "cmake_cache": cache,
        # Keep committed metadata relocatable; the build directory is printed
        # by this script and remains available to run additional tests.
        "binary": binary.name,
        "serial_omp_threads": CI_SERIAL_OMP_THREADS,
        "mpi_omp_threads": 1,
        "mpi_enabled": mpi_enabled,
    }

    env = ci_env.copy()
    env["RSLMTO_REFERENCE_PROVENANCE"] = json.dumps(provenance, sort_keys=True)
    generator = [
        sys.executable, str(ROOT / "tests" / "generate_references.py"),
        "--binary", str(binary),
        "--cases-json", str(cases_json),
        "--references-dir", str(references_dir),
    ]
    if args.scratch_root:
        scratch_root = Path(args.scratch_root)
        if not scratch_root.is_absolute():
            scratch_root = ROOT / scratch_root
        generator.extend(["--scratch-root", str(scratch_root)])
    if selected_cases:
        generator.extend(["--case", *selected_cases])
    run(generator, env=env)

    print("References generated from the existing build.")
    print(f"Binary: {binary}")
    print(f"References: {references_dir}")


if __name__ == "__main__":
    main()
