#!/usr/bin/env python3
"""Build the canonical Linux CI configuration and generate references.

This is intentionally a separate entry point from the generic reference
generator.  It makes the reference provenance explicit and keeps the stored
baselines tied to the production CI toolchain rather than to an arbitrary
developer build directory.

The default configuration mirrors the Linux binary workflow:

* Release build
* Open MPI's ``mpifort`` wrapper (GNU Fortran) with OpenBLAS/LAPACK
* MPI enabled, but no architecture-specific ``-march=native`` flags
* two OpenMP threads for serial test launches, matching ``tests.yml``
* optional MKL Chebyshev kernels disabled

Usage:
  tests/generate_ci_references.py [--case NAME ...]
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
        "--build-dir", default="build-ci-reference",
        help="Build directory (default: build-ci-reference)",
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
        "--fortran-compiler", default=None,
        help="Fortran compiler/wrapper (default: mpifort with MPI, gfortran without MPI)",
    )
    parser.add_argument(
        "--enable-mpi", action=argparse.BooleanOptionalAction, default=True,
        help="Build the reference binary with MPI (default: enabled)",
    )
    parser.add_argument(
        "--profile", choices=("linux-ci-equivalent", "runner-native"),
        default="linux-ci-equivalent",
        help=("Reference build profile. 'linux-ci-equivalent' preserves the "
              "committed Linux/OpenBLAS profile; 'runner-native' mirrors the "
              "current Linux or macOS CI runner."),
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

    # The committed Linux profile uses the repository-supported environment to
    # isolate local MPI installations.  The runner-native profile intentionally
    # leaves the runner environment untouched and lets CMake resolve its tools.
    ci_env = (openmpi_environment() if args.profile == "linux-ci-equivalent"
              else runner_native_environment())
    if args.profile == "linux-ci-equivalent":
        for variable in ("CMAKE_PREFIX_PATH", "CPATH", "C_INCLUDE_PATH", "CPLUS_INCLUDE_PATH"):
            ci_env.pop(variable, None)
    if args.fortran_compiler:
        compiler_name = args.fortran_compiler
    elif args.profile == "linux-ci-equivalent":
        compiler_name = "mpifort" if args.enable_mpi else "gfortran"
    else:
        # This is only used for provenance and fallback lookup.  The
        # runner-native configure below does not force this compiler.
        compiler_name = "gfortran"
    compiler = shutil.which(compiler_name, path=ci_env.get("PATH")) or compiler_name
    if args.enable_mpi:
        ci_env["RSLMTO_MPI_LAUNCHER"] = resolve_mpi_launcher(ci_env)
    ci_env["RSLMTO_OMP_THREADS_SERIAL"] = str(CI_SERIAL_OMP_THREADS)

    configure = [
        "cmake", "-S", str(ROOT), "-B", str(build_dir), "-G", "Ninja",
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DENABLE_MPI={'ON' if args.enable_mpi else 'OFF'}",
        "-DENABLE_MARCH_NATIVE=OFF",
        "-DENABLE_MKL_KERNELS=OFF",
        "-DRUN_REG_TESTS=OFF",
        "-DRUN_EXAMPLE_TESTS=OFF",
        "-DRUN_UNIT_TESTS=OFF",
    ]
    if args.profile == "linux-ci-equivalent" or args.fortran_compiler:
        configure.insert(8, f"-DCMAKE_Fortran_COMPILER={compiler}")
    if args.profile == "linux-ci-equivalent":
        configure.append("-DBLA_VENDOR=OpenBLAS")
    run(configure, env=ci_env)
    run(["cmake", "--build", str(build_dir), "--parallel"], env=ci_env)

    binary = build_dir / "bin" / "rslmto.x"
    if not binary.is_file():
        raise RuntimeError(f"canonical binary was not produced: {binary}")

    cache = cache_values(
        build_dir / "CMakeCache.txt",
        [
            "CMAKE_BUILD_TYPE", "CMAKE_Fortran_COMPILER", "ENABLE_MPI",
            "ENABLE_MARCH_NATIVE", "ENABLE_MKL_KERNELS", "BLA_VENDOR",
            "BLAS_LIBRARIES", "LAPACK_LIBRARIES",
        ],
    )
    effective_compiler = cache.get("CMAKE_Fortran_COMPILER", compiler)
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
        "compiler_version": command_text([effective_compiler, "--version"]),
        "underlying_compiler": command_text([effective_compiler, "--showme:command"]),
        "cmake_version": command_text(["cmake", "--version"]).splitlines()[0],
        "cmake_cache": cache,
        # Keep committed metadata relocatable; the build directory is printed
        # by this script and remains available to run additional tests.
        "binary": binary.name,
        "serial_omp_threads": CI_SERIAL_OMP_THREADS,
        "mpi_omp_threads": 1,
        "mpi_enabled": args.enable_mpi,
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

    print("Canonical CI-equivalent references generated.")
    print(f"Binary: {binary}")
    print(f"References: {references_dir}")


if __name__ == "__main__":
    main()
