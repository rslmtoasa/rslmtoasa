#!/bin/bash
# Thin wrapper for binary invocation.
# Usage: run_binary.sh <binary> [mpi_procs]
# Must be run from the test scratch directory (cwd). Output is written to testrun.log.

set -euo pipefail

binary="$1"
mpi_procs="${2:-1}"

if [ ! -x "$binary" ]; then
    echo "ERROR: binary not executable: $binary" >&2
    exit 1
fi

echo "INFO: run_binary cwd=$(pwd) binary=${binary} mpi_procs=${mpi_procs}" >&2

if [ "$mpi_procs" -gt 1 ]; then
    # Respect an explicitly selected launcher from CI if provided.
    mpi_launcher="${RSLMTO_MPI_LAUNCHER:-}"
    if [ -n "$mpi_launcher" ]; then
        if ! command -v "$mpi_launcher" &>/dev/null; then
            echo "ERROR: requested MPI launcher '$mpi_launcher' not found in PATH" >&2
            exit 3
        fi
    elif command -v mpirun &>/dev/null; then
        mpi_launcher="mpirun"
    elif command -v mpiexec &>/dev/null; then
        mpi_launcher="mpiexec"
    else
        echo "ERROR: MPI run requested (mpi_procs=${mpi_procs}) but no mpirun/mpiexec found on PATH" >&2
        exit 3
    fi

    if [ "$mpi_launcher" = "mpirun" ] || [ "$mpi_launcher" = "mpirun.openmpi" ]; then
        run_cmd=("$mpi_launcher" --oversubscribe -n "$mpi_procs" "$binary")
    else
        run_cmd=("$mpi_launcher" -n "$mpi_procs" "$binary")
    fi

    echo "INFO: launcher=${mpi_launcher}" >&2
    printf 'INFO: command=' >&2
    printf ' %q' "${run_cmd[@]}" >&2
    printf '\n' >&2

    {
        echo "[run_binary] cwd=$(pwd)"
        echo "[run_binary] mpi_procs=${mpi_procs} launcher=${mpi_launcher}"
        printf '[run_binary] command='
        printf ' %q' "${run_cmd[@]}"
        printf '\n'
    } > testrun.log

    set +e
    OMP_NUM_THREADS=1 "${run_cmd[@]}" >> testrun.log 2>&1
    rc=$?
    set -e
    if [ "$rc" -ne 0 ]; then
        echo "ERROR: launcher command failed with exit code ${rc}" >&2
        echo "--- testrun.log (tail -n 120) ---" >&2
        tail -n 120 testrun.log >&2 || true
        echo "----------------------------------" >&2
        exit "$rc"
    fi
else
    set +e
    OMP_NUM_THREADS=1 "$binary" > testrun.log 2>&1
    rc=$?
    set -e
    if [ "$rc" -ne 0 ]; then
        echo "ERROR: serial run failed with exit code ${rc}" >&2
        echo "--- testrun.log (tail -n 120) ---" >&2
        tail -n 120 testrun.log >&2 || true
        echo "----------------------------------" >&2
        exit "$rc"
    fi
fi
