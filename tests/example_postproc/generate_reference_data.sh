#!/bin/bash

set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
manifest="${script_dir}/cases_manifest.csv"
runner="${repo_root}/tests/example_scf/run_example_scf.sh"

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <binary> [--include-cheb] [--scratch-root <path>]"
    exit 2
fi

binary="$1"
shift

include_cheb=0
scratch_root="${repo_root}/build_main/Testing/example_postproc_refgen"
references_root="${script_dir}/references"

while [ "$#" -gt 0 ]; do
    case "$1" in
        --include-cheb)
            include_cheb=1
            shift
            ;;
        --scratch-root)
            scratch_root="$2"
            shift 2
            ;;
        *)
            echo "ERROR: unknown option $1"
            exit 2
            ;;
    esac
done

if [ ! -x "$binary" ]; then
    echo "ERROR: binary not executable: $binary"
    exit 1
fi

mkdir -p "$references_root"

while IFS=',' read -r test_name case_rel nsp recur hoh lld nstep; do
    if [[ -z "${test_name}" ]] || [[ "${test_name:0:1}" == "#" ]]; then
        continue
    fi

    if [[ "$recur" == "chebyshev" && "$include_cheb" -eq 0 ]]; then
        echo "SKIP ${test_name} (chebyshev disabled; use --include-cheb to include)"
        continue
    fi

    case_dir="${script_dir}/cases/${case_rel}"
    if [ ! -d "$case_dir" ]; then
        echo "ERROR: missing test case directory: $case_dir"
        exit 1
    fi
    echo "RUN  ${test_name}"
    /bin/bash "$runner" "$binary" "$case_dir" "$scratch_root" "$test_name" "$nsp" "$recur" "$hoh" "$lld" "$nstep"

    run_dir="${scratch_root}/${test_name}"
    data_file="${run_dir}/data.nml"
    input_file="${run_dir}/input.nml"

    if [ ! -f "$data_file" ]; then
        echo "ERROR: expected output missing: ${data_file}"
        exit 1
    fi

    ref_dir="${references_root}/${test_name}"
    mkdir -p "$ref_dir"
    cp "$data_file" "${ref_dir}/ref.nml"
    cp "$input_file" "${ref_dir}/input.nml"

    cat > "${ref_dir}/meta.txt" <<EOF
name=${test_name}
case_rel=${case_rel}
nsp=${nsp}
recur=${recur}
hoh=${hoh}
lld=${lld}
nstep=${nstep}
EOF

done < "$manifest"

echo "Reference generation complete at ${references_root}"
