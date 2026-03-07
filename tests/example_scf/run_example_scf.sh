#!/bin/bash

set -euo pipefail

if [ "$#" -ne 9 ]; then
    echo "Usage: $0 <binary> <case_dir> <scratch_root> <test_name> <nsp> <recur> <hoh:true|false> <lld> <nstep>"
    exit 2
fi

binary="$1"
binary_dir="$(cd "$(dirname "$binary")" && pwd)"
binary_basename="$(basename "$binary")"
binary="${binary_dir}/${binary_basename}"
case_dir="$2"
scratch_root="$3"
test_name="$4"
nsp="$5"
recur="$6"
hoh_raw="$7"
lld="$8"
nstep="$9"

if [ ! -x "$binary" ]; then
    echo "ERROR: binary not executable: $binary"
    exit 1
fi

if [ ! -d "$case_dir" ]; then
    echo "ERROR: case directory not found: $case_dir"
    exit 1
fi

case "${hoh_raw}" in
    true|TRUE|True) hoh_value=".true." ;;
    false|FALSE|False) hoh_value=".false." ;;
    *)
        echo "ERROR: hoh must be true or false"
        exit 2
        ;;
esac

workdir="${scratch_root}/${test_name}"
rm -rf "$workdir"
mkdir -p "$workdir"
cp -R "${case_dir}/." "$workdir/"
rm -f "$workdir"/*_out.nml "$workdir"/data.nml

input_file="${workdir}/input.nml"
if [ ! -f "$input_file" ]; then
    echo "ERROR: input file not found in copied example: $input_file"
    exit 1
fi

modified_input="${workdir}/input.modified.nml"

awk \
    -v nsp_val="$nsp" \
    -v recur_val="$recur" \
    -v lld_val="$lld" \
    -v nstep_val="$nstep" \
    -v hoh_val="$hoh_value" '
function lower(s,  t) {
    t = s
    gsub(/[A-Z]/, "", t)
    return tolower(s)
}
BEGIN {
    section = ""
    seen_hamiltonian = 0
    ctrl_nsp = 0
    ctrl_recur = 0
    ctrl_lld = 0
    self_nstep = 0
    ham_hoh = 0
}
{
    line = $0
    lline = tolower(line)

    if (lline ~ /^[[:space:]]*&[a-z_][a-z0-9_]*/) {
        if (lline ~ /^[[:space:]]*&control/) section = "control"
        else if (lline ~ /^[[:space:]]*&self/) section = "self"
        else if (lline ~ /^[[:space:]]*&hamiltonian/) {
            section = "hamiltonian"
            seen_hamiltonian = 1
        } else section = "other"
    }

    if (section == "control") {
        if (lline ~ /^[[:space:]]*nsp[[:space:]]*=/) {
            print "nsp = " nsp_val
            ctrl_nsp = 1
            next
        }
        if (lline ~ /^[[:space:]]*recur[[:space:]]*=/) {
            print "recur = '\''" recur_val "'\''"
            ctrl_recur = 1
            next
        }
        if (lline ~ /^[[:space:]]*lld[[:space:]]*=/) {
            print "lld = " lld_val
            ctrl_lld = 1
            next
        }
    }

    if (section == "self") {
        if (lline ~ /^[[:space:]]*nstep[[:space:]]*=/) {
            print "nstep = " nstep_val
            self_nstep = 1
            next
        }
    }

    if (section == "hamiltonian") {
        if (lline ~ /^[[:space:]]*hoh[[:space:]]*=/) {
            print "hoh = " hoh_val
            ham_hoh = 1
            next
        }
    }

    if (lline ~ /^[[:space:]]*\/[[:space:]]*$/) {
        if (section == "control") {
            if (!ctrl_nsp) print "nsp = " nsp_val
            if (!ctrl_recur) print "recur = '\''" recur_val "'\''"
            if (!ctrl_lld) print "lld = " lld_val
        } else if (section == "self") {
            if (!self_nstep) print "nstep = " nstep_val
        } else if (section == "hamiltonian") {
            if (!ham_hoh) print "hoh = " hoh_val
        }
        print line
        section = ""
        next
    }

    print line
}
END {
    if (!seen_hamiltonian) {
        print ""
        print "&hamiltonian"
        print "hoh = " hoh_val
        print "/"
    }
}
' "$input_file" > "$modified_input"

mv "$modified_input" "$input_file"

(
    cd "$workdir"
    OMP_NUM_THREADS=1 "$binary" > testrun.log 2>&1
)

if grep -qi "fatal" "${workdir}/testrun.log"; then
    echo "ERROR: fatal message detected in ${workdir}/testrun.log"
    tail -n 50 "${workdir}/testrun.log"
    exit 1
fi

# Normalize the output so downstream scripts can always read data.nml.
if shopt -q nullglob; then
    nullglob_pre=1
else
    nullglob_pre=0
fi
shopt -s nullglob
out_candidates=("$workdir"/*_out.nml)
if [ "$nullglob_pre" -eq 1 ]; then
    shopt -s nullglob
else
    shopt -u nullglob
fi

selected=""
for candidate in "${out_candidates[@]}"; do
    if [ -f "$candidate" ]; then
        base_candidate="${candidate##*/}"
        if [ -z "$selected" ] || [[ "$base_candidate" < "${selected##*/}" ]]; then
            selected="$candidate"
        fi
    fi
done

if [ -z "$selected" ]; then
    echo "ERROR: expected *_out.nml missing in ${workdir}"
    exit 1
fi

if [ "${#out_candidates[@]}" -gt 1 ]; then
    echo "NOTE: multiple *_out.nml files detected, keeping ${selected##*/} as data.nml"
fi

rm -f "$workdir/data.nml"
mv "$selected" "$workdir/data.nml"

echo "PASS ${test_name}"
