#!/bin/bash
# One SCF per q along Gamma-X (Cartesian z, units of 2pi/alat; X = 1; q/2 is on the 32^3 mesh for every value).
# usage: ./run.sh /path/to/rslmto.x      (QZS="0.0 0.5" ./run.sh ... runs a subset)
BIN=$(cd "$(dirname "${1:-rslmto.x}")" && pwd)/$(basename "${1:-rslmto.x}")
cd "$(dirname "$0")"
for qz in ${QZS:-0.0 0.125 0.25 0.375 0.5 0.625 0.75 0.875 1.0}; do
  mkdir -p q_$qz && cp Ni.nml lattice.nml q_$qz/ && sed "s/@QZ@/$qz/" input.template > q_$qz/input.nml
  (cd q_$qz && "$BIN" input.nml > run.log 2>&1)
  echo "q = $qz: $(grep 'Physical DFT total energy' q_$qz/report.out)"
done
