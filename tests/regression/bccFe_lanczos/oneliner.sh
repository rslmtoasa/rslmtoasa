#!/bin/bash
# usage: oneliner.sh <binary> <scratch dir>; runs on a copy of this directory.

src=$(cd "$(dirname "$0")" && pwd)
binary="$1"
work="$2"

rm -rf "$work"
mkdir -p "$work"
cp -R "$src"/. "$work"/
cp "$src"/../test_comparison.py "$work"/
cd "$work" || exit 1

rm -f Fe_out.nml data.nml ref.nml fort.* clust mad.mat map sbar str.out ves.out view.sbar
$binary > testrun.log 2>&1
mv Fe_out.nml data.nml
cp Fe.nml.ref ref.nml
pytest --tb=line  --no-header test_comparison.py
if [ $? -ne 0 ]; then
    echo "$0: ERROR: $binary failed"
    exit 1
fi
