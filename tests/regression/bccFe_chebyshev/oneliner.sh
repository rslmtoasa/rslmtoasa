#!/bin/bash

if [ $# -eq 1 ]; then
    binary="$1"
else
    binary=../../../build/bin/rslmto.x
fi

rm -f Fe_out.nml data.nml ref.nml fort.* clust mad.mat map sbar str.out ves.out view.sbar
$binary > testrun.log 2>&1
mv Fe_out.nml data.nml
cp Fe.nml.ref ref.nml
pytest --tb=line  --no-header ../test_comparison.py
rm -f data.nml ref.nml fort.* clust mad.mat map sbar str.out ves.out view.sbar
if [ $? -ne 0 ]; then
    echo "$0: ERROR: $binary failed"
    exit 1
fi
