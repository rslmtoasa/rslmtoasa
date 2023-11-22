#!/bin/bash
# Simple script to run the RS-LMTO code and compare outputs for testing
# Usage `run_simple_script.sh <prefix> <RS-LMTO binary path> <test script path>
# i.e. `./run_simple_script.sh Fe /tmp/build/bin /tmp/tests/test-regression

export OMP_NuM_THREADS=1
rm -f ${1}_old.nml
${2}/rslmto.x > testrun.log
cp ${1}_out.nml data.nml
cp ${1}.nml.ref ref.nml
pytest --tb=line  --no-header  ${3}/test_comparison.py
# pytest --no-header  ${3}/test_regression.py
rm -f data.nml ref.nml fort.* clust mad.mat map sbar str.out ves.out view.sbar
