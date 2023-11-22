#!/bin/bash
rm -f Fe_out.nml ; ../../../build/bin/rslmto.x > testrun.log ; mv Fe_out.nml data.nml ; cp Fe.nml.ref ref.nml ;   pytest --tb=line  --no-header ../test_comparison.py ; rm -f data.nml ref.nml fort.* clust mad.mat map sbar str.out ves.out view.sbar
