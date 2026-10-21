#!/usr/bin/env python3
"""PRIVATE summary of completed reciprocal [110]/[-1-10] endpoints.

The original comparison is explicitly the supplied saved output benchmark;
this command does not claim a new unsymmetrized same-executable diagonal run.
"""
import argparse,json
from pathlib import Path
import numpy as np
from causal_summary import complete,directory,parity


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reciprocal',type=Path,required=True)
    p.add_argument('--fixtures',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    cases=('110','1-1-0')
    manifest=complete(args.reciprocal,cases)
    out={'binary_sha256':manifest['binary_sha256'],'cases':{},
         'original_comparison':'Supplied fixture outputs; no new unsymmetrized diagonal calculation',
         'original_benchmark':{}}
    data={}
    for case in cases:
        item=manifest['cases'][case]
        if not item['sbar_override_loaded'] or item['changed_on_copy']!=['sbar.']:
            raise ValueError('Incorrect intervention')
        out['cases'][case],data[case]=directory(args.reciprocal/case,neel='110')
        job=args.fixtures/case/'exchange_jobs/Mn1'
        d=np.loadtxt(job/'dijso.out',ndmin=2)
        parts=np.loadtxt(job/'dijparts.out',ndmin=2)
        j=np.loadtxt(job/'jijso.out',ndmin=2)
        if d.shape[0]!=24 or not np.array_equal(d[:,:5],parts[:,:5]) or not np.array_equal(d[:,:5],j[:,:5]):
            raise ValueError('Benchmark pair coverage/order mismatch')
        rec=out['cases'][case]
        if not np.array_equal(d[:,:5],np.array([b['key'] for b in rec['bonds']])):
            raise ValueError('Benchmark/intervention pair ordering mismatch')
        mask=(d[:,0]==1)&(d[:,1]==2)
        out['original_benchmark'][case]={
            'sum_D_mRy':d[mask,5:8].sum(axis=0).tolist(),
            'sum_CC_mRy':parts[mask,5:8].sum(axis=0).tolist(),
            'sum_SC_mRy':parts[mask,8:11].sum(axis=0).tolist(),
            'D_mRy':d[:,5:8].tolist(),'CC_mRy':parts[:,5:8].tolist(),
            'SC_mRy':parts[:,8:11].tolist(),'Jijso_mRy':j[:,5].tolist()}
    out['occupied_domain_residuals']=parity(data[cases[0]],data[cases[1]])
    for kind in ('D','CC','SC'):
        a,b=[np.array(out['cases'][c]['sum_'+kind+'_mRy']) for c in cases]
        out[kind+'_odd_mRy']=((a-b)/2).tolist()
        out[kind+'_even_mRy']=((a+b)/2).tolist()
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('reciprocal diagonal odd D (mRy):',out['D_odd_mRy'])
    print('reciprocal diagonal even D (mRy):',out['D_even_mRy'])


if __name__=='__main__': main()
