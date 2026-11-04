#!/usr/bin/env python3
"""Summarize benchmark.py runs without loading full moment tensors into RAM."""
import argparse
import csv
from pathlib import Path
import numpy as np


def metrics(a, b):
    delta = a-b
    return float(np.max(np.abs(delta))), float(np.linalg.norm(delta)), float(np.linalg.norm(b))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    with (args.directory/'timings.csv').open() as f:
        timings=list(csv.DictReader(f))
    rows=[]
    curves=[]
    for fixed in timings:
        if fixed['version']!='fixed' or fixed['exit_code']!='0': continue
        case=(fixed['ll'],fixed['kind'],fixed['repeat'])
        run=args.directory/fixed['run_dir']
        for baseline in timings:
            if baseline['version']=='fixed' or baseline['exit_code']!='0': continue
            if (baseline['ll'],baseline['kind'],baseline['repeat'])!=case: continue
            other=args.directory/baseline['run_dir']
            a=np.loadtxt(run/'Pt_cond.out'); b=np.loadtxt(other/'Pt_cond.out')
            np.testing.assert_array_equal(a[:,0],b[:,0])
            assert np.isfinite(a).all() and np.isfinite(b).all()
            for idx,name in [(1,'conductivity_real'),(2,'conductivity_imag')]:
                maximum,numerator,denominator=metrics(a[:,idx],b[:,idx])
                rows.append([*case,baseline['version'],name,maximum,numerator/max(denominator,1e-300)])
            for x,y in zip(a,b):
                curves.append([*case,baseline['version'],x[0],x[1],x[2],y[1],y[2]])
            pa=run/'bench_moments.bin'; pb=other/'bench_moments.bin'
            assert pa.stat().st_size==pb.stat().st_size
            maximum=error2=reference2=0.
            with pa.open('rb') as fa,pb.open('rb') as fb:
                while True:
                    a=np.fromfile(fa,np.complex128,count=65536)
                    if not a.size: break
                    b=np.fromfile(fb,np.complex128,count=a.size)
                    assert np.isfinite(a).all() and np.isfinite(b).all()
                    mx,err,ref=metrics(a,b)
                    maximum=max(maximum,mx); error2+=err**2; reference2+=ref**2
            rows.append([*case,baseline['version'],'moments',maximum,(error2/max(reference2,1e-300))**.5])
    # Split curves into directly usable CSVs; output precision matches Fortran.
    for ll,kind,repeat,baseline in sorted({tuple(row[:4]) for row in curves}):
        selected=[row[4:] for row in curves if row[:4]==[ll,kind,repeat,baseline]]
        np.savetxt(args.output/f'curves_{ll}_{kind}_{repeat}_{baseline}.csv', selected,
                   delimiter=',',fmt='%.8e',
                   header='energy,fixed_real,fixed_imag,baseline_real,baseline_imag',comments='')
    for filename,header,data in [
        ('differences.csv',['ll','kind','repeat','baseline','quantity','max_abs_difference','relative_l2_difference'],rows)]:
        with (args.output/filename).open('w') as f:
            writer=csv.writer(f,lineterminator="\n"); writer.writerow(header); writer.writerows(data)


if __name__=='__main__': main()
