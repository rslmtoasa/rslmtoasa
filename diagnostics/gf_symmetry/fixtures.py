#!/usr/bin/env python3
"""Read-only inventory/regression for external unpublished fixtures.

The small assignment reader handles numeric scalar/vector lines in saved atom files,
not arbitrary Fortran namelists. Calculation execution is handled by run_saved.py.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path
import numpy as np
from symmetry import residual

CASES=('100','1-00','110','1-1-0')
FIELDS=('c(:, 1)','c(:, 2)','srdel(:, 1)','srdel(:, 2)','vmad','xi_p','xi_d','obar(:, 1)','obar(:, 2)')


def assignments(path):
    out={}
    for line in Path(path).read_text().splitlines():
        line=line.split('!')[0]
        if '=' not in line: continue
        key,value=line.split('=',1)
        try: out[re.sub(r'\s+','',key.lower())]=np.array([float(x.strip().replace('D','E').replace('d','e')) for x in value.split(',')])
        except ValueError: continue
    return out


def clust(path):
    lines=Path(path).read_text().splitlines()
    n=int(lines[0].split('=')[1])
    a=np.array(' '.join(lines[1:]).split(),float).reshape(-1,5)
    if len(a)!=n: raise ValueError('Invalid clust atom count')
    return a


def pairs(path):
    return np.array([list(map(int,m)) for m in re.findall(r'ijpair\(\d+,\s*:\)\s*=\s*(\d+),\s*(\d+)',Path(path).read_text(),re.I)])


def keyed_output(path):
    a=np.loadtxt(path,ndmin=2)
    d={tuple(row[:5]):row[5:-1] for row in a}
    if len(d)!=len(a): raise ValueError('Duplicate output bond key')
    return d


def compare_outputs(a,b):
    if set(a)!=set(b): raise ValueError('Directed-bond keys differ')
    aa=np.array([a[k] for k in sorted(a)])
    bb=np.array([b[k] for k in sorted(a)])
    return {'even':residual(aa,bb),'odd':residual(aa,-bb)}


def inventory(root):
    out={'cases':{},'domains':{}}
    lists=[]
    atoms={}
    for case in CASES:
        p=root/case/'exchange_jobs/Mn1'
        cluster=clust(p/'clust')
        pair=pairs(p/'input.nml')
        lists.append((pair,cluster[pair[:,1]-1,:3]-cluster[pair[:,0]-1,:3]))
        d=np.loadtxt(p/'dijso.out'); parts=np.loadtxt(p/'dijparts.out')
        mask=(d[:,0]==1)&(d[:,1]==2)
        atom={label:assignments(p/(label+'_out.nml')) for label in ('Mn1','Mn2','F1','F2')}
        atoms[case]=atom
        mirror={}
        for field in FIELDS:
            key=re.sub(r'\s+','',field)
            mirror[field]=residual(atom['F1'][key],atom['F2'][key])
        out['cases'][case]={'files':sum(f.is_file() for f in (root/case).rglob('*')),
            'pair_count':len(pair),'intersublattice_count':int(mask.sum()),'cluster_count':len(cluster),
            'sum_D_mRy':d[mask,5:8].sum(axis=0).tolist(),
            'sum_CC_mRy':parts[mask,5:8].sum(axis=0).tolist(),
            'sum_SC_mRy':parts[mask,8:11].sum(axis=0).tolist(),
            'parts_sum_max_residual_mRy':float(np.max(abs(parts[:,5:8]+parts[:,8:11]-d[:,5:8]))),
            'F_mirror_parameter_residuals':mirror,
            'moments':{label:a['mom'].tolist() for label,a in atom.items()},
            'lattice_sha256':hashlib.sha256((p/'lattice.nml').read_bytes()).hexdigest(),
            'inputs_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in p.glob('*.nml')},
            'bond_maps':{}}
        for name,diag in [('C2x',[1,-1,-1]),('C2zT',[-1,-1,1]),('I',[-1,-1,-1])]:
            mapping=[]
            for (i,j),r in zip(pair,lists[-1][1]):
                target=r*np.array(diag)
                candidates=[(int(ij[0]),int(ij[1])) for ij,rr in zip(pair,lists[-1][1])
                            if cluster[ij[0]-1,3]==cluster[i-1,3] and cluster[ij[1]-1,3]==cluster[j-1,3]
                            and np.max(abs(rr-target))<1e-7]
                mapping.append({'source_indices':[int(i),int(j)],'r':r.tolist(),'target_indices':candidates})
            out['cases'][case]['bond_maps'][name]=mapping
    out['pair_lists_equal']=all(np.array_equal(lists[0][0],x[0]) for x in lists[1:])
    out['relative_bond_max_difference']=max(float(np.max(abs(lists[0][1]-x[1]))) for x in lists[1:])
    for plus,minus in [('100','1-00'),('110','1-1-0')]:
        p=root/plus/'exchange_jobs/Mn1';m=root/minus/'exchange_jobs/Mn1'
        d=keyed_output(p/'dijso.out');dm=keyed_output(m/'dijso.out')
        key=f'{plus}/{minus}'
        out['domains'][key]={'Dx':compare_outputs({k:v[:1] for k,v in d.items()},{k:v[:1] for k,v in dm.items()}),
                              'Dz':compare_outputs({k:v[2:] for k,v in d.items()},{k:v[2:] for k,v in dm.items()}),
                              'atom_even_residuals':{},'moment_odd_residuals':{}}
        for label in atoms[plus]:
            a,b=atoms[plus][label],atoms[minus][label]
            out['domains'][key]['atom_even_residuals'][label]={field:residual(a[re.sub(r'\s+','',field)],b[re.sub(r'\s+','',field)]) for field in FIELDS}
            out['domains'][key]['moment_odd_residuals'][label]=residual(a['mom'],-b['mom'])
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--regression',type=Path,help='Compare dijso/dijparts in two single job directories, at printed precision')
    args=parser.parse_args()
    if args.regression:
        out={}
        for file in ('dijso.out','dijparts.out'):
            a,b=keyed_output(args.root/file),keyed_output(args.regression/file)
            out[file]=compare_outputs(a,b)['even']
            out[file]['printed_values_identical']=all(np.array_equal(a[k],b[k]) for k in a)
    else: out=inventory(args.root)
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(args.output)


if __name__=='__main__': main()
