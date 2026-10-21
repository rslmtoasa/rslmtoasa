#!/usr/bin/env python3
"""Compact PRIVATE before/after analysis of complete frozen axial exchange runs.

Keeps bond values, peak trace decompositions and norm residuals, not raw matrices.
Only read completed runs: output files alone are not completion evidence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from analyze import (read_dump, stage_report, h_report, bond_key, matrix,
                     same_grid, covariance_report)
from symmetry import SECTORS,BARE_T,MIRROR,P,UT


def complete(root,cases=('100','1-00')):
    manifest=json.loads((root/'manifest.json').read_text())
    for case in cases:
        item=manifest['cases'][case]
        if item['returncode'] or not item['potentials_unchanged']:
            raise ValueError('Incomplete run or changed potentials')
        # Check immutable calculation inputs after execution as well as on copy.
        for name,expected in item['copied_hashes_before_pair_selection'].items():
            if name.endswith('.nml') or name in ('clust','sbar.','view.sbar.','map.'):
                actual=hashlib.sha256((root/case/name).read_bytes()).hexdigest()
                if actual!=expected: raise ValueError(f'Input changed during execution: {case}/{name}')
    if manifest['one_pair'] or manifest['pair_count']:
        raise ValueError('Full 24-pair runs required')
    return manifest


def maximum(store,key,value):
    store[key]=max(store.get(key,0.0),float(value))


def directory(path,neel='100'):
    data=[read_dump(p) for p in sorted(path.glob('gf-audit-gf-*.dat'))]
    if len(data)!=24 or len({bond_key(d) for d in data})!=24:
        raise ValueError('Expected 24 unique directed bonds')
    out={'mirror_fro_max':{},'leakage_fro_max':{},'allowed_fro_max':{},'trace_block_abs_max':{},
         'H_max':{},'H_fro_max':{},'trace_peaks':{},'fixed_probes':{},'reconstruction_max':{},'contraction_error_max':0.0}
    out['full_graph']={k:float(v) for k,v in (line.split() for line in (path/'gf-audit-graph.dat').read_text().splitlines())}
    for d in data:
        s=stage_report(d)
        for e in s['energies']:
            for stage in ('projection','raw_pauli','delta'):
                for v in e[stage].values(): maximum(out['reconstruction_max'],stage,v['max'])
            for k,v in e['mirror'].items(): maximum(out['mirror_fro_max'],k,v['fro'])
            for k,v in e['leakage'].items():
                maximum(out['leakage_fro_max'],k,v['fro'])
                maximum(out['allowed_fro_max'],k,np.sqrt(max(0,e['mirror'][k]['scale']**2-v['fro']**2)))
            if d['pair']==1 and e['index'] in (513,577):
                out['fixed_probes'][str(e['index'])]={'energy':e['energy'],'trace':e['trace']}
            if e['energy']<=d['fermi']:
                for kind,v in e['trace'].items():
                    for k,x in v['blocks'].items(): maximum(out['trace_block_abs_max'],kind+':'+k,abs(x))
                    if kind not in out['trace_peaks'] or abs(v['total'])>abs(out['trace_peaks'][kind]['total']):
                        out['trace_peaks'][kind]={'pair':d['pair'],'index':e['index'],'energy':e['energy'],**v}
                    out['contraction_error_max']=max(out['contraction_error_max'],v['contraction_residual'])
    hfiles=list(path.glob('gf-audit-h-*.dat'))
    if len(hfiles)!=24: raise ValueError('Expected 24 complete H captures')
    for f in hfiles:
        h=h_report(read_dump(f))
        entries={'effective_pair':h['pair_hermiticity'],
                 'effective_onsite':max(h['onsite_hermiticity'].values(),key=lambda v:v['max'])}
        for k,v in h['component_hermiticity'].items():
            entries[k+'_pair']=v['pair']
            entries[k+'_onsite']=max(v['onsites'].values(),key=lambda v:v['max'])
        for k,v in entries.items():
            maximum(out['H_max'],k,v['max']);maximum(out['H_fro_max'],k,v['fro'])
        maximum(out['H_max'],'direct_mzT',max(x['mirror']['max'] for x in h['columns']))
        if h['missing_mirror_atoms']: raise ValueError('Missing H mirror partner')
        if h['pair']==1: out['first_pair_H']=entries
    # C2x does not stabilize the [110] domain; do not mislabel it a symmetry.
    cov=[x for x in covariance_report(data,neel=neel)
         if x.get('status')!='not a symmetry of this domain']
    if any('status' in x for x in cov): raise ValueError('Missing GF spatial partner')
    out['spatial_raw_fro_max']={}
    for x in cov:
        maximum(out['spatial_raw_fro_max'],x['operation'],max(e['raw']['fro'] for e in x['energies']))
    d=np.loadtxt(path/'dijso.out',ndmin=2);parts=np.loadtxt(path/'dijparts.out',ndmin=2)
    if len(d)!=24 or len(parts)!=24 or not np.array_equal(d[:,:5],parts[:,:5]):
        raise ValueError('Output bond ordering/count mismatch')
    mask=(d[:,0]==1)&(d[:,1]==2)
    out['sum_D_mRy']=d[mask,5:8].sum(axis=0).tolist()
    out['sum_CC_mRy']=parts[mask,5:8].sum(axis=0).tolist()
    out['sum_SC_mRy']=parts[mask,8:11].sum(axis=0).tolist()
    out['bonds']=[{'key':row[:5].tolist(),'D':row[5:8].tolist(),
                   'CC':p[5:8].tolist(),'SC':p[8:11].tolist()} for row,p in zip(d,parts)]
    out['Jijso'] = np.loadtxt(path/'jijso.out',ndmin=2)[:,5].tolist()
    return out,data


def parity(plus,minus):
    bykey={bond_key(d):d for d in minus}
    sums={}
    def add(k,a,b):
        v=sums.setdefault(k,np.zeros(5))
        v[:3]+=np.array([np.linalg.norm(b-a)**2,np.linalg.norm(b+a)**2,np.linalg.norm(a)**2])
        v[3]=max(v[3],np.max(abs(b-a)));v[4]=max(v[4],np.max(abs(b+a)))
    for d in plus:
        other=bykey[bond_key(d)];same_grid(d,other)
        if (d['ti'],d['tj'])!=(1,2): continue
        for n in sorted(d['integrands']):
            if d['integrands'][n][0]<=d['fermi']:
                for q,k in enumerate(('CCx','CCy','CCz','SCx','SCy','SCz')):
                    add('integrand_'+k,d['integrands'][n][q+1],other['integrands'][n][q+1])
        for n in sorted({i for i,t in d['matrices']}):
            if d['integrands'][n][0]>d['fermi']: continue
            for direction,reverse in (('ij','ji'),('ji','ij')):
                expected=UT@matrix(d,n,'g'+reverse).T@UT.conj().T
                add('raw_T_'+direction,expected,matrix(other,n,'g'+direction))
            for k in SECTORS:
                for direction in ('ij','ji'):
                    a=matrix(d,n,'g'+k+direction);b=matrix(other,n,'g'+k+direction)
                    add(k+direction,a,BARE_T[k]*b)
                    mask=P[:,None]*P[None,:]==MIRROR[k]
                    for label,m in [('allowed',mask),('forbidden',~mask)]:
                        add(k+direction+'_'+label,a[m],b[m])
    return {k:{'difference_fro':float(np.sqrt(v[0])),'sum_fro':float(np.sqrt(v[1])),
               'reference_fro':float(np.sqrt(v[2])),
               'difference_relative':float(np.sqrt(v[0]/v[2])) if v[2] else None,
               'sum_relative':float(np.sqrt(v[1]/v[2])) if v[2] else None,
               'difference_max':float(v[3]),'sum_max':float(v[4])} for k,v in sums.items()}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',type=Path,required=True)
    p.add_argument('--reciprocal',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();out={}
    for label,root in [('original',args.original),('reciprocal',args.reciprocal)]:
        manifest=complete(root);group={'binary_sha256':manifest['binary_sha256'],'cases':{}};data={}
        for item in manifest['cases'].values():
            if item['sbar_override_loaded']!=(label=='reciprocal'):
                raise ValueError('Incorrect saved-Sbar intervention state')
            if item['changed_on_copy']!=(['sbar.'] if label=='reciprocal' else []):
                raise ValueError('Unexpected input intervention')
        for case in ('100','1-00'): group['cases'][case],data[case]=directory(root/case)
        group['occupied_domain_residuals']=parity(data['100'],data['1-00'])
        for kind in ('D','CC','SC'):
            a=np.array(group['cases']['100']['sum_'+kind+'_mRy']);b=np.array(group['cases']['1-00']['sum_'+kind+'_mRy'])
            group[kind+'_odd_mRy']=((a-b)/2).tolist();group[kind+'_even_mRy']=((a+b)/2).tolist()
        out[label]=group
        print(label,group['D_odd_mRy'],flush=True)
    if out['original']['binary_sha256']!=out['reciprocal']['binary_sha256']:
        raise ValueError('Different executable in before/after runs')
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__=='__main__': main()
