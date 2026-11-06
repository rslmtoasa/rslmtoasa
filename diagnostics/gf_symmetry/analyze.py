#!/usr/bin/env python3
"""Analyze private opt-in dumps. JSON output may contain unpublished results."""
import argparse
import json
from pathlib import Path
import numpy as np
from symmetry import (M, UH, UT, MIRROR, BARE_T, SECTORS, OPS, sectors, pauli,
                      residual, trace_parts, leakage)


def read_dump(path):
    result = {'matrices': {}, 'integrands': {}, 'atoms': {}}
    with Path(path).open() as f:
        for line in f:
            v = line.split()
            if not v or v[0].startswith('#'):
                continue
            if v[0] == 'M':
                index, energy, tag, n = int(v[1]), float(v[2]), v[3], int(v[4])
                if n not in (9,18):
                    raise ValueError(f'Invalid matrix dimension in {path}')
                rows = [np.fromstring(next(f), sep=' ') for _ in range(n)]
                if any(len(row) != 2*n for row in rows):
                    raise ValueError(f'Truncated matrix in {path}')
                a = np.array(rows)
                key = index,tag
                if key in result['matrices']:
                    raise ValueError(f'Duplicate matrix {key} in {path}')
                result['matrices'][key] = (energy,a[:,::2]+1j*a[:,1::2])
            elif v[0] == 'P':
                result['pair'],result['i'],result['j'],result['ti'],result['tj'] = map(int,v[1:6])
                result['r'] = np.array(v[6:9],float)
            elif v[0] == 'I':
                result['integrands'][int(v[1])] = np.array(v[2:],float)
            elif v[0] == 'F':
                result['fermi'],result['local_axis'] = float(v[1]),v[2]
            elif v[0] == 'A':
                result['atoms'][int(v[1])] = (int(v[2]),np.array(v[3:6],float))
            elif v[0] == 'S':
                result['pair'],result['a'],result['b'] = int(v[1]),float(v[2]),float(v[3])
            elif v[0] == 'B':
                result['i'],result['j'] = map(int,v[1:3])
            else:
                raise ValueError(f'Unknown record in {path}: {line}')
    for _,a in result['matrices'].values():
        if not np.isfinite(a).all():
            raise ValueError(f'Nonfinite matrix in {path}')
    return result


def bond_key(d):
    # Only matching keys, not matrix comparison, uses rounded displacement.
    return d['ti'], d['tj'], *np.round(d['r'],9)


def matrix(d,n,tag):
    return d['matrices'][n,tag][1]


def same_grid(a,b):
    if set(a['integrands']) != set(b['integrands']):
        raise ValueError('Energy index sets differ')
    for n in a['integrands']:
        if abs(a['integrands'][n][0]-b['integrands'][n][0]) > 1e-13:
            raise ValueError('Energy values differ; do not compare unmatched grids')
    if abs(a['fermi']-b['fermi']) > 1e-13:
        raise ValueError('Fermi energies differ; integrated parity needs separate analysis')


def stage_report(d,atol=1e-12,rtol=1e-9):
    out = {'pair':d['pair'],'indices':[d['i'],d['j']], 'displacement':d['r'].tolist(),
           'local_axis':d['local_axis'], 'energies':[], 'first_observed_failure':None}
    stages = {}
    for n in sorted({i for i,tag in d['matrices']}):
        energy = d['matrices'][n,'gij'][0]
        if abs(energy-d['integrands'][n][0]) > 1e-13:
            raise ValueError('Snapshot/integrand energy mismatch')
        row = {'index':n,'energy':energy,'mirror':{},'projection':{},'raw_pauli':{},'trace':{},'leakage':{}}
        raw = {direction:matrix(d,n,'g'+direction) for direction in ('ij','ji')}
        for direction,other in [('ij','ji'),('ji','ij')]:
            row['mirror']['raw_'+direction] = residual(raw[direction],UH@raw[other].T@UH.conj().T)
            calc = sectors(raw[direction],raw[other])
            for name in SECTORS:
                g = matrix(d,n,'g'+name+direction)
                row['mirror'][name+direction] = residual(M@g@M,MIRROR[name]*g)
                row['leakage'][name+direction] = leakage(g,MIRROR[name])
                row['projection'][name+direction] = residual(g,calc[name])
            for mu,g in zip('0xyz',pauli(raw[direction])):
                name = {'0':'gnmag','x':'gx','y':'gy','z':'gz'}[mu]
                tag = name[0] + ('i' if direction=='ij' else 'j') + name[1:]
                row['raw_pauli'][mu+direction] = residual(matrix(d,n,tag),g)
        di,dj = matrix(d,n,'Delta_i'),matrix(d,n,'Delta_j')
        row['delta'] = {'i':residual(M@di@M,di),'j':residual(M@dj@M,dj)}
        for k,mu in enumerate('xyz'):
            for kind,left,right in [('CC','01',mu+'0'),('SC','00',mu+'1')]:
                item = trace_parts(di,dj,matrix(d,n,'g'+left+'ij'),matrix(d,n,'g'+right+'ji'),MIRROR[left],MIRROR[right])
                item['production'] = float(d['integrands'][n][1+k+(3 if kind=='SC' else 0)])
                item['contraction_residual'] = abs(item['total']-item['production'])
                row['trace'][kind+mu] = item
        # A downstream failed sector cannot establish a failure at an earlier stage.
        for stage,values in [('raw_GF_mirror',{k:v for k,v in row['mirror'].items() if k.startswith('raw')}),
                             ('Pauli_decomposition',row['raw_pauli']),('sector_projection',row['projection']),
                             ('sector_mirror',{k:v for k,v in row['mirror'].items() if not k.startswith('raw')}),
                             ('Delta_mirror',row['delta'])]:
            for label,r in values.items():
                if r['fro'] > atol+rtol*r['scale']:
                    stages.setdefault(stage,[]).append({'energy_index':n,'quantity':label,'residual':r})
        for label,item in row['trace'].items():
            if item['contraction_residual'] > atol+rtol*max(abs(item['total']),abs(item['production'])):
                stages.setdefault('contraction',[]).append({'energy_index':n,'quantity':label,'absolute':item['contraction_residual']})
        out['energies'].append(row)
    order = ['raw_GF_mirror','Pauli_decomposition','sector_projection','sector_mirror','Delta_mirror','contraction']
    out['failures'] = stages
    out['first_observed_failure'] = next(({'stage':k,'example':stages[k][0]} for k in order if k in stages),None)
    out['scope'] = 'Sampled energies only; a failed raw GF does not identify its upstream cause.'
    # Aggregates supplement, never replace, per-energy matrix residuals.
    out['sector_maxima'] = {label:{metric:max(row['mirror'][label][metric] for row in out['energies'])
                                   for metric in ('max','mean','fro','relative')}
                            for label in out['energies'][0]['mirror']} if out['energies'] else {}
    return out


def domain_report(plus,minus):
    if bond_key(plus) != bond_key(minus):
        raise ValueError('Domain comparison requires matched ordered bond')
    same_grid(plus,minus)
    keys_a,keys_b = set(plus['matrices']),set(minus['matrices'])
    if keys_a != keys_b:
        raise ValueError('Sampled matrix sets differ')
    rows=[]
    for n in sorted({i for i,tag in keys_a}):
        row={'index':n,'raw':{},'sectors':{},'delta':{},
             'real_integrands_even':residual(minus['integrands'][n][1:],plus['integrands'][n][1:])}
        for d,o in [('ij','ji'),('ji','ij')]:
            row['raw'][d]=residual(matrix(minus,n,'g'+d),UT@matrix(plus,n,'g'+o).T@UT.conj().T)
            for name in SECTORS:
                row['sectors'][name+d]=residual(matrix(minus,n,'g'+name+d),BARE_T[name]*matrix(plus,n,'g'+name+d))
        for tag in ('Delta_i','Delta_j'):
            row['delta'][tag]=residual(matrix(minus,n,tag),matrix(plus,n,tag))
        rows.append(row)
    keys_i = plus['integrands']
    return {'pair':plus['pair'],'energies':rows,
            'real_integrands_even':residual(np.array([minus['integrands'][n][1:] for n in sorted(keys_i)]),
                                           np.array([plus['integrands'][n][1:] for n in sorted(keys_i)]))}


def covariance_report(datasets, neel='100'):
    bykey={bond_key(d):d for d in datasets}
    if len(bykey)!=len(datasets):
        raise ValueError('Ambiguous duplicate ordered bonds')
    rows=[]
    for a in datasets:
        for name,(rot,u,anti) in OPS.items():
            if name == 'C2x' and neel != '100':
                rows.append({'pair':a['pair'],'operation':name,'status':'not a symmetry of this domain'})
                continue
            key=(a['ti'],a['tj'],*np.round(rot@a['r'],9))
            if key not in bykey:
                rows.append({'pair':a['pair'],'operation':name,'status':'mapped bond missing'})
                continue
            b=bykey[key]
            same_grid(a,b)
            row={'pair':a['pair'],'mapped_pair':b['pair'],'indices':[a['i'],a['j']],
                 'mapped_indices':[b['i'],b['j']],'operation':name,'energies':[]}
            for n in sorted({i for i,tag in a['matrices']}):
                g=matrix(a,n,'gji' if anti else 'gij')
                if anti: g=g.T
                row['energies'].append({'index':n,'raw':residual(matrix(b,n,'gij'),u@g@u.conj().T)})
            rows.append(row)
    return rows


def h_columns(d):
    zeros=np.zeros((18,18),complex)
    cols={src:{} for src in (d['i'],d['j'])}
    for atom in d['atoms']:
        p=d['matrices'].get((atom,'Hp'),(0,zeros))[1]
        m=d['matrices'].get((atom,'Hm'),(0,zeros))[1]
        for src,sign in [(d['i'],1),(d['j'],-1)]:
            h=d['a']*(p+sign*m)/np.sqrt(2)
            if atom==src: h=h+d['b']*np.eye(18)
            cols[src][atom]=h
    return cols


def h_report(d):
    cols=h_columns(d)
    z=d['atoms'][d['i']][1][2]
    out={'pair':d['pair'],'columns':[],'missing_mirror_atoms':[],
         'scope':'Full output rows of two actual H columns, including HoH/SOC; not the entire Hamiltonian.'}
    atom_ids = list(d['atoms'])
    coords = np.array([d['atoms'][k][1] for k in atom_ids])
    for src,blocks in cols.items():
        for atom,h in blocks.items():
            _,r=d['atoms'][atom]
            reflected=r.copy(); reflected[2]=2*z-r[2]
            hits=[atom_ids[k] for k in np.flatnonzero(np.max(abs(coords-reflected),axis=1)<1e-7)]
            if len(hits)!=1:
                if np.linalg.norm(h)>0: out['missing_mirror_atoms'].append(atom)
                continue
            mapped=hits[0]
            if np.linalg.norm(h)+np.linalg.norm(blocks[mapped])==0: continue
            out['columns'].append({'source':src,'row':atom,'mirror_row':mapped,
                                    'mirror':residual(h,UH@blocks[mapped].conj()@UH.conj().T)})
    i,j=d['i'],d['j']
    out['pair_hermiticity']=residual(cols[j][i],cols[i][j].conj().T)
    if any(tag=='Op' for _,tag in d['matrices']):
        components={name:{src:{} for src in (i,j)} for name in ('hoh','enu','soc','bare')}
        zeros=np.zeros((18,18),complex)
        for atom in d['atoms']:
            for src,sign in [(i,1),(j,-1)]:
                for name,prefix in [('hoh','O'),('enu','E'),('soc','S')]:
                    plus=d['matrices'].get((atom,prefix+'p'),(0,zeros))[1]
                    minus=d['matrices'].get((atom,prefix+'m'),(0,zeros))[1]
                    components[name][src][atom]=(plus+sign*minus)/np.sqrt(2)
                components['bare'][src][atom]=cols[src][atom]+components['hoh'][src][atom]-components['enu'][src][atom]-components['soc'][src][atom]
        out['component_hermiticity']={name:{'pair':residual(c[j][i],c[i][j].conj().T),
                                       'onsites':{str(src):residual(c[src][src],c[src][src].conj().T) for src in c}}
                                       for name,c in components.items()}
    out['onsite_hermiticity']={str(src):residual(cols[src][src],cols[src][src].conj().T) for src in cols}
    return out



def domain_h_report(plus,minus):
    if plus['i']!=minus['i'] or plus['j']!=minus['j']:
        raise ValueError('Hamiltonian source indices differ')
    a,b=h_columns(plus),h_columns(minus)
    rows=[]
    for src,blocks in a.items():
        for atom,h in blocks.items():
            if atom not in b[src]: raise ValueError('Hamiltonian atom sets differ')
            if np.max(abs(plus['atoms'][atom][1]-minus['atoms'][atom][1]))>1e-13:
                raise ValueError('Hamiltonian coordinates differ')
            if np.linalg.norm(h)+np.linalg.norm(b[src][atom])==0: continue
            rows.append({'source':src,'row':atom,'bare_T':residual(b[src][atom],UT@h.conj()@UT.conj().T)})
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--neel',choices=['100','110'],default='100',help='Reference domain direction for surviving operations')
    parser.add_argument('--minus',type=Path,help='Time-reversed domain dump directory')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data=[read_dump(p) for p in sorted(args.directory.glob('gf-audit-gf-*.dat'))]
    if not data: parser.error('No GF dumps found; exchange execution is required')
    output={'pairs':[stage_report(d) for d in data], 'covariance':covariance_report(data,args.neel),
            'hamiltonian':[h_report(read_dump(p)) for p in sorted(args.directory.glob('gf-audit-h-*.dat'))]}
    if args.minus:
        minus={bond_key(d):d for d in map(read_dump,args.minus.glob('gf-audit-gf-*.dat'))}
        if set(minus)!={bond_key(d) for d in data}: parser.error('Domain ordered-bond sets differ')
        output['bare_T']=[domain_report(d,minus[bond_key(d)]) for d in data]
        hp={d['pair']:d for d in map(read_dump,args.directory.glob('gf-audit-h-*.dat'))}
        hm={d['pair']:d for d in map(read_dump,args.minus.glob('gf-audit-h-*.dat'))}
        output['hamiltonian_bare_T']={str(k):domain_h_report(hp[k],hm[k]) for k in hp.keys()&hm.keys()}
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(f'Analyzed {len(data)} pairs; private detailed report written to {args.output}')


if __name__=='__main__':
    main()
