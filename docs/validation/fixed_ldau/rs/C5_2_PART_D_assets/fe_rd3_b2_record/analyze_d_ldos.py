#!/usr/bin/env python3
"""Descriptive max Fe d-LDOS comparison from one-step full RS Green traces.
Requires same absolute energy grid; deliberately does not interpolate or set a gate.
"""
import argparse,csv,json
from decimal import Decimal,getcontext
from pathlib import Path
getcontext().prec=50
PI=Decimal('3.1415926535897932384626433832795028841971693993751')

def rows(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))
def parse(path):
    rs=rows(path)
    maps={}
    greens={}
    pairs_by_grid={}
    for r in rs:
        if r['stage']!='RS_GREEN':continue
        if r['kind']=='basis_map':
            site=int(r['site']); idx=int(r['row'])
            if idx in maps.setdefault(site,{}):raise ValueError(f'duplicate basis map {site} {idx}')
            maps[site][idx]=(int(r['l']),int(r['col']),int(r['spin']))
        elif r['kind']=='green_function':
            site,k=int(r['site']),int(r['k'])
            energy=Decimal(r['energy'])
            i,j=int(r['row']),int(r['col'])
            bucket=greens.setdefault((site,k),{'energy':energy,'diag':{}})
            if bucket['energy']!=energy:raise ValueError(f'energy mismatch inside grid at {site}/{k}')
            pairs=pairs_by_grid.setdefault((site,k),set())
            pair=(i,j)
            if pair in pairs:raise ValueError(f'duplicate Green row at {site}/{k}/{pair}')
            pairs.add(pair)
            if i==j:
                if i in bucket['diag']:raise ValueError(f'duplicate diagonal G at {site}/{k}/{i}')
                bucket['diag'][i]=Decimal(r['im'])
    if not maps or not greens:raise ValueError('missing RS_GREEN basis map or matrix')
    if set(maps)!={0}:raise ValueError(f'expected single impurity/bulk output site 0, got {sorted(maps)}')
    mapping=maps[0]
    n=len(mapping)
    if set(mapping)!=set(range(n)):raise ValueError('basis row indices not contiguous')
    expected={(l,m,s) for l in range(3) for m in range(-l,l+1) for s in (0,1)}
    if set(mapping.values())!=expected:raise ValueError('basis map is not a complete l<=2 spin basis')
    d_by_spin={s:[i for i,(l,m,sp) in mapping.items() if l==2 and sp==s] for s in (0,1)}
    energy_grid={}
    for (site,k),g in greens.items():
        if site!=0:continue
        if set(g['diag'])!=set(range(n)):raise ValueError(f'missing Green diagonal at k={k}')
        if pairs_by_grid[(0,k)]!={(i,j) for i in range(n) for j in range(n)}:raise ValueError(f'incomplete full Green matrix at k={k}')
        energy_grid[k]=g['energy']
    if set(energy_grid)!=set(range(len(energy_grid))):raise ValueError('energy indices not contiguous')
    dldos={}
    for k,energy in energy_grid.items():
        g=greens[(0,k)]['diag']
        dldos[(k,0)]=-sum((g[i] for i in d_by_spin[0]),Decimal(0))/PI
        dldos[(k,1)]=-sum((g[i] for i in d_by_spin[1]),Decimal(0))/PI
    return energy_grid,dldos,{'basis_dimension':n,'energy_points':len(energy_grid),'d_orbitals_per_spin':5,'source_convention':'bands.f90 calculate_dos: dosial=-aimag(Gii)/pi'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('bulk_csv');ap.add_argument('impurity_csv');ap.add_argument('--output',required=True);a=ap.parse_args()
    e0,d0,m0=parse(a.bulk_csv);e1,d1,m1=parse(a.impurity_csv)
    out={'schema':'c52-fe-rd3-dldos-descriptive-v1','bulk_csv':str(Path(a.bulk_csv).resolve()),'impurity_csv':str(Path(a.impurity_csv).resolve()),'energy_grid_equal':e0==e1,'bulk_shape':m0,'impurity_shape':m1,'interpolation':'none'}
    if e0!=e1:
        out.update(status='NOT_COMPUTED_DIFFERENT_ENERGY_GRIDS',bulk_energy_points=len(e0),impurity_energy_points=len(e1))
    else:
        byspin={}
        for s,label in ((0,'spin1'),(1,'spin2')):
            diffs=[abs(d1[(k,s)]-d0[(k,s)]) for k in e0]
            ix=max(range(len(diffs)),key=diffs.__getitem__)
            byspin[label]={'max_abs_difference_per_Ry':str(diffs[ix]),'energy_Ry':str(e0[ix]),'bulk_d_ldos_per_Ry':str(d0[(ix,s)]),'impurity_d_ldos_per_Ry':str(d1[(ix,s)])}
        out.update(status='COMPUTED_DESCRIPTIVE_NO_GATE',max_abs_difference_over_spin_energy=max(Decimal(v['max_abs_difference_per_Ry']) for v in byspin.values()).__str__(),per_spin=byspin,energy_points=len(e0))
    Path(a.output).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':main()
