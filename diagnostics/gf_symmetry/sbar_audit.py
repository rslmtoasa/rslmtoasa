#!/usr/bin/env python3
"""Check screened-structure-constant reciprocity before constructing any GF.

Reads original unformatted sbar. records (little-endian GNU/Intel 4-byte
record markers, nine float64 entries per record), using view.sbar. only for
vector/type labels. The matrices come from the binary, not its rounded preview.
"""
import argparse
import json
from pathlib import Path
import re
import struct
import numpy as np
from fixtures import assignments,clust
from symmetry import V,SIGMA,MIRROR,OPS,UH,sectors,residual,leakage
from analyze import read_dump,h_columns


def read_sbar(job):
    rows=[]
    with (job/'sbar.').open('rb') as f:
        while head:=f.read(4):
            if len(head)!=4: raise ValueError('Truncated record marker')
            n=struct.unpack('<I',head)[0]
            if n!=72: raise ValueError('Unsupported sbar format; expected 9 float64 row')
            chunk=f.read(n)
            if len(chunk)!=n or f.read(4)!=head: raise ValueError('Invalid Fortran record')
            rows.append(np.frombuffer(chunk,dtype='<f8'))
    blocks=np.array(rows).reshape(-1,9,9)
    text=(job/'view.sbar.').read_text()
    labels=re.findall(r'VECTOR=\s*(\S+)\s+(\S+)\s+(\S+)\s+ICLUS=\s*(\d+)',text)
    if len(labels)!=len(blocks): raise ValueError('Preview/binary block counts differ')
    cluster=clust(job/'clust')
    alat=float(re.search(r'^\s*alat\s*=\s*(\S+)',(job/'input.nml').read_text(),re.M|re.I).group(1))
    coords=cluster[:,:3]*alat
    bases={t:coords[np.flatnonzero(cluster[:,3]==t)[0]] for t in np.unique(cluster[:,3]).astype(int)}
    matrices={}
    for block,label in zip(blocks,labels):
        r=np.array(label[:3],float); ti=int(label[3])
        distances=np.linalg.norm(coords-(bases[ti]+r),axis=1)
        hit=int(distances.argmin())
        if distances[hit]>2e-6: raise ValueError('Cannot identify destination type from rounded vector label')
        tj=int(cluster[hit,3])
        key=ti,tj,*r
        if key in matrices: raise ValueError('Duplicate structure-constant key')
        matrices[key]=block
    return matrices


def width_spinor(atom):
    channels=[]
    for spin in (1,2):
        w=atom[f'width_band(:,{spin})']
        channels.append(np.diag(np.repeat(w,[1,3,5])))
    w0=(channels[0]+channels[1])/2
    w1=(channels[0]-channels[1])/2
    return np.kron(SIGMA[0],w0)+sum(np.kron(s*n,w1) for s,n in zip(SIGMA[1:],atom['mom']))


def bare_block(s,ai,aj):
    u=np.kron(np.eye(2),V)
    return u.conj().T@width_spinor(ai)@np.kron(np.eye(2),s.T)@width_spinor(aj)@u


def captured_bare(d,src,row):
    zero=np.zeros((18,18),complex)
    sign=1 if src==d['i'] else -1
    full=h_columns(d)[src][row]
    parts={}
    for name,prefix in [('hoh','O'),('enu','E'),('soc','S')]:
        p=d['matrices'].get((row,prefix+'p'),(0,zero))[1]
        m=d['matrices'].get((row,prefix+'m'),(0,zero))[1]
        parts[name]=(p+sign*m)/np.sqrt(2)
    return full+parts['hoh']-parts['enu']-parts['soc']


def audit(job,h_dump=None,neel="100"):
    matrices=read_sbar(job)
    # This fixture's label/type mapping is explicitly supplied by input.nml.
    text=(job/'input.nml').read_text()
    labelmap={int(i):label for i,label in re.findall(r'label\((\d+)\)\s*=\s*[\'\"]([^\'\"]+)[\'\"]',text,re.I)}
    atoms={t:assignments(job/(label+'_out.nml')) for t,label in labelmap.items()}
    rows=[]
    for key,s in matrices.items():
        ti,tj,*r=key
        reverse=tj,ti,*(-np.array(r))
        if reverse not in matrices:
            rows.append({'key':list(key),'status':'reverse vector missing'})
            continue
        sr=matrices[reverse]
        hi,hj=bare_block(s,atoms[ti],atoms[tj]),bare_block(sr,atoms[tj],atoms[ti])
        item={'key':list(key),'sbar_reciprocity':residual(s,sr.T),
              'bare_h_hermiticity':residual(hi,hj.conj().T)}
        if ti==tj and np.linalg.norm(r)==0:
            # Onsite radial additions are Hermitian and absent from bare_block.
            item['scope']='structure-constant contribution only'
        for name,a in sectors(hi,hj).items():
            item.setdefault('bare_sector_leakage',{})[name]=leakage(a,MIRROR[name])
        rows.append(item)
    out={'blocks':rows,'missing_reverse_count':sum('status' in row for row in rows)}
    operations=dict(OPS)
    operations['mzT']=(np.diag([1,1,-1]),UH,True)
    out['bare_h_spatial']=[]
    for key,s in matrices.items():
        ti,tj,*r=key
        if ti not in (1,2) or tj not in (1,2) or r[2]!=0: continue
        hi=bare_block(s,atoms[ti],atoms[tj])
        for name,(rot,u,anti) in operations.items():
            if name=='C2x' and neel!='100': continue
            target=ti,tj,*np.round(rot@np.array(r),6)
            if target not in matrices: continue
            ht=bare_block(matrices[target],atoms[ti],atoms[tj])
            expected=u@(hi.conj() if anti else hi)@u.conj().T
            item={'source_key':list(key),'target_key':list(target),'operation':name,
                  'direct_H_covariance':residual(ht,expected)}
            rev=tj,ti,*(-np.array(r))
            if anti and rev in matrices:
                hr=bare_block(matrices[rev],atoms[tj],atoms[ti])
                item['transpose_reciprocity']=residual(ht,u@hr.T@u.conj().T)
            out['bare_h_spatial'].append(item)
    # Full four-type geometry map, including F1 <-> F2 under z reflection.
    # Includes onsite radial terms from shifted_band for the supplied HoH jobs.
    out['bare_h_global_spatial']=[]
    def complete_block(k):
        ti,tj,*r=k
        result=bare_block(matrices[k],atoms[ti],atoms[tj])
        if ti==tj and np.linalg.norm(r)==0:
            channels=[np.diag(np.repeat(atoms[ti][f'shifted_band(:,{q})'],[1,3,5])) for q in (1,2)]
            result=result+np.kron(SIGMA[0],(channels[0]+channels[1])/2)
            result=result+sum(np.kron(spin*n,(channels[0]-channels[1])/2) for spin,n in zip(SIGMA[1:],atoms[ti]['mom']))
        return result
    for key in matrices:
        ti,tj,*r=key
        hi=complete_block(key)
        for name,(rot,u,anti) in operations.items():
            if name=='C2x' and neel!='100': continue
            perm={1:1,2:2,3:4,4:3} if rot[2,2]<0 else {t:t for t in atoms}
            target=perm[ti],perm[tj],*np.round(rot@np.array(r),6)
            if target not in matrices:
                out['bare_h_global_spatial'].append({'source_key':list(key),'operation':name,'status':'mapped block missing'})
                continue
            expected=u@(hi.conj() if anti else hi)@u.conj().T
            out['bare_h_global_spatial'].append({'source_key':list(key),'target_key':list(target),'operation':name,
                                                  'direct_H_covariance':residual(complete_block(target),expected)})
    if h_dump:
        d=read_dump(h_dump)
        cluster=clust(job/'clust')
        alat=float(re.search(r'^\s*alat\s*=\s*(\S+)',text,re.M|re.I).group(1))
        i,j=d['i'],d['j'];r=(cluster[j-1,:3]-cluster[i-1,:3])*alat
        ti,tj=int(cluster[i-1,3]),int(cluster[j-1,3])
        keys=[k for k in matrices if k[:2]==(ti,tj) and np.max(abs(np.array(k[2:])-r))<2e-6]
        reverse=[k for k in matrices if k[:2]==(tj,ti) and np.max(abs(np.array(k[2:])+r))<2e-6]
        if len(keys)!=1 or len(reverse)!=1: raise ValueError('Cannot match H dump pair')
        hi=bare_block(matrices[keys[0]],atoms[ti],atoms[tj])
        hj=bare_block(matrices[reverse[0]],atoms[tj],atoms[ti])
        out['captured_bare_agreement']={'ij':residual(hi,captured_bare(d,j,i)),
                                        'ji':residual(hj,captured_bare(d,i,j))}
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('job',type=Path)
    p.add_argument('--h-dump',type=Path)
    p.add_argument('--neel',choices=['100','110'],default='100')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(audit(args.job,args.h_dump,args.neel),indent=2,allow_nan=False)+'\n')
    print(args.output)


if __name__=='__main__': main()
