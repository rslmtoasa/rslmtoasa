#!/usr/bin/env python3
"""Project a COPY of saved real-basis Sbar onto bond reciprocity.

No production code is changed. Input binary and preview ordering are retained;
offsite blocks are paired with the distinct reverse bond, never with themselves.
Only the audited little-endian, 4-byte marker, 9-float64 record format is accepted.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from sbar_audit import read_sbar

TOLERANCE = 2e-6  # absolute Cartesian units of rounded view.sbar. labels


def reverse_map(keys, tolerance=TOLERANCE):
    if not 0 < tolerance <= TOLERANCE:
        raise ValueError('Coordinate tolerance must be positive and at most 2e-6')
    partners=[]
    for ti,tj,*r in keys:
        hits=[q for q,k in enumerate(keys) if k[:2]==(tj,ti)
              and np.max(np.abs(np.array(k[2:])+r))<=tolerance]
        if len(hits)!=1:
            raise ValueError(f'Missing or ambiguous reverse for {(ti,tj,*r)}: {hits}')
        partners.append(hits[0])
    if any(partners[q]!=p for p,q in enumerate(partners)):
        raise ValueError('Reverse map is not an involution')
    if any(p==q and np.any(np.array(keys[p][2:])!=0) for p,q in enumerate(partners)):
        raise ValueError('Only an exactly onsite key may reverse to itself')
    return partners


def project(matrices, tolerance=TOLERANCE):
    keys=list(matrices)
    reverse=reverse_map(keys,tolerance)
    original=np.array(list(matrices.values()))
    result=original.copy()
    for p,q in enumerate(reverse):
        if p>q: continue
        x=(original[p]+original[q].T)*0.5
        result[p]=x
        result[q]=x.T
    return dict(zip(keys,result)),reverse


def encode(original_bytes, matrices):
    """Keep every record marker and offset; replace float64 payloads only."""
    result=bytearray(original_bytes)
    rows=np.array(list(matrices.values()),dtype='<f8').reshape(-1,9)
    if len(result)!=len(rows)*80: raise ValueError('Unexpected binary byte count')
    for q,row in enumerate(rows): result[q*80+4:q*80+76]=row.tobytes()
    return bytes(result)


def norms(a):
    return {'max':float(np.max(np.abs(a))), 'fro':float(np.linalg.norm(a))}


def restore(job, output, binary=None, tolerance=TOLERANCE):
    job=Path(job);output=Path(output);source=Path(binary) if binary else job/'sbar.'
    if output.exists() or output.resolve()==source.resolve():
        raise FileExistsError('Refusing to overwrite any existing file')
    original_bytes=source.read_bytes()
    original=read_sbar(job,source)
    if encode(original_bytes,original)!=original_bytes:
        raise ValueError('Parser roundtrip changed bytes')
    projected,reverse=project(original,tolerance)
    data=encode(original_bytes,projected)
    twice,_=project(projected,tolerance)
    if encode(data,twice)!=data: raise ValueError('Projection not byte-idempotent')
    keys=list(original);a=np.array(list(original.values()));b=np.array(list(projected.values()))
    before=a-np.swapaxes(a[reverse],1,2);after=b-np.swapaxes(b[reverse],1,2)
    report={'input_sha256':hashlib.sha256(original_bytes).hexdigest(),
            'output_sha256':hashlib.sha256(data).hexdigest(),
            'coordinate_tolerance':tolerance,'blocks':len(keys),
            'reverse_coordinate_max':max(float(np.max(np.abs(np.array(keys[p][2:])+keys[q][2:])))
                                         for p,q in enumerate(reverse)),
            'unique_reverse_pairs':sum(p<=q for p,q in enumerate(reverse)),
            'self_reverse_blocks':sum(p==q for p,q in enumerate(reverse)),
            'parser_roundtrip_identical':True,'byte_idempotent':True,
            'before':norms(before),'after':norms(after),'delta_s':norms(b-a),
            'per_block':[{'key':k,'reverse_index':q,'before':norms(before[p]),
                          'after':norms(after[p]),'delta_s':norms(b[p]-a[p])}
                         for p,(k,q) in enumerate(zip(keys,reverse))]}
    # Exclusive create: all validation precedes the first write.
    with output.open('xb') as f: f.write(data)
    reread=read_sbar(job,output)
    if list(reread)!=keys or encode(data,reread)!=data:
        raise ValueError('Written binary failed roundtrip/order validation')
    if source.read_bytes()!=original_bytes: raise ValueError('Original changed during projection')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('job',type=Path)
    p.add_argument('--binary',type=Path,help='Alternate input binary, e.g. idempotence check')
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--report',required=True,type=Path)
    p.add_argument('--tolerance',default=TOLERANCE,type=float)
    args=p.parse_args()
    if args.report.exists() or args.report.resolve() in (args.output.resolve(),(args.binary or args.job/'sbar.').resolve()):
        p.error('Report must be a new file distinct from input and output')
    report=restore(args.job,args.output,args.binary,args.tolerance)
    with args.report.open('x') as f: json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:v for k,v in report.items() if k!='per_block'},indent=2))


if __name__=='__main__': main()
