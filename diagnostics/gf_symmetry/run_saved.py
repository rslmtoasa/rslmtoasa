#!/usr/bin/env python3
"""Run saved exchange jobs in fresh copies; never rerun SCF or overwrite fixtures."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixtures',required=True,type=Path)
    p.add_argument('--binary',required=True,type=Path)
    p.add_argument('--destination',required=True,type=Path)
    p.add_argument('--cases',nargs='+',default=['100','1-00','110','1-1-0'])
    p.add_argument('--one-pair',type=int,help='Smoke run of this original pair only; no SCF/grid changes')
    p.add_argument('--pair-count',type=int,help='Smoke run of first N original pairs only')
    p.add_argument('--audit',default='',help='Global pair index in copied job, or all; empty disables')
    p.add_argument('--h-action',action='store_true')
    p.add_argument('--sbar',type=Path,help='Use this separately projected binary in fresh working copies only')
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--mpi',type=int,default=1)
    p.add_argument('--mpiexec',default='mpiexec')
    args=p.parse_args()
    if args.one_pair and args.pair_count: p.error('Choose one-pair or pair-count')
    binary=args.binary.resolve(strict=True)
    replacement=args.sbar.resolve(strict=True) if args.sbar else None
    dest=args.destination.resolve()
    # No repeated runs into partially populated output directories.
    dest.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy()
    for key in ['RSLMTO_GF_AUDIT','RSLMTO_GF_AUDIT_H','RSLMTO_GF_AUDIT_STRIDE',
                'RSLMTO_AUDIT_SBAR','RSLMTO_AUDIT_SBAR_VIEW']:
        env.pop(key,None)
    env.update(OMP_NUM_THREADS=str(args.threads),OPENBLAS_NUM_THREADS='1')
    if args.audit: env['RSLMTO_GF_AUDIT']=args.audit
    if args.h_action: env['RSLMTO_GF_AUDIT_H']='1'
    manifest={'binary':str(binary),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
              'audit':args.audit,'one_pair':args.one_pair,'pair_count':args.pair_count,
              'threads':args.threads,'mpi_ranks':args.mpi,'mpiexec':args.mpiexec,
              'h_action':args.h_action,'cases':{}}
    for case in args.cases:
        if case not in ('100','1-00','110','1-1-0'): p.error('Unknown fixture case')
        src=args.fixtures.resolve()/case/'exchange_jobs/Mn1'
        target=dest/case
        shutil.copytree(src,target)
        source_hashes={str(f.relative_to(src)):hashlib.sha256(f.read_bytes()).hexdigest()
                       for f in src.rglob('*') if f.is_file()}
        if replacement:
            if not args.audit: p.error('--sbar requires --audit to label the intervention')
            from sbar_audit import read_sbar
            from reciprocal_sbar import project,encode
            raw=replacement.read_bytes()
            projected,_=project(read_sbar(target,replacement))
            if encode(raw,projected)!=raw: raise ValueError('Replacement Sbar is not reciprocal')
            shutil.copyfile(replacement,target/'sbar.')
            env['RSLMTO_AUDIT_SBAR']=str(target/'sbar.')
            env['RSLMTO_AUDIT_SBAR_VIEW']=str(target/'view.sbar.')
        copied_hashes={str(f.relative_to(target)):hashlib.sha256(f.read_bytes()).hexdigest()
                       for f in target.rglob('*') if f.is_file()}
        changed=[name for name,h in source_hashes.items() if copied_hashes.get(name)!=h]
        if changed!=(['sbar.'] if replacement and source_hashes['sbar.']!=copied_hashes['sbar.'] else []):
            raise ValueError(f'Unexpected copied input changes: {changed}')
        input_path=target/'input.nml'
        text=input_path.read_text()
        if not re.search(r"post_processing\s*=\s*['\"]exchange['\"]",text,re.I):
            raise ValueError('Not a saved exchange job')
        if re.search(r'^\s*pre_processing\s*=',text,re.M|re.I):
            raise ValueError('Refusing a pre-processing/SCF job')
        if args.one_pair:
            matches=re.findall(r'^\s*ijpair\((\d+),\s*:\)\s*=\s*(\d+),\s*(\d+).*$',text,re.M|re.I)
            hit=[m for m in matches if int(m[0])==args.one_pair]
            if len(hit)!=1: raise ValueError('Original pair not uniquely found')
            text=re.sub(r'^\s*ijpair\(\d+,\s*:\).*$','',text,flags=re.M|re.I)
            text=re.sub(r'njij\s*=\s*\d+',f'njij = 1\nijpair(1,:) = {hit[0][1]}, {hit[0][2]}',text,flags=re.I)
            input_path.write_text(text)
        if args.pair_count:
            if args.pair_count<1 or args.pair_count>24: p.error('Pair count must be 1..24')
            text=re.sub(r'^\s*ijpair\((\d+),\s*:\).*$',
                        lambda m:m.group(0) if int(m[1])<=args.pair_count else '',text,flags=re.M|re.I)
            text=re.sub(r'njij\s*=\s*\d+',f'njij = {args.pair_count}',text,flags=re.I)
            input_path.write_text(text)
        for old in target.glob('gf-audit-*.dat'): old.unlink()
        start=time.time()
        with (target/'audit-run.log').open('w') as log:
            command=([args.mpiexec,'-n',str(args.mpi)] if args.mpi>1 else [])+[str(binary),'input.nml']
            result=subprocess.run(command,cwd=target,env=env,stdout=log,stderr=subprocess.STDOUT)
        override_loaded='DIAGNOSTIC frozen Sbar override:' in (target/'audit-run.log').read_text()
        manifest['cases'][case]={'returncode':result.returncode,'seconds':time.time()-start,
                                'input_sha256':hashlib.sha256(input_path.read_bytes()).hexdigest(),
                                'source_hashes':source_hashes,'copied_hashes_before_pair_selection':copied_hashes,
                                'changed_on_copy':changed,
                                'sbar_override_loaded':override_loaded,
                                'potentials_unchanged':all(hashlib.sha256((target/name).read_bytes()).hexdigest()==h
                                                          for name,h in source_hashes.items() if name.endswith('_out.nml'))}
        (dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        print(case,{k:v for k,v in manifest['cases'][case].items()
                    if k not in ('source_hashes','copied_hashes_before_pair_selection')},flush=True)
        if result.returncode: raise SystemExit(result.returncode)
        if replacement and not override_loaded:
            raise RuntimeError('Executable did not load saved Sbar; rebuild with the diagnostic hook')


if __name__=='__main__': main()
