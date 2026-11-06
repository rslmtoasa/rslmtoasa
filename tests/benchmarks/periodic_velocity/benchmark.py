#!/usr/bin/env python3
"""Build isolated, identically instrumented versions and benchmark frozen fccPt.
Requires git, CMake, gfortran, f90nml and NumPy. No SCF or repository edits.
See README.md for commands and limits. Output directories must be new.
"""
import argparse
import csv
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[3]
MAIN = '3957b2c8f3a4f4ec449dea93f43ffe185b463200'
LEGACY = 'e604f4029e59d35b6e53551dd2d95919a9baf8e0'
TIMER = '''      block
         integer(kind=8) :: started, finished, rate
         integer :: unit
         call system_clock(started, rate)
         call recursion_obj%compute_moments_stochastic()
         call system_clock(finished)
         write(*,'(a,es24.16)') 'BENCH_MOMENTS_SECONDS=', real(finished-started,8)/real(rate,8)
         open(newunit=unit,file='bench_moments.bin',access='stream',form='unformatted',status='replace')
         write(unit) recursion_obj%mu_nm_stochastic
         close(unit)
      end block'''
BONDS_TIMER = '''      block
         integer(kind=8) :: started, finished, rate
         call system_clock(started, rate)
         call this%setup_kubo_bonds()
         call system_clock(finished)
         write(*,'(a,es24.16)') 'BENCH_BONDS_SECONDS=', real(finished-started,8)/real(rate,8)
      end block'''


def main():
    import f90nml
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--fixed-ref', default='HEAD')
    parser.add_argument('--levels', type=int, nargs='+', default=[50])
    parser.add_argument('--versions', nargs='+', choices=['main','fixed','legacy'],
                        default=['main','fixed','legacy'])
    parser.add_argument('--repeats', type=int, default=1)
    args = parser.parse_args()
    dest = args.workdir.resolve()
    dest.mkdir(parents=True, exist_ok=False)
    refs = dict(main=MAIN, fixed=args.fixed_ref, legacy=LEGACY)
    cmake = os.environ.get('CMAKE', 'cmake')
    env = {**os.environ, 'OMP_NUM_THREADS':'4', 'OMP_DYNAMIC':'FALSE',
           'OPENBLAS_NUM_THREADS':'1', 'MKL_NUM_THREADS':'1', 'BLIS_NUM_THREADS':'1',
           'OMP_STACKSIZE':'64M'}
    # The unchanged integration tensor already exceeds small container limits
    # at ll=500. Avoid starting a run whose unavoidable arrays cannot fit.
    limit_file=Path('/sys/fs/cgroup/memory.max')
    if limit_file.exists() and limit_file.read_text().strip()!='max':
        limit=int(limit_file.read_text())
        accepted=[]
        for ll in args.levels:
            minimum=16*(2510+18*18)*ll*ll
            if minimum>limit:
                print(f'SKIP ll={ll}: gamma+moments alone require {minimum} bytes; limit={limit}',flush=True)
            else:
                accepted.append(ll)
        args.levels=accepted
        if not accepted:
            raise SystemExit('No requested level fits the memory limit')
    for version in args.versions:
        src = dest / ('source-'+version)
        subprocess.run(['git','worktree','add','--detach',str(src),refs[version]],cwd=ROOT,check=True)
        path = src/'source/calculation.f90'
        text = path.read_text()
        needle = '      call recursion_obj%compute_moments_stochastic()'
        assert needle in text
        path.write_text(text.replace(needle,TIMER))
        if version == 'fixed':
            path = src/'source/recursion_transport.f90'
            text = path.read_text()
            needle = '      call this%setup_kubo_bonds()\n      norm_a ='
            assert needle in text
            path.write_text(text.replace(needle,BONDS_TIMER+'\n      norm_a =',1))
        build = dest/('build-'+version)
        commands = [[cmake,'-S',str(src),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release',
                     '-DENABLE_OPENMP=ON','-DENABLE_MPI=OFF','-DBUILD_TESTING=OFF',
                     '-DCMAKE_Fortran_COMPILER='+os.environ.get('FC','gfortran')]
                    + shlex.split(os.environ.get('RS_LMTO_TEST_CMAKE_ARGS','')),
                    [cmake,'--build',str(build),'--target','rslmto.x','-j4']]
        with (dest/(version+'-build.log')).open('w') as log:
            for command in commands:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    fields = ['ll','kind','repeat','version','exit_code','wall_seconds','moments_seconds',
              'bonds_seconds','max_rss_kib','user_seconds','system_seconds','run_dir']
    with (dest/'timings.csv').open('w') as f:
        writer = csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
        writer.writeheader()
        for ll in args.levels:
            for repeat in range(1,args.repeats+1):
                for kind in ['spin','orbital']:
                    versions = args.versions if repeat%2 else args.versions[::-1]
                    for version in versions:
                        run = dest/'runs'/f'{ll}-{kind}-{repeat}-{version}'
                        run.mkdir(parents=True)
                        example = ROOT/'example/conductivity/fccPt'/('sigma_xy_z-spin' if kind=='spin' else 'sigma_xy_z-orb')
                        for path in example.glob('*.nml'):
                            shutil.copy2(path,run/path.name)
                        f90nml.patch(run/'input.nml',{'control':{'cond_ll':ll},
                            'lattice':{'pbc':True,'b1':True,'b2':True,'b3':True,
                                       'n1':12,'n2':12,'n3':12}},run/'patched.nml')
                        (run/'patched.nml').replace(run/'input.nml')
                        print('START',run.name,flush=True)
                        with (run/'run.log').open('w') as log:
                            started = time.monotonic()
                            proc = subprocess.Popen([str(dest/('build-'+version)/'bin/rslmto.x')],
                                cwd=run,env=env,stdout=log,stderr=subprocess.STDOUT)
                            _,status,usage = os.wait4(proc.pid,0)
                            proc.returncode = os.waitstatus_to_exitcode(status)
                            wall = time.monotonic()-started
                        log = (run/'run.log').read_text(errors='replace')
                        def timing(label):
                            match = re.search(label+r'=\s*([\d.Ee+\-]+)',log)
                            return float(match[1]) if match else ''
                        row = dict(ll=ll,kind=kind,repeat=repeat,version=version,
                            exit_code=proc.returncode,wall_seconds=wall,
                            moments_seconds=timing('BENCH_MOMENTS_SECONDS'),
                            bonds_seconds=timing('BENCH_BONDS_SECONDS'),max_rss_kib=usage.ru_maxrss,
                            user_seconds=usage.ru_utime,system_seconds=usage.ru_stime,
                            run_dir=str(run.relative_to(dest)))
                        writer.writerow(row); f.flush()
                        print('DONE',row,flush=True)
                        if proc.returncode or row['moments_seconds']=='':
                            raise RuntimeError('Calculation failed: '+str(run))


if __name__ == '__main__':
    main()
