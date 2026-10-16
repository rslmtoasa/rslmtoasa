#!/usr/bin/env python3
"""Execute preregistered A4.1 lambda-zero captures into durable directories.

Preparation and execution are deliberately separate. This program is a dry-run
unless --execute is supplied, and refuses to reuse any output directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path('/private/tmp/rslmto-ldau-m2-energy-worktree')
PREP = Path('/private/tmp/a4-m21-ks-lambda0-prep')
OUT = Path('/private/tmp/a4-m21-ks-lambda0-execution')
PREREG = ROOT / 'docs/validation/fixed_ldau/evidence/a4_m21_ks_lambda0_calibration_preregistration.json'
LOCK = ROOT / 'docs/validation/fixed_ldau/evidence/a4_m21_ks_lambda0_preparer_lock.json'
EXECUTABLE = Path('/private/tmp/a4-m2-serial/bin/rslmto.x')
EXPECTED = {
    'prereg_sha256': '50f24c646af2e7c1eedc050d6b0b3c14b7cb55b5d6134e02cfa081913e1ba80f',
    'lock_sha256': 'ab454c7b431a2e4ef31206bd42783613f1507d78e30ae030c347f5f97b0408b9',
    'binary_sha256': '536fa80a34c510782013a5608fc1220dca26b52a10189a2334e40ed43e12a0b4',
    'source_sha256': '2184df95a5560d65d730e6074618ec74210e20fb2eee309c8ac8be7bffb3e010',
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, value: object) -> None:
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def preflight() -> tuple[dict, dict]:
    for path, digest in ((PREREG, EXPECTED['prereg_sha256']),
                         (LOCK, EXPECTED['lock_sha256']),
                         (EXECUTABLE, EXPECTED['binary_sha256']),
                         (ROOT / 'source/self.f90', EXPECTED['source_sha256'])):
        if not path.is_file() or sha(path) != digest:
            raise RuntimeError(f'hash binding failed: {path}')
    prereg = json.loads(PREREG.read_text())
    lock = json.loads(LOCK.read_text())
    if prereg['status'] != 'PREPARED_PRE_RESPONSE_NOT_EXECUTED':
        raise RuntimeError('unexpected lambda-zero preregistration status')
    if lock['execution_status'] != 'No process was launched; no response or calibration outcome exists.':
        raise RuntimeError('preparer lock does not state the frozen no-execution status')
    if len(prereg['cases']) != 8 or len(lock['prepared_cases']) != 8:
        raise RuntimeError('expected exactly eight frozen calibration cases')
    for case in prereg['cases']:
        source_dir = PREP / case['run_id']
        if sha(source_dir / 'input.nml.pre') != case['input_nml_pre_sha256']:
            raise RuntimeError(f"PRE input mismatch: {case['run_id']}")
        if sha(source_dir / 'Fe.nml.pre') != case['Fe_nml_pre_sha256']:
            raise RuntimeError(f"PRE seed mismatch: {case['run_id']}")
        for asset, meta in case['static_assets'].items():
            if sha(source_dir / asset) != meta['sha256']:
                raise RuntimeError(f"static asset mismatch: {case['run_id']}/{asset}")
    return prereg, lock


def prepare_case(case: dict, outdir: Path, prereg_sha: str, lock_sha: str) -> dict:
    src = PREP / case['run_id']
    outdir.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(src / 'input.nml.pre', outdir / 'input.nml')
    shutil.copyfile(src / 'Fe.nml.pre', outdir / 'Fe.nml')
    for asset in case['static_assets']:
        shutil.copyfile(src / asset, outdir / asset)
    env = dict(case['environment'])
    env['RSLMTO_FIXED_U_TRACE'] = str(outdir / 'records.csv')
    planned = {
        'status': 'PREPARED_NOT_STARTED', 'run_id': case['run_id'],
        'created_utc_epoch': time.time(), 'cwd': str(outdir),
        'command': [str(EXECUTABLE)], 'environment': env,
        'bindings': {'lambda0_prereg_sha256': prereg_sha, 'preparer_lock_sha256': lock_sha,
                     'binary_sha256': sha(EXECUTABLE), 'input_nml_pre_sha256': sha(src / 'input.nml.pre'),
                     'Fe_nml_pre_sha256': sha(src / 'Fe.nml.pre')},
        'pre_files': {p.name: sha(p) for p in sorted(outdir.iterdir()) if p.is_file()},
        'stdout_path': str(outdir / 'stdout.log'), 'stderr_path': str(outdir / 'stderr.log'),
    }
    write_json(outdir / 'run_manifest.json', planned)
    return planned


def execute_one(case: dict, outdir: Path, manifest: dict) -> dict:
    env = os.environ.copy()
    env.update(manifest['environment'])
    manifest['status'] = 'RUNNING'
    manifest['started_utc_epoch'] = time.time()
    write_json(outdir / 'run_manifest.json', manifest)
    result = None
    failure = None
    try:
        with (outdir / 'stdout.log').open('wb') as stdout, (outdir / 'stderr.log').open('wb') as stderr:
            result = subprocess.run([str(EXECUTABLE)], cwd=outdir, env=env,
                                    stdin=subprocess.DEVNULL, stdout=stdout,
                                    stderr=stderr, check=False)
        manifest['returncode'] = result.returncode
        manifest['status'] = 'PROCESS_EXITED'
    except BaseException as exc:  # persist a receipt even for interruption
        failure = f'{type(exc).__name__}: {exc}'
        manifest['failure'] = failure
        manifest['status'] = 'LAUNCH_OR_CAPTURE_FAILURE'
    finally:
        manifest['finished_utc_epoch'] = time.time()
        manifest['post_files'] = {
            p.name: {'sha256': sha(p), 'size_bytes': p.stat().st_size}
            for p in sorted(outdir.iterdir()) if p.is_file() and p.name != 'run_manifest.json'
        }
        write_json(outdir / 'run_manifest.json', manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--execute', action='store_true', help='launch the frozen SCF cases; default only prepares durable manifests')
    parser.add_argument('--case', action='append', help='run_id to select; may be repeated; default selects all eight')
    parser.add_argument('--out-root', type=Path, default=OUT)
    args = parser.parse_args()
    def interrupted(signum, _frame):
        raise KeyboardInterrupt(f'received signal {signum}')
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    prereg, lock = preflight()
    prereg_sha, lock_sha = sha(PREREG), sha(LOCK)
    selected = [c for c in prereg['cases'] if args.case is None or c['run_id'] in args.case]
    if not selected or (args.case and {c['run_id'] for c in selected} != set(args.case)):
        raise SystemExit('unknown or empty --case selection')
    if args.out_root.exists():
        raise SystemExit(f'refusing to overwrite persistent evidence directory: {args.out_root}')
    args.out_root.mkdir(parents=True)
    batch = {'schema': 'rslmto-a4-m21-lambda0-execution-batch-v1',
             'status': 'PREPARED_NOT_STARTED', 'created_utc_epoch': time.time(),
             'execute_requested': args.execute, 'preregistration_sha256': prereg_sha,
             'preparer_lock_sha256': lock_sha, 'binary_sha256': sha(EXECUTABLE),
             'case_ids': [c['run_id'] for c in selected], 'cases': {}}
    write_json(args.out_root / 'execution_manifest.json', batch)
    for case in selected:
        d = args.out_root / case['run_id']
        m = prepare_case(case, d, prereg_sha, lock_sha)
        batch['cases'][case['run_id']] = {'run_manifest': str(d / 'run_manifest.json'),
                                          'status': m['status']}
        write_json(args.out_root / 'execution_manifest.json', batch)
    if not args.execute:
        print(json.dumps({'status': batch['status'], 'execution_manifest': str(args.out_root / 'execution_manifest.json'),
                          'message': 'dry-run only; no executable was launched'}, sort_keys=True))
        return 0
    batch['status'] = 'EXECUTION_IN_PROGRESS'
    batch['started_utc_epoch'] = time.time()
    write_json(args.out_root / 'execution_manifest.json', batch)
    for case in selected:
        d = args.out_root / case['run_id']
        m = json.loads((d / 'run_manifest.json').read_text())
        done = execute_one(case, d, m)
        batch['cases'][case['run_id']]['status'] = done['status']
        batch['cases'][case['run_id']]['returncode'] = done.get('returncode')
        write_json(args.out_root / 'execution_manifest.json', batch)
    batch['finished_utc_epoch'] = time.time()
    batch['status'] = 'ALL_PROCESSES_EXITED' if all(v['status'] == 'PROCESS_EXITED' for v in batch['cases'].values()) else 'ONE_OR_MORE_LAUNCH_FAILURES'
    write_json(args.out_root / 'execution_manifest.json', batch)
    print(json.dumps({'status': batch['status'], 'execution_manifest': str(args.out_root / 'execution_manifest.json')}, sort_keys=True))
    return 0 if batch['status'] == 'ALL_PROCESSES_EXITED' else 1


if __name__ == '__main__':
    sys.exit(main())
