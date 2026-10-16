#!/usr/bin/env python3
"""Execute only a reviewed, hash-locked A4.1 W1 response registration.

The registration builder prepares immutable case PRE directories. This helper
validates every binding before launch, refuses existing outputs, and runs at
most two independent SCFs concurrently. It does not analyze or classify data.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, obj: object) -> None:
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def resolve(value: str | Path) -> Path:
    return Path(value).expanduser().resolve()


def check_case(case: dict, binary_sha: str) -> tuple[Path, list[Path]]:
    case_id = case['run_id']
    work = resolve(case['pre_dir'])
    if not work.is_dir():
        raise ValueError(f'{case_id}: missing prepared directory')
    if sha(work / 'input.nml.pre') != case['input_nml_pre_sha256']:
        raise ValueError(f'{case_id}: PRE input digest mismatch')
    if sha(work / 'Fe.nml.pre') != case['seed_sha256']:
        raise ValueError(f'{case_id}: PRE seed digest mismatch')
    if sha(work / 'input.nml') != case['input_nml_pre_sha256']:
        raise ValueError(f'{case_id}: launch input differs from frozen PRE input')
    if sha(work / 'Fe.nml') != case['seed_sha256']:
        raise ValueError(f'{case_id}: launch seed differs from frozen PRE seed')
    for name, digest in case['static_assets'].items():
        if sha(work / name) != digest:
            raise ValueError(f'{case_id}: static asset digest mismatch: {name}')

    command = [str(resolve(v)) for v in case['command']]
    if not command or not Path(command[0]).is_file() or sha(Path(command[0])) != binary_sha:
        raise ValueError(f'{case_id}: executable path/hash mismatch')
    outputs = case['outputs']
    required = {'stdout', 'full_scf_log', 'stderr', 'trace', 'final_state'}
    if not required.issubset(outputs):
        raise ValueError(f'{case_id}: output inventory must include {sorted(required)}')
    output_paths = [resolve(v) for v in outputs.values()]
    if resolve(outputs['stdout']) != resolve(outputs['full_scf_log']):
        raise ValueError(f'{case_id}: full_scf_log must alias the captured stdout path')
    unique_paths = {key: resolve(value) for key, value in outputs.items()
                    if key != 'full_scf_log'}
    if len(set(unique_paths.values())) != len(unique_paths):
        raise ValueError(f'{case_id}: output paths collide beyond the declared log alias')
    for path in output_paths:
        if path.exists():
            raise ValueError(f'{case_id}: refusing to overwrite existing output {path}')
        if path.parent != work:
            raise ValueError(f'{case_id}: output must be in its bound PRE directory')
    if case.get('mode') not in ('spectral', 'local') or case.get('operator') != 'w1':
        raise ValueError(f'{case_id}: unexpected mode/operator')
    if case.get('threshold') not in ('standard', 'tight'):
        raise ValueError(f'{case_id}: unexpected threshold label')
    if float(case.get('lambda_ry')) == 0.0:
        raise ValueError(f'{case_id}: response case has zero lambda')
    env = case.get('environment', {})
    if env.get('RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY') != '1':
        raise ValueError(f'{case_id}: trace is not capture-only')
    if env.get('RSLMTO_FIXED_U_TRACE') != str(outputs['trace']):
        raise ValueError(f'{case_id}: trace environment/path mismatch')
    return work, output_paths


def run_case(case: dict, binary_sha: str, reg_sha: str, out_manifest: Path) -> dict:
    work, outpaths = check_case(case, binary_sha)
    output_by_name = {key: resolve(value) for key, value in case['outputs'].items()}
    stdout, stderr = output_by_name['stdout'], output_by_name['stderr']
    executable = resolve(case['command'][0])
    env = os.environ.copy()
    env.update(case['environment'])
    receipt = {
        'schema': 'rslmto-a4-m21-w1-case-run-v1',
        'run_id': case['run_id'], 'status': 'PREPARED',
        'registration_sha256': reg_sha, 'binary_sha256': binary_sha,
        'pre_dir': str(work), 'command': [str(executable)],
        'environment': case['environment'],
        'input_nml_pre_sha256': sha(work / 'input.nml.pre'),
        'seed_sha256': sha(work / 'Fe.nml.pre'),
        'static_assets': {name: sha(work / name) for name in case['static_assets']},
        'outputs': {name: str(path) for name, path in output_by_name.items()},
        'started_utc_epoch': time.time(),
    }
    write_json(out_manifest, receipt)
    receipt['status'] = 'RUNNING'
    write_json(out_manifest, receipt)
    try:
        with stdout.open('xb') as fout, stderr.open('xb') as ferr:
            proc = subprocess.run([str(executable)], cwd=work, env=env,
                                  stdin=subprocess.DEVNULL, stdout=fout, stderr=ferr,
                                  check=False)
        receipt['returncode'] = proc.returncode
        receipt['status'] = 'PROCESS_EXITED'
    except BaseException as exc:
        receipt['status'] = 'LAUNCH_OR_CAPTURE_FAILURE'
        receipt['failure'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        receipt['finished_utc_epoch'] = time.time()
        receipt['post_files'] = {
            key: {'sha256': sha(path), 'size_bytes': path.stat().st_size}
            for key, path in output_by_name.items() if path.is_file()
        }
        missing = set(('stdout', 'full_scf_log', 'stderr', 'trace', 'final_state')) - set(receipt['post_files'])
        if missing:
            receipt['status'] = 'INCOMPLETE_OUTPUTS'
            receipt['missing_outputs'] = sorted(missing)
        write_json(out_manifest, receipt)
    return receipt


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('registration', type=Path)
    ap.add_argument('--registration-sha256', required=True,
                    help='independently reviewed immutable registration digest')
    ap.add_argument('--binary-sha256', required=True,
                    help='expected frozen executable digest from the registration review')
    ap.add_argument('--manifest-dir', type=Path, required=True,
                    help='new persistent directory for execution manifests')
    ap.add_argument('--execute', action='store_true', help='launch; default is preflight only')
    args = ap.parse_args()

    reg_path = resolve(args.registration)
    if sha(reg_path) != args.registration_sha256:
        raise SystemExit('registration digest mismatch')
    reg = json.loads(reg_path.read_text())
    if reg.get('status') != 'PRE_RESPONSE_LOCKED_NO_NONZERO_RUNS':
        raise SystemExit('registration is not locked for response execution')
    if reg.get('schema') != 'rslmto-a4-m21-w1-response-prereg-v1':
        raise SystemExit('unexpected registration schema')
    cases = reg.get('response_cases')
    if not isinstance(cases, list) or len(cases) != 24:
        raise SystemExit('expected exactly 24 signed W1 response cases')
    ids = [case.get('run_id') for case in cases]
    if len(set(ids)) != 24 or any(not x for x in ids):
        raise SystemExit('run identifiers are missing or duplicated')
    if args.manifest_dir.exists():
        raise SystemExit('refusing to overwrite execution manifest directory')
    checked = [(case, *check_case(case, args.binary_sha256)) for case in cases]

    # Reserve a conservative 24 * largest observed calibration footprint,
    # plus 1 GiB free-space margin for filesystem/log growth.
    largest_calibration_bytes = 281_492 * 1024
    required_free = 24 * largest_calibration_bytes + (1 << 30)
    if not args.manifest_dir.parent.is_dir():
        raise SystemExit('manifest directory parent must already exist')
    free = shutil.disk_usage(args.manifest_dir.parent).free
    if free < required_free:
        raise SystemExit(f'insufficient free space for bounded batch: need {required_free} bytes, have {free}')

    args.manifest_dir.mkdir(parents=True)
    batch = {
        'schema': 'rslmto-a4-m21-w1-execution-batch-v1',
        'status': 'PREPARED_NOT_STARTED', 'registration_path': str(reg_path),
        'registration_sha256': args.registration_sha256,
        'binary_sha256': args.binary_sha256, 'case_count': 24,
        'max_concurrency': 2, 'free_bytes_at_preflight': free,
        'minimum_free_bytes': required_free, 'case_manifests': {},
    }
    batch_path = args.manifest_dir / 'execution_manifest.json'
    write_json(batch_path, batch)
    if not args.execute:
        print(json.dumps({'status': batch['status'], 'case_count': 24,
                          'max_concurrency': 2, 'free_bytes': free,
                          'minimum_free_bytes': required_free,
                          'manifest': str(batch_path)}, sort_keys=True))
        return 0

    batch['status'] = 'EXECUTION_IN_PROGRESS'
    batch['started_utc_epoch'] = time.time()
    write_json(batch_path, batch)
    futures = {}
    failures = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        for case, _, _ in checked:
            case_manifest = args.manifest_dir / f"{case['run_id']}.json"
            batch['case_manifests'][case['run_id']] = str(case_manifest)
            futures[pool.submit(run_case, case, args.binary_sha256,
                                args.registration_sha256, case_manifest)] = case['run_id']
        write_json(batch_path, batch)
        for future in as_completed(futures):
            run_id = futures[future]
            try:
                receipt = future.result()
                batch.setdefault('cases', {})[run_id] = {
                    'status': receipt['status'], 'returncode': receipt.get('returncode')}
                if receipt['status'] != 'PROCESS_EXITED' or receipt.get('returncode') != 0:
                    failures.append(run_id)
            except BaseException as exc:
                batch.setdefault('cases', {})[run_id] = {
                    'status': 'LAUNCH_OR_CAPTURE_FAILURE', 'failure': f'{type(exc).__name__}: {exc}'}
                failures.append(run_id)
            write_json(batch_path, batch)
    batch['finished_utc_epoch'] = time.time()
    batch['status'] = 'ALL_PROCESSES_EXITED' if not failures else 'FAILURES_RECORDED'
    batch['failed_run_ids'] = failures
    write_json(batch_path, batch)
    print(json.dumps({'status': batch['status'], 'manifest': str(batch_path),
                      'failed_run_ids': failures}, sort_keys=True))
    return 0 if not failures else 1


if __name__ == '__main__':
    sys.exit(main())
