#!/usr/bin/env python3
"""Run an existing LR integration deck without modifying its source artifacts."""
import argparse
from pathlib import Path
import re
import shutil
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--case', type=Path, required=True)
parser.add_argument('--input', required=True)
parser.add_argument('--scratch', type=Path, required=True)
args = parser.parse_args()
args.case = args.case.resolve()
args.scratch.mkdir(parents=True, exist_ok=True)
text = (args.case / args.input).read_text()
# Preserve the deck's atom database when moving its working directory.
def database(match):
    return "database = '" + str((args.case / match[1]).resolve()) + "/'"
text = re.sub(r"(?im)^\s*database\s*=\s*['\"]([^'\"]+)['\"]", database, text)
(args.scratch / 'input.nml').write_text(text)
if (args.case / 'Fe.nml').exists():
    shutil.copy2(args.case / 'Fe.nml', args.scratch / 'Fe.nml')
with (args.scratch / 'run.log').open('w') as log:
    result = subprocess.run([str(args.binary.resolve()), 'input.nml'], cwd=args.scratch,
                            stdout=log, stderr=subprocess.STDOUT, timeout=720)
if result.returncode:
    raise RuntimeError((args.scratch / 'run.log').read_text()[-6000:])
print(f'PASS: {args.input}; artifacts: {args.scratch}')
