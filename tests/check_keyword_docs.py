#!/usr/bin/env python3
"""Fail when a key of a Fortran namelist group in source/ has no entry in docs/source/keywords/*.rst.

A key counts as documented when a keyword page shows it as ``key``, as a heading line, or as `key =` at the
start of a line (example blocks). Usage: check_keyword_docs.py [REPO_ROOT]
"""
import pathlib
import re
import sys

root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parents[1])
# Include files that no source file includes are never compiled; their namelists cannot be read.
sources = sorted(root.glob('source/**/*.f90'))
included = set()
for path in sources:
    included |= set(re.findall(r"^\s*include\s*['\"]([^'\"]+)['\"]", path.read_text(errors='replace'), re.I | re.M))
keys = {}
for path in sources:
    rel = path.relative_to(root / 'source').as_posix()
    if rel.startswith('include_codes/namelists/') and rel not in included:
        continue
    text = re.sub(r'!.*', '', path.read_text(errors='replace'))
    text = re.sub(r'&\s*\n\s*&?', ' ', text)
    for group, body in re.findall(r'^\s*namelist\s*/\s*(\w+)\s*/\s*([^\n]*)', text, re.I | re.M):
        for key in body.split(','):
            if key.strip():
                keys.setdefault(key.strip().lower(), set()).add(group.lower())

docs = '\n'.join(p.read_text() for p in sorted(root.glob('docs/source/keywords/*.rst'))).lower()
documented = set(re.findall(r'``(\w+)', docs))
documented |= set(re.findall(r'(?m)^(\w+)\s*$', docs))
documented |= set(re.findall(r'(?m)^\s*(\w+)\s*(?:\([^)]*\))?\s*=', docs))

missing = sorted((sorted(groups)[0], key) for key, groups in keys.items() if key not in documented)
for group, key in missing:
    print(f'undocumented: &{group} {key}')
print(f'{len(keys)} namelist keys, {len(missing)} without a keyword entry')
sys.exit(1 if missing else 0)
