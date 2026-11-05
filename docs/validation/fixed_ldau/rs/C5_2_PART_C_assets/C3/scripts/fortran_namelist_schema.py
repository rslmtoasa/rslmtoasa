#!/usr/bin/env python3
"""Source-grounded namelist member validation for staged decks (private C5.2 tooling)."""
from __future__ import annotations
import re
from pathlib import Path
from typing import Iterable
import f90nml


def source_schema(root: Path) -> dict[str, set[str]]:
    schema: dict[str, set[str]] = {}
    for path in sorted((root / 'source/include_codes/namelists').glob('*.f90')):
        lines = path.read_text().splitlines()
        i = 0
        while i < len(lines):
            line = lines[i].split('!')[0]
            match = re.search(r'\bnamelist\s*/\s*([a-z][a-z0-9_]*)\s*/(.*)', line, re.I)
            if not match:
                i += 1
                continue
            group = match.group(1).lower()
            body = match.group(2)
            while body.rstrip().endswith('&') and i + 1 < len(lines):
                body = body.rstrip()[:-1] + ' ' + lines[i + 1].split('!')[0].strip().lstrip('&')
                i += 1
            members = {part.strip().lower() for part in body.replace('&', '').split(',') if part.strip()}
            normalized = set()
            for member in members:
                token = re.match(r'([a-z][a-z0-9_]*)', member, re.I)
                if not token:
                    raise ValueError(f'cannot parse {group} member {member!r} in {path}')
                normalized.add(token.group(1).lower())
            schema.setdefault(group, set()).update(normalized)
            i += 1
    if not schema:
        raise ValueError(f'no source namelist declarations found under {root}')
    return schema


def validate_input(path: Path, schema: dict[str, set[str]]) -> dict[str, object]:
    nml = f90nml.read(str(path))
    groups = {str(group).lower(): values for group, values in nml.items()}
    unknown_groups = sorted(set(groups) - set(schema))
    unknown_members = {group: sorted(set(map(str.lower, values)) - schema[group])
                       for group, values in groups.items() if group in schema}
    unknown_members = {group: bad for group, bad in unknown_members.items() if bad}
    if unknown_groups or unknown_members:
        raise ValueError(f'{path}: unknown groups={unknown_groups}, unknown members={unknown_members}')
    return {'path': str(path.resolve()), 'groups': sorted(groups),
            'member_counts': {group: len(values) for group, values in sorted(groups.items())},
            'schema_valid': True}
