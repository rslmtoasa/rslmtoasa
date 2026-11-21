#!/usr/bin/env python3
"""Write a manifest with its "extends" cases resolved (tests/manifest_defaults.py) to OUT_JSON.

usage: expand_manifest.py CASES_JSON OUT_JSON

Run by tests/register_tests.cmake at configure time; OUT_JSON is in the build tree.
"""

import json
from pathlib import Path
import sys

from manifest_defaults import resolve_extends

with open(sys.argv[1]) as fh:
    manifest = resolve_extends(json.load(fh))
out = Path(sys.argv[2])
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(manifest, indent=2) + "\n")
