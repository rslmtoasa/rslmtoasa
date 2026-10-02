#!/usr/bin/env python3
"""Assert that an API request fails with its method-specific boundary message."""
import pathlib
import subprocess
import sys
import tempfile

binary, case, message = sys.argv[1:]
with tempfile.TemporaryDirectory(prefix='lr-config-reject-') as directory:
    result = subprocess.run([str(pathlib.Path(binary).resolve()), case], cwd=directory,
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=30)
    if result.returncode == 0 or message not in result.stdout:
        raise AssertionError(f'{case}: expected rejection {message!r}; exit={result.returncode}\n{result.stdout}')
print(f'PASS: {case}: {message}')
