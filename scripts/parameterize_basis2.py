#!/usr/bin/env python3
"""
Second-pass basis-dimension parameterisation for RS-LMTO-ASA.

Picks up the remaining literal `18` and `9` occurrences not caught by the
first-pass script, including BLAS leading-dimension arguments, LAPACK workspace
sizes, safe_alloc / allocate argument lists, and array-slice bounds.

Rules applied:
  1.  (/18, 18,  →  (/nb, nb,    (safe_alloc / allocate array-constructor)
  2.  (/18, 18)  →  (/nb, nb)
  3.  (18, kk)  →  (nb, kk)     (local v arrays)
  4.  3*18 - 2  →  3*nb - 2     (LAPACK workspace)
  5.  1:18       →  1:nb         (array slices)
  6.  10:18      →  norb+1:nb   (spin-down block slices)
  7.  1:9        →  1:norb       (spin-up block slices — targeted)
  8.  Broad  \b18\b  →  nb       (catch remaining occurrences outside comments)
  9.  norb+1:nb  already set above; also norb+1:nb alias in spin_off+1:nb ?
     Keep as norb+1:nb for readability.
"""

import re
import sys
from pathlib import Path

FILES = [
    'hamiltonian.f90',
    'recursion.f90',
    'green.f90',
    'bands.f90',
    'density_of_states.f90',
    'mix.f90',
    'conductivity.f90',
    'self.f90',
    'charge.f90',
]


def sub_no_comment(pattern, repl, text, flags=0):
    """Apply regex substitution only on non-comment portions of each line."""
    lines = text.split('\n')
    out = []
    for line in lines:
        # Find comment start (! that is not inside a string literal)
        # Simple heuristic: split at '!' respecting strings
        code_part, _, comment_part = _split_comment(line)
        new_code = re.sub(pattern, repl, code_part, flags=flags)
        if comment_part is not None:
            out.append(new_code + '!' + comment_part)
        else:
            out.append(new_code)
    return '\n'.join(out)


def _split_comment(line):
    """Return (code, '!', rest) or (line, None, None) if no bare comment."""
    in_string = False
    quote_char = None
    for i, ch in enumerate(line):
        if in_string:
            if ch == quote_char:
                in_string = False
        else:
            if ch in ('"', "'"):
                in_string = True
                quote_char = ch
            elif ch == '!':
                return line[:i], '!', line[i+1:]
    return line, None, None


def process_file(src_dir, fname):
    fpath = Path(src_dir) / fname
    if not fpath.exists():
        print(f'  SKIP (not found): {fname}')
        return

    original = fpath.read_text(encoding='utf-8')
    text = original

    # ---- array-constructor (/18, 18, ...) → (/nb, nb, ...)
    text = sub_no_comment(r'\(/\s*18\s*,\s*18\s*,', '(/nb, nb,', text)
    text = sub_no_comment(r'\(/\s*18\s*,\s*18\s*\/', '(/nb, nb/', text)

    # ---- safe_alloc / allocate  (18, kk) or (18, this%...)
    # dimension(18, expr)  already partially done; also catches leading-dim args
    text = sub_no_comment(r'\(18\s*,\s*this%', '(nb, this%', text)
    text = sub_no_comment(r',\s*18\s*,\s*this%', ', nb, this%', text)

    # ---- LAPACK workspace size  3*18 - 2  or  3 * 18 - 2
    text = sub_no_comment(r'3\s*\*\s*18\s*-\s*2', '3*nb - 2', text)

    # ---- array-slice  1:18 → 1:nb
    text = sub_no_comment(r'1:18\b', '1:nb', text)

    # ---- spin-down block slice  10:18 → norb+1:nb
    text = sub_no_comment(r'10:18\b', 'norb+1:nb', text)

    # ---- spin-up block slice  1:9 → 1:norb   (within index context)
    # Only in array slices: preceded by '(' or ',' and followed by ',' or ')'
    text = sub_no_comment(r'(?<=\()1:9(?=[,\)])', '1:norb', text)
    text = sub_no_comment(r'(?<=,\s)1:9(?=[,\)])', '1:norb', text)
    text = sub_no_comment(r'(?<=, )1:9(?=[,\)])', '1:norb', text)
    # Also handle  sum(arr(1:9, ie))  pattern
    text = sub_no_comment(r'\b1:9\b', '1:norb', text)

    # ---- allocate(arr(18, kk))  or  allocate(arr(18,  ))
    text = sub_no_comment(r'\(18\s*,\s*(\w)', r'(nb, \1', text)
    text = sub_no_comment(r',\s*18\s*\)', ', nb)', text)
    text = sub_no_comment(r',\s*18\s*,', ', nb,', text)

    # ---- Broad sweep: any remaining \b18\b outside comments
    text = sub_no_comment(r'\b18\b', 'nb', text)

    if text != original:
        fpath.write_text(text, encoding='utf-8')
        print(f'  UPDATED: {fname}')
    else:
        print(f'  unchanged: {fname}')


def main():
    if len(sys.argv) < 2:
        print('Usage: parameterize_basis2.py <source_dir>')
        sys.exit(1)
    src_dir = sys.argv[1]
    for fname in FILES:
        process_file(src_dir, fname)


if __name__ == '__main__':
    main()
