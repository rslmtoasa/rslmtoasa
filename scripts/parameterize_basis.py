#!/usr/bin/env python3
"""
Mechanical basis-dimension parameterisation for RS-LMTO-ASA.

Replaces hardcoded 18 (nb) and 9 (norb / spin_off) dimension literals with
module variables from basis_mod.  Also adds `use basis_mod` to the USE block
of each file.

Safe replacement rules applied (in order):
  1. dimension(18, 18 ...)  -> dimension(nb, nb ...)
  2. 18*18                  -> nb*nb
  3. ndim   = 18            -> ndim   = nb
  4. do VAR = 1, 18         -> do VAR = 1, nb
  5. Loop copy: do VAR = 1, 18  already covered
  6. +9 spin-offset patterns -> +spin_off
  7. dimension(9, 9 ...)    -> dimension(norb, norb ...)
  8. do VAR = 1, 9          -> do VAR = 1, norb
  9. 9*9                    -> norb*norb
 10. Inject `use basis_mod, only: nb, norb, spin_off` after the last
     non-basis_mod `use` statement before `implicit none`.

Files treated differently:
  recursion.f90: also convert cheb_mom_temp fixed declaration -> allocatable.
"""

import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Per-file configuration
# ---------------------------------------------------------------------------
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
    'exchange.f90',
]

# These files are already parameterised (exchange) or need no change
SKIP_FILES = {'exchange.f90'}

# ---------------------------------------------------------------------------
# Substitution helpers
# ---------------------------------------------------------------------------

def sub(pattern, repl, text, flags=0):
    return re.sub(pattern, repl, text, flags=flags)

def add_basis_use(text):
    """Insert `use basis_mod, only: nb, norb, spin_off` into the USE block."""
    use_line = '   use basis_mod, only: nb, norb, spin_off\n'

    # Skip if already present
    if 'basis_mod' in text:
        return text

    # Find the last 'use ...' line before 'implicit none' in the module body
    # Pattern: look for block of use statements followed by implicit none
    # Insert after the last 'use' line in that block.
    pattern = r'((?:[ \t]*use\s+\S[^\n]*\n)+)([ \t]*implicit\s+none)'
    def inserter(m):
        use_block = m.group(1)
        implicit  = m.group(2)
        return use_block + use_line + implicit
    new_text = re.sub(pattern, inserter, text, count=1, flags=re.IGNORECASE)
    if new_text == text:
        # Fallback: insert just before the first implicit none
        new_text = re.sub(
            r'([ \t]*implicit\s+none)',
            use_line + r'\1',
            text, count=1, flags=re.IGNORECASE)
    return new_text

def replace_nb(text):
    """Replace hardcoded 18 → nb in dimensional/loop contexts."""

    # dimension(18, 18, ...) or dimension(18,18,...)
    text = re.sub(r'\bdimension\s*\(\s*18\s*,\s*18\b', 'dimension(nb, nb', text)

    # allocate(...(18, 18, ...)) - catches the allocation argument
    # Replace ,(18,18,  or (18,18,  or (18,18) in allocate calls
    # We use a lookahead to avoid touching other 18s
    # Safe: replace 18*18
    text = re.sub(r'\b18\s*\*\s*18\b', 'nb*nb', text)

    # ndim = 18  (BLAS leading-dimension scalar)
    text = re.sub(r'\bndim\s*=\s*18\b', 'ndim = nb', text)

    # hblocksize = 18
    text = re.sub(r'\bhblocksize\s*=\s*18\b', 'hblocksize = nb', text)

    # do VAR = 1, 18
    text = re.sub(r'\bdo\s+(\w+)\s*=\s*1\s*,\s*18\b', r'do \1 = 1, nb', text)

    # Standalone dimension size in allocate or local declarations:
    # e.g.  complex(...) :: foo(18, 18)  or  real(...), dimension(18, 18) :: bar
    # Already covered by the first rule above for the 'dimension(18,18' pattern.
    # Also catch allocatable array allocate sizes:
    #   allocate(arr(nb, nb, n)) -- the 18,18 inside allocate(arr(
    # We do a targeted scan: ,(18, or (18, where followed by another 18 or nb or a var
    # Actually the dimension(...) rule above already covers declarations.
    # For allocate() calls the pattern is different; let's replace all remaining
    # standalone (18, 18) patterns:
    text = re.sub(r'\(18\s*,\s*18\)', '(nb, nb)', text)
    text = re.sub(r'\(18\s*,\s*18\s*,', '(nb, nb,', text)
    text = re.sub(r',\s*18\s*,\s*18\)', ', nb, nb)', text)
    text = re.sub(r',\s*18\s*,\s*18\s*,', ', nb, nb,', text)

    return text

def replace_norb(text):
    """Replace hardcoded 9 → norb in orbital-dimension contexts."""

    # dimension(9, 9, ...) or dimension(9,9,...)
    text = re.sub(r'\bdimension\s*\(\s*9\s*,\s*9\b', 'dimension(norb, norb', text)

    # 9*9 → norb*norb
    text = re.sub(r'\b9\s*\*\s*9\b', 'norb*norb', text)

    # do VAR = 1, 9   (orbital loops)
    text = re.sub(r'\bdo\s+(\w+)\s*=\s*1\s*,\s*9\b', r'do \1 = 1, norb', text)

    # Standalone (9, 9) in array index / size context
    text = re.sub(r'\(9\s*,\s*9\)', '(norb, norb)', text)
    text = re.sub(r'\(9\s*,\s*9\s*,', '(norb, norb,', text)

    return text

def replace_spin_off(text):
    """Replace +9 spin-block offset → +spin_off in array indices."""
    # Patterns like: i + 9, j + 9, m + 9, k + 9,  and also +9
    # Only replace in array-index-like context: preceded by identifier or ')'
    # and followed by ',' or ')' or space+operator
    text = re.sub(r'\+\s*9\b', '+spin_off', text)
    return text

# ---------------------------------------------------------------------------
# File-specific fixup: cheb_mom_temp in recursion.f90
# ---------------------------------------------------------------------------

def fix_cheb_mom_temp(text):
    """
    Change:
      complex(rp), dimension(18, 18) :: cheb_mom_temp
    to:
      complex(rp), allocatable :: cheb_mom_temp(:,:)

    Also ensure it is allocated in restore_to_default / wherever needed.
    (The allocation call must be added manually or via a separate pass.)
    """
    text = re.sub(
        r'complex\s*\(\s*rp\s*\)\s*,\s*dimension\s*\(\s*nb\s*,\s*nb\s*\)\s*::\s*cheb_mom_temp',
        'complex(rp), allocatable :: cheb_mom_temp(:,:)',
        text
    )
    return text

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def process_file(src_dir, fname):
    fpath = Path(src_dir) / fname
    if not fpath.exists():
        print(f'  SKIP (not found): {fname}')
        return

    original = fpath.read_text(encoding='utf-8')
    text = original

    text = add_basis_use(text)
    text = replace_nb(text)
    text = replace_norb(text)
    text = replace_spin_off(text)

    if fname == 'recursion.f90':
        text = fix_cheb_mom_temp(text)

    if text != original:
        fpath.write_text(text, encoding='utf-8')
        print(f'  UPDATED: {fname}')
    else:
        print(f'  unchanged: {fname}')


def main():
    if len(sys.argv) < 2:
        print('Usage: parameterize_basis.py <source_dir>')
        sys.exit(1)

    src_dir = sys.argv[1]
    for fname in FILES:
        if fname in SKIP_FILES:
            print(f'  SKIP (excluded): {fname}')
            continue
        process_file(src_dir, fname)

if __name__ == '__main__':
    main()
