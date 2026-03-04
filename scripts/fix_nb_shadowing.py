#!/usr/bin/env python3
"""
Fix local 'nb' variable shadowing in recursion.f90.
Renames local neighbour-loop counter 'nb' to 'ineigh' to avoid shadowing
the module-level 'nb' from basis_mod (the block size).
"""

import re
import sys

def fix_nb_shadowing(filepath):
    """Apply targeted replacements to fix nb shadowing."""

    with open(filepath, 'r') as f:
        content = f.read()

    original_content = content

    # Strategy: Use regex replacements with proper boundary detection
    # These patterns should match the local 'nb' (loop variable/neighbor index)
    # but not the module-level 'nb' (block size in zgemm dimensions, etc.)

    # 1. Replace 'do nb = ' with 'do ineigh = ' (loop variable)
    content = re.sub(r'\bdo\s+nb\s*=', r'do ineigh =', content)

    # 2. Replace 'nn(k, nb)' or 'nn(i, nb)' etc with 'nn(?, ineigh)' (neighbor lookup)
    content = re.sub(r'nn\s*\(\s*([a-z_]\w*)\s*,\s*nb\s*\)', r'nn(\1, ineigh)', content)

    # 3. Replace array 3rd subscript patterns: (:,:,nb, or (:,:,nb) or (:,:,nb ]
    content = re.sub(r'(\(\s*:\s*,\s*:\s*,\s*)nb(\s*[,\)])', r'\1ineigh\2', content)

    # 3b. Replace explicit subscript patterns like (1, 1, nb, or similar
    # This catches cases like ee(1, 1, nb, ih) -> ee(1, 1, ineigh, ih)
    content = re.sub(r'(\([\d\s:,]*,\s*)nb(\s*,)', r'\1ineigh\2', content)

    # 4. Replace in OMP private lists: private(...,nb,...) or private(...,nb)
    # This regex handles: private(k, ih, nr, nb, nnmap) or private(...,nb,...)
    content = re.sub(r'\bprivate\s*\(\s*([^)]*?)\b,\s*nb\s*([,\)])', r'private(\1, ineigh\2', content)
    content = re.sub(r'\bprivate\s*\(\s*nb\s*,', r'private(ineigh,', content)

    # 5. Replace 'integer :: nb,' with 'integer :: ineigh,' in declaration lines
    # This is tricky - we only want to replace nb in the context of integer declarations
    # Look for 'integer :: ... nb ...' patterns
    def replace_in_integer_decl(match):
        """Replace nb with ineigh only in integer declarations."""
        line = match.group(0)
        # Replace nb in the variable list (but not in comments or strings)
        # Handle patterns like 'nb,' 'nb )' ' nb,' ' nb '
        line = re.sub(r'\bnb\b', 'ineigh', line)
        return line

    # Find all lines with 'integer ::' declarations and apply replacement
    content = re.sub(
        r'^(\s*integer\s*::.*?\bnb\b.*?)$',
        replace_in_integer_decl,
        content,
        flags=re.MULTILINE
    )

    # 6. Replace 'end do ! End of the loop in the neighbouring' comment patterns
    # (no change needed, but for reference)

    # 7. Handle remaining patterns: this%lattice%nn(... nb ...) or similar
    # Already handled by the nn(k, nb) replacement

    # 8. Replace patterns like 'hallo(:,:,nb,' or 'hall(:,:,nb,'
    # This is already handled by the (:,:,nb pattern

    if content != original_content:
        with open(filepath, 'w') as f:
            f.write(content)
        print(f"Successfully fixed {filepath}")
        return True
    else:
        print(f"No changes made to {filepath}")
        return False

if __name__ == '__main__':
    filepath = '/Users/andersb/Jobb/rslmto_devel/.claude/worktrees/intelligent-joliot/source/recursion.f90'
    fix_nb_shadowing(filepath)
