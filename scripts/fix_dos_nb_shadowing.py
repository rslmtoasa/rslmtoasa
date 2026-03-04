#!/usr/bin/env python3
"""
Fix local 'nb' variable shadowing in density_of_states.f90.
Renames local 'nb' variable to 'inb' in the specific subroutine where it shadows
the module-level 'nb' from basis_mod (the block size).
"""

import re

def fix_dos_nb_shadowing(filepath):
    """Apply targeted replacements to fix nb shadowing in density_of_states."""

    with open(filepath, 'r') as f:
        content = f.read()

    original_content = content

    # Find the subroutine that contains the problematic local nb variable
    # It's the subroutine around line 255-350 that has "integer :: icode, iii, k, l, ll1, nb, nbp1, ..."

    # Replace in integer declarations: rename 'nb' to 'inb' in lines with "integer :: ... nb ..."
    # Pattern: integer :: ... nb ... (but not in comments or function signatures)

    # 1. Replace 'integer :: ... nb ...' → 'integer :: ... inb ...'
    # Look for the specific pattern with icode and iii (unique to this subroutine)
    content = re.sub(
        r'(integer\s*::\s*icode\s*,\s*iii\s*,\s*k\s*,\s*l\s*,\s*ll1\s*,\s*)nb(\s*,)',
        r'\1inb\2',
        content
    )

    # 2. Replace 'nb = 1' with 'inb = 1' (local assignment)
    content = re.sub(
        r'^\s+nb\s*=\s*1\s*$',
        lambda m: m.group(0).replace('nb = 1', 'inb = 1'),
        content,
        flags=re.MULTILINE
    )

    # 3. Replace 'do k = 1, nb' (where nb was previously set to 1) with 'do k = 1, inb'
    # This is after the 'inb = 1' assignment
    content = re.sub(
        r'(inb\s*=\s*1.*?)(do\s+k\s*=\s*1\s*,\s*)nb\b',
        r'\1\2inb',
        content,
        flags=re.DOTALL,
        count=1  # Only replace the first occurrence after inb = 1
    )

    # 4. Replace 'if (nb > 0)' → 'if (inb > 0)'  and 'nb2(nl) = nb' → 'nb2(nl) = inb'
    content = re.sub(
        r'(inb\s*=\s*1.*?)(if\s*\(\s*)nb(\s*>\s*0\s*\))',
        r'\1\2inb\3',
        content,
        flags=re.DOTALL,
        count=1
    )

    content = re.sub(
        r'(inb\s*=\s*1.*?)(nb2\s*\([^)]*\)\s*=\s*)nb(\b)',
        r'\1\2inb\3',
        content,
        flags=re.DOTALL,
        count=1
    )

    # 5. More robust approach: replace all 'nb' uses that are part of the local variable
    # After "inb = 1" and before "end do" that closes the "do nl" loop
    # Find the section and do targeted replacements

    # Better approach: do multiple passes of specific patterns

    # Replace 'nb2(nl) = nb' with 'nb2(nl) = inb'
    content = re.sub(
        r'nb2\s*\(\s*nl\s*\)\s*=\s*nb\b',
        r'nb2(nl) = inb',
        content
    )

    # Replace 'nb = nb2(nl)' with 'inb = nb2(nl)'
    content = re.sub(
        r'\bnb\s*=\s*nb2\s*\(\s*nl\s*\)',
        r'inb = nb2(nl)',
        content
    )

    # Replace 'do k = 1, nb' (in the context after inb = 1)
    # We'll use a more direct approach
    lines = content.split('\n')
    new_lines = []
    in_target_section = False
    inb_assigned = False

    for i, line in enumerate(lines):
        # Check if we're in the target subroutine
        if 'integer :: icode, iii' in line and 'inb' in line:
            in_target_section = True
        elif in_target_section and line.strip().startswith('end subroutine'):
            in_target_section = False

        # Within the target section, replace patterns
        if in_target_section:
            # After 'inb = 1', replace certain 'nb' with 'inb'
            if 'inb = 1' in line:
                inb_assigned = True

            if inb_assigned and 'do k = 1, nb' in line:
                line = line.replace('do k = 1, nb', 'do k = 1, inb')

            if inb_assigned and 'if (nb > 0)' in line:
                line = line.replace('if (nb > 0)', 'if (inb > 0)')

        new_lines.append(line)

    content = '\n'.join(new_lines)

    if content != original_content:
        with open(filepath, 'w') as f:
            f.write(content)
        print(f"Successfully fixed {filepath}")
        return True
    else:
        print(f"No changes made to {filepath}")
        return False

if __name__ == '__main__':
    filepath = '/Users/andersb/Jobb/rslmto_devel/.claude/worktrees/intelligent-joliot/source/density_of_states.f90'
    fix_dos_nb_shadowing(filepath)
