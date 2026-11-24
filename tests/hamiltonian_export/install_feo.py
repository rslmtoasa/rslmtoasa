#!/usr/bin/env python3
"""Install the archived FeO U4 frozen benchmark into fresh run directories."""
import argparse
import hashlib
from pathlib import Path
import shutil
import tarfile

import f90nml


EXPECTED = 'd6da31bb61056a67583258364dc340411daf79f406912d157dc9e8a2c76a7ea7'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    if hashlib.sha256(args.archive.read_bytes()).hexdigest() != EXPECTED:
        parser.error('Archive SHA-256 differs from the original FeO U4 HOH archive')
    args.destination.mkdir(parents=True, exist_ok=False)
    original = args.destination/'original'
    original.mkdir()
    with tarfile.open(args.archive) as archive:
        files = archive.getmembers()
        if len(files) != 18 or any(not m.isfile() or Path(m.name).name != m.name for m in files):
            parser.error('Unexpected archive contents')
        for member in files:
            with archive.extractfile(member) as src, (original/member.name).open('xb') as dst:
                shutil.copyfileobj(src,dst)

    for hoh in (False,True):
        for scope in ('periodic','kspace'):
            run = args.destination/f'feo-{int(hoh)}-{scope}'
            run.mkdir()
            for path in original.glob('*.nml'):
                if path.name != 'input.nml':
                    shutil.copy2(path,run/path.name)
            nml = f90nml.read(original/'input.nml')
            nml['calculation'].update(pre_processing='none',post_processing='hamiltonian_export')
            nml['self']['cold'] = False
            nml['control'].update(lld=28,cond_ll=28,cheb_backend='legacy',transport_backend='legacy')
            nml['energy']['fix_fermi'] = True
            nml['hamiltonian'].update(hoh=hoh,export='python',export_scope=scope,hubbard_radial_moments='local')
            if scope == 'kspace':
                nml['reciprocal'] = dict(nk1=4,nk2=4,nk3=4,k_offset_x=0.134,
                                         k_offset_y=0.212,k_offset_z=0.163,use_symmetry_reduction=False)
            nml.write(run/'input.nml')
    print(f'Installed original archive and four frozen FeO decks under {args.destination}')


if __name__ == '__main__':
    main()
