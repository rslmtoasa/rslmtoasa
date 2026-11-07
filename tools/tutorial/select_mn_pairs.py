#!/usr/bin/env python3
"""Write a small Mn-only ijpair snippet from a Bravais clust file (morton_sfc=false)."""
import argparse
import re
from pathlib import Path
import numpy as np


def read_cluster(path):
    lines = path.read_text().splitlines()
    match = re.fullmatch(r'\s*II\s*=\s*(\d+)\s*', lines[0])
    if not match:
        raise ValueError('Expected the Bravais clust format: II = count')
    count = int(match[1])
    records = []
    for line in lines[1:]:
        if not line.strip():
            continue
        # Source writer: 3F14.8,2I4,3F14.8,2I4; two sites per line.
        for offset in (0, 50):
            block = line[offset:offset + 50]
            if not block.strip():
                continue
            if len(block) != 50:
                raise ValueError('Incomplete fixed-width clust site record')
            xyz = [float(block[i:i + 14]) for i in (0, 14, 28)]
            kind = int(block[42:46])
            int(block[46:50])  # original-site id, distinct from generated cluster index
            records.append((*xyz, kind))
    if len(records) != count:
        raise ValueError(f'Header says {count} sites; parsed {len(records)}')
    array = np.asarray(records)
    if not np.isfinite(array).all():
        raise ValueError('Cluster contains nonfinite coordinates')
    return array[:, :3], array[:, 3].astype(int)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cluster', type=Path)
    parser.add_argument('--mn-types', type=int, nargs='+', required=True,
                        help='Verified symbolic-type indices for Mn, not a moment cutoff')
    parser.add_argument('--centres', type=int, nargs='+', required=True,
                        help='Verified central Mn cluster indices (1 based)')
    parser.add_argument('--neighbours-per-centre', type=int, default=2)
    parser.add_argument('--output', type=Path, default=Path('mn-pairs.nml.inc'))
    args = parser.parse_args()
    if args.neighbours_per_centre < 1:
        parser.error('--neighbours-per-centre must be positive')
    if args.output.exists():
        parser.error('output already exists; choose a new path to preserve it')
    xyz, kinds = read_cluster(args.cluster)
    pairs = []
    for centre in dict.fromkeys(args.centres):
        if not 1 <= centre <= len(kinds) or kinds[centre - 1] not in args.mn_types:
            parser.error(f'centre {centre} is not a verified Mn site')
        distance = np.linalg.norm(xyz - xyz[centre - 1], axis=1)
        candidates = np.flatnonzero(np.isin(kinds, args.mn_types) & (distance > 1e-8))
        order = candidates[np.lexsort((candidates, distance[candidates]))]
        if len(order) < args.neighbours_per_centre:
            parser.error('not enough noncoincident Mn neighbours in the cluster')
        for neighbour in order[:args.neighbours_per_centre]:
            pairs.append((centre, int(neighbour) + 1))
            print(f'{centre} -> {neighbour + 1}: types {kinds[centre - 1]},{kinds[neighbour]}; '
                  f'r/a={distance[neighbour]:.8f}')
    text = ['! Merge into the existing &lattice; retain Te in the electronic structure.',
            f'njij = {len(pairs)}']
    text += [f'ijpair({n}, :) = {i}, {j}' for n, (i, j) in enumerate(pairs, 1)]
    args.output.write_text('\n'.join(text) + '\n')
    print(f'Wrote {args.output}; demonstration neighbours, not a complete exchange-shell list.')


if __name__ == '__main__':
    main()
