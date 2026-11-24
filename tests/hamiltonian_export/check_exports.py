#!/usr/bin/env python3
"""Compare effective exports with production Fortran H(k) and H/V action captures."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh


def metadata(path):
    lines = path.read_text().splitlines()
    if lines[0] != '# rs-tb-meta-v2':
        raise ValueError('Unexpected metadata version')
    fields = dict(line.split(maxsplit=1) for line in lines if len(line.split()) == 2)
    n = int(fields['n_atoms'])
    start = lines.index('lattice_vectors_cart') + 1
    cell = np.array([list(map(float, line.split())) for line in lines[start:start+3]]).T
    start = lines.index('positions_cart') + 1
    pos = np.array([list(map(float, line.split()[1:4])) for line in lines[start:start+n]])
    return fields, cell, pos


def records(path, fields):
    count = 0
    result = {}
    nb = int(fields['block_size'])
    ry2ev = float(fields['ry_to_eV'])
    with path.open() as stream:
        if next(stream).strip() != '# rs-tb-hoppings-v2':
            raise ValueError('Unexpected hopping record version')
        for line in stream:
            if line.startswith('# end records '):
                if int(line.split()[-1]) != count or count != int(fields['end_records']):
                    raise ValueError('Wrong record count')
                break
            if line.startswith('#'):
                continue
            items = line.split()
            if len(items) != 12:
                raise ValueError('Wrong record width')
            r1, r2, r3, i, j, _, p, q, _, _ = map(int, items[:10])
            re, im = map(float, items[10:])
            key = (i-1,j-1,r1,r2,r3)
            block = result.setdefault(key, np.zeros((nb,nb), dtype=np.complex128))
            if block[p-1,q-1] != 0:
                raise ValueError('Duplicate record')
            block[p-1,q-1] = complex(re, im) / ry2ev
            count += 1
        else:
            raise ValueError('Missing hopping footer')
    return result


def native_matrices(path):
    n = 0
    found = False
    with path.open() as stream:
        for line in stream:
            items = line.split()
            if items[0] == 'matrix_size':
                n = int(items[1])
            elif items[0] == 'k':
                if found:
                    yield k, hk, sk
                k = np.array(items[2:5], dtype=float)
                hk = np.zeros((n,n), dtype=np.complex128)
                sk = None
                found = True
            elif items[0] == 'H':
                i, j = map(int, items[1:3])
                hk[i-1,j-1] = complex(*map(float, items[3:5]))
            elif items[0] == 'S':
                if sk is None:
                    sk = np.zeros_like(hk)
                i, j = map(int, items[1:3])
                sk[i-1,j-1] = complex(*map(float, items[3:5]))
            elif items[0] == 'end_kpoints':
                if found:
                    yield k, hk, sk
                return
    raise ValueError('Missing k-point footer')


def difference(actual, ref):
    delta = np.abs(actual-ref)
    norm = np.linalg.norm(ref)
    return {'max_absolute': float(delta.max()), 'rms': float(np.sqrt(np.mean(delta**2))),
            'relative_norm': float(np.linalg.norm(delta)/norm) if norm > 1e-13 else None}


def check_fourier(periodic, kspace):
    fields, cell, pos = metadata(periodic/'rs_tb.meta')
    nb = int(fields['block_size'])
    blocks = records(periodic/'rs_tb.tb', fields)
    herm = max(np.max(np.abs(b-blocks.get((j,i,-x,-y,-z),0).conjugate().T))
               if (j,i,-x,-y,-z) in blocks else np.max(np.abs(b))
               for (i,j,x,y,z), b in blocks.items())
    errors = []
    band_errors = []
    native_herm = []
    spectra, native_spectra = [], []
    for k, native, _ in native_matrices(kspace/'rs_tb.kspace'):
        k = k - np.floor(k+0.5)
        result = np.zeros_like(native)
        for (i,j,*r), block in blocks.items():
            d = np.linalg.solve(cell, pos[j]-pos[i]) + r
            result[i*nb:(i+1)*nb,j*nb:(j+1)*nb] += block*np.exp(2j*np.pi*(k@d))
        errors.append(difference(result, native))
        native_herm.append(float(np.max(np.abs(native-native.conjugate().T))))
        spectra.append(np.linalg.eigvalsh(result))
        native_spectra.append(np.linalg.eigvalsh(native))
        band_errors.append(float(np.max(np.abs(spectra[-1]-native_spectra[-1]))))
    summary = {name: max(e[name] for e in errors if e[name] is not None) for name in errors[0]}
    summary.update(matrix_units='Ry', kpoints=len(errors), effective_hermiticity_Ry=float(herm),
                   native_hermiticity_Ry=max(native_herm),
                   eigenvalue_max_eV=max(band_errors)*float(fields['ry_to_eV']))
    emin = min(np.min(spectra), np.min(native_spectra))
    emax = max(np.max(spectra), np.max(native_spectra))
    energy = np.linspace(emin-0.05,emax+0.05,1000)
    width = 0.01

    def gaussian_dos(eigenvalues):
        d = (energy[:,None,None]-np.asarray(eigenvalues))/width
        return np.exp(-0.5*d*d).sum(axis=(1,2))/(width*np.sqrt(2*np.pi)*len(eigenvalues))

    summary['gaussian_dos_states_per_Ry'] = difference(gaussian_dos(spectra),gaussian_dos(native_spectra))
    summary['operator_match'] = ('PASS' if max(summary['max_absolute'], summary['relative_norm']) <= 1e-10 else 'FAIL')
    summary['hermiticity_status'] = ('PASS' if max(summary['effective_hermiticity_Ry'],
                                                  summary['native_hermiticity_Ry']) <= 1e-10 else 'FAIL')
    summary['status'] = ('PASS' if summary['operator_match'] == summary['hermiticity_status'] == 'PASS'
                         and summary['eigenvalue_max_eV'] <= 1e-8 else 'FAIL')
    return summary


def action_capture(path):
    with path.open('rb') as stream:
        nb, n = np.fromfile(stream, np.int32, 2)
        arrays = np.fromfile(stream, np.complex128)
    return arrays.reshape((5,n,3,nb))


def check_action(directory):
    native = action_capture(directory/'native.action')
    imported = action_capture(directory/'import.action')
    if not np.array_equal(native[0], imported[0]):
        raise ValueError('Native and imported probes differ')
    results = {}
    for name, actual, ref in zip(('H','Vx','Vy','Vz'), imported[1:], native[1:]):
        metric = difference(actual, ref)
        relative = metric['relative_norm']
        metric['status'] = ('PASS' if (relative <= 1e-10 if relative is not None
                                       else metric['max_absolute'] <= 1e-12) else 'FAIL')
        results[name] = metric
    return results


def check_cluster_file(directory):
    fields, cell, pos = metadata(directory/'rs_tb.meta')
    if fields['scope'] != 'cluster':
        raise ValueError('Independent finite action requires cluster export')
    probe = action_capture(directory/'native.action')
    psi = probe[0]
    applied = np.zeros_like(probe[1:])
    blocks = records(directory/'rs_tb.tb', fields)
    for (i,j,*r), block in blocks.items():
        hpsi = psi[j] @ block.T
        applied[0,i] += hpsi
        displacement = pos[j]-pos[i] + cell@r
        for axis in range(3):
            applied[axis+1,i] += 1j*displacement[axis]*hpsi
    return {name: difference(actual, ref)
            for name, actual, ref in zip(('H','Vx','Vy','Vz'), applied, probe[1:])}


def check_native_eigenpairs(directory, oracle):
    ref_meta = dict(line.split(' = ',1) for line in (oracle/'fermi_surface.meta').read_text().splitlines()
                    if ' = ' in line)
    n = int(ref_meta['n_basis'])
    nk = int(ref_meta['n_kpoints'])
    vectors_end = nk*n*n
    reader = np.fromfile(oracle/'fermi_surface.bin', dtype=np.float64)
    kref = reader[:3*nk].reshape((nk,3))
    eref = reader[4*nk:(4+n)*nk].reshape((nk,n))
    r = reader[(4+n)*nk:(4+n)*nk+vectors_end].reshape((nk,n,n)).transpose(0,2,1)
    im = reader[(4+n)*nk+vectors_end:(4+n)*nk+2*vectors_end].reshape((nk,n,n)).transpose(0,2,1)
    c = r + 1j*im
    errs = []
    s_errors = []
    e_error = 0.0
    n_checked = 0
    for ik, (k, h, s) in enumerate(native_matrices(directory/'rs_tb.kspace')):
        if not np.allclose(k, kref[ik], atol=1e-14, rtol=0):
            raise ValueError('Different k grids')
        if s is None:
            eigvals = eigh(h,eigvals_only=True)
            y = c[ik]
        else:
            eigvals = eigh(h,s,eigvals_only=True)
            y = s @ c[ik]
            s_from_oracle = np.linalg.inv(c[ik] @ c[ik].conjugate().T)
            s_errors.append(difference(s,s_from_oracle))
        native_h = (y*eref[ik]) @ y.conjugate().T
        errs.append(difference(h,native_h))
        e_error = max(e_error,float(np.max(np.abs(eigvals-eref[ik]))))
        n_checked += 1
    if n_checked != nk:
        raise ValueError('Unexpected k-point count')
    summarize = lambda metrics: {key:max(m[key] for m in metrics) for key in metrics[0]}
    h_diff = summarize(errs)
    s_diff = summarize(s_errors) if s_errors else None
    return {'H_Ry':h_diff,'S':s_diff,'kpoints':nk,'eigenvalue_max_Ry':e_error,
            'status':'PASS' if h_diff['max_absolute']<=1e-10 and (not s_diff or s_diff['max_absolute']<=1e-10)
            and e_error<=1e-10 else 'FAIL'}


def check_transport(directory, oracle, baseline=None):
    results = {}
    files = ('response_total.out','response_orbitals_real.out','response_orbitals_imag.out')
    for name in files:
        ref = np.loadtxt(oracle/name)
        actual = np.loadtxt(directory/name)
        if not np.array_equal(ref[:,0],actual[:,0]):
            raise ValueError('Energy grids differ')
        err = difference(actual[:,1:],ref[:,1:])
        scale = float(np.max(np.abs(ref[:,1:])))
        err['reference_scale'] = scale
        err['status'] = ('PASS' if err['max_absolute'] <= max(1e-12,1e-8*scale)
                         and (err['relative_norm'] is None or err['relative_norm'] <= 1e-7) else 'FAIL')
        if baseline:
            err['baseline_byte_identical'] = (oracle/name).read_bytes() == (baseline/name).read_bytes()
        if name == 'response_total.out':
            err['native_max_imaginary'] = float(np.max(np.abs(ref[:,2])))
        results[name] = err
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--kspace', type=Path)
    parser.add_argument('--action', action='store_true')
    parser.add_argument('--cluster-file', action='store_true')
    parser.add_argument('--native-eigenpairs', type=Path)
    parser.add_argument('--transport', type=Path)
    parser.add_argument('--baseline', type=Path)
    args = parser.parse_args()
    if not any((args.kspace,args.action,args.cluster_file,args.native_eigenpairs,args.transport)):
        parser.error('select a reference or action comparison')
    result = {}
    if args.kspace:
        result['fourier'] = check_fourier(args.directory,args.kspace)
    if args.action:
        result['action'] = check_action(args.directory)
    if args.cluster_file:
        result['cluster_file'] = check_cluster_file(args.directory)
    if args.native_eigenpairs:
        result['native_eigenpairs'] = check_native_eigenpairs(args.directory,args.native_eigenpairs)
    if args.transport:
        result['transport'] = check_transport(args.directory,args.transport,args.baseline)
    print(json.dumps(result,indent=2))
    return int('FAIL' in str(result))


if __name__ == '__main__':
    raise SystemExit(main())
