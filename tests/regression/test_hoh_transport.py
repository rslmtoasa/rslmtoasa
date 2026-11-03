"""Compile production transport actions and compare every matrix element to NumPy.
Run with FC=gfortran pytest -q tests/regression/test_hoh_transport.py.
"""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import numpy as np
import pytest
from test_octupole_operators import matrices, PURE

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope='session')
def transport_executable(tmp_path_factory):
    """Build the real recursion module; no copied kernels or mock derived types."""
    fc = shlex.split(os.environ.get('FC', 'gfortran'))
    if not shutil.which(fc[0]):
        pytest.skip('Fortran compiler required (set FC)')
    cmake = os.environ.get('CMAKE', 'cmake')
    if not shutil.which(cmake):
        pytest.skip('CMake required (set CMAKE)')
    build = tmp_path_factory.mktemp('transport_build')
    args = [cmake, '-S', str(ROOT), '-B', str(build),
            '-DCMAKE_Fortran_COMPILER=' + fc[0], '-DCMAKE_BUILD_TYPE=Debug',
            '-DENABLE_OPENMP=ON', '-DENABLE_MPI=OFF', '-DBUILD_TESTING=ON']
    if len(fc) > 1:
        args += ['-DCMAKE_Fortran_FLAGS=' + ' '.join(fc[1:])]
    args += shlex.split(os.environ.get('RS_LMTO_TEST_CMAKE_ARGS', ''))
    log = build / 'build.log'
    with log.open('w') as output:
        for command in [args, [cmake, '--build', str(build), '--target',
                               'transport_operator_test', '-j', '4']]:
            result = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT)
            if result.returncode:
                output.flush()
                pytest.fail('Fortran build failed: ' + str(log) + '\n' + log.read_text()[-8000:])
    return build / 'bin/transport_operator_test'


@pytest.mark.parametrize('zero_overlap', [False, True])
@pytest.mark.parametrize('threads', [1, 2])
def test_compiled_actions(tmp_path, zero_overlap, threads, transport_executable):
    exe = transport_executable
    rng = np.random.default_rng(981)
    b, n, d = 18, 3, 54
    def herm(size):
        a = rng.normal(size=(size,size)) + 1j*rng.normal(size=(size,size))
        return (a+a.conj().T)/20
    dense = herm(d)
    # Chain topology: the HOH correction must generate the absent 1--3 bond.
    dense[:b,2*b:] = 0; dense[2*b:,:b] = 0
    iz = np.array([3,1,2],np.int32)
    nn = np.array([[3,2,0,0],[3,1,3,0],[3,2,0,0]],np.int32)
    h = np.zeros((b,b,4,n),complex)
    for i in range(n):
        for m,j in enumerate([i]+list(nn[i,1:]-1)):
            if j >= 0: h[:,:,m,iz[i]-1] = dense[i*b:(i+1)*b,j*b:(j+1)*b]
    hl = np.zeros_like(h)
    hl[:,:,:,0] = h[:,:,:,iz[0]-1]
    hl[:,:,:,1] = h[:,:,:,iz[1]-1]
    # A different local onsite block tests the impurity override, preserving H=H†.
    hl[:,:,0,0] += herm(b)
    hl[:,:,1,0] *= 1.7
    hl[:,:,1,1] *= 1.7
    o,e,s = [np.stack([herm(b) for _ in range(n)],axis=2) for _ in range(3)]
    if zero_overlap: o[:] = 0
    cr = rng.normal(size=(3,n)); scale = np.array([.4,1.2,2.1])
    direction = np.array([1.,2.,-1.]); direction /= np.linalg.norm(direction)
    psi = np.stack([np.eye(d,dtype=complex)[i*b:(i+1)*b] for i in range(n)],axis=2)
    angular, ops, _ = matrices()
    spin = [np.kron(a,np.eye(9))/2 for a in
            [np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]]
    ls = [np.kron(np.eye(2),a) for a in angular]
    qs = [np.eye(b),*spin,*ls]
    l2 = sum(a@a for a in ls)
    qs += [a@a-l2/3 for a in ls]
    qs += [(ls[i]@ls[j]+ls[j]@ls[i])/2 for i,j in [(0,1),(1,2),(2,0)]]
    qs += [(ls[0]@ls[0]-ls[1]@ls[1])/2,(3*ls[2]@ls[2]-l2)/2]
    qs += [np.kron(np.eye(2),ops[name]) for name in PURE]
    qs = [np.asarray(q, dtype=complex) for q in qs]
    # Exercise the actual slot setup and compiled selectors, not just the kernels.
    cases = [('charge', 'z', np.eye(b))]
    for keyword in ['spin', 'spin_accumulation', 'spin_torque', 'spin_soc_torque', 'soc_spin_torque']:
        cases += [(keyword, axis, q) for axis, q in zip('xyz', spin)]
    for keyword in ['orbital', 'orbital_accumulation', 'orbital_torque']:
        cases += [(keyword, axis, q) for axis, q in zip('xyz', ls)]
    qlabels = ['Qxx', 'Qyy', 'Qzz', 'Qxy', 'Qyz', 'Qzx', 'Qx2y2', 'Q3z2r2']
    for keyword in ['quadrupole', 'oap', 'quadrupole_accumulation', 'oap_accumulation']:
        cases += [(keyword, label, q) for label, q in zip(qlabels, qs[7:15])]
    oct_labels = PURE + ['Mx3y2r2', 'My3z2r2', 'Mz3x2r2', 'Mx3z2r2', 'My3x2r2', 'Mz3y2r2']
    for keyword in ['octupole', 'octupole_current', 'octupole_accumulation']:
        cases += [(keyword, label, np.kron(np.eye(2), ops[label])) for label in oct_labels]
    def write(f,a): np.asarray(a).ravel(order='F').tofile(f)
    with (tmp_path/'fixture.bin').open('wb') as f:
        for a in [nn,iz,cr,scale,h,hl,o,e,s,psi,direction,np.array(len(qs),np.int32),*qs]: write(f,a)
        write(f, np.array(len(cases), np.int32))
        for keyword, label, _ in cases:
            f.write(keyword.ljust(32).encode('ascii'))
            f.write(label.ljust(10).encode('ascii'))
    subprocess.run([str(exe)],cwd=tmp_path,check=True, env={**os.environ, 'OMP_NUM_THREADS': str(threads)})
    raw = np.fromfile(tmp_path/'result.bin',complex)
    offset = 0
    def check(expected):
        nonlocal offset
        actual = raw[offset:offset+b*d*n].reshape((b,d,n),order='F')
        offset += b*d*n
        actual = np.concatenate([actual[:,:,i] for i in range(n)])
        np.testing.assert_allclose(actual,expected,atol=5e-13,rtol=5e-13)
        np.testing.assert_allclose(actual,actual.conj().T,atol=5e-13)
    def local(a):
        out=np.zeros((d,d),complex)
        for i in range(n): out[i*b:(i+1)*b,i*b:(i+1)*b]=a[:,:,iz[i]-1]
        return out
    O,E,S = map(local,(o,e,s))
    x = np.repeat(direction@cr*1.3,b)
    weights=np.repeat(scale[iz-1],b)
    for q in qs:
        Q=np.kron(np.eye(n),q)
        for im in range(2):
            H0=dense.copy()
            if im:
                H0[:b,:b]=hl[:,:,0,0]
                H0[:b,b:2*b]=hl[:,:,1,0]
                H0[b:2*b,:b]=hl[:,:,1,1]
            v=-1j*(x[:,None]-x[None,:])*H0
            for ih in range(2):
                Hns=H0-H0@O@H0+E if ih else H0
                H=Hns+S
                check(H)
                for iw in range(2):
                    v0=v*np.maximum(weights[:,None],weights[None,:]) if iw else v
                    V=v0-v0@O@H0-H0@O@v0 if ih else v0
                    if not iw: np.testing.assert_allclose(V,-1j*(x[:,None]-x[None,:])*H,atol=1e-13)
                    for expected in [V,Q,(Q@V+V@Q)/2,-1j*(Q@Hns-Hns@Q),-1j*(Q@H-H@Q),-1j*(Q@S-S@Q)]: check(expected)
    # Last H0 above is the impurity Hamiltonian; slot b uses the opposite
    # direction and a nonuniform bond weight, while slot a stays unweighted.
    for keyword, label, q in cases:
        Q = np.kron(np.eye(n), q)
        for ih in range(2):
            Hns = H0-H0@O@H0+E if ih else H0
            H = Hns+S
            for iw in range(2):
                v0 = -v*np.maximum(weights[:,None],weights[None,:]) if iw else v
                V = v0-v0@O@H0-H0@O@v0 if ih else v0
                if keyword == 'charge': expected = V
                elif keyword.endswith('_accumulation'): expected = Q
                elif keyword == 'spin_torque': expected = -1j*(Q@Hns-Hns@Q)
                elif keyword == 'orbital_torque': expected = -1j*(Q@H-H@Q)
                elif keyword in ['spin_soc_torque', 'soc_spin_torque']: expected = -1j*(Q@S-S@Q)
                else: expected = (Q@V+V@Q)/2
                check(expected)
    assert offset == raw.size
