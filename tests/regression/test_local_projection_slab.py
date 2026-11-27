"""Two-layer magnetic slab with closed or open lateral boundaries (open z)."""
import csv
import os
import subprocess

import numpy as np
import pytest
from test_hoh_transport import transport_executable

SX = np.array([[0, .5], [.5, 0]])
ACTIVE = [0, 9]


def slab(directory, exe, width, periodic, rep, order, hoh, layout, projection, observable="spin_torque_onsite", pol="x", impurity=False):
    sites = [(x, y, z) for x in range(width) for y in range(width) for z in range(2)]
    count = len(sites)
    types = 2
    refs = [sites.index((rep, rep, z)) for z in range(types)]
    nn = np.zeros((count, 7), np.int32)
    blocks = np.zeros((18, 18, 7, types), complex)
    h = np.zeros((2*count, 2*count), complex)
    v = np.zeros_like(h)
    o = np.zeros((18, 18, types), complex)
    soc = np.zeros_like(o)
    step = [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]
    hopping = np.array([[.14+.05j, .06+.02j], [-.03+.07j, .19-.05j]])
    for i, (x,y,z) in enumerate(sites):
        si = slice(2*i, 2*i+2)
        onsite = np.diag([-.31, .24]) if z else np.eye(2)*.15
        blocks[:,:,0,z][np.ix_(ACTIVE, ACTIVE)] = onsite
        h[si,si] = onsite
        o[:,:,z][np.ix_(ACTIVE, ACTIVE)] = np.diag([.12, .22])
        soc[:,:,z][np.ix_(ACTIVE, ACTIVE)] = np.array([[0,.06j],[-.06j,0]])
        nn[i,0] = 7
        for k,(dx,dy,dz) in enumerate(step):
            coord = (x+dx,y+dy,z+dz)
            if periodic:
                coord = (coord[0]%width, coord[1]%width, coord[2])
            if coord not in sites:
                continue
            j = sites.index(coord)
            nn[i,k+1] = j+1
            hij = hopping*(1+.3*z if dz==0 else 1)*(1+.15*(k//2))
            if k%2:
                hij = hij.conj().T
            blocks[:,:,k+1,z][np.ix_(ACTIVE, ACTIVE)] = hij
            sj = slice(2*j,2*j+2)
            h[si,sj] = hij
            v[si,sj] = 1j*dx*hij
    hall = blocks[:,:,:,:0].copy()
    if impurity:
        refs[0] = 0
        hall = blocks[:,:,:,:1].copy()
        field = np.array([[.14,.09j],[-.09j,-.14]])
        hall[:,:,0,0][np.ix_(ACTIVE,ACTIVE)] += field
        h[:2,:2] += field
    if hoh:
        overlap = np.diag(np.tile([.12,.22], count))
        v = v - v@overlap@h - h@overlap@v
        h = h - h@overlap@h
    a = np.zeros_like(h)
    for i in range(count):
        si = slice(2*i, 2*i+2)
        hii = h[si,si]
        a[si,si] = -1j*(SX@hii-hii@SX)
    h += np.kron(np.eye(count), np.array([[0,.06j],[-.06j,0]]))
    if observable != "spin_torque_onsite":
        # Spread the two active states across p/d orbitals so all tensor ranks
        # can have a finite projected matrix element.
        rng = np.random.default_rng(174)
        columns = rng.normal(size=(18,2)) + 1j*rng.normal(size=(18,2))
        u, _ = np.linalg.qr(columns)
        for z in range(types):
            for k in range(7):
                hij = blocks[:,:,k,z][np.ix_(ACTIVE,ACTIVE)]
                blocks[:,:,k,z] = u @ hij @ u.conj().T
            for matrix in (soc, o):
                matrix[:,:,z] = u @ matrix[:,:,z][np.ix_(ACTIVE,ACTIVE)] @ u.conj().T
    with (directory/'projection.bin').open('wb') as f:
        def write(arr, dtype):
            np.asarray(arr, dtype=dtype).ravel(order='F').tofile(f)
        write([count,types,7,order,int(periodic),layout,int(impurity),int(hoh)], 'i4')
        write(nn, 'i4')
        write([z+1 for _,_,z in sites], 'i4')
        write(np.array(refs)+1, 'i4')
        write(np.array(sites).T, 'f8')
        write([width,width,1], 'i4')
        for arr in (blocks, hall, o, 0*o, soc):
            write(arr, 'c16')
        f.write(observable.encode().ljust(30)+projection.encode().ljust(16)+pol.encode().ljust(10))
    subprocess.run([str(exe), 'projection'], cwd=directory, check=True,
                   env={**os.environ, 'OMP_NUM_THREADS':'1'}, stdout=subprocess.DEVNULL)
    raw = np.fromfile(directory/'moments.bin', complex)
    if observable != "spin_torque_onsite":
        return raw.reshape((order,order,18,2), order='F').sum(axis=2)
    if layout:
        actual = raw.reshape((18,18,order,order,2), order='F')[np.ix_(ACTIVE,ACTIVE)]
    else:
        actual = raw.reshape((order,order,18,2), order='F')[:,:,ACTIVE,:]
    values, U = np.linalg.eigh(h)
    poly = [(U*np.cos(n*np.arccos(values/(8.5/float(np.float32(1.7))))))@U.conj().T for n in range(order)]
    traces = np.zeros((order,order,2), complex)
    error = 0.
    for it, site in enumerate(refs):
        sl = slice(2*site,2*site+2)
        for n, tn in enumerate(poly):
            for m, tm in enumerate(poly):
                if projection == 'local':
                    oracle = a[sl,sl]@tn[sl,:]@v@tm[:,sl]
                else:
                    oracle = tm[sl,:]@a@tn@v[:,sl]
                measured = actual[:,:,n,m,it] if layout else actual[n,m,:,it]
                expected = oracle if layout else np.diag(oracle)
                error = max(error, np.max(abs(expected-measured)))
                traces[n,m,it] = np.trace(oracle)
    assert error < 1e-12, error
    return traces, error


@pytest.mark.parametrize('hoh', [False, True])
@pytest.mark.parametrize('layout', [0, 1])
def test_periodic_slab(tmp_path, transport_executable, hoh, layout):
    def calculate(rep, projection):
        return slab(tmp_path, transport_executable, 3, True, rep, 5, hoh, layout, projection)[0]
    seed = calculate(1, 'seed')
    local = calculate(1, 'local')
    np.testing.assert_allclose(local, calculate(0,'local'), atol=1e-12)
    np.testing.assert_allclose(seed.sum(2), local.sum(2), atol=1e-12, rtol=0)
    assert abs(seed-local).max() > 1e-5
    if not hoh:
        assert abs(local[:,:,0]).max() == 0


def test_open_slab(tmp_path, transport_executable):
    reference, _ = slab(tmp_path, transport_executable, 7, True, 3, 3, False, 0, 'local')
    errors = []
    for width, rep in [(3,0), (3,1), (5,2)]:
        seed, _ = slab(tmp_path, transport_executable, width, False, rep, 3, False, 0, 'seed')
        local, _ = slab(tmp_path, transport_executable, width, False, rep, 3, False, 0, 'local')
        errors.append([abs(local-reference).max(), abs((seed-local).sum(2)).max()])
    assert errors[0][0] > 1e-5
    assert max(errors[-1]) < 1e-12

def response_table(directory):
    with (directory/"local_projection_by_type.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert {r["element"] for r in rows} == {"M"}
    assert len({r["representative_site"] for r in rows}) == 2
    layers = np.array([[complex(float(r["real"]),float(r["imag"])) for r in rows if int(r["type"])==t]
                       for t in (1,2)])
    energy, re, im = np.loadtxt(directory/"response_total.out").T
    total = re+1j*im
    assert abs(total-layers.sum(axis=0)).max() <= 1e-12*abs(layers).max()
    return total, layers


@pytest.mark.parametrize("observable,pol", [
    ("spin_torque_onsite","x"), ("spin_accumulation","x"), ("orbital_accumulation","z"),
    ("quadrupole_accumulation","Qxy"), ("oap_accumulation","Qxy"), ("octupole_accumulation","Mx3y2r2")])
def test_slab_spectra(tmp_path, transport_executable, observable, pol):
    totals = []
    profiles = []
    for projection in ("seed","local"):
        slab(tmp_path, transport_executable, 3, True, 1, 7, True, 0, projection, observable, pol)
        if observable != "spin_torque_onsite" and projection == "seed":
            re, im = np.loadtxt(tmp_path/"response_total.out",usecols=(1,2)).T
            totals.append(re+1j*im)
            # Legacy filenames collide on the deliberately repeated label. The
            # last site's file still supplies one independent profile comparison.
            re, im = np.loadtxt(tmp_path/"M_response.out",usecols=(1,2)).T
            profiles.append(re+1j*im)
        else:
            total, layers = response_table(tmp_path)
            totals.append(total)
            profiles.append(layers[1])
    scale = max(abs(x).max() for x in profiles)
    assert abs(totals[0]-totals[1]).max() < 1e-11*scale
    assert abs(profiles[0]-profiles[1]).max() > 1e-5*scale


def test_nonreciprocal_block_rejected(tmp_path, transport_executable):
    slab(tmp_path, transport_executable, 3, True, 1, 3, True, 0, 'local')
    path = tmp_path/'projection.bin'
    record = bytearray(path.read_bytes())
    sites, types, neighbors = np.frombuffer(record, dtype='i4', count=3)
    # Skip the header, neighbor and type maps, coordinates and supercell counts.
    offset = 4*(8 + sites*neighbors + sites + types + 3) + 8*3*sites
    hopping = np.frombuffer(record, dtype='c16', offset=offset, count=18*18*neighbors*types)
    hopping.reshape((18,18,neighbors,types), order='F')[0,0,1,0] += .001j
    path.write_bytes(record)
    run = subprocess.run([transport_executable,'projection'],cwd=tmp_path,capture_output=True,
                         env={**os.environ,"OMP_NUM_THREADS":"1"},text=True)
    assert run.returncode != 0
    assert "requires Hermitian transport" in run.stdout+run.stderr

@pytest.mark.parametrize("hoh", [False,True])
def test_site_specific_field(tmp_path, transport_executable, hoh):
    # Only one occurrence of type 1 has a transverse local field in hall.
    for projection in ["seed","local"]:
        slab(tmp_path, transport_executable, 3, True, 1, 7, hoh, 1, projection, impurity=True)
