"""Check cyclic KPM block ordering against complete finite-matrix traces."""
import numpy as np
import pytest


@pytest.mark.parametrize("hoh", [False, True])
def test_local_moment_blocks(hoh):
    rng = np.random.default_rng(583)
    def hermitian(size):
        z = rng.normal(size=(size, size)) + 1j*rng.normal(size=(size, size))
        return (z + z.conj().T)/12

    h = hermitian(8)
    h[:2, :2] = np.eye(2)
    if hoh:
        h -= h @ np.diag(np.linspace(.05, .2, 8)) @ h
    position = np.diag(np.repeat(np.arange(4), 2))
    velocity = -1j*(position @ h - h @ position)
    physical = h + hermitian(8)
    energy, states = np.linalg.eigh(physical)
    scale = 1.1*max(abs(energy))
    poly = [(states*np.cos(n*np.arccos(energy/scale))) @ states.conj().T for n in range(7)]
    def propagate(seed):
        vectors = [seed, physical @ seed/scale]
        for _ in range(2, len(poly)):
            vectors.append(2*physical @ vectors[-1]/scale - vectors[-2])
        return vectors

    sx = np.array([[0, .5], [.5, 0]])
    seeds = np.split(np.eye(8), 4, axis=1)
    vertices = []
    for z in seeds:
        block = z.T @ h @ z
        vertices.append(-1j*(sx @ block - block @ sx))
    total = sum(z @ a @ z.T for z, a in zip(seeds, vertices))
    seed = np.zeros((4, 7, 7), complex)
    local = np.zeros_like(seed)
    for i, (z, a) in enumerate(zip(seeds, vertices)):
        x = propagate(z)
        left = propagate(z @ a)
        for n, tn in enumerate(poly):
            for m, tm in enumerate(poly):
                oracle = a @ z.T @ tn @ velocity @ tm @ z
                two = left[n].conj().T @ velocity @ x[m]
                one = (x[m].conj().T @ velocity @ x[n] @ a).conj().T
                np.testing.assert_allclose([one, two], [oracle, oracle], atol=1e-12, rtol=0)
                local[i,n,m] = np.trace(tm @ z @ a @ z.T @ tn @ velocity)
                assert abs(local[i,n,m] - np.trace(one)) < 1e-12
                seed[i,n,m] = np.trace(z.T @ tm @ total @ tn @ velocity @ z)
    np.testing.assert_allclose(local.sum(axis=0), seed.sum(axis=0), atol=1e-12, rtol=0)
    assert abs(local - seed).max() > .001
    assert abs(local - local.swapaxes(1, 2)).max() > .001
    if not hoh:
        assert abs(local[0]).max() == 0
        assert abs(seed[0]).max() > .001


@pytest.mark.parametrize("field", [-.31,.31])
@pytest.mark.parametrize("orbitals", [4,9,16])
def test_field_sign_and_intersite_commutator(field, orbitals):
    pauli = [np.kron(s,np.eye(orbitals)) for s in
             (np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]))]
    hii = field*pauli[2]
    ident = np.eye(2*orbitals)
    def torque(s, h):
        return -.5j*(s@h-h@s)
    for s, expected in zip(pauli, [-field*pauli[1],field*pauli[0],0*hii]):
        np.testing.assert_allclose(torque(s,hii),expected,atol=1e-15)
        onsite = np.kron(np.diag([0.,1.]),hii)
        scalar_hop = np.kron(np.array([[0,1],[1,0]]),ident)
        magnetic_hop = np.kron(np.array([[0,1],[1,0]]), .12*pauli[2])
        q = np.kron(np.eye(2),s)
        np.testing.assert_allclose(torque(q,onsite+scalar_hop),torque(q,onsite))
        delta = torque(q,onsite+scalar_hop+magnetic_hop)-torque(q,onsite)
        np.testing.assert_allclose(delta,torque(q,magnetic_hop))
        field_change = .09*pauli[0]+.05*pauli[1]
        np.testing.assert_allclose(torque(s,hii+field_change)-torque(s,hii),torque(s,field_change))
