"""Complex-spherical orbital diagnostics; no projection of production data."""
import numpy as np

L = np.repeat(np.arange(3), [1, 3, 5])
MQUANT = np.concatenate([np.arange(-l, l + 1) for l in range(3)])
P = (-1.0) ** (L + MQUANT)
M = np.diag(P)
C = np.zeros((9, 9))
for k, (l, m) in enumerate(zip(L, MQUANT)):
    C[k, int(l * (l + 1) - m)] = (-1.0) ** m
INV = np.diag((-1.0) ** L)
SIGMA = np.array([np.eye(2), [[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)
UT = np.kron(1j * SIGMA[2], C)
# Mirror axial spin rotation is C2z. Overall phase is immaterial.
UM = np.kron(SIGMA[3], M)
UH = UM @ UT
SECTORS = ('00', '01', 'x0', 'x1', 'y0', 'y1', 'z0', 'z1')
MIRROR = dict(zip(SECTORS, [1, -1, 1, -1, 1, -1, -1, 1]))
BARE_T = dict(zip(SECTORS, [1, -1, -1, 1, -1, 1, -1, 1]))


def hcpx_v():
    """Numerical transformation entries in math.f90:hcpx (cart = V sph)."""
    v = np.zeros((9, 9), complex)
    q = 1 / np.sqrt(2)
    for i, j, a in [(1,1,1),(2,4,-q),(2,2,q),(3,4,1j*q),(3,2,1j*q),
                    (4,3,1),(5,5,1j*q),(5,9,-1j*q),(6,6,1j*q),
                    (6,8,1j*q),(7,6,q),(7,8,-q),(8,5,q),(8,9,q),(9,7,1)]:
        v[i-1,j-1] = a
    return v


V = hcpx_v()
# Cartesian real harmonics: s,x,y,z,xy,yz,xz,x2-y2,3z2-r2.
C2X = V.conj().T @ np.diag([1,1,-1,-1,-1,1,-1,1,1]) @ V
C2Z = INV @ M
OPS = {
    'C2x': (np.diag([1,-1,-1]), np.kron(SIGMA[1], C2X), False),
    'C2zT': (np.diag([-1,-1,1]), np.kron(SIGMA[3], C2Z) @ UT, True),
    'I': (-np.eye(3), np.kron(np.eye(2), INV), False),
}


def rorb(a):
    """R(A) = C A^t C; transpose, never adjoint/conjugate."""
    return C @ np.asarray(a).T @ C


def project(a, sign):
    return (a + sign * M @ a @ M) / 2


def split(a, b):
    return (a + rorb(b)) / 2, (a - rorb(b)) / 2


def pauli(g):
    a, b, c, d = g[:9,:9], g[:9,9:], g[9:,:9], g[9:,9:]
    return [(a+d)/2, (b+c)/2, 1j*(b-c)/2, (a-d)/2]


def spinor(gs):
    return sum(np.kron(s, g) for s, g in zip(SIGMA, gs))


def sectors(gij, gji):
    a, b = pauli(gij), pauli(gji)
    return {s: x for mu, aa, bb in zip('0xyz', a, b)
            for s, x in zip((mu+'0', mu+'1'), split(aa, bb))}


def residual(a, b):
    a, b = np.asarray(a), np.asarray(b)
    d = np.abs(a-b)
    norm = float(np.linalg.norm(a-b))
    scale = max(float(np.linalg.norm(a)), float(np.linalg.norm(b)))
    return dict(max=float(np.max(d)), mean=float(np.mean(d)), fro=norm,
                scale=scale, relative=norm/scale if scale else 0.0)


def trace_parts(di, dj, a, b, sign_a, sign_b):
    """Resolve both factors by mirror class; each contribution has s/p/d and p_a->p_b detail.

    Diagonal Delta is required: do not silently discard off-diagonal entries.
    """
    if np.linalg.norm(di-np.diag(np.diag(di))) > 1e-12 or np.linalg.norm(dj-np.diag(np.diag(dj))) > 1e-12:
        raise ValueError('trace decomposition requires diagonal Delta')
    weighted = np.diag(di)[:,None] * a * np.diag(dj)[None,:] * b.T
    mask_a = P[:,None]*P[None,:] == sign_a
    mask_b = P[:,None]*P[None,:] == sign_b
    out = {'total': float(weighted.sum().real), 'blocks': {}}
    for qa, ma in [('allowed',mask_a),('forbidden',~mask_a)]:
        for qb, mb in [('allowed',mask_b),('forbidden',~mask_b)]:
            out['blocks'][qa+'/'+qb] = float(weighted[ma & mb].sum().real)
    out['angular_parity'] = {}
    for la, na in enumerate('spd'):
        for lb, nb in enumerate('spd'):
            for p in (-1,1):
                for q in (-1,1):
                    mask = (L[:,None]==la)&(L[None,:]==lb)&(P[:,None]==p)&(P[None,:]==q)
                    out['angular_parity'][f'{na}{p:+}->{nb}{q:+}'] = float(weighted[mask].sum().real)
    return out


def leakage(a, sign):
    forbidden = P[:,None]*P[None,:] != sign
    bad = np.where(forbidden, a, 0)
    row,col = np.unravel_index(np.abs(bad).argmax(),bad.shape)
    return {'fro':float(np.linalg.norm(bad)), 'max':float(abs(bad[row,col])),
            'orbital_indices_1based':[int(row+1),int(col+1)],
            'angular_momenta':[int(L[row]),int(L[col])],
            'mirror_parities':[int(P[row]),int(P[col])],
            'element':[float(bad[row,col].real),float(bad[row,col].imag)]}
