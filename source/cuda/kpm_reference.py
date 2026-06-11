"""Modular Chebyshev-KPM driver for HDF5 tight-binding Hamiltonians.

Architecture
------------
                 read_hdf5_hamiltonian
                          |
                  HoppingModel (R-blocks + n_orb + H0)
                          |
        +--------- make_operator(backend, ...) ----------+
        |            |              |                    |
   CSROperator  BSROperator   FFTOperator         ConvOperator
        |            |              |                    |
        +------------+------+------+--------------------+
                            |
              LinearOperator.apply(v)  <- the ONLY thing KPM sees
                            |
          estimate_bounds -> scaled() -> kpm_moments -> reconstruct_dos

Extending
---------
* New matvec backend: subclass LinearOperator, implement apply() and
  scaled(), register with @register_backend("name").
* New disorder/defect type: subclass Perturbation, implement apply() (and
  diagonal() if it has one, so structured backends can fold it).
* Everything downstream (bounds, KPM, DOS, plotting) is backend-agnostic.

All structured backends are verified against CSR at startup when --verify
is passed (recommended once per new Hamiltonian file).
"""

import argparse
import time
import warnings

import h5py
import numpy as np
import torch
import torch.fft

warnings.filterwarnings("ignore")


# ============================================================================
# Model container + HDF5 I/O
# ============================================================================
class HoppingModel:
    """Translation-invariant part of H as dense n_orb x n_orb blocks per R."""

    def __init__(self, blocks, n_orb, H0=None):
        self.blocks = blocks            # dict {(Rx,Ry,Rz): (n_orb,n_orb) cplx}
        self.n_orb = n_orb
        self.H0 = H0 if H0 is not None else blocks.get((0, 0, 0))

    @property
    def r_max(self):
        return max(max(abs(c) for c in R) for R in self.blocks)

    def is_real(self):
        return all(np.abs(T.imag).max() < 1e-14 for T in self.blocks.values())

    def scaled(self, a, b):
        """Return blocks for H~ = (H - b I)/a (fold KPM rescaling into data)."""
        nb = {R: T / a for R, T in self.blocks.items()}
        z = (0, 0, 0)
        nb[z] = nb.get(z, 0) - (b / a) * np.eye(self.n_orb)
        return HoppingModel(nb, self.n_orb, self.H0)

    def to_coo_lists(self):
        Rs, ois, ojs, vals = [], [], [], []
        for R, T in self.blocks.items():
            ii, jj = np.nonzero(T)
            Rs.append(np.repeat([R], len(ii), axis=0))
            ois.append(ii)
            ojs.append(jj)
            vals.append(T[ii, jj])
        return (np.concatenate(Rs).astype(np.int64), np.concatenate(ois),
                np.concatenate(ojs), np.concatenate(vals).astype(np.complex128))


def read_hdf5_hamiltonian(filename, threshold=1e-8):
    """Read the canonical RS-LMTO HDF5 export into a HoppingModel."""
    print(f"Reading HDF5 Hamiltonian from {filename}...")
    with h5py.File(filename, "r") as h5:
        R = np.asarray(h5["R"][:], dtype=np.int64)
        i_orb = np.asarray(h5["i_orb"][:], dtype=np.int64)
        j_orb = np.asarray(h5["j_orb"][:], dtype=np.int64)
        H_val = np.asarray(h5["H"][:]).astype(np.complex128)

    n_orb = int(max(i_orb.max(), j_orb.max())) + 1

    H0 = np.zeros((n_orb, n_orb), dtype=complex)
    m0 = np.all(R == 0, axis=1)
    H0[i_orb[m0], j_orb[m0]] = H_val[m0]

    keep = np.abs(H_val) > threshold
    R, i_orb, j_orb, H_val = R[keep], i_orb[keep], j_orb[keep], H_val[keep]

    blocks = {}
    Ru, inv = np.unique(R, axis=0, return_inverse=True)
    for k, Rv in enumerate(Ru):
        m = inv == k
        T = np.zeros((n_orb, n_orb), dtype=np.complex128)
        T[i_orb[m], j_orb[m]] = H_val[m]
        blocks[tuple(int(x) for x in Rv)] = T

    print(f"Detected {n_orb} orbitals, {len(blocks)} R-shells, "
          f"{len(H_val)} non-zero hoppings (r_max={max(np.abs(Ru).max(axis=0))}).")
    return HoppingModel(blocks, n_orb, H0)


# ============================================================================
# Perturbations (the symmetry breaking) — extension point
# ============================================================================
class Perturbation:
    """Base class. apply(v) adds the perturbation's contribution to H v."""

    def diagonal(self, n_total):
        """Return a length-n_total complex diagonal, or None if not diagonal.
        Structured backends fold diagonals into a single fused FMA."""
        return None

    def apply(self, v):
        raise NotImplementedError

    def scaled(self, a):
        raise NotImplementedError


class AndersonDisorder(Perturbation):
    """Uncorrelated on-site disorder, uniform in [-W/2, W/2] on every site."""

    def __init__(self, W, n_total, seed=0, device="cpu"):
        rng = np.random.default_rng(seed)
        self.V = torch.from_numpy(
            rng.uniform(-W / 2, W / 2, size=n_total)
        ).to(torch.complex64).to(device)

    def diagonal(self, n_total):
        return self.V

    def apply(self, v):
        return self.V[:, None] * v

    def scaled(self, a):
        p = object.__new__(AndersonDisorder)
        p.V = self.V / a
        return p


class DiagonalFromArray(Perturbation):
    """Arbitrary user-supplied on-site potential (e.g. from .npy file)."""

    def __init__(self, V, device="cpu"):
        self.V = torch.as_tensor(V).to(torch.complex64).to(device)

    def diagonal(self, n_total):
        return self.V

    def apply(self, v):
        return self.V[:, None] * v

    def scaled(self, a):
        return DiagonalFromArray(self.V / a)


class SparseCorrection(Perturbation):
    """Local NON-diagonal defects (e.g. modified hoppings around impurities)
    as a small sparse matrix added on top of a structured backend.
    rows/cols/vals are global site indices into the supercell."""

    def __init__(self, rows, cols, vals, n_total, device="cpu"):
        idx = torch.from_numpy(np.stack([rows, cols]))
        v = torch.as_tensor(vals).to(torch.complex64)
        M = torch.sparse_coo_tensor(idx, v, (n_total, n_total)).coalesce()
        self.M = M.to_sparse_csr().to(device)

    def apply(self, v):
        return torch.sparse.mm(self.M, v)

    def scaled(self, a):
        p = object.__new__(SparseCorrection)
        p.M = (self.M * (1.0 / a))
        return p


# ============================================================================
# Operator backends
# ============================================================================
BACKENDS = {}


def register_backend(name):
    def deco(cls):
        BACKENDS[name] = cls
        cls.backend_name = name
        return cls
    return deco


class LinearOperator:
    """Abstract y = H v. Subclasses own the translation-invariant part;
    perturbations are applied on top (diagonal ones folded where possible)."""

    def __init__(self, model, N, perturbations=(), device="cpu"):
        self.model = model
        self.N = tuple(N)
        self.n_orb = model.n_orb
        self.n_cells = N[0] * N[1] * N[2]
        self.n_total = self.n_cells * model.n_orb
        self.device = device
        # Split perturbations into a single folded diagonal + general rest
        diag = torch.zeros(self.n_total, dtype=torch.complex64, device=device)
        rest = []
        have_diag = False
        for p in perturbations:
            d = p.diagonal(self.n_total)
            if d is not None:
                diag += d.to(device)
                have_diag = True
            else:
                rest.append(p)
        self.diag = diag if have_diag else None
        self.perturbations = rest

    # -- subclass API -------------------------------------------------------
    def apply_periodic(self, v):
        raise NotImplementedError

    # -- common -------------------------------------------------------------
    def apply(self, v):
        y = self.apply_periodic(v)
        if self.diag is not None:
            y = y + self.diag[:, None] * v
        for p in self.perturbations:
            y = y + p.apply(v)
        return y

    def scaled(self, a, b):
        """Operator for H~ = (H - b)/a with everything folded into the data,
        so the KPM loop body contains no extra scaling kernels."""
        new = type(self)(self.model.scaled(a, b), self.N,
                         perturbations=(), device=self.device)
        if self.diag is not None:
            new.diag = self.diag / a
        new.perturbations = [p.scaled(a) for p in self.perturbations]
        return new


def _supercell_coo(model, N):
    """Vectorized COO triplets of the periodic supercell (unique-R factored)."""
    R, oi, oj, val = model.to_coo_lists()
    Nx, Ny, Nz = N
    n_cells = Nx * Ny * Nz
    n_orb = model.n_orb
    gx, gy, gz = np.meshgrid(np.arange(Nx), np.arange(Ny), np.arange(Nz),
                             indexing="ij")
    cx, cy, cz = gx.ravel(), gy.ravel(), gz.ravel()
    c_i = (cx * Ny + cy) * Nz + cz
    Ru, inv = np.unique(R, axis=0, return_inverse=True)
    tx = (cx[:, None] + Ru[None, :, 0]) % Nx
    ty = (cy[:, None] + Ru[None, :, 1]) % Ny
    tz = (cz[:, None] + Ru[None, :, 2]) % Nz
    c_j = ((tx * Ny + ty) * Nz + tz)[:, inv]
    rows = (c_i[:, None] * n_orb + oi[None, :]).ravel()
    cols = (c_j * n_orb + oj[None, :]).ravel()
    vals = np.broadcast_to(val, (n_cells, len(val))).ravel()
    return rows, cols, vals


@register_backend("csr")
class CSROperator(LinearOperator):
    """Reference backend: generic sparse CSR SpMM. Always correct."""

    def __init__(self, model, N, perturbations=(), device="cpu"):
        super().__init__(model, N, perturbations, device)
        rows, cols, vals = _supercell_coo(model, N)
        idx = torch.from_numpy(np.stack([rows, cols]))
        v = torch.from_numpy(np.ascontiguousarray(vals)).to(torch.complex64)
        H = torch.sparse_coo_tensor(idx, v, (self.n_total, self.n_total))
        self.H = H.coalesce().to_sparse_csr().to(device)

    def apply_periodic(self, v):
        return torch.sparse.mm(self.H, v)


@register_backend("bsr")
class BSROperator(LinearOperator):
    """Block-sparse rows with native n_orb x n_orb dense tiles."""

    def __init__(self, model, N, perturbations=(), device="cpu"):
        super().__init__(model, N, perturbations, device)
        rows, cols, vals = _supercell_coo(model, N)
        idx = torch.from_numpy(np.stack([rows, cols]))
        v = torch.from_numpy(np.ascontiguousarray(vals)).to(torch.complex64)
        H = torch.sparse_coo_tensor(idx, v, (self.n_total, self.n_total))
        self.H = (H.coalesce().to_sparse_csr()
                  .to_sparse_bsr(blocksize=(self.n_orb, self.n_orb)).to(device))

    def apply_periodic(self, v):
        return torch.sparse.mm(self.H, v)


@register_backend("fft")
class FFTOperator(LinearOperator):
    """H_periodic v = IFFT[ H(k) @ FFT[v] ]. Cost independent of hopping
    range; operator storage is n_cells * n_orb^2 complex (usually < a few MB)."""

    def __init__(self, model, N, perturbations=(), device="cpu"):
        super().__init__(model, N, perturbations, device)
        Nx, Ny, Nz = N
        n_orb = model.n_orb
        Hk = np.zeros((Nx, Ny, Nz, n_orb, n_orb), dtype=np.complex128)
        kx = 2j * np.pi * np.fft.fftfreq(Nx)
        ky = 2j * np.pi * np.fft.fftfreq(Ny)
        kz = 2j * np.pi * np.fft.fftfreq(Nz)
        for (Rx, Ry, Rz), T in model.blocks.items():
            phase = np.exp(kx[:, None, None] * Rx
                           + ky[None, :, None] * Ry
                           + kz[None, None, :] * Rz)
            Hk += phase[..., None, None] * T
        self.Hk = torch.from_numpy(Hk).to(torch.complex64).to(device)

    def apply_periodic(self, v):
        Nx, Ny, Nz = self.N
        n_rhs = v.shape[1]
        vg = v.reshape(Nx, Ny, Nz, self.n_orb, n_rhs)
        vk = torch.fft.fftn(vg, dim=(0, 1, 2))
        yk = torch.einsum("xyzij,xyzjr->xyzir", self.Hk, vk)
        return torch.fft.ifftn(yk, dim=(0, 1, 2)).reshape(-1, n_rhs)


@register_backend("conv")
class ConvOperator(LinearOperator):
    """H_periodic v as a circular-padded grouped 3D convolution.
    Operator storage ~ KB; the KPM loop becomes a weight-shared deep CNN
    (torch.compile / CUDA-graph / channels_last / multi-GPU halo friendly).
    Complex algebra via 4 real convs; collapses to 2 if hoppings are real."""

    def __init__(self, model, N, perturbations=(), device="cpu"):
        super().__init__(model, N, perturbations, device)
        n_orb = model.n_orb
        r = model.r_max
        K = 2 * r + 1
        self.r_max = r
        Wr = torch.zeros(n_orb, n_orb, K, K, K)
        Wi = torch.zeros(n_orb, n_orb, K, K, K)
        for (Rx, Ry, Rz), T in model.blocks.items():
            Wr[:, :, r + Rx, r + Ry, r + Rz] = torch.from_numpy(T.real).float()
            Wi[:, :, r + Rx, r + Ry, r + Rz] = torch.from_numpy(T.imag).float()
        self.real_hops = model.is_real()
        self.Wr = Wr.to(device)
        self.Wi = None if self.real_hops else Wi.to(device)

    def apply_periodic(self, v):
        Nx, Ny, Nz = self.N
        r = self.r_max
        n_rhs = v.shape[1]
        vg = v.reshape(Nx, Ny, Nz, self.n_orb, n_rhs).permute(4, 3, 0, 1, 2)
        pad = (r,) * 6
        conv = torch.nn.functional.conv3d
        xr = torch.nn.functional.pad(vg.real.contiguous(), pad, mode="circular")
        xi = torch.nn.functional.pad(vg.imag.contiguous(), pad, mode="circular")
        if self.real_hops:
            yr, yi = conv(xr, self.Wr), conv(xi, self.Wr)
        else:
            yr = conv(xr, self.Wr) - conv(xi, self.Wi)
            yi = conv(xr, self.Wi) + conv(xi, self.Wr)
        y = torch.complex(yr, yi).permute(2, 3, 4, 1, 0)
        return y.reshape(-1, n_rhs)


def make_operator(backend, model, N, perturbations=(), device="cpu"):
    if backend not in BACKENDS:
        raise ValueError(f"Unknown backend '{backend}'. "
                         f"Available: {sorted(BACKENDS)}")
    return BACKENDS[backend](model, N, perturbations, device)


# ============================================================================
# Backend-agnostic numerics
# ============================================================================
@torch.no_grad()
def estimate_spectral_bounds(op, num_iters=80, buffer=0.02):
    """Extreme eigenvalues of Hermitian op via Lanczos on op.apply."""
    v = torch.randn(op.n_total, 1, dtype=torch.complex64, device=op.device)
    v /= torch.norm(v)
    v_prev = torch.zeros_like(v)
    alphas, betas = [], []
    beta = 0.0
    for _ in range(num_iters):
        w = op.apply(v)
        alpha = torch.real(torch.sum(v.conj() * w)).item()
        w = w - alpha * v - beta * v_prev
        beta = torch.norm(w).item()
        alphas.append(alpha)
        if beta < 1e-10:
            break
        betas.append(beta)
        v_prev, v = v, w / beta
    T = (np.diag(alphas)
         + np.diag(betas[: len(alphas) - 1], 1)
         + np.diag(betas[: len(alphas) - 1], -1))
    ev = np.linalg.eigvalsh(T)
    half = 0.5 * (ev[-1] - ev[0]) * buffer + 1e-12
    return float(ev[-1] + half), float(ev[0] - half)


@torch.no_grad()
def kpm_moments(op_scaled, target_cell, n_moments):
    """Double-moment KPM. op_scaled must already represent H~ = (H-b)/a,
    so the loop body is exactly one operator application + one fused fma."""
    n_orb = op_scaled.n_orb
    device = op_scaled.device
    v0 = torch.zeros((op_scaled.n_total, n_orb), dtype=torch.complex64,
                     device=device)
    s = target_cell * n_orb
    v0[s:s + n_orb] = torch.eye(n_orb, dtype=torch.complex64, device=device)

    moments = torch.zeros((n_moments, n_orb), dtype=torch.float32, device=device)
    v_pp = v0
    mu0 = torch.sum(v0.real**2 + v0.imag**2, dim=0)
    moments[0] = mu0
    v_p = op_scaled.apply(v0)
    mu1 = torch.sum(v0.conj() * v_p, dim=0).real
    if n_moments > 1:
        moments[1] = mu1

    for n in range(1, (n_moments + 1) // 2):
        if 2 * n < n_moments:
            moments[2 * n] = 2.0 * torch.sum(
                v_p.real**2 + v_p.imag**2, dim=0) - mu0
        if 2 * n + 1 < n_moments:
            v_c = op_scaled.apply(v_p)
            v_c.mul_(2.0).sub_(v_pp)
            moments[2 * n + 1] = 2.0 * torch.sum(
                v_p.conj() * v_c, dim=0).real - mu1
            v_pp, v_p = v_p, v_c
    return moments.cpu().numpy().astype(np.float64)


def reconstruct_dos(moments, a, b, n_points=5000):
    """Jackson-damped Chebyshev sum (Fortran-matching kernel), evaluated as
    one recurrence-built Chebyshev matrix + a single GEMM."""
    if moments.ndim == 1:
        moments = moments[:, None]
    N = moments.shape[0]
    n = np.arange(N, dtype=np.float64)
    scale = N + 1.0
    theta = np.pi * n / scale
    g = ((N - n + 1.0) * np.cos(theta)
         + np.sin(theta) / np.tan(np.pi / scale)) / scale
    coeff = np.full(N, 2.0)
    coeff[0] = 1.0
    damped = moments * (g * coeff)[:, None]

    E_grid = np.linspace(b - a + 1e-4 * a, b + a - 1e-4 * a, n_points)
    x = (E_grid - b) / a
    T = np.empty((N, n_points))
    T[0] = 1.0
    if N > 1:
        T[1] = x
    for k in range(2, N):
        T[k] = 2.0 * x * T[k - 1] - T[k - 2]
    dos = T.T @ damped
    dos *= (1.0 / (np.pi * a * np.sqrt(1.0 - x**2)))[:, None]
    return E_grid, dos


# ============================================================================
# Verification helper
# ============================================================================
def verify_against_csr(op, model, N, perturbations, n_check=4, tol=5e-5):
    """Compare op.apply against the CSR reference on random vectors."""
    if isinstance(op, CSROperator):
        print("Backend is CSR (reference); nothing to verify.")
        return
    ref = make_operator("csr", model, N, perturbations, op.device)
    v = torch.randn(op.n_total, n_check, dtype=torch.complex64,
                    device=op.device)
    y_ref = ref.apply(v)
    err = (torch.norm(op.apply(v) - y_ref) / torch.norm(y_ref)).item()
    status = "OK" if err < tol else "FAILED"
    print(f"Verification vs CSR: rel.err = {err:.2e}  [{status}]")
    if err >= tol:
        raise RuntimeError("Backend disagrees with CSR reference.")
    del ref


# ============================================================================
# Main pipeline
# ============================================================================
def main():
    p = argparse.ArgumentParser(
        description="Modular HDF5 KPM engine with structured matvec backends")
    p.add_argument("--file", type=str, required=True)
    p.add_argument("--backend", type=str, default="fft",
                   choices=sorted(BACKENDS), help="matvec backend")
    p.add_argument("--nx", type=int, default=10)
    p.add_argument("--ny", type=int, default=10)
    p.add_argument("--nz", type=int, default=10)
    p.add_argument("--moments", type=int, default=500)
    p.add_argument("--points", type=int, default=2000)
    p.add_argument("--emin", type=float, default=None)
    p.add_argument("--emax", type=float, default=None)
    p.add_argument("--anderson", type=float, default=0.0,
                   help="on-site Anderson disorder strength W")
    p.add_argument("--disorder-file", type=str, default=None,
                   help=".npy file with a length n_total on-site potential")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--target-cell", type=int, default=0)
    p.add_argument("--verify", action="store_true",
                   help="check backend against CSR before running")
    p.add_argument("--no-plot", action="store_true")
    args = p.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    N = (args.nx, args.ny, args.nz)
    print(f"--- KPM Pipeline | backend={args.backend} | device={device} ---")

    # 1. Model
    model = read_hdf5_hamiltonian(args.file)
    n_total = N[0] * N[1] * N[2] * model.n_orb

    # 2. Perturbations (symmetry breaking)
    perts = []
    if args.anderson > 0:
        print(f"Adding Anderson disorder W={args.anderson} (seed={args.seed})")
        perts.append(AndersonDisorder(args.anderson, n_total,
                                      seed=args.seed, device=device))
    if args.disorder_file:
        print(f"Adding on-site potential from {args.disorder_file}")
        perts.append(DiagonalFromArray(np.load(args.disorder_file),
                                       device=device))

    # 3. Operator
    print(f"Building {args.backend} operator for {N} supercell "
          f"({n_total} sites)...")
    t0 = time.time()
    op = make_operator(args.backend, model, N, perts, device)
    print(f"  Operator built in {time.time() - t0:.2f}s")

    if args.verify:
        verify_against_csr(op, model, N, perts)

    # 4. Bounds
    if args.emin is not None and args.emax is not None:
        E_min, E_max = args.emin, args.emax
        print(f"Using manual spectral bounds: [{E_min:.3f}, {E_max:.3f}]")
    else:
        print("Estimating spectral bounds via Lanczos...")
        E_max, E_min = estimate_spectral_bounds(op)
        print(f"  Estimated Bounds: [{E_min:.3f}, {E_max:.3f}]")
    a, b = (E_max - E_min) / 2.0, (E_max + E_min) / 2.0

    # 5. KPM on the pre-scaled operator
    op_scaled = op.scaled(a, b)
    print(f"Running {args.moments}-moment KPM on all {model.n_orb} orbitals...")
    t0 = time.time()
    moments = kpm_moments(op_scaled, args.target_cell, args.moments)
    dt = time.time() - t0
    print(f"  Loop completed in {dt:.2f}s "
          f"({args.moments / max(dt, 1e-9):.1f} moments/s)")

    # 6. DOS
    energies, dos_all = reconstruct_dos(moments, a, b, n_points=args.points)

    half = model.n_orb // 2
    total_up = np.sum(dos_all[:, :half], axis=1)
    total_dn = np.sum(dos_all[:, half:], axis=1)
    data = np.column_stack((energies, total_up, total_dn, total_up + total_dn))
    np.savetxt("cheb_dos.dat", data, fmt="%.12e",
               header="energy DOS_up DOS_dw DOS_tot")
    print("Wrote cheb_dos.dat")

    if args.no_plot:
        return

    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 2]})
    cax = ax1.matshow(np.abs(model.H0), cmap="magma")
    fig.colorbar(cax, ax=ax1, fraction=0.046, pad=0.04)
    ax1.set_title(r"Intra-cell $|H(R=0)|$", pad=15)
    ax1.set_xlabel("Orbital Index $j$")
    ax1.set_ylabel("Orbital Index $i$")

    for o in range(half):
        ax2.plot(energies, dos_all[:, o], color="dodgerblue",
                 alpha=0.15, linewidth=0.8)
    for o in range(half, model.n_orb):
        ax2.plot(energies, -dos_all[:, o], color="crimson",
                 alpha=0.15, linewidth=0.8)
    ax2.plot(energies, total_up, color="mediumblue", linewidth=2,
             label="Total Spin Up")
    ax2.fill_between(energies, total_up, color="mediumblue", alpha=0.2)
    ax2.plot(energies, -total_dn, color="darkred", linewidth=2,
             label="Total Spin Down")
    ax2.fill_between(energies, -total_dn, color="darkred", alpha=0.2)
    ax2.axhline(0, color="black", linewidth=1.0)
    ax2.axvline(0, color="black", linestyle="--", alpha=0.5, label="$E_F$")
    ax2.set_title(f"Spin/Orbital Resolved LDOS "
                  f"(N={args.moments}, backend={args.backend})")
    ax2.set_xlabel("Energy")
    ax2.set_ylabel("Density of States")
    max_val = max(np.max(total_up), np.max(total_dn)) * 1.1
    ax2.set_ylim(-max_val, max_val)
    ax2.set_xlim(E_min * 0.8, E_max * 0.8)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper right")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()