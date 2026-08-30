#!/usr/bin/env python3
"""Independent fixed-potential covariance gates for GBT-WP03.

The production Fortran helpers and WP02 gate establish the primitive link.  This
test deliberately reconstructs the ordinary one-sublattice reference with
independent dense algebra and checks the composite terms one at a time:

* HOH: eeo(k) h(k) and the second-order Hamiltonian;
* CCOR: scalar and pair-surface combinations of D and Ddot;
* onsite Hubbard-U: q-independent local addition;
* generalized overlap: the incomplete metric remains rejected by source guard.

The ordinary reference never uses a GBT helper or a completed-object phase.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
TOL = 2.0e-12
NEGATIVE_TOL = 1.0e-3
NORB = 3
NB = 2 * NORB
PI = np.pi

# A one-sublattice directed bond table.  The reverse entries are explicitly
# supplied, so every real-space test sees a finite Hermitian operator rather
# than relying on a self-conjugate scalar fixture.
BONDS = np.array(
    [
        [0, 0, 0],
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
    ],
    dtype=float,
)
Q_POINTS = np.array(
    [
        [0.173, -0.117, 0.086],
        [-0.271, 0.314, 0.119],
        [0.223, 0.071, -0.307],
    ]
)
K_POINTS = np.array(
    [
        [0.137, -0.219, 0.083],
        [-0.241, 0.193, 0.157],
        [0.311, 0.047, -0.283],
    ]
)


def relerr(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.linalg.norm(actual - expected) / max(1.0, float(np.linalg.norm(expected))))


def eigerr(actual: np.ndarray, expected: np.ndarray) -> float:
    return float(np.max(np.abs(np.linalg.eigvalsh(actual) - np.linalg.eigvalsh(expected))))


def endpoint_link(q: np.ndarray, displacement: np.ndarray) -> np.ndarray:
    """Independent common-frame endpoint link in the WP02 convention."""
    alpha = 2.0 * PI * float(np.dot(q, displacement))
    return np.diag([np.exp(-0.5j * alpha), np.exp(0.5j * alpha)])


def endpoint_w(w0: np.ndarray, w1: np.ndarray, sign: float = 1.0) -> np.ndarray:
    return np.diag(np.concatenate((w0 + sign * w1, w0 - sign * w1))).astype(np.complex128)


def orbital_spin_diag(values: np.ndarray) -> np.ndarray:
    return np.diag(np.tile(values, 2)).astype(np.complex128)


def primitive(raw: np.ndarray, q: np.ndarray, displacement: np.ndarray, w0: np.ndarray, w1: np.ndarray) -> np.ndarray:
    link = endpoint_link(q, displacement)
    return endpoint_w(w0, w1) @ np.kron(link, raw) @ endpoint_w(w0, w1)


def make_fixture() -> dict[str, object]:
    rng = np.random.default_rng(60303026)
    raw = np.zeros((len(BONDS), NORB, NORB), dtype=np.complex128)
    sdot = np.zeros_like(raw)
    raw[0] = 0.31 * (rng.normal(size=(NORB, NORB)) + 1j * rng.normal(size=(NORB, NORB)))
    raw[0] = 0.5 * (raw[0] + raw[0].conj().T)
    sdot[0] = 0.23 * (rng.normal(size=(NORB, NORB)) + 1j * rng.normal(size=(NORB, NORB)))
    sdot[0] = 0.5 * (sdot[0] + sdot[0].conj().T)
    for forward, reverse in ((1, 2), (3, 4)):
        raw[forward] = 0.12 * (rng.normal(size=(NORB, NORB)) + 1j * rng.normal(size=(NORB, NORB)))
        raw[reverse] = raw[forward].conj().T
        sdot[forward] = 0.07 * (rng.normal(size=(NORB, NORB)) + 1j * rng.normal(size=(NORB, NORB)))
        sdot[reverse] = sdot[forward].conj().T

    return {
        "raw": raw,
        "sdot": sdot,
        "w0": np.array([1.08, 0.74, 0.91]),
        "w1": np.array([0.23, -0.17, 0.11]),
        "o0": np.array([0.031, 0.047, 0.019]),
        "o1": np.array([-0.012, 0.009, 0.015]),
        "e0": np.array([-0.19, 0.27, 0.11]),
        "e1": np.array([0.039, -0.013, 0.021]),
        "a0": np.array([0.017, 0.023, 0.011]),
        "a1": np.array([-0.031, 0.027, 0.019]),
        "lam": 0.37,
        "pair_lam": np.array([0.61, -0.23]),
        "u": np.array(
            [
                [0.081, 0.014, -0.009],
                [0.014, 0.063, 0.012],
                [-0.009, 0.012, 0.097],
            ]
        ),
    }


def gbt_first(fixture: dict[str, object], q: np.ndarray, k: np.ndarray) -> np.ndarray:
    raw = fixture["raw"]
    h = np.zeros((NB, NB), dtype=np.complex128)
    for displacement, block in zip(BONDS, raw):
        phase = np.exp(2.0j * PI * float(np.dot(k, displacement)))
        h += primitive(block, q, displacement, fixture["w0"], fixture["w1"]) * phase
    return h


def gbt_primitive_bonds(fixture: dict[str, object], q: np.ndarray, key: str) -> list[np.ndarray]:
    return [
        primitive(block, q, displacement, fixture["w0"], fixture["w1"])
        for displacement, block in zip(BONDS, fixture[key])
    ]


def local_blocks(fixture: dict[str, object]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    overlap_bar = orbital_spin_diag(fixture["o0"] + fixture["o1"])
    overlap_bar[3:, 3:] = np.diag(fixture["o0"] - fixture["o1"])
    e_nu = orbital_spin_diag(fixture["e0"] + fixture["e1"])
    e_nu[3:, 3:] = np.diag(fixture["e0"] - fixture["e1"])
    u = np.zeros((NB, NB), dtype=np.complex128)
    u[:3, :3] = fixture["u"]
    u[3:, 3:] = fixture["u"] * 0.83
    return overlap_bar, e_nu, u


def fourier(blocks: list[np.ndarray], k: np.ndarray) -> np.ndarray:
    result = np.zeros_like(blocks[0])
    for displacement, block in zip(BONDS, blocks):
        result += block * np.exp(2.0j * PI * float(np.dot(k, displacement)))
    return result


def ordinary_sector(fixture: dict[str, object], k: np.ndarray, spin: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sign = 1.0 if spin == 0 else -1.0
    w = fixture["w0"] + sign * fixture["w1"]
    o = np.diag(fixture["o0"] + sign * fixture["o1"]).astype(np.complex128)
    e = np.diag(fixture["e0"] + sign * fixture["e1"]).astype(np.complex128)
    h_blocks = [w[:, None] * block * w[None, :] for block in fixture["raw"]]
    sdot_blocks = [w[:, None] * block * w[None, :] for block in fixture["sdot"]]
    h = fourier(h_blocks, k)
    sdot = fourier(sdot_blocks, k)
    return h, o, e


def block_diag(up: np.ndarray, down: np.ndarray) -> np.ndarray:
    result = np.zeros((NB, NB), dtype=np.complex128)
    result[:NORB, :NORB] = up
    result[NORB:, NORB:] = down
    return result


def hoh_gbt(fixture: dict[str, object], q: np.ndarray, k: np.ndarray, include_u: bool = False) -> np.ndarray:
    overlap_bar, e_nu, u = local_blocks(fixture)
    h = gbt_first(fixture, q, k)
    eeo = fourier([bond @ overlap_bar for bond in gbt_primitive_bonds(fixture, q, "raw")], k)
    result = e_nu + h - eeo @ h
    if include_u:
        result = result + u
    return result


def hoh_ordinary(fixture: dict[str, object], k: np.ndarray, include_u: bool = False) -> np.ndarray:
    overlap_bar, e_nu, u = local_blocks(fixture)
    h_up, o_up, e_up = ordinary_sector(fixture, k, 0)
    h_down, o_down, e_down = ordinary_sector(fixture, k, 1)
    up = e_up + h_up - (h_up @ o_up) @ h_up
    down = e_down + h_down - (h_down @ o_down) @ h_down
    result = block_diag(up, down)
    if include_u:
        result = result + u
    # Keep this assertion tied to the local construction above: the ordinary
    # reference uses spin sectors, while the GBT side uses the packed onsite
    # matrices from the production HOH contract.
    assert relerr(block_diag(e_up, e_down), e_nu) < TOL
    assert relerr(block_diag(o_up, o_down), overlap_bar) < TOL
    return result


def cc_orbital_blocks(fixture: dict[str, object], q: np.ndarray, pair_surface: bool) -> list[np.ndarray]:
    a1 = orbital_spin_diag(fixture["a1"])
    lam = fixture["lam"]
    pair = fixture["pair_lam"]
    scalar_l = 0.5 * float(pair[0] + pair[1])
    scalar_d = 0.5 * float(pair[0] - pair[1])
    lambda_block = np.diag(
        np.tile(
            [scalar_l + scalar_d, scalar_l + scalar_d, scalar_l + scalar_d],
            2,
        )
    ).astype(np.complex128)
    if pair_surface:
        # Equal endpoints are used here; the source still executes the
        # endpoint-ordered Lambda_i/Lambda_j expression.
        lambda_block = np.diag(
            np.tile([pair[0], pair[0], pair[0], pair[1], pair[1], pair[1]], 1)
        ).astype(np.complex128)
    blocks: list[np.ndarray] = []
    for displacement, s_block, sdot_block in zip(
        BONDS,
        gbt_primitive_bonds(fixture, q, "raw"),
        gbt_primitive_bonds(fixture, q, "sdot"),
    ):
        if pair_surface:
            d = 0.5 * (lambda_block @ sdot_block + sdot_block @ lambda_block)
            d += lambda_block @ a1 @ s_block + s_block @ a1 @ lambda_block
        else:
            d = lam * (sdot_block + a1 @ s_block + s_block @ a1)
        blocks.append(d)

    _, _, _ = local_blocks(fixture)
    onsite_w = endpoint_w(fixture["w0"], fixture["w1"])
    onsite_a0 = orbital_spin_diag(fixture["a0"])
    onsite_a0[NORB:, NORB:] = np.diag(fixture["a0"])
    if pair_surface:
        blocks[0] = blocks[0] + lambda_block @ onsite_a0 @ onsite_w @ onsite_w
    else:
        blocks[0] = blocks[0] + lam * onsite_a0 @ onsite_w @ onsite_w
    return blocks


def ordinary_ccor_sector(fixture: dict[str, object], k: np.ndarray, spin: int, pair_surface: bool) -> np.ndarray:
    sign = 1.0 if spin == 0 else -1.0
    w = fixture["w0"] + sign * fixture["w1"]
    a1 = np.diag(fixture["a1"]).astype(np.complex128)
    a0 = np.diag(fixture["a0"]).astype(np.complex128)
    lam = fixture["lam"]
    pair = fixture["pair_lam"]
    hcc_blocks: list[np.ndarray] = []
    for raw, sdot in zip(fixture["raw"], fixture["sdot"]):
        d = w[:, None] * raw * w[None, :]
        ddot = w[:, None] * sdot * w[None, :]
        if pair_surface:
            pair_lam = pair[spin]
            value = pair_lam * (ddot + a1 @ d + d @ a1)
        else:
            value = lam * (ddot + a1 @ d + d @ a1)
        hcc_blocks.append(value)
    hcc_blocks[0] = hcc_blocks[0] + (pair[spin] if pair_surface else lam) * a0 @ np.diag(w * w)
    return fourier(hcc_blocks, k)


def gbt_ccor(fixture: dict[str, object], q: np.ndarray, k: np.ndarray, pair_surface: bool) -> np.ndarray:
    return fourier(cc_orbital_blocks(fixture, q, pair_surface), k)


def source_contract() -> list[str]:
    failures: list[str] = []
    build = (ROOT / "source/hamiltonian_build.f90").read_text(encoding="utf-8")
    ccor = (ROOT / "source/hamiltonian_ccor.f90").read_text(encoding="utf-8")
    lifecycle = (ROOT / "source/reciprocal_lifecycle.f90").read_text(encoding="utf-8")
    reciprocal = (ROOT / "source/reciprocal_fourier.f90").read_text(encoding="utf-8")
    gbt_start = build.index("module subroutine build_gbt_bulkham")
    gbt_end = build.index("end subroutine build_gbt_bulkham", gbt_start)
    gbt = build[gbt_start:gbt_end]
    ccor_start = ccor.index("module subroutine build_ccor_pair_block_gbt")
    ccor_end = ccor.index("end subroutine build_ccor_pair_block_gbt", ccor_start)
    ccor_gbt = ccor[ccor_start:ccor_end]

    if gbt.count("call gbt_endpoint_link") != 1:
        failures.append("HOH/CCOR GBT builder does not own exactly one primitive endpoint-link calculation")
    if "call this%build_obarm()" not in gbt or "this%eeo(:, :, m, ntype)" not in gbt:
        failures.append("HOH is not derived from linked ee and onsite obarm")
    if "call zgemm('N','N',nmat,nmat,nmat,cone,workspace%eeo" not in reciprocal:
        failures.append("reciprocal HOH is not formed as eeo(k)*h(k)")
    if ccor_gbt.count("call gbt_lift_orbital_block") != 2:
        failures.append("CCOR does not link S and Sdot exactly once each")
    if "raw_sdot(i, j) = normalize_ccor_sdot(this, this%lattice%sdot(j, i, m, ino))" not in gbt:
        failures.append("GBT CCOR Sdot does not use the S-directed neighbor slot/order")
    if "raw_orb(i, j) = this%lattice%sbar(j, i, m, ino)" not in gbt:
        failures.append("GBT CCOR S does not use the S-directed neighbor slot/order")
    if "gbt_endpoint_link" in ccor_gbt or "q_ss" in ccor_gbt:
        failures.append("CCOR applies a second GBT phase after the shared primitive link")
    if "gbt_single_q with intersite Hubbard-V is unsupported." not in gbt:
        failures.append("intersite Hubbard-V is not explicitly rejected")
    if "gbt_single_q with local_axis is not yet audited" not in gbt:
        failures.append("unaudited local-axis GBT path is not explicitly rejected")
    if "gbt_single_q with SOC is unsupported." not in gbt:
        failures.append("SOC is not explicitly rejected in the GBT builder")
    if "if (this%hubbard_u_general_check) call this%calculate_hubbard_u_potential_general()" not in gbt:
        failures.append("onsite Hubbard-U is not built in the GBT path")
    if "G_ii=I" not in gbt:
        failures.append("onsite Hubbard-U lacks the d=0 identity-link contract")
    if "is unsupported: the available generalized-overlap" not in lifecycle:
        failures.append("incomplete generalized-overlap modes are not rejected")
    if "call this%validate_nonzero_q_gbt('reciprocal%build_kspace_hamiltonian')" not in reciprocal:
        failures.append("reciprocal production path does not invoke the GBT feature guard")
    return failures


def run_checks() -> tuple[dict[str, float], list[str]]:
    fixture = make_fixture()
    residuals: dict[str, float] = {}
    failures: list[str] = []

    # HOH is audited alone first.  This covers q=0, shifted-k, reverse bonds,
    # reciprocal Hermiticity, and eigenvalues without any CCOR/U addition.
    for q in [np.zeros(3), Q_POINTS[0], Q_POINTS[1]]:
        for k in K_POINTS:
            actual = hoh_gbt(fixture, q, k)
            expected = hoh_ordinary(fixture, k) if np.allclose(q, 0.0) else block_diag(
                hoh_ordinary(fixture, k - 0.5 * q)[:NORB, :NORB],
                hoh_ordinary(fixture, k + 0.5 * q)[NORB:, NORB:],
            )
            residuals["HOH shifted-k/q0 matrix"] = max(
                residuals.get("HOH shifted-k/q0 matrix", 0.0), relerr(actual, expected)
            )
            residuals["HOH shifted-k/q0 eigenvalues"] = max(
                residuals.get("HOH shifted-k/q0 eigenvalues", 0.0), eigerr(actual, expected)
            )
            residuals["HOH reciprocal Hermiticity"] = max(
                residuals.get("HOH reciprocal Hermiticity", 0.0), relerr(actual, actual.conj().T)
            )

    hoh_bonds = gbt_primitive_bonds(fixture, Q_POINTS[0], "raw")
    residuals["HOH primitive reverse bond"] = max(
        relerr(hoh_bonds[2], hoh_bonds[1].conj().T),
        relerr(hoh_bonds[4], hoh_bonds[3].conj().T),
    )

    # CCOR is then audited independently for both implemented VMT forms.
    for mode in (False, True):
        label = "CCOR pair-surface" if mode else "CCOR scalar"
        for q in [np.zeros(3), Q_POINTS[0], Q_POINTS[2]]:
            for k in K_POINTS:
                actual = gbt_ccor(fixture, q, k, mode)
                expected = (
                    block_diag(
                        ordinary_ccor_sector(fixture, k, 0, mode),
                        ordinary_ccor_sector(fixture, k, 1, mode),
                    )
                    if np.allclose(q, 0.0)
                    else block_diag(
                        ordinary_ccor_sector(fixture, k - 0.5 * q, 0, mode),
                        ordinary_ccor_sector(fixture, k + 0.5 * q, 1, mode),
                    )
                )
                residuals[f"{label} shifted-k/q0 matrix"] = max(
                    residuals.get(f"{label} shifted-k/q0 matrix", 0.0), relerr(actual, expected)
                )
                residuals[f"{label} shifted-k/q0 eigenvalues"] = max(
                    residuals.get(f"{label} shifted-k/q0 eigenvalues", 0.0), eigerr(actual, expected)
                )
                residuals[f"{label} reciprocal Hermiticity"] = max(
                    residuals.get(f"{label} reciprocal Hermiticity", 0.0), relerr(actual, actual.conj().T)
                )
        bonds = cc_orbital_blocks(fixture, Q_POINTS[0], mode)
        residuals[f"{label} reverse bond"] = max(
            residuals.get(f"{label} reverse bond", 0.0), relerr(bonds[2], bonds[1].conj().T), relerr(bonds[4], bonds[3].conj().T)
        )

    # Composite integration is a separate check: the same primitive link may
    # feed HOH and CCOR, but neither completed object receives another phase.
    for q, k in zip(Q_POINTS, K_POINTS):
        actual = hoh_gbt(fixture, q, k) + gbt_ccor(fixture, q, k, False)
        expected = block_diag(
            hoh_ordinary(fixture, k - 0.5 * q)[:NORB, :NORB] + ordinary_ccor_sector(fixture, k - 0.5 * q, 0, False),
            hoh_ordinary(fixture, k + 0.5 * q)[NORB:, NORB:] + ordinary_ccor_sector(fixture, k + 0.5 * q, 1, False),
        )
        residuals["HOH+CCOR composition shifted-k eigenvalues"] = max(
            residuals.get("HOH+CCOR composition shifted-k eigenvalues", 0.0), eigerr(actual, expected)
        )

    # Onsite U is local: it is independent of q and preserves the same
    # shifted-k identity.  This is intentionally fixed-potential evidence,
    # not a claim about self-consistent U feedback.
    for q, k in zip(Q_POINTS, K_POINTS):
        residuals["onsite Hubbard-U shifted-k"] = max(
            residuals.get("onsite Hubbard-U shifted-k", 0.0),
            relerr(hoh_gbt(fixture, q, k, include_u=True) - hoh_gbt(fixture, q, k), local_blocks(fixture)[2]),
        )

    failures.extend(source_contract())
    for label, value in residuals.items():
        if value >= TOL:
            failures.append(f"{label} residual {value:.3e} is not below {TOL:.1e}")

    # Negative controls ensure that the test detects the two most likely
    # composite regressions: a second phase on eeo/CCOR and a missing half q.
    q, k = Q_POINTS[0], K_POINTS[0]
    correct = hoh_gbt(fixture, q, k)
    extra_phase = fourier(
        [bond @ local_blocks(fixture)[0] * np.exp(-1j * PI * float(np.dot(q, displacement)))
         for displacement, bond in zip(BONDS, gbt_primitive_bonds(fixture, q, "raw"))], k
    )
    wrong_hoh = local_blocks(fixture)[1] + gbt_first(fixture, q, k) - extra_phase @ gbt_first(fixture, q, k)
    wrong_half = hoh_gbt(fixture, 2.0 * q, k)
    cc_correct = gbt_ccor(fixture, q, k, False)
    cc_extra_phase = fourier(
        [block * np.exp(-1j * PI * float(np.dot(q, displacement)))
         for displacement, block in zip(BONDS, cc_orbital_blocks(fixture, q, False))], k
    )
    if relerr(cc_extra_phase, cc_correct) <= NEGATIVE_TOL:
        failures.append("negative control: extra CCOR phase was not detected")
    if relerr(wrong_hoh, correct) <= NEGATIVE_TOL:
        failures.append("negative control: extra HOH phase was not detected")
    if relerr(wrong_half, correct) <= NEGATIVE_TOL:
        failures.append("negative control: missing HOH half phase was not detected")
    residuals["negative extra HOH phase"] = relerr(wrong_hoh, correct)
    residuals["negative extra CCOR phase"] = relerr(cc_extra_phase, cc_correct)
    residuals["negative HOH half-phase"] = relerr(wrong_half, correct)
    return residuals, failures


def main() -> int:
    residuals, failures = run_checks()
    for label in sorted(residuals):
        print(f"{label}: {residuals[label]:.6e}")
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    print("WP03 composite covariance oracle PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
