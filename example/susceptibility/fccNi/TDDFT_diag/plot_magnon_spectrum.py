#!/usr/bin/env python3
"""Plot the transverse TDDFT loss spectrum written by the fcc-Ni example.

The program reads one selected `*_q??????_*dyson.dat` product directly, so no
post-processing format conversion is needed.  The horizontal coordinate is accumulated path
length in the *direct reciprocal coordinates* stored in the TDDFT metadata;
it is a path coordinate, not an absolute |q| in inverse Angstrom.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np


RY_TO_MEV = 13_605.693_122_994
# The TDDFT driver writes the legacy Dyson product plus optional raw and
# Goldstone-corrected pair-potential products.  All three share the q index.
Q_PATTERN = re.compile(r"_q(\d+)_(?:(?:pair(?:_corrected)?)_)?dyson\.dat$")
DYSON_SUFFIXES = {
    "corrected_pair": "_pair_corrected_dyson.dat",
    "raw_pair": "_pair_dyson.dat",
    "legacy": "_dyson.dat",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("."),
        help="directory containing fccNi_tddft_q*_dyson.dat (default: current directory)",
    )
    parser.add_argument(
        "--prefix",
        default="fccNi_tddft_sr",
        help="TDDFT output prefix used in *_q??????_*dyson.dat (default: fccNi_tddft_sr)",
    )
    parser.add_argument(
        "--q-path-file",
        type=Path,
        default=None,
        help="optional q-file with ! labels; defaults to fccNi_qpath.dat in --input-dir when present",
    )
    parser.add_argument("--site", type=int, default=1, help="one-based site index for site_loss (default: 1)")
    parser.add_argument(
        "--quantity",
        choices=("site_loss", "trace_loss"),
        default="site_loss",
        help="spectral weight to plot (default: site_loss)",
    )
    parser.add_argument(
        "--dyson-product",
        choices=("auto", *DYSON_SUFFIXES),
        default="auto",
        help=(
            "response product to plot; auto prefers corrected_pair, then raw_pair, then legacy "
            "(default: auto)"
        ),
    )
    parser.add_argument(
        "--energy-max-mev",
        type=float,
        default=None,
        help="discard frequencies above this energy in meV",
    )
    parser.add_argument(
        "--linear-colors",
        action="store_true",
        help="use a linear color scale instead of the default logarithmic scale",
    )
    parser.add_argument(
        "--normalize-each-q",
        action="store_true",
        help="divide every q column by its largest loss before plotting",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("fccNi_magnon_spectrum.png"),
        help="PNG output path (default: fccNi_magnon_spectrum.png)",
    )
    parser.add_argument(
        "--peak-output",
        type=Path,
        default=Path("fccNi_magnon_peaks.dat"),
        help="text peak-track output path (default: fccNi_magnon_peaks.dat)",
    )
    parser.add_argument("--show", action="store_true", help="show the plot interactively after writing it")
    return parser.parse_args()


def q_index(path: Path) -> int:
    match = Q_PATTERN.search(path.name)
    if match is None:
        raise ValueError(f"cannot obtain q index from {path.name}")
    return int(match.group(1))


def read_dyson_file(path: Path, quantity: str, site: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return direct q, frequency (Ry), and selected loss from one Dyson file."""
    q_direct: np.ndarray | None = None
    omega: list[float] = []
    loss: list[float] = []

    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if line.startswith("# q_direct ="):
            q_direct = np.fromstring(line.split("=", maxsplit=1)[1], sep=" ", dtype=float)
            continue
        fields = line.split()
        if quantity == "site_loss" and len(fields) == 4 and fields[0] == "site_loss":
            if int(fields[2]) == site:
                omega.append(float(fields[1]))
                loss.append(float(fields[3]))
        elif quantity == "trace_loss" and len(fields) == 3 and fields[0] == "trace_loss":
            omega.append(float(fields[1]))
            loss.append(float(fields[2]))

    if q_direct is None or q_direct.shape != (3,):
        raise ValueError(f"{path}: missing or malformed '# q_direct =' metadata")
    if not omega:
        suffix = f" for site {site}" if quantity == "site_loss" else ""
        raise ValueError(f"{path}: no {quantity} records{suffix}")

    omega_array = np.asarray(omega)
    loss_array = np.asarray(loss)
    order = np.argsort(omega_array)
    return q_direct, omega_array[order], loss_array[order]


def dyson_paths(input_dir: Path, prefix: str, product: str) -> list[Path]:
    """Return exactly one TDDFT response product, ordered by q index."""
    suffix = DYSON_SUFFIXES[product]
    paths = [path for path in input_dir.glob(f"{prefix}_q*_dyson.dat") if path.name.endswith(suffix)]
    if product == "legacy":
        paths = [path for path in paths if "_pair_" not in path.name]
    return sorted(paths, key=q_index)


def collect_spectrum(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    product = args.dyson_product
    if product == "auto":
        product = next(
            (candidate for candidate in ("corrected_pair", "raw_pair", "legacy")
             if dyson_paths(args.input_dir, args.prefix, candidate)),
            None,
        )
    if product is None:
        raise FileNotFoundError(f"no TDDFT Dyson files found in {args.input_dir.resolve()}")
    paths = dyson_paths(args.input_dir, args.prefix, product)
    if not paths:
        raise FileNotFoundError(
            f"no '{args.prefix}_q*{DYSON_SUFFIXES[product]}' files found in {args.input_dir.resolve()}"
        )

    q_points: list[np.ndarray] = []
    spectra: list[np.ndarray] = []
    reference_omega: np.ndarray | None = None
    for path in paths:
        q_direct, omega, loss = read_dyson_file(path, args.quantity, args.site)
        if reference_omega is None:
            reference_omega = omega
        elif omega.shape != reference_omega.shape or not np.allclose(omega, reference_omega, rtol=0.0, atol=1.0e-12):
            raise ValueError(f"{path}: frequency mesh differs from earlier q files")
        q_points.append(q_direct)
        spectra.append(loss)

    assert reference_omega is not None
    return np.asarray(q_points), reference_omega, np.asarray(spectra)


def path_coordinate(q_points: np.ndarray) -> np.ndarray:
    if len(q_points) == 1:
        return np.zeros(1)
    return np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(q_points, axis=0), axis=1))))


def read_q_labels(path: Path, q_points: np.ndarray, q_path: np.ndarray) -> list[tuple[float, str]]:
    """Read optional ``! LABEL`` markers and validate their q coordinates."""
    records = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not records:
        raise ValueError(f"{path}: empty q-path file")
    try:
        declared_count = int(records[0])
    except ValueError as error:
        raise ValueError(f"{path}: first record must be the q-point count") from error
    if declared_count != len(records) - 1 or declared_count != len(q_points):
        raise ValueError(f"{path}: q-point count does not match the Dyson output")

    labels: list[tuple[float, str]] = []
    for iq, record in enumerate(records[1:]):
        numeric, marker, label = record.partition("!")
        coordinates = np.fromstring(numeric, sep=" ", dtype=float)
        if coordinates.shape != (3,) or not np.allclose(coordinates, q_points[iq], rtol=0.0, atol=1.0e-10):
            raise ValueError(f"{path}: q point {iq + 1} does not match the Dyson metadata")
        if marker and label.strip():
            text = label.strip().upper()
            labels.append((q_path[iq], "Γ" if text == "GAMMA" else text))
    return labels


def write_peak_track(path: Path, q_points: np.ndarray, q_path: np.ndarray, omega_mev: np.ndarray, loss: np.ndarray) -> None:
    peak_indices = np.argmax(loss, axis=1)
    lines = [
        "# q_index q_path q1_direct q2_direct q3_direct peak_energy_meV peak_loss",
    ]
    for iq, peak_index in enumerate(peak_indices, start=1):
        q1, q2, q3 = q_points[iq - 1]
        lines.append(
            f"{iq:d} {q_path[iq - 1]:.12e} {q1:.12e} {q2:.12e} {q3:.12e} "
            f"{omega_mev[peak_index]:.12e} {loss[iq - 1, peak_index]:.12e}"
        )
    path.write_text("\n".join(lines) + "\n")


def plot_spectrum(args: argparse.Namespace) -> None:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
    except ModuleNotFoundError as error:
        raise SystemExit(
            "error: plotting requires matplotlib; install it with "
            "'python3 -m pip install --user matplotlib'"
        ) from error

    q_points, omega_ry, intensity = collect_spectrum(args)
    omega_mev = omega_ry * RY_TO_MEV
    if args.energy_max_mev is not None:
        keep = omega_mev <= args.energy_max_mev
        if not np.any(keep):
            raise ValueError("--energy-max-mev excludes the entire frequency mesh")
        omega_mev = omega_mev[keep]
        intensity = intensity[:, keep]
    peak_intensity = intensity.copy()
    negative_count = int(np.count_nonzero(intensity < 0.0))
    # A finite k mesh can produce a few small negative loss samples even
    # though the exact positive-frequency absorption is non-negative.  Do not
    # reflect them with abs(): that would create fictitious spectral weight.
    # Clipping is solely a plotting treatment; the raw values remain in the
    # Dyson files and the peak track is selected from those raw values.
    plot_intensity = np.maximum(intensity, 0.0)
    if negative_count:
        print(f"warning: clipped {negative_count} negative loss samples to zero for the color map")
    if args.normalize_each_q:
        maxima = np.max(plot_intensity, axis=1)
        nonzero = maxima > 0.0
        plot_intensity[nonzero] /= maxima[nonzero, np.newaxis]
    q_path = path_coordinate(q_points)
    q_label_file = args.q_path_file
    if q_label_file is None:
        candidate = args.input_dir / "fccNi_qpath.dat"
        q_label_file = candidate if candidate.exists() else None
    labels = read_q_labels(q_label_file, q_points, q_path) if q_label_file is not None else []

    fig, axis = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    positive = plot_intensity[plot_intensity > 0.0]
    if positive.size == 0:
        raise ValueError("loss is zero at every q and frequency")
    color_kwargs: dict[str, object] = {"cmap": "magma", "shading": "auto"}
    if args.linear_colors:
        color_kwargs["vmax"] = np.percentile(positive, 99.5)
    else:
        color_kwargs["norm"] = LogNorm(vmin=max(np.percentile(positive, 1.0), np.finfo(float).tiny), vmax=np.max(positive))
    image = axis.pcolormesh(q_path, omega_mev, plot_intensity.T, **color_kwargs)

    peak_indices = np.argmax(peak_intensity, axis=1)
    # axis.plot(q_path, omega_mev[peak_indices], "o-", color="cyan", lw=1.2, ms=4, label="largest loss")
    axis.set_xlabel("q-path coordinate (direct reciprocal coordinates)")
    axis.set_ylabel("energy (meV)")
    axis.set_title(f"fcc Ni transverse TDDFT {args.quantity.replace('_', ' ')}")
    if labels:
        axis.set_xticks([position for position, _ in labels])
        axis.set_xticklabels([label for _, label in labels])
        for position, _ in labels:
            axis.axvline(position, color="white", lw=0.6, alpha=0.45)
    color_label = "normalized loss" if args.normalize_each_q else "loss spectral weight"
    fig.colorbar(image, ax=axis, label=color_label)
    fig.savefig(args.output, dpi=180)
    write_peak_track(args.peak_output, q_points, q_path, omega_mev, peak_intensity)
    print(f"wrote {args.output}")
    print(f"wrote {args.peak_output}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    try:
        plot_spectrum(parse_arguments())
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(f"error: {error}") from error
