#!/usr/bin/env python3
"""Plot a transverse TDDFT dispersion as a q--omega loss imagemap.

The script reads the self-describing production files written by the TDDFT
validation workflow.  It works with both the VAL-18 Fe and VAL-19 Ni response
directories and does not require conversion to an intermediate format.

Examples
--------
Plot the legacy interacting loss for Fe and overlay accepted Lorentzian fits::

    python3 tools/plot_tddft_dispersion.py \
        --response-dir results/validation/VAL-18_bccFe/dispersion_nk16 \
        --product legacy --quantity trace_loss --mode-route legacy \
        --output fe_dispersion.png

Plot the Ni legacy loss.  Since VAL-19 has no accepted fit, this produces the
imagemap without a misleading magnon line::

    python3 tools/plot_tddft_dispersion.py \
        --response-dir results/validation/VAL-19_fccNi/response_nk8 \
        --product legacy --quantity trace_loss --mode-route legacy \
        --output ni_dispersion.png

To inspect the formal Xi crossings, rather than accepted mode fits, use
``--mode-source crossing``.  Rejected or negative-weight crossings are omitted
unless ``--include-rejected`` is supplied.  Signed loss is retained in the
color map; in particular, negative VAL-19 spectral weight is not clipped.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np


RY_TO_MEV = 13_605.693_122_994
Q_INDEX_RE = re.compile(r"_q(\d+)_")
FLOAT_WITHOUT_E_RE = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)[+-]\d+")


class PlotError(ValueError):
    """A malformed or scientifically ambiguous TDDFT plotting input."""


def parse_float(token: str) -> float:
    """Parse Python floats and gfortran's overflow spelling (``1.0+308``)."""
    normalized = token.replace("D", "E").replace("d", "e")
    try:
        return float(normalized)
    except ValueError:
        if not FLOAT_WITHOUT_E_RE.fullmatch(normalized):
            raise
        exponent = re.search(r"[+-]\d+$", normalized)
        if exponent is None:
            raise
        return float(normalized[: exponent.start()] + "e" + exponent.group(0))


def q_index(path: Path) -> int:
    match = Q_INDEX_RE.search(path.name)
    if match is None:
        raise PlotError(f"cannot determine q index from {path.name}")
    return int(match.group(1))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--response-dir",
        type=Path,
        required=True,
        help="TDDFT response directory containing *_manifest.dat and per-q products",
    )
    parser.add_argument("--prefix", help="output prefix; inferred when the directory has one manifest")
    parser.add_argument(
        "--product",
        choices=("auto", "legacy", "raw_pair", "corrected_pair", "chi0"),
        default="legacy",
        help="spectrum family to plot (default: legacy; auto prefers corrected pair, raw pair, then legacy)",
    )
    parser.add_argument(
        "--quantity",
        choices=("trace_loss", "site_loss", "trace", "site_diagonal"),
        default="trace_loss",
        help="record to plot; Dyson products use trace_loss/site_loss, chi0 uses trace/site_diagonal",
    )
    parser.add_argument("--site", type=int, default=1, help="one-based site for site_loss/site_diagonal (default: 1)")
    parser.add_argument(
        "--mode-route",
        choices=("auto", "legacy", "pair"),
        default="auto",
        help="mode file route for the optional overlay (default: auto, preferring legacy)",
    )
    parser.add_argument("--modes", type=Path, help="explicit mode file; overrides --mode-route")
    parser.add_argument(
        "--mode-source",
        choices=("fit", "crossing"),
        default="fit",
        help="overlay accepted fits or Xi crossings (default: fit)",
    )
    parser.add_argument(
        "--include-rejected",
        action="store_true",
        help="also show rejected fits/formal crossings; use only for diagnostics",
    )
    parser.add_argument(
        "--x-axis",
        choices=("path", "qnorm"),
        default="path",
        help="horizontal coordinate: accumulated direct-q path or |q| in direct coordinates",
    )
    parser.add_argument("--energy-unit", choices=("meV", "Ry"), default="meV")
    parser.add_argument("--energy-max", type=float, help="discard frequencies above this value in --energy-unit")
    parser.add_argument(
        "--linear-colors",
        action="store_true",
        help="use a linear signed color scale; otherwise use log scale for positive data and symmetric-log for signed data",
    )
    parser.add_argument("--title", help="plot title (default: derived from response directory and product)")
    parser.add_argument("--output", type=Path, default=Path("tddft_dispersion.png"))
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--show", action="store_true", help="show the plot after writing it")
    return parser.parse_args()


def find_manifest(response_dir: Path, prefix: str | None) -> tuple[str, Path]:
    if prefix is not None:
        path = response_dir / f"{prefix}_manifest.dat"
        if not path.exists():
            raise FileNotFoundError(path)
        return prefix, path
    manifests = sorted(response_dir.glob("*_manifest.dat"))
    if len(manifests) != 1:
        names = ", ".join(path.name for path in manifests) or "none"
        raise PlotError(f"expected one manifest in {response_dir}, found: {names}; pass --prefix")
    path = manifests[0]
    return path.name[: -len("_manifest.dat")], path


def read_manifest(path: Path) -> tuple[np.ndarray, list[dict[str, object]]]:
    q_points: list[list[float]] = []
    records: list[dict[str, object]] = []
    for raw_line in path.read_text().splitlines():
        fields = raw_line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) < 7:
            raise PlotError(f"{path}: malformed manifest row: {raw_line}")
        q_points.append([parse_float(value) for value in fields[1:4]])
        records.append(
            {
                "index": int(fields[0]),
                "chi0": fields[4],
                "legacy": fields[5],
                "raw_pair": fields[6],
            }
        )
    if not records:
        raise PlotError(f"{path}: no q-point records")
    return np.asarray(q_points, dtype=float), records


def product_paths(response_dir: Path, prefix: str, records: list[dict[str, object]], product: str) -> list[Path]:
    if product == "corrected_pair":
        paths = sorted(response_dir.glob(f"{prefix}_q*_pair_corrected_dyson.dat"), key=q_index)
    else:
        key = {"legacy": "legacy", "raw_pair": "raw_pair", "chi0": "chi0"}[product]
        paths = [response_dir / str(record[key]) for record in records]
    if not paths or not all(path.exists() for path in paths):
        missing = [str(path) for path in paths if not path.exists()]
        detail = ", ".join(missing[:3]) or "no matching files"
        raise FileNotFoundError(f"{product} products unavailable in {response_dir}: {detail}")
    return paths


def choose_product(response_dir: Path, prefix: str, records: list[dict[str, object]], requested: str) -> tuple[str, list[Path]]:
    candidates = ("corrected_pair", "raw_pair", "legacy", "chi0") if requested == "auto" else (requested,)
    errors: list[str] = []
    for product in candidates:
        try:
            return product, product_paths(response_dir, prefix, records, product)
        except FileNotFoundError as error:
            errors.append(str(error))
    raise FileNotFoundError("; ".join(errors))


def read_spectrum(path: Path, product: str, quantity: str, site: int) -> tuple[np.ndarray, np.ndarray]:
    omega: list[float] = []
    value: list[float] = []
    is_chi0 = product == "chi0"
    if is_chi0 and quantity not in {"trace", "site_diagonal"}:
        raise PlotError("chi0 products require --quantity trace or site_diagonal")
    if not is_chi0 and quantity not in {"trace_loss", "site_loss"}:
        raise PlotError("Dyson products require --quantity trace_loss or site_loss")

    for raw_line in path.read_text().splitlines():
        fields = raw_line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if is_chi0:
            if len(fields) < 7 or fields[1] != quantity:
                continue
            if quantity == "site_diagonal" and int(fields[2]) != site:
                continue
            omega.append(parse_float(fields[0]))
            value.append(parse_float(fields[6]))
        elif quantity == "trace_loss":
            if len(fields) >= 3 and fields[0] == "trace_loss":
                omega.append(parse_float(fields[1]))
                value.append(parse_float(fields[2]))
        elif len(fields) >= 4 and fields[0] == "site_loss" and int(fields[2]) == site:
            omega.append(parse_float(fields[1]))
            value.append(parse_float(fields[3]))

    if not omega:
        raise PlotError(f"{path}: no {quantity} records{f' for site {site}' if 'site' in quantity else ''}")
    order = np.argsort(omega)
    return np.asarray(omega, dtype=float)[order], np.asarray(value, dtype=float)[order]


def collect_spectrum(
    paths: list[Path], q_points: np.ndarray, product: str, quantity: str, site: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reference_omega: np.ndarray | None = None
    spectra: list[np.ndarray] = []
    for path in paths:
        omega, values = read_spectrum(path, product, quantity, site)
        if reference_omega is None:
            reference_omega = omega
        elif omega.shape != reference_omega.shape or not np.allclose(omega, reference_omega, rtol=0.0, atol=1.0e-12):
            raise PlotError(f"{path}: frequency mesh differs from earlier q points")
        spectra.append(values)
    assert reference_omega is not None
    if len(q_points) != len(spectra):
        raise PlotError(f"manifest has {len(q_points)} q points but products contain {len(spectra)} spectra")
    return q_points, reference_omega, np.asarray(spectra, dtype=float)


def x_coordinates(q_points: np.ndarray, axis: str) -> np.ndarray:
    if axis == "qnorm":
        return np.linalg.norm(q_points, axis=1)
    if len(q_points) == 1:
        return np.zeros(1)
    steps = np.linalg.norm(np.diff(q_points, axis=0), axis=1)
    return np.concatenate(([0.0], np.cumsum(steps)))


def mode_path(response_dir: Path, prefix: str, route: str, explicit: Path | None) -> Path | None:
    if explicit is not None:
        if explicit.is_absolute():
            path = explicit
        else:
            cwd_path = Path.cwd() / explicit
            path = cwd_path if cwd_path.exists() else response_dir / explicit
        if not path.exists():
            raise FileNotFoundError(path)
        return path
    routes = ("legacy", "pair") if route == "auto" else (route,)
    for candidate in routes:
        path = response_dir / f"{prefix}_{candidate}_modes.dat"
        if path.exists():
            return path
    return None


def read_modes(path: Path, source: str, include_rejected: bool) -> dict[int, float]:
    modes: dict[int, float] = {}
    for raw_line in path.read_text().splitlines():
        fields = raw_line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if source == "fit" and fields[0] == "fit" and len(fields) >= 7:
            accepted = fields[2].upper().startswith("T")
            center = parse_float(fields[3])
            if (accepted or include_rejected) and math.isfinite(center) and center > 0.0:
                modes[int(fields[1])] = center
        elif source == "crossing" and fields[0] == "crossing" and len(fields) >= 10:
            present = fields[2].upper().startswith("T")
            omega = parse_float(fields[3])
            weight = parse_float(fields[7])
            usable = present and math.isfinite(omega) and omega >= 0.0 and weight > 0.0
            if include_rejected:
                usable = present and math.isfinite(omega) and omega >= 0.0
            if usable:
                modes[int(fields[1])] = omega
    return modes


def make_plot(args: argparse.Namespace) -> None:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm, Normalize, SymLogNorm
    except ModuleNotFoundError as error:
        raise SystemExit("plotting requires matplotlib (and numpy); install matplotlib in the active Python environment") from error

    response_dir = args.response_dir.resolve()
    if not response_dir.is_dir():
        raise FileNotFoundError(response_dir)
    prefix, manifest = find_manifest(response_dir, args.prefix)
    q_points, records = read_manifest(manifest)
    product, paths = choose_product(response_dir, prefix, records, args.product)
    q_points, omega_ry, values = collect_spectrum(paths, q_points, product, args.quantity, args.site)

    energy = omega_ry if args.energy_unit == "Ry" else omega_ry * RY_TO_MEV
    if args.energy_max is not None:
        keep = energy <= args.energy_max
        if not np.any(keep):
            raise PlotError("--energy-max excludes the entire frequency grid")
        energy = energy[keep]
        values = values[:, keep]

    x = x_coordinates(q_points, args.x_axis)
    finite = np.isfinite(values)
    if not np.all(finite):
        raise PlotError(f"spectrum contains {int(np.count_nonzero(~finite))} non-finite values")
    if np.allclose(values, 0.0):
        raise PlotError("spectrum is zero everywhere")

    fig, axis = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)
    minimum = float(np.min(values))
    maximum = float(np.max(values))
    signed = minimum < 0.0
    if signed:
        scale = max(abs(minimum), abs(maximum))
        if args.linear_colors:
            norm = Normalize(vmin=-scale, vmax=scale)
        else:
            norm = SymLogNorm(linthresh=max(scale * 1.0e-3, np.finfo(float).tiny), vmin=-scale, vmax=scale)
        cmap = "RdBu_r"
    else:
        positive = values[values > 0.0]
        if args.linear_colors:
            norm = Normalize(vmin=0.0, vmax=maximum)
        else:
            norm = LogNorm(vmin=max(float(np.percentile(positive, 1.0)), np.finfo(float).tiny), vmax=maximum)
        cmap = "magma"

    image = axis.pcolormesh(x, energy, values.T, shading="auto", cmap=cmap, norm=norm)
    color_label = "signed spectral weight" if signed else "spectral weight"
    fig.colorbar(image, ax=axis, label=color_label)

    mode_route = args.mode_route
    if mode_route == "auto" and product in {"raw_pair", "corrected_pair"}:
        mode_route = "pair"
    elif mode_route == "auto" and product == "legacy":
        mode_route = "legacy"
    mode_file = mode_path(response_dir, prefix, mode_route, args.modes)
    if mode_file is not None:
        modes = read_modes(mode_file, args.mode_source, args.include_rejected)
        indices = sorted(index for index in modes if 1 <= index <= len(x))
        if indices:
            mode_x = np.asarray([x[index - 1] for index in indices])
            mode_omega = np.asarray([modes[index] for index in indices])
            if args.energy_unit == "meV":
                mode_omega *= RY_TO_MEV
            label = f"{mode_file.stem}: {args.mode_source}"
            if args.include_rejected:
                label += " (incl. rejected)"
            axis.plot(mode_x, mode_omega, "o-", color="cyan", lw=1.5, ms=4.0, label=label)
            axis.legend(loc="best")
        else:
            print(f"warning: {mode_file.name} has no overlayable {args.mode_source} records")
    elif args.mode_route != "auto" or args.modes is not None:
        print("warning: requested mode overlay file was not found")

    axis.set_xlabel("|q| in direct reciprocal coordinates" if args.x_axis == "qnorm" else "q-path coordinate in direct reciprocal coordinates")
    axis.set_ylabel(f"energy ({args.energy_unit})")
    axis.set_title(args.title or f"{response_dir.name}: {product} {args.quantity}")
    axis.set_xticks(x)
    axis.set_xticklabels([f"{point[0]:.3g},{point[1]:.3g},{point[2]:.3g}" for point in q_points], rotation=35, ha="right")
    axis.grid(axis="y", alpha=0.18)
    if signed:
        print(f"note: preserving signed spectral weight, range [{minimum:.6e}, {maximum:.6e}]")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi)
    print(f"wrote {args.output}")
    if mode_file is not None:
        print(f"mode source: {mode_file}")
    if args.show:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    args = parse_args()
    try:
        make_plot(args)
    except (FileNotFoundError, OSError, PlotError) as error:
        raise SystemExit(f"error: {error}") from error


if __name__ == "__main__":
    main()
