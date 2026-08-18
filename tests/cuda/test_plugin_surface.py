#!/usr/bin/env python3
"""Static coverage contract for the existing RS CUDA plugin.

This test deliberately does not load CUDA.  It catches interface drift between
the Fortran wrapper, the public C ABI, the C++ dispatcher, and the production
CPU/GPU route call sites.  Numerical equivalence remains a hardware test.
"""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FORTRAN = ROOT / "source/rsrec_cuda_plugin.f90"
HEADER = ROOT / "source/cuda/rsrec_cuda.h"
DISPATCHER = ROOT / "source/cuda/rsrec_cuda.cpp"
CUDA = ROOT / "source/cuda/rsrec_gpu.cu"
CORE = ROOT / "source/recursion_core.f90"


ROUTES = {
    "scalar Lanczos": (
        "recur()",
        "scalar_lanczos",
        "rsrec_cuda_scalar_lanczos",
        "rsrec_scalar_lanczos",
    ),
    "Block Lanczos": (
        "recur_b()",
        "block_lanczos",
        "rsrec_cuda_block_lanczos",
        "rsrec_block_lanczos",
    ),
    "Block intersite Lanczos": (
        "recur_b_ij()",
        "block_lanczos",
        "rsrec_cuda_block_lanczos",
        "rsrec_block_lanczos",
    ),
    "Chebyshev moments": (
        "chebyshev_recur()",
        "chebyshev_moments",
        "rsrec_cuda_chebyshev_moments",
        "rsrec_chebyshev_moments",
    ),
    "Chebyshev intersite moments": (
        "chebyshev_recur_ij()",
        "chebyshev_moments",
        "rsrec_cuda_chebyshev_moments",
        "rsrec_chebyshev_moments",
    ),
    "stochastic moments": (
        "compute_moments_stochastic()",
        "stochastic_moments",
        "rsrec_cuda_stochastic_moments",
        "rsrec_stochastic_moments",
    ),
    "orbital moments": (
        "chebyshev_orbital_mod()",
        "orbital_moments",
        "rsrec_cuda_orbital_moments",
        "rsrec_orbital_moments",
    ),
    "Chebyshev DOS": (
        "chebyshev_green_core",
        "chebyshev_dos",
        "rsrec_cuda_chebyshev_dos",
        "rsrec_chebyshev_dos",
    ),
    "Chebyshev eta GF": (
        "chebyshev_green_core",
        "chebyshev_gf_eta",
        "rsrec_cuda_chebyshev_gf_eta",
        "rsrec_chebyshev_gf_eta",
    ),
    "Block DOS": (
        "block_green_core",
        "block_dos",
        "rsrec_cuda_block_dos",
        "rsrec_block_dos",
    ),
    "Block eta GF": (
        "calculate_intersite_gf_eta_gpu",
        "block_gf_eta",
        "rsrec_cuda_block_gf_eta",
        "rsrec_block_gf_eta",
    ),
}


def functions(text: str, prefix: str) -> set[str]:
    return set(re.findall(rf"\b({re.escape(prefix)}[A-Za-z0-9_]*)\s*\(", text))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> int:
    for path in (FORTRAN, HEADER, DISPATCHER, CUDA, CORE):
        require(path.is_file(), f"missing plugin source: {path}")

    fortran = FORTRAN.read_text()
    header = HEADER.read_text()
    dispatcher = DISPATCHER.read_text()
    cuda = CUDA.read_text()
    core = CORE.read_text()
    production_sources = "\n".join(
        path.read_text()
        for path in (ROOT / "source").glob("*.f90")
    )

    bindings = set(re.findall(r"bind\(C,\s*name='([^']+)'\)", fortran))
    header_api = functions(header, "rsrec_cuda_")
    dispatcher_api = functions(dispatcher, "rsrec_cuda_")
    cuda_api = functions(cuda, "rsrec_")

    require(bindings <= header_api, "Fortran C bindings missing from rsrec_cuda.h: "
            + repr(sorted(bindings - header_api)))
    require(bindings <= dispatcher_api, "Fortran C bindings missing from C++ dispatcher: "
            + repr(sorted(bindings - dispatcher_api)))

    for name, (cpu_route, method, cuda_name, cpu_name) in ROUTES.items():
        require(cpu_route in production_sources,
                f"{name}: CPU route marker {cpu_route!r} is missing")
        require(f"procedure :: {method}" in fortran or f"%{method}" in production_sources,
                f"{name}: Fortran method {method!r} is missing")
        require(cuda_name in dispatcher_api, f"{name}: C++ dispatcher {cuda_name!r} is missing")
        require(cpu_name in cuda_api, f"{name}: CUDA implementation {cpu_name!r} is missing")

    require("does not support orbital moments" in core,
            "structured orbital-moment rejection must remain explicit")
    require("allow_hoh=.true." in (ROOT / "source/recursion_haydock.f90").read_text(),
            "Block Lanczos HOH capability guard is missing")
    require("allow_hoh=.true." in (ROOT / "source/recursion_transport.f90").read_text(),
            "transport HOH capability guard is missing")

    print(f"CUDA plugin surface: {len(bindings)} Fortran C bindings, "
          f"{len(ROUTES)} production/reconstruction routes: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
