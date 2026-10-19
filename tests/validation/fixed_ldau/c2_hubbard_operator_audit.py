#!/usr/bin/env python3
"""Compare opt-in C2 records for the first bulk and local Hubbard insertions.

Set RSLMTO_C2_HUBBARD_AUDIT=1 and optionally
RSLMTO_C2_HUBBARD_AUDIT_FILE=<path> when running a one-step D1 deck. The
production run also needs the normal full C1 trace to record P_OUT matrices.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def read_records(path: Path) -> dict[tuple[str, int, int], dict[str, object]]:
    records: dict[tuple[str, int, int], dict[str, object]] = {}
    for line_number, raw in enumerate(path.read_text().splitlines(), 1):
        fields = raw.split()
        if not fields or fields[0] not in {"META", "OCC", "BLOCK"}:
            continue
        try:
            tag, assembler = fields[0], fields[1]
            build_id, site, ptype = map(int, fields[2:5])
            key = (assembler, build_id, site)
            record = records.setdefault(key, {"ptype": ptype, "occ": {}, "block": {}})
            if int(record["ptype"]) != ptype:
                raise ValueError("potential type changes within one site/build record")
            if tag == "META":
                if len(fields) != 9:
                    raise ValueError(f"META expected 9 fields, got {len(fields)}")
                record.update(symbol=fields[5], u_ry=float(fields[6]), j_ry=float(fields[7]), count=int(fields[8]))
            elif tag == "OCC":
                if len(fields) != 10:
                    raise ValueError(f"OCC expected 10 fields, got {len(fields)}")
                spin, i, j = map(int, fields[5:8])
                record["occ"][spin, i, j] = complex(float(fields[8]), float(fields[9]))
            else:
                if len(fields) != 15:
                    raise ValueError(f"BLOCK expected 15 fields, got {len(fields)}")
                spin, i, j, count = map(int, fields[5:9])
                record["block"][spin, i, j] = {
                    "count": count,
                    "before": complex(float(fields[9]), float(fields[10])),
                    "correction": complex(float(fields[11]), float(fields[12])),
                    "after": complex(float(fields[13]), float(fields[14])),
                }
        except (IndexError, ValueError) as exc:
            raise ValueError(f"{path}:{line_number}: malformed C2 record: {exc}") from exc
    return records


def choose_record(records, assembler: str, site: int, build: int | None):
    matches = [
        value for (kind, build_id, site_id), value in records.items()
        if kind == assembler and site_id == site and (build is None or build_id == build)
    ]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {assembler} record for site={site}, build={build}; found {len(matches)}")
    record = matches[0]
    if not all(key in record for key in ("u_ry", "j_ry", "count", "occ", "block")):
        raise ValueError(f"incomplete {assembler} record for site={site}")
    if len(record["occ"]) != 50 or len(record["block"]) != 50:
        raise ValueError(f"{assembler} site={site} must contain 50 OCC and 50 BLOCK entries")
    return record


def compare_records(bulk, local, tolerance: float) -> dict[str, float | int]:
    def max_difference(a, b):
        if set(a) != set(b):
            raise ValueError("bulk and local records have different matrix-entry keys")
        return max(abs(a[key] - b[key]) for key in a)

    u_diff = abs(float(bulk["u_ry"]) - float(local["u_ry"]))
    j_diff = abs(float(bulk["j_ry"]) - float(local["j_ry"]))
    occ_diff = max_difference(bulk["occ"], local["occ"])
    correction_diff = max_difference(
        {key: value["correction"] for key, value in bulk["block"].items()},
        {key: value["correction"] for key, value in local["block"].items()},
    )
    symbols_match = bulk["symbol"] == local["symbol"]
    insertion_counts_match = all(
        value["count"] == 1 for record in (bulk, local) for value in record["block"].values()
    ) and int(bulk["count"]) == 1 and int(local["count"]) == 1
    arithmetic_error = max(
        abs(value["after"] - value["before"] - value["correction"])
        for record in (bulk, local) for value in record["block"].values()
    )
    result = {
        "u_diff_ry": u_diff,
        "j_diff_ry": j_diff,
        "occupation_max_abs_diff": occ_diff,
        "correction_max_abs_diff": correction_diff,
        "symbols_match": int(symbols_match),
        "insertion_counts_match": int(insertion_counts_match),
        "post_minus_pre_max_abs_error": arithmetic_error,
    }
    if max(u_diff, j_diff, occ_diff, correction_diff, arithmetic_error) > tolerance or \
            not symbols_match or not insertion_counts_match:
        raise AssertionError(f"bulk/local C2 records differ beyond {tolerance:g}: {result}")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit_file", type=Path)
    parser.add_argument("--bulk-site", type=int, required=True, help="site_index from BULK_EE META record")
    parser.add_argument("--local-site", type=int, required=True, help="site_index from LOCAL_HALL META record")
    parser.add_argument("--bulk-build", type=int)
    parser.add_argument("--local-build", type=int)
    parser.add_argument("--tolerance", type=float, default=1.0e-12)
    args = parser.parse_args(argv)
    records = read_records(args.audit_file)
    bulk = choose_record(records, "BULK_EE", args.bulk_site, args.bulk_build)
    local = choose_record(records, "LOCAL_HALL", args.local_site, args.local_build)
    result = compare_records(bulk, local, args.tolerance)
    for name, value in result.items():
        print(f"{name}={value}")
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError) as exc:
        print(f"RESULT: FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
