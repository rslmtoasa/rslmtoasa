#!/usr/bin/env python3
"""Compare full-grid d-orbital DOS outputs without interpolation or tolerance."""
from __future__ import annotations
import argparse, hashlib, json
from decimal import Decimal
from pathlib import Path

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def read(path: Path):
    rows=[]
    for line_no, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        vals=line.split()
        if len(vals) != 19:
            raise ValueError(f"{path}:{line_no}: expected energy + 18 orbital columns, got {len(vals)}")
        try: row=[Decimal(v.replace("D","E").replace("d","e")) for v in vals]
        except Exception as e: raise ValueError(f"{path}:{line_no}: invalid numeric value") from e
        if not all(v.is_finite() for v in row):
            raise ValueError(f"{path}:{line_no}: nonfinite value")
        rows.append(row)
    if not rows: raise ValueError(f"{path}: no numeric rows")
    return rows

def analyze(base: Path, cand: Path):
    a,b=read(base),read(cand)
    grid_a=[r[0] for r in a]; grid_b=[r[0] for r in b]
    if len(a)!=len(b) or grid_a != grid_b:
        report={"status":"NOT_COMPUTED_DIFFERENT_ENERGY_GRIDS",
                "definition":"At each existing orbital-DOS grid point in [-1,0.5] Ry relative to that run's EF, sum the five l=2 orbitals separately for each spin and take the maximum absolute candidate-minus-baseline difference. No interpolation or acceptance tolerance.",
                "baseline":{"path":str(base),"sha256":sha(base),"rows":len(a),"energy_relative_to_EF_range_Ry":[str(min(grid_a)),str(max(grid_a))]},
                "candidate":{"path":str(cand),"sha256":sha(cand),"rows":len(b),"energy_relative_to_EF_range_Ry":[str(min(grid_b)),str(max(grid_b))]},
                "energy_grid":{"equal":False,"interpolation":"none","comparison":"not computed","precision":"F16.5 emitted energy values; no matching energies are fabricated"},
                "first_row_difference":None}
        for i in range(min(len(a),len(b))):
            if grid_a[i] != grid_b[i]:
                report["first_row_difference"]={"row_index_zero_based":i,"baseline_energy_relative_to_EF_Ry":str(grid_a[i]),"candidate_energy_relative_to_EF_Ry":str(grid_b[i]),"candidate_minus_baseline_Ry":str(grid_b[i]-grid_a[i])}
                break
        return report
    # File orbital columns are one-based 1..18 after the energy value;
    # l=2 occupies 5..9 up and 14..18 down (bands.f90:492-496,587-588).
    output={}
    for spin, cols in (("up",range(5,10)),("down",range(14,19))):
        diffs=[]
        for i,(ra,rb) in enumerate(zip(a,b)):
            energy=ra[0]
            if Decimal("-1") <= energy <= Decimal("0.5"):
                da=sum((rb[c]-ra[c] for c in cols),Decimal(0))
                diffs.append((abs(da), energy, da, i))
        if not diffs: raise ValueError("no points inside registered energy window [-1,0.5] Ry")
        maxrow=max(diffs,key=lambda x:x[0])
        output[spin]={"max_abs_sum_d_DOS_difference_states_per_Ry":str(maxrow[0]),
                      "energy_relative_to_EF_Ry":str(maxrow[1]),
                      "signed_candidate_minus_baseline":str(maxrow[2]),
                      "row_index_zero_based":maxrow[3],
                      "window_points":len(diffs)}
    return {"status":"DESCRIPTIVE_RECORD_ONLY",
            "definition":"At each existing orbital-DOS grid point in [-1,0.5] Ry relative to that run's EF, sum the five l=2 orbitals separately by spin and take the maximum absolute candidate-minus-baseline difference. No interpolation or acceptance tolerance.",
            "baseline":{"path":str(base),"sha256":sha(base),"rows":len(a)},
            "candidate":{"path":str(cand),"sha256":sha(cand),"rows":len(b)},
            "energy_grid":{"equal":True,"comparison":"exact Decimal equality","precision":"F16.5 emitted ene(i)-fermi; cross-run differences below output resolution are not resolved"},
            "source_mapping":{"path":"source/bands.f90","lines":"492-496, 587-588","columns_one_based_after_energy":{"up_d":[5,9],"down_d":[14,18]},"basis_order":"s(1), p(2:4), d(5:9), spin-down offset norb=9"},
            "results":output}

def main():
    p=argparse.ArgumentParser();p.add_argument("baseline",type=Path);p.add_argument("candidate",type=Path);p.add_argument("--output",type=Path,required=True);a=p.parse_args()
    result=analyze(a.baseline,a.candidate)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":result["status"],"output":str(a.output),"output_sha256":sha(a.output)},sort_keys=True))
if __name__=="__main__": main()
