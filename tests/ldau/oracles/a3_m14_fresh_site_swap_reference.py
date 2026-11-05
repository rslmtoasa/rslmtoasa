#!/usr/bin/env python3
"""Generate a fresh independent finite-H site/type mapping control fixture.

This is intentionally separate from the already measured M1.4 fixture. It
uses a new deterministic seed, computes the target state from exact eigenpairs,
then predicts the wrong same-state result obtained by substituting the other
site's J=0 V while holding target P/M1 and radial moments fixed.
"""
from __future__ import annotations

import argparse
import json
import hashlib
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import a3_m14_independent_reference as reference  # noqa: E402

SEED = 20261016
EPS = np.finfo(np.float64).eps


def gamma(n: int) -> float:
    return n * EPS/(1.0-n*EPS)


def trace(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.trace(a @ b).real)


def build() -> dict:
    rng = np.random.default_rng(SEED)
    cases = [reference.one_case(rng, site, spin) for site in (1, 2) for spin in (1, 2)]
    by_key = {(x['site'], x['spin']): x for x in cases}
    records = []
    for case in cases:
        other = by_key[(3-case['site'], case['spin'])]
        p = np.asarray(case['P_re']) + 1j*np.asarray(case['P_im'])
        m1 = np.asarray(case['M1_re']) + 1j*np.asarray(case['M1_im'])
        vwrong = np.asarray(other['Vin_re']) + 1j*np.asarray(other['Vin_im'])
        q0, c, q2 = case['q0'], case['c_ry'], case['q2_ry2']
        t_wrong = trace(vwrong, p)
        c_wrong = c-t_wrong/q0
        m1c = m1-c*p
        q2_wrong = q2-2.0*trace(vwrong,m1c)+trace(vwrong@vwrong,p)-t_wrong*t_wrong/q0

        h0 = np.asarray(case['H0_re'])+1j*np.asarray(case['H0_im'])
        u = np.asarray(case['eigenvectors_re'])+1j*np.asarray(case['eigenvectors_im'])
        e = np.asarray(case['eigenvalues_ry'])
        f = np.asarray(case['occupations'])
        a = u[:5,:]
        # This is a held-P mapping control, not a second finite-H solve:
        # target P/M1 are identical in both calculations. Thus finite-H
        # residual uncertainty is already common-mode and cancels; bound the
        # two independently evaluated trace expressions by their arithmetic
        # forward errors only, using the same gamma operation-count model.
        wrong_c_scale=(abs(c)+abs(t_wrong/q0)+abs(case['cprime_direct_ry'])
            +float(np.sum(f[None,:]*np.abs(a)*(np.abs(h0)@np.abs(u))[:5,:])/q0)
            +float(np.sum(np.abs(vwrong)*np.abs(p.T))/q0)
            +float(np.sum(np.abs(m1))/q0))
        wrong_c_bound=gamma(4096)*wrong_c_scale
        wrong_q2_scale=(abs(q2)
            +2*float(np.sum(np.abs(vwrong)*(np.abs(m1.T)+abs(c)*np.abs(p.T))))
            +float(np.sum((np.abs(vwrong)@np.abs(vwrong))*np.abs(p.T)))
            +abs(t_wrong*t_wrong/q0)
            +float(np.sum(np.abs(vwrong)*np.abs(p.T)))*(1+abs(t_wrong/q0)))
        wrong_q2_bound=gamma(16384)*wrong_q2_scale
        c_sep=abs(c_wrong-case['cprime_direct_ry'])
        q2_sep=abs(q2_wrong-case['q2prime_direct_ry2'])
        target_c_bound=case['c_error_bound_ry']; target_q2_bound=case['q2_error_bound_ry2']
        records.append({
            'site':case['site'],'spin':case['spin'],'q0':q0,'c_ry':c,'q2_ry2':q2,
            'cprime_direct_ry':case['cprime_direct_ry'],'q2prime_direct_ry2':case['q2prime_direct_ry2'],
            'c_error_bound_ry':target_c_bound,'q2_error_bound_ry2':target_q2_bound,
            'cprime_wrong_v_ry':c_wrong,'q2prime_wrong_v_ry2':q2_wrong,'t_wrong_ry':t_wrong,
            'wrong_v_source_site':other['site'],'center_separation_ry':c_sep,
            'q2_separation_ry2':q2_sep,'wrong_center_error_bound_ry':wrong_c_bound,
            'wrong_q2_error_bound_ry2':wrong_q2_bound,
            'center_combined_bound_ry':target_c_bound+wrong_c_bound,
            'q2_combined_bound_ry2':target_q2_bound+wrong_q2_bound,
            'center_eligible':bool(c_sep>target_c_bound+wrong_c_bound),
            'q2_eligible':bool(q2_sep>target_q2_bound+wrong_q2_bound),
            'P_re':case['P_re'],'P_im':case['P_im'],'M1_re':case['M1_re'],'M1_im':case['M1_im'],
            'Vin_re':case['Vin_re'],'Vin_im':case['Vin_im'],
            'Vwrong_re':other['Vin_re'],'Vwrong_im':other['Vin_im'],
        })
    return {'schema':'rslmto-a3-m14-fresh-held-p-site-swap-v1','seed':SEED,
        'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_oracle_sha256':hashlib.sha256((HERE/'a3_m14_independent_reference.py').read_bytes()).hexdigest(),
        'method':'Fresh finite-H cases use a new seed. Correct target comes from direct H0 eigenpair sums; wrong prediction holds target P/M1/q0/c/q2 fixed and substitutes the other site same-spin V. Eligibility is separation greater than correct plus wrong error bounds.',
        'error_model':'Correct prediction uses the independent finite-H reference bound (gamma operation count plus explicit eigensolver-residual propagation). The held-P alternate-V control shares exactly the same P/M1 inputs, so finite-H residual terms are common-mode; its independent trace-formula forward error is bounded by gamma_4096/16384 times the recomputed absolute intermediate scales for alternate V. Combined admissibility bound is the sum of correct and alternate-V bounds. Units are Ry/Ry^2.',
        'records':records,'all_controls_eligible':all(r['center_eligible'] or r['q2_eligible'] for r in records)}


def write_fortran(obj: dict, path: Path) -> None:
    out=['A3M14SWAP 1',str(len(obj['records']))]
    for r in obj['records']:
        fields=['site','spin','q0','c_ry','q2_ry2','cprime_direct_ry','q2prime_direct_ry2',
                'c_error_bound_ry','q2_error_bound_ry2','cprime_wrong_v_ry','q2prime_wrong_v_ry2',
                'center_combined_bound_ry','q2_combined_bound_ry2']
        out.append(' '.join(f"{r[k]:.17g}" for k in fields))
        for name in ('P','M1','Vin','Vwrong'):
            mat=np.asarray(r[f'{name}_re'])+1j*np.asarray(r[f'{name}_im'])
            for j in range(5):
                for i in range(5): out.append(f'{mat[i,j].real:.17g} {mat[i,j].imag:.17g}')
    path.write_text('\n'.join(out)+'\n')


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument('--json',type=Path,required=True)
    parser.add_argument('--fortran',type=Path,required=True)
    args=parser.parse_args()
    obj=build(); args.json.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
    write_fortran(obj,args.fortran)


if __name__=='__main__': main()
