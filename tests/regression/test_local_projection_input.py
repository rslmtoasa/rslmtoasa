"""Exercise namelist validation through the normal material entry point."""
import os
from pathlib import Path
import shutil
import subprocess

import f90nml
import pytest

EXE = Path(os.environ.get("RSLMTO_EXE", "build/bin/rslmto.x")).resolve()
FIXTURE = Path(__file__).resolve().parents[1]/"postproc/cases/conductivity/fccPt"


@pytest.mark.skipif(not EXE.exists(), reason="set RSLMTO_EXE to run input tests")
@pytest.mark.parametrize("override,diagnostic", [
    ({"response_projection":"sideways"}, "response_projection must be"),
    ({"cond_calctype":"random_vec"}, "per_type and charge"),
    ({"recur":"block"}, "requires Chebyshev"),
    ({"constraints_enable":True,"constraints_i_cons":2}, "restore SCF constraining fields"),
    ({"linear_in":"spin_torque_onsite"}, "output operator"),
    ({"linear_in":"spin"}, "per_type and charge"),
    ({"linear_out":"spin_torque"}, "onsite observable"),
    ({"linear_out":"orbital"}, "onsite observable"),
])
def test_invalid_request(tmp_path, override, diagnostic):
    nml = f90nml.read(FIXTURE/"input.nml")
    nml["control"].update(response_projection="local", linear_out="spin_torque_onsite")
    if "constraints_enable" in override:
        nml["constraints"] = override
    else:
        nml["control"].update(override)
    nml.write(tmp_path/"input.nml")
    shutil.copy(FIXTURE/"Pt.nml", tmp_path)
    run = subprocess.run([EXE], cwd=tmp_path, capture_output=True, text=True,
                         env={**os.environ,"OMP_NUM_THREADS":"1"})
    assert run.returncode != 0
    assert diagnostic in run.stdout+run.stderr

@pytest.mark.skipif(not EXE.exists(), reason="set RSLMTO_EXE to run input tests")
@pytest.mark.parametrize("post,route,message", [
    ("conductivity","lehmann","requires gf_route=recursion"),
    ("conductivity","dyson","requires gf_route=recursion"),
    ("conductivity_p2rs","recursion","native transport Hamiltonian"),
    ("conductivity_import","recursion","native transport Hamiltonian"),
    ("conductivity_effective","recursion","native transport Hamiltonian"),
])
def test_unsupported_source(tmp_path, post, route, message):
    nml = f90nml.read(FIXTURE/"input.nml")
    nml["calculation"].update(post_processing=post, gf_route=route)
    nml["control"].update(response_projection="local", linear_out="spin_torque_onsite")
    nml.write(tmp_path/"input.nml")
    run = subprocess.run([EXE], cwd=tmp_path, capture_output=True, text=True)
    assert run.returncode != 0
    assert message in run.stdout+run.stderr
