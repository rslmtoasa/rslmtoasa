"""Compiled tests of output units, trace normalization and bounded quadrature."""
import os
from pathlib import Path
import subprocess
import numpy as np
import pytest
from test_hoh_transport import transport_executable

E=1.602176634e-19
HBAR=1.054571817e-34
RY=2.179874099e-18
OPS=[('charge',1,1,0),('spin',1,1,1),('orbital',1,1,1),
     ('quadrupole',1,1,2),('octupole',1,1,3),
     ('spin_accumulation',0,0,1),('orbital_accumulation',0,0,1),
     ('quadrupole_accumulation',0,0,2),('octupole_accumulation',0,0,3),
     ('oap',1,1,2),('oap_accumulation',0,0,2),('octupole_current',1,1,3),
     ('soc_spin_torque',1,0,1),('spin_torque',1,0,1),('orbital_torque',1,0,1),('spin_soc_torque',1,0,1)]


def legacy_integral(y,x):
    weights=np.zeros(len(x))
    for i in range(1,len(x)-1,2):
        weights[i-1:i+2]+=[1,4,1]
    return np.array([np.sum(weights*y*np.where(x<ef,1.,np.where(x==ef,.5,0.)))
                     *(x[1]-x[0])/3 for ef in x])


def run(executable,path,geometry,operator='charge',input_operator='charge',height=2,extra='',open_sample=False,replicas=1):
    path.mkdir(exist_ok=True,parents=True)
    (path/'input.nml').write_text(f"&conductivity\n response_geometry='{geometry}'\n {extra}\n/\n")
    return subprocess.run([str(executable.with_name('response_units_test')),operator,input_operator,str(height),'open' if open_sample else 'periodic',str(replicas) if replicas!=1 else ''],
                          cwd=path,text=True,capture_output=True,
                          env={**os.environ,'OMP_NUM_THREADS':'2'})


@pytest.mark.parametrize('geometry', ['bulk','sheet'])
@pytest.mark.parametrize('operator,pe,pl,rank',OPS)
def test_output_pipeline(tmp_path,transport_executable,geometry,operator,pe,pl,rank):
    result=run(transport_executable,tmp_path,geometry,operator)
    assert result.returncode==0,result.stdout+result.stderr
    label,nums,powers=(tmp_path/'units.txt').read_text().splitlines()
    weights=np.array([float(s) for s in nums.split()]); assert tuple(map(int,powers.split()))==(pe,pl,rank)
    measure={'bulk':16.,'sheet':4.,'per_atom':5.}[geometry]
    np.testing.assert_allclose(weights[:2],np.array([1,1])/measure,rtol=1e-14)
    if pl:
        conversion={'bulk':E**2/HBAR*1e8,'sheet':2*np.pi,'per_atom':E**2/HBAR*1e-20}[geometry]
    else:
        conversion=E*1e-10/(RY if pe==0 else HBAR)*{'bulk':1e30,'sheet':1e20,'per_atom':1}[geometry]
    np.testing.assert_allclose(weights[2],conversion,rtol=2e-14)
    if rank: assert f'hbar^{rank}' in label
    elif geometry=='sheet': assert label=='e^2/h'
    else: assert 'S' in label
    x=np.linspace(-.7,.7,12)
    primitive=lambda t:t+t*t+t**3
    integral=legacy_integral(1+2*x+3*x*x,x)
    # Legacy prefactor and seed average retained; only bounded Simpson panels.
    raw=np.loadtxt(tmp_path/'cond_total.out')
    np.testing.assert_allclose(raw[:,1],2*16/(np.pi*3.4**2)*legacy_integral(1+2*x+3*x*x,x*float(np.float32(2)-np.float32(.3))/1.7),atol=6e-7,rtol=6e-7)
    np.testing.assert_allclose(raw[:,2],0,atol=1e-30)
    physical=np.loadtxt(tmp_path/'response_total.out')
    np.testing.assert_allclose(physical[:,1],2*16/(np.pi*3.4**2)*legacy_integral(1+2*x+3*x*x,x*float(np.float32(2)-np.float32(.3))/1.7)/measure*conversion,rtol=3e-14,atol=1e-30)
    np.testing.assert_allclose(physical[:,0],(2*x+.1)*13.605703976,atol=2e-14)
    resolved=np.loadtxt(tmp_path/'response_orbitals_real.out')
    np.testing.assert_allclose(resolved[:,1:].sum(axis=1),physical[:,1],rtol=2e-14)
    # Each atomic subspace enters once, regardless of cluster populations.
    types=[np.loadtxt(tmp_path/f'Pt{i}_response.out') for i in (1,2)]
    np.testing.assert_allclose(types[0][:,1],1/measure*integral*conversion,rtol=3e-14,atol=1e-30)
    np.testing.assert_allclose(types[1][:,1],3/measure*integral*conversion,rtol=3e-14,atol=1e-30)
    cell=np.loadtxt(tmp_path/'cell_total.out')
    np.testing.assert_allclose(cell[:,1],types[0][:,1]+types[1][:,1],rtol=3e-14)
    legacy=np.loadtxt(tmp_path/'legacy_comparison.txt')
    np.testing.assert_allclose(legacy[:,0],legacy[:,1],rtol=3e-14,atol=1e-14)
    # Compare with independently assembled legacy Simpson/Fermi weights.
    quad=np.loadtxt(tmp_path/'quadrature.txt')
    z=(x+1)**2
    for power in range(3):
        expected=legacy_integral(z**power,x)
        np.testing.assert_allclose(quad[power*12:(power+1)*12,2],expected,atol=2e-14,rtol=2e-14)
    # Actual production gamma: keep Lorentz damping and n=0 half weights.
    actual=np.fromfile(tmp_path/'gamma.bin',complex).reshape((12,4,4),order='F')
    x=x*float(np.float32(2)-np.float32(.3))/1.7
    n=np.arange(4); damping=np.sinh(6*(1-n/4))/np.sinh(6); damping[0]*=.5
    t=np.cos(np.arccos(x[:,None])*n)
    cn=(x[:,None]-1j*n*np.sqrt(1-x[:,None]**2))*np.exp(1j*n*np.arccos(x[:,None]))
    expected=(cn[:,:,None]*t[:,None,:]+cn.conj()[:,None,:]*t[:,:,None])
    expected*=damping[None,:,None]*damping[None,None,:]/(1-x[:,None,None]**2)**2
    np.testing.assert_allclose(actual,expected,atol=2e-14,rtol=2e-14)


def test_sheet_is_explicit_and_vacuum_independent(tmp_path,transport_executable):
    curves={}
    for geom in ['bulk','sheet']:
        for height in [2,20]:
            path=tmp_path/f'{geom}-{height}'
            r=run(transport_executable,path,geom,height=height)
            assert r.returncode==0,r.stderr
            curves[geom,height]=np.loadtxt(path/'response_total.out')[:,1]
    np.testing.assert_allclose(curves['sheet',2],curves['sheet',20],rtol=1e-14)
    np.testing.assert_allclose(curves['bulk',2],10*curves['bulk',20],rtol=1e-14)
    # A strictly planar lattice (zero third vector) is valid only for sheet units.
    assert run(transport_executable,tmp_path/'planar','sheet',height=0).returncode==0
    assert run(transport_executable,tmp_path/'badplanar','bulk',height=0).returncode!=0


def test_generalized_input_and_explicit_measure(tmp_path,transport_executable):
    r=run(transport_executable,tmp_path,'bulk','octupole_accumulation','spin_accumulation',extra='response_measure=80')
    assert r.returncode==0,r.stdout+r.stderr
    lines=(tmp_path/'units.txt').read_text().splitlines()
    assert 'raw Ry^-2 Angstrom^0' in lines[0]
    np.testing.assert_allclose(np.fromstring(lines[1],sep=' '),[1/80,1/80,1])
    assert 'S/cm' not in (tmp_path/'response_units.txt').read_text()


@pytest.mark.parametrize('geometry,extra',[('auto',''),('bulk','response_measure=-1'),('per_atom','response_measure=2'),('bulk','unknown_input=1')])
def test_reject_bad_units_input(tmp_path,transport_executable,geometry,extra):
    r=run(transport_executable,tmp_path,geometry,extra=extra)
    assert r.returncode!=0


def test_open_measure_and_default_bulk(tmp_path,transport_executable):
    path=tmp_path/'default';path.mkdir()
    (path/'input.nml').write_text('&conductivity /\n')
    r=subprocess.run([str(transport_executable.with_name('response_units_test'))],cwd=path,capture_output=True)
    assert r.returncode==0
    assert (path/'units.txt').read_text().splitlines()[0]=='S/cm'
    assert run(transport_executable,tmp_path/'missing','bulk',open_sample=True).returncode==0
    assert run(transport_executable,tmp_path/'explicit','bulk',open_sample=True,extra='response_measure=16').returncode==0
    assert run(transport_executable,tmp_path/'atom','per_atom',open_sample=True).returncode!=0


def test_cell_sum_independent_of_cluster(tmp_path,transport_executable):
    for geometry in ['bulk','sheet']:
        first=tmp_path/(geometry+'-first')
        second=tmp_path/(geometry+'-second')
        assert run(transport_executable,first,geometry).returncode==0
        assert run(transport_executable,second,geometry,replicas=4,open_sample=True).returncode==0
        for name in ['units.txt','response_total.out','Pt1_response.out','Pt2_response.out']:
            assert (first/name).read_bytes()==(second/name).read_bytes()
