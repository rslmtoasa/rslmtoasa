#!/usr/bin/env python3
"""Read-only Ward decomposition for the accepted 05R Fe gate.

Convention derivation (no empirical factor): sigma+ has its up/down entry 1.
The repository circular bubble is 2*(f_up-f_down)/(E_up-E_down+i*eta).
For D=H_down-H_up, <u|D|d>=(E_down-E_up)<u|d>. A y-axis
rigid rotation has T_up,down=(H_up-H_down)/2=-D/2. Consequently
chi_repo*(-D/2) gives rho_up-rho_down. The D-source coefficient oracle
has prefactor -(f_up-f_down)/(E_up-E_down+i*eta), i.e. -chi_repo/2.
The factor and sign follow from the commutator and the band spectral identity.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
ETA = np.array([.01,.005,.0025,.00125,.000625,.0003125,0.0])
ALGEBRA_TOL = 1.e-10  # double-precision finite-matrix identity; not a material closure tolerance

def coefficient_ward(h, ev, vectors, occupations, weights, etas=ETA):
    n = h.shape[0]//2
    names = ['full_D','half_D','minus_D','diagonal_D','d_sector_D','onsite_diagonal_D']
    onsite_diff = np.einsum('ijk,k->ij',h[n:,n:,:]-h[:n,:n,:],weights/weights.sum())
    sums = np.zeros((len(etas), len(names)))
    target2 = 0.0
    for k, wk in enumerate(weights/weights.sum()):
        u, d = vectors[:n,:,k], vectors[n:,:,k]
        f = occupations[:,k]
        target = (u*f)@u.conj().T-(d*f)@d.conj().T
        target2 += wk*np.linalg.norm(target)**2
        diff = h[n:,n:,k]-h[:n,:n,k]
        only_d = np.zeros_like(diff); only_d[4:9,4:9]=diff[4:9,4:9]
        controls = [diff,diff/2,-diff,np.diag(np.diag(diff)),only_d,np.diag(np.diag(onsite_diff))]
        delta = ev[:,k,None]-ev[None,:,k]
        df = f[:,None]-f[None,:]
        for a, eta in enumerate(etas):
            ratio = np.zeros_like(delta,dtype=complex)
            np.divide(-df,delta+1j*eta,out=ratio,where=(delta+1j*eta)!=0)
            for j, source in enumerate(controls):
                source_band = u.conj().T@source@d
                action = u@(ratio*source_band)@d.conj().T
                sums[a,j] += wk*np.linalg.norm(action-target)**2
    return dict(zip(names, (np.sqrt(sums/target2)).T.tolist()))

def self_test():
    # Noncommuting Hermitian spin blocks, independent occupied density matrices.
    rng = np.random.default_rng(501)
    h = np.zeros((18,18,1),complex); vec = np.zeros_like(h); ev=np.zeros((18,1))
    for spin in range(2):
        a=rng.normal(size=(9,9))+1j*rng.normal(size=(9,9))
        block=(a+a.conj().T)/8+spin*.13*np.eye(9)
        e,v=np.linalg.eigh(block)
        h[9*spin:9*(spin+1),9*spin:9*(spin+1),0]=block
        vec[9*spin:9*(spin+1),9*spin:9*(spin+1),0]=v
        ev[9*spin:9*(spin+1),0]=e
    f=1/(1+np.exp(ev/.07))
    result=coefficient_ward(h,ev,vec,f,np.ones(1),np.array([0.0]))
    assert result['full_D'][0]<ALGEBRA_TOL, result
    assert all(result[name][0]>ALGEBRA_TOL for name in result if name!='full_D'), result
    print('coefficient Ward independent synthetic oracle and five negative controls: PASS')

def read_dump(path):
    with path.open('rb') as stream:
        def read(shape, dtype='<f8'):
            size=int(np.prod(shape)); data=np.fromfile(stream,dtype=dtype,count=size)
            if len(data)!=size: raise RuntimeError('truncated diagnostic export')
            return data.reshape(shape,order='F')
        magic,version,nk,nb,nmat,nr,nc,lmax,npair=read((9,),'<i4')
        assert magic==501 and version==1
        ef,temp,a,b,moment,adapter=read((6,))
        ev=read((nb,nk)); vectors=read((nmat,nb,nk),'<c16'); occ=read((nb,nk))
        weights=read((nk,)); kpoints=read((3,nk))
        h=read((nmat,nmat,nk),'<c16'); torque=read(h.shape,'<c16'); torque_base=read(h.shape,'<c16')
        radius=read((nr,)); w=read((nr,)); sr=read((nr,)); bxc=read((nr,))
        total=read((nr,)); val=read((nr,)); core=read((nr,))
        projected=read((nc,4),'<c16'); branches=read((nr,lmax+1,6,2))
        indices=read((4,nc),'<i4'); modes=read((nr,nc),'<c16')
        dtype=np.dtype([('index','<i4',(3,)),('transition','<c16',(nc,))])
        pairs=np.fromfile(stream,dtype=dtype,count=npair)
        assert len(pairs)==npair and not stream.read(1)
    return locals()

def analyze(path, output):
    z=read_dump(path)
    ev,vec,occ,weights,h=[z[x] for x in ['ev','vectors','occ','weights','h']]
    nr,nc,norb=z['nr'],z['nc'],z['nmat']//2
    k,ib,jb=(z['pairs']['index']-1).T
    up=vec[:norb,ib,k]; down=vec[norb:,jb,k]
    el,er=ev[ib,k],ev[jb,k]
    powers=np.array([np.ones_like(el),el,er,el*er,el*el,er*er])
    coefficients=np.array([np.sum(up[(l*l):((l+1)**2)].conj()*down[(l*l):((l+1)**2)],axis=0)
                           for l in range(z['lmax']+1)])
    point=[]
    for case in range(2):
        point.append(np.einsum('rlb,lp,bp->rp',z['branches'][:,:,:,case],coefficients,powers,optimize=True)/np.sqrt(4*np.pi))
    tc=z['pairs']['transition'].T
    b00=np.sqrt(4*np.pi)*z['bxc']
    sources=[t.T@(z['w']*b00) for t in point]
    source_compact=z['projected'][:,0].conj()@tc
    projection_source_error=np.linalg.norm(source_compact-sources[0])/np.linalg.norm(sources[0])
    tband=[]; tbaseband=[]
    for p in range(len(k)):
        tband.append(np.vdot(vec[:,ib[p],k[p]],z['torque'][:,:,k[p]]@vec[:,jb[p],k[p]]))
        tbaseband.append(np.vdot(vec[:,ib[p],k[p]],z['torque_base'][:,:,k[p]]@vec[:,jb[p],k[p]]))
    tband=np.array(tband);tbaseband=np.array(tbaseband)
    wp=weights[k]/weights.sum()
    norm=lambda x: float(np.sqrt(2*np.sum(wp*np.abs(x)**2)))
    torque_norm=norm(tband); source_norm=norm(sources[0]); torque_relative=norm(tband-sources[0])/torque_norm
    coeff=coefficient_ward(h,ev,vec,occ,weights)
    direct_rotation_error=np.max(np.abs(z['torque'][:norb,norb:,:]+(h[norb:,norb:,:]-h[:norb,:norb,:])/2))
    eigen_error=max(np.max(np.abs(h[:,:,q]@vec[:,:,q]-vec[:,:,q]*ev[:,q])) for q in range(z['nk']))
    rho_moment=sum(weights[q]*np.sum(occ[:,q]*(np.sum(abs(vec[:norb,:,q])**2,axis=0)-np.sum(abs(vec[norb:,:,q])**2,axis=0)))
                   for q in range(z['nk']))/weights.sum()
    def reconstruct(c):
        x=np.zeros(((z['lmax']*2+1)**2,nr),complex)
        # Positive-measure point values only; origin is never assigned a floor.
        for i,(_,l,m,mode) in enumerate(z['indices'].T):
            flat=l*l+l+m
            x[flat,1:]+=z['modes'][1:,i]*c[i]/np.sqrt(z['w'][1:])
        return x
    def point_stats(res,target):
        absnorm=np.sqrt(np.sum(z['w']*abs(res)**2))
        targetnorm=np.sqrt(np.sum(z['w']*abs(target)**2))
        if res.ndim==1: radial_abs=abs(res)
        else: radial_abs=np.max(abs(res),axis=0)
        ix=1+np.argmax(radial_abs[1:])
        return {'absolute':float(absnorm),'relative':float(absnorm/targetnorm),
                'infinity_point_Ylm':float(radial_abs[ix]),'maximum_radius_bohr':float(z['radius'][ix]),'maximum_radial_index':int(ix+1)}
    targets={name:np.sqrt(4*np.pi)*z[key] for name,key in [('total','total'),('valence','val'),('core','core')]}
    projected_targets={'total':z['projected'][:,1],'valence':z['projected'][:,2],'core':z['projected'][:,3]}
    point_norms={name:float(np.sqrt(np.sum(z['w']*abs(v)**2))) for name,v in targets.items()}
    compact_norms={name:float(np.linalg.norm(v)) for name,v in projected_targets.items()}
    core_weight=4*np.pi*z['w']*abs(z['core']); core90=z['radius'][np.searchsorted(np.cumsum(core_weight),.9*core_weight.sum())]
    core_max_index=1+np.argmax(abs(z['core'][1:])); eta_rows=[]; profiles={}
    for eta in ETA:
        factor=np.zeros(len(k),complex)
        np.divide(2*wp*(occ[ib,k]-occ[jb,k]),el-er+1j*eta,out=factor,where=(el-er+1j*eta)!=0)
        compact_response=tc@(factor*sources[0].conj())
        compact_full_source=tc@(factor*tband.conj())
        row={'eta_Ry':float(eta),'compact':{},'point_SR':{},'radial_production':{},'full_H_source':{}}
        for name in ['total','valence']:
            target=projected_targets[name];res=compact_response-target
            cp=point_stats(reconstruct(res),reconstruct(target))
            cp.update(absolute=float(np.linalg.norm(res)),relative=float(np.linalg.norm(res)/np.linalg.norm(target)),
                      infinity_compact=float(np.max(abs(res))),overlap=[float(x) for x in [np.vdot(target,compact_response).real/np.vdot(target,target).real,
                                                                                          np.vdot(target,compact_response).imag/np.vdot(target,target).real]])
            row['compact'][name]=cp
            for label,t,source in [('point_SR',point[0],sources[0]),('radial_production',point[1],sources[1])]:
                response=t@(factor*source.conj())
                row[label][name]=point_stats(response-targets[name],targets[name])
        target=projected_targets['valence'];res=compact_full_source-target
        row['full_H_source']['compact_valence']=point_stats(reconstruct(res),reconstruct(target))
        row['full_H_source']['compact_valence']['relative']=float(np.linalg.norm(res)/np.linalg.norm(target))
        row['full_H_source']['point_SR_valence']=point_stats(point[0]@(factor*tband.conj())-targets['valence'],targets['valence'])
        row['full_H_source']['point_Pauli_valence']=point_stats(point[1]@(factor*tband.conj())-targets['valence'],targets['valence'])
        eta_rows.append(row)
        if eta==0:
            profiles['compact_response']=reconstruct(compact_response)
            profiles['point_SR_response']=point[0]@(factor*sources[0].conj())
            profiles['radial_response']=point[1]@(factor*sources[1].conj())
    # Finite-mesh Fourier R=0 component is the weighted k average; all
    # remaining translations are grouped as hopping (subject to mesh aliasing).
    tb_coeff=np.zeros((norb,norb,z['nk']),complex)
    for q in range(z['nk']):
        ps=np.where(k==q)[0]
        for p in ps:
            tb_coeff[:,:,q]+=np.outer(up[:,p],down[:,p].conj())*sources[0][p]
    rot=z['torque'][:norb,norb:,:]
    def weighted_matrix_norm(x):return float(np.sqrt(2*np.einsum('ijk,ijk,k->',x.conj(),x,weights/weights.sum()).real))
    sector={}
    for name,mask in [('d_sector',(np.arange(norb)[:,None]>=4)&(np.arange(norb)[None,:]>=4)),
                      ('non_d_or_cross',~((np.arange(norb)[:,None]>=4)&(np.arange(norb)[None,:]>=4))),
                      ('orbital_diagonal',np.eye(norb,dtype=bool)),('orbital_offdiagonal',~np.eye(norb,dtype=bool))]:
        rr=rot*mask[:,:,None];bb=tb_coeff*mask[:,:,None]
        sector[name]={'rotation_norm':weighted_matrix_norm(rr),'local_source_norm':weighted_matrix_norm(bb),
                      'difference_norm':weighted_matrix_norm(rr-bb)}
    onsite_rot=np.einsum('ijk,k->ij',rot,weights/weights.sum())
    onsite_b=np.einsum('ijk,k->ij',tb_coeff,weights/weights.sum())
    sector['R0_onsite']={'rotation_norm':float(np.sqrt(2)*np.linalg.norm(onsite_rot)),
                         'local_source_norm':float(np.sqrt(2)*np.linalg.norm(onsite_b)),
                         'difference_norm':float(np.sqrt(2)*np.linalg.norm(onsite_rot-onsite_b))}
    sector['R_nonzero_hopping']={'rotation_norm':weighted_matrix_norm(rot-onsite_rot[:,:,None]),
                                'local_source_norm':weighted_matrix_norm(tb_coeff-onsite_b[:,:,None]),
                                'difference_norm':weighted_matrix_norm(rot-tb_coeff-(onsite_rot-onsite_b)[:,:,None])}
    core_ratio=point_norms['core']/point_norms['valence']
    # Classifications use the actual residuals. No parameter is fitted or fed
    # back into any production object; the raw numbers remain the evidence.
    core_policy='CORE-CONTRIBUTES-BUT-NOT-DOMINANT' if core_ratio>1.e-10 else 'CORE-NEGLIGIBLE'
    classification='BLOCKED — LOCAL BXC TO LMTO TRANSVERSE-VERTEX MAPPING NOT CLOSED'
    if coeff['full_D'][-1]>ALGEBRA_TOL:classification='FAIL — COEFFICIENT-SPACE SPIN-FLIP WARD IDENTITY BROKEN'
    elif eta_rows[-1]['compact']['valence']['relative'] < ALGEBRA_TOL and eta_rows[-1]['compact']['total']['relative'] > ALGEBRA_TOL:
        classification='BLOCKED — FROZEN-CORE WARD CONTRACT NOT CLOSED';core_policy='CORE-DOMINANT-WARD-MISMATCH'
    elif torque_relative<ALGEBRA_TOL:
        classification='BLOCKED — RESPONSE BASIS / RADIAL REPRESENTATION BREAKS WARD IDENTITY'
        if eta_rows[-1]['compact']['total']['relative']<ALGEBRA_TOL and eta_rows[-1]['radial_production']['valence']['relative']<ALGEBRA_TOL:
            classification='PASS — ALSDA RESPONSE REPRESENTATION CLOSURE ESTABLISHED'
    report={'EF_Ry':float(z['ef']),'temperature_K':float(z['temp']),'k_mesh':[4,4,4],'k_count':int(z['nk']),
            'SCF_moment_muB':float(z['moment']),'coefficient_valence_moment':float(rho_moment),
            'M_SR':float(4*np.pi*np.sum(z['w']*z['sr'])),
            'M_P':{name:float(4*np.pi*np.sum(z['w']*z[key])) for name,key in [('total','total'),('valence','val'),('core','core')]},
            'point_metric_norms':point_norms,'compact_metric_norms':compact_norms,'core_to_valence_norm_ratio':float(core_ratio),
            'core_90_percent_abs_moment_radius_bohr':float(core90),'core_maximum_density_radius_bohr':float(z['radius'][core_max_index]),
            'core_maximum_density':float(z['core'][core_max_index]),'core_additivity_max':float(np.max(abs(z['total']-z['val']-z['core']))),
            'core_difference_statistics':point_stats(targets['core'],targets['valence']),
            'core_difference_relative_norm_denominator':'valence target; independent of eta',
            'coefficient_eta_Ry':ETA.tolist(),'coefficient_relative_residuals':coeff,
            'rotation_norm':torque_norm,'local_Bxc_source_norm':source_norm,'torque_source_relative_mismatch':torque_relative,
            'radial_Pauli_source_norm':norm(sources[1]),'torque_Pauli_source_relative_mismatch':norm(tband-sources[1])/torque_norm,
            'HOH_derivative_norm':norm(tband-tbaseband),'no_HOH_derivative_norm':norm(tbaseband),
            'site_offdiagonal_norm':0.0,'site_offdiagonal_note':'one-site primitive cell; translated bonds are reported separately',
            'torque_sectors':sector,'local_field_compact_projection_relative_error':float(projection_source_error),
            'rotation_commutator_max_error':float(direct_rotation_error),'Hamiltonian_fixture_max_error':float(z['adapter']),
            'accepted_eigenpair_max_error':float(eigen_error),'eta_rows':eta_rows,'core_policy':core_policy,'classification':classification,
            'production_physics_modified':False,'Goldstone_correction':'OFF'}
    output.mkdir(parents=True,exist_ok=True)
    (output/'diagnostics.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(output/'profiles.npz',radius=z['radius'],weights=z['w'],m_SR=z['sr'],m_P_total=z['total'],
             m_P_valence=z['val'],m_P_core=z['core'],bxc=z['bxc'],**profiles)
    print(json.dumps({key:value for key,value in report.items() if key!='eta_rows'},indent=2))
    print('eta, coeff, compact_total, compact_valence, radial_valence, fullH_compact_valence')
    for i,row in enumerate(eta_rows):
        print(row['eta_Ry'],coeff['full_D'][i],row['compact']['total']['relative'],row['compact']['valence']['relative'],
              row['radial_production']['valence']['relative'],row['full_H_source']['compact_valence']['relative'])
    assert projection_source_error<ALGEBRA_TOL
    assert eigen_error<ALGEBRA_TOL and z['adapter']<ALGEBRA_TOL and direct_rotation_error<ALGEBRA_TOL
    assert coeff['full_D'][-1]<ALGEBRA_TOL
    assert all(coeff[name][-1]>ALGEBRA_TOL for name in coeff if name!='full_D')
    return report

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--binary',type=Path)
    parser.add_argument('--scratch',type=Path)
    parser.add_argument('--dump',type=Path)
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    self_test()
    if args.self_test:return
    if args.scratch is None:parser.error('--scratch required')
    args.scratch=args.scratch.resolve();args.scratch.mkdir(parents=True,exist_ok=True)
    path=args.dump
    if path is None:
        if args.binary is None:parser.error('--binary or --dump required')
        case=ROOT/'tests/integration/tddft_driver_smoke'
        text=(case/'input_compact_alsda.nml').read_text()
        text=re.sub(r"database\s*=\s*'([^']+)'",lambda m:"database = '"+str((case/m[1]).resolve())+"/'",text)
        (args.scratch/'input.nml').write_text(text)
        path=args.scratch/'ward_inputs.bin'
        env={**os.environ,'OMP_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1',
             'RSLMTO_LR_WARD_DUMP':str(path)}
        with (args.scratch/'run.log').open('w') as log:
            subprocess.run([str(args.binary.resolve()),'input.nml'],cwd=args.scratch,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1800)
        if 'Converged!' not in (args.scratch/'run.log').read_text():raise RuntimeError('accepted Fe gate did not converge')
    analyze(path.resolve(),args.scratch)

if __name__=='__main__':main()
