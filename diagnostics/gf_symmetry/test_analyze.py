import tempfile
from pathlib import Path
import unittest
import numpy as np
from symmetry import UH, UT, MIRROR, pauli, sectors, trace_parts
from analyze import read_dump, stage_report, domain_report, covariance_report, h_report


def synthetic():
    rng=np.random.default_rng(44)
    h=rng.normal(size=(36,36))+1j*rng.normal(size=(36,36))
    h=(h+h.conj().T)/2
    u=np.kron(np.eye(2),UH)
    h=(h+u@h.conj()@u.conj().T)/2
    return h


def dataset(h):
    d=dict(pair=1,i=1,j=2,ti=1,tj=2,r=np.array([.5,.5,0]),
           fermi=0.,local_axis='F',matrices={},integrands={})
    g=np.linalg.inv((.37+.19j)*np.eye(36)-h)
    for direction,raw,rev in [('ij',g[:18,18:],g[18:,:18]),('ji',g[18:,:18],g[:18,18:])]:
        d['matrices'][1,'g'+direction]=(.37,raw)
        for name,a in sectors(raw,rev).items():
            d['matrices'][1,'g'+name+direction]=(.37,a)
        for tag,a in zip(['gnmag','gx','gy','gz'],pauli(raw)):
            d['matrices'][1,tag[0]+('i' if direction=='ij' else 'j')+tag[1:]]=(.37,a)
    di,dj=np.diag(np.arange(1,10)),np.diag(np.arange(2,11))
    d['matrices'][1,'Delta_i']=(.37,di)
    d['matrices'][1,'Delta_j']=(.37,dj)
    values=[]
    for a,b in [('01','0'),('00','1')]:
        for mu in 'xyz':
            values.append(np.trace(di@d['matrices'][1,'g'+a+'ij'][1]@dj@d['matrices'][1,'g'+mu+b+'ji'][1]).real)
    d['integrands'][1]=np.array([.37,*values])
    return d


class AnalyzerTests(unittest.TestCase):
    def test_end_to_end_and_deliberate_failure(self):
        d=dataset(synthetic())
        report=stage_report(d)
        self.assertIsNone(report['first_observed_failure'])
        self.assertLess(report['sector_maxima']['x0ij']['fro'],1e-13)
        d['matrices'][1,'gx1ji'][1][0,0]+=1e-4
        report=stage_report(d)
        self.assertIn('sector_projection',report['failures'])
        self.assertIn('sector_mirror',report['failures'])
        self.assertIn('contraction',report['failures'])
        contribution=report['energies'][0]['trace']['SCx']['blocks']['allowed/forbidden']
        self.assertNotEqual(contribution,0)
        # Independent CC leakage, on an opposite-parity index.
        d=dataset(synthetic())
        d['matrices'][1,'gx0ji'][1][2,0]+=1e-4
        result=stage_report(d)['energies'][0]['trace']['CCx']
        self.assertNotEqual(result['blocks']['allowed/forbidden'],0)

    def test_domain_and_mismatched_grid(self):
        h=synthetic()
        u=np.kron(np.eye(2),UT)
        a,b=dataset(h),dataset(u@h.conj()@u.conj().T)
        report=domain_report(a,b)
        self.assertLess(report['real_integrands_even']['fro'],1e-11)
        self.assertLess(report['energies'][0]['sectors']['x1ij']['fro'],1e-13)
        b['integrands'][1][0]+=.1
        with self.assertRaises(ValueError): domain_report(a,b)

    def test_dump_roundtrip(self):
        d=dataset(synthetic())
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'dump.dat'
            with p.open('w') as f:
                f.write('P 1 1 2 1 2 .5 .5 0\nF 0 F\n')
                f.write('I 1 '+' '.join(map(str,d['integrands'][1]))+'\n')
                for (n,tag),(e,a) in d['matrices'].items():
                    f.write(f'M {n} {e} {tag} {len(a)}\n')
                    for row in a:
                        f.write(' '.join(f'{x.real:.17e} {x.imag:.17e}' for x in row)+'\n')
            parsed=read_dump(p)
            self.assertIsNone(stage_report(parsed)['first_observed_failure'])
            self.assertEqual(covariance_report([parsed])[0]['status'],'mapped bond missing')
            self.assertEqual(covariance_report([parsed],neel='110')[0]['status'],'not a symmetry of this domain')

    def test_actual_h_action_reconstruction(self):
        h=synthetic()
        a,b=3.,.7
        scaled=(h-b*np.eye(36))/a
        wplus=(scaled[:,:18]+scaled[:,18:])/np.sqrt(2)
        wminus=(scaled[:,:18]-scaled[:,18:])/np.sqrt(2)
        d=dict(i=1,j=2,pair=1,a=a,b=b,atoms={1:(1,np.array([0.,0.,0.])),2:(2,np.array([.5,.5,0.]))},matrices={})
        for atom in (1,2):
            s=slice((atom-1)*18,atom*18)
            d['matrices'][atom,'Hp']=(0,wplus[s,:])
            d['matrices'][atom,'Hm']=(0,wminus[s,:])
        r=h_report(d)
        self.assertLess(r['pair_hermiticity']['fro'],1e-13)
        self.assertLess(max(x['mirror']['fro'] for x in r['columns']),1e-13)
        d['matrices'][1,'Hp'][1][0,0]+=1e-3j
        self.assertGreater(max(x['mirror']['fro'] for x in h_report(d)['columns']),1e-4)


if __name__=='__main__': unittest.main(verbosity=2)
