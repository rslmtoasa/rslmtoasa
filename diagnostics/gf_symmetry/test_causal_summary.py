import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from causal_summary import complete,parity
from symmetry import SECTORS,BARE_T,MIRROR,P,UT


class CausalSummaryTests(unittest.TestCase):
    def test_diagonal_input_integrity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);cases=('110','1-1-0')
            manifest={'one_pair':None,'pair_count':None,'cases':{}}
            for case in cases:
                (root/case).mkdir();f=root/case/'input.nml';f.write_text('saved diagonal state\n')
                manifest['cases'][case]={'returncode':0,'potentials_unchanged':True,
                    'copied_hashes_before_pair_selection':{'input.nml':hashlib.sha256(f.read_bytes()).hexdigest()}}
            (root/'manifest.json').write_text(json.dumps(manifest))
            self.assertEqual(complete(root,cases),manifest)
            (root/'1-1-0/input.nml').write_text('changed diagonal state\n')
            with self.assertRaisesRegex(ValueError,'Input changed'): complete(root,cases)

    def test_completed_manifest_and_mutated_input_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            manifest={'one_pair':None,'pair_count':None,'cases':{}}
            for case in ('100','1-00'):
                (p/case).mkdir();f=p/case/'input.nml';f.write_text('frozen input\n')
                manifest['cases'][case]={'returncode':0,'potentials_unchanged':True,
                    'copied_hashes_before_pair_selection':{'input.nml':hashlib.sha256(f.read_bytes()).hexdigest()}}
            path=p/'manifest.json';path.write_text(json.dumps(manifest))
            self.assertEqual(complete(p),manifest)
            (p/'100/input.nml').write_text('changed input\n')
            with self.assertRaisesRegex(ValueError,'Input changed'): complete(p)
            (p/'100/input.nml').write_text('frozen input\n')
            manifest['cases']['100']['returncode']=-9;path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,'Incomplete run'): complete(p)
            manifest['cases']['100']['returncode']=0;manifest['one_pair']=1
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,'Full 24-pair'): complete(p)

    def test_known_even_scalar_and_odd_vector_leakage(self):
        a={'ti':1,'tj':2,'r':np.array([-.5,.5,0]),'fermi':0.,
           'matrices':{},'integrands':{1:np.array([-.1,1.,0.,3.,2.,0.,4.])}}
        b={**a,'matrices':{},'integrands':{1:np.array([-.1,-1.,0.,3.,-2.,0.,4.])}}
        for d in (a,b):
            for direction in ('ij','ji'):
                d['matrices'][1,'g'+direction]=(-.1,np.eye(18,dtype=complex))
        for k in SECTORS:
            allowed=P[:,None]*P[None,:]==MIRROR[k]
            for direction in ('ij','ji'):
                x=np.where(allowed,1.,.01).astype(complex)
                y=BARE_T[k]*x
                if k=='01': y=np.where(allowed,y,x)  # measured even forbidden scalar
                if k=='x1': y=np.where(allowed,y,-x) # measured odd forbidden vector
                a['matrices'][1,'g'+k+direction]=(-.1,x)
                b['matrices'][1,'g'+k+direction]=(-.1,y)
        r=parity([a],[b])
        self.assertEqual(r['01ij_forbidden']['difference_fro'],0)
        self.assertEqual(r['x1ji_forbidden']['sum_fro'],0)
        self.assertEqual(r['x0ji_allowed']['sum_fro'],0)
        self.assertEqual(r['integrand_CCx']['sum_fro'],0)
        self.assertEqual(r['integrand_SCz']['difference_fro'],0)
        self.assertGreater(r['01ij']['difference_fro'],0) # nominal sector parity fails
        for k in SECTORS:
            for direction in ('ij','ji'):
                tag='g'+k+direction
                b['matrices'][1,tag]=(-.1,BARE_T[k]*a['matrices'][1,tag][1])
        r=parity([a],[b])
        self.assertEqual(r['01ij']['difference_fro'],0)
        self.assertEqual(r['x1ji']['difference_fro'],0)

    def test_raw_retarded_reciprocal_time_reversal(self):
        rng=np.random.default_rng(42)
        a={'ti':1,'tj':2,'r':np.array([-.5,.5,0]),'fermi':0.,
           'matrices':{},'integrands':{1:np.zeros(7)}}
        b={**a,'matrices':{}}
        for direction in ('ij','ji'):
            raw=rng.normal(size=(18,18))+1j*rng.normal(size=(18,18))
            a['matrices'][1,'g'+direction]=(-.1,raw)
        for direction,reverse in (('ij','ji'),('ji','ij')):
            raw=a['matrices'][1,'g'+reverse][1]
            b['matrices'][1,'g'+direction]=(-.1,UT@raw.T@UT.conj().T)
            for k in SECTORS:
                a['matrices'][1,'g'+k+direction]=(-.1,np.zeros((9,9),complex))
                b['matrices'][1,'g'+k+direction]=(-.1,np.zeros((9,9),complex))
        result=parity([a],[b])
        self.assertEqual(result['raw_T_ij']['difference_fro'],0.)
        self.assertEqual(result['raw_T_ji']['difference_fro'],0.)
        b['matrices'][1,'gij'][1][0,0]+=1e-4
        result=parity([a],[b])
        self.assertAlmostEqual(result['raw_T_ij']['difference_fro'],1e-4)
        self.assertEqual(result['raw_T_ji']['difference_fro'],0.)


if __name__=='__main__': unittest.main()
