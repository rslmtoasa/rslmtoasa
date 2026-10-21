import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from causal_summary import complete,parity
from symmetry import SECTORS,BARE_T,MIRROR,P


class CausalSummaryTests(unittest.TestCase):
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


if __name__=='__main__': unittest.main()
