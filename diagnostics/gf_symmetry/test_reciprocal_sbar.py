import struct
import tempfile
import unittest
from pathlib import Path
import numpy as np
from reciprocal_sbar import project,encode,restore,reverse_map
from sbar_audit import read_sbar


class ReciprocalTests(unittest.TestCase):
    def setUp(self):
        self.rng=np.random.default_rng(170)

    def matrices(self):
        # Mn-Mn, Mn-F and F-F, with onsite and distinct offsite reverse blocks.
        keys=[(1,2,1.,0.,0.),(2,1,-1.,0.,0.),(1,3,0.,1.,0.),
              (3,1,0.,-1.,0.),(3,4,1.,-1.,0.),(4,3,-1.,1.,0.),(1,1,0.,0.,0.)]
        return {k:self.rng.normal(size=(9,9)) for k in keys}

    def test_pair_projection_shuffle_and_idempotence(self):
        a=self.matrices();b,rev=project(a)
        for p,(k,s) in enumerate(b.items()):
            np.testing.assert_array_equal(s,list(b.values())[rev[p]].T)
            np.testing.assert_array_equal(s,(a[k]+list(a.values())[rev[p]].T)/2)
        shuffled={k:a[k] for k in reversed(a)}
        c,_=project(shuffled)
        for k in b: np.testing.assert_array_equal(b[k],c[k])
        d,_=project(b)
        for k in b: np.testing.assert_array_equal(b[k],d[k])
        self.assertFalse(np.allclose(b[(1,2,1.,0.,0.)],b[(1,2,1.,0.,0.)].T))

    def test_missing_ambiguous_and_invalid_tolerance(self):
        a=self.matrices();a.pop((2,1,-1.,0.,0.))
        with self.assertRaisesRegex(ValueError,'Missing or ambiguous'): project(a)
        keys=list(self.matrices());keys.append((2,1,-1.+1e-7,0.,0.))
        with self.assertRaisesRegex(ValueError,'Missing or ambiguous'): reverse_map(keys)
        for tol in (0,-1,1e-3):
            with self.assertRaises(ValueError): project(self.matrices(),tol)

    def fixture(self,p):
        a=self.rng.normal(size=(9,9));b=self.rng.normal(size=(9,9))
        (p/'clust').write_text('II = 2\n0 0 0 1 1  1 0 0 2 2\n')
        (p/'view.sbar.').write_text(' VECTOR= 1 0 0 ICLUS= 1\n VECTOR= -1 0 0 ICLUS= 2\n')
        (p/'input.nml').write_text('alat = 1\n')
        marker=struct.pack('<I',72)
        data=b''.join(marker+np.asarray(row,dtype='<f8').tobytes()+marker for s in (a,b) for row in s)
        (p/'sbar.').write_bytes(data)
        return data

    def test_binary_roundtrip_exclusive_write_and_marker_preservation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);data=self.fixture(p)
            self.assertEqual(encode(data,read_sbar(p)),data)
            r=restore(p,p/'sbar.reciprocal')
            self.assertEqual((p/'sbar.').read_bytes(),data)
            self.assertEqual(r['after']['max'],0)
            self.assertTrue(r['byte_idempotent'])
            out=(p/'sbar.reciprocal').read_bytes()
            for q in range(18):
                self.assertEqual(data[q*80:q*80+4],out[q*80:q*80+4])
                self.assertEqual(data[q*80+76:q*80+80],out[q*80+76:q*80+80])
            restore(p,p/'twice',p/'sbar.reciprocal')
            self.assertEqual((p/'twice').read_bytes(),out)
            for target in (p/'sbar.',p/'sbar.reciprocal'):
                with self.assertRaises(FileExistsError): restore(p,target)

    def test_format_errors_and_ambiguous_geometry(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);data=self.fixture(p)
            invalid=[b'',data[:2],data[:-1],data[:80],struct.pack('>I',72)+data[4:],
                     data[:76]+b'xxxx'+data[80:]]
            for bad in invalid:
                (p/'sbar.').write_bytes(bad)
                with self.assertRaises(ValueError): read_sbar(p)
            (p/'sbar.').write_bytes(data)
            (p/'clust').write_text('II = 3\n0 0 0 1 1  1 0 0 2 2  1 0 0 3 3\n')
            with self.assertRaisesRegex(ValueError,'ambiguous destination'): read_sbar(p)


if __name__=='__main__': unittest.main()
