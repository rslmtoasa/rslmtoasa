import struct
import tempfile
from pathlib import Path
import unittest
import numpy as np
from sbar_audit import read_sbar,audit


class SbarTests(unittest.TestCase):
    def test_binary_reciprocity_and_known_perturbation(self):
        rng=np.random.default_rng(9)
        a=rng.normal(size=(9,9))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/'clust').write_text('II = 2\n0 0 0 1 1  1 0 0 2 2\n')
            (p/'view.sbar.').write_text(' VECTOR= 1 0 0 ICLUS= 1\n VECTOR= -1 0 0 ICLUS= 2\n')
            (p/'input.nml').write_text("alat = 1\nlabel(1) = 'A'\nlabel(2) = 'B'\n")
            for label,mom in [('A','1,0,0'),('B','-1,0,0')]:
                (p/(label+'_out.nml')).write_text('width_band(:,1) = .2,.3,.4\nwidth_band(:,2) = .1,.2,.3\nmom = '+mom+'\n')
            for eps in [0,1e-5]:
                b=a.T.copy();b[0,1]+=eps
                with (p/'sbar.').open('wb') as f:
                    for block in (a,b):
                        for row in block:
                            chunk=np.asarray(row,dtype='<f8').tobytes()
                            mark=struct.pack('<I',len(chunk));f.write(mark+chunk+mark)
                self.assertEqual(len(read_sbar(p)),2)
                r=audit(p)['blocks'][0]
                self.assertAlmostEqual(r['sbar_reciprocity']['max'],eps)
                if eps: self.assertGreater(r['bare_h_hermiticity']['max'],1e-8)
                else: self.assertLess(r['bare_h_hermiticity']['max'],1e-16)


if __name__=='__main__': unittest.main()
