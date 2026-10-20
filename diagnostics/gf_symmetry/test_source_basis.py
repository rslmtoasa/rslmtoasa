"""Tie the abstract tests to actual numerical basis/SOC conventions in source."""
from pathlib import Path
import re
import unittest
import numpy as np
from symmetry import V, C, M, UH, UT, SIGMA


class SourceBasisTests(unittest.TestCase):
    def test_actual_orbital_operators_and_soc(self):
        text=(Path(__file__).resolve().parents[2]/'source/math.f90').read_text()
        operators=[]
        for name,sign in [('L_x',-1),('L_y',-1),('L_z',1)]:
            raw=re.search(r'parameter :: '+name+r'.*?reshape\(\(/(.*?)/\)',text,re.S).group(1)
            raw=raw.replace('&','').replace('_rp','').replace('sqrt_three',str(np.sqrt(3)))
            values=np.fromstring(raw.replace('\n',''),sep=',')
            self.assertEqual(len(values),81)
            a=V.conj().T@(values.reshape((9,9),order='F')*-1j)@V
            operators.append(a)
            np.testing.assert_allclose(M@a@M,sign*a,atol=5e-16)
            np.testing.assert_allclose(C@a.T@C,-a,atol=5e-16)
        soc=sum(np.kron(s,l)/2 for s,l in zip(SIGMA[1:],operators))
        np.testing.assert_allclose(soc,soc.conj().T,atol=5e-16)
        np.testing.assert_allclose(soc,UH@soc.conj()@UH.conj().T,atol=5e-16)
        np.testing.assert_allclose(soc,UT@soc.conj()@UT.conj().T,atol=5e-16)
        # Signs of printed m labels and eigenvalues differ; never infer R from labels alone.
        np.testing.assert_allclose(np.diag(operators[2]).real,[0,1,0,-1,2,1,0,-1,-2],atol=5e-16)


if __name__=='__main__': unittest.main()
