import unittest
import numpy as np
from symmetry import (C, M, P, V, UT, UH, UM, SIGMA, OPS, MIRROR, BARE_T,
                      project, rorb, split, pauli, spinor, sectors, residual, trace_parts)


class SymmetryTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(1729)

    def matrix(self, n=9):
        return self.rng.normal(size=(n,n))+1j*self.rng.normal(size=(n,n))

    def test_source_indexed_reciprocity(self):
        a = self.matrix()
        b = np.empty_like(a)
        for k in range(1,10):
            for j in range(1,10):
                l1, l2 = int(np.sqrt(k-.9)), int(np.sqrt(j-.9))
                k0, j0 = l1*(l1+1)+1, l2*(l2+1)+1
                b[k-1,j-1] = (-1)**(k+j)*a[2*j0-j-1,2*k0-k-1]
        np.testing.assert_allclose(rorb(a), b)
        np.testing.assert_allclose(rorb(rorb(a)), a)
        # R is an anti-automorphism, not a conjugate-linear map.
        b = self.matrix()
        np.testing.assert_allclose(rorb(a@b), rorb(b)@rorb(a), atol=2e-14)
        np.testing.assert_allclose(rorb(1j*a), 1j*rorb(a))

    def test_basis_and_group(self):
        np.testing.assert_allclose(V.conj().T@V, np.eye(9), atol=3e-16)
        np.testing.assert_allclose(V.conj().T@V.conj(), C, atol=3e-16)
        cart_mirror = np.diag([1,1,1,-1,1,-1,-1,1,1])
        np.testing.assert_allclose(V.conj().T@cart_mirror@V, M, atol=3e-16)
        np.testing.assert_allclose(C@C, np.eye(9))
        np.testing.assert_allclose(C@M, M@C)
        np.testing.assert_allclose(UT@UT.conj(), -np.eye(18))
        np.testing.assert_allclose(UH@UH.conj(), np.eye(18))
        for _,u,anti in OPS.values():
            np.testing.assert_allclose(u@u.conj().T, np.eye(18), atol=5e-16)

    def test_nonzero_gx0_gx1_zero_traces(self):
        di, dj = np.diag(np.arange(1,10)), np.diag(np.arange(2,11))
        gx0, gx1 = project(self.matrix(),1), project(self.matrix(),-1)
        g00, g01 = project(self.matrix(),1), project(self.matrix(),-1)
        self.assertGreater(np.linalg.norm(gx0), 1)
        self.assertGreater(np.linalg.norm(gx1), 1)
        self.assertAlmostEqual(np.trace(di@g00@dj@gx1).real, 0)
        self.assertAlmostEqual(np.trace(di@g01@dj@gx0).real, 0)
        # A known same-parity leakage in Gx1 produces a detectable real trace.
        g00[0,0] = 1
        gx1[0,0] = 1e-5
        result = trace_parts(di,dj,g00,gx1,1,-1)
        self.assertAlmostEqual(result['total'], 2e-5)
        self.assertAlmostEqual(result['blocks']['allowed/forbidden'], 2e-5)
        self.assertGreater(residual(M@gx1@M,-gx1)['fro'], 1e-5)

    def test_allowed_z_traces(self):
        a = project(self.matrix(),1)
        b = project(self.matrix(),-1)
        self.assertGreater(np.trace(a@a.conj().T).real, 0)
        self.assertGreater(np.trace(b@b.conj().T).real, 0)

    def test_fixed_bond_mirror_for_both_directions(self):
        raw = [self.matrix() for _ in range(4)]
        signs = [1,1,1,-1]
        reverse = [s*rorb(M@g@M) for s,g in zip(signs,raw)]
        for aa,bb in [(raw,reverse),(reverse,raw)]:
            for mu,a,b,s in zip('0xyz',aa,bb,signs):
                even,odd = split(a,b)
                np.testing.assert_allclose(M@even@M,s*even)
                np.testing.assert_allclose(M@odd@M,-s*odd)

    def test_retarded_transpose_and_domain_parity(self):
        # A Hermitian 2-site H invariant under m_z*T, including SOC-like terms.
        h = self.matrix(36)
        h = (h+h.conj().T)/2
        site_u = np.kron(np.eye(2),UH)
        h = (h+site_u@h.conj()@site_u.conj().T)/2
        z = .37 + .19j
        g = np.linalg.inv(z*np.eye(36)-h)
        a,b = g[:18,18:], g[18:,:18]
        np.testing.assert_allclose(a,UH@b.T@UH.conj().T,atol=3e-15)
        self.assertGreater(np.linalg.norm(a-UH@b.conj()@UH.conj().T), .1)
        for name,mat in sectors(a,b).items():
            np.testing.assert_allclose(M@mat@M,MIRROR[name]*mat,atol=3e-15)
        u = np.kron(np.eye(2),UT)
        hm = u@h.conj()@u.conj().T
        gm = np.linalg.inv(z*np.eye(36)-hm)
        np.testing.assert_allclose(gm[:18,18:],UT@b.T@UT.conj().T,atol=3e-15)
        plus,minus = sectors(a,b),sectors(gm[:18,18:],gm[18:,:18])
        for name in plus:
            np.testing.assert_allclose(minus[name],BARE_T[name]*plus[name],atol=3e-15)

    def test_nonhermitian_mirror_does_not_imply_retarded_transpose(self):
        # A controlled counterexample to the missing simplifying assumption.
        # This is synthetic algebra, not a modification of a production Hamiltonian.
        h = self.matrix(36)
        u = np.kron(np.eye(2),UH)
        h = (h+u@h.conj()@u.conj().T)/2
        self.assertGreater(np.linalg.norm(h-h.conj().T),1)
        z = .37+.19j
        g = np.linalg.inv(z*np.eye(36)-h)
        lower = np.linalg.inv(z.conjugate()*np.eye(36)-h)
        np.testing.assert_allclose(g,u@lower.conj()@u.conj().T,atol=5e-13)
        self.assertGreater(np.linalg.norm(g-u@g.T@u.conj().T),1)

    def test_numerical_nonreciprocity_can_produce_odd_traces(self):
        # Synthetic counterexample, not an SOC-off material calculation.
        cart_p = np.diag([1,1,1,-1,1,-1,-1,1,1])
        def cart_block():
            a=self.rng.normal(size=(9,9))
            return (a+cart_p@a@cart_p)/2
        a,b=cart_block(),cart_block()
        a=(a+a.T)/2; b=(b+b.T)/2
        hopping=cart_block()
        orbital=np.block([[a,hopping],[hopping.T+1e-3*cart_block(),b]])
        field=np.diag(np.r_[np.linspace(.1,.9,9),-np.linspace(.2,1.0,9)])
        # Construct site-major spin blocks, each 18 x 18.
        transform=np.kron(np.eye(2),V)
        def hamiltonian(sign):
            h=np.zeros((36,36),complex)
            for i in range(2):
                for j in range(2):
                    sl=slice(i*9,(i+1)*9); sr=slice(j*9,(j+1)*9)
                    block=np.kron(SIGMA[0],orbital[sl,sr])+sign*np.kron(SIGMA[1],field[sl,sr])
                    h[i*18:(i+1)*18,j*18:(j+1)*18]=transform.conj().T@block@transform
            return h
        vals=[]
        for sign in [1,-1]:
            h=hamiltonian(sign)
            u=np.kron(np.eye(2),UH)
            np.testing.assert_allclose(h,u@h.conj()@u.conj().T,atol=1e-15)
            g=np.linalg.inv((.37+.19j)*np.eye(36)-h)
            aa,bb=sectors(g[:18,18:],g[18:,:18]),sectors(g[18:,:18],g[:18,18:])
            vals.append([np.trace(aa['00']@bb['x1']).real,np.trace(aa['01']@bb['x0']).real])
        self.assertGreater(np.linalg.norm(vals[0]),1e-6)
        np.testing.assert_allclose(vals[0],-np.array(vals[1]),atol=1e-14)

    def test_pauli_roundtrip_and_operation_squares(self):
        a = self.matrix(18)
        np.testing.assert_allclose(spinor(pauli(a)),a,atol=5e-16)
        for _,u,anti in OPS.values():
            transformed = u@(a.T if anti else a)@u.conj().T
            transformed = u@(transformed.T if anti else transformed)@u.conj().T
            np.testing.assert_allclose(transformed,a,atol=4e-15)


if __name__ == '__main__':
    unittest.main(verbosity=2)
