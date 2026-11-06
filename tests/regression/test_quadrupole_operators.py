"""Algebraic checks of the quadrupole selector's actual matrix expressions.

Run: python tests/regression/test_quadrupole_operators.py
Requires NumPy. This evaluates the source expressions, not compiled Fortran.
Checks use Cartesian L matrices: polynomial identities and traces are invariant
under the unitary cart2sph transformation used by the production selector.
"""
from pathlib import Path
import re

import numpy as np
import unittest

ROOT = Path(__file__).resolve().parents[2]


def operators():
    math_source = (ROOT / 'source/math.f90').read_text()
    env = {'matmul': np.matmul}
    angular = []
    for axis in 'xyz':
        block = math_source.split('L_' + axis + ' = reshape((/', 1)[1]
        values = block.split('/),', 1)[0]
        values = values.replace('&', '').replace('_rp', '')
        values = values.replace('sqrt_three', str(np.sqrt(3)))
        matrix = -1j * np.array([float(v) for v in values.split(',')]).reshape((9, 9), order='F')
        angular.append(matrix)
        env['mL' + axis] = matrix
    source = (ROOT / 'source/hamiltonian_build.f90').read_text()
    selector = source.split('   module subroutine select_quadrupole_operator(', 1)[1]
    selector = selector.split('   end subroutine select_quadrupole_operator', 1)[0]
    for name, expr in re.findall(r'^\s*(L[xyz]2|L2) = (.+)$', selector, re.M):
        env[name] = eval(expr.replace('_rp', ''), {'__builtins__': {}}, env)
    result = {}
    for aliases, expr in re.findall(r"case\(([^\n]+)\)\s*\n\s*Q9 = ([^\n]+)", selector):
        matrix = eval(expr.replace('_rp', ''), {'__builtins__': {}}, env)
        for alias in re.findall(r"'([^']+)'", aliases):
            result[alias] = matrix
    return angular, result


def check_pure_rank_two(label):
    angular, q = operators()
    matrix = q[label]
    np.testing.assert_allclose(matrix, matrix.conj().T, atol=1e-13)
    for block in (slice(0, 1), slice(1, 4), slice(4, 9)):
        np.testing.assert_allclose(np.trace(matrix[block, block]), 0, atol=1e-13)
    # A rank-k operator obeys sum_a [La,[La,Q]] = k(k+1) Q.
    casimir = sum(a @ a @ matrix - 2 * a @ matrix @ a + matrix @ a @ a for a in angular)
    np.testing.assert_allclose(casimir, 6 * matrix, atol=1e-13)


def check_cartesian_tensor_and_pdf_normalization():
    angular, q = operators()
    np.testing.assert_allclose(q['Qxx'] + q['Qyy'] + q['Qzz'], 0, atol=1e-13)
    np.testing.assert_allclose(q['Qx2y2'], (q['Qxx'] - q['Qyy']) / 2, atol=1e-13)
    np.testing.assert_allclose(q['Q3z2r2'], 1.5 * q['Qzz'], atol=1e-13)
    # PDF Q_z2 eigenvalues in p and d shells, independent of basis.
    np.testing.assert_allclose(np.linalg.eigvalsh(q['Q3z2r2'][1:4, 1:4]), [-1, .5, .5], atol=1e-13)
    np.testing.assert_allclose(np.linalg.eigvalsh(q['Q3z2r2'][4:9, 4:9]), [-3, -1.5, -1.5, 3, 3], atol=1e-13)
    for i, j, label in [(0, 1, 'Qxy'), (1, 2, 'Qyz'), (2, 0, 'Qzx')]:
        np.testing.assert_allclose(q[label], (angular[i] @ angular[j] + angular[j] @ angular[i]) / 2, atol=1e-13)
    for label in ['Qxx', 'Qyy', 'Qzz', 'Qxy', 'Qyz', 'Qzx', 'Qx2y2', 'Q3z2r2']:
        np.testing.assert_array_equal(q[label], q[label.lower()])
        np.testing.assert_array_equal(q[label], q[label[1:]])
    np.testing.assert_array_equal(q['Qxz'], q['Qzx'])


class QuadrupoleTests(unittest.TestCase):
    def test_rank_two_components(self):
        for label in ['Qxx', 'Qyy', 'Qzz', 'Qxy', 'Qyz', 'Qzx', 'Qx2y2', 'Q3z2r2']:
            with self.subTest(component=label):
                check_pure_rank_two(label)

    def test_normalization_and_aliases(self):
        check_cartesian_tensor_and_pdf_normalization()


if __name__ == '__main__':
    unittest.main()
