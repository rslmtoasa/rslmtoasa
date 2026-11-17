"""Source-extracted algebra/routing checks, NOT a compiled Fortran solver test.
Run: python tests/regression/test_octupole_operators.py
Requires numpy and f90nml (repository requirements).
"""
from pathlib import Path
import re
import unittest

import f90nml
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
H = (ROOT / 'source/hamiltonian_build.f90').read_text()
R = (ROOT / 'source/recursion_transport.f90').read_text()
M = (ROOT / 'source/math.f90').read_text()
PAIRS = [('x', 'y'), ('y', 'z'), ('z', 'x'), ('x', 'z'), ('y', 'x'), ('z', 'y')]
PURE = ['O' + a + b + b for a, b in PAIRS] + ['Oxyz']


def subroutine(source, name):
    start = re.search(r'^   (?:module )?subroutine ' + name + r'\(', source, re.M)
    return source[start.end():].split('   end subroutine ' + name, 1)[0]


def evaluate(expr, env):
    return eval(expr.replace('_rp', ''), {'__builtins__': {}}, env)


def matrices():
    """Read actual L constants and actual hcpx transform assignments."""
    cart = []
    for axis in 'xyz':
        block = M.split('L_' + axis + ' = reshape((/', 1)[1]
        values = block.split('/),', 1)[0]
        values = values.replace('&', '').replace('_rp', '').replace('sqrt_three', str(np.sqrt(3)))
        cart.append(-1j * np.array([float(v) for v in values.split(',')]).reshape((9, 9), order='F'))
    v = np.zeros((9, 9), complex)
    for i, j, expr in re.findall(r'^\s*v\((\d+), (\d+)\) = ([^\n]+)', subroutine(M, 'hcpx'), re.M):
        v[int(i)-1, int(j)-1] = evaluate(expr, {'cone': 1., 'const': 1/np.sqrt(2), 'i_unit': 1j})
    np.testing.assert_allclose(v.conj().T @ v, np.eye(9), atol=1e-14)
    angular = [v.conj().T @ a @ v for a in cart]
    env = {'matmul': np.matmul, **{'mL'+a: l for a, l in zip('xyz', angular)}}
    selector = subroutine(H, 'select_octupole_operator')
    for name, expr in re.findall(r'^\s*(L[xyz]2|L2)\s*= (.+)$', selector, re.M):
        env[name] = evaluate(expr, env)
    # Join Fortran continuation lines without changing the expressions.
    selector = re.sub(r'&\s*\n\s*', '', selector)
    ops, rejected = {}, set()
    for aliases, body in re.findall(r"case\(([^\n]+)\)(.*?)(?=\n\s*case |\n\s*case\(|\n\s*end select)", selector, re.S):
        names = re.findall(r"'([^']+)'", aliases)
        expr = re.search(r'^\s*O9 = ([^\n]+)', body, re.M)
        if expr:
            for name in names:
                ops[name] = evaluate(expr.group(1), env)
        elif 'g_logger%fatal' in body:
            rejected.update(names)
    return angular, ops, rejected


def comm(a, b):
    return a @ b - b @ a


class OctupoleTests(unittest.TestCase):
    def test_all_pure_components(self):
        angular, ops, _ = matrices()
        for name in PURE:
            with self.subTest(component=name):
                o = ops[name]
                np.testing.assert_allclose(o, o.conj().T, atol=1e-13)
                np.testing.assert_allclose(o[:4, :], 0, atol=1e-13)
                np.testing.assert_allclose(o[:, :4], 0, atol=1e-13)
                for shell in [slice(0, 1), slice(1, 4), slice(4, 9)]:
                    np.testing.assert_allclose(np.trace(o[shell, shell]), 0, atol=1e-13)
                    for a in angular:
                        np.testing.assert_allclose(np.vdot(a[shell, shell], o[shell, shell]), 0, atol=1e-13)
                casimir = sum(comm(a, comm(a, o)) for a in angular)
                np.testing.assert_allclose(casimir, 12 * o, atol=1e-12)
        self.assertEqual(np.linalg.matrix_rank(np.array([ops[n].reshape(-1) for n in PURE]), tol=1e-12), 7)

    def test_aliases_and_legacy_migration(self):
        angular, ops, rejected = matrices()
        L2 = sum(a @ a for a in angular)
        for a, b in PAIRS:
            name = 'O'+a+b+b
            for alias in [name.lower(), a+b+b, 'O'+a+'5'+b+'2r2', 'o'+a+'5'+b+'2r2']:
                np.testing.assert_array_equal(ops[name], ops[alias])
            mixed = 'M'+a+'3'+b+'2r2'
            np.testing.assert_array_equal(ops[mixed], ops[mixed.lower()])
            La = angular['xyz'.index(a)]
            np.testing.assert_allclose(ops[name] - ops[mixed], (2*L2/15 - np.eye(9)/10) @ La, atol=1e-13)
            self.assertIn('O'+a+'3'+b+'2r2', rejected)
            self.assertIn('o'+a+'3'+b+'2r2', rejected)
        for alias in ['oxyz', 'xyz']:
            np.testing.assert_array_equal(ops['Oxyz'], ops[alias])
        self.assertNotIn('z', ops)
        self.assertTrue(all(len(name) <= 10 for name in ops))

    def test_routing_and_guards(self):
        setup = subroutine(R, 'setup_kubo_operators')
        self.assertIn("set_kubo_operator_slot(linear_in, 'b')", setup)
        self.assertIn("set_kubo_operator_slot(linear_out, 'a')", setup)
        dispatch = subroutine(R, 'set_kubo_operator_slot')
        self.assertRegex(dispatch, r"case\('a'\)\s+pol = this%hamiltonian%pol_alpha")
        self.assertRegex(dispatch, r"case\('b'\)\s+pol = this%hamiltonian%pol_beta")
        self.assertIn("case('octupole', 'octupole_current', 'octupole_accumulation')", dispatch)
        self.assertIn('call select_octupole_operator(pol, q)', dispatch)
        self.assertIn("kind = 'local'", dispatch)
        apply = subroutine(R, 'apply_kubo_operator')
        self.assertIn('this%transport_dir_a, .false., psi, out, work)', apply)
        self.assertIn('this%transport_dir_b, .true., psi, out, work)', apply)
        moments = subroutine(R, 'compute_moments_stochastic')
        self.assertEqual(moments.count('call this%kubo_ham_vec_matmul('), 4)
        self.assertEqual(moments.count('call this%apply_kubo_operator('), 2)
        self.assertIn('Unknown linear-response operator type', dispatch)
        for name, filename in [('hamiltonian', 'hamiltonian_build'), ('control', 'control')]:
            source = (ROOT / ('source/'+filename+'.f90')).read_text()
            self.assertIn("g_logger%fatal('Invalid &"+name+' namelist', source)

    def test_example_inputs(self):
        _, ops, _ = matrices()
        for filename, output in [('octupole_hall.nml', 'octupole'), ('octupole_accumulation.nml', 'octupole_accumulation')]:
            nml = f90nml.read(ROOT / 'example/transport' / filename)
            self.assertEqual(nml['control']['linear_in'], 'charge')
            self.assertEqual(nml['control']['linear_out'], output)
            self.assertIn(nml['hamiltonian']['pol_alpha'], ops)
            self.assertFalse(nml['hamiltonian']['hoh'])
            for key in ['v_alpha', 'v_beta']:
                self.assertGreater(np.linalg.norm(nml['hamiltonian'][key]), 0)
            for group, values in nml.items():
                declarations = (ROOT / f'source/include_codes/namelists/{group}.f90').read_text()
                for name in values:
                    self.assertRegex(declarations, r'\b'+name+r'\b')


if __name__ == '__main__':
    unittest.main()
