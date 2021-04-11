import unittest

from ..mathematics.quaternions import *


q_case_1 = Quaternion(3.2, 17.4, 9.25, -8.47)


class TestQuaternions(unittest.TestCase):

    def test_sqrd_norm(self):

        self.assertAlmostEqual(Quaternion.zero().sqrdNorm(), 0.0)
        self.assertAlmostEqual(Quaternion.identity().sqrdNorm(), 1.0)
        self.assertAlmostEqual(Quaternion.i().sqrdNorm(), 1.0)
        self.assertAlmostEqual(Quaternion.j().sqrdNorm(), 1.0)
        self.assertAlmostEqual(Quaternion.k().sqrdNorm(), 1.0)

    def test_normalized(self):

        norm_quat = q_case_1.normalize()

        self.assertAlmostEqual(norm_quat.sqrdNorm(), 1.0)

        cnj_norm = norm_quat.conjugate
        inv_norm = norm_quat.inverse
        assert cnj_norm.isCloseTo(inv_norm)

        quat_1 = Quaternion(0.696, 0.322, -0.152, 0.624)
        assert areClose(quat_1.normalize().norm, 1.0)


if __name__ == '__main__':

    unittest.main()


