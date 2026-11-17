import unittest,numpy as np
from calibration import draw
class Tests(unittest.TestCase):
 def test_reproducible_and_finite(self):
  self.assertEqual(draw(123),draw(123));self.assertTrue(np.isfinite(draw(456,'narrow')))
 def test_generator_validation(self):
  self.assertGreaterEqual(draw(789,inject=.2),0.)
if __name__=='__main__':unittest.main()
