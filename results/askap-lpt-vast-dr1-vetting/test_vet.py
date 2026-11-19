import unittest,importlib.util
from pathlib import Path
spec=importlib.util.spec_from_file_location("vet",Path("results/askap-lpt-vast-dr1-vetting/vet.py"))
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
class VetTests(unittest.TestCase):
 def test_dist(self):
  self.assertAlmostEqual(v.angular(1,2,1,2),0)
  self.assertTrue(359<v.angular(0,0,.1,0)<361)
 def test_frozen(self):
  s,m,ids,on=v.shortlist()
  self.assertEqual((len(ids),len(m),len(on)),(60,693,14))
 def test_no_inferred_counterpart(self):
  self.assertEqual(v.classify("x",{}, {},{}),"unresolved_incomplete_native_images")
if __name__=="__main__":unittest.main()
