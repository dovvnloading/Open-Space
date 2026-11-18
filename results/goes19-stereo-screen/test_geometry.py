import unittest
import numpy as np
from geometry import A,B,surface,intersect,triangulate,altitude_family
class GeometryTests(unittest.TestCase):
 def test_offset_ellipsoid_and_range(self):
  s=surface(-75,0,A+35786.0234,B+35786.0234); p=surface(-100,10,A+14,B+6); q=intersect(s,p,A+80,B+80)
  self.assertAlmostEqual((q[0]**2+q[1]**2)/(A+80)**2+q[2]**2/(B+80)**2,1,places=11)
  family=altitude_family(s,-100,10,[60,80,100]); self.assertGreater(np.linalg.norm(family[0]-family[2]),30)
 def test_known_2023_peak(self):
  s1=surface(-75.1999969,0,A+35786.0234,B+35786.0234);s2=surface(-137,0,A+35786.0234,B+35786.0234)
  p,miss=triangulate(s1,(-111.04264605,-11.90889803),s2,(-110.05316415,-11.91404661))
  np.testing.assert_allclose(p,[-2212.183436,-5929.958195,-1304.251066],atol=0.001);self.assertLess(miss,1.1)
 def test_reject_miss_and_bad_height(self):
  with self.assertRaises(ValueError): intersect([42164,0,0],[42164,1,0],A+80,B+80)
  with self.assertRaises(ValueError): altitude_family([42164,0,0],-100,10,[121])
if __name__=='__main__': unittest.main()
