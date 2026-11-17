import unittest,math
from guard import gate,predict_ray
class Tests(unittest.TestCase):
 def test_phase_seal(self):
  for sat in (18,19):
   with self.assertRaises(ValueError): gate('development',sat,[0],[0])
 def test_time_quality(self):
  self.assertEqual(gate('synthetic',19,[0,10,119.9],[0,1,0]),[0,119.9])
  for t,q in [([120],[0]),([2,1],[0,0]),([1],[2])]:
   with self.assertRaises(ValueError): gate('synthetic',19,t,q)
 def test_shell_geometry(self):
  for height in (60.,80.,100.):
   point=predict_ray((42164.,0.,0.),(6371.,0.,0.),height)
   self.assertAlmostEqual(math.sqrt(sum(v*v for v in point)),6371.+height,places=7)
  with self.assertRaises(ValueError): predict_ray((42164.,0.,0.),(42164.,1.,0.),80.)
  with self.assertRaises(ValueError): predict_ray((42164.,0.,0.),(6371.,0.,0.),70.)
if __name__=='__main__': unittest.main()
