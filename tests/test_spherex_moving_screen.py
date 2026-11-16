import unittest
import numpy as np
from examples.spherex_moving_screen import aperture
class ApertureTests(unittest.TestCase):
 def test_constant_background(self):
  im=np.ones((60,60));v=np.ones_like(im);f=np.zeros(im.shape,dtype=int)
  self.assertEqual(aperture(im,v,f,30,30,2.5)['sum_mjysr'],0)
 def test_outlier_gate_is_explicit(self):
  im=np.ones((60,60));im[30,30]=10;v=np.ones_like(im);f=np.zeros(im.shape,dtype=int);f[30,30]=524288
  self.assertIsNone(aperture(im,v,f,30,30,2.5))
  self.assertEqual(aperture(im,v,f,30,30,2.5,True)['sum_mjysr'],9)
 def test_edge(self):
  a=np.ones((60,60));self.assertIsNone(aperture(a,a,np.zeros(a.shape,dtype=int),0,0,2.5))
if __name__=='__main__':unittest.main()
