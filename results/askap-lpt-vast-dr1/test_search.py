import unittest,pandas as pd
from search import screen,dist
class Tests(unittest.TestCase):
 def test_wrap(self):self.assertLess(dist(359.999,0,.001,0),8)
 def test_frozen_cuts(self):
  base=dict(avg_compactness=1.1,n_selavy=3,n_siblings=0,n_measurements=10,
    max_snr=25,max_flux_peak=20,min_flux_peak=2,eta_peak=150,v_peak=.9)
  rows=[dict(source=1,**base),dict(source=2,**{**base,"n_siblings":1}),
        dict(source=3,**{**base,"avg_compactness":2.1}),
        dict(source=4,**{**base,"min_flux_peak":19,"eta_peak":1,"v_peak":.1}),
        dict(source=5,**{**base,"min_flux_peak":-1,"eta_peak":1,"v_peak":.1})]
  self.assertEqual(list(screen(pd.DataFrame(rows)).source),[1])
if __name__=="__main__":unittest.main()
