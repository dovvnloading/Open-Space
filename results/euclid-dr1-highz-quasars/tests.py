import unittest
from run import mag,sep_arcsec,reasons,toy_p,toy_tracks
class Verify(unittest.TestCase):
 def test_ab(self):self.assertAlmostEqual(mag(1),23.9)
 def test_sep(self):self.assertAlmostEqual(sep_arcsec(0,0,.001,0),3.6,places=3)
 def test_no_vis_not_non_detection(self):
  a={"flux_y_3fwhm_aper":"0.3","fluxerr_y_3fwhm_aper":"0.1",
   "flux_j_3fwhm_aper":"1.5","fluxerr_j_3fwhm_aper":"0.1",
   "flux_h_3fwhm_aper":"1.3","fluxerr_h_3fwhm_aper":"0.1",
   "ellipticity":"0.1","flag_j":"0","flag_h":"0"}
  self.assertIn("VIS_MISSING_OR_DETECTION",reasons(a))
  a.update({"flux_vis_3fwhm_aper":"0","fluxerr_vis_3fwhm_aper":"0.1"})
  self.assertEqual(reasons(a),[])
 def test_toy_probability_nonclaim(self):
  self.assertIsNone(toy_p(None,0,toy_tracks()))
  self.assertTrue(0<=toy_p(1,0,toy_tracks())<=1)
if __name__=="__main__":unittest.main()
