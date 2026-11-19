"""Network-free guardrails for Euclid wide-survey selection and model semantics."""
import unittest
import math
import numpy as np
import pyarrow as pa
import scan
import vet

def make_tab(**override):
    x={"ra":193.285625,"dec":70.908972,"object_id":"control",
       "vis_det":0,"point_like_prob":None,"spurious_prob":.04,
       "ellipticity":.1,"mumax_minus_mag":-2.9,
       "flag_y":0,"flag_j":0,"flag_h":0,"flag_vis":0,
       "flux_vis_3fwhm_aper":0.0,"fluxerr_vis_3fwhm_aper":.10,
       "flux_y_3fwhm_aper":.70,"fluxerr_y_3fwhm_aper":.10,
       "flux_j_3fwhm_aper":2.00,"fluxerr_j_3fwhm_aper":.10,
       "flux_h_3fwhm_aper":2.00,"fluxerr_h_3fwhm_aper":.10}
    x.update(override)
    return pa.Table.from_pylist([x])

class SelectionTests(unittest.TestCase):
    def count(self,**kw):
        b,s,y,j,l=scan.mask_rows(make_tab(**kw))
        return bool(b[0]),bool(s[0]),float(y[0]),float(j[0])
    def test_standard_dropout_passes(self):
        b,s,y,j=self.count()
        self.assertTrue(b)
        self.assertTrue(s)
        self.assertGreater(y,.4)
        self.assertAlmostEqual(j,0)
    def test_vis_detection_fails_strict(self):
        self.assertEqual(self.count(flux_vis_3fwhm_aper=.22)[:2],(True,False))
    def test_pointlike_prob_none_does_not_deselect(self):
        self.assertTrue(self.count(point_like_prob=None)[1])
    def test_ellipticity_cut(self):
        self.assertFalse(self.count(ellipticity=.31)[1])
    def test_j_flag(self):
        self.assertFalse(self.count(flag_j=1)[1])
    def test_y_colour(self):
        self.assertEqual(self.count(flux_y_3fwhm_aper=1.6)[:2],(True,False))
    def test_negative_y_upper_limit(self):
        self.assertTrue(self.count(flux_y_3fwhm_aper=-.05)[1])
    def test_pathological_blue_remains_in_frozen_cut(self):
        b,s,y,j=self.count(flux_h_3fwhm_aper=.45)
        self.assertTrue(b and s)
        self.assertLess(j,-.5)
    def test_high_jh_rejected(self):
        self.assertFalse(self.count(flux_h_3fwhm_aper=4.0)[1])
    def test_low_j_snr_rejected(self):
        self.assertFalse(self.count(fluxerr_j_3fwhm_aper=.4)[1])

class VetTests(unittest.TestCase):
    def test_bayesian_prior_sensitivity(self):
        q=np.array([[.85,.1],[.95,.0],[1.1,.2],[1.3,.15],[1.5,.1]])
        d=np.array([[.0,-.5],[.3,-.1],[.4,-.35],[.6,-.5],[.75,-.3]])
        a=vet.p_q_given_q_or_d((1,.1),q,d,prior=.001)
        b=vet.p_q_given_q_or_d((1,.1),q,d,prior=.01)
        self.assertTrue(0<a<b<1)
    def test_too_few_controls_withholds_probability(self):
        r=vet.diagnostics(np.array([[1,0]]),np.array([[0,0]]))
        self.assertEqual(r["model_status"],"INSUFFICIENT_MATCHED_LABELED_CONTROLS")
    def test_gaia_movement_rejects_dwarf(self):
        d=vet.motion_and_wise({
            "Gaia_DR3":{"catalog_match":{"pmra":"20","pmra_error":"1"}},
            "catwise_2020":{"catalog_match":{"w1mpro":"16","w2mpro":"15"}}})
        self.assertTrue(d["foreground_motion_diagnostic"].startswith("REJECT"))
        self.assertAlmostEqual(d["W1_minus_W2_vega_mag"],1)

if __name__=="__main__":
    unittest.main()
