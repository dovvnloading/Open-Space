#!/usr/bin/env python3
"""Offline falsification tests for colour limits, epoch propagation and counts."""
import unittest
import sys
from pathlib import Path
import numpy as np
import pyarrow as pa
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from selection import cold_mask,cold_row
from bootstrap import propagate,angular_arcsec,subtype,cold_type
from referee import validate_row

def tab(**changes):
    row={"vis_det":0,"flag_j":0,"flag_h":0,"spurious_prob":.1,
         "ellipticity":.1,"flux_j_3fwhm_aper":10.,"fluxerr_j_3fwhm_aper":1.,
         "flux_y_3fwhm_aper":3.,"fluxerr_y_3fwhm_aper":1.,
         "flux_h_3fwhm_aper":5.,"fluxerr_h_3fwhm_aper":1.,
         "flux_vis_3fwhm_aper":0.,"fluxerr_vis_3fwhm_aper":1.,
         "object_id":123,"ra":45.,"dec":-30.}
    row.update(changes)
    return pa.Table.from_pylist([row])

class ColdDwarfTest(unittest.TestCase):
    def test_measured_colour_pass(self):
        b,s,y,h,l=cold_mask(tab())
        self.assertTrue(s[0])
    def test_h_nondetection_is_blue_upper_limit(self):
        t=tab(flux_h_3fwhm_aper=-1.)
        b,s,y,h,l=cold_mask(t)
        self.assertTrue(s[0])
        r=cold_row(t,0,y,h,l)
        self.assertTrue(r["J_minus_H_is_upper_limit"])
        self.assertLess(r["J_minus_H_for_selection"],0)
    def test_y_nondetection_lower_limit(self):
        t=tab(flux_y_3fwhm_aper=-1.)
        b,s,y,h,l=cold_mask(t)
        self.assertTrue(s[0])
        r=cold_row(t,0,y,h,l)
        self.assertTrue(r["Y_minus_J_is_lower_limit"])
    def test_optical_veto(self):
        self.assertFalse(cold_mask(tab(flux_vis_3fwhm_aper=5.))[1][0])
    def test_extended_or_wrong_shape_rejected(self):
        self.assertFalse(cold_mask(tab(ellipticity=.35))[1][0])
    def test_epoch_and_types(self):
        a,d,has_pm=propagate(180.,0.,2000.,1000.,0.)
        self.assertTrue(has_pm)
        self.assertAlmostEqual(angular_arcsec(a,d,180,0),24.25,places=3)
        self.assertTrue(cold_type("Y0.5 pec"))
        self.assertFalse(cold_type("T7.0"))
    def test_independent_mask(self):
        r=tab().to_pylist()[0]
        r["Y_minus_J_for_selection"]=-2.5*np.log10(.3)
        r["J_minus_H_for_selection"]=-2.5*np.log10(2)
        self.assertTrue(validate_row(r))
        r["flux_vis_3fwhm_aper"]=9
        self.assertFalse(validate_row(r))

if __name__=="__main__":unittest.main(verbosity=2)
