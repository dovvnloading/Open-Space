"""Synthetic CSC ellipse, original row identity, nonoverlap and sparse-count regressions."""
from unittest import TestCase

import numpy as np
from astropy.wcs import WCS

import deblend
from diagnostics import DTYPE


class DeblendTest(TestCase):
    def test_rotated_csc_ellipse_and_two_distinct_chandra_centers(self):
        self.assertAlmostEqual(deblend.nearest_ellipse(5,0,2,1,90),3)
        self.assertAlmostEqual(deblend.nearest_ellipse(5,0,2,1,0),4)
        self.assertEqual(deblend.nearest_ellipse(.9,0,2,1,0),0.)
        wcs=WCS(naxis=2)
        wcs.wcs.crpix=[100.,100.]
        wcs.wcs.crval=[250.,30.]
        wcs.wcs.cdelt=[-.5/3600,.5/3600]
        wcs.wcs.ctype=["RA---TAN","DEC--TAN"]
        ra,dec=wcs.all_pix2world([[108,100]],1)[0]
        xra,xdec=wcs.all_pix2world([[113,100]],1)[0]
        event={"aperture_x":100,"aperture_y":100,"aperture_radius_px":10,"ra_deg":xra,"dec_deg":xdec,
               "position_error_arcsec":.4}
        master={"ra":ra,"dec":dec,"err_ellipse_r0":1.,"err_ellipse_r1":.5,"err_ellipse_ang":90.}
        geom=deblend.geometry(wcs,event,master)
        self.assertAlmostEqual(geom["separation_px"],8.)
        self.assertAlmostEqual(geom["core_radius_px"],3.6)
        self.assertTrue(geom["extraction_circle"]["csc_center_in_circle"])
        self.assertFalse(geom["localization_circle"]["csc_ellipse_intersects_circle"])

    def test_phases_union_shares_and_unsupported_geometry(self):
        # Centers 0 and 10; r=4.5. A photon in the inter-core gap remains an original trial photon.
        pos=[(1,10),(2,9),(3,10),(10,0),(11,1),(12,3),(13,5),(14,10),(15,20),(20,10),(21,-10)]
        data=np.zeros(len(pos),dtype=DTYPE)
        for i,(t,x) in enumerate(pos):
            data[i]["row"],data[i]["t"],data[i]["x"]=50+i,t,x
            data[i]["y"]=0
        geom={"trial_xy":[0,0],"csc_xy":[10,0],"midpoint_xy":[5,0],"original_radius_px":8,
              "core_radius_px":4.5,"annulus_px":[14.,23.]}
        ev={"live_fraction":.99,"diagnostics":{"gti_met_s":[0.,30.]},"start_met_s":10.,"stop_met_s":20.,
            "before_counts":0,"on_counts":4,"after_counts":0}
        phases,tags,tests=deblend.measure(data,geom,ev,0.,30.)
        self.assertEqual([p["phase"] for p in tags if p["events_row_zero_based"]==53],["trigger"])
        self.assertEqual([p["region"] for p in tags if p["events_row_zero_based"]==56],["original_only"])
        self.assertEqual(len(tags),len({r["events_row_zero_based"] for r in tags}))
        self.assertEqual(tests["all_recorded"]["trigger_then_pooled_flanks_by_trial_and_csc"],[[3,1],[0,4]])
        self.assertIsNotNone(tests["all_recorded"]["fisher_two_sided_descriptive_p"])
        self.assertIsNone(tests["common_geometry_supported"]["fisher_two_sided_descriptive_p"])
        self.assertEqual(phases["trigger"]["common_geometry"]["supported_live_s"],0)
        self.assertEqual(phases["trigger"]["shared_annulus"],1)
        self.assertEqual(phases["trigger"]["original_photons_outside_both_cores"],1)
        self.assertEqual(phases["trigger"]["trial_fraction_among_both_cores"]["n"],4)
        self.assertEqual(deblend.fraction(0,0)["cp95"],[None,None])
        with self.assertRaisesRegex(ValueError,"original photometry"):
            deblend.measure(data,geom,{**ev,"on_counts":5},0.,30.)
