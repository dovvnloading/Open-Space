import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest


@unittest.skipUnless(importlib.util.find_spec("astropy"), "archive extra required")
class FFIImageTests(unittest.TestCase):
    def test_reflected_shift_and_injected_star(self):
        import numpy as np
        from astropy.wcs import WCS
        from open_space.ffi import detect
        n = 17
        x,y = 8.0,7.0
        yy,xx = np.indices((n,n))
        wcs = WCS(naxis=2)
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs.wcs.crval, wcs.wcs.crpix = [124.,22.], [9.,9.]
        wcs.wcs.cdelt = [-0.00583,0.00583]
        time = np.arange(0,11,.04)
        rng = np.random.default_rng(2026)
        star = np.exp(-((x-xx)**2+(y-yy)**2)/(2*.83**2))
        ref = np.ones((n,n))*100.
        for cx,cy in ((3,3),(6,12),(14,5),(13,13)):
            ref += 80*np.exp(-((cx-xx)**2+(cy-yy)**2)/4)
        gy,gx = np.gradient(ref)
        shift = 0.12*np.sin(2*np.pi*time/0.33)
        cube = ref + shift[:,None,None]*gx + 0.5*shift[:,None,None]*gy
        cube += 4*np.cos(2*np.pi*time/1.39)[:,None,None]
        cube += rng.normal(0,.25,size=cube.shape)
        error = np.full_like(cube,.25)
        quiet_peaks,_ = detect((time,cube,error,wcs))
        self.assertEqual([p for p in quiet_peaks if not p["artifact_flags"]],[])
        flux = cube + (6*np.cos(2*np.pi*time/.721))[:,None,None]*star
        measured, counts = detect((time,flux,error,wcs))
        self.assertEqual(counts["eligible_pixels"], (n-4)**2)
        nearest = min(measured,key=lambda s:(s["centroid_x"]-x)**2+(s["centroid_y"]-y)**2)
        self.assertAlmostEqual(nearest["period_days"],.721,delta=.008)
        self.assertLess(np.hypot(nearest["centroid_x"]-x,nearest["centroid_y"]-y),0.5)
        self.assertEqual(nearest["artifact_flags"],[])

    def test_matching_is_one_to_one_and_vetoes_survive(self):
        from open_space.ffi import reconcile
        d = {"ra_deg": 0.001, "dec_deg":1.2,"period_days":1.0,"artifact_flags":[]}
        far = dict(d,ra_deg=70.)
        later = dict(d,artifact_flags=["bipolar_residual"])
        match, unmatched = reconcile([d,far],[later])
        self.assertEqual(unmatched,[])
        self.assertEqual([m["disposition"] for m in match],["image_artifact_veto","single_sector_image_peak"])
        match,_ = reconcile([d],[dict(d,period_days=1.1)])
        self.assertEqual(match[0]["disposition"],"frequency_inconsistent")

    def test_pilot_is_not_dependent_on_header_creation_time(self):
        from open_space.ffi_archive import request_url
        root = Path(__file__).resolve().parents[1]
        plan=json.loads((root/"data/ffi-fields.json").read_bytes())
        result=json.loads((root/"results/ffi-fields.json").read_bytes())
        for field, report in zip(plan["fields"],result["fields"]):
            for product, measured in zip(field["products"],report["products"]):
                self.assertEqual(product["science_sha256"],measured["science_sha256"])
                self.assertEqual(request_url(field,product["sector"]),measured["url"])
                self.assertEqual(len(measured["parent_ffi_list_sha256"]),64)
        self.assertTrue(result["fields"][0]["sources"][0]["candidate"])
