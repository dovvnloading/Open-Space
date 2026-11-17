"""Synthetic regression cases for bookkeeping and the instrument-control tests."""
import hashlib
import tempfile
import unittest
from pathlib import Path

import numpy as np
from vetting import (REST, background_hypotheses, forced_line, reference_mask,
                     select_controls, verified_review, vet)


def epoch(expid, bary=5018.7):
    return dict(expid=expid, fiber=100, barycentric_A=bary, lsf_fwhm_A=2.0,
                star_to_bary=1.0, topo_to_bary=1.00006, rv_error_kms=2.0, target_flux=12.,
                expected_sky_fibers=[112, 185, 301], trace_fibers_measured=[99, 101],
                sky_controls=[dict(fiber=f, flux=0.0, sigma=1.0) for f in (112,185,301)],
                good_bary_intervals_A=[[3600.,5800.]], other_peaks_bary_A=[])


class VettingTests(unittest.TestCase):
    def test_sites_and_cameras(self):
        c = select_controls([0,1,2,3,4], [False,False,True,True,True],
                            [359.999,90,0.0,359.995,31.], [0]*5, 0)
        self.assertEqual(c["sky_fibers"], [2,3,4])
        self.assertEqual(c["trace_fibers"], [1])

    def test_wavelength_factors_and_completeness(self):
        # An artificial large separation makes a double application of a correction detectable.
        line, h = REST["halpha"], 1.001
        a = epoch(1, line * h)
        a["star_to_bary"] = h
        result = vet({"expected_expids":[1,2], "exposures":[a,epoch(2)]})
        self.assertEqual(result["status"], "excluded")
        self.assertIn("halpha", result["measurements"][0]["stellar_window_hits"])
        self.assertTrue(reference_mask(np.array([line*h]), 2., a["topo_to_bary"], h, 2.)[0])

    def test_unmeasured_and_weak_sky(self):
        a = epoch(1)
        a["sky_controls"][0]["sigma"] = 8.
        result = vet({"expected_expids":[1,2], "exposures":[a]})
        self.assertEqual(result["status"], "pending")
        self.assertIn("complete_distinct_exposure_inventory", result["pending_checks"])
        self.assertIn("1/insensitive_sky_fiber_112", result["pending_checks"])
        a["sky_controls"][0].update(flux=12., sigma=2.)
        self.assertEqual(vet({"expected_expids":[1,2], "exposures":[a]})["status"], "excluded")

    def test_clipped_controls_stay_pending(self):
        a = epoch(1)
        a["sky_controls"] = a["sky_controls"][:2]
        a["trace_fibers_measured"] = [99]
        result = vet({"expected_expids":[1,2], "exposures":[a,epoch(2)]})
        self.assertIn("1/three_nearest_distinct_sky_fibers", result["pending_checks"])
        self.assertIn("1/adjacent_physical_traces", result["pending_checks"])

    def test_forced_photometry_and_bad_core(self):
        wave = np.linspace(5000,5020,41)
        profile = np.exp(-0.5*((wave-5010)/1.2)**2)
        profile /= np.sum(profile * 0.5)
        data = 12 + 0.01*(wave-5010) + 7*profile
        ivar, mask = np.full(len(wave), 100.0), np.zeros(len(wave), int)
        fit = forced_line(wave, data, ivar, mask, profile)
        self.assertAlmostEqual(fit["flux"], 7.0)
        mask[20] = 1
        self.assertIsNone(forced_line(wave, data, ivar, mask, profile))

    def test_background_companions(self):
        red = 1.1
        peaks = [REST["oiii_4959"]*red, REST["oiii_5007"]*red]
        hypo = background_hypotheses(peaks[1], 2, peaks, [(5000,8000)])
        oiii = next(h for h in hypo if h["anchor"]=="oiii_5007")
        self.assertAlmostEqual(oiii["z"], .1)
        self.assertEqual(oiii["companions"]["oiii_4959"]["peaks_within_one_fwhm_A"], [peaks[0]])
        self.assertNotIn("oiii_5007", oiii["companions"])

    def test_untracked_review_stays_pending(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root/"report.txt").write_text("synthetic test only")
            info = {"file":"report.txt","sha256":"bad","status":"clear","reason":"test"}
            self.assertIsNone(verified_review(info,root))
            info["sha256"] = hashlib.sha256((root/"report.txt").read_bytes()).hexdigest()
            self.assertEqual(verified_review(info,root),"clear")


if __name__ == "__main__":
    unittest.main()
