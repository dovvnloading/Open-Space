"""Deterministic synthetic controls. These are NOT ACT observational detections."""
import math
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mf_core import (gaussian_peak_psf, inject, matched_filter, peaks,
                     negative_tail_count, poisson_zero_rate_upper_95)


class MatchedFilterSyntheticTests(unittest.TestCase):
    def setUp(self):
        self.psf = gaussian_peak_psf(1.7)
        self.base = np.zeros((96, 96), dtype=float)
        self.weights = np.ones_like(self.base)

    def run_filter(self, observed, reference=None, io=None, ir=None, **kw):
        return matched_filter(observed, self.base if reference is None else reference,
                              self.weights if io is None else io,
                              self.weights if ir is None else ir,
                              self.psf, independent_reference=True, **kw)

    def test_rejects_nonindependent_reference(self):
        with self.assertRaisesRegex(ValueError, "leave-one-out"):
            matched_filter(self.base, self.base, self.weights, self.weights, self.psf)

    def test_rejects_pixel_mismatch_and_bad_weights(self):
        with self.assertRaisesRegex(ValueError, "identical"):
            self.run_filter(self.base, ir=np.ones((94, 94)))
        result = self.run_filter(self.base, io=np.zeros_like(self.base))
        self.assertFalse(np.any(result["eligible"]))
        self.assertTrue(np.all(np.isnan(result["snr_white"])))

    def test_exact_synthetic_amplitude(self):
        observed = inject(self.base, self.psf, 10.0, 48, 48)
        result = self.run_filter(observed)
        self.assertAlmostEqual(result["amplitude_native"][48, 48], 10., places=5)
        self.assertGreater(result["snr_white"][48, 48], 10.)
        self.assertLess(result["sigma_white_native"][48, 48], 1)
        self.assertEqual(len(peaks(result["snr_white"], threshold=8)), 1)

    def test_reference_noise_is_included(self):
        result = self.run_filter(self.base)
        high_reference_ivar = self.run_filter(self.base, ir=100 * self.weights)
        self.assertGreater(result["sigma_white_native"][48, 48],
                           high_reference_ivar["sigma_white_native"][48, 48])

    def test_edge_and_scan_holes_vetoed(self):
        obs = self.base.copy()
        obs[:16] = np.nan
        obs[30:40, 30:40] = np.nan
        result = self.run_filter(obs)
        self.assertFalse(result["eligible"][2, 48])
        self.assertFalse(result["eligible"][35, 35])
        self.assertTrue(result["eligible"][75, 75])
        self.assertTrue(np.isnan(result["amplitude_native"][35, 35]))

    def test_noise_calibration_is_white_noise_only(self):
        rng = np.random.default_rng(413)
        observed = rng.normal(size=self.base.shape)
        reference = rng.normal(size=self.base.shape)
        result = self.run_filter(observed, reference)
        good = result["snr_white"][16:-16, 16:-16]
        self.assertGreater(np.std(good), 0.75)
        self.assertLess(np.std(good), 1.25)
        self.assertEqual(len(peaks(result["snr_white"], threshold=7)), 0)
        self.assertEqual(negative_tail_count(result["snr_white"], threshold=7), 0)

    def test_synthetic_injection_recovery(self):
        rng = np.random.default_rng(631)
        recovered = 0
        for y, x in ((y, x) for y in (20, 48, 74) for x in (20, 48, 74)):
            observed = inject(rng.normal(size=self.base.shape), self.psf, 8.0, y, x)
            reference = rng.normal(size=self.base.shape)
            result = self.run_filter(observed, reference)
            recovered += bool(result["snr_white"][y, x] > 7)
        self.assertEqual(recovered, 9)  # synthetic test only, not survey completeness

    def test_rate_requires_validated_footprint(self):
        with self.assertRaisesRegex(ValueError, "Unvalidated"):
            poisson_zero_rate_upper_95(1000)
        with self.assertRaises(ValueError):
            poisson_zero_rate_upper_95(0, validated=True)
        self.assertAlmostEqual(poisson_zero_rate_upper_95(100, validated=True),
                               -math.log(0.05)/100)


if __name__ == "__main__":
    unittest.main()
