"""Regression checks for the PR #48 held-out radio result.

This deliberately runs the whole native-pixel FITS analysis rather than only
checking that the script imports or that a checked-in CSV looks plausible.
Only the repository's existing 'archive' extra (NumPy + Astropy) is required.
"""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "results" / "wise-j0458-heldout" / "recheck.py"


class WiseJ0458HeldoutReplay(unittest.TestCase):
    def test_four_native_fits_blends_and_controls(self):
        spec = importlib.util.spec_from_file_location("wise_j0458_heldout", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertEqual(len(list((module.PARENT / "inputs").glob("VLASS[1-4].fits"))), 4)

        # Avoid leaving an untracked measurements.json in the working tree.
        with tempfile.TemporaryDirectory() as temp_dir:
            module.HERE = Path(temp_dir)
            module.main()
            data = json.loads((Path(temp_dir) / "measurements.json").read_text())

        self.assertEqual(data["field_controls"], 64)
        epochs = data["epochs"]
        self.assertEqual([row["epoch"] for row in epochs], [1, 2, 3, 4])

        expected = (
            (0.79867, -0.09506),
            (-0.25660, 0.10901),
            (0.14823, 0.31591),
            (0.48848, -0.00889),
        )
        for row, (moving, stationary) in zip(epochs, expected):
            actual = row["moving_plus_static"]["amplitudes"]
            self.assertAlmostEqual(actual[0], moving, delta=0.004)
            self.assertAlmostEqual(actual[1], stationary, delta=0.004)
            self.assertGreater(row["moving_plus_static"]["n_pixels"], 250)
            self.assertTrue(row["moving_plus_static"]["sse_mjy2"] >= 0)
            self.assertIn("15", row["window_sensitivity"])
            self.assertIn("19", row["window_sensitivity"])
            self.assertIn("23", row["window_sensitivity"])

        # The 2017 and 2025 moving-only fits should outperform the
        # independently specified single stationary LoTSS catalog position.
        for idx in (0, 3):
            row = epochs[idx]
            self.assertGreater(row["sse_static_minus_moving"], 0)
            self.assertEqual(row["contrast_control"]["rank_ge_including_target"], 1)
            self.assertLess(abs(row["moving_plus_static"]["amplitudes"][1]), 0.16)

        self.assertEqual(data["endpoint_joint"]["n_controls"], 64)
        self.assertEqual(data["endpoint_joint"]["rank_ge_including_target"], 1)

        # Adversarial hypothesis: two fixed endpoint sites can mimic the
        # intermediate 2023 flux without either site moving.
        both_2023 = epochs[2]["two_stationary_endpoints"]["both"]
        self.assertAlmostEqual(both_2023["amplitudes"][0], 0.350, delta=0.006)
        self.assertAlmostEqual(both_2023["amplitudes"][1], 0.302, delta=0.006)

        morph = data["morphology"]
        self.assertEqual(morph["local_pointings_n"], 5)
        self.assertAlmostEqual(morph["maximum_catalog_pm_smear_arcsec"], 0.238,
                               delta=0.005)
        self.assertGreater(morph["major_over_smear"], 19)
        self.assertLess(morph["major_over_smear"], 21)

        self.assertAlmostEqual(
            data["variability_constant_flux_diagnostic"]["1-2-3-4"]["chi2"],
            22.6793424, delta=0.015)
        ra, dec = module.predicted("2028-09-06T14:10:22.8")
        self.assertAlmostEqual(ra, 74.72724859, delta=0.000002)
        self.assertAlmostEqual(dec, 64.58276496, delta=0.000002)


if __name__ == "__main__":
    unittest.main()
