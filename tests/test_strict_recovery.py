import math
import unittest
from pathlib import Path
from unittest.mock import patch

from open_space import followup, spectra
from open_space.dataio import load_observation
from open_space.followup import measure


def request(index):
    return {"check": "strict_mask_recovery", "feature_id": f"F{index:04d}"}


class StrictRecoveryTests(unittest.TestCase):
    def test_rerun_is_shared_and_baseline_ids_are_preserved(self):
        rows = [dict(row, or_mask=0) for row in spectra.synthetic_rows()]
        redshift = 0.07
        original = spectra.analyze(rows, redshift)["features"]
        rows[original[0]["center_pixel"]]["or_mask"] = 4
        with patch.object(followup, "screen", wraps=followup.screen) as rerun:
            first, second = measure([request(1), request(2)], original, rows, redshift)
        self.assertEqual(rerun.call_count, 1)
        self.assertEqual(rerun.call_args.args[1], redshift)
        self.assertEqual(first["strict_search"]["outcome"], "unsearchable")
        self.assertEqual(second["feature_id"], "F0002")
        strict = second["strict_search"]
        self.assertEqual(strict["outcome"], "hit")
        self.assertEqual(strict["required_direction"], original[1]["direction"])
        self.assertAlmostEqual(strict["observed_angstrom"], original[1]["wavelength"])
        self.assertNotIn("feature_id", strict)
        self.assertAlmostEqual(math.expm1(strict["velocity_half_window_kms"] / followup.C_KMS)
                               * original[1]["wavelength"], original[1]["match_tolerance_angstrom"])

    def test_missing_or_partial_mask_is_unavailable(self):
        rows = [dict(row, or_mask=0) for row in spectra.synthetic_rows()]
        del rows[-1]["or_mask"]
        features = spectra.analyze(rows)["features"]
        with patch.object(followup, "screen") as rerun:
            result = measure([request(1)], features, rows)[0]
        rerun.assert_not_called()
        self.assertEqual(result["strict_search"], {"outcome": "unavailable", "reason": "missing_or_mask"})

    def test_covered_window_without_same_sign_match_is_a_miss(self):
        rows = [dict(row, or_mask=0) for row in spectra.synthetic_rows()]
        features = spectra.analyze(rows)["features"]
        # Replay the recorded query against a spectrum with a flat flux window.
        for row in rows:
            row["flux"] = 20.0
        result = measure([request(1)], features, rows)[0]["strict_search"]
        self.assertEqual(result["outcome"], "miss")

    def test_original_galaxy_oxygen_has_no_strict_coverage(self):
        example = Path(__file__).resolve().parents[1] / "data/sdss-6497-56329-0139.json"
        meta, rows = load_observation(example)
        z = meta["redshift"]
        features = spectra.analyze(rows, z)["features"]
        index = min(range(len(features)), key=lambda i: abs(features[i]["rest_wavelength"] - 5008.239))
        self.assertLess(abs(features[index]["rest_wavelength"] - 5008.239), 2)
        result = measure([request(index + 1)], features, rows, z)[0]
        self.assertEqual(result["strict_search"]["outcome"], "unsearchable")
