import math
import unittest

from examples.residual_search import covered, feature_evidence, group_record, residual_pixels, searchable_waves
from open_space import spectra


class ResidualSearchTests(unittest.TestCase):
    def test_subtraction_replays_and_does_not_modify_source(self):
        rows = spectra.synthetic_rows(0)
        logs = [math.log10(r["wavelength"]) for r in rows]
        residuals = residual_pixels(rows, [r["flux"] for r in rows], logs)
        self.assertEqual(spectra.analyze(residuals)["features"], [])
        self.assertEqual(len(spectra.analyze(rows)["features"]), 3)
        self.assertEqual(searchable_waves(residuals), searchable_waves(rows))
        with self.assertRaisesRegex(ValueError, "grid"):
            residual_pixels(rows, [r["flux"] for r in rows], list(reversed(logs)))
        with self.assertRaises(ValueError):
            residual_pixels(rows, [1], logs)

    def test_two_visits_must_include_two_plates_and_two_nights(self):
        ids = ("a", "b", "c")
        members = [dict(observation_id=i, feature_id="F0001", direction="emission", score=9,
                        wavelength=5100, rest_wavelength=5000) for i in ids]
        group = dict(members=members, direction="emission", target_count=3, median_rest_angstrom=5000.,
                     median_observed_angstrom=5100., median_score=9., observed_span_angstrom=0.,
                     frame_span_kms=0., visits=[{"plate": 2, "mjd": 9}, {"plate": 5, "mjd": 9}])
        coverage = {i: (.02, [5100.]) for i in ids}
        coverage["other"] = (.03, [])
        self.assertFalse(group_record(group, 1, coverage, "rest")["passes_triage"])
        group["visits"][1]["mjd"] += 1
        out = group_record(group, 1, coverage, "rest")
        self.assertTrue(out["passes_triage"])
        self.assertEqual(out["searchable_targets_at_median"], 3)
        self.assertFalse(covered([], 5000.))
        self.assertFalse(covered([6000.], 5000.))

    def test_strict_quality_control_absence_and_measured_sky_are_distinct(self):
        rows = [dict(r, or_mask=r["mask"]) for r in spectra.synthetic_rows()]
        feature = spectra.analyze(rows)["features"][0]
        i = feature["center_pixel"]
        for row in rows[i-1:i+2]:
            row["or_mask"] = 7
        sky = [1.] * len(rows)
        sky[i] = 10.
        result = feature_evidence(feature, rows, sky, ([], []))
        self.assertEqual(result["core_or_masks"], [7, 7, 7])
        self.assertEqual(result["sky_at_center"], 10.)
        self.assertEqual(result["sky_sideband_median"], 1.)
        self.assertFalse(result["strict_or_window_covered"])
        self.assertFalse(result["strict_or_same_sign_match"])


if __name__ == "__main__":
    unittest.main()
