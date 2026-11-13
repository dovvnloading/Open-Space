import copy
import json
import unittest

from open_space.calibrate import _recover, calibrate


def background():
    return [{"wavelength": 4400.0 + i, "flux": i * 0.01,
             "ivar": 1.0, "mask": int(45 <= i < 55)} for i in range(100)]


class CalibrationTests(unittest.TestCase):
    manifest = {"observation_id": "control", "redshift": 0.0,
                "csv_sha256": "a" * 64, "source_sha256": "b" * 64}

    def test_paired_coverage_replay_and_no_mutation(self):
        rows = background()
        before = copy.deepcopy(rows)
        result = calibrate(self.manifest, rows, trials=3, seed=7)
        self.assertEqual(rows, before)
        self.assertEqual(result, json.loads(json.dumps(calibrate(self.manifest, rows, 3, 7))))
        self.assertEqual(result["data_role"], "calibration")
        self.assertEqual(result["searchable_centers"], 42)
        self.assertEqual(result["injected_spectra"], 18)
        self.assertEqual(result["baseline_features"], 0)
        for position in result["positions"]:
            self.assertTrue(12 <= position["pixel"] <= 32 or 67 <= position["pixel"] <= 87)
            self.assertEqual(len(position["outcomes"]), 6)
        for trial in result["summary"]:
            self.assertEqual(trial["attempted"], trial["assessable"] + trial["baseline_confounded"])
            self.assertEqual(trial["assessable"], trial["missed"] + trial["recovered"])
            if trial["peak_noise_units"] == 10:
                self.assertEqual(trial["recovered"], 3)
            if trial["peak_noise_units"] == 3:
                self.assertEqual(trial["recovered"], 0)

    def test_existing_features_and_wrong_sign(self):
        baseline = [{"wavelength": 5000.0, "direction": "emission"}]
        self.assertEqual(_recover(5000.0, 2.0, "emission", baseline, baseline)[0],
                         "baseline_confounded")
        self.assertEqual(_recover(5000.0, 2.0, "absorption", [], baseline)[0], "missed")
        shifted = [{"wavelength": 5001.5, "direction": "absorption"}]
        self.assertEqual(_recover(5003, 2.0, "absorption", baseline, shifted)[0],
                         "baseline_confounded")

    def test_no_coverage_and_bad_trials(self):
        with self.assertRaisesRegex(ValueError, "no searchable"):
            calibrate(self.manifest, background()[:12], 1)
        with self.assertRaisesRegex(ValueError, "positive integer"):
            calibrate(self.manifest, background(), 0)
