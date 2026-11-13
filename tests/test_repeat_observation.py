import unittest

from examples.repeat_observation import exposure_ids, fixed_window, paired_targets, ROOT


class RepeatObservationTests(unittest.TestCase):
    def test_real_catalog_selects_one_repeat_and_includes_null_matches(self):
        anchors, pairs = paired_targets((ROOT / "data/repeat-catalog.csv").read_bytes())
        self.assertEqual(len(anchors), 12)
        self.assertEqual(pairs, [((6497, 56329, 139), (6012, 56101, 585),
                                 "1237665180057927873", 0.0)])

    def test_distinct_camera_rows_cannot_double_count_the_integration(self):
        header = {"NEXP": 2, "EXPID001": "b1-00000012-00000010-00000011",
                  "EXPID002": "r1-00000012-00000010-00000011"}
        self.assertEqual(exposure_ids(header)["science_ids"], ["00000012"])
        with self.assertRaisesRegex(ValueError, "incomplete"):
            exposure_ids({**header, "NEXP": 3})

    def test_forced_window_flux_sign_and_missing_coverage(self):
        rows = [{"wavelength": float(i), "flux": 8 if abs(i - 100) <= 3 else 10,
                 "ivar": 4, "and_mask": 0, "or_mask": int(abs(i - 100) <= 3)}
                for i in range(75, 126)]
        index = fixed_window(rows, 0, 100, "and_mask")
        self.assertAlmostEqual(index["fractional_depth"], 0.2)
        self.assertEqual(index["pixels_blue_core_red"], [8, 7, 8])
        self.assertEqual(fixed_window(rows, 0, 100, "or_mask")["status"], "insufficient_pixels")
