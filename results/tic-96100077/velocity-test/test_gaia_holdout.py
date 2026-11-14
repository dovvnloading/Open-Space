"""Synthetic unit tests for mission timing, release partitioning, flags and velocity sign."""

import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np
from gaia_holdout import DR3_END_OBMT, REQUIRED, TIME_ZERO, compare, read_rows
from orbit import predict_jd


class GaiaHoldoutTest(unittest.TestCase):
    def make_record(self, rev, valid="true", status=0):
        # Encode an artificial focal-plane transit using the published bit convention.
        tid = (int(rev * 21600e9 / 204800) << 17) + 7
        jd = 2457023.75 + (rev - 1717.6256) / 4
        return {
            "source_id": "1234567890123456789",
            "transit_id": str(tid),
            "rv_validity_flag": valid,
            "radial_velocity_error": "1.5",
            "radial_velocity": str(predict_jd(jd)),
            "obs_time_rv": str(jd - TIME_ZERO),
            "radial_velocity_status": str(status),
            "obs_time_bary_corr": "300000000000",
        }

    def parse(self, records, clock):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder) / "sample.csv"
            with p.open("w", newline="") as s:
                writer = csv.DictWriter(
                    s, fieldnames=sorted(REQUIRED | {"obs_time_bary_corr"})
                )
                writer.writeheader()
                writer.writerows(records)
            return read_rows(p, "1234567890123456789", clock)

    def test_partition_and_clock(self):
        data = [
            self.make_record(DR3_END_OBMT - 0.0001),
            self.make_record(DR3_END_OBMT + 0.0001, status=1 << 55),
            self.make_record(7000, valid="false"),
        ]
        rows, counts = self.parse(data, "barycentric")
        shifted, _ = self.parse(data, "spacecraft")
        self.assertEqual(
            counts,
            {
                "archive_rows": 3,
                "overlapping_dr3": 1,
                "scored": 1,
                "invalidated_by_archive": 1,
            },
        )
        self.assertTrue(rows[0][4])
        self.assertAlmostEqual((shifted[0][1] - rows[0][1]) * 86400, 300.0, places=4)

    def test_malformed_and_duplicate(self):
        base = self.make_record(7000)
        for record in (
            dict(base, transit_id="0"),
            dict(base, source_id="9"),
            dict(base, obs_time_rv="58000"),
            dict(base, rv_validity_flag=""),
        ):
            with self.subTest(record=record), self.assertRaises(ValueError):
                self.parse([record], "barycentric")
        with self.assertRaises(ValueError):
            self.parse([base, base], "barycentric")

    def test_injected_orbit_and_constant(self):
        rows, _ = self.parse(
            [self.make_record(r) for r in (6000, 6200, 6400, 6600, 6800)], "barycentric"
        )
        fixed = compare(rows)
        self.assertLess(fixed["weighted_rms_fixed_km_s"], 1e-8)
        self.assertAlmostEqual(fixed["diagnostic_free_slope"], 1.0)
        flat = compare([(t, d, 9.0, e, b) for t, d, v, e, b in rows])
        self.assertLess(flat["chisq_flat_minus_fixed"], 0)
        self.assertTrue(
            np.isfinite(
                flat["quadratic_form_with_joint_orbit_and_measurement_covariance"]
            )
        )
        self.assertEqual(compare([])["n"], 0)
        self.assertNotIn("diagnostic_free_slope", compare(rows[:1]))
