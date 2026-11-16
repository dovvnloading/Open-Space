#!/usr/bin/env python3
"""Self-contained controls for the PDS Mars gravity-field frame comparison."""
import importlib.util
import math
import pathlib
import unittest

PATH = pathlib.Path(__file__).with_name("analyze.py")
SPEC = importlib.util.spec_from_file_location("mars_gravity_analyze", PATH)
grav = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(grav)


def synthetic_pair(yaw):
    """Degree-dependent shapes with a known pure longitude relabelling."""
    old, new = {}, {}
    for n in range(1, 41):
        for m in range(n + 1):
            value = complex(math.cos(0.3 * n + 0.71 * m) / (1 + n*n),
                            math.sin(0.2*n - 0.81*m) / (1 + n*n) if m else 0.0)
            shifted = value * complex(math.cos(m*yaw), math.sin(m*yaw))
            old[n, m] = (value.real, value.imag, 1e-12, 1e-12)
            new[n, m] = (shifted.real, shifted.imag, 1e-12, 1e-12)
    return old, new


class TestMarsGravityAudit(unittest.TestCase):
    def test_known_synthetic_yaw_and_disjoint_degree_holdouts(self):
        yaw = 6.067 / grav.ARCSEC_PER_RAD
        a, b = synthetic_pair(yaw)
        fitted = grav.fit_yaw(a, b, 2, 10)
        self.assertAlmostEqual(fitted, yaw, delta=1e-13)
        self.assertAlmostEqual(grav.fit_yaw(a, b, 2, 10, True), yaw, delta=1e-13)
        self.assertAlmostEqual(grav.fit_yaw(a, b, 3, 10), yaw, delta=1e-13)
        for lo, hi in [(11, 20), (21, 40)]:
            stats = grav.band_stats(a, b, lo, hi, fitted, 42828.375, 400)
            self.assertEqual(stats["coefficient_rows"], sum(n + 1 for n in range(lo, hi + 1)))
            self.assertGreater(stats["fraction_raw_power_explained_by_rotation"], 0.999999999)
            self.assertLess(stats["aligned_coefficient_l2"], 1e-13)

    def test_unchanged_models_cannot_imply_spin(self):
        a, _ = synthetic_pair(0)
        self.assertAlmostEqual(grav.fit_yaw(a, a, 2, 10), 0.0, delta=1e-14)

    def test_shadr_parser_rejects_duplicate_coefficient(self):
        header = "3396.0,42828.375,0.0004,120,120,1,0,0"
        rows = [f"{n},{m},1e-10,0.0,1e-12,1e-12"
                for n in range(1, 121) for m in range(n + 1)]
        data = (header + "\n" + "\n".join(rows)).encode("ascii")
        parsed = grav.parse_table(data)
        self.assertEqual(len(parsed["coefficients"]), 7380)
        self.assertEqual(parsed["reference_radius_km"], 3396.0)
        with self.assertRaises(AssertionError):
            grav.parse_table(data + b"\n" + rows[0].encode("ascii"))

    def test_checksum_constants_are_valid_sha256_digests(self):
        for model, receipts in grav.EXPECTED_SHA256.items():
            self.assertIn(model, grav.FILES)
            for digest in receipts.values():
                self.assertEqual(len(digest), 64)
                int(digest, 16)


if __name__ == "__main__":
    unittest.main()
