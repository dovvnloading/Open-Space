import csv
import json
import math
import unittest
from pathlib import Path

from open_space.spectra import analyze, benchmark, synthetic_rows


def artificial_line(center=5000.0, amplitude=20.0):
    return [{"wavelength": center + dx, "flux": 0.3 * dx + amplitude * math.exp(-0.5 * dx ** 2),
             "ivar": 1.0, "mask": 0} for dx in range(-30, 31)]


class SpectralScreenTests(unittest.TestCase):
    def test_recovery_and_replay(self):
        report = analyze(synthetic_rows())
        self.assertEqual(report["counts"], {"catalog_match": 1, "sky_match": 1, "unidentified": 1})
        self.assertEqual(report["quality"]["searched_centers"], 784)
        self.assertEqual(report, json.loads(json.dumps(analyze(synthetic_rows()), allow_nan=False)))
        found = {f["classification"]: f for f in report["features"]}
        self.assertAlmostEqual(found["catalog_match"]["wavelength"], 4862.683, delta=1.1)
        self.assertEqual(found["unidentified"]["direction"], "absorption")
        self.assertEqual(found["unidentified"]["matches"], [])

    def test_source_redshift_and_observer_airglow(self):
        z = 0.08
        report = analyze(artificial_line(5008.24 * (1 + z)), redshift=z)
        hit, = report["features"]
        self.assertAlmostEqual(hit["rest_wavelength"], 5008.24, places=3)
        self.assertEqual(hit["matches"][0]["name"], "[O III]")
        sky, = analyze(artificial_line(5578.9), redshift=z)["features"]
        self.assertEqual(sky["classification"], "sky_match")

    def test_sky_and_source_coincidence_keeps_both_alternatives(self):
        hit, = analyze(artificial_line(6302.046))["features"]
        self.assertEqual(hit["classification"], "sky_match")
        self.assertEqual({m["frame"] for m in hit["matches"]}, {"source", "observer"})

    def test_coverage_is_distinct_from_a_null_search(self):
        flat = artificial_line(amplitude=0)
        self.assertEqual(analyze(flat)["quality"]["status"], "searched")
        self.assertEqual(analyze(flat)["features"], [])
        for short in ([], flat[:24], [{**row, "mask": 1} for row in flat]):
            report = analyze(short)
            self.assertEqual(report["quality"]["status"], "no_searchable_centers")
            self.assertEqual(report["quality"]["searched_centers"], 0)
            self.assertEqual(report["features"], [])

    def test_masks_invalid_pixels_and_gaps_break_the_continuum_window(self):
        rows = artificial_line()
        for key, value in (("mask", 4), ("ivar", 0), ("ivar", -1),
                           ("flux", math.nan), ("ivar", math.inf)):
            with self.subTest(key=key, value=value):
                split = list(rows)
                split[30] = {**split[30], key: value}
                result = analyze(split)
                self.assertEqual(result["quality"]["searched_centers"], 12)
                self.assertEqual(result["features"], [])
        # Removing the center leaves a two-Angstrom discontinuity in a one-Angstrom grid.
        result = analyze(rows[:30] + rows[31:])
        self.assertEqual(result["quality"]["searched_centers"], 12)
        self.assertEqual(result["features"], [])

    def test_isolated_spike_is_rejected(self):
        rows = artificial_line(amplitude=0)
        rows[30]["flux"] += 1000
        self.assertEqual(analyze(rows)["features"], [])

    def test_malformed_metadata_rejected(self):
        rows = artificial_line()
        for bad_z in (math.nan, math.inf, -1, -1.5):
            with self.assertRaises(ValueError):
                analyze(rows, redshift=bad_z)
        for change in ({"mask": 0.7}, {"wavelength": math.nan}, {"wavelength": 0},
                       {"wavelength": rows[-2]["wavelength"]}, {"flux": "missing"}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                analyze(rows[:-1] + [{**rows[-1], **change}])

    def test_labelled_benchmark_and_independent_null(self):
        report = benchmark(3)
        self.assertEqual(report["recovered"], {"catalog_match": 3, "sky_match": 3, "unidentified": 3})
        self.assertEqual(report["null_features"], 0)
        self.assertEqual(report["masked_artifact_features"], 0)
        self.assertEqual(report["null_searched_centers"], 3 * 784)
        with self.assertRaises(ValueError):
            benchmark(0)

    def test_pinned_sdss_positive_control(self):
        path = Path(__file__).resolve().parents[1] / "data/sdss-6497-56329-0139.csv"
        with path.open(encoding="utf-8", newline="") as spectrum:
            report = analyze(csv.DictReader(spectrum), redshift=0.08566851168870926)
        self.assertEqual(report["quality"]["searched_centers"], 3289)
        found = {(match["name"], match["reference_angstrom"])
                 for feat in report["features"] for match in feat["matches"]}
        for line in (("H beta", 4862.683), ("[O III]", 4960.295), ("[O III]", 5008.24)):
            self.assertIn(line, found)
        self.assertEqual(report["counts"], {"catalog_match": 12, "sky_match": 0, "unidentified": 0})


if __name__ == "__main__":
    unittest.main()
