"""Offline regressions for the frozen, candidate NEA linkage evidence.

These check published inputs and numerics, not a new orbit fit. Reproducing
Find_Orb and ASSIST trajectories requires the pinned external dependencies
documented under results/nea-pair-search.
"""
import csv
import gzip
import hashlib
import json
import math
import unittest
from collections import Counter
from pathlib import Path


PACKET = Path(__file__).resolve().parents[1] / "results" / "nea-pair-search"


def records(path):
    with (PACKET / path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def frozen_json(path):
    return json.loads((PACKET / path).read_text(encoding="utf-8"))


class NearEarthPairEvidenceTests(unittest.TestCase):
    def test_snapshot_denominator_and_candidate_ranking(self):
        source = frozen_json("screen.json")
        frozen = frozen_json("population.json")
        ranks = records("ranking.csv")
        n = source["objects"]
        self.assertEqual(n, 43098)
        self.assertEqual(source["possible_pairs"], n * (n - 1) // 2)
        self.assertEqual(source["pairs_below_threshold"], len(ranks))
        self.assertEqual(len(ranks), 730)
        self.assertEqual([int(row["rank"]) for row in ranks], list(range(1, 731)))
        self.assertEqual([float(row["score"]) for row in ranks],
                         sorted(float(row["score"]) for row in ranks))
        self.assertEqual((ranks[5]["first"], ranks[5]["second"]),
                         ("2007 VE138", "2026 VV15"))
        self.assertEqual((ranks[6]["first"], ranks[6]["second"]),
                         ("2017 RB15", "2026 RG39"))
        self.assertEqual(frozen["objects"], n)
        with (PACKET / "population.csv.gz").open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else None
        if digest is None:
            digest = hashlib.sha256((PACKET / "population.csv.gz").read_bytes()).hexdigest()
        self.assertEqual(digest, frozen["input_sha256"])
        with gzip.open(PACKET / "population.csv.gz", "rt", newline="") as handle:
            self.assertEqual(sum(1 for _ in csv.DictReader(handle)), n)

    def test_optical_counts_and_cfht_out_of_fit_positions(self):
        astrometry = records("orbits/astrometry.csv")
        counts = Counter(row["designation"] for row in astrometry)
        self.assertEqual(len(astrometry), 260)
        self.assertEqual(len({row["obsid"] for row in astrometry}), 260)
        self.assertEqual(counts, {
            "2007 VE138": 43, "2026 VV15": 67,
            "2017 RB15": 19, "2026 RG39": 43,
            "2026 SA12": 28, "2026 SW11": 20,
            "2026 VO": 21, "2026 VR1": 19,
        })
        for a, b, sites in [
            ("2007 VE138", "2026 VV15", 25),
            ("2017 RB15", "2026 RG39", 13),
        ]:
            self.assertEqual(len({r["site"] for r in astrometry
                                  if r["designation"] in (a, b)}), sites)
        original_times = {row["utc"] for row in astrometry}
        cfht = records("cfht/measurements.csv")
        self.assertEqual(len(cfht), 3)
        self.assertEqual(len({row["utc_mid"] for row in cfht}), 3)
        for row in cfht:
            self.assertNotIn(row["utc_mid"], original_times)
            self.assertTrue(row["utc_mid"].startswith("2017-09-15"))
            self.assertAlmostEqual(
                math.hypot(float(row["east_arcsec"]), float(row["north_arcsec"])),
                float(row["separation_arcsec"]), delta=0.001)

    def test_outward_fit_and_unfavorable_temporal_holdout(self):
        force = frozen_json("pressure/fits/VE138.json")
        fit = force["all/A12"]
        self.assertEqual(fit["training"]["n"], 110)
        self.assertGreater(fit["state"][6], 0)
        self.assertAlmostEqual(fit["state"][6] / 1e-8, 1.424, delta=0.001)
        self.assertAlmostEqual(
            math.sqrt(fit["covariance"]["site_date"][6][6]) / 1e-8,
            0.094, delta=0.002)
        self.assertLess(fit["training"]["radial_rms_arcsec"],
                        force["all/grav"]["training"]["radial_rms_arcsec"])
        residual = records("pressure/fits/VE138-residuals.csv")
        self.assertEqual(len(residual), 110)
        reconstructed = math.sqrt(sum(float(r["ra_cosdec_arcsec"]) ** 2
                                      + float(r["dec_arcsec"]) ** 2
                                      for r in residual) / len(residual))
        self.assertAlmostEqual(reconstructed, fit["training"]["radial_rms_arcsec"],
                               delta=0.00001)
        gravity = force["first_new8/grav"]["withheld"]
        acceleration = force["first_new8/A12"]["withheld"]
        self.assertEqual(gravity["n"], 59)
        self.assertEqual(acceleration["n"], 59)
        self.assertAlmostEqual(gravity["radial_rms_arcsec"], 2.506, delta=0.001)
        self.assertAlmostEqual(acceleration["radial_rms_arcsec"], 3.103, delta=0.001)
        self.assertLess(gravity["radial_rms_arcsec"],
                        acceleration["radial_rms_arcsec"])
        self.assertEqual(force["sites/A12"]["withheld"]["n"], 14)

    def test_distinct_implementations_not_distinct_astrometry(self):
        pressure = frozen_json("pressure/fits/VE138.json")
        independent = frozen_json("independent/independent.json")
        corrected = frozen_json("independent/corrected.json")
        primary_a1 = pressure["all/A1"]["state"][6]
        reproduced_a1 = independent["radial"]["heliocentric_ecliptic_j2000_state_and_A1"][-1]
        self.assertEqual(independent["radial"]["positions"], 110)
        self.assertLess(abs(primary_a1 - reproduced_a1) / primary_a1, 0.0001)
        self.assertEqual(corrected["radial"]["positions"], 110)
        self.assertEqual(len(corrected["catalog_corrections"]), 53)
        self.assertAlmostEqual(
            corrected["radial"]["heliocentric_ecliptic_j2000_state_and_A1"][-1]
            / 1e-8, 1.442, delta=0.001)
        result_readme = (PACKET / "README.md").read_text(encoding="utf-8")
        root_readme = (PACKET.parents[1] / "README.md").read_text(encoding="utf-8")
        for expected in ("59 positions", "2.506", "3.103",
                         "share the same optical data and DE440"):
            self.assertIn(expected, result_readme)
        self.assertIn("same optical observations", root_readme)
        self.assertIn("later 59-position holdout favors gravity", root_readme)


if __name__ == "__main__":
    unittest.main()
