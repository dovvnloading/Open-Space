"""Regression checks on archived Inner Niger Delta hydraulic measurements.

These guard the reported conclusions using checked-in numeric source extracts.
They do not substitute for refetching products, checking gauge datum history,
or independently validating SWOT discharge.
"""
import csv
import unittest
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "results/mission05-swot-rivers/inner-niger"


def records(relative_path):
    with (ROOT / relative_path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


class InnerNigerDeltaEvidenceTests(unittest.TestCase):
    def test_same_track_stage_matches_are_two_and_within_coverage_sensitivity(self):
        pairs = records("analysis/matched-limbs.csv")
        self.assertEqual([int(p["season"]) for p in pairs], [2024, 2025])
        for pair, expected_change in zip(pairs, (4.6361, 3.5493)):
            self.assertEqual(pair["version"], "D")
            self.assertEqual(float(pair["floor"]), .75)
            self.assertEqual(int(pair["track"]), 363)
            self.assertEqual(pair["crid_rise"], pair["crid_fall"])
            self.assertLessEqual(abs(float(pair["mopti_change_m"])), .20)
            self.assertGreater(float(pair["min_frac"]), .75)
            self.assertAlmostEqual(float(pair["dire_change_m"]), expected_change, places=3)
            self.assertLess(abs(float(pair["dire_slope_change_cm_km"])),
                            float(pair["sigma_slope_change_cm_km"]))

    def test_fixed_nodes_corroborate_each_downstream_change(self):
        node_rows = records("analysis/matched-nodes.csv")
        pairs = {int(p["season"]): float(p["dire_change_m"])
                 for p in records("analysis/matched-limbs.csv")}
        for year, change in pairs.items():
            usable = [float(n["change_m"]) for n in node_rows
                      if int(n["season"]) == year and n["usable"] == "True"]
            self.assertGreaterEqual(len(usable), 1)
            for node_change in usable:
                self.assertLess(abs(node_change - change), .03)

    def test_gauges_agree_without_absolute_datum_transfer(self):
        staff = records("external/hydrometry/staff_stage.csv")
        lookup = {(r["station"], r["start"], r["end"]): r for r in staff}
        dates = (
            (2024, "2024-07-21", "2024-07-31", "2025-01-26", "2025-01-26", 4.43, -.26),
            (2025, "2025-07-13", "2025-07-13", "2026-01-19", "2026-01-19", 3.49, .03),
        )
        swot = {int(p["season"]): float(p["dire_change_m"])
                for p in records("analysis/matched-limbs.csv")}
        for year, j0, j1, f0, f1, expected_dire, expected_mopti in dates:
            for site, expected in (("Diré", expected_dire), ("Mopti", expected_mopti)):
                earlier, later = lookup[(site, j0, j1)], lookup[(site, f0, f1)]
                diff = (int(later["stage_cm"]) - int(earlier["stage_cm"])) / 100
                self.assertAlmostEqual(diff, expected, places=2)
                self.assertNotEqual(earlier["source_id"], "")
                self.assertNotEqual(later["source_id"], "")
                self.assertEqual(later["statistic"], "daily_stage")
                if year == 2024:
                    self.assertEqual(earlier["statistic"], "dekadal_mean")
                else:
                    self.assertEqual(earlier["statistic"], "daily_stage")
            self.assertLess(abs(swot[year] - expected_dire), .35)

    def test_quality_screen_counts_and_no_good_reach_flags(self):
        all_rows = records("swot/swot_coverage.csv")
        groups = defaultdict(list)
        for row in all_rows:
            groups[row["sample"]].append(row)
        self.assertEqual(len(groups["all"]), 9)
        self.assertEqual(sum(int(x["n_rows"]) for x in groups["all"]), 543)
        self.assertEqual(sum(int(x["n_reach_q_0"]) for x in groups["all"]), 0)
        self.assertEqual(sum(int(x["n_rows"]) for x in groups["primary"]), 300)
        self.assertEqual(sum(int(x["n_sos_consensus_q"]) for x in groups["primary"]), 193)
        self.assertEqual(sum(int(x["n_rows"]) for x in groups["strict_geolocation"]), 5)

    def test_primary_multiple_comparison_result_is_not_significant(self):
        rows = records("analysis/anniversary-tests.csv")
        primary = [r for r in rows if r["version"] == "2.0"
                   and r["floor"] == "0.8" and r["season"] == "2023"]
        self.assertEqual(len(primary), 63)
        adjusted = [float(r["primary_bonferroni_p"]) for r in primary
                    if r["primary_bonferroni_p"]]
        self.assertTrue(adjusted)
        self.assertTrue(all(p >= .05 for p in adjusted))

    def test_independent_basin_precipitation_is_context_not_flow_validation(self):
        rows = records("external/precipitation/chirps_seasonal.csv")
        self.assertEqual(len(rows), 35)
        peak = max(rows, key=lambda r: float(r["july_to_september_mm"]))
        self.assertEqual(peak["flood_year"], "2024")
        self.assertAlmostEqual(float(peak["july_to_september_mm"]), 871.9, places=1)


if __name__ == "__main__":
    unittest.main()
