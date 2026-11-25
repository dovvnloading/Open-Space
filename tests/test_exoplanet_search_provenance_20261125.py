"""Offline regression checks for frozen-screen provenance and denominators."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1] / "results/exoplanet-search-20261125"
sys.path.insert(0, str(ROOT))
try:
    import search_provenance as p
finally:
    sys.path.pop(0)


class FrozenScreenGuards(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.screen = self.root / "screen"
        self.screen.mkdir()
        targets = [{"tic": tic, "sector": 108,
                    "filename": f"tess2026250081459-s0108-{tic:016d}-0310-s_lc.fits",
                    "selection": "eligible"} for tic in (17, 29, 42)]
        targets.sort(key=lambda r: p.digest(f"open-space-20261125:{r['tic']}".encode()))
        for rank, row in enumerate(targets):
            row["rank"] = rank
        self.plan = {"sector": 108, "eligible_targets": 3,
                     "inventory_targets": 3, "inventory": {"n": 3}, "targets": targets}
        (self.root / "plan.json").write_text(json.dumps(self.plan))
        with (self.root / "selection.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(targets[0]))
            writer.writeheader()
            writer.writerows(targets)
        _, self.plan_sha = p.load_plan(self.root)
        self.config = "a" * 64

    def record(self, row, status="searched"):
        product = {"filename": row["filename"], "mast_uri": "mast:TESS/product/" + row["filename"],
                   "sha256": "b" * 64, "columns_sha256": "c" * 64, "bytes": 10}
        peaks = [{"peak_rank": rank, "preliminary_pass": rank == 1,
                  "period_days": 2.5, "epoch_btjd": 10., "duration_days": .1,
                  "white_snr": 9.} for rank in range(1, 6)] if status == "searched" else []
        rec = {"tic": row["tic"], "sector": row["sector"], "rank": row["rank"],
               "status": status, "peaks": peaks,
               "provenance": p.binding(row, self.plan_sha, self.config, product)}
        if status == "error":
            rec["error"] = "retrieval failed"
        if status == "unsearchable":
            rec["reason"] = "coverage below threshold"
        return rec

    def save(self, row, rec):
        (self.screen / p.record_name(row)).write_text(json.dumps(rec))

    def census(self):
        return p.census(self.plan, self.plan_sha, self.config, self.screen)

    def test_incomplete_catalog_fails_closed(self):
        index = self.root / "index.json"
        value = {"status": "complete", "required_failures": [], "parse_issues": [],
                 "tic_exclusions": [7, 8], "unique_tic_exclusion_count": 2}
        index.write_text(json.dumps(value))
        self.assertEqual(p.checked_catalog(index), value)
        for change in ({"status": "incomplete"}, {"required_failures": ["toi"]},
                       {"parse_issues": ["bad row"]}, {"tic_exclusions": [7, 7]},
                       {"unique_tic_exclusion_count": 1}):
            index.write_text(json.dumps({**value, **change}))
            with self.subTest(change=change), self.assertRaises(ValueError):
                p.checked_catalog(index)

    def test_frozen_ledger_mismatch_rejected(self):
        path = self.root / "selection.csv"
        path.write_text(path.read_text().replace("eligible", "excluded", 1))
        with self.assertRaises(ValueError):
            p.load_plan(self.root)

    def test_complete_census_preserves_coverage_exclusions(self):
        for row, status in zip(self.plan["targets"], ("searched", "unsearchable", "searched")):
            self.save(row, self.record(row, status))
        summary, ledger, shortlist = self.census()
        self.assertEqual(summary["status"], "complete_screening_census")
        self.assertEqual(summary["counts"]["searched"], 2)
        self.assertEqual(summary["counts"]["unsearchable"], 1)
        self.assertEqual(len(ledger), 3)
        self.assertEqual(len(shortlist), 2)

    def test_failed_or_partial_screen_cannot_enter_shortlist(self):
        first, second, third = self.plan["targets"]
        self.save(first, self.record(first, "error"))
        rec = self.record(second)
        rec.update(status="error", error="interrupted after a peak")
        self.save(second, rec)
        summary, ledger, shortlist = self.census()
        self.assertEqual(summary["status"], "incomplete")
        self.assertEqual(summary["counts"]["error"], 1)
        self.assertEqual(summary["counts"]["invalid"], 1)
        self.assertEqual(summary["counts"]["missing"], 1)
        self.assertEqual(shortlist, [])

    def test_stale_plan_or_configuration_rejected(self):
        row = self.plan["targets"][0]
        for key in ("plan_sha256", "configuration_sha256", "target_sha256"):
            rec = self.record(row)
            rec["provenance"][key] = "f" * 64
            self.save(row, rec)
            summary, _, shortlist = self.census()
            self.assertEqual(summary["counts"]["invalid"], 1)
            self.assertEqual(shortlist, [])

    def test_foreign_record_is_reported_not_ranked(self):
        row = self.plan["targets"][0]
        self.save(row, self.record(row))
        (self.screen / "tic999_s108.json").write_text(json.dumps(self.record(row)))
        summary, _, shortlist = self.census()
        self.assertEqual(summary["foreign_records"], ["tic999_s108.json"])
        self.assertEqual(len(shortlist), 1)
        self.assertEqual(summary["status"], "incomplete")

    def test_changed_product_not_reusable(self):
        row = self.plan["targets"][0]
        rec = self.record(row)
        changed = {**rec["provenance"]["product"], "sha256": "e" * 64}
        with self.assertRaisesRegex(ValueError, "photometry changed"):
            p.check_record(row, rec, self.plan_sha, self.config, changed)

    def test_wrong_target_and_legacy_records_rejected(self):
        first, second, _ = self.plan["targets"]
        rec = self.record(first)
        self.save(second, rec)
        rec.pop("provenance")
        self.save(first, rec)
        summary, _, shortlist = self.census()
        self.assertEqual(summary["counts"]["invalid"], 2)
        self.assertFalse(shortlist)


if __name__ == "__main__":
    unittest.main()
