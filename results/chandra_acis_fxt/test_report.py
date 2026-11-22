"""Coverage, control and catalog regression fixtures for the manual Actions run."""
import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import report


def product(name, spans, replicas=(0, 1)):
    return {"complete": True, "input": {"filename": name, "sha256": name},
            "search": {"start_support_by_duration": {"60": spans}, "triggers": [],
                       "controls": [{"replica": n, "triggers": []} for n in replicas],
                       "injections": []}}


class TestSurveyReport(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.base = Path(temp.name)
        (self.base / "inventory").mkdir()
        (self.base / "observations").mkdir()
        self.enterContext(patch.object(report, "EVIDENCE", self.base))
        self.enterContext(patch.object(report, "N_CONTROL", 2))
        (self.base / "inventory/inventory_exceptions.json").write_text("[]")

    def inventory(self, *obsids):
        entries = [{"obsid": oid} for oid in obsids]
        (self.base / "inventory/inventory.json").write_text(json.dumps(entries))

    def ledger(self, oid, status, products, acquired):
        row = {"status": status, "errors": [], "products": products,
               "acquired": [r["input"] for r in acquired], "catalog_cache": "fixture"}
        path = self.base / f"observations/{oid}.json.gz"
        path.write_bytes(gzip.compress(json.dumps(row).encode()))

    def test_only_fully_scanned_inputs_enter_the_union(self):
        self.inventory(100, 200, 300, 400)
        a = product("stream-a", [[0, 10, .5], [0, 10, .5]])
        b = product("stream-b", [[5, 15, .75]])
        missing = product("undone", [[100, 900, 1]])
        # Overlaps in two CCDs and an alternating exposure stream.
        self.ledger(100, "complete", [a, b], [a, b])
        self.ledger(200, "error", [a], [a, missing])
        # A stale status may outlive an interrupted multi-product write.
        self.ledger(300, "complete", [a], [a, missing])
        # ObsID 400 has no worker ledger at all.
        _, _, rows, _, _, support, receipts, _, complete, searched = report.collect()
        self.assertEqual((complete, searched), (1, 1))
        self.assertEqual([r["obsid"] for r in rows if r["scanned_all_products"]], [100])
        self.assertEqual(rows[-1]["status"], "unattempted")
        self.assertAlmostEqual(support["60"], 5 * .5 + 10 * .75)
        self.assertEqual(sum(r["included_in_inference"] for r in receipts), 2)

    def test_missing_or_repeated_scrambles_cannot_supply_a_complete_ensemble(self):
        self.inventory(100)
        for replicas in [(0,), (0, 0, 1)]:
            with self.subTest(replicas=replicas):
                row = product("science", [[0, 10, 1]], replicas)
                self.ledger(100, "complete", [row], [row])
                with self.assertRaisesRegex(ValueError, "control"):
                    report.collect()

    def test_catalog_relabel_is_idempotent(self):
        original = "\n" * 8 + (
            "The catalog has 1 entries, 0 new measurements, 1 candidate discoveries, "
            "0 forecasts, 0 corrections and 0 carefully recorded null results.\n"
            "- **Fixture.** An unrelated lead. *Candidate*\n"
            f"- **Chandra.** Search pending. [Details]({report.TAG})\n"
        )
        root_readme = self.base / "README.md"
        root_readme.write_text(original)
        with patch.object(report, "ROOT", self.base):
            report.update_catalog("One search completed.", "Measured")
            measured = root_readme.read_text()
            report.update_catalog("One search completed.", "Measured")
            self.assertEqual(root_readme.read_text(), measured)
            self.assertIn("2 entries, 1 new measurements, 1 candidate discoveries", measured.splitlines()[8])
            self.assertEqual(measured.count(report.TAG), 1)
            report.update_catalog("No flare survived.", "Null result")
            final = root_readme.read_text()
            self.assertIn("2 entries, 0 new measurements, 1 candidate discoveries", final.splitlines()[8])
            self.assertIn("1 carefully recorded null", final.splitlines()[8])
            self.assertEqual(final.count("*Null result*"), 1)


if __name__ == "__main__":
    unittest.main()
