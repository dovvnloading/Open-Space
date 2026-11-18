"""Displayed observing readouts must survive changes in ordering, epochs and source data."""

import csv
import io
import json
from pathlib import Path
import unittest

from observatory.verify_spotlights import verify_spotlights


ROOT = Path(__file__).resolve().parents[1]
SPOTLIGHTS = ROOT / "observatory/data/spotlights.json"


def read_blob(path):
    return (ROOT / path).read_bytes()


class SpotlightReadoutsTests(unittest.TestCase):
    def test_all_printed_values_have_checked_source_data(self):
        items = json.loads(SPOTLIGHTS.read_bytes())
        verify_spotlights(items, read_blob)
        self.assertEqual(sum(len(item["measurements"]) for item in items), 12)
        for item in items:
            for metric in item["measurements"]:
                with self.subTest(target=item["id"], metric=metric["key"]):
                    self.assertTrue(metric["method"])
                    self.assertTrue(metric["sourcePaths"])
                    for path in metric["sourcePaths"]:
                        self.assertTrue((ROOT / path).is_file())

    def test_mislabeled_dossier_set_cannot_skip_the_audit(self):
        items = json.loads(SPOTLIGHTS.read_bytes())
        with self.assertRaisesRegex(ValueError, "each audited target"):
            verify_spotlights(items[:-1] + [items[0]], read_blob)

    def test_visible_tables_preserve_observing_units_and_nondetections(self):
        items = json.loads(SPOTLIGHTS.read_bytes())
        verify_spotlights(items, read_blob)
        tables = {item["id"]: item["table"] for item in items}
        self.assertEqual({key: len(value["rows"]) for key, value in tables.items()},
                         {"mars": 31, "sun": 30, "quasar": 76, "dwarf": 4})
        for table in tables.values():
            self.assertTrue((ROOT / table["sourcePath"]).is_file())
            self.assertTrue(all(len(row) == len(table["columns"]) for row in table["rows"]))
        peak_day = next(row for row in tables["mars"]["rows"] if row[0] == "2025-08-24")
        self.assertEqual([float(v) for v in peak_day[-2:]], [157000, 616500])
        self.assertEqual(max(tables["sun"]["rows"], key=lambda r: float(r[1])), ["18:27:38.220", "286.350"])
        gaia = tables["quasar"]["rows"]
        peaks = [min((row for row in gaia if row[0] == image), key=lambda row: float(row[2]))
                 for image in ("A", "B")]
        self.assertAlmostEqual(float(peaks[0][1]) - float(peaks[1][1]), 202.90, places=2)
        self.assertEqual(tables["dwarf"]["rows"][1], ["2020-08-25", "-0.150", "0.148", "No centroid recovered"])

    def test_mars_peak_must_match_in_both_deposits(self):
        def altered(path):
            if path.endswith("/maven-2025-sep-events/summary.json"):
                data = json.loads(read_blob(path))
                data["ranked_events"][0]["peak_hourly_median_ion"] *= 2
                return json.dumps(data).encode()
            return read_blob(path)
        with self.assertRaisesRegex(ValueError, "peak disagrees"):
            verify_spotlights(json.loads(SPOTLIGHTS.read_bytes()), altered)

    def test_delayed_gaia_image_changes_the_displayed_measurement(self):
        def delayed(path):
            if not path.endswith("/j0826-lens-delay/gaia-nights.csv"):
                return read_blob(path)
            rows = list(csv.DictReader(io.StringIO(read_blob(path).decode())))
            time = "bjd_tcb_minus_2400000p5"
            for row in rows:
                if row["image"] == "A":
                    row[time] = str(float(row[time]) + 50)
            buf = io.StringIO()
            writer = csv.DictWriter(buf, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
            return buf.getvalue().encode()
        with self.assertRaisesRegex(ValueError, r"quasar\.gaia: displayed"):
            verify_spotlights(json.loads(SPOTLIGHTS.read_bytes()), delayed)


if __name__ == "__main__":
    unittest.main()
