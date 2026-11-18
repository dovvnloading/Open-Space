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
