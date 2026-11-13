import csv
import json
import math
from pathlib import Path
import unittest
from tempfile import TemporaryDirectory

from open_space import holdout, spectra
from open_space.dataio import digest


class FrozenHoldoutTests(unittest.TestCase):
    def test_protocol_pins_published_group_list_and_detector(self):
        root = Path(__file__).resolve().parents[1]
        protocol = json.loads((root / "protocols/holdout-v1.json").read_bytes())
        self.assertEqual(protocol["detector_sha256"], digest(Path(spectra.__file__).read_bytes()))
        self.assertEqual(protocol["discovery_groups_csv_sha256"],
                         digest((root / "results/first-visit-groups.csv").read_bytes()))
        self.assertEqual(len(protocol["predictions"]), 9)

    def test_sample_lock_checks_roster_bytes_and_metadata(self):
        root = Path(__file__).resolve().parents[1]
        protocol = root / "protocols/holdout-v1.json"
        roster, lock = root / "data/holdout-targets.csv", root / "protocols/holdout-sample.json"
        with TemporaryDirectory() as tmp:
            new = Path(tmp)
            local_roster = new / roster.name
            local_roster.write_bytes(roster.read_bytes() + b"\n")
            with self.assertRaisesRegex(ValueError, "sample lock"):
                holdout.run(protocol, local_roster, lock, new)
            with roster.open(newline="") as stream:
                first = next(csv.DictReader(stream))
            name = f"sdss-{int(first['plate']):04d}-{first['mjd']}-{int(first['fiber']):04d}.json"
            (new / name).write_text("{}")
            with self.assertRaisesRegex(ValueError, "metadata differs"):
                holdout.run(protocol, roster, lock, new)

    def test_masked_prediction_has_no_vote_even_when_a_detection_is_nearby(self):
        rows = spectra.synthetic_rows()
        result = holdout.window(*holdout.screen(rows, 0), 5488.0, "emission", 150)
        self.assertEqual(result["outcome"], "unsearchable")

    def test_polarity_catalog_and_velocity_window(self):
        data = holdout.screen(spectra.synthetic_rows(), 0)
        positive = holdout.window(*data, 4862.68, "emission", 150)
        wrong_sign = holdout.window(*data, 4862.68, "absorption", 150)
        off_line = holdout.window(*data, 4862.68 * math.exp(1200 / holdout.C_KMS), "emission", 150)
        self.assertEqual(positive["catalog_classification"], "catalog_match")
        self.assertEqual([wrong_sign["outcome"], off_line["outcome"]], ["miss", "miss"])

    def test_strict_mask_and_unknown_alternate_fits(self):
        rows = [dict(r, or_mask=0) for r in spectra.synthetic_rows()]
        rows[63]["or_mask"] = 256
        rows[63 - 16]["mask"] = 1
        protocol = {"predictions": [{"id":"X", "direction":"emission", "vacuum_rest_angstrom":4862.68}],
                    "velocity_half_window_kms":150, "side_control_offsets_kms":[-1200,1200],
                    "minimum_targets":3, "minimum_visits":2}
        measured = holdout.measure_spectrum(protocol, {"redshift": 0}, rows)
        self.assertEqual(measured[0]["primary"]["outcome"], "hit")
        self.assertEqual(measured[0]["strict_mask"]["outcome"], "unsearchable")
        summary = holdout.summarize(protocol, [{"visit": "1-2", "measurements": measured}])[0]
        self.assertEqual(summary["wide_continuum_evaluable_hits"], 0)
        self.assertFalse(summary["meets_recurrence_rule"])
