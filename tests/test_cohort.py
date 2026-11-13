import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from open_space import cohort

DATA = Path(__file__).resolve().parents[1] / "data"


class CohortTests(unittest.TestCase):
    def test_live_query_is_visit_bounded_and_sorted(self):
        raw = b"#Table1\nplate,mjd,fiber\n6497,56329,2\n6497,56329,1\n"
        with tempfile.TemporaryDirectory() as directory, patch.object(
            cohort, "urlopen", return_value=io.BytesIO(raw)
        ) as remote:
            path = Path(directory) / "plan.json"
            doc = cohort.plan(6497, 56329, 0, 3, None, path)
            self.assertEqual([t["fiber"] for t in doc["targets"]], [1, 2])
            for part in ["TOP 3", "s.plate = 6497", "s.mjd = 56329", "ORDER BY s.fiberID"]:
                self.assertIn(part, doc["query"])
            self.assertIn("s.sciencePrimary = 1", doc["query"])
            self.assertEqual(json.loads(path.read_bytes()), doc)
            self.assertTrue(remote.call_args.args[0].startswith(cohort.CATALOG + "?"))

    def test_no_data_and_bad_payloads(self):
        header = b"plate,mjd,fiber\n"
        self.assertEqual(cohort.target_list(header, 6497, 56329, 0, 2), [])
        for raw in [b"Server error", header + b"42,56329,1\n", header + b"6497,56329,1\n" * 2,
                    header + b"6497,56329,1\n6497,56329,2\n6497,56329,3\n"]:
            with self.assertRaises(ValueError):
                cohort.target_list(raw, 6497, 56329, 0, 2)
        with self.assertRaises(ValueError):
            cohort.selection_query(6497, 56329, 0, 101)

    def test_operator_plan_and_replay_receipt(self):
        pinned = DATA / "cohort-6497-56329.csv"
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            source = json.loads((DATA / "sdss-6497-56329-0139.json").read_bytes())
            manifest = base / "plan.json"
            plan = cohort.plan(6497, 56329, 0, 3, pinned, manifest)
            self.assertEqual(plan["origin"], "operator_csv")
            self.assertEqual([t["fiber"] for t in plan["targets"]], [88, 96, 100])

            def ingest(obs, output, expected_hash):
                key, _ = cohort.archive_source(*obs)
                name = f"sdss-{key}"
                m = dict(source, observation_id=name, csv_file=f"{name}.csv", source_sha256=expected_hash)
                path = output / f"{name}.json"
                path.write_text(json.dumps(m))
                return path

            with patch.object(cohort, "ingest", side_effect=ingest) as load:
                out = base / "out"
                receipt = cohort.fetch(manifest, out)
                published = json.loads(receipt.read_bytes())
                self.assertEqual(len(published["observations"]), 3)
                self.assertEqual(published["observations"][0]["provenance"], "sdss-6497-56329-0088.json")
                self.assertEqual(cohort.fetch(manifest, out), receipt)
                self.assertEqual(load.call_count, 6)

    def test_plan_fails_before_download_when_target_is_out_of_order(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            manifest = base / "plan.json"
            doc = cohort.plan(6497, 56329, 0, 3, DATA / "cohort-6497-56329.csv", manifest)
            doc["targets"].reverse()
            manifest.write_text(json.dumps(doc))
            with patch.object(cohort, "ingest") as fetch:
                with self.assertRaisesRegex(ValueError, "ordered unique"):
                    cohort.fetch(manifest, base)
                fetch.assert_not_called()
