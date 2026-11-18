"""Visitor catalog regressions for grouped findings, source links and new reports."""

import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from observatory.build_catalog import build_records, build_spotlights, source_link
from open_space.verify.tree import RepoTree


class ObservatoryCatalogTests(unittest.TestCase):
    def setUp(self):
        readme = """# Observatory

### Stars

- **A moving source.** A blend remains possible. *Candidate* · [Images](results/radio/README.md#epochs) · [Test](results/blend/followup.md)
- **A second test.** *Null result* · [Data](results/radio/README.md)
- **A control.** *Measured* · [Control study](docs/CONTROL.md)
"""
        self.tree = RepoTree({
            "README.md": readme.encode(),
            "results/radio/README.md": b"# Radio images\n\n## Epochs\n",
            "results/blend/followup.md": b"# Blend check\n",
            "results/plan/SEARCH_PLAN.md": b"# New search plan\n\nNo images have been measured.\n",
            "docs/CONTROL.md": b"# Calibration data\n",
        }, commit="1" * 40)
        self.repo = "https://github.com/example/project"

    def test_grouped_results_keep_status_limits_and_all_evidence(self):
        updates = {"results/radio/README.md": "2026-01-04T13:00:00+03:00",
                   "results/blend/followup.md": "2026-01-04T11:30:00+00:00"}
        records = build_records(self.tree, self.repo, updates)
        self.assertEqual(len(records), 4)
        self.assertEqual(len({r["id"] for r in records}), len(records))
        first = records[0]
        self.assertEqual(first["status"], "Candidate")
        self.assertEqual(first["topic"], "Stars")
        self.assertEqual(first["summary"], "A blend remains possible.")
        self.assertEqual(first["updatedAt"], updates["results/blend/followup.md"])
        self.assertEqual(first["url"], f"{self.repo}/blob/{self.tree.commit}/results/radio/README.md#epochs")
        self.assertEqual(first["links"][1]["label"], "Test")
        self.assertEqual(records[1]["status"], "Null result")
        self.assertEqual(records[2]["summary"], "Calibration data")
        supplemental = records[3]
        self.assertEqual(supplemental["title"], "New search plan")
        self.assertEqual(supplemental["status"], "Research record")
        self.assertEqual(supplemental["summary"], "No images have been measured.")
        self.assertIsNone(supplemental["updatedAt"])

    def test_directory_and_missing_links(self):
        directory = source_link(self.tree, self.repo, "results/radio/", "Data")
        self.assertIn(f"/tree/{self.tree.commit}/results/radio", directory["url"])
        with self.assertRaisesRegex(ValueError, "Missing evidence"):
            source_link(self.tree, self.repo, "results/unknown/README.md", "Data")

    def test_spotlights_and_rendered_figure_come_from_requested_commit(self):
        revision = self.tree.commit
        raw = b"source table at revision"
        image = b"rendered image with source digest " + hashlib.sha256(raw).hexdigest().encode()
        curated = [{"path": "results/radio/README.md", "evidence": [],
                    "figure": {"path": "results/radio/table.csv", "src": "data/figures/radio.svg"}}]
        tree = RepoTree({
            "results/radio/README.md": b"# Radio images",
            "results/radio/table.csv": raw,
        }, commit=revision)
        records = [{"id": "radio", "status": "Candidate", "updatedAt": None,
                    "links": [{"path": "results/radio/README.md"}]}]
        blobs = {
            ("show", f"{revision}:observatory/data/spotlights.json"): json.dumps(curated).encode(),
            ("show", f"{revision}:observatory/data/figures/radio.svg"): image,
            ("show", f"{revision}:results/radio/table.csv"): raw,
        }
        with TemporaryDirectory() as tmp, patch("observatory.build_catalog.git") as git:
            git.side_effect = lambda *arguments: blobs[arguments]
            result = build_spotlights(tree, self.repo, records, Path(tmp) / "catalog.json")
            self.assertEqual((Path(tmp) / "figures/radio.svg").read_bytes(), image)
            self.assertEqual(result[0]["status"], "Candidate")
            self.assertEqual(git.call_count, 3)

    def test_packaged_dossiers_have_figures_and_exact_evidence_links(self):
        site = Path(__file__).resolve().parents[1] / "observatory"
        catalog = json.loads((site / "data/catalog.json").read_text(encoding="utf8"))
        self.assertRegex(catalog["revision"], r"^[0-9a-f]{40}$")
        records = {r["id"]: r for r in catalog["records"]}
        for item in catalog["spotlights"]:
            with self.subTest(target=item["target"]):
                self.assertEqual(item["status"], records[item["recordId"]]["status"])
                self.assertTrue(item["limitation"])
                self.assertTrue(item["measurements"])
                self.assertTrue((site / item["figure"]["src"]).is_file())
                self.assertIn(catalog["revision"], item["figure"]["url"])
                for link in item["links"]:
                    self.assertIn(catalog["revision"], link["url"])


if __name__ == "__main__":
    unittest.main()
