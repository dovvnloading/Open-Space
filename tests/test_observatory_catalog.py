"""Visitor catalog regressions for grouped findings, source links and new reports."""

import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from observatory.build_catalog import archive_audit, build_records, build_spotlights, source_link
from open_space.verify.runner import run_checks
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
        self.assertEqual(first["provenance"], {"kind": "catalog",
                         "url": f"{self.repo}/blob/{self.tree.commit}/README.md?plain=1#L5"})
        self.assertEqual(first["summary"], "A blend remains possible.")
        self.assertEqual(first["updatedAt"], updates["results/blend/followup.md"])
        self.assertEqual(first["url"], f"{self.repo}/blob/{self.tree.commit}/results/radio/README.md#epochs")
        self.assertEqual(first["links"][1]["label"], "Test")
        self.assertEqual(records[1]["status"], "Null result")
        self.assertEqual(records[2]["summary"], "Calibration data")
        supplemental = records[3]
        self.assertEqual(supplemental["title"], "New search plan")
        self.assertEqual(supplemental["status"], "Research record")
        self.assertEqual(supplemental["provenance"], {"kind": "report", "url": supplemental["url"]})
        self.assertEqual(supplemental["summary"], "No images have been measured.")
        self.assertIsNone(supplemental["updatedAt"])

    def test_directory_and_missing_links(self):
        directory = source_link(self.tree, self.repo, "results/radio/", "Data")
        self.assertIn(f"/tree/{self.tree.commit}/results/radio", directory["url"])
        with self.assertRaisesRegex(ValueError, "Missing evidence"):
            source_link(self.tree, self.repo, "results/unknown/README.md", "Data")
        with self.assertRaisesRegex(ValueError, "Missing evidence anchor"):
            source_link(self.tree, self.repo, "results/radio/README.md#other-epoch", "Data")
        self.assertEqual(source_link(self.tree, self.repo, "results/radio/#epochs", "Test")["url"],
                         f"{self.repo}/tree/{self.tree.commit}/results/radio#epochs")

    def test_inconsistent_catalog_totals_stop_publication(self):
        tree = RepoTree({"README.md": b"The catalog has 2 entries.\n- **One.** *Measured* [A](a.md)\n",
                         "a.md": b"# A"}, commit="2" * 40)
        with self.assertRaisesRegex(ValueError, "states 2 entries, catalog has 1"):
            build_records(tree, self.repo, {})

    def test_scoped_ledger_preserves_checks_and_rejects_wrong_identity_or_coverage(self):
        manifest = self.tree.read("README.md").replace(b"docs/CONTROL.md", b"results/README.md")
        tree = RepoTree({p: self.tree.read(p) for p in self.tree.paths} | {
            "README.md": manifest, "docs/VERIFICATION.md": b"# Method", "results/README.md": b"# Results"
        }, origin="commit", commit=self.tree.commit)
        records = build_records(tree, self.repo, {})
        before = archive_audit(tree, self.repo, records, None)
        self.assertEqual((before["catalogEntries"], before["supplementalRecords"]), (3, 1))
        self.assertEqual((before["recordLinks"], before["uniqueEvidencePaths"]), (5, 4))
        self.assertIsNone(before["repositoryChecks"])
        self.assertEqual(records[0]["evidenceChecks"], [])
        ledger = run_checks(tree)
        measured = archive_audit(tree, self.repo, records, ledger)
        self.assertEqual(measured["repositoryChecks"]["resultFolders"], 3)
        self.assertEqual(measured["repositoryChecks"]["declaredClaims"], 0)
        self.assertEqual([c["folder"] for c in records[0]["evidenceChecks"]], ["blend", "radio"])
        self.assertEqual(records[2]["evidenceChecks"], [])  # a top-level index, outside result folders
        for wrong in [dict(ledger, scope=["radio"]), dict(ledger, checks=["catalog"]),
                      dict(ledger, tree=dict(ledger["tree"], commit="2" * 40)),
                      dict(ledger, tree=dict(ledger["tree"], origin="worktree")),
                      dict(ledger, results=ledger["results"][1:])]:
            with self.subTest(scope=wrong["scope"]), self.assertRaises(ValueError):
                archive_audit(tree, self.repo, records, wrong)

    def test_published_shared_report_ids_survive_qualification_changes(self):
        primary = "results/neowise-nuclear-flares/README.md"
        tree = RepoTree({"README.md": (
            f"- **New title for flares.** *Candidate* [Data]({primary})\n"
            f"- **Scoped fading result.** *Null result* [Sample]({primary}#no-vanishing-supergiant)\n"
        ).encode(), primary: b"# Report\n\n## No vanishing supergiant\n"}, commit="a" * 40)
        records = build_records(tree, self.repo, {})
        self.assertEqual([r["id"] for r in records], ["results-neowise-nuclear-flares-28374e29",
                                                      "results-neowise-nuclear-flares-3dc81494"])

    def test_spotlights_and_rendered_figure_come_from_requested_commit(self):
        revision = self.tree.commit
        raw = b"source table at revision"
        image = b"rendered image with source digest " + hashlib.sha256(raw).hexdigest().encode()
        legacy_metric = {"label": "Survey band", "value": "Radio"}
        curated = [{"id": "radio", "path": "results/radio/README.md", "evidence": [], "measurements": [legacy_metric],
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
            self.assertEqual(result[0]["measurements"], [legacy_metric])
            self.assertEqual(git.call_count, 3)
            current_tree = RepoTree({p: tree.read(p) for p in tree.paths} |
                                    {"observatory/verify_spotlights.py": b""}, commit=revision)
            with self.assertRaisesRegex(ValueError, "each audited target"):
                build_spotlights(current_tree, self.repo, records, Path(tmp) / "current.json")

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
