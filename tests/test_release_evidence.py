"""Check archived task lineage, all hash-pinned members, and the portable CSV index."""

import csv
from hashlib import sha256
import json
from pathlib import Path
import unittest
from zipfile import ZipFile

from examples.image_survey import plan_grid

ROOT = Path(__file__).resolve().parents[1]


class PublishedEvidenceTest(unittest.TestCase):
    def test_export_is_self_contained(self):
        results = ROOT / "results"
        manifest = json.loads((results / "discovery-run.json").read_bytes())
        aggregate = results / "discovery-evidence.zip"
        listing = results / "discovery-registry.csv"
        self.assertEqual(
            sha256(aggregate.read_bytes()).hexdigest(), manifest["bundle_sha256"]
        )
        self.assertEqual(
            sha256(listing.read_bytes()).hexdigest(), manifest["csv_sha256"]
        )
        with ZipFile(aggregate) as bundle:
            names = set(bundle.namelist())
            with listing.open(encoding="utf-8", newline="") as stream:
                public = list(csv.DictReader(stream))
            self.assertTrue(all(row["evidence"] in names for row in public))
            for campaign in manifest["campaigns"]:
                name = campaign["name"]
                raw = bundle.read(name + "/registry.json")
                self.assertEqual(
                    sha256(raw).hexdigest(), campaign["registry_json_sha256"]
                )
                data = json.loads(raw)
                ids = {task["id"] for task in data["tasks"]}
                self.assertEqual(len(ids), len(data["tasks"]))
                self.assertEqual(
                    campaign["measurement_rows"],
                    len([r for r in public if r["campaign"] == name]),
                )
                project = data["campaign"]
                for digest in (
                    project["tess_plan_sha256"]
                    + project["ffi_plan_sha256"]
                    + project["wise_plan_sha256"]
                    + [project["catalog_sha256"]]
                ):
                    self.assertEqual(
                        sha256(bundle.read(f"{name}/inputs/{digest}.json")).hexdigest(),
                        digest,
                    )
                for task in data["tasks"]:
                    self.assertTrue(
                        ids.issuperset(task["parents"] + task["prerequisites"])
                    )
                    self.assertEqual(
                        sha256(bundle.read(name + "/" + task["evidence"])).hexdigest(),
                        task["evidence_sha256"],
                    )
                self.assertTrue(
                    ids.issuperset(row["task_id"] for row in data["findings"])
                )

    def test_image_grid_pages_cover_complete_selection(self):
        def plan(start, limit):
            return plan_grid(359.8, -69.4, 3, 2, 0.4, 0.2, [1, 2], 21, start, limit)[
                "fields"
            ]

        grid = plan(0, 6)
        self.assertEqual(grid, plan(0, 4) + plan(4, 4))
        self.assertEqual(len({row["field_id"] for row in grid}), 6)
        self.assertEqual(grid[1]["ra_deg"], 0.2)
        self.assertEqual(grid[3]["dec_deg"], -69.2)
        with self.assertRaises(ValueError):
            plan(6, 1)
