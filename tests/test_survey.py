"""Exercise survey restart and budget boundaries with a real local worker queue."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from open_space import discovery, survey
from open_space.tess_archive import CONTROL_TIC


class SurveyContinuationTest(unittest.TestCase):
    def patch_attribute(self, obj, attribute, **kwargs):
        patcher = patch.object(obj, attribute, **kwargs)
        self.addCleanup(patcher.stop)
        return patcher.start()

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.base = Path(tmp.name)
        self.root, self.cache = self.base / "survey", self.base / "cache"
        self.stars = {110, 121, 134, 205, 416}
        self.populate_cache(self.cache, self.stars)
        self.calls = []
        self.fetch = self.patch_attribute(survey, "fetch_plan", side_effect=lambda plan, *_: plan)
        self.patch_attribute(survey, "snapshot", side_effect=lambda ids, _: {
            "targets": {str(tic): [] for tic in ids}})
        self.patch_attribute(discovery, "execute", side_effect=self.measure)
        self.patch_attribute(survey.discovery_agents, "settings", return_value={"fixture": True})
        self.provider = self.patch_attribute(survey.discovery_agents, "decide",
                                             side_effect=AssertionError("Zero-call budget exceeded"))

    @staticmethod
    def populate_cache(cache, stars):
        inventory = cache / "inventories"
        inventory.mkdir(parents=True, exist_ok=True)
        for sector in (1, 2):
            lines = []
            for tic in sorted(stars | {CONTROL_TIC}):
                name = f"tess2018206045859-s{sector:04}-{tic:016}-0120-s_lc.fits"
                lines.append(f"curl -o {name} https://archive.invalid/{name}")
            (inventory / f"s{sector}.sh").write_text("\n".join(lines) + "\n")

    def measure(self, job, _paths):
        kind, payload = job["kind"], job["payload"]
        self.calls.append((kind, job["subject"]))
        if kind == "counterparts":
            return {"records": [], "radius_arcsec": 42, "variable_catalog": {}}
        target = payload["target"]
        core = {"tic_id": target["tic_id"], "role": target["role"], "ra_deg": 2, "dec_deg": -5,
                "sources": target["products"]}
        if kind == "transit":
            return {**core, "disposition": "below_threshold", "catalog": {}}
        self.assertEqual(kind, "variability")
        # No epoch has been supplied in this fixture. Its one action is a sky lookup.
        return {**core, "status": "two_sector_coherent_variability", "period_days": .8}

    def advance(self, **kwargs):
        values = dict(root=self.root, search=1, validation=2, start=0, size=2, pages=1,
                      workers=3, jobs=80, model_calls=0, models=True, seed="stable", cache_root=self.cache)
        return survey.run(**{**values, **kwargs})

    def page_tics(self):
        return [[t["tic_id"] for t in json.loads(p.read_bytes())["targets"]]
                for p in sorted((self.root / "pages").glob("*/targets.json"))]

    def test_restart_page_budget_and_no_duplicate_targets(self):
        first = self.advance(jobs=1)
        self.assertEqual(first["status"], "in_progress")
        self.assertEqual(first["completed_pages"], 0)
        self.assertEqual(len(self.page_tics()), 1)
        self.assertEqual(len(self.calls), 1)
        partial_plan = self.page_tics()[0]
        page1 = self.advance()
        self.assertEqual((page1["status"], page1["completed_pages"], page1["next_offset"]),
                         ("page_budget_reached", 1, 2))
        self.assertEqual(self.page_tics()[0], partial_plan)
        self.assertEqual(self.fetch.call_count, 1)
        final = self.advance(pages=8)
        self.assertEqual((final["status"], final["explored_stars"], final["completed_pages"]),
                         ("exhausted", 5, 3))
        pages = self.page_tics()
        self.assertEqual([len(p) for p in pages], [2, 2, 1])
        self.assertEqual({tic for page in pages for tic in page}, self.stars)
        self.assertEqual(self.fetch.call_count, 3)
        self.assertEqual(len(self.calls), len(set(self.calls)))
        calls = list(self.calls)
        self.assertEqual(self.advance(), final)
        self.assertEqual(self.calls, calls)
        self.assertEqual(final["completed_model_calls"], 0)
        self.provider.assert_not_called()
        self.assertEqual(sum(kind == "counterparts" for kind, _ in calls), len(self.stars))
        for page in (self.root / "pages").iterdir():
            registry = json.loads((page / "export" / "registry.json").read_bytes())
            self.assertEqual(registry["counts"]["model_calls_reserved"], 0)
            plans = [t for t in registry["tasks"] if t["kind"] == "plan"]
            self.assertTrue(plans)
            for t in plans:
                decision = json.loads((page / "export" / t["evidence"]).read_bytes())
                self.assertTrue(decision["model_budget_fallback"])

    def test_survey_uses_same_population_after_cache_move(self):
        self.advance()
        other_cache = self.base / "new-cache"
        self.populate_cache(other_cache, {999, 1001, 1002, 1003, 1004, 1005})
        resumed = self.advance(pages=8, cache_root=other_cache)
        self.assertEqual(resumed["status"], "exhausted")
        pages = self.page_tics()
        flattened = [tic for page in pages for tic in page]
        self.assertEqual(sorted(flattened), sorted(self.stars))
        self.assertEqual(resumed["explored_stars"], len(self.stars))

    def test_export_then_checkpoint_interruption_does_not_repeat_work(self):
        publish = survey.atomic_bytes

        def interrupt_progress(path, data):
            if path.name == "progress.json":
                raise OSError("interrupted after page export")
            publish(path, data)

        with patch.object(survey, "atomic_bytes", side_effect=interrupt_progress):
            with self.assertRaisesRegex(OSError, "after page export"):
                self.advance()
        self.assertFalse((self.root / "progress.json").exists())
        self.assertEqual(len(self.page_tics()), 1)
        completed_work = list(self.calls)
        checkpoint = self.advance()
        self.assertEqual(checkpoint["completed_pages"], 1)
        self.assertEqual(checkpoint["explored_stars"], 2)
        self.assertEqual(self.calls, completed_work)
        self.assertEqual(self.fetch.call_count, 1)
