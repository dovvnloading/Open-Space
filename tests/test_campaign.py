import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from open_space import campaign, cohort
from open_space.runner import run_work

DATA = Path(__file__).resolve().parents[1] / "data"
PILOT = "sdss-6497-56329-0139"
CODE = "c" * 64


class CampaignTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        fixture = self.root / "ids.csv"
        fixture.write_text("plate,mjd,fiber\n6497,56329,139\n")
        self.original_plan = cohort.plan

        def plan(plate, mjd, cursor, limit, unused, output):
            return self.original_plan(plate, mjd, cursor, limit, fixture, output)

        def ingest(obs, output, expected_hash):
            self.assertEqual(obs, (6497, 56329, 139))
            for ext in ("csv", "json"):
                shutil.copyfile(DATA / f"{PILOT}.{ext}", output / f"{PILOT}.{ext}")
            return output / f"{PILOT}.json"

        plan_mock = patch.object(cohort, "plan", side_effect=plan)
        fetch_mock = patch.object(cohort, "ingest", side_effect=ingest)
        self.planner = plan_mock.start()
        self.ingest = fetch_mock.start()
        self.addCleanup(plan_mock.stop)
        self.addCleanup(fetch_mock.stop)

    def advance(self, config=None, calls=0, code=CODE, limit=3):
        return campaign.run_campaign(self.root, 6497, 56329, limit, 1, calls, config, code)

    def test_completed_page_exports_and_rerun_does_no_work(self):
        result = self.advance()
        self.assertEqual(result["status"], "exhausted")
        self.assertEqual(result["after_fiber"], 139)
        packet = json.loads((self.root / "pages/after-0000/reports" / f"{PILOT}.json").read_bytes())
        self.assertEqual(packet["state"], "completed")
        self.assertEqual(packet["observation"]["observation_id"], PILOT)
        self.assertEqual(packet["provider_calls_reserved"], 0)
        comparison = json.loads((self.root / "recurrence.json").read_bytes())
        self.assertEqual(comparison["distinct_targets"], 1)
        self.assertEqual(comparison["rest_frame_groups"], [])
        self.assertEqual(self.advance()["pages_committed"], 0)
        self.planner.assert_called_once()
        self.ingest.assert_called_once()
        for changes in ({"code": "changed"}, {"config": {"model": "changed"}}, {"limit": 4}):
            with self.assertRaisesRegex(ValueError, "differ"):
                self.advance(**changes)

    def test_budget_pause_resumes_current_page_without_new_archive_or_analyst_calls(self):
        roles = []

        def role(stage, scan, previous, config):
            roles.append(stage)
            return {"requests": [], "summary": "test agent"}

        def work(*args, **kwargs):
            return run_work(*args, **kwargs, role=role)

        with patch.object(campaign, "run_work", side_effect=work):
            result = self.advance(config={"test": True}, calls=1)
            self.assertEqual(result["status"], "paused_for_call_budget")
            self.assertEqual(result["after_fiber"], 0)
            self.assertEqual(result["provider_calls_reserved_this_invocation"], 1)
            result = self.advance(config={"test": True}, calls=1)
            self.assertEqual(result["status"], "exhausted")
            self.assertEqual(result["provider_calls_reserved_this_invocation"], 1)
        self.assertEqual(roles, ["analyst", "reviewer"])
        self.planner.assert_called_once()
        self.ingest.assert_called_once()

    def test_export_interruption_replays_completed_job_and_preserves_cursor(self):
        original = campaign.checkpoint

        def interrupt(path, data):
            if path.parent.name == "reports":
                raise OSError("simulated full disk")
            return original(path, data)

        with patch.object(campaign, "checkpoint", side_effect=interrupt):
            with self.assertRaisesRegex(OSError, "simulated full disk"):
                self.advance()
        self.assertEqual(json.loads((self.root / "state.json").read_bytes())["after_fiber"], 0)
        resumed = self.advance()
        self.assertEqual(resumed["jobs_completed_this_invocation"], 0)
        self.assertEqual(resumed["after_fiber"], 139)
        self.planner.assert_called_once()
        self.ingest.assert_called_once()

    def test_exclusive_lock_prevents_overlapping_intake(self):
        with campaign.coordinator(self.root):
            with self.assertRaisesRegex(RuntimeError, "Another coordinator"):
                self.advance()
        self.planner.assert_not_called()
