import json
from pathlib import Path
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor

from open_space import ledger
from open_space.dataio import receipt_paths
from open_space.runner import run_work


PILOT = Path(__file__).resolve().parents[1] / "data/sdss-6497-56329-0139.json"


class DurableWorkerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "jobs.sqlite"
        self.db = ledger.connect(self.path)
        self.addCleanup(self.db.close)
        self.code = ledger.fingerprint()
        self.config = {"endpoint": "http://localhost/v1/chat/completions", "analyst_model": "test"}

    def add(self, config=None):
        return ledger.enqueue(self.db, PILOT, config, self.code)

    def test_duplicate_enqueue_and_empty_work_are_idempotent(self):
        job = self.add()
        self.assertEqual(job, self.add())
        stats = run_work(self.db, self.code, 2, 0)
        self.assertEqual(stats["completed"], 1)
        report = ledger.evidence(self.db, job)
        self.assertEqual(report["stages"]["detection"]["counts"]["catalog_match"], 12)
        self.assertEqual(report["provider_calls_reserved"], 0)
        self.assertEqual(run_work(self.db, self.code, 1, 0)["completed"], 0)

    def test_checkpoints_preserve_analyst_and_budget(self):
        job = self.add(self.config)
        calls = []

        def provider(stage, scan, previous, config):
            calls.append(stage)
            if stage == "analyst":
                self.assertIsNone(previous)
                return {"requests": [{"feature_id": "F0001", "check": "mask_neighborhood"},
                                     {"feature_id": "F0001", "check": "continuum_sensitivity"}]}
            self.assertEqual(len(previous["followup_measurements"]), 2)
            self.assertEqual(previous["followup_measurements"][1]["fits"][0]["status"], "measured")
            if calls.count("reviewer") == 1:
                raise RuntimeError("endpoint unavailable")
            return {"findings": []}

        self.assertEqual(run_work(self.db, self.code, 3, 0, provider)["paused_for_call_budget"], 1)
        self.assertEqual(calls, [])
        self.assertEqual(run_work(self.db, self.code, 3, 1, provider)["paused_for_call_budget"], 1)
        self.assertEqual(run_work(self.db, self.code, 3, 1, provider)["failed"], 1)
        self.assertEqual(ledger.evidence(self.db, job)["execution_error"], "endpoint unavailable")
        ledger.retry(self.db, job, self.code)
        self.assertEqual(run_work(self.db, self.code, 1, 1, provider)["completed"], 1)
        report = ledger.evidence(self.db, job)
        self.assertEqual(calls, ["analyst", "reviewer", "reviewer"])
        self.assertEqual(report["provider_calls_reserved"], 3)
        self.assertEqual(set(report["stages"]), {"detection", "analyst", "reviewer", "measurements"})

    def test_competing_workers_and_lease_recovery(self):
        job = self.add()

        def attempt(_):
            db = ledger.connect(self.path)
            try:
                return ledger.claim(db, self.code)
            finally:
                db.close()

        with ThreadPoolExecutor(max_workers=3) as pool:
            claimed = [x for x in pool.map(attempt, range(3)) if x]
        self.assertEqual(len(claimed), 1)
        self.db.execute("UPDATE jobs SET lease=0 WHERE id=?", (job,))
        with self.assertRaises(ledger.LostLease):
            ledger.checkpoint(self.db, claimed[0], "analyst", {})
        self.assertEqual(run_work(self.db, self.code, 1, 0)["completed"], 1)
        with self.assertRaises(ledger.LostLease):
            ledger.finish(self.db, claimed[0], "failed", "stale")

    def test_offline_target_never_runs_pending_model_job(self):
        model_job = self.add(self.config)
        pilot = self.add()
        run_work(self.db, self.code, 1, 0, only_job=pilot)
        self.assertEqual(ledger.evidence(self.db, model_job)["stages"], {})
        self.assertEqual(ledger.evidence(self.db, pilot)["state"], "completed")

    def test_changed_revision_is_separate_experiment(self):
        first = self.add()
        second = ledger.enqueue(self.db, PILOT, None, "different")
        self.assertNotEqual(first, second)
        self.assertEqual(run_work(self.db, "different", 2, 0)["completed"], 1)
        self.assertEqual(ledger.evidence(self.db, first)["state"], "queued")

    def test_receipt_rejects_inconsistent_provenance(self):
        manifest = json.loads(PILOT.read_text())
        temp = self.path.parent
        (temp / PILOT.name).write_text(json.dumps(manifest))
        entry = {k: manifest[k] for k in ("csv_sha256", "source_sha256", "csv_file", "redshift", "observation_id")}
        entry["provenance"] = PILOT.name
        receipt = temp / "receipt.json"
        receipt.write_text(json.dumps({"observations": [entry]}))
        self.assertEqual(receipt_paths(receipt), [temp / PILOT.name])
        entry["source_sha256"] = "0" * 64
        receipt.write_text(json.dumps({"observations": [entry]}))
        with self.assertRaisesRegex(ValueError, "source_sha256"):
            receipt_paths(receipt)

    def test_inconsistent_csv_is_never_queued(self):
        manifest = json.loads(PILOT.read_text())
        temp = self.path.parent
        (temp / manifest["csv_file"]).write_text("corrupt")
        path = temp / PILOT.name
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "checksum"):
            ledger.enqueue(self.db, path, None, self.code)
        self.assertEqual(ledger.summary(self.db), [])
