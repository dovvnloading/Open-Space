"""Offline fixtures for the distinct primary-upload recovery producer."""
import hashlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import ingest_primary as ingest
import upload_recovery as recovery


class UploadRecoveryTests(unittest.TestCase):
    def test_failed_upload_requires_completed_correct_science_step(self):
        job = {"id": 1001, "run_id": ingest.RUN, "run_attempt": 1, "head_sha": ingest.CODE,
               "status": "completed", "conclusion": "failure", "name": "observation (72)", "html_url": "example",
               "steps": [{"number": 5, "name": "Search assigned observations and record all outcomes",
                          "conclusion": "success"},
                         {"number": 6, "name": "Run actions/upload-artifact@frozen", "conclusion": "failure"},
                         {"number": 7, "name": "Run actions/upload-artifact@frozen", "conclusion": "success"}]}
        log = "Shard 72 has 34 assigned observations\nObsID 12345 complete: []\nCreateArtifact timed out"
        self.assertEqual(recovery.failed_job(job, 72, 1001, 34, log)["original_worker_status_lines"],
                         [("12345", "complete")])
        for wrong in (dict(job, name="observation (71)"), dict(job, head_sha="changed"),
                      dict(job, steps=[dict(job["steps"][0], conclusion="failure"), *job["steps"][1:]])):
            with self.assertRaises(ValueError):
                recovery.failed_job(wrong, 72, 1001, 34, log)
        with self.assertRaises(ValueError):
            recovery.failed_job(job, 72, 1001, 34, log.replace("CreateArtifact", "Photometry"))

    def test_witnessed_stream_exact_while_other_stream_requires_live_archive(self):
        import archive
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            original = tmp / "original"
            original.write_bytes(b"exact-original-evt2")
            oid = 12345
            old, new, alt = (f"acisf{oid}N001_evt2.fits.gz", f"acisf{oid}N004_evt2.fits.gz",
                             f"acisf{oid}_aN003_evt2.fits.gz")
            rec = {"filename": old, "sha256": hashlib.sha256(original.read_bytes()).hexdigest(),
                   "bytes": original.stat().st_size, "fits_sha256": "x", "member": f"{oid}/{old}"}
            saved = {(oid, recovery.stream_key(old)): (original, rec)}
            info = {"id": 9, "digest": "sha256:fixture"}
            url = f"{archive.ROOTS[0]}/5/{oid}/primary/"
            html = f'<a href="{new}"></a><a href="{alt}"></a>'
            with patch.object(archive, "download", return_value={"source": "fresh"}) as fresh:
                with recovery.exact_retained(oid, saved, info):
                    self.assertEqual(set(archive.product_names(html, oid)), {old, alt})
                    dest = tmp / "retrieved"
                    self.assertEqual(archive.download(url+old, dest)["source_member"], rec["member"])
                    self.assertEqual(dest.read_bytes(), original.read_bytes())
                    self.assertEqual(archive.download(url+alt, dest), {"source": "fresh"})
                    with self.assertRaises(ValueError):
                        archive.download(url+new, dest)
                    with self.assertRaises(ValueError):
                        archive.product_names(f'<a href="{alt}"></a>', oid)
                fresh.assert_called_once()
            self.assertFalse(recovery.witness_match(oid, {"acquired": []}, saved))
            self.assertTrue(recovery.witness_match(oid, {"acquired": [rec]}, saved))

    def test_terminal_and_original_heldout_gate(self):
        import heldout_original
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            row = {"artifact_status": "missing_artifact", "artifact_id": None, "job": {"id": 22, "conclusion": "failure"},
                   "assigned_without_verified_worker_ledger": [12345]}
            coverage = {"groups": {"72": row}}
            source = {"artifacts": {}, "jobs": {"72": {"id": 22}}, "run": {"status": "completed"}}
            with self.assertRaises(ValueError):
                recovery.adoption_state(root, 72, 22, [12345], source, coverage)
            marker = root/"held_out"/heldout_original.RECEIPT
            ingest.save(marker, {"source_run": ingest.RUN})
            recovery.adoption_state(root, 72, 22, [12345], source, coverage)
            for bad in (dict(source, run={"status": "in_progress"}),
                        dict(source, artifacts={"chandra-evidence-72": {"id": 8}})):
                with self.assertRaises(ValueError):
                    recovery.adoption_state(root, 72, 22, [12345], bad, coverage)

    def test_stage_manifest_detects_modified_and_unexpected_members(self):
        with tempfile.TemporaryDirectory() as directory:
            stage = Path(directory)
            ingest.save(stage/"original"/"receipt.json", {"source": 1})
            ingest.save(stage/"seal.json", {"schema": 1, "files": recovery.filemap(stage)})
            recovery.verify_stage(stage)
            ingest.save(stage/"science"/"extra.json", {})
            with self.assertRaises(ValueError):
                recovery.verify_stage(stage)
