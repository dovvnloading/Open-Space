"""Process and API contract of the JSON evidence report."""
import io
import json
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch

import agent_check
import verify


class AgentCheckTests(unittest.TestCase):
    def invoke(self, *args):
        result = subprocess.run([sys.executable, str(verify.HERE / "agent_check.py"), *args],
                                text=True, capture_output=True)
        return result, json.loads(result.stdout)

    def test_bad_arguments_return_json_on_stdout_and_exit_two(self):
        result, message = self.invoke("--unknown")
        self.assertEqual(result.returncode, 2)
        self.assertEqual(message["status"], "error")
        self.assertEqual(message["error"]["type"], "ArgumentError")
        self.assertNotIn("coverage", message)

    def test_wrong_catalogue_returns_error_and_no_scientific_result(self):
        with tempfile.TemporaryDirectory() as location:
            bad = Path(location) / "table.html"
            bad.write_bytes(b"revised catalogue")
            result, message = self.invoke("--catalog-html", str(bad))
        self.assertEqual(result.returncode, 1)
        self.assertIn("differs from the pinned input", message["error"]["message"])
        self.assertNotIn("external", message)

    def test_requested_bad_proof_fails_instead_of_falling_back(self):
        with tempfile.TemporaryDirectory() as location:
            bad = Path(location) / "proof.pdf"
            bad.write_bytes(b"incorrect proof")
            result, message = self.invoke("--paper-pdf", str(bad))
        self.assertEqual(result.returncode, 1)
        self.assertIn("pinned PDF proof", message["error"]["message"])
        self.assertNotIn("parkes_proof", message)

    def test_success_reports_values_returned_by_replay_and_separates_logs(self):
        def replay_without_pdf(*args):
            print("progress")
            return {"coverage": {"eligible_intervals": 81},
                    "parkes_proof": {"status": "transcription_only", "verified_pdf_sha256": None}}
        out, err = io.StringIO(), io.StringIO()
        with patch.object(agent_check, "replay", replay_without_pdf), redirect_stdout(out), redirect_stderr(err):
            self.assertEqual(agent_check.main([]), 0)
        report = json.loads(out.getvalue())
        self.assertEqual(report["coverage"]["eligible_intervals"], 81)
        self.assertEqual(report["parkes_proof"]["status"], "transcription_only")
        self.assertIsNone(report["parkes_proof"]["verified_pdf_sha256"])
        self.assertEqual(err.getvalue().strip(), "progress")

    def test_nonfinite_numbers_fail_without_partial_json(self):
        with patch.object(agent_check, "replay", return_value={"external": {"skill": float("nan")}}):
            out = io.StringIO()
            with redirect_stdout(out):
                self.assertEqual(agent_check.main([]), 1)
        message = json.loads(out.getvalue())
        self.assertEqual(message["status"], "error")
        self.assertNotIn("external", message)


if __name__ == "__main__":
    unittest.main()
