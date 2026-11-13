import json
import os
import threading
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from unittest.mock import patch

from open_space import review
from open_space.spectra import analyze, synthetic_rows


class ReviewTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.analysis = analyze(synthetic_rows())

    def setUp(self):
        env = patch.dict(os.environ, {"OPEN_SPACE_MODEL_BASE_URL": "https://example.org/v1",
                                      "OPEN_SPACE_ANALYST_MODEL": "test/analyst",
                                      "OPEN_SPACE_REVIEWER_MODEL": "test/critic"}, clear=True)
        env.start()
        self.addCleanup(env.stop)
        self.config = review.config_from_env()
        self.analyst = {"summary": "Check the local baseline.", "hypotheses": [
            {"feature_id": "F0002", "explanation": "A local baseline error could mimic absorption.",
             "prediction": "A different sideband fit will reduce the signed contrast.",
             "falsifier": "The contrast persists under independent continuum fits."}],
            "requests": [{"feature_id": "F0002", "check": "continuum_sensitivity"}]}
        self.reviewer = {"summary": "One coadd cannot exclude a detector effect.", "findings": [
            {"feature_id": "F0002", "assessment": "inconclusive",
             "challenge": "The follow-up reuses the original exposure.",
             "next_test": "Measure the line in an independently acquired spectrum."}]}

    def test_http_roundtrip_and_audit(self):
        calls = []
        scripts = iter((self.analyst, self.reviewer))

        class Server(BaseHTTPRequestHandler):
            def do_POST(self):
                request = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
                calls.append((self.path, self.headers["Authorization"], request))
                result = {"id": "local-id", "model": request["model"], "usage": {"completion_tokens": 80},
                          "choices": [{"finish_reason": "stop",
                                       "message": {"content": json.dumps(next(scripts))}}]}
                self.send_response(200)
                self.end_headers()
                self.wfile.write(json.dumps(result).encode("utf-8"))

            def log_message(self, *args):
                pass

        server = ThreadingHTTPServer(("127.0.0.1", 0), Server)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        self.addCleanup(server.server_close)
        self.addCleanup(thread.join)
        self.addCleanup(server.shutdown)
        self.config["endpoint"] = f"http://127.0.0.1:{server.server_port}/v1/chat/completions"
        with patch.dict(os.environ, {"OPEN_SPACE_API_KEY": "private-placeholder",
                                     "no_proxy": "127.0.0.1"}):
            sample = {**self.analysis, "source_context": {"kind": "synthetic", "observation_id": "fixture-0"}}
            analyst = review.run_role("analyst", sample, None, self.config)
            checks = [{"feature_id": "F0002", "check": "continuum_sensitivity", "contrast": -8.5}]
            reviewer = review.run_role("reviewer", sample, {
                "analyst": analyst, "followup_measurements": checks}, self.config)
        self.assertEqual(analyst["requests"], self.analyst["requests"])
        self.assertEqual(reviewer["findings"], self.reviewer["findings"])
        self.assertEqual(len(calls), 2)
        for call, result in zip(calls, (analyst, reviewer)):
            path, auth, payload = call
            self.assertEqual(path, "/v1/chat/completions")
            self.assertEqual(auth, "Bearer private-placeholder")
            self.assertNotIn("private-placeholder", json.dumps(result))
            self.assertNotIn("private-placeholder", json.dumps(self.config))
            self.assertEqual(payload["max_completion_tokens"], 1600)
            self.assertEqual(result["review"]["messages"], payload["messages"])
            self.assertEqual(result["review"]["request_sha256"], review._digest(payload))
        packet = json.loads(calls[1][2]["messages"][1]["content"])
        self.assertEqual(packet["followup_measurements"], checks)
        self.assertEqual(packet["source_context"]["kind"], "synthetic")
        self.assertNotIn("review", packet["analyst"])

    def test_unknown_id_and_unimplemented_check_rejected(self):
        for entry in ({"feature_id": "F9999", "check": "mask_neighborhood"},
                      {"feature_id": "F0001", "check": "web_search"}):
            with self.subTest(entry=entry), self.assertRaises(ValueError):
                review._validate("analyst", {**self.analyst, "requests": [entry]}, {"F0001", "F0002"})

    def test_packet_truncation_preserves_original_feature_ids(self):
        source = self.analysis["features"][0]
        sample = [{**source, "score": float(i), "classification": "unidentified"} for i in range(30)]
        packet = review._evidence({**self.analysis, "features": sample}, None)
        self.assertEqual(packet["omitted_features"], 6)
        ids = [row["feature_id"] for row in packet["features"]]
        self.assertEqual(ids, [f"F{i:04d}" for i in range(30, 6, -1)])
        with self.assertRaises(ValueError):
            review._validate("analyst", self.analyst, set(ids))

    def test_malformed_truncated_and_oversize_output_fail_after_one_call(self):
        for content, reason in (('{"summary":"x","summary":"y"}', "stop"),
                                (json.dumps(self.analyst), "length"),
                                (json.dumps({**self.analyst, "summary": "x" * 701}), "stop"),
                                (json.dumps({**self.analyst, "requests": self.analyst["requests"] * 2}), "stop")):
            result = {"choices": [{"finish_reason": reason, "message": {"content": content}}]}
            with self.subTest(content=content), patch.object(review, "_post", return_value=result) as post:
                with self.assertRaises(ValueError):
                    review.run_role("analyst", self.analysis, None, self.config)
                post.assert_called_once()

    def test_no_endpoint_or_invalid_token_budget(self):
        with patch.dict(os.environ, {"OPEN_SPACE_MAX_TOKENS": "9000"}), self.assertRaises(ValueError):
            review.config_from_env()
        with patch.dict(os.environ, {}, clear=True), self.assertRaises(ValueError):
            review.config_from_env()

    def test_api_base_rejects_embedded_credentials(self):
        for base in ("https://user:secret@example.org/v1", "https://example.org/v1?token=secret",
                     "https://example.org/v1#secret"):
            with self.subTest(base=base), patch.dict(os.environ, {"OPEN_SPACE_MODEL_BASE_URL": base}):
                with self.assertRaises(ValueError):
                    review.config_from_env()

    def test_reported_usage_cannot_exceed_request(self):
        response = {"choices": [{"finish_reason": "stop", "message": {"content": json.dumps(self.analyst)}}],
                    "usage": {"completion_tokens": self.config["max_tokens"] + 1}}
        with patch.object(review, "_post", return_value=response), self.assertRaises(ValueError):
            review.run_role("analyst", self.analysis, None, self.config)

    def test_empty_evidence_is_valid_for_abstention(self):
        empty = {**self.analysis, "features": []}
        content = json.dumps({"summary": "No features are available.", "hypotheses": [], "requests": []})
        with patch.object(review, "_post", return_value={"choices": [
            {"finish_reason": "stop", "message": {"content": content}}]}):
            self.assertEqual(review.run_role("analyst", empty, None, self.config)["requests"], [])


if __name__ == "__main__":
    unittest.main()
