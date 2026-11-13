import json
from contextlib import closing
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch

from open_space import discovery, discovery_agents, discovery_queue as queue
from open_space.archive_index import parse_votable
from open_space.ledger import fingerprint, transaction
from open_space.planet_catalog import assess


class ArchiveDiscoveryTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.db_path = self.root / discovery.DB

    def test_catalog_coverage_and_alias(self):
        cat = {"targets": {"1": [], "2": [{"period_days": 2}]}}
        self.assertEqual(assess("3", 1, cat)["status"], "not_queried")
        self.assertEqual(
            assess("1", 1, cat)["status"], "no_exact_tic_entry_in_snapshot"
        )
        self.assertEqual(
            assess("2", 1, cat)["matches"][0]["period_relation"], "half_period"
        )
        self.assertEqual(
            assess("2", 4, cat)["matches"][0]["period_relation"], "double_period"
        )
        with self.assertRaisesRegex(ValueError, "truncated"):
            parse_votable(
                b'<VOTABLE><RESOURCE><INFO name="QUERY_STATUS" value="OVERFLOW"/></RESOURCE></VOTABLE>'
            )

    def test_expired_workers_and_transitive_review(self):
        queue.initialize(self.db_path, {"test": True})
        with closing(queue.connect(self.db_path)) as db:
            pixel = queue.enqueue(db, "pixels", "field", {})
            lease0 = queue.claim(db, 2)
            db.execute("UPDATE task SET expires=0")
            lease1 = queue.claim(db, 2)
            self.assertEqual(lease0["id"], lease1["id"])
            with self.assertRaises(queue.LostLease):
                queue.complete(db, lease0, {}, [])
            reviewer = queue.enqueue(db, "review", "field", {"children": [pixel]})
            queue.complete(db, lease1, {}, [("counterparts", "field", {})])
            follow = queue.claim(db, 2)
            self.assertEqual(follow["kind"], "counterparts")
            self.assertIsNone(queue.claim(db, 2))
            queue.complete(db, follow, {}, [])
            self.assertEqual(queue.claim(db, 2)["id"], reviewer)
            self.assertEqual(
                db.execute(
                    "SELECT COUNT(*) FROM dependency WHERE task_id=?", (reviewer,)
                ).fetchone()[0],
                2,
            )

    def test_model_controls_real_http_and_bounded_parallel_followup(self):
        calls = []

        class Endpoint(BaseHTTPRequestHandler):
            def do_POST(s):
                body = json.loads(s.rfile.read(int(s.headers["Content-Length"])))
                packet = json.loads(body["messages"][1]["content"])
                calls.append(packet)
                answer = {
                    "assessment": "The control needs a repeat.",
                    "falsifier": "Test another epoch.",
                    "actions": packet["available_actions"],
                }
                response = {
                    "choices": [
                        {
                            "finish_reason": "stop",
                            "message": {"content": json.dumps(answer)},
                        }
                    ],
                    "model": body["model"],
                    "usage": {"completion_tokens": 24},
                }
                s.send_response(200)
                s.end_headers()
                s.wfile.write(json.dumps(response).encode())

            def log_message(s, *args):
                pass

        server = ThreadingHTTPServer(("127.0.0.1", 0), Endpoint)
        t = threading.Thread(target=server.serve_forever, daemon=True)
        t.start()
        self.addCleanup(server.server_close)
        self.addCleanup(t.join)
        self.addCleanup(server.shutdown)
        model = {
            "endpoint": f"http://127.0.0.1:{server.server_port}/v1/chat/completions",
            "analyst_model": "local-analyst",
            "reviewer_model": "local-critic",
            "max_tokens": 256,
            "version": discovery_agents.VERSION,
        }
        queue.initialize(self.db_path, {"model": model, "code_sha256": fingerprint()})
        field = {"id": "test"}
        inputs = {
            "field": field,
            "counts": {},
            "candidates": [{"ra_deg": 1, "dec_deg": 2}],
            "sources": [],
        }
        counterpart = {
            "records": [],
            "radius_arcsec": 45,
            "variable_catalog": {"status": "complete"},
        }
        with closing(queue.connect(self.db_path)) as db:
            queue.enqueue(db, "motion", "test", {"field": field})

        def measure(job, _):
            return inputs if job["kind"] == "motion" else counterpart

        with patch.object(discovery, "execute", side_effect=measure), patch.dict(
            "os.environ", {"no_proxy": "127.0.0.1"}
        ):
            one = discovery.work(self.root, {}, 3, 1, 2, 2)
            self.assertEqual(one["tasks"], {"done": 1, "queued": 1})
            result = discovery.work(self.root, {}, 3, 8, 2, 2)
        self.assertEqual(result["tasks"], {"done": 4})
        self.assertEqual(result["model_calls_reserved"], 2)
        self.assertEqual(len(calls), 2)
        self.assertEqual(calls[0]["available_actions"], ["counterpart_1"])
        self.assertEqual(
            {r["kind"] for r in calls[1]["measurements"]},
            {"motion", "plan", "counterparts"},
        )
        with closing(queue.connect(self.db_path)) as db:
            packet = json.loads(
                db.execute("SELECT result FROM task WHERE kind='review'").fetchone()[0]
            )
        self.assertEqual(packet["requested_model"], "local-critic")

    def test_shards_partition_and_resume(self):
        plan = self.root / "plan.json"
        catalog = self.root / "catalog.json"
        plan.write_text(
            json.dumps(
                {
                    "targets": [
                        {"tic_id": i, "products": [{"sector": 1}, {"sector": 2}]}
                        for i in range(1, 16)
                    ]
                }
            )
        )
        catalog.write_text('{"targets":{}}')
        sets = []
        for i in range(3):
            page = self.root / str(i)
            for _ in range(2):
                discovery.initialize(page, [plan], [], catalog, shard=i, shards=3)
            with closing(queue.connect(page / discovery.DB)) as db:
                sets.append({r[0] for r in db.execute("SELECT id FROM task")})
        self.assertEqual(sum(map(len, sets)), 30)
        self.assertEqual(len(set.union(*sets)), 30)

    def test_coherent_stellar_signal_has_predictive_minimum(self):
        from open_space.stellar import search
        import numpy as np

        t = np.arange(0, 26, 0.02) + 1300

        def f(tt):
            return 0.016 * np.cos(2 * np.pi * (tt - 1300) / 1.36 - 0.4)

        measurement = search([(t, f(t)), (t + 29, f(t + 29))])
        self.assertEqual(measurement["status"], "two_sector_coherent_variability")
        self.assertLess(f(measurement["epoch_min_btjd"]), -0.015)
        self.assertAlmostEqual(measurement["period_days"], 1.36, delta=0.01)
