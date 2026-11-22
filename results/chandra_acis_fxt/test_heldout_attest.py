import json
import unittest

import heldout_attest as att


class HeldoutScheduler(unittest.TestCase):
    def test_a_queued_timestamp_is_not_a_worker_start(self):
        j={"id":10,"name":"observation (96)","status":"queued","runner_id":0,"steps":[],
           "started_at":"2026-11-22T05:02:00Z","run_id":att.ingest.RUN,"head_sha":att.ingest.CODE,
           "run_attempt":att.ingest.ATTEMPT}
        data=json.dumps({"jobs":[j]}).encode()
        jobs={}
        self.assertEqual(att.read_jobs(data,jobs),1)
        self.assertTrue(att.queued(jobs[96]))
        self.assertFalse(att.queued({**jobs[96],"runner_id":12345}))
        self.assertFalse(att.queued({**jobs[96],"step_count":1}))
        with self.assertRaises(ValueError):
            att.read_jobs(data,jobs)
