"""Check that separated spectral injections reproduce isolated real-pipeline trials."""
from unittest import TestCase

from laser_search import inject_line, scan_target
from recovery import _matches
from survey_injections import experiment
from test_laser_search import fixture


class BatchTests(TestCase):
    def test_batch_recovers_same_lines_as_one_at_a_time(self):
        records = [fixture(), fixture(7, 1.2, 0.05)]
        experiment_result = experiment(records, (1,7), 31, 31041)
        self.assertTrue(any(t["eligible"] for t in experiment_result["trials"]))
        for t in experiment_result["trials"]:
            if not t["eligible"]:
                self.assertEqual(set(t["reasons"]), {"outside_joint_coverage"})
                continue
            w = t["rest_A"]
            for level, f in enumerate(experiment_result["flux17"]):
                search = scan_target([inject_line(e,w,f) for e in records], (1,7))
                near = [c for c in search.candidates if _matches(c,w)]
                self.assertEqual(t["reasons"][level] == "recovered", bool(near))
                if near:
                    for measured, h in zip(t["fitted"][level], near[0].hits):
                        self.assertAlmostEqual(measured[0],h.observed_A,places=8)
