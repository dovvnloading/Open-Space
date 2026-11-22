"""Fixed secondary boundaries and symmetric handling of incomplete classifications."""
import unittest

from heldout import criteria, fap


class TestHeldOut(unittest.TestCase):
    def event(self):
        return {"instrument_pass":True, "vet":{"status":"no_counterpart","errors":{},
                "position":{"galactic_latitude_deg":-30}}, "trial_duration_s":60, "duration_s":500,
                "before_counts":1, "after_counts":1, "before_s":40, "after_s":80,
                "seed_score":8.4, "ra_deg":20, "dec_deg":20, "start_met_s":100, "stop_met_s":160,
                "position_error_arcsec":1}

    def test_both_equalities_pass_using_trial_duration(self):
        self.assertEqual(criteria(self.event())[::4], ("pass",1))
        t=self.event()
        t["after_counts"]=2
        self.assertEqual(criteria(t)[0],"fail")

    def test_unknown_remains_unknown_unless_another_rule_definitively_fails(self):
        t=self.event()
        t["vet"]["status"]="catalog_incomplete"
        t["vet"]["errors"]={"ned":"unavailable"}
        self.assertEqual(criteria(t)[0],"unknown")
        t["instrument_pass"]=False
        self.assertEqual(criteria(t)[0],"fail")
        t=self.event()
        t["before_s"]=0
        self.assertEqual(criteria(t)[0],"fail")
        t=self.event()
        del t["before_s"]
        self.assertEqual(criteria(t)[0],"unknown")

    def test_tied_maxima_count_and_two_test_floor(self):
        r=fap(9,[9,10]+[None]*22)
        self.assertEqual((r["k"],r["survey_fap"],r["two_test_bonferroni"]),(2,3/25,6/25))
        self.assertEqual(fap(9,[None]*24)["two_test_bonferroni"],.08)
        self.assertEqual(fap(9,[10]*24)["two_test_bonferroni"],1)
