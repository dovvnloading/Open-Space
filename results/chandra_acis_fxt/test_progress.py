"""Offline fixture for partial-coverage reporting and overlapping ACIS streams."""
from copy import deepcopy
import unittest

from progress import LENGTHS, describe


class ProgressTest(unittest.TestCase):
    def test_shared_start_windows_are_unioned_and_partial_rows_are_separate(self):
        event={"obsid":11,"seed_score":20,"start_met_s":100,"stop_met_s":115,"ra_deg":12,"dec_deg":23,
               "position_error_arcsec":2,"instrument_pass":True}
        products=[]
        for n,interval in enumerate(((0.,10.,.5),(5.,15.,1.))):
            products.append({"input":{"filename":f"event-{n}","sha256":str(n)},"complete":True,
                             "search":{"triggers":[deepcopy(event),{"instrument_pass":False}],
                                       "start_support_by_duration":{d:[interval] for d in LENGTHS}}})
        obs={"obsid":11,"products":products,"acquired":[p["input"] for p in products],"status":"followup_incomplete"}
        measured=describe(obs)
        self.assertTrue(measured["scanned_all_products"])
        self.assertFalse(measured["complete_followup"])
        self.assertEqual(measured["eligible_start_s_60"],12.5)
        self.assertEqual(measured["real_rows_included"],4)
        self.assertEqual(measured["detector_pass_rows_included"],2)
        self.assertEqual(measured["distinct_detector_passes_included"],1)
        products[1]["complete"]=False
        obs["status"]="error"
        partial=describe(obs)
        self.assertFalse(partial["scanned_all_products"])
        self.assertEqual(partial["eligible_start_s_60"],0)
        self.assertEqual(partial["real_rows_included"],0)
        self.assertEqual(partial["real_rows_recorded"],4)
