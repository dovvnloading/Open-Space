import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import campaign
import lead_optical


class LeadCatalogAudit(unittest.TestCase):
    def test_duplicate_trigger_numbers_in_distinct_products_keep_their_receipts(self):
        with tempfile.TemporaryDirectory() as tmp:
            evidence=Path(tmp)
            folder=evidence/"observations"
            folder.mkdir()
            base={"trigger_key":"real/0","ra_deg":11.3,"dec_deg":-23.4,
                  "position_error_arcsec":1.2,"event_mjd":58001.}
            event={**base,"instrument_pass":True,"vet":{"status":"ambiguous",
                   "matches":[{"catalog":"gaia","name":"123"}]}}
            products=[{"input":{"sha256":s*64},"complete":True,"search":{"triggers":[event]}} for s in "ab"]
            campaign.save_obs(folder/"100.json.gz",{"obsid":100,"status":"complete","products":products})
            self.assertEqual(len(list(lead_optical.inputs(evidence,None))),2)
            choice=Path(tmp)/"leads.json"
            row={**base,"obsid":100,"input_sha256":"b"*64}
            choice.write_text(json.dumps([row]))
            self.assertEqual(list(lead_optical.inputs(evidence,choice)),[(row,["123"])])
            choice.write_text(json.dumps([{**row,"ra_deg":11.2}]))
            with self.assertRaisesRegex(ValueError,"disagrees"):
                list(lead_optical.inputs(evidence,choice))

    def test_incomplete_cone_remains_distinct_from_empty_validated_cone(self):
        with patch.object(lead_optical.optical,"ps1_service",side_effect=ValueError("missing columns")):
            self.assertEqual(lead_optical.run_catalog("ps1_mean",{},None,True)["status"],"incomplete")
        with patch.object(lead_optical.optical,"ps1_service"):
            self.assertEqual(lead_optical.run_catalog("ps1_mean",{},None,True)["status"],"empty_validated_cone")
