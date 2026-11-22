"""The exploratory follow-up must select the original instrument-passing real event."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import field_followup as field
import ingest_primary as ingest


class FixedPositionInput(unittest.TestCase):
    def test_checked_source_and_post_ingestion_change(self):
        with TemporaryDirectory() as work:
            root = Path(work)
            event = {"obsid":27089, "trigger_key":"real/0", "instrument_pass":True,
                     "ra_deg":197., "dec_deg":-23., "event_mjd":59968., "position_error_arcsec":2.6,
                     "trial_duration_s":3600., "duration_s":2215., "diagnostics":{"files":{
                         "diagnostics.png":"032c0c6469ca3780f603158cb047b451988489df2d54fa1f097537c6ae16bc17"}}}
            source = {"sha256":field.PRODUCT, "filename":"test.fits.gz", "url":"https://example.org/test.fits.gz"}
            product = {"complete":True, "input":source, "search":{"triggers":[event]}}
            rel = "observations/27089.json.gz"
            ingest.save(root/rel, {"obsid":27089, "products":[product]}, zipped=True)
            members = root/ingest.HISTORY/"members/80.json.gz"
            ingest.save(members, [{"source":rel, "destination":rel, "sha256":ingest.sha(root/rel)}], zipped=True)
            ingest.save(ingest.receipt_file(root,"80"), {"group":"80", "artifact_id":field.ARTIFACT,
                        "zip_digest":field.ZIP_SHA, "code_commit":ingest.CODE,
                        "member_manifest_sha256":ingest.sha(members)})
            original, chosen, lead = field.select(root)
            self.assertEqual(original, source)
            self.assertEqual(chosen, event)
            self.assertEqual(lead["input_sha256"], field.PRODUCT)
            event["ra_deg"] += .001
            ingest.save(root/rel, {"obsid":27089, "products":[product]}, zipped=True)
            with self.assertRaisesRegex(ValueError, "Original observation"):
                field.select(root)
