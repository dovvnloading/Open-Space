"""Synthetic stage and adjudication regressions for the final workflow."""
from copy import deepcopy
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import final_assembly as manual
import finalize


class FinalCheck(unittest.TestCase):
    def test_legacy_reviews_must_still_describe_the_same_measurement(self):
        k=(11,"a"*64,"real/2")
        row={"decision":"unresolved","reason":"Sky source with undetermined origin.",
             "evidence":["fixture.md"],"obsid":k[0],"input_sha256":k[1],"trigger_key":k[2]}
        legacy=[("PRIMARY_ODD_DECISIONS.json",row)]
        current="b"*64
        self.assertIsNone(manual.choices(k,current,[],True)[0])
        accepted,issue=manual.choices(k,current,legacy,True)
        self.assertIsNone(issue)
        self.assertEqual(accepted["decision"],"unresolved")
        self.assertEqual(manual.choices(k,current,legacy,False)[1],"reviewed_event_or_catalog_changed")
        old_candidate={**row,"decision":"candidate"}
        self.assertEqual(manual.choices(k,current,[("old",old_candidate)],True)[1],"candidate_needs_final_measurement_review")
        override={**row,"review_sha256":current,"decision":"rejected"}
        updated,issue=manual.choices(k,current,legacy+[("new",override)],False)
        self.assertIsNone(issue)
        self.assertEqual(updated["decision"],"rejected")
        self.assertEqual(manual.choices(k,current,[("a",override),("b",{**override,"decision":"candidate"})],True)[1],
                         "conflicting_current_manual_decisions")

    def test_fingerprint_carries_request_bytes_but_ignores_local_cache_location(self):
        row={"instrument_pass":True,"catalog_cache":"old","vet":{"status":"ambiguous","queries":{"ned":{
            "ok":True,"sha256":"x","url":"https://example.invalid","attempts":[]}}}}
        other=deepcopy(row)
        other["catalog_cache"]="temporary-retry"
        self.assertEqual(manual.event_state(row),manual.event_state(other))
        other["vet"]["queries"]["ned"]["sha256"]="y"
        self.assertNotEqual(manual.event_state(row),manual.event_state(other))

    def test_stages_reject_changed_bytes_foreign_observations_and_shared_symlinks(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ,{"GITHUB_RUN_ID":"1","GITHUB_RUN_ATTEMPT":"1"}):
            stage=Path(tmp)
            f=stage/"observations/11.json.gz"
            f.parent.mkdir()
            f.write_bytes(b"fixture")
            plan={"groups":[[11]]}
            finalize.worker_ownership(finalize.tree(stage),0,plan,"1","1")
            with self.assertRaisesRegex(ValueError,"Foreign"):
                finalize.worker_ownership({"observations/12.json.gz":None},0,plan,"1","1")
            with self.assertRaisesRegex(ValueError,"Unexpected"):
                finalize.worker_ownership({"catalogs/followup-2-1/11/file.bin":None},0,plan,"1","1")
            original=finalize.stamped(stage,"workers",{"group":0})
            self.assertEqual(original,finalize.verified(stage,"workers"))
            f.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError,"sealed bytes"):
                finalize.verified(stage,"workers")
            (stage/"alias").symlink_to(f)
            with self.assertRaisesRegex(ValueError,"symlink"):
                finalize.tree(stage)
