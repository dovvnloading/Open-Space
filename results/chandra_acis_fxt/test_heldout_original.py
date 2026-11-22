"""The original held-out analysis cannot acquire later products or classifications."""
import gzip
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import campaign
import heldout_original as orig


class TestOriginalHeldOut(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root=Path(tmp.name)/"evidence"
        self.view=Path(tmp.name)/"view"
        (self.view/"observations").mkdir(parents=True)
        self.root.mkdir()
        self.name="observations/98765.json.gz"
        self.blob=gzip.compress(b'{"obsid":98765,"status":"error"}',mtime=0)
        self.member={"source":self.name,"destination":self.name,"bytes":len(self.blob),
                     "sha256":orig.digest(self.blob)}
        self.selection={"group":96,"obsid":98765,"exclusions":[]}

    def file(self,path,data):
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(data)
        return path

    def test_current_corrected_and_intermediate_backups_are_not_original(self):
        first=self.root/"history/ned-untyped-v1/10-1"
        second=self.root/"history/ned-untyped-v1/11-1"
        self.file(self.root/self.name,b"corrected")
        self.file(second/"before"/self.name,b"intermediate")
        expected=self.file(first/"before"/self.name,self.blob)
        plan=self.file(first/"plan.json",orig.encoded({"obsids":{"98765":{"primary_sha256":self.member["sha256"]}}}))
        self.file(first/"result.json",orig.encoded({"plan_sha256":orig.sha(plan),"changed_or_attempted":
                   {"98765":{"before_sha256":self.member["sha256"]}}}))
        rows=orig.mount_group(self.view,self.root,[self.selection],{self.name:self.member},{},[second,first])
        self.assertEqual((self.view/self.name).resolve(),expected)
        orig.coverage_check(rows,[{**self.selection,"ledger_sha256":self.member["sha256"],"disposition":"incomplete"}])
        with self.assertRaisesRegex(ValueError,"another version"):
            orig.coverage_check(rows,[{**self.selection,"ledger_sha256":"e"*64,"disposition":"evaluable"}])
        expected.unlink()
        with self.assertRaisesRegex(ValueError,"cannot be verified"):
            orig.locate(self.root,self.member,[first,second],{})

    def test_no_worker_member_omits_later_ledger_from_every_stream(self):
        self.file(self.root/self.name,self.blob)
        rows=orig.mount_group(self.view,self.root,[self.selection],{},{},[])
        self.assertFalse((self.view/self.name).exists())
        orig.coverage_check(rows,[{**self.selection,"ledger_sha256":None,"disposition":"unattempted"}])
        with self.assertRaisesRegex(ValueError,"later ledger"):
            orig.coverage_check(rows,[{**self.selection,"ledger_sha256":None,"disposition":"evaluable"}])

    def test_exact_interrupted_worker_wins_over_better_selected_checkpoint(self):
        logical=(orig.ingest.HISTORY/"interrupted_worker/96"/self.name).as_posix()
        self.file(self.root/self.name,b"preferred primary checkpoint")
        producer=self.file(self.root/logical,self.blob)
        got=orig.locate(self.root,{**self.member,"destination":logical},[],{})
        self.assertEqual(got,producer)

    def test_sealed_evaluation_is_read_only_and_detects_changed_output(self):
        dest=self.root/"held_out"
        dest.mkdir()
        projection={"source_implementation":{"vet.py":"original"},"selection_sha256":"a"*64,"ledgers":[]}
        def inputs(root,view,work):
            work.mkdir()
            return [],projection
        def evaluate(selection,view,work):
            for name in orig.OUTPUTS:
                self.file(work/name,b"fixture")
            self.file(work/"coverage.json",b"[]")
            self.file(work/"summary.json",orig.encoded({"primary_implementation":projection["source_implementation"],
                                                      "selection_sha256":projection["selection_sha256"]}))
        env={"GITHUB_RUN_ID":"1","GITHUB_RUN_ATTEMPT":"1","GITHUB_SHA":"f"*40}
        with (patch.object(orig,"inputs",side_effect=inputs),patch.object(orig.ingest,"git",return_value=b"f"*40),
              patch.dict(os.environ,env),patch.object(orig.heldout,"evaluate",side_effect=evaluate) as evaluator):
            orig.run(self.root,"evaluate")
            receipt=orig.sha(dest/orig.RECEIPT)
            orig.run(self.root,"verify")
            orig.run(self.root,"evaluate")
            self.assertEqual(evaluator.call_count,1)
            self.assertEqual(receipt,orig.sha(dest/orig.RECEIPT))
            (dest/"controls.json").write_text("altered")
            with self.assertRaisesRegex(ValueError,"result changed"):
                orig.run(self.root,"verify")


if __name__=="__main__":
    unittest.main()
