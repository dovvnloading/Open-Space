"""Synthetic artifact ingress and resumption; no GitHub or science data are used."""
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from zipfile import ZipFile

import ingest_primary as recover


def raw(data):
    return json.dumps(data,sort_keys=True).encode()


def compressed(data):
    return gzip.compress(raw(data),mtime=0)


def observation(oid):
    h=hashlib.sha256(str(oid).encode()).hexdigest()
    rec={"filename":f"acisf{oid:05d}N001_evt2.fits.gz","obsid":oid,"sha256":h,"bytes":16}
    support={k:[] for k in ("60","180","600","1800","3600")}
    stream={"obsid":oid,"seed":int(hashlib.sha256(("chandra-acis-fxt-v1:"+h).encode()).hexdigest()[:16],16)>>1,
            "triggers":[],"controls":[{"replica":i,"triggers":[]} for i in range(24)],
            "injections":[{"trial":i,"attempted":False,"recovered":None} for i in range(12)],
            "start_support_by_duration":support,"exposure_by_duration_s":{k:0 for k in support},"gti_livetime_s":0}
    return {"obsid":oid,"status":"complete","errors":[],"acquired":[rec],
            "products":[{"input":rec,"complete":True,"search":stream}]}


class FrozenFixture:
    def __init__(self):
        self.groups=[[11,12]]+[[] for _ in range(191)]
        self.population={11,12}
        self.assigned={11,12}
        self.sources={}
        self.supplement=set()
        self.blobs={}
        self.prior={}

    def baseline(self,name):
        return self.prior[name]


class IngressTests(unittest.TestCase):
    def setUp(self):
        area=tempfile.TemporaryDirectory()
        self.addCleanup(area.cleanup)
        self.area=Path(area.name)
        self.root=self.area/"root"
        self.root.mkdir()
        self.frozen=FrozenFixture()
        self.name="chandra-evidence-0"
        self.ns=f"{recover.RUN}-1-0"
        self.z=self.area/(self.name+".zip")
        self.meta={"run":{"status":"completed"},"jobs":{"0":{"conclusion":"success"}},"artifacts":{}}
        self.files={"environment-0.json":raw({"run_id":str(recover.RUN),"git_commit":recover.CODE}),
                    "shard-0.json":raw({"group":0,"assigned":[11,12],"commit":recover.CODE,"completed_worker":True}),
                    "observations/11.json.gz":compressed(observation(11)),
                    f"catalogs/{self.ns}/fixture.bin":b"original catalog bytes"}

    def archive(self):
        with ZipFile(self.z,"w") as z:
            for name,data in self.files.items():
                z.writestr(name,data)
        self.meta["artifacts"][self.name]={"id":10,"name":self.name,"digest":"sha256:"+recover.sha(self.z),
                                           "size_in_bytes":self.z.stat().st_size,"expired":False}

    def ingest(self,packed=False):
        with patch.object(recover.catalog_pack,"protected",return_value=set()), patch.object(
                recover,"api",side_effect=AssertionError("unexpected network access")):
            return recover.import_one(self.root,self.frozen,self.meta,"0",packed,self.area)

    def prior(self,row):
        name="observations/11.json.gz"
        self.frozen.prior[name]=compressed(row)
        self.frozen.blobs[name]="test-original-git-blob"
        p=self.root/name
        p.parent.mkdir(exist_ok=True)
        p.write_bytes(self.frozen.prior[name])

    def test_packed_import_resumes_without_source_zip_and_reports_unattempted(self):
        previous={"obsid":11,"status":"error","errors":["archive service"],"products":[]}
        self.prior(previous)
        self.archive()
        rec=self.ingest(True)
        self.assertEqual("complete",rec["worker_ledger_status"]["11"]["selected_status"])
        self.assertEqual("unattempted",rec["worker_ledger_status"]["12"]["worker_status"])
        retained=self.root/recover.HISTORY/"checkpoint_before/0/observations/11.json.gz"
        self.assertEqual(self.frozen.prior["observations/11.json.gz"],retained.read_bytes())
        self.assertFalse((self.root/"catalogs"/self.ns/"fixture.bin").exists())
        self.z.unlink()
        self.assertEqual(rec,self.ingest(True))
        result=recover.audit(self.root,self.frozen,self.meta)
        self.assertFalse(result["complete_artifact_transport"])
        self.assertEqual([12],result["groups"]["0"]["assigned_without_verified_worker_ledger"])
        self.assertEqual([12],result["dated_residual_obsids"])

    def test_regressive_retry_is_preserved_separately_from_preexisting_scan(self):
        old=observation(11)
        old["status"]="followup_incomplete"
        self.prior(old)
        latest={"obsid":11,"status":"error","errors":["worker deadline"],"products":[]}
        self.files["observations/11.json.gz"]=compressed(latest)
        self.archive()
        rec=self.ingest()
        self.assertEqual("checkpoint_after_interrupted_retry",rec["worker_ledger_status"]["11"]["selected_source"])
        self.assertEqual(compressed(old),(self.root/"observations/11.json.gz").read_bytes())
        self.assertEqual(compressed(latest),(self.root/recover.HISTORY/"interrupted_worker/0/observations/11.json.gz").read_bytes())

    def test_conflicting_scans_and_foreign_paths_are_rejected(self):
        old=observation(11)
        new=deepcopy(old)
        new["products"][0]["search"]["controls"][2]["triggers"].append({"new_summary":True})
        with self.assertRaisesRegex(ValueError,"frozen real scan"):
            recover.revisit(old,new)
        self.files["observations/875.json.gz"]=compressed(observation(875))
        self.archive()
        with self.assertRaisesRegex(ValueError,"Foreign artifact path"):
            self.ingest()
        self.assertFalse(recover.receipt_file(self.root,"0").exists())

    def test_changed_zip_and_missing_control_stream_are_not_certified(self):
        self.archive()
        with self.z.open("ab") as f:
            f.write(b"extra")
        with self.assertRaisesRegex(ValueError,"SHA-256 mismatch"):
            self.ingest()
        row=observation(11)
        row["products"][0]["search"]["controls"].pop()
        self.files["observations/11.json.gz"]=compressed(row)
        self.archive()
        with self.assertRaisesRegex(ValueError,"controls or injections"):
            self.ingest()
        self.assertFalse(recover.receipt_file(self.root,"0").exists())
