"""Synthetic transport-free comparisons with native visit reduction."""
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import campaign
import ingest_primary as ingest
import visits
import visit_workers as w


class VisitWorkerTests(unittest.TestCase):
    def test_native_product_and_full_multiplicity(self):
        leads=[{"obsid":100+i,"input_sha256":f"{i+1:064x}","trigger_key":"real/0","ra_deg":20.+i,
                "dec_deg":10.,"position_error_arcsec":2.,"event_mjd":59000.} for i in range(2)]
        def rec(i):
            return {"obsid":i,"filename":f"{i}.gz","bytes":i*1024,"sha256":f"{i:064x}",
                    "fits_sha256":f"{i+2:064x}","url":f"https://example.invalid/{i}"}
        matches=[(8,rec(8),leads),(9,rec(9),leads[1:])]
        env={"GITHUB_ACTIONS":"true","GITHUB_REF":"refs/heads/"+ingest.BRANCH,"GITHUB_EVENT_NAME":"push",
             "GITHUB_RUN_ID":"22","GITHUB_RUN_ATTEMPT":"1"}
        def download(manifest,directory):
            [r]=campaign.read(manifest)
            (directory/r["filename"]).write_bytes(b"synthetic")
        def measure(path,oid,rows):
            self.assertEqual(path.read_bytes(),b"synthetic")
            return {visits._lead_id(r):{"status":"measured_aperture","temporal_relation":"before",
                    "local_excess_p":.01,"conditional_aperture_rate_upper95_ct_s":.05} for r in rows}
        with tempfile.TemporaryDirectory() as d, patch.dict(os.environ,env):
            root=Path(d)
            out,other=root/"independent_visits",root/"reference"
            ingest.save(root/"leads.json",leads)
            ingest.save(root/"inventory/inventory.json",[{"obsid":8},{"obsid":9}])
            with patch.object(visits,"_receipt_plan",return_value=([8,9],matches,[])), patch.object(
                    visits,"replay",download),patch.object(visits,"_measure",measure):
                plan=w.plan(root,2)
                visits.run_visits(root,leads,other)
                self.assertEqual(sorted(i for g in plan["groups"] for i in g),[0,1])
                for j in reversed(plan["jobs"]):
                    w.perform(root,out,out/"exec",j,leads,[])
                    self.assertEqual((out/j["record"]).read_bytes(),(other/j["record"]).read_bytes())
                with patch.object(visits,"replay",side_effect=AssertionError("duplicate fetch")):
                    for j in plan["jobs"]:
                        w.perform(root,out,out/"exec",j,leads,[])
                    s=visits.run_visits(root,leads,out)
                self.assertEqual((out/"independent_visits.json").read_bytes(),(other/"independent_visits.json").read_bytes())
                self.assertEqual([r["per_lead_bonferroni_p"] for r in s["visits"][visits._lead_id(leads[1])]],[.02,.02])
                j=plan["jobs"][0]
                p=out/j["record"]
                corrupt=campaign.read(p)
                corrupt["lead_ids"]=[]
                ingest.save(p,corrupt)
                before=ingest.sha(p)
                self.assertEqual(w.state(out,j,leads),"invalid_checkpoint")
                w.perform(root,out,out/"exec",j,leads,[])
                self.assertEqual(w.state(out,j,leads),"complete")
                self.assertEqual(ingest.sha(out/"exec/previous"/j["record"]/(before+".json")),before)

    def test_zero_and_disjoint_assignment(self):
        self.assertEqual(w.partition([],4),[])
        jobs=[{"receipt":{"bytes":n},"selected":[1]*i} for n,i in [(300000000,1),(16,3),(120000000,4),(90,1)]]
        groups=w.partition(jobs,2)
        self.assertEqual(groups,w.partition(jobs,2))
        self.assertCountEqual([k for g in groups for k in g],range(len(jobs)))
