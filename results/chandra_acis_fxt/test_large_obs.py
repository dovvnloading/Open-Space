"""Offline illustrations of group-1 recovery, run in the Chandra Actions workflow."""
from copy import deepcopy
import os
from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import campaign
import diagnostics
import large_obs as big
import vet


class TestLargeProduct(unittest.TestCase):
    def setUp(self):
        temp = TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.temp = Path(temp.name)
        self.primary = self.temp/"primary"
        self.staged = self.temp/"stage"
        self.primary.mkdir()
        self.staged.mkdir()

    def test_real_control_assignment_is_exhaustive_balanced_and_keeps_classified_rows(self):
        raw = {"ra_deg":1.0, "dec_deg":1.0, "instrument_pass":False}
        prior = {**raw,"vet":{"queries":{"ned":{"ok":True}}}}
        toy = {"search":{"triggers":[raw,raw,prior],
                         "controls":[{"replica":i,"triggers":[dict(raw)]} for i in range(24)]}}
        rows = big.assignments(toy,4)
        self.assertEqual(len(rows),27)
        self.assertIsNone(rows["real/2"]["group"])
        counts = [sum(r["group"] == g for r in rows.values()) for g in range(4)]
        self.assertEqual(sum(counts),26)
        self.assertLessEqual(max(counts),7)
        self.assertEqual(rows["control-23/0"]["raw_sha256"],big.fix.value_sha(raw))

    def test_new_row_recovers_exact_original_http_reply_without_network(self):
        rel = Path(big.SOURCE_CACHE)
        source = self.primary/rel
        source.mkdir(parents=True)
        cache = self.staged/"catalogs"/"part-0"
        cache.mkdir(parents=True)
        url = "https://fixture.invalid/ok"
        stem = "ned-" + vet.digest(url.encode())[:24]
        payload = b"primary source reply bytes"
        body = stem+".bin"
        (source/body).write_bytes(payload)
        memo = {"url":url,"ok":True,"attempts":[],"sha256":vet.digest(payload),
                "size_bytes":len(payload),"body_file":body}
        big.save(source/(stem+".json"),memo)
        name = str(rel/(stem+".json"))
        members = {name:{"destination":name,"sha256":big.sha(self.primary/name)}}
        issues = []
        read = big.old_transport(self.primary,members,vet.read_input,issues)
        records = {}
        self.assertEqual(read("ned",url,cache,True,records),payload)
        self.assertEqual(records["ned"],memo)
        self.assertEqual(issues,[])
        (cache/(stem+".json")).unlink()
        (source/body).write_bytes(b"altered")
        with self.assertRaisesRegex(ValueError,"missing or altered"):
            read("ned",url,cache,True,records)
        self.assertTrue(issues)

    def test_group_receipts_require_all_assigned_rows_and_raw_source_bytes(self):
        ph="p"*64
        key="control-0/0"
        event={"ra_deg":1.5,"dec_deg":1.0,"instrument_pass":False}
        h=big.fix.value_sha(event)
        plan={"groups":1,"rows":{key:{"group":0,"raw_sha256":h}}}
        cache=self.staged/"catalogs"/"part-0"
        cache.mkdir(parents=True)
        (cache/"example.bin").write_bytes(b"preserved unsuccessful response")
        producer={"run":"1","attempt":"1","commit":"0"*40,"job":"catalogs"}
        packet={"producer":producer,"plan_sha256":ph,"group":0,"row":key,"raw_sha256":h,"trigger":
                {**event,"obsid":big.OID,"trigger_key":key},"verified":False}
        path=big.packet_path(self.staged,key)
        big.save(path,packet,True)
        seal={"producer":producer,"plan_sha256":ph,"group":0,"rows":{key:big.sha(path)},"source_files":big.cache_files(cache)}
        big.save(self.staged/"group.json",seal)
        records,groups=big.load_groups(plan,ph,self.staged)
        self.assertEqual(campaign.read_obs(records[key][1]),packet)
        self.assertEqual(groups,{0:self.staged})
        (cache/"example.bin").write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError,"Sealed response"):
            big.load_groups(plan,ph,self.staged)
        (self.staged/"group.json").unlink()
        with self.assertRaisesRegex(ValueError,"Missing catalog work group"):
            big.load_groups(plan,ph,self.staged)

    def test_raw_unit_result_and_distributed_catalogs_use_stock_end_to_end_processor(self):
        # A preassembled raw search and real/control catalog responses feed the
        # ordinary process. Only the synthetic acquisition and CIAO calls are mocked.
        repo,bundle,prep,out=self.primary,self.temp/"bundle",self.temp/"prep",self.temp/"worker"
        producer={"run":"1234","attempt":"1","commit":"0"*40,"job":"unit-test"}
        event=self.temp/big.EVENT[0]
        event.write_bytes(b"fixture event archive, recovery integration")
        h=big.sha(event)
        receipt={"filename":event.name,"bytes":event.stat().st_size,"sha256":h,"fits_sha256":"2"*64,
                 "obsid":big.OID,"url":"https://fixture.invalid/input"}
        raw={"obsid":big.OID,"seed":int(vet.digest(("chandra-acis-fxt-v1:"+h).encode())[:16],16)>>1,
             "triggers":[{"ra_deg":1.0,"dec_deg":2.0,"instrument_pass":True}],
             "controls":[{"replica":i,"triggers":[]} for i in range(24)],
             "injections":[{"trial":i,"attempted":True,"recovered":i==0} for i in range(12)]}
        raw["controls"][0]["triggers"]=deepcopy(raw["triggers"])
        prior={"obsid":big.OID,"status":"error","acquired":[receipt],"products":[{"input":receipt}],
               "catalog_cache":"catalogs/primary","errors":["time limit"]}
        old=big.fix.old_path(repo,big.OID)
        big.save(old,prior,True)
        p={"input":receipt,"search":raw}
        plan={"source":{"ledger_sha256":big.sha(old)},"raw_science_sha256":big.fix.raw_scan(p),
              "primary_scan":"absent","assembly_sha256":None,"groups":1,"rows":big.assignments(p,1)}
        big.save(prep/"plan.json",plan)
        ph=big.sha(prep/"plan.json")
        big.save(bundle/"plan.json",{"obsids":{str(big.OID):{"primary_sha256":big.sha(old),"reasons":["error"]}},
                 "repair_sha256":big.sha(Path(big.fix.__file__)),"inventory_sha256":big.fix.INVENTORY,
                 "catalog_reader_sha256":big.sha(Path(big.catalog_pack.__file__))})
        cache=self.staged/"catalogs"/"part-0"
        cache.mkdir(parents=True)
        queries={}
        for label in ("gaia","simbad","ned","quirola2023","dillmann2025","cudidi2026"):
            url="https://fixture.invalid/"+label
            stem=label+"-"+vet.digest(url.encode())[:24]
            body=stem+".bin"
            content=("synthetic recorded table "+label).encode()
            (cache/body).write_bytes(content)
            q={"url":url,"ok":True,"attempts":[],"body_file":body,"size_bytes":len(content),"sha256":vet.digest(content)}
            big.save(cache/(stem+".json"),q)
            queries[label]={**q,"row_count":0}
        v={"status":"no_counterpart","errors":{},"queries":queries}
        seals={}
        for stream,i,subject in big.fix.triggers(p):
            key=f"{stream}/{i}"
            file=big.packet_path(self.staged,key)
            big.save(file,{"plan_sha256":ph,"producer":producer,"group":0,"raw_sha256":plan["rows"][key]["raw_sha256"],
                          "trigger":{**subject,"obsid":big.OID,"trigger_key":key},"row":key,"verified":True,"vet":v},True)
            seals[key]=big.sha(file)
        big.save(self.staged/"group.json",{"group":0,"producer":producer,"plan_sha256":ph,"rows":seals,
                                          "source_files":big.cache_files(cache)})
        def exact_input(pin,oid,directory,retained,records):
            self.assertEqual(pin,receipt)
            shutil.copyfile(event,directory/event.name)
            records.append({"filename":event.name,"fits_sha256":pin["fits_sha256"]})
        with (patch.object(big,"open_plan",return_value=(plan,p,{})),patch.object(big.fix,"checked"),
              patch.object(big.fix,"event_bytes",side_effect=exact_input),patch.object(big.fix,"load_vet",return_value=vet),
              patch.dict(sys.modules,{"vet":vet}),patch.object(vet,"vet_trigger",side_effect=AssertionError("external request")),
              patch.object(campaign,"response",return_value={"peak_flux_erg_cm2_s":1.0}) as flux,
              patch.object(diagnostics,"save_diagnostics",return_value={"files":{"lightcurve":"fixture"}}) as diagnostic,
              patch.dict(os.environ,{"GITHUB_RUN_ID":"1234","GITHUB_RUN_ATTEMPT":"1","RETAIN_EVENTS":str(self.temp/"kept")})):
            result=big.finish(repo,prep,self.staged,bundle,out,[])
        self.assertEqual(result["outcome"],"complete",result)
        self.assertEqual(campaign.read_obs(old),prior)
        revised=campaign.read_obs(big.fix.old_path(out,big.OID))["products"][0]
        self.assertEqual(big.fix.raw_scan(revised),big.fix.raw_scan(p))
        self.assertEqual(revised["search"]["injections"],raw["injections"])
        self.assertEqual(revised["search"]["controls"][0]["triggers"][0]["vet"],v)
        self.assertNotIn("flux",revised["search"]["controls"][0]["triggers"][0])
        self.assertEqual(flux.call_count,1)
        self.assertEqual(diagnostic.call_count,1)


if __name__ == "__main__":
    unittest.main()
