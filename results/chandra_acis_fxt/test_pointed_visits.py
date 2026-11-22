"""Synthetic active-workflow fixtures for a pinned, multi-band discovery GTI."""
import csv
import gzip
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import ingest_primary as ingest
import pointed_bands as bands
import pointed_visits as pointed


class PointedBandTests(unittest.TestCase):
    def test_energy_boundaries_clipping_counts_and_exact_intervals(self):
        with TemporaryDirectory() as d:
            file=Path(d)/"rows.csv.gz"
            photons=[(-3600,"aperture",2000),(-180,"annulus",500),(0,"aperture",1999),
                     (300,"aperture",2000),(600,"aperture",7000),(899,"annulus",2000),(900,"aperture",500),
                     (1080,"annulus",7000),(4500,"aperture",4000)]
            with gzip.open(file,"wt",newline="") as f:
                w=csv.writer(f)
                w.writerow(["events_row_zero_based","region","phase","time_met_s","energy_eV","grade","in_peak"])
                for i,(t,region,e) in enumerate(photons):
                    phase="before" if t<0 else "trigger" if t<900 else "after"
                    w.writerow([i,region,phase,t,e,0,int(0<=t<60)])
            diag={"gti_met_s":[-3600.,4600.],"files":{"event_rows.csv.gz":ingest.sha(file)},
                  "sky_center_px":[1.,2.],"apertures_sky_px":[3.,6.,12.],"input_sha256":"a"*64,"event_table_rows":len(photons)}
            event={"diagnostics":diag,"start_met_s":0,"stop_met_s":900,"before_s":3600,"after_s":3700,
                   "live_fraction":.8,"annulus_area_ratio":12.,"before_counts":1,"on_counts":3,"after_counts":2,
                   "annulus_counts":1,"annulus_off_counts":2,"peak_start_met_s":0,"peak_window_s":60,"peak_counts":1,"peak_annulus_counts":0}
            result=bands.profile(file,event)
            rows={r["name"]:r for r in result["bins"]}
            self.assertEqual(len(rows),13)
            self.assertEqual(rows["earlier_baseline"]["live_s"],0)
            self.assertEqual(rows["earlier_baseline"]["aperture"]["gross_hard_fraction_cp95"],[None,None])
            self.assertEqual(rows["trial_first_third"]["aperture"]["soft"],1)
            self.assertEqual(rows["trial_middle_third"]["aperture"]["hard"],1)
            self.assertEqual(rows["trial_last_third"]["aperture"]["hard"],1)
            self.assertEqual(rows["after_0_180"]["aperture"]["soft"],1)
            self.assertEqual(rows["after_180_600"]["annulus"]["hard"],1)
            self.assertEqual(rows["later_baseline"]["aperture"]["hard"],1)
            tri=result["counts_by_phase"]["trigger"]
            self.assertAlmostEqual(tri["aperture"]["gross_hard_fraction"],2/3)
            self.assertEqual(tri["net_aperture_counts"]["hard"],2.-1/12)
            self.assertAlmostEqual(rows["trial_middle_third"]["aperture"]["gross_hard_fraction_cp95"][0],.025)
            self.assertEqual(rows["trial_middle_third"]["aperture"]["gross_hard_fraction_cp95"][1],1.)
            self.assertAlmostEqual(sum(r["live_s"] for r in rows.values()),8200*.8)
            with self.assertRaisesRegex(ValueError,"tally"):
                bands.profile(file,{**event,"before_counts":2})
            with self.assertRaisesRegex(ValueError,"SHA"):
                bands.profile(file,{**event,"diagnostics":{**diag,"files":{"event_rows.csv.gz":"0"*64}}})

    def test_nonzero_trigger_key_and_original_producer_identity(self):
        with TemporaryDirectory() as d:
            root=Path(d); oid=42; group=5; art=123; h="b"*64; png="c"*64
            rel=f"observations/{oid}.json.gz"
            events=[{"trigger_key":"real/0","obsid":oid,"instrument_pass":True,"diagnostics":{"files":{"diagnostics.png":"x"}}},
                    {"trigger_key":"real/1","obsid":oid,"instrument_pass":True,"diagnostics":{"files":{"diagnostics.png":png}}}]
            ingest.save(root/rel,{"obsid":oid,"products":[{"input":{"sha256":h},"complete":True,"search":{"triggers":events}}]},zipped=True)
            stamp={"source":rel,"destination":rel,"bytes":(root/rel).stat().st_size,"sha256":ingest.sha(root/rel)}
            members=root/ingest.HISTORY/"members/5.json.gz"
            ingest.save(members,[stamp],zipped=True)
            ingest.save(ingest.receipt_file(root,"5"),{"artifact_id":art,"code_commit":ingest.CODE,
                        "artifact_name":"chandra-evidence-5","plan_commit":ingest.PLAN,"member_manifest_sha256":ingest.sha(members)})
            _,selected,_,_=pointed.original(root,group,art,oid,h,png,"real/1")
            self.assertEqual(selected["trigger_key"],"real/1")
            with self.assertRaisesRegex(ValueError,"trigger key"):
                pointed.original(root,group,art,oid,h,png,"real/01")
            self.assertEqual([r[2] for r in pointed.COHORTS["high-latitude"]],[22929,23779])
            self.assertEqual([r[2] for r in pointed.COHORTS[pointed.SECOND]],[22558])
