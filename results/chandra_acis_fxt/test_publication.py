"""Publish a complete mini-survey using synthetic archived measurements and written decisions."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch

import publication as pub


def digest(data):
    return hashlib.sha256(data).hexdigest()


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    data=(json.dumps(value,sort_keys=True)+"\n").encode()
    path.write_bytes(gzip.compress(data,mtime=0) if path.suffix==".gz" else data)
    return path


def table(path,fields,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def source(oid):
    raw=f"Synthetic EVT2 payload for fixture observation {oid}".encode()
    wire=gzip.compress(raw,mtime=0)
    name=f"acisf{oid}N001_evt2.fits.gz"
    return dict(obsid=oid,filename=name,sha256=digest(wire),fits_sha256=digest(raw),bytes=len(wire),
                url=f"https://cxc.cfa.harvard.edu/cdaftp/byobsid/{oid%10}/{oid}/primary/{name}")


def fixture_lead(evidence,receipt,disposition,index=0):
    oid,h=receipt["obsid"],receipt["sha256"]
    candidate=disposition=="candidate"
    tag=f"real/{index}"
    folder=evidence/"photometry"/str(oid)/receipt["filename"]/tag
    weights=save(folder/"peak_weights.fits.gz",{"fixture":"mock calibrated photon weights"})
    diag=save(folder/"diagnostics.json",{"fixture":"simulated aperture diagnostics"})
    host={"name":"Fixture Galaxy","redshift":.12,"separation_arcsec":.45} if candidate else None
    event=dict(ra_deg=120.25,dec_deg=30.5,position_error_arcsec=.9,event_mjd=59500.5,
               duration_s=46.5,trial_duration_s=60,trigger_key=tag,instrument_pass=True,
               on_counts=21,annulus_counts=4,annulus_area_ratio=12.,before_counts=1,before_s=1500.,
               after_counts=0,after_s=2400.,seed_score=9. if candidate else 8.,seed_local_p=1e-9 if candidate else 1e-8,
               flux={"peak_flux_erg_cm2_s":2e-12,"flux_counting_sigma_erg_cm2_s":3e-13,
                     "gross_flux_erg_cm2_s":2.3e-12,"scaled_background_flux_erg_cm2_s":.3e-12,
                     "event_sha256":h,"calibrated_fits_sha256":digest(weights.read_bytes())},
               diagnostics={"input_sha256":h,"files":{"diagnostics.json":digest(diag.read_bytes())}},
               vet={"status":"host_candidate" if candidate else "no_counterpart","errors":{},"host":host,
                    "matches":[{"catalog":"fixture catalog","name":host["name"],"separation_arcsec":.45}] if host else []})
    if disposition=="unresolved":
        event["flux"]={"error":"Missing fixture calibration"}
    if disposition=="rejected":
        event["flux"].update(peak_flux_erg_cm2_s=-2e-12,gross_flux_erg_cm2_s=1e-12,scaled_background_flux_erg_cm2_s=3e-12)
    p=evidence/"observations"/f"{oid}.json.gz"
    if index:
        ledger=json.loads(gzip.decompress(p.read_bytes()))
        ledger["products"][0]["search"]["triggers"].append(event)
    else:
        ledger={"obsid":oid,"status":"followup_incomplete" if disposition=="unresolved" else "complete","products":[
                {"input":receipt,"complete":True,"search":{"seed":int(h[:12],16),"triggers":[event]}}]}
    save(p,ledger)
    lo,hi=pub.report.bounds(1,24)
    return {**{k:event[k] for k in ("ra_deg","dec_deg","position_error_arcsec","event_mjd","duration_s","trial_duration_s")},
            "obsid":oid,"trigger_key":tag,"input_sha256":h,"input_name":receipt["filename"],"input_url":receipt["url"],
            "survey_fap":2/25,"control_exceedances":1,"null_reference":"instrument_maximum","status":event["vet"]["status"],
            "peak_flux_erg_cm2_s":event["flux"].get("peak_flux_erg_cm2_s"),"fap_lo95":lo,"fap_hi95":hi}


def build_supplement(supplement):
    meta={}
    for name in pub.GROUPS:
        measured=name=="unresolved_date_sky"
        directory=supplement/name
        oid=92001
        cov=[{"obsid":oid,"status":"complete","scanned_all_products":True}] if measured else []
        table(directory/"coverage.csv",["obsid","status","scanned_all_products"],cov)
        src=[{**source(oid),"included_in_inference":True}] if measured else []
        table(directory/"event_inputs.csv",["obsid","filename","sha256","fits_sha256","bytes","url","included_in_inference"],src)
        save(directory/"control_counts.json",[{"replica":i,"retained_or_unresolved":0} for i in range(24)])
        save(directory/"injections.json",[{"duration_s":60,"planned":2,"attempted":2,"recovered":1}] if measured else [])
        save(directory/"leads.json",[])
        meta[name]=dict(eligible_start_support_s={"60":500000 if measured else 0},obsids=[oid] if measured else [],
                        fully_scanned=int(measured),followup_complete=int(measured),mean_control_retained_or_unresolved=0 if measured else None,
                        attempted_injections=2 if measured else 0,recovered_injections=int(measured),
                        surviving_or_unresolved=0,real_raw_triggers=0)
    save(supplement/"summary.json",{"groups":meta})
    save(supplement/"release_bounds.json",{"obsids":[{"obsid":92001,"sky_status":"nominal_imaging","proven_in_window":False}]})


class TestPublication(unittest.TestCase):
    def test_negative_signed_flux_and_component_disagreement(self):
        f = dict(peak_flux_erg_cm2_s=-4.22e-11, flux_counting_sigma_erg_cm2_s=4.41e-11,
                 gross_flux_erg_cm2_s=1.9e-12, scaled_background_flux_erg_cm2_s=4.41e-11)
        self.assertEqual(pub.flux_components(f), (1.9e-12,4.41e-11))
        f["peak_flux_erg_cm2_s"] *= -1
        with self.assertRaisesRegex(ValueError, "inconsistent aperture"):
            pub.flux_components(f)

    def test_candidate_unresolved_separate_date_bounds_and_catalog_counts(self):
        with tempfile.TemporaryDirectory() as d:
            repo=Path(d)
            here=repo/"results/chandra_acis_fxt"
            evidence=here/"evidence"
            supplement=evidence/"supplemental_search/synthetic"
            corrected={"vet.py":"9f3615c635e42d60425be4738dd3ffd62dfe8aee74bc717736df951af2a996e8"}
            receipt=[source(91001),source(91002)]
            leads=[fixture_lead(evidence,receipt[0],"candidate"),fixture_lead(evidence,receipt[1],"unresolved"),
                   fixture_lead(evidence,receipt[0],"rejected",1)]
            save(evidence/"leads.json",leads)
            save(evidence/"inventory/campaign.json",{"implementation":corrected,"followup_revision":"history/fixture"})
            save(evidence/"history/fixture/result.json",{"changed_or_attempted":{}})
            table(evidence/"obsid_coverage.csv",["obsid","status","scanned_all_products"],
                  [{"obsid":oid,"status":status,"scanned_all_products":i<2}
                   for i,(oid,status) in enumerate(((91001,"complete"),(91002,"followup_incomplete"),(91003,"unattempted")))])
            table(evidence/"event_receipts.csv",[*receipt[0],"included_in_inference"],
                  [{**r,"included_in_inference":True} for r in receipt])
            table(evidence/"control_trials.csv",["replica","surviving_or_unresolved","scanned_obsids","maximum_score","instrument_maximum"],
                  [{"replica":i,"surviving_or_unresolved":int(i==0),"scanned_obsids":2,
                    "maximum_score":10 if i==0 else "","instrument_maximum":10 if i==0 else ""} for i in range(24)])
            grid=[{"duration_s":length,"input_counts":counts,"psf_sigma_sky_px":1.,"scheduled":planned,
                   "attempted":attempted,"recovered":recovered} for length,counts,planned,attempted,recovered in
                  [(60,8,3,3,2),(180,16,1,1,1),(600,8,1,0,0)]]
            table(evidence/"completeness.csv",list(grid[0]),grid)
            save(evidence/"summary.json",dict(exposure_by_duration_s={"60":2000000,"180":800000},scanned_obsids=2,
                 completed_obsids=1,inventory_obsids=3,scheduled_injections=5,attempted_injections=4,recovered_injections=3,
                 control_replicas=24,mean_control_survivors=1/24,scientific_leads=3,input_products=2,real_raw_triggers=9,
                 real_instrument_pass=3))
            build_supplement(supplement)
            note=here/"REVIEW.md"
            note.write_text("# Fixture review\n"+"".join(
                f"\n## ObsID {r['obsid']}\n\nObsID {r['obsid']} real/0 and real/1 EVT2 {r['sha256']}\n" for r in receipt))
            field=evidence/"field_followup/27089/45008126305-1"
            save(field/"history/history.json",{})
            frame=evidence/"frame_forensics/24912"/("0"*64)/"real/0/fixture"
            save(frame/"result.json",{"status":"complete","measurement":{"instrument_ccd":8},
                                      "target":{"obsid":24912,"trigger":"real/0"}})
            (frame/"README.md").write_text("Synthetic S4 result\n")
            (here/"METHODS.md").write_text("# Synthetic fixture\n")
            old_entries=["- Another measurement. *Measured*","- Another correction. *Correction*"]
            (repo/"README.md").write_text("\n"*8+"3 entries, 2 new measurements, 0 candidate discoveries, 0 forecasts, "
                 "1 corrections and 0 carefully recorded null results.\n"+"\n".join(old_entries)+
                 f"\n- Pending survey. *Measured* · [Details]({pub.report.TAG})\n")
            self.addCleanup(pub.file_sha.cache_clear)
            self.addCleanup(pub.obs.cache_clear)
            with (patch.object(pub,"HERE",here),patch.object(pub.report,"ROOT",repo),patch.object(pub,"identity",return_value=corrected)):
                groups,measured,signature,sources=pub.load(supplement)
                written=[{"obsid":r["obsid"],"input_sha256":r["input_sha256"],"trigger_key":r["trigger_key"],
                          "decision":decision,"reason":reason,"evidence":[f"REVIEW.md#obsid-{r['obsid']}"]}
                         for r,decision,reason in zip(leads,["candidate","unresolved","rejected"],
                             ["A plausible galaxy association survived review.","The source distance remains unknown.",
                              "The CCD rows share a detector column."],strict=True)]
                manifest=save(evidence/"decisions.json",{"reviewed_inputs_sha256":signature,"decisions":written})
                decisions,cited=pub.adjudicate(campaign_read(manifest),measured,signature)
                result=pub.publish(groups,measured,decisions,sources,cited,signature,manifest,supplement)
            self.assertEqual(campaign_read(evidence/"publication.json"),result)
            self.assertEqual(result["dispositions"]["dated"],{"candidate":1,"unresolved":1,"rejected":1})
            self.assertEqual(result["label"],"Candidate")
            self.assertEqual(result["reviewed_inputs_sha256"],signature)
            self.assertEqual(result["decisions_sha256"],digest(manifest.read_bytes()))
            with (evidence/"fxt_candidates.csv").open(newline="") as f:
                rows=list(csv.DictReader(f))
            self.assertEqual(len(rows),1)
            row=rows[0]
            self.assertEqual((row["obsid"],row["trigger_key"],row["host"],row["nearest_name"]),
                             ("91001","real/0","Fixture Galaxy","Fixture Galaxy"))
            self.assertEqual((float(row["ra_deg"]),float(row["dec_deg"]),float(row["position_error_arcsec"]),
                              float(row["peak_flux_erg_cm2_s"]),float(row["duration_s"])),(120.25,30.5,.9,2e-12,46.5))
            self.assertEqual((row["event_sha256"],row["decompressed_fits_sha256"],row["event_url"],int(row["event_bytes"])),
                             (receipt[0]["sha256"],receipt[0]["fits_sha256"],receipt[0]["url"],receipt[0]["bytes"]))
            self.assertEqual((int(row["control_exceedances_of_24"]),float(row["survey_fap"])),(1,2/25))
            self.assertLess(float(row["exceedance_lo95"]),1/24)
            self.assertGreater(float(row["exceedance_hi95"]),1/24)
            self.assertEqual(result["candidate_table_sha256"],digest((evidence/"fxt_candidates.csv").read_bytes()))
            catalog=evidence/"reviewed_fxt_triggers.csv"
            with catalog.open(newline="") as f:
                every=list(csv.DictReader(f))
            self.assertEqual(result["reviewed_rows"],3)
            self.assertEqual(result["candidate_rows"],1)
            self.assertEqual(result["reviewed_table_sha256"],digest(catalog.read_bytes()))
            self.assertEqual([r for r in every if r["disposition"]=="candidate"],rows)
            by_state={r["disposition"]:r for r in every}
            self.assertEqual(set(by_state),{"candidate","unresolved","rejected"})
            unavailable=by_state["unresolved"]
            self.assertEqual((unavailable["photometry_status"],unavailable["photometry_error"]),
                             ("unavailable","Missing fixture calibration"))
            self.assertTrue(all(unavailable[k]=="" for k in ("peak_flux_erg_cm2_s","counting_sigma_erg_cm2_s",
                          "gross_peak_flux_erg_cm2_s","scaled_peak_background_erg_cm2_s","calibrated_weights_sha256")))
            rejected=by_state["rejected"]
            self.assertEqual((rejected["trigger_key"],rejected["event_sha256"]),("real/1",receipt[0]["sha256"]))
            self.assertEqual(tuple(float(rejected[k]) for k in ("peak_flux_erg_cm2_s","gross_peak_flux_erg_cm2_s",
                             "scaled_peak_background_erg_cm2_s","counting_sigma_erg_cm2_s")),(-2e-12,1e-12,3e-12,3e-13))
            self.assertEqual(float(rejected["seed_score"]),8.)
            self.assertEqual(float(rejected["seed_local_p"]),1e-8)
            self.assertEqual(int(rejected["search_random_seed"]),int(receipt[0]["sha256"][:12],16))
            self.assertEqual(rejected["control_reference"],"instrument_maximum")
            self.assertEqual(rejected["source_ledger_sha256"],digest((here/rejected["source_ledger"]).read_bytes()))
            self.assertEqual(json.loads(unavailable["human_evidence"]),["REVIEW.md#obsid-91002"])
            plot=(evidence/"injection_recovery.png").read_bytes()
            self.assertEqual(plot[:8],b"\x89PNG\r\n\x1a\n")
            width,height=struct.unpack(">II",plot[16:24])
            self.assertGreater(width,height)
            self.assertGreater(height,400)
            self.assertEqual(result["recovery_figure_sha256"],digest(plot))
            bounds=result["conditional_upper95_detected_events_per_Ms"]
            self.assertTrue(all(v is None for v in bounds["dated"].values()))
            self.assertAlmostEqual(bounds["unresolved_date_sky"]["60"],-2*math.log(.05))
            self.assertIsNone(bounds["unresolved_date_sky"]["180"])
            conservative=result["instrument_pass_conservative_upper95_detected_events_per_Ms"]
            review_rates=result["review_selection_upper95_events_per_Ms"]
            self.assertEqual(review_rates,{})
            self.assertIn("Incomplete counterpart or photometric follow-up prevents an upper bound",
                          (here/"README.md").read_text())
            fixture_exposure={str(d):2_000_000 if d==60 else 0 for d in pub.DURATIONS}
            accepted_mu=2*pub.instrument_rate_bounds(1,fixture_exposure)["60"]
            possible_mu=2*pub.instrument_rate_bounds(2,fixture_exposure)["60"]
            self.assertAlmostEqual(math.exp(-accepted_mu)*(1+accepted_mu),.05)
            self.assertAlmostEqual(math.exp(-possible_mu)*(1+possible_mu+possible_mu**2/2),.05)
            self.assertLess(accepted_mu,possible_mu)
            mu=2*conservative["60"]
            self.assertAlmostEqual(math.exp(-mu)*(1+mu+mu**2/2+mu**3/6),.05)
            self.assertIsNone(conservative["600"])
            readme=(here/"README.md").read_text()
            for evidence_text in ("in 2 dated ObsIDs covering 2.000000 Ms","3 assigned ObsIDs, 2 fully scanned",
                                  "3 of 4 attempted (75.0%)","evidence/injection_recovery.png","evidence/fxt_candidates.csv",
                                  "evidence/reviewed_fxt_triggers.csv","GW170817 field · all-years visits and exclusions",
                                  "ACIS-S4 frame measurement · ObsID 24912, real/0","Unresolved · ObsID 91002, real/0.","The source distance remains unknown.",
                                  "2.000000–2.500000 Ms",receipt[0]["sha256"]):
                self.assertIn(evidence_text,readme)
            self.assertNotIn("GW170817 field · time-averaged rates",readme)
            lines=(repo/"README.md").read_text().splitlines()
            self.assertEqual(lines[8],"3 entries, 1 new measurements, 1 candidate discoveries, 0 forecasts, "
                             "1 corrections and 0 carefully recorded null results.")
            self.assertEqual(lines[9:11],old_entries)
            self.assertEqual(sum(pub.report.TAG in line for line in lines),1)
            self.assertIn("*Candidate*",lines[11])


def campaign_read(path):
    return json.loads(path.read_text())


if __name__=="__main__":
    unittest.main()
