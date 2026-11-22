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


class OfflineFieldReplay(unittest.TestCase):
    def test_archive_members_positions_dates_and_unlocalized(self):
        import io,gzip,json,tarfile
        from datetime import datetime,timezone
        from astropy.time import Time
        import alerts,field_replay as replay
        peak=datetime(2024,2,1,12,tzinfo=timezone.utc)
        lead={"obsid":1,"input_sha256":"a"*64,"trigger_key":"real/0","ra_deg":13.5,"dec_deg":-2.25,
              "event_mjd":float(Time(peak).mjd),"trial_duration_s":60}
        docs=[{"circularId":1,"subject":"GRB 240201A, location","body":"RA = 13.5 deg, DEC = -2.25 deg. T0=2024-02-01T12:00:03 UTC",
               "createdOn":peak.timestamp()*1000},{"circularId":2,"subject":"Report without a location","body":"No labelled coordinates.",
               "createdOn":peak.timestamp()*1000}]
        with TemporaryDirectory() as directory:
            root=Path(directory)
            src,out=root/"source",root/"out"
            src.mkdir()
            blob=io.BytesIO()
            matches=[]
            with tarfile.open(fileobj=blob,mode="w:gz") as z:
                for doc in docs:
                    body=json.dumps(doc).encode()
                    h=replay.digest(body)
                    ident=doc["circularId"]
                    name=f"2024/{ident}.json"
                    f=tarfile.TarInfo(name)
                    f.size=len(body)
                    z.addfile(f,io.BytesIO(body))
                    row=alerts.review(lead,peak,alerts.dates_in_name(" "+doc["subject"]),peak,False,
                                      list(alerts.positions(doc["body"])),alerts.onset_times(doc["body"]))
                    matches.append({**row,"file":f"members/{ident}-{h[:16]}.json","archive_member":name,"sha256":h,
                                    "id":ident,"subject":doc["subject"],"event":"","published_utc":peak.isoformat()})
            archive_url="https://gcn.nasa.gov/circulars/fixture.tar.gz"
            memo={}
            for label,url,body in (("archive",archive_url,blob.getvalue()),("index",alerts.INDEX,
                                  f'<a href="{archive_url}">JSON (31 MB)</a>'.encode())):
                stem=label+"-"+replay.digest(url.encode())[:24]
                h=replay.digest(body)
                name=stem+"-"+h[:20]+".bin"
                (src/name).write_bytes(body)
                attempt={"http_status":200,"body_file":name,"sha256":h,"size_bytes":len(body)}
                memo[label]={**attempt,"ok":True,"url":url,"attempts":[attempt]}
                ingest.save(src/(stem+".json"),memo[label])
            review={"receipts":memo,"leads":[{**lead,"matches":matches,"status":"review_shortlist"}]}
            bodies=replay.restore(src,out,review,lambda p:p,memo["archive"]["sha256"])
            derived=replay.circulars(review,lead,bodies)
            self.assertEqual([r["position_status"] for r in derived["notices"]],["parsed","unlocalized_by_parser"])
            r=derived["notices"][0]
            self.assertAlmostEqual(r["positions"][0]["offset_arcsec"],0)
            self.assertAlmostEqual(r["explicit_onsets"][0]["seconds_from_peak"],3,places=4)
            self.assertEqual(r["event_name_days"][0]["calendar_days_from_peak"],0)
            replay.copy_members(out,root/"future",review)
            for k,b in bodies.items():
                self.assertEqual((root/"future"/k).read_bytes(),b)
            (out/matches[0]["file"]).write_bytes(b"different")
            with self.assertRaisesRegex(ValueError,"Different bytes"):
                replay.restore(src,out,review,lambda p:p,memo["archive"]["sha256"])

    def test_field_visit_receipts_epoch_sums_multiplicity_and_exclusion(self):
        import gzip
        from scipy.stats import chi2
        import campaign,field_replay as r,visits
        lead={"obsid":1,"input_sha256":"b"*64,"trigger_key":"real/0","ra_deg":13.5,"dec_deg":0.,
              "position_error_arcsec":2.,"event_mjd":60000.}
        key=visits._lead_id(lead)
        with TemporaryDirectory() as d:
            base=Path(d); history=base/"history"
            supported=[]; gaps=[]; sources=[]
            for oid,epoch,lo,gross,bg in [(2,"before",59990.,3,12),(3,"after",60001.,0,24),(4,"gap",59999.,0,0)]:
                start=0.;stop=100.
                rel=f"visits/{oid}/visit.json"
                p=history/rel
                p.parent.mkdir(parents=True)
                name=f"acisf{oid:05d}N001_evt2.fits.gz"
                url=f"https://cda.harvard.edu/byobsid/{oid%10}/{oid}/primary/"
                html=f'<a href="{name}">file</a>'.encode()
                (p.parent/"primary_0.html.gz").write_bytes(gzip.compress(html))
                rec={"filename":name,"obsid":oid,"revision":1,"sha256":f"{oid:064x}","fits_sha256":f"{oid+10:064x}",
                     "url":url+name,"index_url":url,"index_sha256":r.digest(html),"bytes":420}
                ref={"obsid":oid,"record":rel,"filename":name,"input_sha256":rec["sha256"]}
                retrieved={k:rec[k] for k in ("filename","url","bytes","sha256")}
                if epoch=="gap":
                    m={"status":"coverage_unconfirmed"}
                    gaps.append({**ref,**m})
                else:
                    m={"status":"measured_aperture","accepted_live_s":90.,"annulus_area_ratio":12.,"annulus_counts":bg,
                       "gross_counts":gross,"event_mjd_utc":60000.,"visit_mjd_utc":[lo,lo+.001],"temporal_relation":epoch,
                       "net_aperture_rate_ct_s":(gross-bg/12)/90.,"local_excess_p":.2,
                       "conditional_aperture_rate_upper95_ct_s":float(.5*chi2.ppf(.95,2*(gross+1)))/90.,
                       "visit_iso_utc":["start","end"],"segments":[{"accepted_bins":[[start,stop]],"duty":.9,
                                                                    "gross_counts":gross,"annulus_counts":bg}]}
                    stamp={k:m[k] for k in ("status","visit_iso_utc","temporal_relation","accepted_live_s","gross_counts",
                                    "annulus_counts","net_aperture_rate_ct_s","local_excess_p","conditional_aperture_rate_upper95_ct_s")}
                    supported.append({**ref,**stamp,"per_target_bonferroni_p":.4})
                sources.append({"obsid":oid,"record":rel,"status":"acquired"})
                ingest.save(p,{"obsid":oid,"lead_ids":[key],"status":"acquired","products":[{"input":rec,
                            "status":"measured_or_geometry_excluded","measurements":{key:m}}],"archive_attempts":[
                            {"index_url":url,"index":{"sha256":r.digest(html),"bytes":len(html)},"selected_files":[name],"downloaded":[retrieved]}]})
            ids=[2,3,4]
            ingest.save(history/"target-00/selection.json",{"included":{str(i):[] for i in ids},"excluded":[{"obsid":1,"reason":"own trigger"}]})
            ingest.save(history/"history.json",{"selected_leads":[lead],"targets":{key:{"folder":"target-00","pointing_query_status":"ok",
                       "nearby_obsids":ids,"supported":supported,"gaps":gaps}},"visits":sources})
            execution={"visit_coverage":{"pointing_query_status":"ok","nearby_obsids":3,"supported_products":2,"coverage_gaps":1}}
            result=r.visit_records(base,lead,execution,lambda p:p)
            self.assertEqual(result["coverage_gaps"],gaps)
            self.assertEqual(result["per_target_trials"],2)
            self.assertEqual(result["pointing_exclusions"][0]["obsid"],1)
            b,a=result["summary_by_epoch"]["before"],result["summary_by_epoch"]["after"]
            self.assertEqual((b["summed_live_s"],b["gross_counts"],b["scaled_annulus_counts"]),(90.,3,1.))
            self.assertAlmostEqual(a["mean_net_aperture_rate_ct_s"],-2/90)
            self.assertAlmostEqual(a["gross_poisson_upper95_ct_s"],2.9957322735539895/90)
            self.assertIsNone(result["summary_by_epoch"]["overlapping_epoch"]["gross_poisson_upper95_ct_s"])
            p=history/"visits/2/primary_0.html.gz"
            p.write_bytes(gzip.compress(b"another official-looking index"))
            with self.assertRaisesRegex(ValueError,"index changed"):
                r.visit_records(base,lead,execution,lambda p:p)
