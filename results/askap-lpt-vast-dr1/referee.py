#!/usr/bin/env python3
"""Independent integrity / light-curve gates. NOT a pixel-level transient referee."""
import csv,hashlib,json,sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

def sha(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
 return h.hexdigest()
def csvrows(path):
 with path.open(newline="") as f:return list(csv.DictReader(f))
def num(x):return float(x)
def threshold(r):
 c=num(r["avg_compactness"]);f=num(r["max_flux_peak"]);low=num(r["min_flux_peak"])
 return (.75<=c<=1.4 and num(r["n_measurements"])>=5 and num(r["n_selavy"])>=2
  and num(r["n_siblings"])<=0 and num(r["max_snr"])>=8 and f>=5
  and ((low>0 and f/max(low,.5)>=5) or
       (num(r["eta_peak"])>=100 and num(r["v_peak"])>=.8)))
def referee(out):
 out=Path(out);s=json.loads((out/"summary.json").read_text());checks={}
 checks["science_incomplete_not_discovery"]=s.get("claim")=="NO VALIDATED NEW LONG-PERIOD TRANSIENT" and s.get("not_a_period_search") is True
 checks["screen_completed"]=s.get("state")=="catalogue_screen_only_no_periodicity"
 checks["protocol_digest"]=sha(out/"frozen_protocol.json")==s.get("protocol_sha256")
 p=out/"frozen_images.jsonl"
 checks["image_digest"]=p.exists() and sha(p)==s.get("frozen_image_sha256")
 images=[json.loads(x) for x in p.read_text().splitlines()] if p.exists() else []
 image_ids={str(i["image_id"]) for i in images}
 checks["actual_images"]=len(images)==2945 and len(image_ids)==s.get("distinct_images")
 checks["actual_sources"]=s.get("source_rows")==548113
 checks["one_row_publication_discrepancy_documented"]=(s.get("measurement_rows")==6449203
  and s.get("paper_reported_measurement_rows")==6449202 and s.get("measured_minus_paper_count")==1)
 all_sources=csvrows(out/"all_screened_sources.csv");top=csvrows(out/"ranked_sources.csv")
 checks["all_149_preserved"]=len(all_sources)==149 and len(all_sources)==s.get("screen_all")
 checks["top_60_preserved"]=len(top)==60 and len(top)==s.get("ranked_rows")
 checks["distinct_screen_ids"]=len({r["source"] for r in all_sources})==len(all_sources)
 checks["all_thresholds_reproduced"]=all(threshold(r) for r in all_sources)
 checks["top60_same_prefix"]=[r["source"] for r in top]==[r["source"] for r in all_sources[:60]]
 checks["ranking_monotonic"]=all(num(x["screen_score"])>=num(y["screen_score"])
                                   for x,y in zip(all_sources,all_sources[1:]))
 m=csvrows(out/"ranked_epoch_measurements.csv")
 checks["epoch_693_rows"]=len(m)==693 and len(m)==s.get("epoch_review",{}).get("rows")
 checks["epoch_60_sources"]=len({r["source"] for r in m})==60 and {r["source"] for r in m}=={r["source"] for r in top}
 checks["all_epochs_in_frozen_images"]=all(str(r["image_id"]) in image_ids for r in m)
 by=defaultdict(list)
 for r in m:by[r["source"]].append(r)
 patterns=[];hourpairs=[]
 for sid,g in by.items():
  high=[r for r in g if r["forced"]=="False" and r["flag_c4"]=="False"
        and r["has_siblings"]=="False" and num(r["snr"])>=8]
  quiet=[r for r in g if r["forced"]=="True" and r["flag_c4"]=="False"
         and num(r["snr"])<3]
  if len(high)>=2 and quiet:
   patterns.append({"source":sid,"high_detected_epochs":len(high),
     "forced_below_3sigma_epochs":len(quiet),
     "brightest_detected_mjy_beam":round(max(num(r["flux_peak"]) for r in high),3)})
  dates=sorted(datetime.fromisoformat(r["time"]) for r in g)
  if len(dates)>1:
   dt=min((y-x).total_seconds()/60 for x,y in zip(dates,dates[1:]))
   if dt<60:hourpairs.append({"source":sid,"min_separation_minutes":round(dt,3)})
 checks["14_disappearance_patterns"]=len(patterns)==14
 checks["no_false_positive_rate_invented"]=s.get("contaminant_prior",{}).get("false_positive_rate","").startswith("NOT CALIBRATED")
 checks["known_controls_not_covered"]=s.get("controls_in_approximate_footprints")==0
 status=s.get("atnf",{}).get("state","unavailable")
 outdoc={"pass":all(checks.values()),"checks":checks,
  "lightcurve_disappearance_patterns":sorted(patterns,key=lambda r:-r["brightest_detected_mjy_beam"]),
  "under_one_hour_epoch_pairs":hourpairs,
  "external_vetting":{"gaia_queried":sum(z.get("state")=="queried" for z in s.get("gaia",{}).get("checked",[])),
                       "atnf_status":status},
  "independent_review_scope":"Independent CSV/epoch/field consistency and threshold reimplementation; no native image or pulse-time proof",
  "disposition":"No new validated long-period transient; 14 are unvetted light-curve patterns",
  "calibrated_all_sky_null":False,"pending":s.get("pending",[])}
 (out/"referee.json").write_text(json.dumps(outdoc,indent=2,sort_keys=True)+"\n")
 print(json.dumps({"pass":outdoc["pass"],"checks":checks,"patterns":len(patterns),"minute_pairs":len(hourpairs)},indent=2))
 return outdoc["pass"]
if __name__=="__main__":sys.exit(0 if referee(sys.argv[1]) else 1)
