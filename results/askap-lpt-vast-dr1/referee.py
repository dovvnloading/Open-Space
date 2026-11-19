#!/usr/bin/env python3
"""Independent count/hash assertions; not native-image scientific referee."""
import csv,hashlib,json,sys
from pathlib import Path
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(out):
 out=Path(out);s=json.loads((out/"summary.json").read_text());a={}
 a["complete"]=s.get("state")=="catalogue_screen_only_no_periodicity"
 a["no_discovery"]=s.get("claim")=="NO VALIDATED NEW LONG-PERIOD TRANSIENT"
 a["no_period"]=s.get("not_a_period_search") is True
 a["protocol_hash"]=digest(out/"frozen_protocol.json")==s.get("protocol_sha256")
 a["image_hash"]=(out/"frozen_images.jsonl").exists() and digest(out/"frozen_images.jsonl")==s.get("frozen_image_sha256")
 a["source_rows"]=s.get("source_rows")==548113
 a["measurement_rows"]=s.get("measurement_rows")==6449202
 a["image_rows"]=s.get("frozen_images")==2945
 p=out/"ranked_sources.csv"
 if p.exists():
  with p.open() as f:rows=list(csv.DictReader(f))
  scores=[float(x["screen_score"]) for x in rows]
  a["rank_rows"]=len(rows)==s.get("ranked_rows")
  a["rank_sorted"]=all(x>=y for x,y in zip(scores,scores[1:]))
  a["rank_unique"]=len({x["source"] for x in rows})==len(rows)
 else:a["rank_rows"]=False
 a["no_fpr_claim"]=s.get("contaminant_prior",{}).get("false_positive_rate","").startswith("NOT CALIBRATED")
 r={"pass":all(a.values()),"checks":a,
    "independent_review_scope":"audit of counts, hashes and claims only; native pixel re-reduction NOT done",
    "disposition":"No validated new LPT: variability is not a measured minute-to-hour period",
    "pending":s.get("pending",["data retrieval"])}
 (out/"referee.json").write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
 print(json.dumps({"pass":r["pass"],"checks":a},indent=2))
 return r["pass"]
if __name__=="__main__":sys.exit(0 if check(sys.argv[1]) else 1)
