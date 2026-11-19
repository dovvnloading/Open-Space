#!/usr/bin/env python3
"""Offline accounting referee, NOT a substitute for independent astronomer pixel review."""
import csv,json,sys
from collections import Counter
from pathlib import Path
d=Path("results/askap-lpt-vast-dr1-vetting/output")
if not (d/"summary.json").exists():
 print("No completed evidence; cannot certify.",file=sys.stderr);sys.exit(1)
s=json.loads((d/"summary.json").read_text())
def read(n):return list(csv.DictReader((d/n).open()))
a=read("classifications.csv");e=read("epoch_image_audit.csv")
assert len(a)==60 and len(set(x["source"] for x in a))==60
assert len(e)==693
assert sum(x["onoff"]=="True" for x in a)==14
assert sum(bool(x["psr"]) for x in a)==15
assert all(x["RACS_V"]=="not_measured" and x["VAST_V"]=="not_measured" for x in a)
assert not s["discovery_validated"]
assert s["native_images"]==sum(x["image_state"]=="audited" for x in e)
assert dict(Counter(x["class"] for x in a))==s["classes"]
out={"status":"PASS accounting/anti-overclaim ONLY","n":60,"epochs":693,"onoff":14,
     "known_ATNF_matches":15,"native_epoch_audits":s["native_images"],
     "classes":s["classes"],"stokes_V_tested":False,
     "periodic_LPT_confirmed":False,"independent_visual_referee":"not established"}
(d/"referee.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print(json.dumps(out,indent=2))
