"""Verify the bounded samples and their principal science counts from the committed tables."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

BASE=Path(__file__).resolve().parent
def rows(path):
    with (BASE/path).open() as f:return list(csv.DictReader(f))
def info(path):return json.loads((BASE/path).read_text())
def digest(data):return hashlib.sha256(data).hexdigest()
def check(test,message):
    if not test:raise ValueError(message)
def main():
    pilot=rows("pilot/shortlist.csv")
    frozen=info("pilot/frozen_inputs.json")
    for packet in frozen:
        blob=(BASE/"pilot"/packet["file"]).read_bytes()
        check(digest(blob)==packet["gzip_sha256"] and digest(gzip.decompress(blob))==packet["sha256"],packet["file"])
    for r in info("pilot/catwise_manifest.json"):
        p=BASE/"pilot"/r["file"]
        raw=p.read_bytes() if p.exists() else gzip.decompress(Path(str(p)+".gz").read_bytes())
        check(digest(raw)==r["sha256"],r["file"])
    check(len(pilot)==658 and sum(int(r["bright_tier"]) for r in pilot)==3,"pilot selection")
    counters=[sum(int(rows("pilot/inputs/"+f+"_"+k+"_count.csv")[0]["n"]) for f in ("north","south","fornax")) for k in ("all","snr","filtered")]
    check(counters==[3601560,1905431,11299],"pilot counters")
    sc=rows("pilot/untimely/track_scores.csv")
    check({r["source"] for r in sc}=={r["source_name"] for r in pilot},"pilot epoch coverage")
    supported=[]
    for r in sc:
        other=max(int(r["epochs_"+k]) for k in ["stationary","rot90","rot180","rot270"])
        if int(r["epochs_moving"])>=5 and float(r["moving_match_baseline_yr"] or 0)>=5 and int(r["epochs_moving"])-other>=3:supported.append(r["source"])
    check(set(supported)=={"J041358.20-475038.1","J181125.34+665806.4"},"pilot corroboration")
    controls=rows("controls/epoch_matched.csv")
    check(len(controls)==41 and sum(int(r["rows_2p5"])>0 for r in controls)==31 and sum(int(r["bright_rows_2p5"])>0 for r in controls)==30,"control recovery")
    wide=rows("wide/selected.csv");caps=info("wide/inputs_manifest.json")
    check({r["cap"] for r in caps}==set(range(24)) and sum(r["rows"] for r in caps)==2829,"cap coverage")
    for r in caps:
        p=BASE/"wide"/"inputs"/f"cap{r['cap']:02d}.csv"
        check(digest(p.read_bytes())==r["sha256"],p.name)
    check(len(wide)==len({r["source_id"] for r in wide})==1017,"wide rows")
    literature={r["source_name"] for r in rows("catalogues/wide_matches.csv")}
    check(len(literature)==29,"wide literature alerts")
    for p,src in (("counterparts/gaia_pilot/manifest.json",pilot),("counterparts/gaia_wide/manifest.json",wide),("fast/gaia/manifest.json",rows("fast/tracks.csv"))):
        targets=[n for r in info(p) for n in r["targets"]]
        check(len(targets)==len(src) and set(targets)=={r["source_name"] for r in src},p)
    w=rows("fast/w1/track_colors.csv");fast=rows("fast/tracks.csv")
    check(len(fast)==138 and len(w)==138 and sum(r["tile"]!="1342m076" for r in fast)==131,"fast fields")
    check([r["track"] for r in w if r["color_status"]=="sensitive_catalog_dropout"]==["track1342m076_1"],"W1 dropout control")
    f=rows("fast/l1b/summary.csv")
    check([r["source"] for r in f if r["triage"].startswith("track-supported")]==["track1342m076_1"],"NEOWISE fast holdout")
    spectral=[r for name in ("wise0855","shortlist","wide","additional") for r in rows(f"spherex/{name}_fits.csv")]
    check(len(spectral)==1115,"SPHEREx frame count")
    print("Verified 3,601,560 pilot catalogue entries, 24 wide caps, 41 control objects, 138 W2 links, and 1,115 spectral frame measurements.")
if __name__=="__main__":main()
