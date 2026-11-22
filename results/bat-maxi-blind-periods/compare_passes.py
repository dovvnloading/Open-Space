"""Summarize coverage changes relative to a previous frozen-input search."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path


def read(folder):
    found = {}
    for f in folder.glob("search-*.jsonl.gz"):
        with gzip.open(f,"rt") as handle:
            for line in handle:
                r=json.loads(line)
                found[r["roster_id"]]=r
    return found


def passes(previous,current):
    old,new=read(previous),read(current)
    if not old:
        return
    if old.keys()!=new.keys():
        raise ValueError("Cannot compare passes on different rosters")
    summary={"previous_status":Counter(),"current_status":Counter(),"source_transitions":Counter(),"band_transitions":{}}
    changes=[]
    for sid in sorted(old):
        x,y=old[sid],new[sid]
        summary["previous_status"][x["status"]]+=1
        summary["current_status"][y["status"]]+=1
        summary["source_transitions"][x["status"]+" / "+y["status"]]+=1
        for band in ("bat","maxi"):
            b1=x.get("instruments",{}).get(band,{})
            b2=y.get("instruments",{}).get(band,{})
            label = str(b1.get("status"))+" / "+str(b2.get("status"))
            summary["band_transitions"].setdefault(band,Counter())[label]+=1
            if b1.get("status")!=b2.get("status"):
                changes.append(dict(id=sid,band=band,before=b1.get("status"),after=b2.get("status"),
                                    samples_before=b1.get("samples"),samples_after=b2.get("samples"),
                                    effective_before=b1.get("effective_samples"),effective_after=b2.get("effective_samples")))
    summary["earlier_stable"]=[(s,c["period_days"]) for s,r in old.items() for c in r.get("candidates",[])]
    summary["current_stable"]=[(s,c["period_days"]) for s,r in new.items() for c in r.get("candidates",[])]
    (current/"pass-comparison.json").write_text(json.dumps(summary,indent=2)+"\n")
    with (current/"changed-coverage.jsonl").open("w") as target:
        for row in changes:
            target.write(json.dumps(row)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("previous",type=Path)
    parser.add_argument("current",type=Path)
    args=parser.parse_args()
    passes(args.previous,args.current)
