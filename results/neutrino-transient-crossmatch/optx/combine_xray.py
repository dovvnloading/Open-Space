"""Combine independent EP, Swift-BAT and Swift-XRT episodes. Four joint GRBs are named in the cited EP circulars."""
import csv
from datetime import datetime
from pathlib import Path
from optical import write
from swift import distance

DIR=Path(__file__).resolve().parent
CONFIRMED={"EP250430a":("GRB250430A",40299),"EP250615a":("GRB250615A",40743),
           "EP250706a":("GRB250706A",41005),"EP251230a":("GRB251230A",43259)}
def read(name):
    return list(csv.DictReader((DIR/name).open()))

def main():
    ep=read("ep_events.csv");bat=read("bat_events.csv");swift=read("swift_events.csv")
    aliases=[];reviewed=set()
    for left,right in ((ep,bat),(ep,swift),(bat,swift)):
        for a in left:
            for b in right:
                dt=abs((datetime.fromisoformat(a["utc"])-datetime.fromisoformat(b["utc"])).total_seconds())
                if dt>7*86400:continue
                sep=distance(a["ra_deg"],a["dec_deg"],b["ra_deg"],b["dec_deg"])
                if sep>0.25:continue
                assert CONFIRMED.get(a["source_id"])==(b["source_id"],int(a.get("gcn") or 0)),(a["source_id"],b["source_id"])
                reviewed.add(a["source_id"])
                aliases.append(dict(ep_id=a["source_id"],bat_id=b["source_id"],ep_gcn=a["gcn"],delta_s=round(dt,3),separation_deg=round(sep,6)))
    assert reviewed==set(CONFIRMED)
    write("xray_aliases.csv",aliases)
    combined=[]
    for r in ep+bat+swift:
        if r["source_id"] in CONFIRMED:continue
        combined.append({k:r[k] for k in ("source_id","ra_deg","dec_deg","utc","family","anchor","catalog")})
    assert len(combined)==len({r["source_id"]for r in combined})
    write("xray_events.csv",sorted(combined,key=lambda r:(r["utc"],r["source_id"])))
if __name__=="__main__":main()
