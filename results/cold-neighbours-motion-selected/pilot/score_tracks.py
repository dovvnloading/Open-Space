"""Compare predicted motion with a static source and three rotated-motion controls.

These are association diagnostics, not fits to an independent blind detection
sample. Flux/uncertainty and morphology cuts are made before comparing tracks.
"""
import csv
import gzip
import hashlib
import io
import math
import json
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def load(path):
    with path.open() as f: return list(csv.DictReader(f))

def main():
    cats=load(ROOT/"shortlist.csv")
    combined=ROOT/"untimely"/"detections.csv"
    if combined.exists():
        detections=load(combined)
    else:
        packets=ROOT/"untimely"/"neighbourhoods"
        detections=[]
        for packet in json.loads((packets/"manifest.json").read_text()):
            compressed=(packets/(packet["tile"]+".csv.gz")).read_bytes()
            if hashlib.sha256(compressed).hexdigest()!=packet["gzip_sha256"]:raise ValueError(packet["tile"])
            detections.extend(csv.DictReader(io.StringIO(gzip.decompress(compressed).decode())))
    by_source=defaultdict(list)
    for d in detections:
        by_source[d["source"]].append(d)
    manifest=json.loads((ROOT/"untimely"/"manifest.json").read_text())
    observed_tiles={item["url"].rsplit("/",1)[-1][:8] for item in manifest if "url" in item}
    result=[]
    for s in cats:
        if s["source_id"][:8] not in observed_tiles: continue
        ra0=float(s["ra2015_deg"]); dec0=float(s["dec2015_deg"])
        ma=float(s["pmra"]); md=float(s["pmdec"])
        names=["moving","stationary","rot90","rot180","rot270"]
        velocities=[(ma,md),(0,0),(-md,ma),(-ma,-md),(md,-ma)]
        eligible=[x for x in by_source[s["source_name"]] if float(x["dflux"])>0 and
                  float(x["flux"])/float(x["dflux"])>=4 and float(x["qf"])>0.85]
        data={k:{} for k in names}
        for d in eligible:
            mjd=float(d["MJDMEAN"]);t=(mjd-57170)/365.25
            a=(float(d["ra"])-ra0)*3600*math.cos(math.radians(dec0))
            b=(float(d["dec"])-dec0)*3600
            for key,(mua,mud) in zip(names,velocities):
                sep=math.hypot(a-t*mua,b-t*mud)
                if sep < data[key].get(mjd,999): data[key][mjd]=sep
        local=defaultdict(list)
        for d in by_source[s["source_name"]]:
            # Local dflux in the neighborhood is a depth proxy, not a
            # completeness test; retain the original 25-arcsec aperture.
            da=(float(d["ra"])-ra0)*3600*math.cos(math.radians(dec0))
            dd=(float(d["dec"])-dec0)*3600
            if (math.hypot(da,dd)<=25 and float(d["dflux"])>0 and
                float(d["qf"])>.85 and int(d["nm"])>=3):
                local[d["EPOCH"]].append(22.5-2.5*math.log10(4*float(d["dflux"])))
        depths=[statistics.median(v) for v in local.values()]
        result_row=dict(source=s["source_name"],field=s["field"],w2=s["w2"],
                        local_depth_epochs=len(depths),
                        median_local_four_sigma_w2=round(statistics.median(depths),3) if depths else "",
                        local_depth_reaches_w2=sum(v>=float(s["w2"]) for v in depths),
                        bright_tier=s["bright_tier"],mu=s["mu"])
        for key in names:
            times=[mjd for mjd,v in data[key].items() if v<=1.5]
            result_row["epochs_"+key]=len(times)
            if key=="moving":
                result_row["closest_arcsec"]=round(min(data[key].values(),default=math.inf),3)
                result_row["earliest_match_mjd"]=round(min(times),3) if times else ""
                result_row["latest_match_mjd"]=round(max(times),3) if times else ""
                result_row["moving_match_baseline_yr"]=round((max(times)-min(times))/365.25,3) if times else 0
                result_row["matches_after_mjd58465"]=sum(mjd>58465 for mjd in times)
        result.append(result_row)
    with (ROOT/"untimely"/"track_scores.csv").open("w") as out:
        w=csv.DictWriter(out,fieldnames=result[0]);w.writeheader();w.writerows(result)
    for r in result:
        n=r["epochs_moving"];ctrl=max(r["epochs_"+x] for x in ["stationary","rot90","rot180","rot270"])
        if n>=5 and n-ctrl>=3 and r["moving_match_baseline_yr"]>=5:
            print(r["source"],r["w2"],r["mu"],n,ctrl,r["moving_match_baseline_yr"])
if __name__=="__main__":main()
