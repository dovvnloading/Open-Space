"""Freeze converged detections, boundaries and receipts for an offline flare replay."""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import shutil
from statistics import median

from fetch_lcr import HERE, eligible

EPOCH=dt.datetime(2001,1,1,tzinfo=dt.timezone.utc)
MJD0=dt.datetime(1858,11,17,tzinfo=dt.timezone.utc)
CUTOFF=dt.datetime.fromisoformat("2026-11-17T12:09:02+00:00")
WEEK=7*86400
# UTC effective dates of the leap seconds elapsed since the Fermi epoch.
LEAP_BOUNDARIES=[(dt.datetime(y,m,1,tzinfo=dt.timezone.utc)-EPOCH).total_seconds()+i+1
                 for i,(y,m) in enumerate([(2006,1),(2009,1),(2012,7),(2015,7),(2017,1)])]
PEAK_FIELDS=["source_id","ra_deg","dec_deg","class","association","peak_mjd","peak_utc","first_bin_center_utc","last_bin_center_utc","bins","photon_flux","median_flux","ts_at_peak","error_at_peak"]

def instant(met):
    leaps=sum(met>=boundary for boundary in LEAP_BOUNDARIES)
    return EPOCH+dt.timedelta(seconds=met-leaps)

def iso(value):
    return value.strftime("%Y-%m-%dT%H:%M:%SZ")

def save_csv(path, rows, fields):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

def keyed(rows, n):
    valid={r[0]:r[1:] for r in rows if len(r)==n}
    if len(valid)!=len(rows):
        raise ValueError("duplicate timestamps or invalid API row length")
    return valid

def extract(raw):
    ts=keyed(raw["ts"],2)
    flux=keyed(raw["flux"],2)
    errors=keyed(raw["flux_error"],3)
    conv=keyed(raw["fit_convergence"],2)
    bins=sorted(t for t in ts if instant(t+WEEK//2)<=CUTOFF)
    selected=[]
    for t in bins:
        if t not in flux:
            continue
        f=flux[t][0]
        lo,hi=errors[t]
        if (conv[t][0]==0 and all(isinstance(v,(int,float)) and math.isfinite(v) for v in (ts[t][0],f,lo,hi))
            and ts[t][0]>=9 and f>0 and lo<f<hi and f/((hi-lo)/2)>=3):
            selected.append((t,ts[t][0],f,lo,hi))
    return bins,selected

def freeze(cache):
    folder=HERE/"weekly_inputs"
    folder.mkdir(exist_ok=True)
    cat=cache/"4FGL.json.gz"
    shutil.copyfile(cat,HERE/"lcr_4fgl.json.gz")
    sources=[s for s in json.loads(gzip.decompress(cat.read_bytes())) if eligible(s)]
    receipts=[{"source_id":"4FGL",**json.loads(cat.with_suffix(".receipt").read_text())}]
    manifest=[]
    for shard,offset in enumerate(range(0,len(sources),25)):
        stream=io.StringIO()
        w=csv.writer(stream,lineterminator="\n")
        w.writerow(["source_id","met_mid","ts","flux","lo","hi"])
        for s in sources[offset:offset+25]:
            name=s["Source_Name"]
            file=cache/(name.replace(" ","_")+".json.gz")
            raw=gzip.decompress(file.read_bytes())
            receipt=json.loads(file.with_suffix(".receipt").read_text())
            if hashlib.sha256(raw).hexdigest()!=receipt["sha256"]:
                raise ValueError(f"response checksum mismatch for {name}")
            receipts.append({"source_id":name,**receipt})
            bins,selected=extract(json.loads(raw))
            w.writerows((name,*row) for row in selected)
            manifest.append(dict(source_id=name,bins_before_cutoff=len(bins),valid_bins=len(selected),
                    first_week_mid_utc=iso(instant(bins[0])) if bins else "",
                    last_week_mid_utc=iso(instant(bins[-1])) if bins else ""))
        (folder/f"valid_{shard:02d}.csv.gz").write_bytes(gzip.compress(stream.getvalue().encode(),mtime=0))
    save_csv(HERE/"lcr_source_counts.csv",manifest,list(manifest[0]))
    save_csv(HERE/"lcr_receipts.csv",receipts,list(receipts[0]))

def flare_groups(bins):
    groups=[]
    for row in bins:
        if not groups or row[0]-groups[-1][-1][0]!=WEEK:
            groups.append([])
        groups[-1].append(row)
    return groups

def replay():
    catalogue=json.loads(gzip.decompress((HERE/"lcr_4fgl.json.gz").read_bytes()))
    sources={s["Source_Name"]:s for s in catalogue if eligible(s)}
    counts={r["source_id"]:r for r in csv.DictReader((HERE/"lcr_source_counts.csv").read_text().splitlines())}
    valid={name:[] for name in sources}
    for path in sorted((HERE/"weekly_inputs").glob("valid_*.csv.gz")):
        with gzip.open(path,"rt") as f:
            for r in csv.DictReader(f):
                valid[r["source_id"]].append((int(r["met_mid"]),float(r["ts"]),float(r["flux"]),float(r["lo"]),float(r["hi"])))
    peaks=[]
    for name,bins in valid.items():
        s=sources[name]
        if len(bins)!=int(counts[name]["valid_bins"]):
            raise ValueError(f"valid-bin denominator mismatch for {name}")
        if len(bins)<10:
            continue
        med=median(row[2] for row in bins)
        for group in flare_groups([r for r in bins if r[2]>=3*med]):
            t,ts,flux,lo,hi=max(group,key=lambda r:r[2])
            epoch=instant(t)
            peaks.append(dict(source_id=name,ra_deg=s["RAJ2000"],dec_deg=s["DEJ2000"],association=s["ASSOC1"],
                        **{"class":s["CLASS1"]},peak_mjd=round((epoch-MJD0).total_seconds()/86400,8),
                        peak_utc=iso(epoch),first_bin_center_utc=iso(instant(group[0][0])),
                        last_bin_center_utc=iso(instant(group[-1][0])),bins=len(group),
                        photon_flux=flux,median_flux=med,ts_at_peak=ts,error_at_peak=(hi-lo)/2))
    save_csv(HERE/"lcr_flare_peaks.csv",peaks,PEAK_FIELDS)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,help="freeze from fetched responses before replay")
    config=p.parse_args()
    if config.cache:
        freeze(config.cache)
    replay()
