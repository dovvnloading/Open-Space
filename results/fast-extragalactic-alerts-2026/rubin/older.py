"""Supplement: screen older broker IDs then test for any earlier positive Rubin detection.

Downloads every positive S/N >=5 row in the primary UTC interval for older IDs.
Screen eligible pairs before looking for pre-window positives in their full histories.
"""
import argparse
import csv
import gzip
import json
import math
import hashlib
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from retrieve import fetch, OFFSET
from screen import QUALITY

START, END = 61345 + OFFSET, 61361 + OFFSET
ROOT = Path(__file__).resolve().parent
FIELDS = ("d.oid,d.measurement_id,d.mjd,d.band,p.visit,d.ra,d.dec,p.psfflux,p.psffluxerr,"
          "p.ssobjectid,p.diaobjectid," + ",".join("p." + key for key in QUALITY))
JOIN = (" FROM alerce_tap.detection AS d JOIN alerce_tap.lsst_detection AS p "
        "ON d.measurement_id=p.measurement_id AND p.oid=d.oid AND p.sid=d.sid")

def fetch_days(directory, days):
    def get_day(day):
        sql = (f"SELECT {FIELDS}{JOIN} JOIN alerce_tap.object AS o ON o.oid=d.oid AND o.sid=d.sid "
               f"WHERE d.sid=1 AND d.mjd>={day+OFFSET:.12f} AND d.mjd<{day+1+OFFSET:.12f} "
               f"AND o.firstmjd<{START:.12f} AND o.n_det>=3 AND p.psffluxerr>0 "
               "AND p.psfflux>=5*p.psffluxerr AND p.isnegative='False'")
        key = f"older_day_{day}"
        if (directory/(key+".csv.gz")).exists():
            return
        fetch((key,sql), directory)
    directory.mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(2) as ex:
        list(ex.map(get_day, days))

def load_rows(directory, key):
    body=gzip.decompress((directory/(key+".csv.gz")).read_bytes())
    receipt=json.loads((directory/(key+".json")).read_text())
    if hashlib.sha256(body).hexdigest()!=receipt["sha256"]:
        raise ValueError(f"digest mismatch for {key}")
    rows=list(csv.DictReader(body.decode().splitlines()))
    if len(rows)!=receipt['rows']:
        raise ValueError(f"row-count mismatch for {key}")
    return rows

def pair_screen(directory):
    groups=defaultdict(dict)
    for day in range(61345,61361):
        for r in load_rows(directory,f"older_day_{day}"):
            group=groups[r["oid"]]
            if group.setdefault(r["measurement_id"],r)!=r:
                raise ValueError("conflicting duplicate")
    pairs=[]
    eligible=0
    sampled=set()
    for oid, measurements in groups.items():
        bands=defaultdict(list)
        for r in measurements.values():
            bands[r["band"]].append(r)
        for rows in bands.values():
            rows.sort(key=lambda r:float(r["mjd"]))
            for i,a in enumerate(rows):
                for b in rows[i+1:]:
                    dt=float(b["mjd"])-float(a["mjd"])
                    if dt>1:break
                    if dt<0.02 or a["visit"]==b["visit"]:continue
                    eligible+=1; sampled.add(oid)
                    f1,f2=float(a["psfflux"]),float(b["psfflux"])
                    dm=2.5*math.log10(f1/f2)
                    err=2.5/math.log(10)*math.hypot(float(a["psffluxerr"])/f1,float(b["psffluxerr"])/f2)
                    if abs(dm)<0.30 or abs(dm)/dt<=1 or abs(dm)<3*err:continue
                    flags=[q for q in QUALITY if any(r[q]=="True" for r in (a,b))]
                    if any(r["ssobjectid"] not in ("","0") or r["diaobjectid"] in ("","0") for r in (a,b)):
                        flags.append("ss_association")
                    pairs.append(dict(oid=oid,band=a["band"],ra=a["ra"],dec=a["dec"],measurement1=a["measurement_id"],
                                 measurement2=b["measurement_id"],t1=a["mjd"],t2=b["mjd"],dm_mag=dm,sigma_mag=err,flags="|".join(flags)))
    counts=dict(positive_rows=sum(map(len,groups.values())),positive_objects=len(groups),
                eligible_pairs=eligible, eligible_objects=len(sampled), passing_pairs=len(pairs),
                passing_objects=len({r["oid"] for r in pairs}))
    (directory/"screen_counts.json").write_text(json.dumps(counts,indent=2)+"\n")
    with (directory/"passing.csv").open("w") as f:
        w=csv.DictWriter(f,lineterminator="\n",fieldnames=["oid","band","ra","dec","measurement1","measurement2","t1","t2","dm_mag","sigma_mag","flags"])
        w.writeheader();w.writerows(sorted(pairs,key=lambda r:(r["oid"],r["t1"],r["t2"])))
    return pairs

def prior(directory):
    pairs=pair_screen(directory)
    ids=sorted({r["oid"] for r in pairs})
    history={}
    for i in range(0,len(ids),100):
        key=f"prior_{i:04d}"
        sql=(f"SELECT d.oid,COUNT(*) AS n,MIN(d.mjd) AS first_positive,MAX(d.mjd) AS last_positive{JOIN} "
            f"WHERE d.sid=1 AND d.oid IN ({','.join(ids[i:i+100])}) AND d.mjd<{START:.12f} "
            "AND p.psfflux>0 AND p.isnegative='False' GROUP BY d.oid")
        if not (directory/(key+".csv.gz")).exists():
            fetch((key,sql), directory)
        for r in load_rows(directory,key): history[r["oid"]]=r
    summary=[]
    for oid in ids:
        source=[r for r in pairs if r["oid"]==oid]
        h=history.get(oid,{})
        summary.append(dict(oid=oid,passing_pairs=len(source),unflagged_pairs=sum(not p["flags"] for p in source),
                       prior_positives=int(h.get("n",0)),first_prior_tai=h.get("first_positive",""),
                       last_prior_tai=h.get("last_positive","")))
    with (directory/"prior_summary.csv").open("w") as f:
        w=csv.DictWriter(f,lineterminator="\n",fieldnames=["oid","passing_pairs","unflagged_pairs","prior_positives","first_prior_tai","last_prior_tai"])
        w.writeheader();w.writerows(summary)
    print(len(summary),"passing IDs;",sum(x["prior_positives"]==0 for x in summary),"first positive in window.")

if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path); p.add_argument("days",nargs="*",type=int)
    p.add_argument("--history",action="store_true")
    a=p.parse_args()
    if a.history:
        prior(a.directory)
    else:
        fetch_days(a.directory, a.days or range(61345,61361))
