"""Final-release single-exposure follow-up of every literature-unmatched bright wide source.
python screen.py --selected ../../wide/selected.csv --matches ../../catalogues/wide_matches.csv
Run again as further caps arrive. Each target's native archive response is cached.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.coordinates import get_body_barycentric
from astropy.table import Table
from astropy.time import Time

BASE=Path(__file__).resolve().parent
COLS="ra,dec,mjd,w2mpro,w2sigmpro,w2snr,w2rchi2,w1mpro,w1snr,sso_flg,cc_flags,scan_id,frame_num,sigra,sigdec,qual_frame"
CONTROL="J085715.68+560405.8"

def fit_epochs(hit,a,d):
    # Assign points by the input track. Centroids and this diagnostic fit are consequently prior-selected.
    hit=hit.sort_values("mjd").copy()
    hit["visit"]=(hit.mjd.diff()>30).cumsum()
    e=[]
    for _,g in hit.groupby("visit"):
        if len(g)<2:continue
        item=dict(mjd=float(g.mjd.median()),n=len(g))
        for col in ("x","y"):
            med=float(g[col].median())
            item[col]=med
            item["s"+col]=float(np.hypot(.15,1.2533*1.4826*np.median(abs(g[col]-med))/np.sqrt(len(g))))
        e.append(item)
    if len(e)<6 or e[-1]["mjd"]-e[0]["mjd"]<4*365.25:return e,None
    g=pd.DataFrame(e)
    xx,yy,zz=get_body_barycentric("earth",Time(g.mjd,format="mjd",scale="utc")).xyz.to_value("AU")
    ar,dr=np.deg2rad([a,d])
    px=xx*np.sin(ar)-yy*np.cos(ar)
    py=xx*np.cos(ar)*np.sin(dr)+yy*np.sin(ar)*np.sin(dr)-zz*np.cos(dr)
    n=len(g)
    A=np.zeros((2*n,5));A[:n,0]=A[n:,1]=1
    A[:n,2]=A[n:,3]=(g.mjd-58500.)/365.25
    A[:,4]=np.r_[px,py]
    sig=np.r_[g.sx,g.sy];xy=np.r_[g.x,g.y];aw=A/sig[:,None]
    coef,_,rank,_=np.linalg.lstsq(aw,xy/sig,rcond=None)
    if rank<5:return e,None
    chi=np.sum(((xy-A@coef)/sig)**2)
    err=np.sqrt(np.diag(np.linalg.inv(aw.T@aw))*max(1,chi/(2*n-5)))
    return e,dict(parameters_arcsec=coef.tolist(),formal_errors_arcsec=err.tolist(),chi2=float(chi),dof=2*n-5)

def one(row, heldout_mjd=58484.0):
    name=row["source_name"]
    a,d,mx,my,w2=[float(row[k]) for k in ("ra2015_deg","dec2015_deg","pmra","pmdec","w2")]
    radius=math.ceil(math.hypot(mx,my)*10+14)
    url="https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?"+urlencode(dict(
        catalog="neowiser_p1bs_psd",spatial="Cone",objstr=f"{a:.7f} {d:.7f}",radius=radius,
        radunits="arcsec",selcols=COLS,outfmt=1))
    rawpath=BASE/"cache"/(name+".tbl.gz")
    result=BASE/"retained"/(name+".json")
    code_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if result.exists():
        previous=json.loads(result.read_text())
        if (previous.get("code_sha256")==code_hash and previous.get("heldout_mjd")==heldout_mjd and previous["ra"]==a and previous["dec"]==d
            and previous["pmra"]==mx and previous["pmdec"]==my and previous["w2_catwise"]==w2):
            return previous["summary"]
    if not rawpath.exists():
        for i in range(3):
            try:
                b=urlopen(url,timeout=150).read()
                if b"RowsRetrieved" not in b:raise ValueError(b[:200])
                rawpath.write_bytes(gzip.compress(b,mtime=0))
                break
            except (TimeoutError,OSError,ValueError):
                if i==2:raise
                time.sleep(2**i)
    b=gzip.decompress(rawpath.read_bytes())
    t=Table.read(b.decode().splitlines(),format="ascii.ipac").to_pandas()
    advertised=int(next(l.split("=")[1] for l in b.decode().splitlines() if l.startswith("\\RowsRetrieved")))
    if len(t)!=advertised:raise ValueError("Row count mismatch")
    t["x"]=(t.ra-a)*np.cos(np.deg2rad(d))*3600;t["y"]=(t.dec-d)*3600
    dt=(t.mjd-57170.)/365.25
    clean=(t.w2snr>=5)&(t.w2rchi2<3)&(t.sso_flg==0)&(t.qual_frame>=5)&t.cc_flags.astype(str).str.zfill(4).str[1].eq("0")
    late=t.mjd>=heldout_mjd
    masks={}
    for label,u,v in [("forward",mx,my),("stationary",0,0),("rot90",-my,mx),("reverse",-mx,-my),("rot270",my,-mx)]:
        dist=np.hypot(t.x-u*dt,t.y-v*dt)
        masks[label]=dist
    counts={}
    for label,dist in masks.items():
        z=t[clean&late&(dist<=2)].sort_values("mjd")
        counts[label]=dict(n=len(z),visits=int((z.mjd.diff()>30).sum()+1) if len(z) else 0,
                           w2median=float(z.w2mpro.median()) if len(z) else None)
    hit=t[clean&(masks["forward"]<=2)]
    epochs,fit=fit_epochs(hit,a,d)
    near=(masks["forward"]<=5)&late
    n3=int(np.sum(near&(t.w2snr>=3)))
    n5=int(np.sum(near&(t.w2snr>=5)))
    f=counts["forward"]
    compare=max(counts[z]["n"] for z in ("rot90","reverse","rot270","stationary"))
    # Operational triage, insufficient for discovery. Offset tracks have correlated footprints.
    if f["visits"]>=4 and f["n"]>=8 and f["n"]>2*max(1,compare):
        outcome="track-supported; inspect images"
    elif f["visits"]>=3 and f["n"]>=6:
        outcome="ambiguous track or stationary alternative"
    elif n3==0:
        outcome="zero held-out W2 S/N>=3 within 5 arcsec"
    else:outcome="sparse or inconclusive"
    retained=t[clean&np.logical_or.reduce([dist<5. for dist in masks.values()])]
    table_path=BASE/"retained"/(name+".csv.gz")
    retained.to_csv(table_path,index=False,compression={"method":"gzip","mtime":0})
    record=dict(source=name,role="known T8 control" if name==CONTROL else "unmatched selected source",
                track_ref_mjd=57170.,heldout_mjd=heldout_mjd,ra=a,dec=d,pmra=mx,pmdec=my,w2_catwise=w2,
                archive_url=url,cone_arcsec=radius,n_field=len(t),n_clean=int(clean.sum()),n_field_late=int(late.sum()),
                field_mjd=[float(t.mjd.min()),float(t.mjd.max())] if len(t) else [],
                trace_n3_5arcsec_heldout=n3,trace_n5_5arcsec_heldout=n5,heldout_2arcsec=counts,
                quality_selected_epochs=epochs,prior_selected_five_parameter_fit=fit,
                raw_ipac_sha256=hashlib.sha256(b).hexdigest(),
                raw_science_rows_sha256=hashlib.sha256(b"\n".join(l for l in b.splitlines() if not l.startswith(b"\\"))).hexdigest(),
                retained_csv_sha256=hashlib.sha256(gzip.decompress(table_path.read_bytes())).hexdigest(),triage=outcome)
    summary=dict(source=name,role=record["role"],w2=w2,n_field_late=record["n_field_late"],
                n_w2_3_within5=n3,n_w2_5_within5=n5,n_clean_within2=f["n"],visits=f["visits"],
                n_stationary=counts["stationary"]["n"],max_alt=compare,triage=outcome)
    record.update(summary=summary,code_sha256=code_hash)
    result.write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    print(name,outcome,f["n"],f["visits"],compare,flush=True)
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--selected",type=Path,required=True)
    p.add_argument("--heldout-mjd",type=float,default=58484.0);p.add_argument("--matches",type=Path,required=True);p.add_argument("--control",action="store_true");args=p.parse_args()
    (BASE/"retained").mkdir(exist_ok=True);(BASE/"cache").mkdir(exist_ok=True)
    seen={r["source_name"] for r in csv.DictReader(args.matches.open())}
    targets=[r for r in csv.DictReader(args.selected.open()) if r["source_name"] not in seen or (args.control and r["source_name"]==CONTROL)]
    with ThreadPoolExecutor(4) as ex:records=list(ex.map(lambda row: one(row,args.heldout_mjd),targets))
    pd.DataFrame(records).to_csv(BASE/"summary.csv",index=False)
if __name__=="__main__":main()
