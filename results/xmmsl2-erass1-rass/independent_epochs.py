"""Look for independent pointed observations at the most plausible galaxy positions."""
import argparse
import hashlib
import io
import json
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse

from full_only_review import inspect, plain

TAP = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
RADIUS = 90.0


def tap(query, receipts, limit=2500):
    params = dict(REQUEST="doQuery", LANG="ADQL", FORMAT="votable", MAXREC=limit, QUERY=query)
    r = requests.get(TAP, params=params, timeout=(15, 80))
    receipt = {"endpoint":TAP, "query":query, "http":r.status_code,
               "sha256":hashlib.sha256(r.content).hexdigest(), "bytes":len(r.content)}
    receipts.append(receipt)
    r.raise_for_status()
    result = parse(io.BytesIO(r.content), verify="ignore")
    states = [(i.value,i.content) for res in result.resources for i in res.infos if i.name=="QUERY_STATUS"]
    if any(state in ("ERROR","OVERFLOW") for state, _ in states):
        raise ValueError(str(states)[:250])
    table = result.get_first_table().to_table(use_names_over_ids=True)
    if len(table)>=limit:
        raise ValueError("Output cap reached")
    receipt["rows"] = len(table)
    return table


def extract(table):
    def clean(x):
        if np.ma.is_masked(x):
            return None
        if isinstance(x,bytes):
            return x.decode(errors="replace")
        if isinstance(x,np.generic):
            x=x.item()
        if isinstance(x,float) and not np.isfinite(x):
            return None
        return x
    return [{n:clean(row[n]) for n in table.colnames} for row in table]


def archive(name, positions, receipts):
    cols = extract(tap("SELECT column_name, unit, description FROM TAP_SCHEMA.columns "+
                       f"WHERE table_name='{name}' OR table_name='\"{name}\"'", receipts))
    fields = {r["column_name"]:r for r in cols}
    if not {"ra","dec"}.issubset(fields):
        return {"status":"table unavailable or unrecognized", "columns":cols}, []
    keep = [k for k in fields if k in {"ra","dec","name","time","end_time","obsid","obs_id",
            "observation_id","source_id","detection_id","summary_flag","detection_flag","det_flag",
            "source_number","source_name","pos_err","position_error","extent","exposure","exposure_time",
            "spectrum_flag","var_flag","flux","flux_error","count_rate","count_rate_error"}
            or (k.startswith(("ep_","soft_","med_","medium_","broad_","peak_","mean_","band"))
                and any(s in k for s in ("flux","flag","rate","err","like","extent","expos")))]
    records=[]
    by_field = {}
    projection = ",".join('"'+k+'"' for k in keep)
    for uid, ra, dec in positions:
        sky = SkyCoord(ra*u.deg,dec*u.deg)
        within = (f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},{RADIUS/3600:.8f}))")
        try:
            rows=extract(tap(f'SELECT {projection} FROM "{name}" WHERE {within}',receipts))
            for row in rows:
                pos=SkyCoord(float(row["ra"])*u.deg,float(row["dec"])*u.deg)
                row.update(table=name,detection_uid=uid,separation_arcsec=float(sky.separation(pos).arcsec))
            records.extend(rows)
            by_field[uid]={"status":"ok","rows":len(rows)}
        except Exception as exc:
            by_field[uid]={"status":"incomplete","reason":str(exc)[:300]}
    return {"columns":[fields[k] for k in keep], "cones":by_field},records


def main(out, work):
    reg=pd.read_csv(out/"review_register.csv.gz",low_memory=False)
    cats=pd.read_csv(out/"catalog_sources.csv.gz",dtype={"xmm_obs":str},low_memory=False).set_index("detection_uid",drop=False)
    work=work/"independent"
    work.mkdir(parents=True,exist_ok=True)
    picks=reg.loc[reg.host_type.fillna("").str.fullmatch(r"G|Galaxy|Sy1|Sy2|AGN|QSO|LINER",case=False)
                  & reg.robust_model_floor_ratio.ge(10)
                  & reg.local_galaxy_pcc.le(.01)
                  & reg.stellar_evidence.fillna("").eq("")]
    picks=picks.sort_values("guarded_ratio",ascending=False).drop_duplicates("source_group")
    receipts=[]
    epochs=[]
    status={}
    positions=[(r.detection_uid,r.ra,r.dec) for r in picks.itertuples()]
    for name in ("xmmssc","swift2sxps","swift1sxps"):
        try:
            info, rows=archive(name,positions,receipts)
            status[name]=info
            epochs.extend(rows)
        except Exception as exc:
            status[name]={"status":"incomplete","reason":str(exc)[:350]}
    images=[]
    inputs=json.loads((out/"catalog_inputs.json").read_text())
    for r in picks.itertuples():
        if r.high!="xmm":
            continue
        src=cats.loc[r.detection_uid]
        img=inspect(src,work,receipts,inputs)
        images.append({"name":r.name,"detection_uid":r.detection_uid,"xmm_ra":src.ra,"xmm_dec":src.dec,
                       "net_soft_counts":src.xmm_counts,"exposure_s":src.xmm_exp,"band6_detection_likelihood":src.xmm_ml,
                       "instrument_checks":img})
    result={"selection":"Grouped strong-galaxy follow-up: galaxy/AGN in host register, original model floor at least 10, local galaxy chance probability <=0.01, no decisive stellar veto. Cone centres are refined X-ray positions.",
            "radius_arcsec":RADIUS,"sources":[{"name":r.name,"detection_uid":r.detection_uid,"ra":r.ra,"dec":r.dec}
            for r in picks.itertuples()],"archive_tables":status,"pointed_detections":epochs,"slew_images":images,
            "interpretation":"Every row inside the wide search cone is retained as a possible neighbour. Catalogue flux columns keep their original bands, flags, and units. Empty cones provide no flux limit. Earlier high/low detections alone do not establish prior publication of the change.",
            "input_receipts":receipts}
    (out/"review_independent_epochs.json").write_text(json.dumps(plain(result),indent=2,allow_nan=False,default=str)+"\n")
    print(f"Checked {len(positions)} plausible galaxy fields; {len(epochs)} pointed catalogue rows in search cones")

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    p.add_argument("--work",required=True,type=Path)
    a=p.parse_args()
    main(a.out,a.work)
