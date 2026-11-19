#!/usr/bin/env python3
"""Freeze published known-object controls and literature comparison on Actions."""
import csv
import hashlib
import io
import json
import math
import re
import statistics
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
EPOCH = 2024.25  # nominal NISP epoch, not an object-specific observation MJD
SOURCES = {
    "UltracoolSheet_v2.1": (
        "https://zenodo.org/api/records/15802304/files/UltracoolSheet%20-%20Main.csv/content",
        "2307303072ebf729be6e7af6ec67c0a9feed6c2ab6f0adb35444bfab5724a184"),
    "Kirkpatrick2024_20pc": (
        "https://exoplanetarchive.ipac.caltech.edu/data/Contributed/20pcCensus/Kirkpatrick2024_20pcCensus_Table4_Parameters.csv",
        "376b2ad0a681f230043deece0f7e85d3325c0ef2da80a712d7f74361a185941f"),
    "BackyardWorlds_Schneider2026": (
        "https://content.cld.iop.org/journals/1538-3881/171/6/341/revision1/ajae5b60t1_mrt.txt",
        "992b9e186161215499c47c3b39bb7799e35295af0598e025ee603db86de6d0a2")
}

def finite(v):
    try:
        n = float(str(v).strip())
        return n if math.isfinite(n) else None
    except (ValueError, TypeError):
        return None

def number(row, *keys):
    for k in keys:
        v = finite(row.get(k))
        if v is not None:
            return v
    return None

def anycolumn(row, pattern):
    for key in row:
        if re.search(pattern, key, re.I):
            v = finite(row[key])
            if v is not None:
                return v
    return None

def subtype(value):
    m = re.search(r'(?<![A-Za-z])(T[0-9](?:\.[0-9])?|Y[0-9](?:\.[0-9])?)', str(value or "").upper())
    return m.group(1) if m else ""

def cold_type(value):
    s = subtype(value)
    return s.startswith("Y") or s.startswith(("T8", "T9"))

def extended_type(value):
    s = subtype(value)
    return cold_type(s) or s.startswith(("T6", "T7"))

def propagate(ra, dec, epoch, pmra, pmdec):
    if epoch is None or pmra is None or pmdec is None:
        return ra, dec, False
    d = dec + (EPOCH-epoch)*pmdec/3600000
    a = (ra + (EPOCH-epoch)*pmra/(3600000*max(.05,math.cos(math.radians(dec)))))%360
    return a, d, True

def angular_arcsec(a,d,b,e):
    dx=((a-b+180)%360-180)*math.cos(math.radians((d+e)/2))
    return math.hypot(dx,d-e)*3600

def make(source,name,sp,ra,dec,epoch,pmra=None,pmdec=None,plx=None,w2=None,teff=None):
    a,d,proper=propagate(ra,dec,epoch,pmra,pmdec)
    return dict(name=name, source=source, spectral_type=sp, ra=a, dec=d,
                reference_ra=ra, reference_dec=dec, reference_epoch=epoch,
                motion_propagated=proper, pmra_masyr=pmra, pmdec_masyr=pmdec,
                parallax_mas=plx, w2_vega=w2, teff_published_K=teff)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    essential=[OUT/"known_population.json.gz",OUT/"control_catalog.json",OUT/"literature_receipts.json"]
    if all(x.exists() for x in essential):
        json.loads((OUT/"control_catalog.json").read_text())
        print("CONTROLS_FROZEN",flush=True)
        return
    raw,receipts={},{}
    for key,(url,expected) in SOURCES.items():
        last=None
        for attempt in range(3):
            try:
                response=requests.get(url,timeout=(12,120),
                    headers={"User-Agent":"Open-Space-Euclid-DR1F-TY-controls/1.0"})
                response.raise_for_status()
                data=response.content
                digest=hashlib.sha256(data).hexdigest()
                if digest!=expected:
                    raise ValueError(f"Frozen {key} changed: {digest}")
                raw[key]=data
                receipts[key]={"url":url,"sha256":digest,"bytes":len(data),
                    "retrieved_utc":datetime.now(timezone.utc).isoformat()}
                break
            except requests.RequestException as ex:
                last=ex
        if key not in raw:
            raise RuntimeError(f"Cannot download frozen control source {key}: {last}")
    import gzip
    literature=[]
    uc=list(csv.DictReader(io.StringIO(raw["UltracoolSheet_v2.1"].decode("utf-8-sig"))))
    for row in uc:
        ra=number(row,"ra_j2000_formula");dec=number(row,"dec_j2000_formula")
        if ra is None or dec is None:continue
        name=row.get("name","").strip()
        sp=subtype(row.get("spt_ir") or row.get("spt_opt"))
        pmra=number(row,"pmra_formula")
        pmdec=number(row,"pmdec_formula")
        plx=number(row,"plx_formula")
        w2=number(row,"w2mpro_formula","w2_formula","w2mpro","W2")
        if w2 is None:w2=anycolumn(row,r"^(?:w2|w2mpro).*formula$")
        literature.append(make("UltracoolSheet_v2.1",name,sp,ra,dec,2000.0,
                               pmra,pmdec,plx,w2))
    census=list(csv.DictReader(io.StringIO(raw["Kirkpatrick2024_20pc"].decode("utf-8-sig"))))
    for row in census:
        ra=number(row,"RA (deg)");dec=number(row,"Dec (deg)")
        if ra is None or dec is None:continue
        ep=number(row,"Epoch")
        pmra=number(row,"PMRA (mas/yr)")
        pmdec=number(row,"PMDec (mas/yr)")
        plx=number(row,"Parallax (mas)")
        sp=""
        for key in row:
            if "spectral" in key.lower() or "spt" in key.lower():
                sp=subtype(row[key])
                if sp:break
        teff=anycolumn(row,r"(?:^teff|effective.*temp|t_eff)")
        w2=number(row,"W2 (mag)","W2")
        if w2 is None: w2=anycolumn(row,r"^W2(?:\s|$|\()")
        literature.append(make("Kirkpatrick2024_20pc",row.get("DefaultName",""),sp,
                               ra,dec,ep,pmra,pmdec,plx,w2,teff))
    # Schneider et al. published motion-selected IDs; name-derived coordinates
    # have *no adopted proper motion* and are positional alerts, not exclusions.
    p=re.compile(r'^J(\d\d)(\d\d)(\d\d\.\d+)([+-])(\d\d)(\d\d)(\d\d\.\d+)')
    for line in raw["BackyardWorlds_Schneider2026"].decode("utf-8","replace").splitlines():
        m=p.match(line)
        if not m:continue
        h,mi,s,sign,d,dm,ds=m.groups()
        ra=15*(int(h)+int(mi)/60+float(s)/3600)
        dec=(1 if sign=="+" else -1)*(int(d)+int(dm)/60+float(ds)/3600)
        literature.append(make("BackyardWorlds_Schneider2026",m.group(), "",ra,dec,2015.405))
    # Prior cold-neighbour work did not survey Euclid DR1-F HATS; its dated
    # publications are reference sets only, not measurements reused here.
    tcontrols=[r for r in literature if extended_type(r["spectral_type"]) and
               r["source"]=="UltracoolSheet_v2.1"]
    for r in literature:
        if r["source"]!="Kirkpatrick2024_20pc" or not extended_type(r["spectral_type"]):
            continue
        if any(angular_arcsec(r["ra"],r["dec"],x["ra"],x["dec"])<3.0 for x in tcontrols):
            continue
        tcontrols.append(r)
    controls=[]
    for i,r in enumerate(tcontrols):
        controls.append({**r,"name":"COLDCTRL:"+str(i)+":"+r["name"],"class":"cold_dwarf" if cold_type(r["spectral_type"]) else "warm_T_control",
                         "z":None})
    # Frozen spectroscopic quasars are an explicit contaminant control.
    import sys
    sys.path.insert(0,str(ROOT.parent/"euclid-dr1-wide-highz-vetting"))
    import scan as upstream
    for r in upstream.load_controls():
        if r["class"]=="quasar":
            controls.append({**r,"name":"QSOCTRL:"+r["name"],"source":r["source"],
                             "spectral_type":"","motion_propagated":False})
    # Parallax based empirical M_W2 references. Do not manufacture missing
    # photometry or extrapolate beyond published classified controls.
    anchors={}
    for group in ("T8","T9","Y0","Y1","Y2"):
        rows=[]
        for r in literature:
            if not r["spectral_type"].startswith(group):continue
            pmas,w2=r["parallax_mas"],r["w2_vega"]
            if pmas is None or not 10<pmas<1000 or w2 is None or not 5<w2<23:continue
            mw2=w2-5*math.log10(1000/pmas/10)
            if 5<mw2<24:rows.append((mw2,r.get("teff_published_K")))
        if len(rows)>=3:
            vals=[x[0] for x in rows]
            anchors[group]={"n":len(vals),"median_M_W2":round(statistics.median(vals),3),
                 "scatter_M_W2_MAD":round(statistics.median(abs(v-statistics.median(vals)) for v in vals)*1.4826,3),
                 "median_published_Teff_K":statistics.median([x[1] for x in rows if x[1] and 150<x[1]<1800]) if any(x[1] and 150<x[1]<1800 for x in rows) else None}
    with gzip.open(OUT/"known_population.json.gz","wt",encoding="utf8") as handle:
        json.dump(literature,handle,separators=(",",":"),allow_nan=False)
    (OUT/"control_catalog.json").write_text(json.dumps(controls,indent=2,allow_nan=False)+"\n")
    (OUT/"reference_calibration.json").write_text(json.dumps({"anchors":anchors,
        "method":"Median published-parallax-derived absolute W2 Vega magnitudes by literature subtype; no invented relation",
        "no_Euclid_bandpass_Teff_relation_claimed":True},indent=2)+"\n")
    receipts["counts"]={"ultracool":len(uc),"census":len(census),
        "literature_positions":len(literature),"controls_cold":sum(x["class"]=="cold_dwarf" for x in controls),
        "controls_warm":sum(x["class"]=="warm_T_control" for x in controls),
        "controls_quasar":sum(x["class"]=="quasar" for x in controls)}
    (OUT/"literature_receipts.json").write_text(json.dumps(receipts,indent=2)+"\n")
    print("CONTROLS",json.dumps(receipts["counts"]),flush=True)

if __name__=="__main__":
    main()
