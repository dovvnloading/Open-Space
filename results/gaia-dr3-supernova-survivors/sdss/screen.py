"""Crossmatch SDSS-V SnowWhite star products with Gaia DR3, preserving release overlap."""
import argparse
import csv
import gzip
import hashlib
import json
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.io import fits
from zero_point import zpt

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"desi"))
from fetch import BATCH, digest, fetch_block

SOURCES = {
    "DR19": ("0.6.0","a7489a2f7a0a155fe127badc7c21a5b21eea5791d097300e2f683e0ad6a923a9"),
    "DR20": ("0.8.1","a19c7aa645459d83aabe58915693fa95312bd59b9c9559bd4f0459a57b819099"),
}


def read_survey(cache, rel):
    version, checksum = SOURCES[rel]
    name=f"astraAllStarSnowWhite-{version}.fits.gz"
    url=f"https://{rel.lower()}.sdss.org/sas/{rel.lower()}/spectro/astra/{version}/summary/{name}"
    path=cache/name
    if not path.exists():
        with urlopen(url,timeout=240) as source, path.open("wb") as out:
            shutil.copyfileobj(source,out)
    with path.open("rb") as f:
        if hashlib.file_digest(f,"sha256").hexdigest()!=checksum:
            raise ValueError(f"{rel} digest differs")
    plain=cache/f"{rel.lower()}_snow.fits"
    # Always extract the verified compressed bytes before reading the memory map.
    with gzip.open(path,"rb") as source, plain.open("wb") as out:
        shutil.copyfileobj(source,out,1024**2)
    all_rows = classified = 0
    parts=[]
    with fits.open(plain,memmap=True) as hdus:
        for hdu in hdus[1:]:
            d=hdu.data
            if d is None or len(d)==0:
                continue
            all_rows+=len(d)
            types=np.char.strip(d["classification"])
            good=(types!="")
            classified+=int(good.sum())
            keep=good & (d["gaia_dr3_source_id"]>0)
            parts.append(pd.DataFrame({k:np.asarray(d[k][keep]).astype(str) for k in
                ("gaia_dr3_source_id","classification","sdss_id","telescope","min_mjd","max_mjd","snr")}))
    cat=pd.concat(parts,ignore_index=True).rename(columns={"gaia_dr3_source_id":"source_id"})
    meta=dict(url=url,sha256=checksum,bytes=path.stat().st_size,rows=all_rows,classified_rows=classified,
              identified_rows=len(cat),unique_ids=cat.source_id.nunique())
    return cat,meta


def run(cache):
    cache.mkdir(parents=True,exist_ok=True)
    cats,meta={},{}
    for release in SOURCES:
        cats[release],meta[release]=read_survey(cache,release)
    ids=set().union(*(set(c.source_id) for c in cats.values()))
    prior=json.loads((HERE.parent/"desi"/"receipt.json").read_text())
    gaia=[]
    for b in prior["batches"]:
        p=cache/f"gaia_{b['batch']:03}.csv"
        if digest(p.read_bytes())!=b["response_sha256"]:
            raise ValueError("prior DESI census cache mismatch")
        gaia.append(pd.read_csv(p,dtype={"source_id":str}))
    g=pd.concat(gaia,ignore_index=True)
    missing=sorted(ids-set(g.source_id),key=int)
    fresh=cache/"sdss_gaia";fresh.mkdir(exist_ok=True)
    jobs=[(fresh,i//BATCH,missing[i:i+BATCH]) for i in range(0,len(missing),BATCH)]
    with ThreadPoolExecutor(6) as ex:
        receipts=list(ex.map(fetch_block,jobs))
    if jobs:
        g=pd.concat([g,*[pd.read_csv(fresh/f"gaia_{b['batch']:03}.csv",dtype={"source_id":str})
                         for b in receipts]],ignore_index=True)
    g=g[g.source_id.isin(ids)].copy()
    if len(g)!=len(ids) or g.source_id.duplicated().any():
        raise ValueError("Gaia exact-ID census incomplete")
    zpt.load_tables()
    ok=g.astrometric_params_solved.isin([31,95])&(g.parallax_error>0)&g.pmra.notna()&g.pmdec.notna()
    fields=["phot_g_mean_mag","nu_eff_used_in_astrometry","pseudocolour","ecl_lat","astrometric_params_solved"]
    g["zp_mas"]=np.nan
    g.loc[ok,"zp_mas"]=zpt.get_zpt(*(g.loc[ok,c].to_numpy() for c in fields),_warnings=False)
    p=g.parallax-g.zp_mas
    g["plx_snr"]=p/g.parallax_error
    g["vt_kms"]=4.74047*np.hypot(g.pmra,g.pmdec)/p.where(p>0)
    g["quality"]=(g.ruwe>0)&(g.ruwe<1.4)&(g.visibility_periods_used>=9)&g.duplicated_source.eq(False)&(g.ipd_frac_multi_peak<=2)&(g.ipd_gof_harmonic_amplitude<.1)
    ctrl=pd.read_csv(HERE.parent/"control"/"controls.csv",dtype={"source_id":str})
    reports=[];matches=[];overall=[]
    for rel,c in cats.items():
        c["snr"]=pd.to_numeric(c.snr)
        # Per-source display uses the classification of its highest-S/N spectrum.
        cc=c.sort_values(["snr","telescope"],ascending=[False,True]).drop_duplicates("source_id")
        joined=cc.merge(g,on="source_id",validate="one_to_one")
        matches.append(joined.merge(ctrl,on="source_id").assign(release=rel))
        for sn in [3,5]:
            sub=joined[joined.quality & (joined.plx_snr>=sn)].sort_values("vt_kms",ascending=False)
            reports.append(dict(release=rel,snr_min=sn,eligible=len(sub),over_600=int((sub.vt_kms>600).sum()),
                over_1000=int((sub.vt_kms>1000).sum()),max_vt=round(sub.vt_kms.max(),2),max_id=sub.iloc[0].source_id))
            if sn==3:overall.append(sub.head(max(5, int((sub.vt_kms>600).sum()))).assign(release=rel))
    cols=["release","source_id","sdss_id","classification","snr","telescope","min_mjd","max_mjd",
          "parallax","parallax_error","zp_mas","plx_snr","vt_kms","ruwe","quality"]
    pd.concat(matches)[["name",*cols]].to_csv(HERE/"controls.csv",index=False)
    pd.concat(overall)[cols].to_csv(HERE/"tail.csv",index=False)
    overlap=set(cats["DR19"].source_id)&set(cats["DR20"].source_id)
    summary=dict(sources=meta,dr19_dr20_shared_ids=len(overlap),distinct_ids=len(ids),
        matched_from_desi=len(ids)-len(missing),additional_queried=len(missing),gates=reports)
    (HERE/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (HERE/"gaia_receipt.json").write_text(json.dumps(receipts,indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("cache",type=Path)
    run(p.parse_args().cache)
