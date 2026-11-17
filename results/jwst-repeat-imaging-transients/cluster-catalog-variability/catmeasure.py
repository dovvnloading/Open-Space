"""Bright persistent compact-object variability in epoch-specific STScI catalogs.
Both epochs must have detections. This selection does not test dropouts.
"""
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from catfetch import obs

HERE=Path(__file__).resolve().parent

def sources(epoch,band):
    path=HERE/"catalog-cache"/(epoch+"_nircam_clear-"+band+".ecsv")
    d=pd.read_csv(path,sep=r"\s+",comment="#",usecols=["label","sky_centroid.ra","sky_centroid.dec","aper_total_flux","aper_total_flux_err","is_extended"])
    return d[np.isfinite(d["sky_centroid.ra"])&np.isfinite(d["sky_centroid.dec"])].copy().reset_index(drop=True)

def pairs(a,b,radius):
    t=cKDTree(b); dist,j=t.query(a); i=cKDTree(a).query(b)[1]
    ii=np.where((dist<radius)&(i[j]==np.arange(len(a))))[0]
    return ii,j[ii],dist[ii]

def mad(x):
    return float(1.4826*np.median(np.abs(x-np.median(x))))

def xy(d,origin):
    a=d[["sky_centroid.ra","sky_centroid.dec"]].to_numpy()-origin
    return a*np.array([np.cos(np.deg2rad(origin[1])),1])*3600

def measure(field,epochs,band):
    a,b=[sources(ep,band) for ep in epochs]
    origin=a[["sky_centroid.ra","sky_centroid.dec"]].median().to_numpy()
    x,y=xy(a,origin),xy(b,origin)
    ia,ib,_=pairs(x,y,.5)
    offset=np.median(y[ib]-x[ia],axis=0)
    ia,ib,sep=pairs(x,y-offset,.12)
    aa,bb=a.iloc[ia],b.iloc[ib]
    f,g=[d.aper_total_flux.to_numpy()*1e9 for d in (aa,bb)]
    ef,eg=[d.aper_total_flux_err.to_numpy()*1e9 for d in (aa,bb)]
    comp=(~aa.is_extended.to_numpy())&(~bb.is_extended.to_numpy())
    bright=(ef>0)&(eg>0)&(f>=20*ef)&(g>=20*eg)&comp
    usable=(ef>0)&(eg>0)&(f>=5*ef)&(g>=5*eg)&comp
    if bright.sum()<5:raise ValueError("Insufficient compact controls")
    ratio=g[bright]/f[bright]
    gain=float(np.median(ratio)); floor=mad(np.log(ratio))
    g=g/gain;eg=eg/gain
    sigma=(g-f)/np.sqrt(eg**2+ef**2+(floor*np.minimum(f,g))**2)
    cut=usable&(np.abs(sigma)>=7)&((f>=2*g)|(g>=2*f))
    rows=[]
    for k in np.flatnonzero(cut):
        rows.append(dict(field=field,band=band,label0=int(aa.iloc[k].label),label1=int(bb.iloc[k].label),ra=float(aa.iloc[k]["sky_centroid.ra"]),dec=float(aa.iloc[k]["sky_centroid.dec"]),sep_arcsec=float(sep[k]),epoch0_njy=float(f[k]),epoch1_gain_corrected_njy=float(g[k]),formal_err0_njy=float(ef[k]),formal_err1_corrected_njy=float(eg[k]),sigma=float(sigma[k])))
    mag=lambda n:float(31.4-2.5*np.log10(n))
    re=dict(field=field,band=band,epoch0=epochs[0],epoch1=epochs[1],catalog_counts=[len(a),len(b)],offset_east_north_arcsec=offset.tolist(),mutual_matches=len(ia),median_pair_separation_arcsec=float(np.median(sep)),both_epoch_compact_5sigma=int(usable.sum()),compact_20sigma_controls=int(bright.sum()),relative_gain=gain,empirical_ln_flux_ratio_scatter=floor,flags=len(rows),median_catalog_5sigma_ab_at_eligible_positions=[mag(5*np.median(ef[usable])),mag(5*np.median(eg[usable]))])
    return rows,re

def main():
    rows=[];summary=[]
    for field,(epochs,bands) in obs.items():
        for b in bands:
            r,s=measure(field,epochs,b);rows+=r;summary.append(s)
    (HERE/"catalog_screen.json").write_text(json.dumps(summary,indent=2)+"\n")
    fields="field band label0 label1 ra dec sep_arcsec epoch0_njy epoch1_gain_corrected_njy formal_err0_njy formal_err1_corrected_njy sigma".split()
    with (HERE/"catalog_flags.csv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=fields);w.writeheader();w.writerows(rows)
    print([(s['field'],s['band'],s['both_epoch_compact_5sigma'],s['flags']) for s in summary])
if __name__=="__main__":main()
