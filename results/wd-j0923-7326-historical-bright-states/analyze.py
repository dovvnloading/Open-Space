"""Rebuild the WD J0923+7326 measurements from public calibrated catalog rows."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd

PIN = "83ecde5b286d824f644cea04f18c5015663d483c"
TARGET = 196121409642518764
RA, DEC = 140.96433333333331, 73.43991111111112
PSCOLS = "objID,detectID,imageID,ra,dec,obsTime,filterID,psfFlux,psfFluxErr,apFlux,apFluxErr,psfQfPerfect,psfMajorFWHM,psfChiSq,infoFlag,infoFlag2,infoFlag3"
URLS = {
    "ps1field.csv": "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/detection.csv?" + urlencode(dict(ra=RA, dec=DEC, radius=.04, pagesize=25000, columns="[" + PSCOLS + "]")),
    "irsa_nei.csv": "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?" + urlencode(dict(POS=f"CIRCLE {RA} {DEC} 0.008333", FORMAT="CSV")),
    "forced.json": f"https://raw.githubusercontent.com/janvanroestel/WD-transitingdebris-alerts/{PIN}/data/forcedphot/WDJ0923%2B7326.json",
}
NEIGHBORS = [196131409579720621, 196121409582756231]
# Fixed field-specific OIDs keep reductions separate and exclude the other two resolved stars.
ZTF_IDS = [[859101300004227, 859201300000131], [859101300004224, 859201300000119], [859101300000058, 859201300000144]]


def ps1(raw):
    x = pd.read_csv(raw)
    good = (x.psfFlux / x.psfFluxErr > 10) & (x.psfQfPerfect >= .90) & (x.psfChiSq < 3) & (x.psfMajorFWHM < 3) & (x.apFlux > 0) & ((x.infoFlag & 8) == 0)
    x = x[good].copy()
    x["mjd"] = x.obsTime
    x["night"] = np.floor(x.mjd).astype(int)
    x["mag"] = -2.5 * np.log10(x.psfFlux / 3631)
    x["apmag"] = -2.5 * np.log10(x.apFlux / 3631)
    group = x.groupby(["objID", "filterID"])
    x["baseline"] = group.mag.transform("median")
    x["n_nights"] = group.night.transform("nunique")
    x["residual"] = x.mag - x.baseline
    separation = 3600 * np.hypot((x.ra - RA) * np.cos(np.deg2rad(DEC)), x.dec - DEC)
    refs = x[(separation > 4) & (x.n_nights >= 3) & x.baseline.between(15.5, 20.5)]
    zp = refs.groupby("imageID").residual.agg(offset="median", n_comp="count")
    x = x.join(zp, on="imageID")
    x = x[x.n_comp >= 5].copy()
    x["dmag"] = x.mag - x.offset
    x["apmag"] -= x.offset
    x["band"] = x.filterID.map(dict(enumerate("grizy", start=1)))
    t = x[x.objID.isin([TARGET] + NEIGHBORS)].copy()
    t["source"] = t.objID.map(dict(zip([TARGET] + NEIGHBORS, ["target", "neighbor_8.6arcsec", "neighbor_9.6arcsec"])))
    t["offset_arcsec"] = 3600 * np.hypot((t.ra - RA) * np.cos(np.deg2rad(DEC)), t.dec - DEC)
    columns = ["source", "objID", "imageID", "detectID", "mjd", "band", "mag", "dmag", "apmag", "offset", "n_comp", "psfQfPerfect", "psfMajorFWHM", "offset_arcsec"]
    return t[columns].sort_values(["mjd","source"]), x


def epoch_contrasts(x):
    # Dates selected from the light curve. Relative ranks are descriptive, not discovery p-values.
    rows = []
    for band in "griz":
        for series in ["psf", "aperture"]:
            col = "dmag" if series == "psf" else "apmag"
            s = x[x.band == band]
            early = (55260, 55280) if band in "gr" else (56280, 56380)
            def window(bounds):
                q = s[s.night.between(*bounds)]
                return q.groupby(["objID", "night"])[col].median().groupby("objID").agg(["median", "size"])
            bright, faint = window(early), window((56640, 56745))
            tab = bright.join(faint, how="inner", lsuffix="_bright", rsuffix="_faint")
            tab["delta"] = tab.median_faint - tab.median_bright
            if TARGET not in tab.index:
                continue
            a = tab.loc[TARGET]
            # Relative comparator brightness and repeated-night criterion identical for peers.
            peers = tab[(tab.median_faint.between(a.median_faint - 1.5, a.median_faint + 1.5)) & (tab.size_bright >= a.size_bright) & (tab.size_faint >= a.size_faint)].drop(index=TARGET)
            rows.append(dict(band=band, estimator=series, bright_mjd=list(early), faint_mjd=[56640,56745], bright_mag=a.median_bright, faint_mag=a.median_faint, delta_mag=a.delta, bright_nights=int(a.size_bright), faint_nights=int(a.size_faint), bright_to_faint_flux_ratio=10**(.4*a.delta), peer_sources=len(peers), peers_as_large=int((peers.delta.abs() >= abs(a.delta)).sum()), neighbors={str(n):float(tab.loc[n,"delta"]) for n in NEIGHBORS if n in tab.index}))
    return rows


def ztf(path):
    d = json.loads(path.read_text())
    x = pd.DataFrame({key:d[key] for key in ("jd","filter","mag","mag_unc")}).drop_duplicates(["jd","filter"])
    x = x[x.mag_unc.between(0, .20, inclusive="neither") & x["filter"].isin(["ZTF_g","ZTF_r"])].copy()
    x["mjd"] = x.jd - 2400000.5
    x["night"] = np.floor(x.mjd)
    return x.groupby(["night","filter"]).agg(mjd=("mjd","median"),mag=("mag","median"),n_exp=("mag","size")).reset_index()


def bins(nights, band, width=14, origin=0):
    sample = nights[nights["filter"]==band]
    group = sample.groupby(np.floor((sample.mjd-57999.5-origin)/width))
    return group.agg(t=("mjd","median"), m=("mag","median"), n=("night","size")).query("n >= 2")


def phase_check(g, epoch, period, width):
    results, draws = [], []
    rng = np.random.default_rng(5407)
    for k in range(2,6):
        e = epoch + period*k
        z = g[np.abs(g.t-e) <= period/2]
        phase, mag = (z.t-e).to_numpy(), z.m.to_numpy()
        def score(offset):
            inside = np.abs((phase-offset+period/2) % period-period/2) < width
            if inside.sum() < 2 or (~inside).sum() < 2:
                return np.nan
            return float(np.median(mag[~inside])-np.median(mag[inside]))
        val=score(0)
        rotating = np.array([score(i) for i in np.arange(0, period, 1)])
        rotating = rotating[np.isfinite(rotating)]
        draws.append(rng.choice(rotating, size=100000))
        results.append(dict(cycle=k, forecast_mjd=e, inside_bins=int((np.abs(phase)<width).sum()), available_bins=len(z), brightening_mag=val))
    scores = [r["brightening_mag"] for r in results]
    tail = float((1+(np.mean(draws,axis=0)>=np.mean(scores)).sum())/100001) if np.all(np.isfinite(scores)) else None
    for r in results:
        if not np.isfinite(r["brightening_mag"]): r["brightening_mag"]=None
    return dict(half_width_days=width, cycles=results, conditional_phase_tail=tail)


def irsa(path):
    x = pd.read_csv(path)
    x = x.query("catflags == 0 and magerr < 0.1")
    rows = []
    for source, ids in zip(["target","neighbor_8.6arcsec","neighbor_9.6arcsec"],ZTF_IDS):
        s=x[x.oid.isin(ids)].copy()
        s["source"]=source
        s["night"]=np.floor(s.mjd).astype(int)
        rows.append(s[["source","oid","filtercode","expid","mjd","mag","magerr","night"]])
    return pd.concat(rows,ignore_index=True)



def anchor_sensitivity(nights):
    rows = []
    for width in [7, 10, 14, 21, 28]:
        values = []
        for origin in range(width):
            g = bins(nights, "ZTF_g", width, origin)
            a = g[g.t < 58349.5].nsmallest(1, "m").iloc[0].t
            b = g[g.t.between(58679.5, 59029.5)].nsmallest(1, "m").iloc[0].t
            p = b - a
            def difference(t):
                return abs((t - a + p / 2) % p - p / 2)
            values.append([p, difference(55268.3), difference(56289.6)])
        arr = np.array(values)
        rows.append({"bin_width": width, "origins_tested": width,
                     "period_range_days": [arr[:,0].min(), arr[:,0].max()],
                     "2010_offset_range_days": [arr[:,1].min(), arr[:,1].max()],
                     "2012_offset_range_days": [arr[:,2].min(), arr[:,2].max()]})
    return rows


def direct_crosscheck(i, epoch, period):
    rows = []
    for (source, band), z in i[i.n_nights >= 2].groupby(["source", "filtercode"]):
        for k in range(2,6):
            center = epoch + k*period
            a = z[np.abs(z.mjd - center) <= period/2]
            mask = np.abs(a.mjd - center) < 56
            value = None
            if mask.sum() >= 2 and (~mask).sum() >= 2:
                value = a.loc[~mask, "mag"].median() - a.loc[mask, "mag"].median()
            rows.append(dict(source=source, band=band, cycle=k, inside_bins=int(mask.sum()),
                             outside_bins=int((~mask).sum()), brightening_mag=value))
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    args.out.mkdir(parents=True,exist_ok=True)
    receipt={}
    for file,url in URLS.items():
        path=args.cache/file
        if not path.exists():
            with urlopen(url,timeout=150) as resp:
                data=resp.read()
            path.write_bytes(data)
        receipt[file]={"url":url,"sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"bytes":path.stat().st_size}
    t,x=ps1(args.cache/"ps1field.csv")
    t.to_csv(args.out/"ps1-measurements.csv",index=False,float_format="%.7f")
    n=ztf(args.cache/"forced.json"); g=bins(n,"ZTF_g")
    a=g[g.t<58349.5].nsmallest(1,"m").iloc[0]
    b=g[g.t.between(58679.5,59029.5)].nsmallest(1,"m").iloc[0]
    epoch,period=a.t,b.t-a.t
    catalogs=irsa(args.cache/"irsa_nei.csv")
    catalogs["block"]=np.floor((catalogs.mjd-57999.5)/14).astype(int)
    i=catalogs.groupby(["source","filtercode","block"]).agg(mjd=("mjd","median"),mag=("mag","median"),n_nights=("night","nunique")).reset_index()
    i.to_csv(args.out/"ztf-direct-blocks.csv",index=False,float_format="%.6f")
    both=pd.concat([bins(n,band).assign(band=band) for band in ["ZTF_g","ZTF_r"]])
    both.to_csv(args.out/"ztf-forced-blocks.csv",index=False,float_format="%.6f")
    out=dict(receipts=receipt,ps1_rows=len(pd.read_csv(args.cache/"ps1field.csv")),ps1_measurements=len(t[t.source=="target"]),ps1_contrasts=epoch_contrasts(x),anchor_mjd=[a.t,b.t],period_days=period,anchor_sensitivity=anchor_sensitivity(n),direct_crosscheck=direct_crosscheck(i,epoch,period),next_window_center_mjd=epoch+period*6,phase_checks={band:[phase_check(bins(n,band),epoch,period,w) for w in [45,56,70,84]] for band in ["ZTF_g","ZTF_r"]})
    (args.out/"measurement.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
