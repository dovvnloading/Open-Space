"""Rebuild the bulge pulsation measurements from OGLE, VizieR and ESA Gaia.

Usage: python replay.py OUTPUT_DIRECTORY [--screen]
Requires numpy, pandas and scipy>=1.15. Successful responses are cached locally.
"""
import argparse
import csv
import hashlib
import io
import json
from math import log
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.signal import lombscargle

OGLE = "https://www.astrouw.edu.pl/ogle/"
III = OGLE + "ogle3/OIII-CVS/blg/lpv/"
IV = OGLE + "ogle4/OCVS/blg/lpv/"
STARS = {
    "046289": (267.16262500000005, -22.26075, "cg1183l11257", 4068994807280968448),
    "020182": (265.660625, -24.239916666666666, "cg1147k8195", 4068323353622502912),
    "151705": (269.645875, -32.025416666666665, "cg0800n26843", None),
    "225280": (273.00895833333334, -28.21172222222222, None, 4050863349471479936),
    "178146": (270.31283333333334, -26.825444444444443, None, 4063972513037118336),
}
CONTROLS = {"225280": ["225131","225753","228233","228672"],
            "178146": ["178788","175061","173732","175508","174372"]}
CATALOGS = {
    III + "Miras.dat": "b117629491b092ca3e7148cca16013ed2709a26b8cabc5bb8c7ce9b6b2766340",
    III + "SRVs.dat": "2921522f6f89abe53e9293f67a2b6fcc2f0696d0b85c48e6dae7c33cfe5d9ee7",
    IV + "Miras.dat": "48acc7f5c2618a0ea24c94fea95a72e7261ea375f37c740613870626b5a26493",
}


class Archive:
    def __init__(self, root):
        self.root = root
        self.raw = root / "raw"
        self.raw.mkdir(parents=True, exist_ok=True)
        self.sources = {}

    def get(self, name, url):
        path = self.raw / name
        if not path.exists():
            for attempt in range(3):
                try:
                    request_url = url.replace("www.astrouw.edu.pl", "ogle.astrouw.edu.pl") if attempt == 1 and url.startswith(OGLE) else url
                    with urlopen(request_url, timeout=90) as response:
                        data = response.read()
                    break
                except HTTPError as exc:
                    if exc.code == 404 or attempt == 2:
                        raise
                except (TimeoutError, URLError):
                    if attempt == 2:
                        raise
                time.sleep(2 * (attempt + 1))
            if not data or data.lstrip().lower().startswith((b"<html",b"<!doctype",b"<?xml")):
                raise ValueError(f"Unexpected data format at {url}")
            path.write_bytes(data)
        raw = path.read_bytes()
        self.sources[name] = {"url":url,"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest()}
        return raw

    def ogle(self, k, era, band="I", days=1):
        name = f"OGLE-BLG-LPV-{k}.dat"
        raw = self.get(f"{k}-{era}-{band}.dat", IV + f"phot_ogle{era}/{band}/" + name)
        t,y,e = np.loadtxt(io.BytesIO(raw), ndmin=2).T
        good = (e>0) & (e<.35)
        return bin_time(t[good], y[good], days)

    def eros(self,k,band="R",neighbor=None):
        ra,dec,source,_ = STARS[k]
        query = urlencode({"-source":"II/390/xtime","-c":f"{ra} {dec}","-c.rs":"5","-out.max":"20000","-out":"**"})
        raw = self.get(k+"-eros.tsv","https://vizier.idia.ac.za/viz-bin/asu-tsv?"+query)
        lines = [line for line in raw.decode().splitlines() if line and not line.startswith("#")]
        if len(lines)<4:
            raise ValueError(f"No EROS measurements received for {k}")
        a = pd.read_csv(io.StringIO("\n".join(lines[:1]+lines[3:])),sep="\t",skipinitialspace=True)
        err, y = a["e_"+band+"mag"], a[band+"mag"]
        use = (a.EROS.str.strip()==(neighbor or source)) & (err>0) & (err<.18) & (y>8) & (y<23)
        return bin_time(a.HJD[use], y[use])

    def gaia(self,k):
        ra,dec,_,source = STARS[k]
        if source is None:
            query = f"SELECT source_id,ra,dec FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},0.0015))"
            url = "https://gea.esac.esa.int/tap-server/tap/sync?"+urlencode(dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",QUERY=query))
            rows = csv.DictReader(io.StringIO(self.get(k+"-gaia-field.csv",url).decode()))
            source = min(rows,key=lambda v: ((float(v["ra"])-ra)*np.cos(np.deg2rad(dec)))**2+(float(v["dec"])-dec)**2)["source_id"]
        url = f"https://gea.esac.esa.int/data-server/data?RETRIEVAL_TYPE=EPOCH_PHOTOMETRY&ID={source}&RELEASE=Gaia+DR3&BAND=G&VALID_DATA=TRUE&FORMAT=ECSV"
        a = pd.read_csv(io.BytesIO(self.get(k+"-gaia.ecsv",url)),comment="#")
        good = ~a.variability_flag_g_reject & ~a.rejected_by_photometry
        return bin_time(a.g_transit_time[good]+5197.5,a.g_transit_mag[good])


def bin_time(t,y,width=1):
    a = pd.DataFrame({"t":t,"y":y}).dropna()
    grouped = a.groupby(np.floor(a.t/width)).median()
    return grouped.t.to_numpy(),grouped.y.to_numpy()


def design(t,periods):
    cols = [np.ones(len(t))]
    for p in periods:
        phase = 2*np.pi*t/p
        cols.extend([np.cos(phase),np.sin(phase)])
    return np.column_stack(cols)


def measure(t,y):
    f = np.linspace(1/1200,1/70,1401)
    t0 = t-t.mean()
    power = lombscargle(t0,y,2*np.pi*f,normalize=True,floating_mean=True)
    i = np.argmax(power)
    def loss(p):
        x = design(t0,[p])
        return np.mean((y-x@np.linalg.lstsq(x,y,rcond=None)[0])**2)
    bounds = (1/f[min(i+1,len(f)-1)],1/f[max(i-1,0)])
    best = minimize_scalar(loss,method="bounded",bounds=bounds)
    return {"period_d":round(best.x,3),"variance_fraction":round(1-best.fun/np.var(y),3),
            "days":len(y),"start":round(t.min(),2),"end":round(t.max(),2)}


def blocks(archive,k,width=1,drop_season=False):
    rows = []
    for era in (3,4):
        t,y = archive.ogle(k,era,days=width)
        edges = np.linspace(t.min()-.01,t.max()+.01,max(1,int(np.ptp(t)//950))+1)
        for left,right in zip(edges[:-1],edges[1:]):
            use = (t>=left)&(t<right)
            tt,yy = t[use],y[use]
            row = measure(tt,yy)
            if drop_season:
                seasons = np.floor((tt-95)/365.25)
                vals = []
                for season in np.unique(seasons):
                    keep = seasons!=season
                    if keep.sum()>=25 and np.ptp(tt[keep])>=700:
                        vals.append(measure(tt[keep],yy[keep])["period_d"])
                row["season_deletion_range_d"] = [min(vals),max(vals)]
            rows.append(row)
    return rows


def harmonic(archive):
    sets = [archive.eros("046289"),archive.ogle("046289",3),archive.ogle("046289",4)]
    grid = np.arange(238.,306.,.25)
    scores = []
    for t,y in sets:
        seasons = np.floor((t-95)/365.25)
        def rms(p):
            x = design(t,[p,p/2])
            residual = []
            for s in np.unique(seasons):
                test = seasons==s
                if test.sum()>=6 and (~test).sum()>=20:
                    b = np.linalg.lstsq(x[~test],y[~test],rcond=None)[0]
                    residual.extend((y[test]-x[test]@b)**2)
            return np.mean(residual)
        scores.append([rms(p) for p in grid])
    scores = np.array(scores)
    norm = np.array([np.var(y) for t,y in sets])
    i = np.argmin((scores/norm[:,None]).sum(axis=0))
    return {"common_period_d":grid[i],"separate_d":grid[np.argmin(scores,axis=1)].tolist(),
            "common_rms_mag":np.sqrt(scores[:,i]).round(4).tolist(),
            "separate_rms_mag":np.sqrt(scores.min(axis=1)).round(4).tolist()}


def screen(archive):
    catalogs = []
    for i,(url,expected) in enumerate(CATALOGS.items()):
        raw = archive.get(f"catalog-{i}.dat",url)
        if hashlib.sha256(raw).hexdigest()!=expected:
            raise ValueError(f"Original catalog changed: {url}")
        catalogs.append({a[0]:float(a[3]) for a in (s.split() for s in raw.decode().splitlines())})
    old,late = catalogs[0]|catalogs[1],catalogs[2]
    shared = old.keys() & late.keys()
    hits = sorted((k for k in shared if abs(log(late[k]/old[k]))>.07),key=lambda k:-abs(log(late[k]/old[k])))
    assert len(shared)==6539 and len(hits)==52
    with (archive.root/"screen.csv").open("w",newline="") as out:
        w = csv.writer(out)
        w.writerow(["id","catalog_III","catalog_IV","days_III","fit_III","days_IV","fit_IV"])
        for name in hits:
            row = [name,old[name],late[name]]
            for era in (3,4):
                try:
                    t,y = archive.ogle(name[-6:],era)
                    row.extend([len(t),measure(t,y)["period_d"] if len(t)>=10 else ""])
                except HTTPError as e:
                    if e.code!=404: raise
                    row.extend([0,"404"])
            w.writerow(row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--screen",action="store_true")
    args = parser.parse_args()
    a = Archive(args.output)
    j = {"dominance":{},"shortening":{},"controls":{}}
    for k,(_,_,eros_id,_) in STARS.items():
        if eros_id:
            raw = {"OGLE-III I":a.ogle(k,3),"OGLE-IV I":a.ogle(k,4),"OGLE-IV V":a.ogle(k,4,"V"),
                   "EROS R":a.eros(k),"EROS B":a.eros(k,"B"),"Gaia G":a.gaia(k)}
            t,y = raw["EROS R"]
            mask = t<raw["OGLE-III I"][0].min()
            raw["EROS before OGLE-III"]=(t[mask],y[mask])
            j["dominance"][k]={s:measure(*v) for s,v in raw.items()}
            ps=[j["dominance"][k][s]["period_d"] for s in ("OGLE-III I","OGLE-IV I")]
            for label in ("OGLE-III I","OGLE-IV I"):
                t,y=raw[label]
                b=np.linalg.lstsq(design(t,ps),y,rcond=None)[0]
                j["dominance"][k][label]["joint_peak_to_peak_mag"]=(2*np.hypot(b[1::2],b[2::2])).round(3).tolist()
        else:
            j["shortening"][k]={"days":blocks(a,k,drop_season=True),"three_days":blocks(a,k,3),
                                "Gaia G":measure(*a.gaia(k)),"OGLE-IV V":measure(*a.ogle(k,4,"V"))}
    for k,others in CONTROLS.items():
        j["controls"][k]={v:blocks(a,v) for v in others}
    j["harmonic_046289"]=harmonic(a)
    j["EROS_neighbor_046289"]=measure(*a.eros("046289",neighbor="cg1183l11405"))
    if args.screen:
        screen(a)
    (args.output/"measurements.json").write_text(json.dumps(j,indent=2)+"\n")
    (args.output/"sources.json").write_text(json.dumps(a.sources,indent=2)+"\n")


if __name__ == "__main__":
    main()
