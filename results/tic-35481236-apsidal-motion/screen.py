"""Replay the 43-object exploratory QLP screen. Output includes unavailable products."""
import argparse
import csv
import urllib.error
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
from tess_stars2px import tess_stars2px_function_entry

CATALOG = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J%2FApJ%2F912%2F123&-out.all&-out.max=20000"


def select_rows(catalog):
    results = []
    # Two CDS tables share this response. Table 2 has the empirical timings.
    for line in catalog.split("#Table")[1].splitlines():
        v = line.split("\t")
        if len(v)!=12 or not v[0].strip().isdigit():
            continue
        _, tic, period, t1, t2, ecosw, d1, d2, mag, _, ra, dec = v
        if .02<float(ecosw)<.6 and float(d2)>.02 and float(period)<7 and float(mag)<13.5:
            results.append((int(tic),*map(float,(period,t1,t2,ecosw,d1,d2,mag,ra,dec))))
    return results


def center(t,f,ref,other,period):
    x = (t-ref+period/2)%period-period/2
    window = min(.16*period,.40*min((other-ref)%period,(ref-other)%period))
    keep = abs(x)<window
    if keep.sum()<20:
        return
    x, y = x[keep], f[keep]
    level = np.percentile(y,80)
    edges = np.linspace(-window,window,41)
    bins = [y[(x>=a)&(x<b)] for a,b in zip(edges,edges[1:])]
    med = np.array([np.median(values) if len(values) else np.nan for values in bins])
    initial = np.clip(edges[np.nanargmin(med)]+window/40,-.69*window,.69*window)
    depth = max(.001,level-np.nanmin(med))
    def residual(p):
        c,w,k,d,b,s = p
        return b+s*x-d*np.exp(-.5*abs((x-c)/w)**k)-y
    scale = max(.00025,np.median(abs(y-np.median(y)))*.25)
    try:
        return least_squares(residual,[initial,window/5,2,depth,level,0],
                             bounds=([-.85*window,.004,1,0,level-.30,-2],
                                     [.85*window,window,8,1.7,level+.30,2]),
                             loss="soft_l1",f_scale=scale).x
    except ValueError:
        return


def measure(task):
    row, sector, cache = task
    tic, period, t1,t2,*_ = row
    key = f"{tic:016d}"
    folder = "/".join(key[i:i+4] for i in range(0,16,4))
    url = (f"https://archive.stsci.edu/hlsps/qlp/s{sector:04d}/{folder}/"
           f"hlsp_qlp_tess_ffi_s{sector:04d}-{key}_tess_v01_llc.fits")
    dest = cache/f"{tic}_{sector}.fits"
    output = dict(tic=tic,sector=int(sector),period_days=period,requested_url=url,
                  status="no fit",median_btjd=None,primary_offset_days=None,secondary_offset_days=None,
                  differential_minutes=None,depth1=None,depth2=None)
    if not dest.exists():
        try:
            with urlopen(url,timeout=30) as stream:
                data = stream.read()
        except urllib.error.HTTPError as exc:
            output["status"] = f"HTTP {exc.code}"
            return output
        if data[:6]!=b"SIMPLE":
            output["status"] = "not FITS"
            return output
        dest.write_bytes(data)
    tab = fits.getdata(dest)
    t = tab["TIME"]
    f = tab["KSPSAP_FLUX" if "KSPSAP_FLUX" in tab.names else "DET_FLUX"].astype(float)
    use = (tab["QUALITY"]==0)&np.isfinite(t+f)
    if use.sum()<80:
        return output
    t, f = t[use], f[use]/np.median(f[use])
    a,b = center(t,f,t1,t2,period),center(t,f,t2,t1,period)
    if a is not None and b is not None:
        output.update(status="fit",median_btjd=np.median(t),primary_offset_days=a[0],
                      secondary_offset_days=b[0],differential_minutes=(b[0]-a[0])*1440,
                      depth1=a[3],depth2=b[3])
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,default=Path("screen-cache"))
    parser.add_argument("--output",type=Path,default=Path("screen.csv"))
    args = parser.parse_args()
    args.cache.mkdir(exist_ok=True,parents=True)
    catalog = args.cache/"justesen.tsv"
    if not catalog.exists():
        with urlopen(CATALOG,timeout=60) as r:
            catalog.write_bytes(r.read())
    stars = select_rows(catalog.read_text())
    tic, ra, dec = zip(*[(r[0],r[-2],r[-1]) for r in stars])
    projected = tess_stars2px_function_entry(list(tic),list(ra),list(dec))
    tasks = []
    for row in stars:
        sectors = sorted({int(s) for s in projected[3][projected[0]==row[0]] if 1<=s<=104})
        if sectors:
            pick = {sectors[0], sectors[-1], sectors[len(sectors)//2], sectors[max(0,len(sectors)-2)]}
            tasks.extend((row,sector,args.cache) for sector in sorted(pick))
    with ThreadPoolExecutor(8) as pool, args.output.open("w",newline="") as stream:
        for index,result in enumerate(pool.map(measure,tasks)):
            if index==0:
                writer = csv.DictWriter(stream,fieldnames=result.keys(),lineterminator="\n")
                writer.writeheader()
            writer.writerow(result)


if __name__=="__main__":
    main()
