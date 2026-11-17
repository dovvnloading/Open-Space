"""Measure signed flux changes at the Warhol flag and all blind peak positions."""
import argparse
import csv
from pathlib import Path

import numpy as np

import macs0416 as obs

HERE = Path(__file__).resolve().parent
EXTRAS = (
    ("jw01176-o211_t009", "f444w"),
    ("jw01176-o212_t009", "f444w"),
    ("jw01208-o004_t002", "f444w"),
    ("jw03538-o019_t003", "f300m"),
    ("jw06882-o054_t054", "f300m"),
)


def extract_pair(model, grid, cache, early, late, band, positions):
    paths = [cache / f"{e}_nircam_clear-{band}.fits" for e in (early, late)]
    d = model.analyze(band, paths, grid)
    receipts = d["meta"]["epoch"]
    rows = []
    for peak in positions:
        x, y = int(peak["x"]), int(peak["y"])
        empirical_err = d["noise"][y, x] * d["to_njy"] / d["resp"]
        rows.append(dict(
            flag_band=peak["band"], x=x, y=y, ra=peak["ra"], dec=peak["dec"],
            compared_band=band, earlier_product=paths[0].stem, later_product=paths[1].stem,
            earlier_mjd=receipts[0]["mjd"], later_mjd=receipts[1]["mjd"],
            valid=bool(d["good"][y, x]), delta_njy=float(d["score"][y, x]*empirical_err),
            local_sigma_njy=float(empirical_err), signed_score=float(d["score"][y, x]),
            astrometry_rms_arcsec=d["meta"]["astrom_rms_arcsec"],
        ))
    return rows


def main():
    import json
    import sys
    from concurrent.futures import ThreadPoolExecutor
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("action",choices=["fetch","measure"])
    p.add_argument("cache",type=Path)
    p.add_argument("--base",type=Path,default=obs.LEGACY)
    args=p.parse_args()
    if args.action=="fetch":
        sys.argv=[sys.argv[0],str(args.cache)]
        m=obs.module("fetch_elgordo",args.base)
        m.CENTER=obs.CENTER
        m.SIZE_ARCSEC=39
        with ThreadPoolExecutor(2) as pool:
            receipts=list(pool.map(m.fetch,EXTRAS))
        (HERE/"extra_ranges.json").write_text(json.dumps(receipts,indent=1)+"\n")
        return
    m=obs.module("analyze_elgordo",args.base)
    m.CENTER=obs.CENTER
    m.SIDE=601
    grid=m.getgrid()
    with (HERE/"single_band_peaks.csv").open() as h:
        peaks=list(csv.DictReader(h))
    rows=[]
    for band in obs.BANDS:
        rows.extend(extract_pair(m,grid,args.cache,*obs.EPOCHS,band,peaks))
    warhol=[p for p in peaks if int(p["x"])==345 and int(p["y"])==378]
    assert len(warhol)==1
    for ep,band in EXTRAS[:-1]:
        rows.extend(extract_pair(m,grid,args.cache,ep,obs.EPOCHS[1],band,warhol))
    with (HERE/"forced_differences.csv").open("w",newline="") as h:
        writer=csv.DictWriter(h,fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)

if __name__=="__main__":main()
