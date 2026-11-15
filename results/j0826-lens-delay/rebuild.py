"""Rebuild the reduced tables from archive responses. Requires NumPy and pandas."""
import argparse
import csv
import hashlib
import json
import urllib.request
from pathlib import Path
from xml.etree import ElementTree
import numpy as np
import pandas as pd
from replay import ROOT, TARGET, read


def write(folder, name, headers, rows):
    with (folder / name).open("w") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(rows)


def distance(ra, dec, best_ra, best_dec):
    return 3600 * np.hypot((ra - best_ra) * np.cos(np.deg2rad(best_dec)), dec - best_dec)


def ztf_quality(df):
    good = (df.catflags == 0) & (df.magerr > 0) & (df.magerr < .28)
    return df[good & (df.mag < df.limitmag - .05)]


def make_ztf(cache, out, targets):
    nights = []
    for name, group in targets.groupby("system", sort=False):
        src = cache / "lcs" / (name.replace(" ", "_").replace("/", "_") + ".csv")
        df = pd.read_csv(src)
        if not len(df):
            continue
        df = ztf_quality(df)
        objects = df.groupby("oid")[["ra", "dec"]].median()
        for row in group.itertuples():
            ids = objects.index[distance(objects.ra, objects.dec, row.ra_deg, row.dec_deg) < .85]
            for band, selected in df[df.oid.isin(ids)].groupby("filtercode"):
                if len(selected) < 8 or (name != TARGET and band == "zi"):
                    continue
                grouped = selected.groupby(selected.mjd.astype(int))[["mjd", "mag", "magerr"]].median()
                for t, m, e in grouped.to_numpy():
                    nights.append([name, row.image, band, round(t, 7), round(m, 6), round(e, 6)])
    write(out, "ztf-nights.csv", ["system","image","band","mjd","magnitude","median_catalog_error"], nights)


def make_gaia(cache, out):
    rows = []
    for image, file in [("A", "aip0.csv"), ("B", "aipB.xml")]:
        doc = ElementTree.parse(cache / file)
        names = [f.attrib["name"] for f in doc.findall(".//{*}FIELD")]
        values = dict(zip(names, [e.text.split() for e in doc.findall(".//{*}TD")]))
        for band, prefix in [("g", "g_transit"), ("bp", "bp"), ("rp", "rp")]:
            t = "g_transit_time" if band == "g" else band + "_obs_time"
            col = lambda k: pd.to_numeric(pd.Series(values[k]), errors="coerce").to_numpy()
            obs = pd.DataFrame(dict(time=55197 + col(t), mag=col(prefix + "_mag"),
                                    err=1.086 / col(prefix + "_flux_over_error")))
            ok = np.array(values["variability_flag_" + band + "_reject"]) == "False"
            obs = obs[ok & np.isfinite(obs.sum(1,skipna=False)) & (obs.err > .00001) & (obs.err < .3)]
            for t, m, e in obs.groupby(obs.time.astype(int)).median().to_numpy():
                rows.append([image, band, round(t,7), round(m,7), round(e,7)])
    write(out,"gaia-nights.csv",["image","band","bjd_tcb_minus_2400000p5","magnitude","median_catalog_error"],rows)


def make_followup(cache, out, positions):
    ps = pd.read_csv(cache / "ps.csv")
    rows = []
    for row in positions.itertuples():
        sub = ps[distance(ps.ra, ps.dec, row.ra_deg, row.dec_deg) < .8]
        sub = sub[(sub.psfQfPerfect > .8) & (sub.psfFlux > 5*sub.psfFluxErr) & (sub.psfFluxErr > 0)]
        for p in sub.itertuples():
            rows.append([row.image,p.objID,p.detectID,p.filterID,p.obsTime,p.psfFlux,p.psfFluxErr])
    write(out,"ps1-detections.csv",["image","object_id","detection_id","filter_id","mjd","psf_jy","error_jy"],rows)

    a = positions[positions.image == "A"].iloc[0]
    alerts = []
    for file in ["ZTF18achddwe.json","ZTF18adpabfq.json"]:
        alerts.extend(json.loads((cache / file).read_text()))
    df = pd.DataFrame(alerts).drop_duplicates("candid")
    df["off"] = distance(df.ra,df.dec,a.ra_deg,a.dec_deg)
    selected = df[(df.rb > .55) & (~df.dubious) & (df.sigmapsf < .18) & (df.off < 1.2)]
    cols = ["candid","fid","mjd","ra","dec","off","isdiffpos","magpsf","rb"]
    selected[cols].to_csv(out / "difference-detections.csv",index=False,float_format="%.7f")

    field = ztf_quality(pd.read_csv(cache / "field.csv"))
    rows = []
    for (band, oid), group in field[field.field==840].groupby(["filtercode","oid"]):
        ra, dec = group.ra.median(), group.dec.median()
        label = next((p.image for p in positions.itertuples()
                      if distance(ra,dec,p.ra_deg,p.dec_deg)<.85),"field")
        stats = []
        for lo, hi in [(58900,59500),(60350,60970)]:
            sub = group[(group.mjd>=lo)&(group.mjd<hi)]
            med = sub.groupby(sub.mjd.astype(int)).mag.median()
            stats.extend([len(med),round(med.median(),6)])
        if min(stats[::2])>10:
            rows.append([band,oid,label,round(ra,6),round(dec,6),*stats,round(stats[3]-stats[1],6)])
    write(out,"field-controls.csv",["band","oid","object","ra","dec","early_nights","early_median",
                                   "late_nights","late_median","late_minus_early_mag"],rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    if args.out.resolve() == ROOT:
        parser.error("choose a different output directory so the frozen products remain available for comparison")
    args.out.mkdir(parents=True,exist_ok=True)
    for source in read("sources.csv"):
        path = args.cache / source["source_filename"]
        if not path.exists():
            path.parent.mkdir(parents=True,exist_ok=True)
            with urllib.request.urlopen(source["url"], timeout=180) as response:
                path.write_bytes(response.read())
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != source["sha256"]:
            print(f"{path.name}: archive response differs from frozen retrieval ({digest})")
    targets = pd.read_csv(ROOT / "targets.csv")
    make_ztf(args.cache,args.out,targets)
    make_gaia(args.cache,args.out)
    make_followup(args.cache,args.out,targets[targets.system==TARGET])


if __name__ == "__main__":
    main()
