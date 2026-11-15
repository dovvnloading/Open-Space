"""Measure late optical rebrightening and centroid shifts in the SN 2023alo field."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / "inputs"
OBJECTS = {"zg": 375103300019510, "zr": 375203300006172}
WINDOWS = {"pre": (58100, 59850), "first": (59900, 60190),
           "trough": (60420, 60650), "late": (60725, 60855)}


def interval(table, bounds):
    return table.loc[(table.index >= bounds[0]) & (table.index <= bounds[1])]


def measure_catalog():
    data = pd.read_csv(INPUTS / "ztf_field.csv.gz")
    # Same primary field, detector, quadrant and filter as the target.
    data = data[(data.catflags == 0) & (data.magerr < .15)]
    results = {}
    for band, target in OBJECTS.items():
        source = data[data.oid == target].iloc[0]
        selected = data[(data.filtercode == band) & (data.field == source.field) &
                        (data.ccdid == source.ccdid) & (data.qid == source.qid)].copy()
        selected["night"] = np.floor(selected.mjd).astype(int)
        nightly = selected.pivot_table(index="night", columns="oid", values="mag", aggfunc="median")
        nightly = nightly[nightly[target].notna()]
        epochs = {name: interval(nightly,bounds) for name,bounds in WINDOWS.items()}
        medians = pd.DataFrame({name:x.median().where(x.count() >= max(5, .55*len(x)))
                                for name,x in epochs.items()}).T
        ids = medians.columns[(medians.notna().all()) & (medians.columns != target) &
                              (abs(medians.loc["pre"] - medians.loc["pre", target]) < 2.5)]
        changes = {str(oid):float(medians.at["late",oid]-medians.at["trough",oid]) for oid in ids}
        # Seven-day blocks sample the observed cadence. These errors describe the mean change
        # across the selected seasons and do not include a search-trials or host-model term.
        rng = np.random.default_rng(1984)
        estimated = []
        for key in ["trough","late"]:
            s = epochs[key][target].dropna()
            blocks = s.groupby(s.index // 7).median().values
            sampled = rng.choice(blocks, (10000,len(blocks)), replace=True)
            estimated.append(np.median(sampled,axis=1))
        results[band] = {"catalog_object":target,
                         "nights":{k:len(v) for k,v in epochs.items()},
                         "median_magnitude":{k:float(medians.at[k,target]) for k in WINDOWS},
                         "difference_late_minus_trough":float(medians.at["late",target]-medians.at["trough",target]),
                         "weekly_block_difference_percentiles":np.percentile(estimated[1]-estimated[0],[2.5,50,97.5]).tolist(),
                         "controls_late_minus_trough":changes}
    return results


def measure_alerts():
    alerts = pd.read_json(INPUTS / "ztf_alerts.json.gz")
    alerts = alerts[(alerts.isdiffpos == 1) & (alerts.rb > .55)].copy()
    host = pd.read_csv(INPUTS / "panstarrs.csv", dtype={"objID":str}).sort_values("nDetections").iloc[-1]
    host_center = np.array([host.raMean, host.decMean])
    scale = 3600 * np.array([np.cos(np.radians(host.decMean)),1.0])
    astrometry = {}
    for select in ["all",1,2]:
        data = alerts if select == "all" else alerts[alerts.fid == select]
        rng = np.random.default_rng(2026)
        bootstrap = []
        positions = {}
        for key in ["first","late"]:
            low,high = WINDOWS[key]
            g = data[data.mjd.between(low,high)]
            points = g.groupby(np.floor(g.mjd))[["ra","dec"]].median().values
            draws = rng.integers(0,len(points),(5000,len(points)))
            bootstrap.append(np.median(points[draws],axis=1))
            positions[key] = {"nights":len(points),"alerts":len(g),
                              "host_offset_east_north_arcsec":((np.median(points,axis=0)-host_center)*scale).tolist()}
        shift = (bootstrap[1]-bootstrap[0])*scale
        astrometry[str(select)] = {"positions":positions,
                                  "shift_east_north_percentiles_arcsec":np.percentile(shift,[2.5,50,97.5],axis=0).tolist(),
                                  "shift_length_percentiles_arcsec":np.percentile(np.linalg.norm(shift,axis=1),[2.5,50,97.5]).tolist()}
    light = {}
    for fid,g in alerts.groupby("fid"):
        q=g[g.mjd.between(*WINDOWS["late"])]
        light[str(fid)] = {"template_ids":sorted(int(k) for k in g.rfid.dropna().unique()),
                           "late_nights":len(np.unique(np.floor(q.mjd))),
                           "brightest_late_diffmag":float(q.magpsf.min())}
    return {"host_Gaia_DR3":6300335462633830528,"host_PS1_objID":host.objID,
            "photometry":light,"astrometry":astrometry}


def asassn():
    data = pd.read_parquet(INPUTS / "asassn.parquet")
    data = data[(data.quality == "G") & (data.flux_err < .15)].copy()
    data["mjd"] = data.jd - 2400000.5
    samples = data.set_index("mjd")
    return {"usable_rows_by_window":{k:len(interval(samples,v)) for k,v in WINDOWS.items()},
            "late_visits":interval(samples,WINDOWS["late"]).reset_index()[["mjd","flux","flux_err","camera"]].to_dict("records")}


if __name__ == "__main__":
    for path,record in json.loads((ROOT / "provenance.json").read_text()).items():
        if hashlib.sha256((INPUTS/path).read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError(f"Input checksum changed for {path}")
    print(json.dumps({"windows_mjd":WINDOWS,"ztf_catalog":measure_catalog(),
                      "ztf_difference":measure_alerts(),"asassn":asassn()},indent=2,allow_nan=False))
