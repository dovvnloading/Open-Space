"""Inspect the identified tracks and compare their radiants to the IAU MDC."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from survey import COLS, unit, xyz

HERE = Path(__file__).resolve().parent
MDC_URL = "https://iaumdc.astro.sk/downloads/lists_shw_data/streamfulldata.json"
DETAIL = {**COLS, 5: "solar_lon", 8: "sigma_ra", 10: "sigma_dec", 16: "sigma_vg",
          23: "a_au", 25: "ecc", 27: "incl_deg", 37: "q_au", 59: "v_init",
          64: "sigma_lat", 66: "sigma_lon", 68: "sigma_h", 75: "duration_s",
          76: "abs_mag", 81: "residual_arcsec", 85: "stations"}


def collect(cache):
    groups = pd.read_csv(HERE/"primary/groups.csv")
    ids = {x for x in groups.iloc[:, 1:].to_numpy().ravel() if isinstance(x, str)}
    frames = []
    for month in sorted({key[:6] for key in ids}):
        df = pd.read_csv(cache/(month+".txt"), sep=";", header=None, comment="#",
                         skipinitialspace=True, usecols=sorted(DETAIL)).rename(columns=DETAIL)
        df = df.loc[df.id.str.strip().isin(ids)].copy()
        for key in ["id", "utc", "shower", "stations"]:
            df[key] = df[key].str.strip()
        frames.append(df)
    pd.concat(frames).sort_values("utc").to_csv(HERE/"tracks.csv", index=False)


def measure(frame):
    pos, direction = xyz(frame)
    rad = unit(frame.ra, frame.dec)
    start = unit(frame.lon, frame.lat)
    time = pd.to_datetime(frame.utc)
    rows = []
    for i, j in itertools.combinations(range(len(frame)), 2):
        off = pos[i]-pos[j]
        chord = np.linalg.norm(start[i]-start[j])
        rows.append({
            "id_i": frame.id.iat[i], "id_j": frame.id.iat[j],
            "dt_s": (time.iat[j]-time.iat[i]).total_seconds(),
            "radiant_deg": np.rad2deg(2*np.arcsin(np.linalg.norm(rad[i]-rad[j])/2)),
            "speed_difference_kms": abs(frame.vg.iat[i]-frame.vg.iat[j]),
            "entry_ground_arc_km": 6371*2*np.arcsin(chord/2),
            "cross_track_km": min(np.linalg.norm(np.cross(off,direction[k])) for k in (i,j)),
        })
    return pd.DataFrame(rows)


def catalog_matches(shower_catalog, representatives):
    rows = []
    for date, meteors in representatives.groupby(representatives.utc.str[:10]):
        ra, dec, v, sol = meteors[["ra","dec","vg","solar_lon"]].median()
        rad = unit([ra], [dec])[0]
        candidates = []
        for shower in shower_catalog["data"]:
            for s in shower["solution"]:
                if any(s.get(k) is None for k in ["LoS", "Ra", "De", "Vg"]):
                    continue
                delta = (sol - s["LoS"] + 180) % 360 - 180
                if abs(delta) > 15:
                    continue
                a = s["Ra"] + (s.get("dRa") or 0)*delta
                b = s["De"] + (s.get("dDe") or 0)*delta
                if abs(b) > 90:
                    continue
                angle = np.rad2deg(2*np.arcsin(np.linalg.norm(unit([a],[b])[0]-rad)/2))
                dv = s["Vg"]-v
                candidates.append(dict(date=date, code=shower["Code"], name=shower["Name"],
                                       solution=s["AdNo"], flag=s["s"], d_solar_lon_deg=delta,
                                       radiant_offset_deg=angle, speed_offset_kms=dv,
                                       scaled_distance=np.hypot(angle/4, dv/2.2)))
        rows.extend(sorted(candidates, key=lambda x: x["scaled_distance"])[:8])
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", type=Path, help="extract track rows from monthly input directory")
    parser.add_argument("--mdc", type=Path, help="compare all IAU MDC solutions in this JSON release")
    args = parser.parse_args()
    if args.extract:
        collect(args.extract)
    d = pd.read_csv(HERE/"tracks.csv")
    groups = pd.read_csv(HERE/"geographic/groups.csv")
    ids = {x for x in groups.iloc[:, 1:].to_numpy().ravel() if isinstance(x, str)}
    representatives = d[d.id.isin(ids)]
    for date, subset in representatives.groupby(representatives.utc.str[:10]):
        print(date, len(subset), "paths")
        print(measure(subset).round(3).to_string(index=False))
    august = representatives[representatives.utc.str.startswith("2026-08-24")]
    assert len(august) == 4
    assert len(set.union(*(set(s.split(",")) for s in august.stations))) == august.nst.sum()
    pairs = measure(august)
    assert 8.94 < pairs.dt_s.max() < 8.95 and pairs.cross_track_km.min() > 200
    if args.mdc:
        with args.mdc.open("rb") as f:
            print("MDC SHA-256", hashlib.file_digest(f, "sha256").hexdigest())
        results = catalog_matches(json.loads(args.mdc.read_bytes()), representatives)
        results.to_csv(HERE/"shower-comparison.csv", index=False, float_format="%.5f")


if __name__ == "__main__":
    main()
