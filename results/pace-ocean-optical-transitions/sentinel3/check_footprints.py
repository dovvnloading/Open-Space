"""Collocate independently gridded OLCI Rrs with each pair of matched PACE coordinates."""

import argparse
import csv
import datetime as dt
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from olci import WAVES, curvature, distance, save
from olci_zarr import RECEIPTS, Grid


def compare(path):
    with open(path) as f:
        records = list(csv.DictReader(f))
    grid = Grid("timeChunked")
    lat = np.array([float(r["lat"]) for r in records])
    lon = np.array([float(r["lon"]) for r in records])
    y, x = grid.bounds(
        lat.min() - 0.1, lat.max() + 0.1, lon.min() - 0.1, lon.max() + 0.1
    )
    ys, xs = grid.lat[y], grid.lon[x]
    dates = sorted({r[k] for r in records for k in ["start", "end"]})
    jobs = [(d, b) for d in dates for b in WAVES]

    def load(job):
        day, wave = job
        return grid.extract(wave, dt.date.fromisoformat(day), y, x)

    with ThreadPoolExecutor(max_workers=12) as pool:
        bands = dict(zip(jobs, pool.map(load, jobs)))
    output = []
    grouped = defaultdict(list)
    for r, la, lo in zip(records, lat, lon):
        i = np.abs(ys - la).argmin()
        j = np.abs(xs - lo).argmin()
        oci, olci, red_ok, raw = [], [], [], []
        for k in ["start", "end"]:
            blue, cyan, green = [float(r[f"{k}_Rrs_{b}"]) for b in (442, 490, 560)]
            oci.append(np.log(cyan) - (70 * np.log(blue) + 48 * np.log(green)) / 118)
            native = [bands[r[k], b][i, j] for b in WAVES]
            raw.extend(native)
            blue, cyan, green, red = native
            olci.append(float(curvature(native)))
            red_ok.append(
                all(0 < v < 0.05 for v in native[:3])
                and green > 0.00005
                and 0 <= red < 0.0015
                and red / green < 0.5
            )
        key = (r["region"], r["start"], r["end"])
        output.append(
            [
                *key,
                la,
                lo,
                r["depth_m"],
                ys[i],
                xs[j],
                distance(la, lo, ys[i], xs[j]),
                *oci,
                *olci,
                *raw,
                all(red_ok),
            ]
        )
        grouped[key].append((oci[1] - oci[0], olci[1] - olci[0], all(red_ok)))
    hdr = "region,start,end,oci_lat,oci_lon,depth_m,olci_lat,olci_lon,distance_km,oci_C_before,oci_C_after,olci_C_before,olci_C_after,"
    hdr += (
        ",".join(f"{when}_olci_{w}" for when in ["start", "end"] for w in WAVES)
        + ",olci_screen_pass"
    )
    save("olci_footprints.csv", hdr, output)
    summaries = []
    for key, changes in sorted(grouped.items()):
        p, o, q = np.array(changes).T
        both = np.isfinite(p) & np.isfinite(o)
        for label, valid in [
            ("all_positive_bands", both),
            ("plus_red_and_range_checks", both & q.astype(bool)),
        ]:
            a, b = p[valid], o[valid]
            qs = np.percentile(b, [25, 75]) if len(b) else ["", ""]
            rho = (
                np.corrcoef(a, b)[0, 1]
                if len(a) > 1 and a.std() > 0 and b.std() > 0
                else ""
            )
            summaries.append(
                [
                    *key,
                    label,
                    len(p),
                    len(a),
                    np.median(a) if len(a) else "",
                    np.median(b) if len(a) else "",
                    *qs,
                    np.mean(np.sign(a) == np.sign(b)) if len(a) else "",
                    rho,
                ]
            )
    save(
        "olci_footprints_summary.csv",
        "region,start,end,selection,oci_pairs,collocations,oci_median_delta_C,olci_median_delta_C,olci_q25,olci_q75,paired_sign_agreement,pearson_r",
        summaries,
    )
    save(
        "footprint_receipts.csv",
        "url,bytes,sha256",
        [(u, *v) for u, v in sorted(RECEIPTS.items())],
    )


if __name__ == "__main__":
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument(
        "footprints",
        help="CSV of paired PACE coordinates and native 442/490/560 nm reflectances",
    )
    compare(args.parse_args().footprints)
