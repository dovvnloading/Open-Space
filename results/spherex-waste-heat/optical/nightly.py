"""Recalculate nightly ZTF r residuals across field/quadrant-specific light curves."""
from collections import defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from statistics import median
import sys
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "wise"))
from fetch import position, separation, targets


def collect(target, rows, catalog_positions=False):
    fields = defaultdict(lambda: defaultdict(list))
    for row in rows:
        mjd = float(row["mjd"])
        if int(row["catflags"]) or not 0 < float(row["magerr"]) <= 0.1:
            continue
        epoch = 2000 + (mjd - 51544.5) / 365.25
        if not catalog_positions and separation(float(row["ra"]), float(row["dec"]),
                                                *position(target, epoch)) > 1.5:
            continue
        fields[row["oid"]][int(mjd)].append(row)
    dates = defaultdict(list)
    for oid, nights in fields.items():
        if len(nights) < 5:
            continue
        means = {day: median([float(row["mag"]) for row in points]) for day, points in nights.items()}
        center = median(means.values())
        for day, points in nights.items():
            dates[day].append((oid, means[day] - center,
                               median([float(row["magerr"]) for row in points]),
                               [row["mjd"] for row in points]))
    nights = []
    for day, groups in sorted(dates.items()):
        nights.append(dict(label=target["label"], mjd_day=day, r_residual=median(x[1] for x in groups),
                           err_proxy=median(x[2] for x in groups), oids="|".join(sorted(x[0] for x in groups)),
                           unique_times=len({time for x in groups for time in x[3]})))
    if nights:
        center = median(row["r_residual"] for row in nights)
        scatter = 1.4826 * median(abs(row["r_residual"] - center) for row in nights)
        for row in nights:
            row["r_residual"] -= center
            row["scatter_proxy"] = scatter
            row["local_score"] = row["r_residual"] / max(scatter, row["err_proxy"])
    return nights


def write(path, keys, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, keys.split(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    receipts = {r["label"]: r for r in json.loads((ROOT / "query_receipts.json").read_text())}
    hats_path = ROOT / "hats_receipts.json"
    hats = {r["label"]: r for r in json.loads(hats_path.read_text())["targets"]} if hats_path.exists() else {}
    all_nights = []
    summary = []
    with zipfile.ZipFile(ROOT / "lightcurve_inputs.zip") as archive:
        for target in targets():
            label = target["label"]
            receipt = receipts[label]
            source = "IRSA API"
            if "sha256" in receipt:
                data = gzip.decompress(archive.read(f"cache/{label}.csv.gz"))
                if hashlib.sha256(data).hexdigest() != receipt["sha256"]:
                    raise ValueError(f"ZTF input hash differs for {label}")
                rows = list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))
                nights = collect(target, rows)
            else:
                fallback = ROOT / f"hats_{label}.csv.gz"
                if fallback.is_file():
                    source = "IRSA HATS"
                    data = gzip.decompress(fallback.read_bytes())
                    if hashlib.sha256(data).hexdigest() != hats[label]["sha256"]:
                        raise ValueError(f"HATS input hash differs for {label}")
                    rows = list(csv.DictReader(io.StringIO(data.decode())))
                    nights = collect(target, rows, catalog_positions=True)
                else:
                    rows, nights, source = [], [], "API query failed"
            rec = dict(label=label, source=source, raw=len(rows), nights=len(nights))
            if nights:
                values = [x["r_residual"] for x in nights]
                rec.update(mjd_first=nights[0]["mjd_day"], mjd_last=nights[-1]["mjd_day"],
                           p95=round(np.percentile(values, 95), 5), faintest=round(max(values), 5),
                           n_0p2=sum(v >= .2 for v in values),
                           n_0p2_score3=sum(x["r_residual"] >= .2 and x["local_score"] >= 3 for x in nights),
                           n_0p2_score5=sum(x["r_residual"] >= .2 and x["local_score"] >= 5 for x in nights))
            summary.append(rec)
            all_nights.extend(nights)
    write(ROOT / "nightly.csv", "label mjd_day r_residual err_proxy oids unique_times scatter_proxy local_score", all_nights)
    write(ROOT / "nightly_summary.csv", "label source raw nights mjd_first mjd_last p95 faintest n_0p2 n_0p2_score3 n_0p2_score5", summary)


if __name__ == "__main__":
    main()
