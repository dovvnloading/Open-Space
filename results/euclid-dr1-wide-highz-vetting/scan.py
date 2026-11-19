#!/usr/bin/env python3
"""Full wide HATS Parquet streaming screen; executed on networked GitHub Actions only."""
import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
from xml.etree import ElementTree as ET

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.fs as pafs
import requests
from astropy import units as u
from astropy_healpix import HEALPix

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
BATCHES = OUT / "batches"
BASE = "https://nasa-irsa-euclid-dr1.s3.us-east-1.amazonaws.com/"
BUCKET = "nasa-irsa-euclid-dr1"
PREFIX = ("contributed/dr1/hats/reformatted_catalogs/mer_catalogue_wide_survey/"
          "euclid_dr1_reformatted_catalogs_mer_catalogue_wide_survey-hats/")
INDEX_URL = BASE + PREFIX + "partition_info.csv"
FROZEN_INDEX_SHA256 = "f7940bcaed19b6f67295c5592cc81783d798c0f858062d7293e5e3dfbfa57c2a"
PUBLISHED_ROWS = 888152165
BATCH_SIZE = 24
BASE_COLUMNS = ("ra dec object_id tileid vis_det point_like_prob spurious_prob "
                "ellipticity mumax_minus_mag det_quality_flag flag_y flag_j flag_h "
                "flag_vis flux_vis_3fwhm_aper fluxerr_vis_3fwhm_aper "
                "flux_y_3fwhm_aper fluxerr_y_3fwhm_aper "
                "flux_j_3fwhm_aper fluxerr_j_3fwhm_aper "
                "flux_h_3fwhm_aper fluxerr_h_3fwhm_aper").split()
EXT_BANDS = ("g_ext_decam r_ext_decam i_ext_decam z_ext_decam "
             "i_ext_panstarrs z_ext_panstarrs").split()
COLUMNS = BASE_COLUMNS + [f"{pre}_{band}_3fwhm_aper" for band in EXT_BANDS
                          for pre in ("flux", "fluxerr")]
REQUIRED = set(("ra dec object_id vis_det ellipticity flag_j flag_h "
                "flux_vis_3fwhm_aper fluxerr_vis_3fwhm_aper "
                "flux_y_3fwhm_aper fluxerr_y_3fwhm_aper "
                "flux_j_3fwhm_aper fluxerr_j_3fwhm_aper "
                "flux_h_3fwhm_aper fluxerr_h_3fwhm_aper spurious_prob").split())
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "Open-Space-DR1-wide-source-vetting/1.0"})
HP_CACHE = {}

def utc():
    return datetime.now(timezone.utc).isoformat()

def sha(data):
    return hashlib.sha256(data).hexdigest()

def fetch(url, timeout=65):
    for k in range(4):
        try:
            r = SESSION.get(url, timeout=timeout)
            r.raise_for_status()
            return r
        except requests.RequestException:
            if k == 3:
                raise
            time.sleep(2 ** k)

def catalogue_index():
    OUT.mkdir(exist_ok=True)
    response = fetch(INDEX_URL)
    raw = response.content
    digest = sha(raw)
    # Fail closed when the publisher's catalogue changes under a frozen scan.
    if digest != FROZEN_INDEX_SHA256:
        raise RuntimeError("Publisher partition index changed: " + digest)
    (OUT / "partition_info.csv").write_bytes(raw)
    items = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    fields = {"Norder", "Npix", "Dir", "num_rows"}
    if not items or not fields.issubset(items[0]):
        raise RuntimeError("Unexpected publisher HATS index schema")
    for x in items:
        for key in ("Norder", "Npix", "num_rows"):
            x[key] = int(x[key])
    items.sort(key=lambda x: (x["Norder"], x["Npix"]))
    n = sum(x["num_rows"] for x in items)
    if n != PUBLISHED_ROWS:
        raise RuntimeError(f"Index total rows changed: {n} != {PUBLISHED_ROWS}")
    return items

def load_controls():
    out = []
    for filename, typ in (("controls.tsv", "quasar"), ("dwarfs.tsv", "dwarf")):
        with (ROOT / filename).open() as f:
            for row in csv.DictReader(f, delimiter="\t"):
                name = row["name"]
                m = re.search(r"J(\d{2})(\d{2})(\d{2}\.\d+)([+-])(\d{2})(\d{2})(\d{2}\.\d+)", name)
                if not m:
                    raise ValueError("Invalid control name " + name)
                hh, mm, ss, sign, dd, dm, ds = m.groups()
                ra = (int(hh) + int(mm)/60 + float(ss)/3600)*15
                dec = (int(dd) + int(dm)/60 + float(ds)/3600)*(1 if sign == "+" else -1)
                out.append({"name": name, "class": typ, "ra": ra, "dec": dec,
                            "z": float(row.get("z", "nan")) if typ == "quasar" else None})
    return out

def hp_for(order, controls):
    if order not in HP_CACHE:
        hp = HEALPix(nside=2**order, order="nested", frame="icrs")
        ra = np.array([c["ra"] for c in controls])
        dec = np.array([c["dec"] for c in controls])
        HP_CACHE[order] = hp.lonlat_to_healpix(ra*u.deg, dec*u.deg)
    return HP_CACHE[order]

def as_float(tab, key):
    if key not in tab.column_names:
        return np.full(tab.num_rows, np.nan, dtype=np.float64)
    return np.asarray(tab[key].to_numpy(zero_copy_only=False), dtype=np.float64)

def mask_rows(tab):
    v = as_float(tab, "flux_vis_3fwhm_aper")
    ve = as_float(tab, "fluxerr_vis_3fwhm_aper")
    y = as_float(tab, "flux_y_3fwhm_aper")
    ye = as_float(tab, "fluxerr_y_3fwhm_aper")
    j = as_float(tab, "flux_j_3fwhm_aper")
    je = as_float(tab, "fluxerr_j_3fwhm_aper")
    h = as_float(tab, "flux_h_3fwhm_aper")
    he = as_float(tab, "fluxerr_h_3fwhm_aper")
    spur = as_float(tab, "spurious_prob")
    visdet = as_float(tab, "vis_det")
    ell = as_float(tab, "ellipticity")
    fj = as_float(tab, "flag_j")
    fh = as_float(tab, "flag_h")
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        broad = ((visdet == 0) & (j > .65) & (h > 0) &
                 (y < .85*j) & (j > .65*h) & (je > 0) & (he > 0) &
                 (~np.isfinite(spur) | (spur < .7)))
        yj = -2.5*np.log10(y/j)
        jh = -2.5*np.log10(j/h)
        y_lower_limit = -2.5*np.log10(np.maximum(2*ye, y)/j)
        # The pilot uses a detected-flux colour for y>0, and a 2-sigma
        # lower limit only when y<=0. Keep that selection unchanged.
        ypass = np.where(y > 0, yj >= .40, y_lower_limit >= .40)
        strict = (broad & (j/je >= 6) & (h/he >= 4) &
                  (ve > 0) & (np.abs(v/ve) < 2) & ypass & (jh <= .45) &
                  np.isfinite(ell) & (ell < .30) &
                  np.isfinite(fj) & (fj == 0) &
                  np.isfinite(fh) & (fh == 0))
    return broad, strict, yj, jh, y_lower_limit

def nearest(ra, dec, co, radius_arcsec=3):
    ddec = np.abs(dec - co["dec"])
    close = ddec <= radius_arcsec/3600
    dra = np.abs((ra-co["ra"]+180) % 360 - 180)
    close &= dra <= radius_arcsec/3600 / max(.01, math.cos(math.radians(co["dec"])))
    ii = np.flatnonzero(close)
    if len(ii) == 0:
        return []
    lat1 = np.radians(dec[ii])
    lat2 = math.radians(co["dec"])
    dlon = np.radians(ra[ii]-co["ra"])
    a = np.sin((lat1-lat2)/2)**2 + np.cos(lat1)*math.cos(lat2)*np.sin(dlon/2)**2
    dist = np.degrees(2*np.arcsin(np.sqrt(np.clip(a, 0, 1))))*3600
    return sorted([(int(i), float(s)) for i, s in zip(ii, dist) if s <= radius_arcsec],
                  key=lambda x: x[1])

def scalar(v):
    if v is None:
        return None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    if isinstance(v, float) and not math.isfinite(v):
        return None
    return v

def row_dict(tab, i, yj, jh, yl):
    r = {k: scalar(v) for k,v in tab.slice(i, 1).to_pylist()[0].items()}
    r.update({"Y_minus_J_AB": scalar(yj[i]), "J_minus_H_AB": scalar(jh[i]),
              "Y_minus_J_2sig_limit": scalar(yl[i])})
    r["id"] = str(r["object_id"])
    return r

def get_partition_path(record):
    folder = (PREFIX + f"dataset/Norder={record['Norder']}/"
              f"Dir={record['Dir']}/Npix={record['Npix']}/")
    listing = fetch(BASE + "?list-type=2&prefix=" + quote(folder, safe="/"))
    xml = ET.fromstring(listing.content)
    keys = [x.find("{*}Key").text for x in xml.findall("{*}Contents")
            if x.find("{*}Key") is not None and x.find("{*}Key").text.endswith(".parquet")]
    if len(keys) != 1:
        raise ValueError("Expected one Parquet object, found " + str(keys[:5]))
    key = keys[0]
    head = SESSION.head(BASE + key, timeout=45)
    head.raise_for_status()
    return key, {"key": key, "bytes": int(head.headers.get("Content-Length", "0")),
                 "etag": head.headers.get("ETag"), "list_sha256": sha(listing.content)}

def scan_partition(record, fs, controls):
    order, pix = record["Norder"], record["Npix"]
    key, receipt = get_partition_path(record)
    with fs.open_input_file(BUCKET + "/" + key) as stream:
        pf = pq.ParquetFile(stream)
        cols = [x for x in COLUMNS if x in pf.schema_arrow.names]
        missing = REQUIRED - set(cols)
        if missing:
            raise RuntimeError("Missing required MER columns: " + ",".join(sorted(missing)))
        if pf.metadata.num_rows != record["num_rows"]:
            raise RuntimeError("Publisher partition rows don't match Parquet footer")
        pixcontrols = [c for c,q in zip(controls, hp_for(order, controls)) if q == pix]
        result_rows, result_controls = [], []
        count_broad = count_strict = count_blue = 0
        for rg in range(pf.num_row_groups):
            tab = pf.read_row_group(rg, columns=cols, use_threads=False)
            broad, strict, yj, jh, yl = mask_rows(tab)
            count_broad += int(broad.sum())
            count_strict += int(strict.sum())
            count_blue += int((strict & (jh < -.5)).sum())
            selected = np.flatnonzero(strict)
            for i in selected:
                r = row_dict(tab, int(i), yj, jh, yl)
                r["Norder"] = order
                r["Npix"] = pix
                r["morphology_diagnostic"] = bool(
                    r.get("mumax_minus_mag") is not None and
                    -3.25 <= r["mumax_minus_mag"] <= -2.65 and
                    r.get("ellipticity") is not None and r["ellipticity"] < .2)
                result_rows.append(r)
            if pixcontrols:
                ra = as_float(tab, "ra")
                dec = as_float(tab, "dec")
                for c in pixcontrols:
                    for i, dist in nearest(ra, dec, c):
                        r = row_dict(tab, i, yj, jh, yl)
                        result_controls.append({
                            "name": c["name"], "class": c["class"],
                            "z": c["z"], "object_id": str(r["object_id"]),
                            "sep_arcsec": round(dist, 5),
                            "broad_pass": bool(broad[i]),
                            "strict_pass": bool(strict[i]),
                            "ra": r["ra"], "dec": r["dec"],
                            "J_minus_H_AB": r["J_minus_H_AB"],
                            "Y_minus_J_AB": r["Y_minus_J_AB"],
                            "mumax_minus_mag": r.get("mumax_minus_mag"),
                            "ellipticity": r["ellipticity"]})
        receipt.update({"Norder": order, "Npix": pix, "rows": record["num_rows"],
                        "row_groups": pf.num_row_groups, "column_count": len(cols),
                        "columns_sha256": sha("\n".join(cols).encode()),
                        "broad": count_broad, "strict": count_strict,
                        "blue_JH_below_minus_0p5": count_blue,
                        "controls": result_controls})
    return receipt, result_rows

def main():
    BATCHES.mkdir(parents=True, exist_ok=True)
    items = catalogue_index()
    controls = load_controls()
    n_batches = math.ceil(len(items) / BATCH_SIZE)
    done = []
    for n in range(n_batches):
        f = BATCHES / f"batch_{n:03d}_manifest.json"
        if f.exists():
            m = json.loads(f.read_text())
            if m["status"] != "COMPLETE" or m["batch"] != n:
                raise RuntimeError("Invalid saved batch manifest")
            data = BATCHES / f"batch_{n:03d}_candidates.jsonl.gz"
            if sha(data.read_bytes()) != m["candidates_sha256"]:
                raise RuntimeError("Batch candidate output hash mismatch " + str(n))
            done.append(n)
    n = next((x for x in range(n_batches) if x not in done), None)
    if n is None:
        state = {"phase": "SCANNED_ALL_PARTITIONS", "time_utc": utc(),
                 "published_rows": PUBLISHED_ROWS, "total_partitions": len(items),
                 "completed_batches": len(done), "n_batches": n_batches}
        (OUT / "status.json").write_text(json.dumps(state, indent=2)+"\n")
        print("COMPLETE_SCAN", json.dumps(state), flush=True)
        return
    selected = items[n*BATCH_SIZE:(n+1)*BATCH_SIZE]
    fs = pafs.S3FileSystem(anonymous=True, region="us-east-1", request_timeout=90,
                          connect_timeout=15, proxy_options=os.environ.get("HTTPS_PROXY") or None)
    rows, receipts = [], []
    try:
        for k, entry in enumerate(selected):
            receipt, part_rows = scan_partition(entry, fs, controls)
            rows.extend(part_rows)
            receipts.append(receipt)
            print(f"SCAN batch {n} partition {k+1}/{len(selected)} "
                  f"({entry['Norder']},{entry['Npix']}) broad={receipt['broad']} "
                  f"strict={receipt['strict']}", flush=True)
    except Exception as e:
        state = {"phase": "BLOCKED_PARTITION_FAILURE", "batch": n,
                 "read_partitions_before_failure": len(receipts),
                 "failed": {"Norder": entry["Norder"], "Npix": entry["Npix"],
                            "error": type(e).__name__ + ": " + str(e)[:500]},
                 "time_utc": utc()}
        (OUT / "status.json").write_text(json.dumps(state, indent=2)+"\n")
        print("BLOCKED_SCAN", json.dumps(state), flush=True)
        raise
    data = BATCHES / f"batch_{n:03d}_candidates.jsonl.gz"
    with gzip.open(data, "wt", encoding="utf-8", compresslevel=6) as f:
        for r in rows:
            f.write(json.dumps(r, sort_keys=True, allow_nan=False)+"\n")
    manifest = {"batch": n, "status": "COMPLETE", "time_utc": utc(),
                "partitions": receipts, "candidate_count": len(rows),
                "advertised_rows": sum(x["rows"] for x in receipts),
                "candidates_sha256": sha(data.read_bytes()),
                "index_sha256": FROZEN_INDEX_SHA256,
                "pipeline_sha256": sha(Path(__file__).read_bytes())}
    (BATCHES / f"batch_{n:03d}_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True)+"\n")
    state = {"phase": "SCANNING", "time_utc": utc(), "batch_just_finished": n,
             "published_rows": PUBLISHED_ROWS, "total_partitions": len(items),
             "completed_batches": len(done)+1, "n_batches": n_batches}
    (OUT / "status.json").write_text(json.dumps(state, indent=2)+"\n")
    print("BATCH_COMPLETE", json.dumps({"batch": n, "candidate_count": len(rows),
          "broad": sum(x["broad"] for x in receipts),
          "rows": manifest["advertised_rows"]}), flush=True)

if __name__ == "__main__":
    main()
