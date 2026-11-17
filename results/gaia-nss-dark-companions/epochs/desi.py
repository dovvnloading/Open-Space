"""Scan the per-exposure DR1 Gaia crossmatches; retrieve RVs and times by row."""

import argparse
import csv
import hashlib
import io
import json
import math
import re
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

from catalog_io import byte_range, read_header, table_dtype, target_ids

BASE = "https://data.desi.lbl.gov/public/dr1/vac/dr1/mws/iron/v1.0/rv_output/240521/"
CORRECTION = "https://zenodo.org/records/15469272/files/backup_correction.fits"
EXPECTED_ROWS = 10012925


def extensions(url):
    found = {}
    position, size = 0, None
    while size is None or position < size:
        block, size = byte_range(url, position, position + 65535)
        buffer = io.BytesIO(block)
        header = read_header(buffer)
        start = position + buffer.tell()
        length = int(header.get("NAXIS1", 0)) * int(header.get("NAXIS2", 0))
        if length:
            found[header["EXTNAME"]] = (header, start)
        position = start + ((length + 2879) // 2880) * 2880
    return found


def scan_product(name, ids):
    url = BASE + name + ".fits"
    extensions_by_name = extensions(url)
    gaia, start = extensions_by_name["GAIA"]
    gaia_type = table_dtype(gaia)
    n = int(gaia["NAXIS2"])
    step = max(1, 8 * 1024**2 // gaia_type.itemsize)

    def scan(first):
        stop = min(first + step, n)
        data, _ = byte_range(
            url,
            start + first * gaia_type.itemsize,
            start + stop * gaia_type.itemsize - 1,
        )
        entries = np.frombuffer(data, gaia_type)
        selected = np.flatnonzero(np.isin(entries["SOURCE_ID"], ids))
        return hashlib.sha256(data).hexdigest(), [
            (first + int(i), int(entries["SOURCE_ID"][i])) for i in selected
        ]

    hashes, hits = [], []
    with ThreadPoolExecutor(6) as workers:
        for digest, matches in workers.map(scan, range(0, n, step)):
            hashes.append(digest)
            hits.extend(matches)

    def fetch_entry(ext, index):
        header, offset = extensions_by_name[ext]
        if int(header["NAXIS2"]) != n:
            raise ValueError(f"Unaligned {ext} table")
        dtype = table_dtype(header)
        data, _ = byte_range(
            url,
            offset + index * dtype.itemsize,
            offset + (index + 1) * dtype.itemsize - 1,
        )
        return np.frombuffer(data, dtype)[0], hashlib.sha256(data).hexdigest()

    def fields(hit):
        index, sid = hit
        r, rhash = fetch_entry("RVTAB", index)
        f, fhash = fetch_entry("FIBERMAP", index)
        if r["TARGETID"] != f["TARGETID"] or r["EXPID"] != f["EXPID"]:
            raise ValueError("RVTAB–FIBERMAP identity mismatch")
        return dict(
            source_id=str(sid),
            row=index,
            targetid=str(r["TARGETID"]),
            expid=int(r["EXPID"]),
            mjd=float(f["MJD"]),
            exptime=float(f["EXPTIME"]),
            rv=float(r["VRAD"]),
            err=float(r["VRAD_ERR"]),
            success=r["SUCCESS"].decode(),
            warn=int(r["RVS_WARN"]),
            fiberstatus=int(f["FIBERSTATUS"]),
            ref_id=str(r["REF_ID"]),
            ref_cat=r["REF_CAT"].decode().strip(),
            ra=float(r["TARGET_RA"]),
            dec=float(r["TARGET_DEC"]),
            rvs_sha256=rhash,
            fibermap_sha256=fhash,
        )

    with ThreadPoolExecutor(4) as workers:
        rows = list(workers.map(fields, hits))
    return dict(
        url=url,
        gaia_rows=n,
        gaia_row_width=gaia_type.itemsize,
        gaia_chunk_size_rows=step,
        gaia_hdu_start=start,
        gaia_chunk_sha256=hashes,
        matches=rows,
    )


def correction_table(data):
    digest = hashlib.md5(data).hexdigest()
    if digest != "f48a4b21b541e94d61f4372f4c555f12":
        raise ValueError("Correction file changed")
    f = io.BytesIO(data)
    read_header(f)
    header = read_header(f)
    return np.frombuffer(
        data, dtype=table_dtype(header), count=int(header["NAXIS2"]), offset=f.tell()
    )


def write_results(records, targets, corrections, out):
    rows = []
    corr = correction_table(corrections)
    for record in records:
        name = record["url"].rsplit("/", 1)[1].removesuffix(".fits")
        program = name.split("-")[-1]
        for r in record["matches"]:
            offset = 0.0
            if name == "rvpix_exp-main-backup":
                matches = corr[corr["TARGETID"] == int(r["targetid"])]
                if len(matches) != 1:
                    raise ValueError("Missing or duplicate DESI backup correction")
                offset = float(matches["VRAD_BIAS"][0])
            floor = {"backup": 2.0, "bright": 1.0, "dark": 1.6}.get(program, 2.0)
            rv = r["rv"] - offset
            quality = int(
                r["success"] == "T"
                and r["warn"] == 0
                and r["fiberstatus"] == 0
                and r["err"] > 0
                and math.isfinite(rv)
            )
            rows.append(
                dict(
                    source_id=r["source_id"],
                    survey="DESI DR1",
                    instrument="DESI-" + program,
                    time_jd=round(r["mjd"] + 2400000.5 + r["exptime"] / 172800, 8),
                    time_scale="UTC",
                    rv_kms=rv,
                    err_kms=math.hypot(r["err"], floor),
                    quality=quality,
                    product=name,
                    targetid=r["targetid"],
                    expid=r["expid"],
                    raw_rv_kms=r["rv"],
                    formal_err_kms=r["err"],
                    subtracted_offset_kms=offset,
                    floor_kms=floor,
                    rvs_warn=r["warn"],
                    fiberstatus=r["fiberstatus"],
                )
            )
    rows.sort(key=lambda r: (int(r["source_id"]), r["time_jd"], r["expid"]))
    if len({(r["source_id"], r["expid"]) for r in rows}) != len(rows):
        raise ValueError("Duplicate target/exposure")
    scanned = sum(r["gaia_rows"] for r in records)
    if scanned != EXPECTED_ROWS:
        raise ValueError("Incomplete exposure search")
    keys = [
        "source_id",
        "survey",
        "instrument",
        "time_jd",
        "time_scale",
        "rv_kms",
        "err_kms",
        "quality",
        "product",
        "targetid",
        "expid",
        "raw_rv_kms",
        "formal_err_kms",
        "subtracted_offset_kms",
        "floor_kms",
        "rvs_warn",
        "fiberstatus",
    ]
    with (out / "desi_visits.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    ledger = dict(
        products=records,
        target_count=len(target_ids(targets)),
        target_file_sha256=hashlib.sha256(targets.read_bytes()).hexdigest(),
        n_scanned=scanned,
        n_matched=len(rows),
        n_sources=len({r["source_id"] for r in rows}),
        corrections_url=CORRECTION,
        corrections_sha256=hashlib.sha256(corrections).hexdigest(),
    )
    (out / "desi_depth.json").write_text(json.dumps(ledger, indent=2) + "\n")


def main(targets, out):
    out.mkdir(exist_ok=True, parents=True)
    ids = target_ids(targets)
    with urllib.request.urlopen(BASE, timeout=60) as response:
        names = sorted(
            re.findall(r'href="(rvpix_exp-[^"]+)\.fits"', response.read().decode())
        )
    if len(names) != 18:
        raise ValueError("DESI program file inventory changed")
    with ThreadPoolExecutor(3) as workers:
        records = list(workers.map(lambda name: scan_product(name, ids), names))
    with urllib.request.urlopen(CORRECTION, timeout=120) as response:
        corrections = response.read()
    write_results(records, targets, corrections, out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("targets", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    main(args.targets, args.out)
