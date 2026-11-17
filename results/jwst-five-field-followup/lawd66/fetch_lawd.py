"""Freeze or replay LAWD 66 STScI Level 3 catalogues and image products."""

import argparse
import concurrent.futures
import datetime
import hashlib
import json
import re
import sys
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("repo", type=Path)
p.add_argument("cache", type=Path)
p.add_argument("--freeze", action="store_true")
args = p.parse_args()
sys.path.insert(0, str(args.repo / "results/jwst-nearby-companions/miri"))
from fetch_inputs import download, product_url

out = Path(__file__).resolve().parent
args.cache.mkdir(exist_ok=True, parents=True)
manifest = out / "inputs.json"


def inventory(pid):
    prefix = f"jwst/public/jw{pid}/L3/"
    url = "https://s3.amazonaws.com/stpubdata?" + urlencode(
        {"list-type": 2, "max-keys": 1000, "prefix": prefix}
    )
    with urlopen(url, timeout=70) as res:
        raw = res.read()
    root = ET.fromstring(raw)
    if (
        root.findtext("{*}IsTruncated") != "false"
        or root.findtext("{*}Prefix") != prefix
    ):
        raise ValueError(pid)
    filename = f"stpubdata_{pid}.xml"
    (out / filename).write_bytes(raw)
    rows = []
    for item in root.findall("{*}Contents"):
        name = item.findtext("{*}Key").rsplit("/", 1)[-1]
        m = re.fullmatch(
            r"(jw\d+-o\d+_t001_nircam_[a-z0-9-]+)_(cat.ecsv|i2d.fits)", name
        )
        if not m:
            continue
        rows.append(
            dict(
                obs_id=m[1],
                filename=name,
                url=product_url(m[1], m[2]),
                bytes=int(item.findtext("{*}Size")),
                etag=item.findtext("{*}ETag"),
                product_updated_utc=item.findtext("{*}LastModified"),
            )
        )
    return rows, dict(
        url=url,
        file=filename,
        sha256=hashlib.sha256(raw).hexdigest(),
        retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    )


if args.freeze:
    if manifest.exists():
        raise FileExistsError(manifest)
    with concurrent.futures.ThreadPoolExecutor(3) as pool:
        result = list(pool.map(inventory, ("02692", "04519", "06584")))
    rows = sorted(
        [r for part, _ in result for r in part],
        key=lambda r: (r["filename"].endswith(".fits"), r["filename"]),
    )
    if len(rows) != 32:
        raise ValueError(len(rows))
    mast_receipts = []
    mast = {}
    from fetch_inputs import ENDPOINT

    for pid in ("2692", "4519", "6584"):
        columns = "obsid,obs_id,proposal_id,target_name,instrument_name,t_min,t_max,filters,dataURL,dataRights"
        terms = [
            {"paramName": "proposal_id", "values": [pid]},
            {"paramName": "obs_collection", "values": ["JWST"]},
        ]
        query = {
            "service": "Mast.Caom.Filtered",
            "params": {"columns": columns, "filters": terms},
            "format": "json",
            "pagesize": 5000,
            "page": 1,
        }
        url = ENDPOINT + urlencode({"request": json.dumps(query)})
        fname = f"mast_{pid}.json"
        if not (out / fname).exists():
            (out / fname).write_bytes(urlopen(url, timeout=180).read())
        raw = (out / fname).read_bytes()
        reply = json.loads(raw)
        assert reply["status"] == "COMPLETE" and len(reply["data"]) < 5000
        mast.update(
            {
                r["obs_id"]: r
                for r in reply["data"]
                if r["dataRights"] == "PUBLIC" and r["target_name"] == "LAWD-66-SOURCE"
            }
        )
        mast_receipts.append(
            dict(url=url, file=fname, sha256=hashlib.sha256(raw).hexdigest())
        )
    for r in rows:
        m = mast[r["obs_id"]]
        r.update(
            obsid=m["obsid"],
            t_min=m["t_min"],
            t_max=m["t_max"],
            mast_filters=m["filters"],
            data_rights=m["dataRights"],
            mast_image_uri=m["dataURL"],
        )
    data = dict(
        mast_queries=mast_receipts,
        source="STScI public cloud distribution verified against MAST program and visit metadata",
        receipts=[s for _, s in result],
    )
else:
    data = json.loads(manifest.read_text())
    rows = data.pop("inputs")


def one(r):
    ret = download(r, args.cache)
    if ret["bytes"] != r["bytes"]:
        raise ValueError(f"truncated input {r['filename']}")
    # These STScI S3 multipart objects use 50 MiB parts. This independently
    # checks the frozen bytes against the service's recorded ETag.
    with (args.cache / r["filename"]).open("rb") as f:
        parts = [hashlib.md5(b).digest() for b in iter(lambda: f.read(50 * 2**20), b"")]
    etag = hashlib.md5(b"".join(parts)).hexdigest() + "-" + str(len(parts))
    if etag != r["etag"].strip('"'):
        raise ValueError(f"Etag mismatch for {r['filename']}")
    print(r["filename"], r["bytes"], flush=True)
    return ret


with concurrent.futures.ThreadPoolExecutor(5) as pool:
    results = list(pool.map(one, rows))
if args.freeze:
    manifest.write_text(json.dumps(dict(data, inputs=results), indent=2) + "\n")
