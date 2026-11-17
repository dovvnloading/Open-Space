"""Recover the three missing warm-sample ZTF DR24 cones and compute unique-night fading."""
from concurrent.futures import ThreadPoolExecutor
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
from urllib.request import urlopen

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "optical"))
import hats_fallback as hats
import nightly


def retry_api(receipt):
    meta = dict(label=receipt["label"], url=receipt["url"], utc=dt.datetime.now(dt.timezone.utc).isoformat())
    try:
        with urlopen(receipt["url"], timeout=90) as f:
            data = f.read()
        reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
        if "mjd" not in (reader.fieldnames or []):
            raise ValueError(data[:200])
        rows = list(reader)
        (ROOT / f"retry_{receipt['label']}.csv.gz").write_bytes(gzip.compress(data, mtime=0))
        meta.update(rows=len(rows), sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))
    except Exception as exc:
        meta["error"] = str(exc)
    return meta


def selections(t, rows, method):
    points = nightly.collect(t, rows, catalog_positions=(method == "HATS"))
    result = dict(source_id=t["source_id"], method=method, raw=len(rows), nights=len(points))
    if points:
        vals = [v["r_residual"] for v in points]
        result.update(mjd_first=points[0]["mjd_day"], mjd_last=points[-1]["mjd_day"],
                      p95=round(np.percentile(vals,95),5), faintest=round(max(vals),5),
                      n_0p2=sum(v >= .2 for v in vals),
                      n_0p2_score3=sum(v["r_residual"] >= .2 and v["local_score"] >= 3 for v in points))
    return result


def main():
    failed = [r for r in json.loads((ROOT / "ztf_receipts.json").read_text()) if "error" in r]
    assert len(failed) == 3
    ids = {r["label"] for r in failed}
    with (ROOT.parent / "literature" / "fgk_astrometry.csv").open() as stream:
        targets = {r["source_id"]: dict(r,label=r["source_id"]) for r in csv.DictReader(stream) if r["source_id"] in ids}
    hats.ROOT = ROOT
    H, results = [], []
    with ThreadPoolExecutor(3) as pool:
        retry = [pool.submit(retry_api,r) for r in failed]
        try:
            remote = hats.fs.S3FileSystem(anonymous=True, region=hats.fs.resolve_s3_region(hats.BUCKET))
            index_path = f"{hats.CAT}/partition_info.csv"
            with remote.open_input_stream(index_path) as stream:
                index = stream.read()
            leaves = {(int(r["Norder"]),int(r["Npix"])) for r in csv.DictReader(io.StringIO(index.decode()))}
            H_index = dict(url="s3://"+index_path, sha256=hashlib.sha256(index).hexdigest())
            for sid,t in targets.items():
                try:
                    H.append(hats.retrieve(t,remote,leaves))
                except Exception as exc:
                    H.append(dict(label=sid,error=str(exc)))
        except Exception as exc:
            H_index = dict(error=str(exc))
            H = [dict(label=s,error="Partition-index retrieval failed") for s in ids]
        primary = [f.result() for f in retry]
    audit = dict(utc=dt.datetime.now(dt.timezone.utc).isoformat(),hats_index=H_index,
                 hats=H,primary_api=primary)
    (ROOT/"recovery_receipts.json").write_text(json.dumps(audit,indent=2)+"\n")
    for method, receipts, prefix in [("API",primary,"retry"),("HATS",H,"hats")]:
        for r in receipts:
            if "error" in r:
                print(method, r["label"], r["error"], flush=True)
                continue
            data=gzip.decompress((ROOT/f"{prefix}_{r['label']}.csv.gz").read_bytes())
            assert hashlib.sha256(data).hexdigest()==r["sha256"]
            rows = list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))
            result = selections(targets[r["label"]],rows,method)
            results.append(result)
            print(result,flush=True)
    nightly.write(ROOT/"recovery_summary.csv","source_id method raw nights mjd_first mjd_last p95 faintest n_0p2 n_0p2_score3",results)
    if {r["source_id"] for r in results} != ids:
        raise ValueError("Some original ZTF failures remain unavailable in both archives")


if __name__ == "__main__":
    main()
