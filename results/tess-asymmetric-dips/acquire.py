#!/usr/bin/env python3
"""Rebuild the frozen QLP sample and verify FITS bytes against the run receipts."""
import argparse
import csv
import hashlib
import heapq
import io
import json
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

ROOT = "https://archive.stsci.edu/hlsps/qlp"


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def fetch(url):
    # A partial HTTP response can start with a valid FITS header too.
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=50) as response:
                if response.status != 200:
                    raise IOError(f"expected HTTP 200, got {response.status} from {url}")
                payload = response.read()
                advertised = response.headers.get("Content-Length")
                if advertised and int(advertised) != len(payload):
                    raise IOError(f"incomplete response from {url}")
                return payload
        except urllib.error.HTTPError as exc:
            if exc.code < 500 and exc.code != 429:
                raise
            if attempt == 2:
                raise
        except (TimeoutError, OSError):
            if attempt == 2:
                raise
        time.sleep(2 ** attempt)


def select(payload, n):
    rows = csv.DictReader(io.TextIOWrapper(io.BytesIO(payload), encoding="ascii"))
    seen = set()
    total = 0

    def ranking():
        nonlocal total
        for row in rows:
            tic = int(row["#TIC_ID"])
            if tic <= 0:
                raise ValueError(f"invalid TIC {tic}")
            total += 1
            if tic not in seen:
                seen.add(tic)
                yield (sha256(str(tic).encode("ascii")), tic, float(row["RA"]), float(row["DEC"]))

    targets = [dict(tic=tic, ra=ra, dec=dec) for _, tic, ra, dec in heapq.nsmallest(n, ranking())]
    if len(targets) != n:
        raise ValueError("target list is smaller than the registered sample")
    return total, len(seen), targets


def product_url(sector, tic):
    sector_tag = f"s{sector:04d}"
    padded = f"{tic:016d}"
    folders = "/".join(padded[i:i+4] for i in range(0, 16, 4))
    return f"{ROOT}/{sector_tag}/{folders}/hlsp_qlp_tess_ffi_{sector_tag}-{padded}_tess_v01_llc.fits"


def write_json(path, value):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    tmp.replace(path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path, help="fresh working directory for rerun, with space for full FITS")
    p.add_argument("--sample", required=True, type=Path, help="frozen sample.json")
    p.add_argument("--products", required=True, type=Path, help="original products.csv with TIC and SHA-256")
    p.add_argument("--workers", type=int, default=4)
    args = p.parse_args()
    sample = json.loads(args.sample.read_text())
    expected = {int(r["tic"]): r for r in csv.DictReader(args.products.open())}
    args.out.mkdir(parents=True, exist_ok=True)
    fits_dir = args.out / "fits"
    fits_dir.mkdir(exist_ok=True)
    report = {"started_utc": datetime.now(timezone.utc).isoformat(), "requests": []}
    report_path = args.out / "replay.json"
    url = sample["url"]
    try:
        payload = fetch(url)
        digest = sha256(payload)
        count, unique, targets = select(payload, len(sample["targets"]))
        report["requests"].append(dict(url=url, sha256=digest, bytes=len(payload), raw_rows=count, unique=unique))
        if (digest, len(payload), count, unique, targets) != (
            sample["sha256"], sample["bytes"], sample["raw_rows"], sample["unique"], sample["targets"]
        ):
            raise ValueError("archive target list or sample differs from frozen inputs")
        if set(expected) != {t["tic"] for t in targets}:
            raise ValueError("product receipt and sample differ")
        sector = sample["sector"]

        def replay(target):
            tic = target["tic"]
            url = expected[tic].get("url") or product_url(sector, tic)
            receipt = {"tic": tic, "url": url}
            try:
                raw = fetch(url)
                digest = sha256(raw)
                receipt.update(bytes=len(raw), sha256=digest)
                if not raw.startswith(b"SIMPLE  =") or len(raw) % 2880:
                    raise ValueError("invalid FITS file framing")
                if not expected[tic]["sha256"]:
                    raise ValueError("no frozen product hash")
                if digest != expected[tic]["sha256"]:
                    raise ValueError("FITS bytes differ from original run")
                # Publish only verified products. Never overwrite a file with a partial transfer.
                dest = fits_dir / Path(product_url(sector, tic)).name
                tmp = dest.with_suffix(".part")
                tmp.write_bytes(raw)
                tmp.replace(dest)
                receipt["status"] = "verified"
            except (OSError, ValueError) as err:
                receipt.update(status="error", detail=str(err))
            return receipt

        with ThreadPoolExecutor(args.workers) as pool:
            for receipt in pool.map(replay, targets):
                report["requests"].append(receipt)
                write_json(report_path, report)
        if any(r.get("status") == "error" for r in report["requests"]):
            raise ValueError("one or more frozen FITS products could not be verified")
        with (args.out / "scan_manifest.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, ["path", "tic", "sector", "sha256"])
            writer.writeheader()
            for item in targets:
                tic = item["tic"]
                path = fits_dir / Path(product_url(sector, tic)).name
                writer.writerow(dict(path=str(path.resolve()), tic=tic, sector=sector,
                                     sha256=expected[tic]["sha256"]))
    except Exception as exc:
        report["error"] = str(exc)
        raise
    finally:
        write_json(report_path, report)


if __name__ == "__main__":
    main()
