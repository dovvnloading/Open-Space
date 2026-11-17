"""Freeze program 8142 mosaics and the two F2550W standards reused from PR 142."""

import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import parse_qs, urlencode, urlparse
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parents[1] / "jwst-nearby-companions" / "g29"


def sha(blob):
    return hashlib.sha256(blob).hexdigest()


def fetch(row, cache):
    path = cache / row["filename"]
    if not path.exists():
        with urlopen(row["url"], timeout=120) as response:
            data = response.read()
        if not data.startswith(b"SIMPLE  =") or len(data) % 2880:
            raise ValueError(f"Incomplete FITS file: {path.name}")
        path.write_bytes(data)
    data = path.read_bytes()
    digest = sha(data)
    if row.get("sha256", digest) != digest:
        raise ValueError(f"SHA-256 mismatch: {path.name}")
    return dict(row, sha256=digest, bytes=len(data))


def main(cache, freeze):
    cache.mkdir(parents=True, exist_ok=True)
    manifest = HERE / "inputs.json"
    if freeze:
        template = json.loads((PREVIOUS / "query_urls.json").read_text())["science"]
        request = json.loads(parse_qs(urlparse(template).query)["request"][0])
        for condition in request["params"]["filters"]:
            if condition["paramName"] == "proposal_id":
                condition["values"] = ["8142"]
        url = "https://mast.stsci.edu/api/v0/invoke?" + urlencode(
            {"request": json.dumps(request)}
        )
        with urlopen(url, timeout=120) as response:
            raw = response.read()
        result = json.loads(raw)
        if result["status"] != "COMPLETE" or len(result["data"]) >= request["pagesize"]:
            raise ValueError("Incomplete CAOM response")
        (HERE / "receipts/caom.json").write_bytes(raw)
        receipt = dict(
            url=url,
            rows=len(result["data"]),
            sha256=sha(raw),
            retrieved_utc=datetime.now(timezone.utc).isoformat(),
        )
        (HERE / "receipts/query.json").write_text(json.dumps(receipt, indent=2) + "\n")
        science = []
        for row in result["data"]:
            if (
                row["target_name"].upper().replace("-", "").replace(" ", "")
                != "HD53143"
            ):
                continue
            if row["filters"] != "F2550W" or not row["dataURL"].endswith("_i2d.fits"):
                continue
            name = row["dataURL"].split("/")[-1]
            url = "https://mast.stsci.edu/api/v0.1/Download/file?uri=" + row["dataURL"]
            role = "combined" if "-c" in name else "roll"
            science.append(dict(row, filename=name, url=url, role=role))
        if sorted(row["role"] for row in science) != ["combined", "roll", "roll"]:
            raise ValueError(
                f"Expected two rolls and their combined mosaic, found {len(science)}"
            )
        standards = [
            dict(r, role="standard")
            for r in json.loads((PREVIOUS / "inputs.json").read_text())
            if r["target_name"] == "BD+60-1753" and r["filters"] == "F2550W"
        ]
        assert len(standards) == 2
        inputs = sorted(science, key=lambda r: r["t_min"]) + sorted(
            standards, key=lambda r: r["t_min"]
        )
    else:
        inputs = json.loads(manifest.read_text())
    with ThreadPoolExecutor(max_workers=4) as pool:
        checked = list(pool.map(lambda row: fetch(row, cache), inputs))
    if freeze:
        manifest.write_text(json.dumps(checked, indent=2) + "\n")
    return checked


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    main(args.cache, args.freeze)
