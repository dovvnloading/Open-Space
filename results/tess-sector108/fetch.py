"""Download the official SPOC manifest, intersect with a frozen catalog, and fetch FITS."""

import csv
import hashlib
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

CATALOG_URL = "https://raw.githubusercontent.com/eumemic/exoplanets/e2eee928fd621f75c1621f6b4206d34bb222ae59/candidates/ephemerides.csv"
MANIFEST_URL = "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_108_lc.sh"
ROOT = Path(__file__).resolve().parent


def get(url, path):
    if not path.exists():
        data = urlopen(url, timeout=90).read()
        path.write_bytes(data)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    expected = json.loads((ROOT / "input_hashes.json").read_text()).get(url)
    if expected is not None and digest != expected:
        raise ValueError(
            f"Source revision changed for {url}. Expected {expected}, obtained {digest}."
        )
    return digest


def download(row):
    tic, name = row
    padded = f"{tic:016}"
    directory = "/".join(padded[i : i + 4] for i in range(0, 16, 4))
    url = f"https://stpubdata.s3.amazonaws.com/tess/public/tid/s0108/{directory}/{name}"
    path = ROOT / "fits" / name
    return {"tic": tic, "name": name, "url": url, "sha256": get(url, path)}


def main():
    (ROOT / "fits").mkdir(exist_ok=True)
    catalog_hash = get(CATALOG_URL, ROOT / "ephemerides.csv")
    archive_hash = get(MANIFEST_URL, ROOT / "targets.sh.txt")
    text = (ROOT / "targets.sh.txt").read_text()
    names = re.findall(r"-o (tess\d+-s0108-(\d{16})-\d+-s_lc\.fits) ", text)
    by_id = {int(tic): name for name, tic in names}
    if len(by_id) != len(names):
        raise ValueError("Duplicate 2-minute products; inspect the archive manifest.")
    rows = list(csv.DictReader((ROOT / "ephemerides.csv").read_text().splitlines()))
    ids = sorted({int(r["TIC"]) for r in rows} & by_id.keys())
    with ThreadPoolExecutor(4) as pool:
        files = list(pool.map(download, [(i, by_id[i]) for i in ids]))
    provenance = {
        "catalog_url": CATALOG_URL,
        "catalog_sha256": catalog_hash,
        "n_catalog_rows": len(rows),
        "n_catalog_stars": len({r["TIC"] for r in rows}),
        "target_manifest_url": MANIFEST_URL,
        "target_manifest_sha256": archive_hash,
        "n_2min_targets": len(names),
        "selected": files,
    }
    (ROOT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print([(f["tic"], f["sha256"][:8]) for f in files])


if __name__ == "__main__":
    main()
