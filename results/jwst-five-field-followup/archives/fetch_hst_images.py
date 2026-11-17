"""Freeze selected HST images for direct checks of JWST peaks."""
import argparse
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
# WEDGEA1.0, five subexposures per CRJ. Faramaz et al. 2019, Table 3.
STIS = ["ocjc" + suffix + "040_crj.fits" for suffix in ("01", "02", "04", "05", "06", "08")]


def products():
    for name in STIS:
        yield dict(name=name, field="HD202628", url="https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/" + name)
    with (HERE / "hst_footprints.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    selected = set()
    for r in rows:
        target, obs = r["target"], r["obs_id"]
        if float(r["t_exptime"] or 0) < 1350:
            continue
        if target == "hd53143" and obs.startswith("oed8"):
            selected.add((target, obs, obs + "_crj.fits"))
        if target in ("hd21997", "hd38206") and obs.startswith("j9f9"):
            selected.add((target, obs, obs[:-1] + "1_drz.fits"))
    for target, obs, name in sorted(selected):
        yield dict(name=name, field=target.upper(), url=f"https://stpubdata.s3.amazonaws.com/hst/public/{obs[:4]}/{obs}/{name}")
    for prop, target, band in [(10599, "hd-53143", "f110w"), (7897, "hd53143", "f160w")]:
        name = f"hlsp_alice_hst_nicmos_{prop}-{target}_{band}_v1_multi-roll.fits"
        yield dict(name=name, field="HD53143",
                   url=f"https://archive.stsci.edu/missions/hlsp/alice/{prop}/{target}/{band}/{name}")


def main(cache, freeze=False):
    cache.mkdir(parents=True, exist_ok=True)
    manifest = HERE / "hst_inputs.json"
    rows = list(products()) if freeze else json.loads(manifest.read_text())
    def one(row):
        path = cache / row["name"]
        if not path.exists():
            with urlopen(row["url"], timeout=90) as response:
                raw = response.read()
            if raw[:8] != b"SIMPLE  " or len(raw) % 2880:
                raise ValueError(f"Invalid FITS {path.name}")
            path.write_bytes(raw)
        raw = path.read_bytes()
        sha = hashlib.sha256(raw).hexdigest()
        if sha != row.get("sha256", sha):
            raise ValueError(f"Changed bytes: {row['name']}")
        return dict(row, sha256=sha, bytes=len(raw))
    with ThreadPoolExecutor(5) as pool:
        result = list(pool.map(one, rows))
    if freeze:
        manifest.write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    main(args.cache, args.freeze)
