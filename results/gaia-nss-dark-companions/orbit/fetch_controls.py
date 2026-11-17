"""Fetch the four benchmark data products, preserving raw-byte hashes."""
import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
ESA_SHA = "c5c7165bd13ca7bdad0454719b10c5637b4b2d32"
PAPERS = {"BH1": ("2209.06833v3", 41), "BH2": ("2302.07880v2", 48),
          "NS1": ("2402.06722v2", 59)}
QUERY = ("SELECT * FROM gaiadr3.nss_two_body_orbit WHERE source_id IN "
         "(4373465352415301632,5870569352746779008,6328149636482597888)")
GAIA_URL = "https://gea.esac.esa.int/tap-server/tap/sync?" + urlencode(
    dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=QUERY))
ROW = re.compile(r"^\s*(24\d{5}\.\d+)\s+([−\-\d.]+)\s*±\s*([\d.]+)\s+([\w-]+)")


def fetch(cache):
    cache.mkdir(parents=True, exist_ok=True)
    manifest = {}
    def retrieve(name, url):
        path = cache / name
        if not path.exists():
            with urlopen(url, timeout=90) as response:
                path.write_bytes(response.read())
        manifest[name] = dict(url=url, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        return path
    gaia = retrieve("dr3.csv", GAIA_URL)
    # ESA's archive can return rows in any order.
    with gaia.open() as source:
        reader = csv.DictReader(source)
        records = sorted(reader, key=lambda r: r["source_id"])
        if len(records) != 3:
            raise ValueError("Incomplete Gaia control set")
        with (HERE / "nss_controls.csv").open("w") as dest:
            writer = csv.DictWriter(dest, fieldnames=reader.fieldnames)
            writer.writeheader()
            writer.writerows(records)
    all_rows = []
    for name, (paper, expected) in PAPERS.items():
        pdf = retrieve(name.lower()+".pdf", f"https://arxiv.org/pdf/{paper}")
        lines = subprocess.check_output(
            ["pdftotext", "-layout", str(pdf), "-"], text=True).splitlines()
        rows = [[name, t, v.replace("−", "-"), err, instrument]
                for line in lines if (m := ROW.match(line))
                for t, v, err, instrument in [m.groups()]]
        if len(rows) != expected or len({r[1] for r in rows}) != expected:
            raise ValueError(f"Unexpected RV row count for {name}: {len(rows)}")
        all_rows.extend(rows)
    with (HERE / "published_rv.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(["control","hjd_utc","rv_kms","err_kms","instrument"])
        writer.writerows(all_rows)
    for name in ["epochast.dat", "epochrv.dat"]:
        retrieve(name, f"https://raw.githubusercontent.com/esa/gaia-bhthree/{ESA_SHA}/{name}")
    (HERE / "sources.json").write_text(json.dumps(manifest, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cache", type=Path, help="Directory for immutable source files")
    fetch(parser.parse_args().cache)
