"""Retrieve the CDS I/361 field catalogue and its input-catalogue crosswalk."""
import concurrent.futures
import gzip
import hashlib
import json
import os
from pathlib import Path
from urllib.request import Request, urlopen

BASES = [
    "https://cdsarc.cds.unistra.fr/ftp/I/361/",
    "https://cdsarc.cds.unistra.fr/ftp/cats/I/361/",
    "https://vizier.cfa.harvard.edu/ftp/cats/I/361/",
    "https://cdsarc.u-strasbg.fr/ftp/I/361/",
]
FILES = ("ReadMe", "lenscand.dat.gz", "lenscatn.dat.gz")


def retrieve(filename):
    cache = Path(os.environ["GAIA_CACHE"])
    cache.mkdir(parents=True, exist_ok=True)
    dest = cache / filename
    receipt = cache / (filename + ".json")
    if dest.exists() and receipt.exists():
        prior = json.loads(receipt.read_text())
        with dest.open("rb") as src:
            if hashlib.file_digest(src, "sha256").hexdigest() == prior["sha256"]:
                return prior
    errors = []
    for base in BASES:
        url = base + filename
        temp = dest.with_suffix(".part")
        try:
            digest = hashlib.sha256()
            with urlopen(Request(url, headers={"User-Agent": "Open-Space-research/1.0"}), timeout=75) as response, temp.open("wb") as output:
                while chunk := response.read(1024 * 1024):
                    digest.update(chunk)
                    output.write(chunk)
            if filename.endswith(".gz"):
                with gzip.open(temp) as data:
                    if not data.readline(2000).strip():
                        raise ValueError("Empty compressed data")
            elif "lenscand" not in temp.read_text(errors="replace"):
                raise ValueError("Unexpected catalogue dictionary")
            temp.replace(dest)
            row = dict(file=filename, url=url, bytes=dest.stat().st_size, sha256=digest.hexdigest(), errors=errors)
            receipt.write_text(json.dumps(row, indent=2))
            print(filename, row["bytes"], flush=True)
            return row
        except Exception as exc:
            errors.append(dict(url=url, error=str(exc)))
    return dict(file=filename, errors=errors)


if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(3) as pool:
        rows = list(pool.map(retrieve, FILES))
    out = Path(__file__).parent / "evidence"
    out.mkdir(exist_ok=True)
    (out / "fpr-downloads.json").write_text(json.dumps(rows, indent=2) + "\n")
    if any("sha256" not in row for row in rows):
        raise SystemExit("Some CDS files unavailable; inspect input receipts.")
