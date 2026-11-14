"""Fetch the pinned OGLE inputs for population.py into a separate directory."""
import argparse
import hashlib
import json
import shutil
import tarfile
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("destination", type=Path)
    args = ap.parse_args()
    args.destination.mkdir(parents=True, exist_ok=True)
    sources = json.loads((ROOT / "sources.json").read_text())
    for name in ("iv.tar.gz", "iii.tar.gz", "iv_ecl.dat", "iv_remarks.txt", "iv_ident.dat", "iii_ident.dat"):
        record = sources[name]
        target = args.destination / name
        if not target.exists():
            partial = target.with_suffix(target.suffix + ".part")
            with urlopen(record["url"], timeout=180) as response, partial.open("wb") as out:
                shutil.copyfileobj(response, out)
            if hashlib.sha256(partial.read_bytes()).hexdigest() != record["sha256"]:
                raise ValueError(f"Changed upstream input: {name}")
            partial.replace(target)
        if hashlib.sha256(target.read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError(f"Changed input: {name}")
        if name.endswith(".tar.gz"):
            destination = args.destination / name.split(".")[0]
            with tarfile.open(target) as archive:
                for entry in archive:
                    path = Path(entry.name)
                    if (entry.isfile() and len(path.parts) == 3 and path.parts[:2] in
                            (("phot", "I"), ("phot", "V")) and path.name.startswith("OGLE-SMC-") and
                            path.suffix == ".dat"):
                        output = destination / path
                        output.parent.mkdir(exist_ok=True, parents=True)
                        with archive.extractfile(entry) as src, output.open("wb") as dst:
                            shutil.copyfileobj(src, dst)


if __name__ == "__main__":
    main()
