"""Fetch the three public OGLE cutouts and check each original SHA-256."""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent


def main():
    inputs = json.loads((ROOT / "provenance.json").read_text())["inputs"]
    for name in ("science.fits", "reference.fits", "difference.fits"):
        path = ROOT / "inputs" / name
        info = inputs[name]
        if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() == info["sha256"]:
            continue
        with urlopen(info["url"], timeout=90) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != info["sha256"]:
            raise ValueError(f"Archive content differs for {name}")
        temp = path.with_suffix(".part")
        temp.write_bytes(data)
        temp.replace(path)


if __name__ == "__main__":
    main()
