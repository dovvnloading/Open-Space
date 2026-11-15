"""Download the five catalog inputs used by analyze.py, or reuse cached files."""
from pathlib import Path
from urllib.request import urlopen
import hashlib
import json
import sys

BASE = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
MAST = "https://archive.stsci.edu/hlsps/phangs-cat/dr5/hlsp_phangs-cat_hst_acs-uvis_dr5-targets_v1_"
URLS = {
    "m31.tsv": BASE + "-source=J/ApJ/824/42&-out.max=unlimited",
    "pne.tsv": BASE + "-out.all&-source=J/MNRAS/511/6087&-out.max=unlimited",
    "congiu.tsv": BASE + "-out.all&-source=J/A%2bA/672/A148&-out.max=unlimited",
    "human.fits": MAST + "human-class1+2-primary.fits",
    "machine.fits": MAST + "machine-class1+2-primary.fits",
}

def digest(path):
    data = path.read_bytes()
    if path.suffix == ".tsv":
        # VizieR adds the query date and service metadata to the response.
        data = b"\n".join(line for line in data.splitlines() if line and not line.startswith(b"#"))
    return hashlib.sha256(data).hexdigest()

if __name__ == "__main__":
    directory = Path(sys.argv[1])
    directory.mkdir(exist_ok=True, parents=True)
    expected = json.loads(Path(__file__).with_name("inputs.json").read_text())
    for name, url in URLS.items():
        path = directory / name
        if not path.exists():
            with urlopen(url, timeout=120) as response:
                path.write_bytes(response.read())
        if digest(path) != expected[name]["sha256"]:
            raise ValueError(f"Input changed: {path}. Compare catalog revisions.")
        print(name, "verified")
