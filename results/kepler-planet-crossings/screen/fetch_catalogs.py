"""Fetch the two public inputs for the exploratory Kepler screen."""

import hashlib
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

columns = (
    "kepid,kepoi_name,kepler_name,koi_disposition,koi_period,koi_time0bk,"
    "koi_duration,koi_depth,koi_impact,koi_model_snr,koi_kepmag"
)
urls = {
    "koi.csv": "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"
    + urlencode({"query": "select " + columns + " from cumulative", "format": "csv"}),
    "times.tsv": "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
    "-source=J/ApJS/225/9/table3&-out.all&-out.max=unlimited",
}
for name, url in urls.items():
    with urlopen(url, timeout=180) as response:
        data = response.read()
    if name == "koi.csv" and not data.startswith(b"kepid,"):
        raise ValueError("Expected NASA Exoplanet Archive KOI table")
    if name == "times.tsv" and b"295187" not in data[-1200:]:
        raise ValueError("Incomplete Holczer timing table")
    Path(name).write_bytes(data)
    print(name, hashlib.sha256(data).hexdigest())
