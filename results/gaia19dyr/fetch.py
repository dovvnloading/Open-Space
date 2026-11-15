"""Download public source bytes, check hashes, and rebuild the nightly data."""
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import urllib.request

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def reduce_photometry():
    rows, objects = [], []
    for item in json.loads((HERE / "sources.json").read_text()):
        with urllib.request.urlopen(item["url"], timeout=180) as response:
            data = response.read()
        if sha256(data).hexdigest() != item["sha256"]:
            raise ValueError(f"Archive response has changed for {item['star']} {item['archive']}")
        star = item["star"]
        if item["archive"] == "Gaia":
            table = pd.read_csv(BytesIO(data), skiprows=1)
            table["mag"] = pd.to_numeric(table.averagemag, errors="coerce")
            table = table.dropna(subset="mag").copy()
            table["night"] = np.floor(table["JD(TCB)"] - 2400000.5).astype(int)
            baseline = table.mag.quantile(.2)
            for night, group in table.groupby("night"):
                rows.append((star, "G", night, group["JD(TCB)"].median(),
                             group.mag.median(), group.mag.median() - baseline, len(group)))
            continue
        table = pd.read_csv(BytesIO(data)).query("catflags == 0 and magerr < .15").copy()
        for (band, oid), group in table.groupby(["filtercode", "oid"]):
            objects.append((star, band, str(oid), len(group), group.ra.median(),
                            group.dec.median(), group.mag.quantile(.2)))
        table["delta"] = table.mag - table.groupby("oid").mag.transform(lambda m: m.quantile(.2))
        table["night"] = np.floor(table.mjd).astype(int)
        for (band, night), group in table.groupby(["filtercode", "night"]):
            rows.append((star, band, night, group.hjd.median(), group.mag.median(),
                         group.delta.median(), len(group)))
    columns = ["star", "filter", "mjd_night", "jd", "mag", "delta_mag", "n"]
    pd.DataFrame(rows,columns=columns).to_csv(HERE / "nightly.csv", index=False, float_format="%.8f")
    columns = ["star", "filter", "oid", "n", "ra_median", "dec_median", "mag_p20"]
    pd.DataFrame(objects,columns=columns).to_csv(HERE / "objects.csv", index=False, float_format="%.8f")


if __name__ == "__main__":
    reduce_photometry()
