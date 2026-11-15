"""Retrieve Kepler-18's short-cadence products and bin isolated transit windows."""

from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from astropy.io import fits

HERE = Path(__file__).resolve().parent
BASE = "https://stpubdata.s3.amazonaws.com/"
PREFIX = "kepler/public/lightcurves/0086/008644288/"
epochs = pd.read_csv(HERE / "epochs.csv")


def extract(key):
    rows = []
    with urlopen(BASE + key, timeout=90) as stream, fits.open(BytesIO(stream.read())) as hdus:
        data = hdus[1].data
        time, flux = data["TIME"], data["PDCSAP_FLUX"]
        keep = np.isfinite(time + flux) & (data["SAP_QUALITY"] == 0)
        time, flux = time[keep], flux[keep]
        if len(time) == 0:
            return rows
        candidates = epochs[epochs.time.between(time.min(), time.max()) & (epochs.planet != "b")]
        for event in candidates.itertuples():
            other = epochs[epochs.planet != event.planet]
            if (abs(other.time - event.time) < 0.4).any():
                continue
            select = abs(time - event.time) < 0.27
            dt, y = time[select] - event.time, flux[select]
            out = abs(dt) > 0.117
            if select.sum() < 600 or out.sum() < 120:
                continue
            trend = np.polyfit(dt[out], y[out], 2)
            y = (y / np.polyval(trend, dt) - 1) * 1e6
            sigma = 1.4826 * np.median(abs(y[out] - np.median(y[out])))
            groups = np.round(dt * 1440 / 3).astype(int)
            for group in np.unique(groups):
                mask = groups == group
                if mask.sum() >= 2:
                    rows.append(
                        [
                            event.planet,
                            event.time,
                            dt[mask].mean(),
                            y[mask].mean(),
                            sigma / np.sqrt(mask.sum()),
                        ]
                    )
    return rows


if __name__ == "__main__":
    with urlopen(BASE + "?" + urlencode({"prefix": PREFIX}), timeout=90) as response:
        root = ET.fromstring(response.read())
    keys = sorted(e.text for e in root.findall(".//{*}Key") if e.text.endswith("_slc.fits"))
    with ThreadPoolExecutor(5) as pool:
        rows = [row for group in pool.map(extract, keys) for row in group]
    pd.DataFrame(rows, columns=["planet", "epoch", "dt", "flux", "error"]).to_csv(
        HERE / "isolated.csv.gz", index=False, float_format="%.9g"
    )
