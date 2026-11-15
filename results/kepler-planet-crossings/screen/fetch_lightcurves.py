from pathlib import Path
import concurrent.futures, urllib.request as net, urllib.parse as parse, xml.etree.ElementTree as ET
import pandas as pd, time
from astropy.io import fits
from io import BytesIO
import numpy as np

threads = 14
url = "https://stpubdata.s3.amazonaws.com/"
Path("cache").mkdir(exist_ok=True)


def read(u):
    for n in range(3):
        try:
            return net.urlopen(u, timeout=16).read()
        except Exception:
            if n == 2:
                raise
            time.sleep(1 + n)


def listing(k):
    kd = f"{int(k):09d}"
    prefix = f"kepler/public/lightcurves/{kd[:4]}/{kd}/"
    root = ET.fromstring(read(url + "?" + parse.urlencode(dict(prefix=prefix))))
    return [n.text for n in root.findall(".//{*}Key") if n.text.endswith("_llc.fits")]


def get(k):
    p = Path("cache") / (k.split("/")[-1].replace(".fits", ".npz"))
    if p.exists():
        return
    try:
        with fits.open(BytesIO(read(url + k))) as h:
            d = h[1].data
            cols = [
                "TIME",
                "PDCSAP_FLUX",
                "PDCSAP_FLUX_ERR",
                "SAP_FLUX",
                "SAP_QUALITY",
                "MOM_CENTR1",
                "MOM_CENTR2",
            ]
            np.savez_compressed(p, **{c: np.copy(d[c]) for c in cols})
    except Exception as e:
        raise RuntimeError(f"Failed to fetch {k}") from e


if __name__ == "__main__":
    import sys

    pairs = pd.read_csv("pairs.csv")
    stars = pairs[pairs.snr > 6].k.drop_duplicates()
    with concurrent.futures.ThreadPoolExecutor(threads) as ex:
        keys = [k for ks in ex.map(listing, stars) for k in ks]
        print(len(stars), len(keys), flush=True)
        list(ex.map(get, keys))
