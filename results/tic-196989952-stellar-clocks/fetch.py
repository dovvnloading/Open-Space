"""Download the exact mission and survey products used by analyze.py."""
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode, quote
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.io import fits

HERE = Path(__file__).resolve().parent

def retrieve(url):
    with urlopen(url, timeout=100) as response:
        return response.read()

def main():
    inputs = []
    def fetch(url):
        data = retrieve(url)
        inputs.append({"url": url, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
        return data

    rows = []
    for sector in [20, 60, 73]:
        key = f"s{sector:04d}"
        product = f"mast:HLSP/qlp/{key}/0000/0001/9698/9952/hlsp_qlp_tess_ffi_{key}-0000000196989952_tess_v01_llc.fits"
        data = fetch("https://mast.stsci.edu/api/v0.1/Download/file?uri=" + product)
        with fits.open(io.BytesIO(data)) as hdus:
            events = hdus[1].data
            good = (events["QUALITY"] == 0) & np.isfinite(events["TIME"] + events["SAP_FLUX"])
            rows.append(pd.DataFrame({"sector": sector, "btjd": events["TIME"][good].astype(float),
                                      "flux": events["SAP_FLUX"][good].astype(float)}))
    pd.concat(rows, ignore_index=True).to_csv(HERE / "tess.csv.gz", index=False)

    address = "https://wasp.cerit-sc.cz/csv?object=" + quote("1SWASP J062202.73+402104.6", safe="")
    wasp = pd.read_csv(io.BytesIO(fetch(address)))
    wasp.to_csv(HERE / "wasp.csv.gz", index=False)

    query = {"-source": "I/355/epphot", "Source": "956785024104626432", "-out.max": 9000, "-out.all": ""}
    raw = fetch("https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?" + urlencode(query))
    lines = [line for line in raw.decode().splitlines() if line and not line.startswith("#")]
    gaia = pd.read_csv(io.StringIO("\n".join([lines[0]] + lines[3:])), sep="\t")
    gaia[["Source", "TransitID", "TimeG", "FG", "e_FG", "Gmag", "GrVFlag", "noisyFlag"]].to_csv(HERE / "gaia.csv", index=False)
    (HERE / "inputs.json").write_text(json.dumps(inputs, indent=2) + "\n")

if __name__ == "__main__":
    main()
