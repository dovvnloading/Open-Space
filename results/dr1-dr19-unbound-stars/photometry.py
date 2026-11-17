"""Sparse photometric variability check at the position of the strongest RV lead."""

import hashlib, json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import numpy as np
from make_queue import rows

if __name__ == "__main__":
    url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?" + urlencode(
        dict(
            POS="CIRCLE 180.19617815 78.79471224 0.000417",
            FORMAT="csv",
            BAD_CATFLAGS_MASK=32768,
        )
    )
    p = Path("ztf.csv")
    if not p.exists():
        p.write_bytes(urlopen(url, timeout=120).read())
    data = rows(p)
    assert data and "catflags" in data[0]
    stats = []
    for oid in sorted({r["oid"] for r in data}):
        group = [r for r in data if r["oid"] == oid and int(r["catflags"]) == 0]
        mag = np.array([float(r["mag"]) for r in group])
        err = [float(r["magerr"]) for r in group]
        mjd = np.array([float(r["mjd"]) for r in group])
        stats.append(
            dict(
                oid=oid,
                filter=group[0]["filtercode"],
                n=len(group),
                distinct_mjd_floors=len(set(np.floor(mjd))),
                span_days=float(np.ptp(mjd)),
                mag_med=float(np.median(mag)),
                mag_q05=float(np.quantile(mag, 0.05)),
                mag_q95=float(np.quantile(mag, 0.95)),
                median_formal_error=float(np.median(err)),
            )
        )
    Path("ztf_summary.json").write_text(
        json.dumps(
            dict(
                url=url, sha256=hashlib.sha256(p.read_bytes()).hexdigest(), groups=stats
            ),
            indent=2,
        )
        + "\n"
    )
