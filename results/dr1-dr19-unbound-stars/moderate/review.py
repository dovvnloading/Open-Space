"""Follow the additional 550–900 km/s SpecDis screen in this directory."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from screen_desi import cartesian
from make_queue import rows, save, number


def build():
    b = {x["source_id"]: x for x in rows("bj.csv")}
    s = {x["SOURCE_ID"]: x for x in rows("desi_screen.csv")}
    checks = rows("desi_checks.csv")
    q = []
    for c in checks:
        sid = c["source_id"]
        a = s[sid]
        j = b.get(sid, {})
        d = number(j.get("r_med_photogeo")) / 1000
        v = (
            np.linalg.norm(
                cartesian(
                    *[number(a[k]) for k in ("RA", "DEC")],
                    d,
                    *[number(a[k]) for k in ("PMRA", "PMDEC", "VRAD")],
                )[1]
            )
            if d > 0
            else np.nan
        )
        c["bj_photogeo_speed"] = v
        if c["rv_ok"] != "1" or not v >= 550:
            continue
        q.append(
            dict(
                source_id=sid,
                targetid=c["targetid"],
                channel="moderate",
                rv=c["rv"],
                e_rv=c["rv_err"],
                spec_d=c["spec_d"],
                spec_e=c["spec_sigma"],
                sp_d=c["spdist_d"],
                sp_e=c["spdist_sigma"],
                zp=a["PARALLAX_ZPC"],
                **{k: v for k, v in j.items() if k != "source_id"},
            )
        )
    from astropy.io import fits
    import hashlib

    with open(sys.argv[1], "rb") as f:
        assert (
            hashlib.file_digest(f, "sha256").hexdigest()
            == "eb4da91267db39f285a277989489f991ea48371336efedd28bd07d2b58e4a400"
        )
    backup = {
        r["targetid"]
        for r in rows("desi/queries/coadds.csv")
        if r["program"] == "backup"
    }
    with fits.open(sys.argv[1]) as f:
        a = f[1].data
        lookup = {
            str(r["TARGETID"]): r["VRAD_BIAS"]
            for r in a[np.isin(a["TARGETID"], [int(s["targetid"]) for s in q])]
        }
    for s in q:
        if s["targetid"] in backup and s["targetid"] in lookup:
            s["rv_raw"] = s["rv"]
            s["rv_bias"] = lookup[s["targetid"]]
            s["rv"] = number(s["rv"]) - s["rv_bias"]
            s["rv_floor"] = 2.0
    save("queue.csv", q)
    save("screen_review.csv", checks)


if __name__ == "__main__":
    build()
