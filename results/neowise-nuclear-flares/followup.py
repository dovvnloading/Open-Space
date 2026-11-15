"""Independent checks on one flare: ZTF nuclear photometry and WISE exposures.

Usage: python followup.py BATCHDIR ID RA DEC OUT.json

ZTF season medians come from the public IRSA light-curve service (PSF
photometry within 2.2 arcsec, unflagged epochs, objects with 40 or more
points). The WISE part lists, for each half-year, the number of exposures,
median magnitudes, the reduced chi-square of the point-source fit and the
artifact flags.
"""
import io
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import requests

import epochs

ZTF = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"


def ztf_seasons(ra, dec):
    out = {}
    for band in ("g", "r"):
        response = requests.get(ZTF, params={"POS": f"CIRCLE {ra} {dec} 0.0006", "BANDNAME": band,
                                             "FORMAT": "csv", "BAD_CATFLAGS_MASK": 32768}, timeout=300)
        table = pd.read_csv(io.StringIO(response.text))
        if len(table) == 0:
            continue
        table = table[table.catflags == 0]
        for oid, group in table.groupby("oid"):
            if len(group) < 40:
                continue
            season = np.floor(2000 + (group.mjd - 51544.5) / 365.25 + 0.4).astype(int)
            rows = []
            for year, part in group.groupby(season):
                mad = 1.4826 * np.median(np.abs(part.mag - part.mag.median()))
                rows.append({"season": int(year), "n": int(len(part)), "median_mag": round(float(part.mag.median()), 3),
                             "error": round(float(mad / np.sqrt(len(part))), 3)})
            out[f"{band}:{oid}"] = rows
    return out


def wise_exposures(batchdir, name):
    for path in sorted(pathlib.Path(batchdir).glob("batch_*.ids")):
        ids = path.read_text().splitlines()
        if name in ids:
            table = epochs.read_batch(path.with_name(path.name.replace(".ids", ".tbl.gz")))
            g = table[(table.cntr == ids.index(name) + 1) & (table.dist <= 3)].copy()
            g["half_year"] = (2000 + (g.mjd - 51544.5) / 365.25).round(1)
            rows = [{"half_year": float(k), "n": int(len(v)), "w1": round(float(v.w1.median()), 3),
                     "w2": round(float(v.w2.median()), 3), "chi2_w1": round(float(v.chi1.median()), 1),
                     "chi2_w2": round(float(v.chi2.median()), 1)} for k, v in g.groupby("half_year")]
            return {"exposures": int(len(g)), "cc_flags": g.cc_flags.value_counts().to_dict(), "by_half_year": rows}
    return None


def main():
    batchdir, name, ra, dec, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
    result = {"id": name, "ra": ra, "dec": dec, "ztf_season_medians": ztf_seasons(ra, dec) if dec > -28 else {},
              "wise_exposures": wise_exposures(batchdir, name)}
    json.dump(result, open(out, "w"), indent=1)
    print(name, "ZTF objects:", list(result["ztf_season_medians"]), "WISE exposures:", result["wise_exposures"]["exposures"])


if __name__ == "__main__":
    main()
