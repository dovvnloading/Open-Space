"""Follow up the product dependence for TIC 257433294, including earlier visits and pixels."""

import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from measure import ROOT, compile_result, measure_events, write_rows
from fetch import get

TIC = 257433294
SECTORS = {
    4: "2018292075959-s0004-0000000257433294-0124",
    105: "2026164183000-s0105-0000000257433294-0307",
    106: "2026192185000-s0106-0000000257433294-0308",
    107: "2026221103959-s0107-0000000257433294-0309",
    108: "2026250081459-s0108-0000000257433294-0310",
}


def inputs():
    manifest = []
    for sector, stem in SECTORS.items():
        for kind in (["lc", "tp"] if sector == 108 else ["lc"]):
            name = f"tess{stem}-s_{kind}.fits"
            url = f"https://stpubdata.s3.amazonaws.com/tess/public/tid/s{sector:04}/0000/0002/5743/3294/{name}"
            path = ROOT / "fits" / name
            manifest.append(
                {
                    "sector": sector,
                    "kind": kind,
                    "name": name,
                    "url": url,
                    "sha256": get(url, path),
                }
            )
    return manifest


def measure(t, flux, mask, cadence, c):
    m = mask & np.isfinite(t) & np.isfinite(flux) & (flux > 0)
    e = float(c["T0_BJD_TDB"]) - 2457000
    period = float(c["period_d"])
    duration = float(c["duration_h"]) / 24
    result = measure_events(t[m], flux[m], e, period, duration, cadence, 2.0)
    return result, compile_result(result, (t[m].min() + t[m].max()) / 2)


def main():
    manifest = inputs()
    (ROOT / "followup_provenance.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    with (ROOT / "ephemerides.csv").open() as inp:
        c = next(r for r in csv.DictReader(inp) if int(r["TIC"]) == TIC)
    all_summaries = []
    all_events = []
    for item in manifest:
        with fits.open(ROOT / "fits" / item["name"]) as h:
            tab = h[1].data
            t = np.array(tab["TIME"])
            cadence = h[1].header["TIMEDEL"]
            mask = tab["QUALITY"] == 0
            if item["kind"] == "lc":
                raw = np.asarray(tab["SAP_FLUX"], dtype=float)
                pdc = np.asarray(tab["PDCSAP_FLUX"], dtype=float)
                common = mask & np.isfinite(raw) & np.isfinite(pdc)
                trials = [
                    ("SAP", raw, mask),
                    ("PDCSAP", pdc, mask),
                    ("SAP_common", raw, common),
                    ("PDCSAP_common", pdc, common),
                ]
            else:
                aper = h[2].data
                y, x = np.indices(aper.shape)
                tx, ty = WCS(h[2].header, naxis=2).world_to_pixel_values(
                    h[0].header["RA_OBJ"], h[0].header["DEC_OBJ"]
                )
                dist = np.hypot(x - tx, y - ty)
                cube = tab["FLUX"]
                trials = []
                for name, region in [
                    ("SPOC_mask", (aper & 2) != 0),
                    ("one_pixel", dist == dist.min()),
                    (
                        "core_3x3",
                        (abs(x - round(float(tx))) <= 1)
                        & (abs(y - round(float(ty))) <= 1),
                    ),
                    ("annulus_3_to_5", (dist >= 3) & (dist <= 5)),
                ]:
                    # Sum background-subtracted calibrated pixels. Keep missing individual pixels missing.
                    flux = np.sum(cube[:, region], axis=1, dtype=np.float64)
                    trials.append((name, flux, mask))
            for name, flux, selected in trials:
                events, summary = measure(t, flux, selected, cadence, c)
                prefix = {"tic": TIC, "sector": item["sector"], "extraction": name}
                all_summaries.append({**prefix, **summary})
                all_events.extend(
                    {
                        **prefix,
                        "cycle": n,
                        "centre_btjd": round(centre, 7),
                        "depth_ppm": d,
                    }
                    for n, centre, d in events
                )
    write_rows("followup_summary.csv", all_summaries)
    write_rows("followup_events.csv", all_events)


if __name__ == "__main__":
    main()
