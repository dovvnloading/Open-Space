"""Test fixed ephemerides on TESSCut FFIs, with local background subtraction."""

import argparse
import csv
import hashlib
import io
import json
import warnings
import zipfile
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from measure import ROOT, compile_result, measure_events, write_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sector", choices=[108, 109], type=int, default=109)
    sector = parser.parse_args().sector
    source = json.loads((ROOT / "provenance.json").read_text())
    ids = {r["tic"] for r in source["selected"]}
    with (ROOT / "ephemerides.csv").open() as f:
        candidates = [r for r in csv.DictReader(f) if int(r["TIC"]) in ids]
    rec = []
    sums = []
    tests = []
    epochs = []
    for c in candidates:
        tic = int(c["TIC"])
        ra = float(c["ra"])
        dec = float(c["dec"])
        params = {"ra": ra, "dec": dec, "x": 9, "y": 9, "units": "px", "sector": sector}
        url = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?" + urlencode(params)
        path = ROOT / "fits" / f"tesscut{sector}_{tic}.zip"
        if not path.exists():
            path.write_bytes(urlopen(url, timeout=120).read())
        b = path.read_bytes()
        info = {"tic": tic, "query_url": url}
        if not zipfile.is_zipfile(io.BytesIO(b)):
            response = json.loads(b)
            if response.get(
                "results"
            ) != [] or "Missing sector information" not in response.get("msg", ""):
                raise ValueError(response)
            rec.append({**info, "status": response["msg"]})
            continue
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            names = [n for n in z.namelist() if n.endswith(".fits")]
            assert len(names) == 1
            blob = z.read(names[0])
            info["filename"] = names[0]
        with fits.open(io.BytesIO(blob)) as h:
            assert (
                h[0].header["SECTOR"] == sector
                and h[1].header["TIMEREF"] == "SOLARSYSTEM"
            )
            tab = h[1].data
            t = tab["TIME"].astype(float)
            cube = tab["FLUX"].astype(float)
            flags = tab["QUALITY"]
            d = h[1].header["TIMEDEL"]
            yy, xx = np.indices(cube.shape[1:])
            tx, ty = WCS(h[2].header, naxis=2).world_to_pixel_values(ra, dec)
            dist = np.hypot(xx - tx, yy - ty)
            # TESSCut may regenerate packaging metadata. Hash the numerical inputs instead.
            digest = hashlib.sha256()
            for arr, dtype in [
                (t, "<f8"),
                (cube, "<f8"),
                (flags, "<i8"),
                (dist, "<f8"),
            ]:
                digest.update(np.ascontiguousarray(arr, dtype=dtype).tobytes())
            rec.append(
                {
                    **info,
                    "status": "covered",
                    "numeric_sha256": digest.hexdigest(),
                    "fits_sha256": hashlib.sha256(blob).hexdigest(),
                    "camera": h[0].header["CAMERA"],
                    "ccd": h[0].header["CCD"],
                    "cadences": len(t),
                }
            )
            edge = (abs(xx - 4) == 4) | (abs(yy - 4) == 4)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                median = np.nanmedian(cube, axis=0)
                chosen = edge & (median <= np.median(median[edge]))
                sky = np.nanmedian(cube[:, chosen], axis=1)
            for label, aperture in [
                (
                    "3x3",
                    (abs(xx - round(float(tx))) <= 1)
                    & (abs(yy - round(float(ty))) <= 1),
                ),
                ("2x2", dist <= np.sort(dist.ravel())[3]),
            ]:
                flux = np.sum(cube[:, aperture], axis=1) - aperture.sum() * sky
                m = (flags == 0) & np.isfinite(t) & np.isfinite(flux) & (flux > 0)
                tm = t[m]
                f = flux[m]
                e = float(c["T0_BJD_TDB"]) - 2457000
                p = float(c["period_d"])
                dur = float(c["duration_h"]) / 24
                ev = measure_events(tm, f, e, p, dur, d, 2.0)
                s = compile_result(ev, (tm[0] + tm[-1]) / 2)
                depths = []
                for k in range(1, 200):
                    off = k / 200
                    if min(off, 1 - off) * p < 2.5 * dur:
                        continue
                    trial = measure_events(tm, f, e + off * p, p, dur, d, 2.0)
                    ctl = compile_result(trial, (tm[0] + tm[-1]) / 2)
                    if ctl["covered"] == 0 or abs(ctl["covered"] - s["covered"]) > max(
                        1, round(0.2 * s["covered"])
                    ):
                        continue
                    depths.append(ctl["depth_ppm"])
                    tests.append(
                        {
                            "tic": tic,
                            "aperture": label,
                            "phase": off,
                            "depth_ppm": ctl["depth_ppm"],
                            "covered": ctl["covered"],
                        }
                    )
                sums.append(
                    {
                        "tic": tic,
                        "aperture": label,
                        **s,
                        "offphase_count": len(depths),
                        "offphase_equal_deeper": (
                            sum(v >= s["depth_ppm"] for v in depths)
                            if s["depth_ppm"] is not None
                            else None
                        ),
                    }
                )
                epochs.extend(
                    {
                        "tic": tic,
                        "aperture": label,
                        "cycle": n,
                        "centre_btjd": round(t0, 7),
                        "depth_ppm": depth,
                    }
                    for n, t0, depth in ev
                )
    (ROOT / f"ffi{sector}_provenance.json").write_text(json.dumps(rec, indent=2) + "\n")
    write_rows(f"ffi{sector}_results.csv", sums)
    write_rows(f"ffi{sector}_events.csv", epochs)
    write_rows(f"ffi{sector}_phase_controls.csv", tests)


if __name__ == "__main__":
    main()
