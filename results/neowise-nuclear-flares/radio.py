"""Radio emission at the flaring nuclei: VLASS quick-look images, NVSS and FIRST.

Usage: python radio.py OUT.json

For each of the seven unreported flares north of declination -40 deg (VLASS
coverage), this asks the CADC archive for every VLASS quick-look Stokes I
image that covers the nucleus, cuts out a 3 arcmin field around it and
measures the brightest pixel within 2.5 arcsec of the nucleus against the
local noise (1.4826 x MAD beyond 10 arcsec). A detection needs 5 sigma;
otherwise 3 sigma is the limit. NVSS (1993-1996) and FIRST (1993-2011)
catalogue matches give the 1.4 GHz state before the flares.
"""
import csv
import io
import json
import math
import sys
import time
from urllib.parse import urlencode

import numpy as np
import requests
from astropy.io import fits
from astropy.wcs import WCS

TAP = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync"
SODA = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/sync"
VIZIER = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"
FLARES = ["02265317+4141211", "06440345-5554310", "11514275-3116072", "21172918-8048085",
          "19265088-1147169", "03183820+3736290", "01110461-4558435"]
QUERY = ("SELECT o.observationID, a.uri, p.time_bounds_lower FROM caom2.Observation o "
         "JOIN caom2.Plane p ON o.obsID = p.obsID JOIN caom2.Artifact a ON p.planeID = a.planeID "
         "WHERE o.collection = 'VLASS' AND p.productID LIKE '%quicklook' "
         "AND a.uri LIKE '%.I.iter1.image.pbcor.tt0.subim.fits' "
         "AND INTERSECTS(POINT('ICRS', {ra}, {dec}), p.position_bounds) = 1")
MPC_CM = 3.0857e24
# Other radio sources worth measuring in the same cutouts: SN 2007D lies 6.3 arcsec from the UGC 2653 nucleus.
EXTRA = {"03183820+3736290": ("SN 2007D", 49.66129167, 37.60733333)}


def get(url, params=None, tries=5, timeout=300):
    for attempt in range(tries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if r.status_code == 200:
                return r
            problem = f"HTTP {r.status_code}: {r.text[:150]!r}"
        except requests.RequestException as error:
            problem = repr(error)
        time.sleep(10 * (attempt + 1))
    raise IOError(problem)


def images(ra, dec):
    r = get(TAP, {"QUERY": QUERY.format(ra=ra, dec=dec), "LANG": "ADQL", "REQUEST": "doQuery", "FORMAT": "csv"})
    return sorted(({"observation": row["observationID"], "uri": row["uri"],
                    "mjd_start": float(row["time_bounds_lower"]) if row["time_bounds_lower"] else None}
                   for row in csv.DictReader(io.StringIO(r.text))), key=lambda d: d["observation"])


def measure(uri, ra, dec):
    r = get(SODA, {"ID": uri, "CIRCLE": f"{ra} {dec} 0.025"})
    with fits.open(io.BytesIO(r.content)) as h:
        hdr, data = h[0].header, np.squeeze(h[0].data).astype(float)
    wcs = WCS(hdr).celestial
    x, y = wcs.all_world2pix([[ra, dec]], 0)[0]
    yy, xx = np.indices(data.shape)
    scale = abs(wcs.pixel_scale_matrix[1, 1]) * 3600.0
    dist = np.hypot(xx - x, yy - y) * scale
    good = np.isfinite(data)
    inner, outer = good & (dist <= 2.5), good & (dist > 10.0)
    if inner.sum() == 0 or outer.sum() < 200:
        return {"status": "no valid pixels"}
    rms = 1.4826 * np.median(np.abs(data[outer] - np.median(data[outer])))
    peak = data[inner].max()
    return {"date_obs": hdr.get("DATE-OBS"), "peak_mjy_beam": 1e3 * peak, "rms_mjy_beam": 1e3 * rms,
            "snr": peak / rms, "detected": bool(peak >= 5 * rms),
            "limit_3sigma_mjy": None if peak >= 5 * rms else 3e3 * rms, "pixel_arcsec": scale}


def catalogue(source, ra, dec, radius):
    r = get(VIZIER, {"-source": source, "-c": f"{ra} {dec}", "-c.rs": radius, "-out.all": "", "-out.max": 10})
    lines = [ln for ln in r.text.splitlines() if ln and not ln.startswith("#")]
    if len(lines) < 4:
        return []
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[3:]]


def main():
    flares = {r["id"]: r for r in csv.DictReader(open("data/flares.csv", encoding="utf-8"))}
    out = []
    for fid in FLARES:
        t = flares[fid]
        ra, dec, dist = float(t["ra"]), float(t["dec"]), float(t["distance_mpc"])
        entry = {"id": fid, "host": t["host"], "ra": ra, "dec": dec, "distance_mpc": dist,
                 "onset_year": float(t["onset_year"])}
        if dec < -40:
            entry["vlass"] = "outside VLASS (declination below -40)"
        else:
            epochs = []
            for img in images(ra, dec):
                m = measure(img["uri"], ra, dec)
                m.update(img)
                if fid in EXTRA:
                    name, xra, xdec = EXTRA[fid]
                    other = measure(img["uri"], xra, xdec)
                    m["other_source"] = {"name": name, **{k: other.get(k) for k in ("peak_mjy_beam", "rms_mjy_beam",
                                                                                     "snr", "detected")}}
                if m.get("limit_3sigma_mjy") or m.get("peak_mjy_beam"):
                    flux = m["peak_mjy_beam"] if m.get("detected") else m.get("limit_3sigma_mjy")
                    m["nu_L_nu_erg_s"] = 3e9 * 4 * math.pi * (dist * MPC_CM) ** 2 * flux * 1e-26
                epochs.append(m)
                print(fid, img["observation"], {k: m.get(k) for k in ("date_obs", "peak_mjy_beam", "rms_mjy_beam",
                                                                       "detected")}, flush=True)
            entry["vlass"] = epochs
        entry["nvss"] = catalogue("VIII/65/nvss", ra, dec, 30)
        entry["first"] = catalogue("VIII/92/first14", ra, dec, 5) if dec > -40 else []
        out.append(entry)
    with open(sys.argv[1], "w", newline="\n") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
