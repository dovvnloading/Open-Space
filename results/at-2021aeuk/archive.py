"""Verify the pinned Swift catalogue and photon measurements against the public archives.

Requires NumPy and Astropy. Hashes cover the archived compressed FITS bytes.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parent
POSITION = SkyCoord(243.24934167, 42.32784722, unit="deg")


def check_lsxps():
    request = dict(APIFunc="getSXPSLC", APIVersion="1.0.5", whichCat="LSXPS",
                   sourceID=1248139, bands="all", binning="Observation", timeFormat="MJD",
                   showObs=True, getAllTypes=True)
    payload = Request("https://www.swift.ac.uk/API/main.php", data=json.dumps(request).encode(),
                      headers={"Content-Type": "application/json"})
    with urlopen(payload, timeout=90) as response:
        data = json.load(response)
    if not data.get("OK"):
        raise ValueError(data)
    live = {}
    for band in ("Total", "Soft", "Medium", "Hard"):
        for kind in ("blind", "retro", "nondet"):
            block = data[f"{band}_{kind}_ratesData"]
            for values in block["data"]:
                row = dict(zip(block["columns"], values))
                live[band, kind, row["ObsID"]] = row
    with (ROOT / "lsxps.csv").open() as stream:
        pinned = list(csv.DictReader(stream))
    for row in pinned:
        source = live[row["band"], row["status"], row["obsid"]]
        for field, column in (("mjd","Time"), ("source_counts","CtsInSrc"),
                              ("expected_background","BGInSrc"), ("correction","CorrFact"),
                              ("exposure_s","Exposure")):
            np.testing.assert_allclose(float(row[field]), source[column], rtol=0, atol=5.1e-9)
    print(f"Verified {len(pinned)} pinned band/observation rows. Current catalogue has {len(live)}.")


def xrt_measure(hdu):
    hdr = hdu.header
    # EVENTS sky columns 2 and 3 use a TAN tangent plane in degrees.
    wcs = WCS(naxis=2)
    for attr, prefix in (("crpix","TCRPX"), ("crval","TCRVL"), ("cdelt","TCDLT"), ("ctype","TCTYP")):
        setattr(wcs.wcs, attr, [hdr[prefix+"2"], hdr[prefix+"3"]])
    x,y = wcs.all_world2pix([[POSITION.ra.deg, POSITION.dec.deg]], 1)[0]
    t = hdu.data
    good = (t["PI"] >= 30) & (t["PI"] <= 1000) & (t["GRADE"] <= 12)
    d = np.hypot(t["X"][good] - x, t["Y"][good] - y)
    return [dict(exp=hdr["EXPOSURE"], n10=int(sum(d<=10)), n15=int(sum(d<=15)),
                 n20=int(sum(d<=20)), annulus=int(sum((d>30)&(d<=60))))]


def exposure_measure(hdu):
    x,y = WCS(hdu.header, naxis=2).all_world2pix([[POSITION.ra.deg, POSITION.dec.deg]], 0)[0]
    data = hdu.data
    yy,xx = np.ogrid[:data.shape[0],:data.shape[1]]
    d2 = (xx-x)**2 + (yy-y)**2
    return [dict(center=float(data[round(y),round(x)]), median10=float(np.median(data[d2<100])),
                 mean10=float(np.mean(data[d2<100])), medianann=float(np.median(data[(d2>900)&(d2<=3600)])))]


def uv_measure(hdul):
    sources = []
    controls = []
    for hdu in hdul:
        if not isinstance(hdu, fits.BinTableHDU) or "FILTER" not in hdu.header or not len(hdu.data):
            continue
        band = hdu.header["FILTER"]
        t = hdu.data
        d = separation(t["RA"], t["DEC"], POSITION.ra.deg, POSITION.dec.deg)
        i = np.argmin(d)
        row = t[i]
        sources.append(dict(band=band,sep=round(float(d[i]),4),mag=float(row["MAG"]),
                            magerr=float(row["MAG_ERR"]),flags=int(row["FLAGS"]),
                            rate=float(row["SENSCORR_RATE"]),rateerr=float(row["SENSCORR_RATE_ERR"]),
                            exp=hdu.header["EXPOSURE"]))
        if band == "UVW1":
            for row in t[(t["FLAGS"] == 0) & (t["MAG_ERR"] < .15) & np.isfinite(t["MAG"])]:
                controls.append([round(float(row[key]),7) for key in ("RA","DEC","MAG","MAG_ERR")])
    return sources, controls


def separation(ra, dec, ra0, dec0):
    return 3600*np.hypot((ra-ra0)*np.cos(np.deg2rad(dec0)), dec-dec0)


def check_product(entry):
    with urlopen(entry["url"], timeout=90) as response:
        blob = response.read()
    if hashlib.sha256(blob).hexdigest() != entry["sha256"]:
        raise ValueError(f"Product revision changed at {entry['url']}")
    controls = []
    with fits.open(io.BytesIO(gzip.decompress(blob)), memmap=False) as hdul:
        if entry["mode"] == "XRT":
            measured = xrt_measure(hdul["EVENTS"])
        elif entry["mode"] == "EXPO":
            image = next(h for h in hdul if h.header.get("NAXIS") == 2 and h.data is not None)
            measured = exposure_measure(image)
        else:
            measured, controls = uv_measure(hdul)
    if len(measured) != len(entry["rows"]):
        raise ValueError(f"Extension count differs in {entry['url']}")
    for actual,expected in zip(measured, entry["rows"]):
        for key,value in expected.items():
            if isinstance(value, str):
                assert actual[key] == value
            else:
                np.testing.assert_allclose(actual[key], value, rtol=2e-6, atol=1e-6)
    return entry["obs"], controls


def select_controls(catalogues, early_ids):
    tracks = []
    for ra, dec, mag, error in catalogues["00049502090"]:
        if not (20 < separation(ra,dec,POSITION.ra.deg,POSITION.dec.deg) < 420 and 16 < mag < 19):
            continue
        visits = []
        for obs, catalogue in sorted(catalogues.items()):
            matches = [r for r in catalogue if separation(r[0],r[1],ra,dec) <= 1.5]
            if len(matches) == 1:
                visits.append([obs,matches[0][2]])
        n_early = sum(obs in early_ids for obs,mag in visits)
        if min(n_early,len(visits)-n_early) >= 6:
            tracks.append(dict(ra=ra,dec=dec,visits=visits))
    return tracks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lsxps-only", action="store_true")
    args = parser.parse_args()
    check_lsxps()
    if not args.lsxps_only:
        packet = json.loads((ROOT / "swift-products.json").read_text())
        with ThreadPoolExecutor(max_workers=4) as pool:
            catalogues = {obs:field for obs,field in pool.map(check_product,packet["products"]) if field}
        with (ROOT / "lsxps.csv").open() as stream:
            early_ids = {r["obsid"] for r in csv.DictReader(stream) if 60800 <= float(r["mjd"]) < 61000}
        assert select_controls(catalogues,early_ids) == packet["controls"]
        print(f"Verified {len(packet['products'])} FITS products and the full UVW1 field-control selection.")
