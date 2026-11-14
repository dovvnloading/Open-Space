"""Extract the two MAXI photon apertures and convert event times to timing bins."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
from astropy import units as u
from astropy.coordinates import get_body_barycentric
from astropy.io import fits
from astropy.time import Time
from scipy.interpolate import CubicSpline

ROOT = Path(__file__).resolve().parent
TARGET_RA, TARGET_DEC = np.deg2rad([235.597083, -52.386111])
START, STOP, SECONDS = 61132, 61325, 16


def light_time_grid():
    direction = np.array([np.cos(TARGET_DEC)*np.cos(TARGET_RA),
                          np.cos(TARGET_DEC)*np.sin(TARGET_RA), np.sin(TARGET_DEC)])
    grid = np.arange(START-5.,STOP+5.)
    earth = get_body_barycentric("earth",Time(grid,format="mjd",scale="tt"))
    return CubicSpline(grid,earth.xyz.to_value(u.m).T@direction/299792458.0)


def emission_times(local_mjd, light_time):
    bary = local_mjd + light_time(local_mjd)/86400
    eccentricity = 0.174
    omega = np.deg2rad(64.)
    anomaly = 2*np.pi*(bary-51441.989)/3.7283749 + np.pi/2-omega
    ecc = anomaly.copy()
    for _ in range(5):
        ecc -= (ecc-eccentricity*np.sin(ecc)-anomaly)/(1-eccentricity*np.cos(ecc))
    delay = 56.6*((np.cos(ecc)-eccentricity)*np.sin(omega)
                  + np.sqrt(1-eccentricity**2)*np.sin(ecc)*np.cos(omega))
    return bary-delay/86400


def process(cache, download, output):
    refs = json.loads((ROOT/"sources.json").read_text())
    light_time = light_time_grid()
    count = int((STOP-START)*86400/SECONDS)
    source, annulus = np.zeros(count,dtype=int), np.zeros(count,dtype=int)
    rows = []
    for item in refs["maxi"]:
        path = cache/item["url"].split("/")[-1]
        if download and not path.exists():
            cache.mkdir(parents=True,exist_ok=True)
            urlretrieve(item["url"],path)
        if hashlib.sha256(path.read_bytes()).hexdigest()!=item["sha256"]:
            raise ValueError(f"File differs from original archive bytes: {path}")
        with fits.open(path) as hdus:
            header, table = hdus[1].header, hdus[1].data
            ra, dec = np.deg2rad(table["RA"].astype(float)), np.deg2rad(table["DEC"].astype(float))
            cos_dist = np.sin(TARGET_DEC)*np.sin(dec)+np.cos(TARGET_DEC)*np.cos(dec)*np.cos(ra-TARGET_RA)
            radius = np.rad2deg(np.arccos(np.clip(cos_dist,-1,1)))
            good = ((table["PI"]>=60)&(table["PI"]<300)&(table["QUAL_POS"]==0)
                    &(table["QUAL_PI"]==0)&(table["DFLAG"]==0))
            reference = header["MJDREFI"]+header["MJDREFF"]
            for result,region in [(source,radius<1.5),(annulus,(radius>2.5)&(radius<4.0))]:
                local = reference+table["TIME"][good & region]/86400
                ts = emission_times(local,light_time)
                bins = np.floor((ts-START)*86400/SECONDS).astype(int)
                bins = bins[(bins>=0)&(bins<count)]
                np.add.at(result,bins,1)
    for idx in np.flatnonzero(source+annulus):
        rows.append((idx,source[idx],annulus[idx]))
    # mtime=0 also makes the compressed representation deterministic.
    import io
    with output.open("wb") as raw, gzip.GzipFile(fileobj=raw,filename="",mode="wb",mtime=0) as zipped:
        with io.TextIOWrapper(zipped,encoding="ascii",newline="") as stream:
            writer = csv.writer(stream,lineterminator="\n")
            writer.writerow(["bin16_since_mjd61132","source","annulus"])
            writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--download",action="store_true")
    args = parser.parse_args()
    process(args.cache,args.download,args.output)
