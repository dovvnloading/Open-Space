#!/usr/bin/env python3
"""Independent CODE IONEX + IGRF-14 thin-shell RM check against JPL/ionFR values.

Research diagnostic, NOT ionFR/RMextract bit-for-bit emulation:
- CODE final GNSS IONEX (not the JPL GIM in the released MeerKAT tables)
- IGRF-14 via ppigrf
- 450-km-class IONEX ionospheric shell, spherical ray/geometry
- bilinear spatial and linear UTC time interpolation; Faraday sign for the ray
  propagating *toward* the telescope
- conventional 2.63e-13 * integral n_e[m^-3] B[T] dl[m] -> rad/m².

No external data are bundled. --fetch makes explicit network requests to the
independent CODE archive and writes source URL plus downloaded-file SHA256.
A failed/inaccessible product cannot silently be replaced by a fallback model.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import re
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np

from diagnostics import LAT, LON, PRIOR, solar_altitude

HERE = Path(__file__).resolve().parent
BASE = "https://download.aiub.unibe.ch/CODE"
R_E = 6371.0
C_RM = 2.63e-13


def parse_ionex(data):
    """Return time-stamped maps (lat,lon,TEC) and their shell height.

    IONEX has a 60-column header/tag convention and 5-column integer TEC
    values. Values 9999 denote missing points.
    """
    text = data.decode("ascii", "replace")
    exponent, height = None, None
    lines = text.splitlines()
    for line in lines:
        tag = line[60:].strip() if len(line) >= 60 else ""
        if tag == "EXPONENT":
            exponent = int(line[:60].split()[0])
        elif tag == "HGT1 / HGT2 / DHGT":
            height = float(line[:60].split()[0])
        elif tag == "END OF HEADER":
            break
    if exponent is None:
        raise ValueError("IONEX header lacks EXPONENT")
    maps, current, stamp = [], [], None
    i = 0
    while i < len(lines):
        line = lines[i]
        tag = line[60:].strip() if len(line) >= 60 else ""
        if tag == "START OF TEC MAP":
            current, stamp = [], None
        elif tag == "EPOCH OF CURRENT MAP":
            yy,mm,dd,hh,mi,ss = map(int, line[:60].split()[:6])
            stamp = dt.datetime(yy,mm,dd,hh,mi,ss,tzinfo=dt.timezone.utc)
        elif tag == "LAT/LON1/LON2/DLON/H":
            z = [float(t) for t in re.findall(r"[-+]?\d+(?:\.\d+)?", line[:60])]
            if len(z) != 5:
                raise ValueError("Invalid IONEX latitude row header")
            lat, lon0, lon1, dlon, shell = z
            count = int(round((lon1-lon0)/dlon)) + 1
            vals = []
            while len(vals) < count:
                i += 1
                if i >= len(lines):
                    raise ValueError("Incomplete IONEX latitude TEC data")
                vals.extend(map(float, lines[i].split()))
            vals = np.array(vals[:count]) * 10**exponent
            vals[vals >= 999] = np.nan
            current.append((lat, lon0, dlon, vals))
            height = shell
        elif tag == "END OF TEC MAP":
            if stamp is not None and current:
                maps.append((stamp, current))
        i += 1
    if not maps or height is None:
        raise ValueError("IONEX contains no valid TEC map")
    return maps, float(height)


def map_tec(grid, lat, lon):
    """Bilinear VTEC; neighbor cells must all be valid."""
    rows = sorted(grid, key=lambda x:x[0])
    lats = [v[0] for v in rows]
    i = max(0, min(len(rows)-2, int(np.searchsorted(lats, lat)-1)))
    def sample(row):
        _, start, step, tec = row
        x = ((lon-start) % 360) + start
        if x >= start + step*(len(tec)-1):
            x -= 360
        f = (x-start)/step
        j = max(0, min(len(tec)-2, int(np.floor(f))))
        factor = f-j
        if factor < -1e-6 or factor > 1+1e-6:
            raise ValueError("IONEX longitude outside grid")
        return tec[j]*(1-factor) + tec[j+1]*factor
    top,bot = sample(rows[i]),sample(rows[i+1])
    f = (lat-lats[i])/(lats[i+1]-lats[i])
    return top*(1-f) + bot*f


def time_interp(maps, stamp, lat, lon):
    epochs = np.array([z[0].timestamp() for z in maps])
    moment = stamp.timestamp()
    if not epochs[0] <= moment <= epochs[-1]:
        raise ValueError("Observation outside TEC map coverage")
    i = max(0, min(len(maps)-2, int(np.searchsorted(epochs, moment)-1)))
    a = map_tec(maps[i][1],lat,lon)
    b = map_tec(maps[i+1][1],lat,lon)
    f = (moment-epochs[i])/(epochs[i+1]-epochs[i])
    return (1-f)*a+f*b


def faraday_at(stamp, ra, dec, maps, height):
    """Spherical pierce-point and IGRF thin-layer approximation at the ray."""
    import ppigrf
    mjd = (stamp - dt.datetime(1858,11,17,tzinfo=dt.timezone.utc)).total_seconds()/86400
    alt, solarhour, ha = solar_altitude(mjd, ra, dec)
    if alt < 10:
        raise ValueError("Low elevation outside safe model geometry")
    lat, lon, ha, de = np.deg2rad([LAT, LON, ha*15, dec])
    sin_el = np.sin(np.deg2rad(alt))
    east = -np.cos(de)*np.sin(ha)
    north = np.sin(de)*np.cos(lat) - np.cos(de)*np.sin(lat)*np.cos(ha)
    up = sin_el
    east0 = np.array([-np.sin(lon),np.cos(lon),0])
    north0 = np.array([-np.sin(lat)*np.cos(lon),-np.sin(lat)*np.sin(lon),np.cos(lat)])
    up0 = np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    ray = east*east0+north*north0+up*up0
    t = -R_E*sin_el + np.sqrt((R_E*sin_el)**2 + (R_E+height)**2 - R_E**2)
    pp = R_E*up0 + t*ray
    plat = np.rad2deg(np.arcsin(pp[2]/np.linalg.norm(pp)))
    plon = np.rad2deg(np.arctan2(pp[1],pp[0]))
    vtec = time_interp(maps,stamp,plat,plon)
    b_r,b_theta,b_phi = ppigrf.igrf_gc(
        R_E+height,90-plat,plon,stamp.replace(tzinfo=None))
    br,bn,be = float(np.ravel(b_r)[0]),-float(np.ravel(b_theta)[0]),float(np.ravel(b_phi)[0])
    platr,plonr = np.deg2rad([plat,plon])
    eu = np.array([np.cos(platr)*np.cos(plonr),np.cos(platr)*np.sin(plonr),np.sin(platr)])
    en = np.array([-np.sin(platr)*np.cos(plonr),-np.sin(platr)*np.sin(plonr),np.cos(platr)])
    ee = np.array([-np.sin(plonr),np.cos(plonr),0])
    b_toward_nT = -np.dot(ray, br*eu+bn*en+be*ee)
    slant = 1/np.sqrt(1-(R_E/(R_E+height)*np.cos(np.deg2rad(alt)))**2)
    rm = C_RM * (vtec*1e16*slant) * (b_toward_nT*1e-9)
    return {"code_rm":rm,"vtec":vtec,"slant":slant,"pp_lat":plat,
            "pp_lon":plon,"altitude":alt,"solar_hour":solarhour}


def archive_filename(date):
    year, doy = date.year, date.timetuple().tm_yday
    # AIUB finalized long naming starts 2022 day 331; earlier files are .Z.
    if year > 2022 or (year == 2022 and doy >= 331):
        return f"COD0OPSFIN_{year}{doy:03d}0000_01D_01H_GIM.INX.gz"
    return f"CODG{doy:03d}0.{year%100:02d}I.Z"


def load_day(date, cache, fetch):
    filename = archive_filename(date)
    path = cache / filename
    url = f"{BASE}/{date.year}/{filename}"
    if not path.exists():
        if not fetch:
            raise FileNotFoundError(f"Missing archived product {path}; pass --fetch")
        req = urllib.request.Request(url,headers={"User-Agent":"MeerKAT-CODE-IONEX-research/1.0"})
        with urllib.request.urlopen(req,timeout=45) as r:
            payload = r.read()
        if len(payload) < 30_000:
            raise ValueError("Unexpectedly short IONEX download")
        cache.mkdir(parents=True,exist_ok=True)
        path.write_bytes(payload)
    payload=path.read_bytes()
    if filename.endswith(".gz"):
        decoded=gzip.decompress(payload)
    else:
        import unlzw3
        decoded=unlzw3.unlzw(payload)
    maps, height = parse_ionex(decoded)
    return (maps,height),{"file":filename,"url":url,
                          "sha256":hashlib.sha256(payload).hexdigest(),
                          "n_maps":len(maps),"shell_height_km":height}


def read_data():
    stars = defaultdict(list)
    with (PRIOR/"measurements.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            z = {k:float(v) if k!="psr" else v for k,v in row.items()}
            z["e"] = float(np.hypot(z["rm_err"],z["iono_err"]))
            z["y"] = z["rm"]-z["iono"]
            stars[z["psr"]].append(z)
    baselines = {n:float(np.median([z["y"] for z in v if z["e"]<2]))
                 for n,v in stars.items()
                 if len(v)>=18 and sum(z["e"]<2 for z in v)>=12}
    return stars,baselines


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--dates",default="2022-10-08,2022-11-12,2022-12-16,2023-01-31",
                    help="Comma-separated UTC dates; preserve date selection in receipts")
    ap.add_argument("--cache",type=Path,default=HERE/".ionex_cache")
    ap.add_argument("--fetch",action="store_true")
    ap.add_argument("--write",action="store_true")
    args=ap.parse_args()
    dates=[dt.date.fromisoformat(s.strip()) for s in args.dates.split(",")]
    stars,baselines=read_data()
    allrows=[z for n,v in stars.items() if n in baselines for z in v
             if z["e"]<1.8 and dt.date.fromordinal(
                 dt.date(1858,11,17).toordinal()+int(z["mjd"])) in dates]
    summaries=[];rows=[];receipts=[]
    for day in dates:
        (maps,height),receipt=load_day(day,args.cache,args.fetch)
        receipts.append(receipt)
        start=len(rows)
        for z in allrows:
            epoch=dt.datetime(1858,11,17,tzinfo=dt.timezone.utc)+dt.timedelta(days=z["mjd"])
            if epoch.date()!=day:
                continue
            try:
                calc=faraday_at(epoch,z["ra_deg"],z["dec_deg"],maps,height)
            except ValueError:
                continue
            if not np.isfinite(calc["code_rm"]):
                continue
            rows.append({"date":day.isoformat(),"psr":z["psr"],"mjd":z["mjd"],
                         "ra_deg":z["ra_deg"],"dec_deg":z["dec_deg"],
                         "ionfr_jpl_rm":z["iono"],"ionex_code_rm":calc["code_rm"],
                         "code_minus_jpl":calc["code_rm"]-z["iono"],
                         "jpl_corrected_departure":z["y"]-baselines[z["psr"]],
                         **calc})
        these=rows[start:]
        if not these:
            raise ValueError(f"No finite CODE independent RM for {day}")
        def med(k):return float(np.median([z[k] for z in these]))
        corr=float(np.corrcoef([z["ionfr_jpl_rm"] for z in these],
                                [z["ionex_code_rm"] for z in these])[0,1])
        summaries.append({"date":day.isoformat(),"n_samples":len(these),
                          "n_pulsars":len({z["psr"] for z in these}),
                          "median_code_rm":med("code_rm"),
                          "median_jpl_rm":med("ionfr_jpl_rm"),
                          "median_code_minus_jpl":med("code_minus_jpl"),
                          "median_baseline_departure_jpl":med("jpl_corrected_departure"),
                          "jpl_code_across_sources_r":corr})
    result={"method":"CODE GNSS IONEX + IGRF-14, spherical thin shell, not ionFR",
            "dates_preselected": [d.isoformat() for d in dates],
            "sample_counts":len(rows),
            "daily":summaries,"receipts":receipts,
            "caution":"Pilot selected dates; no independent full-series baseline or validated RMextract reproduction; correlation does not identify an origin."}
    if args.write:
        (HERE/"code_ionex_probe_checks.json").write_text(json.dumps(result,indent=2)+"\n")
        with (HERE/"code_ionex_probe_rows.csv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]))
            w.writeheader();w.writerows(rows)
    print(json.dumps({"n":len(rows),"daily":summaries}))


if __name__=="__main__":
    main()
