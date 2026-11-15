"""Recalculate barycentric times of the two Figure 3 flashes with JPL DE430.

Run after extract.py with the same input directory. See README for the
independent all-catalog time-label check. Requires Astropy and jplephem.
"""
import argparse
import csv
import math
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric, get_body_barycentric_posvel, solar_system_ephemeris
from astropy.io import fits
from astropy.time import Time, TimeDelta

from extract import HERE, inputs


def barycentric(event_time, header, orbit_file, ra_deg, dec_deg):
    """Barycentric time using Earth and spacecraft Roemer, Einstein, solar Shapiro terms."""
    mission_time = event_time + header.get("TIMEZERO",0)
    tt_epoch = Time(header["MJDREFI"],header["MJDREFF"],format="mjd",scale="tt")
    t = tt_epoch + TimeDelta(mission_time,format="sec")
    with fits.open(orbit_file) as hdus:
        orbit = hdus["ORBIT"]
        track = orbit.data["TIME"] + orbit.header.get("TIMEZERO",0)
        if not track[0] <= mission_time <= track[-1]:
            raise ValueError("No spacecraft orbit coverage")
        sat = np.array([np.interp(mission_time,track,orbit.data[k]) for k in "XYZ"])
    ra, dec = map(math.radians,(ra_deg,dec_deg))
    source = [math.cos(ra)*math.cos(dec),math.sin(ra)*math.cos(dec),math.sin(dec)]
    earth, speed = get_body_barycentric_posvel("earth",t)
    obs = earth.xyz.to_value("m") + sat
    from_sun = obs - get_body_barycentric("sun",t).xyz.to_value("m")
    c = 299792458
    roemer = np.dot(obs,source)/c
    einstein = np.dot(speed.xyz.to_value("m/s"),sat)/c**2
    shapiro = 9.8509819e-6*math.log1p(np.dot(from_sun,source)/np.linalg.norm(from_sun))
    return t, t.tdb + TimeDelta(roemer+einstein+shapiro,format="sec")


def main(cache):
    labels = {"1020240141":"2018-11-22T11:58:14.425441",
              "5505010502":"2022-07-24T23:59:07.273212"}
    with (HERE/"data/pointing.csv").open() as file:
        entries = {r["obsid"]:r for r in csv.DictReader(file)}
    with solar_system_ephemeris.set(str(cache.resolve()/"de430.bsp")):
        for obs,label in labels.items():
            info=entries[obs]
            h=fits.getheader(cache/f"ni{obs}_0mpu7_cl.evt.gz","EVENTS")
            t,bary=barycentric(float(info["spacecraft_event_time"]),h,cache/f"ni{obs}.orb.gz",
                               float(info["source_ra"]),float(info["source_dec"]))
            bary.precision=6
            delta=(bary-Time(label,scale="tdb")).sec
            print(f"{obs}  spacecraft UTC {t.utc.isot}  SSB TDB {bary.isot}  "
                  f"offset from Figure 3 TDB label {delta:.6f} s")
            if abs(delta-69.184)>0.003:
                raise ValueError("Expected offset differs by more than 3 milliseconds")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    parser.add_argument("--download",action="store_true")
    args=parser.parse_args()
    inputs(args.cache,args.download)
    main(args.cache)
