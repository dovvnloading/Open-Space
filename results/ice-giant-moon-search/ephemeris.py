"""Cache a short, densely sampled spacecraft ephemeris for each recognized moon."""
import argparse
from pathlib import Path
import numpy as np
from geometry import AU, C, MOONS, PLANET, direction, query


def build(planet, observer, start_mjd, end_mjd, cache, output, spk=None):
    from astropy.time import Time
    dates = np.linspace(start_mjd, end_mjd, max(4, int(np.ceil((end_mjd-start_mjd)*1440/6))+1))
    if not 0 < end_mjd-start_mjd < 2:
        raise ValueError("use one continuous visit group of less than two days")
    ids = [PLANET[planet], *MOONS[planet]]
    values = []
    for body in ids:
        if body == 75052 and spk is not None:
            from moon_spk import offsets
            planet_data = values[0]
            tdb = Time(dates, format="mjd", scale="utc").tdb.jd
            xyz = offsets(spk, tdb, planet_data)
            delta = np.linalg.norm(xyz, axis=1)
            ra = np.rad2deg(np.arctan2(xyz[:,1], xyz[:,0])) % 360
            dec = np.rad2deg(np.arcsin(xyz[:,2]/delta))
            values.append(np.column_stack((ra, dec, np.full_like(ra,np.nan),planet_data[:,3],delta/AU)))
        else:
            values.append(query(body, observer, dates+2400000.5, cache))
        print(f"{observer} {body} {len(dates)} epochs", flush=True)
    np.savez_compressed(output, mjd=dates, ids=ids, values=values, planet=planet, observer=observer)


def interpolate(filename, dates):
    from scipy.interpolate import CubicSpline
    with np.load(filename) as f:
        ids, times, src = f["ids"], f["mjd"], f["values"]
        planet, observer = str(f["planet"]), str(f["observer"])
    if list(ids) != [PLANET[planet], *MOONS[planet]]:
        raise ValueError("catalog incomplete or outdated")
    if np.min(dates)<times[0] or np.max(dates)>times[-1]:
        raise ValueError("ephemeris does not bracket exposures")
    src[:,:,0] = np.rad2deg(np.unwrap(np.deg2rad(src[:,:,0]),axis=1))
    diam = src[:,:,2]
    if not np.isfinite(diam[0]).all():
        raise ValueError("missing planetary angular diameter")
    if np.any(np.isfinite(diam).any(axis=1) & ~np.isfinite(diam).all(axis=1)):
        raise ValueError("intermittent missing angular diameter")
    src[:,:,2] = np.nan_to_num(diam)  # moons without a measured diameter get the PSF mask
    return planet, observer, ids, CubicSpline(times, src, axis=1)(dates)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("planet", choices=PLANET)
    p.add_argument("observer", choices=["JWST","HST"])
    p.add_argument("start_mjd", type=float)
    p.add_argument("end_mjd", type=float)
    p.add_argument("cache")
    p.add_argument("output")
    p.add_argument("--ura184", help="local excerpt or full ura184.bsp for S/2025 U1")
    a = p.parse_args()
    build(a.planet, a.observer, a.start_mjd, a.end_mjd, a.cache, a.output, a.ura184)
