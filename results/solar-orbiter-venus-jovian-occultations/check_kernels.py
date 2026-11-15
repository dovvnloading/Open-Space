"""Check Solar Orbiter occultation contacts against two ESA trajectory kernels."""
import argparse
import hashlib
import urllib.request
from pathlib import Path

import numpy as np
import spiceypy as spice
from scipy.optimize import brentq

from reproduce import CASES, VENUS_RADIUS_KM, write_csv

MIRROR = "https://datacenter.stix.i4ds.net/pub/spice/latest/kernels/"
OCTOBER = "solo_ANC_soc-orbit-stp_20200210-20301120_432_V1_00588_V01.bsp"
NOVEMBER = "solo_ANC_soc-orbit-stp_20200210-20301120_436_V1_00593_V01.bsp"
KERNELS = {
    "lsk/naif0012.tls": "678e32bdb5a744117a467cd9601cd6b373f0e9bc9bbde1371d5eee39600a039b",
    "spk/de421.bsp": "a20a7139da04cbc462454634918e9a9ca69127044e2cc9d4f9c16e238d2deedc",
    "spk/" + OCTOBER: "faf21183107b18f57e27a1789db7a4ec75a9a75705122eebaf4b0a201076e5f1",
    "spk/" + NOVEMBER: "ef0ce3146c9df4ed78975450233135f21e8207e9c98ee4521c88b6f3f8d198e4",
}


def ensure_kernels(cache):
    cache.mkdir(parents=True, exist_ok=True)
    for name, digest in KERNELS.items():
        path = cache / Path(name).name
        if not path.exists():
            with urllib.request.urlopen(MIRROR + name, timeout=120) as response:
                content = response.read()
            if hashlib.sha256(content).hexdigest() != digest:
                raise ValueError(f"Version mismatch for {name}")
            path.write_bytes(content)
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Version mismatch for {path}")


def contacts(start, stop, body):
    def limb(seconds):
        venus, _ = spice.spkpos("299", start + seconds, "J2000", "CN+S", "-144")
        background, _ = spice.spkpos(body, start + seconds, "J2000", "CN+S", "-144")
        return spice.vsep(venus, background) - np.arcsin(VENUS_RADIUS_KM / spice.vnorm(venus))

    grid = np.arange(0, stop - start + 1, 10)
    values = np.array([limb(t) for t in grid])
    roots = [brentq(limb, grid[i], grid[i + 1]) for i in np.where(np.diff(np.signbit(values)))[0]]
    if len(roots) not in (0, 2):
        raise ValueError("Window does not bracket both contacts")
    return [spice.et2utc(start + s, "ISOC", 3) + "Z" for s in roots] if roots else ["", ""]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True, help="Directory for approximately 54 MB of ESA kernels")
    cache = parser.parse_args().cache
    ensure_kernels(cache)
    rows = []
    for orbit in (OCTOBER, NOVEMBER):
        spice.kclear()
        for name in ["naif0012.tls", "de421.bsp", orbit]:
            spice.furnsh(str(cache / name))
        for day, (_, observer, first, last) in CASES.items():
            if observer != "-144":
                continue
            start, stop = [spice.str2et(day + "T" + minute) for minute in (first, last)]
            for name, body in [("earth", "399"), ("jupiter_barycenter", "5"), ("sun", "10")]:
                rows.append([orbit, day, name, *contacts(start, stop, body)])
    write_csv("kernel_contacts.csv", ["orbit_kernel", "date", "direction", "ingress_utc", "egress_utc"], rows)

    # Earth and spacecraft arrival times for the same outgoing Jupiter wavefront.
    representative = {"2021-08-09": "04:33:00", "2025-02-18": "20:49:00",
                      "2026-12-24": "23:04:00", "2028-03-18": "08:21:00"}
    views = []
    for day, clock in representative.items():
        reception = spice.str2et(day + "T" + clock)
        _, light_time = spice.spkpos("5", reception, "J2000", "CN", "-144")
        emission = reception - light_time
        solo, _ = spice.spkpos("-144", emission, "J2000", "XCN", "5")
        earth, light_time_to_earth = spice.spkpos("399", emission, "J2000", "XCN", "5")
        views.append([day, clock, round(float(np.degrees(spice.vsep(solo, earth))), 5),
                      round(light_time_to_earth - light_time, 3), round(light_time, 3)])
    write_csv("observer_geometry.csv", ["date", "reception_utc_on_solo", "view_separation_at_jupiter_deg",
                                        "earth_minus_solo_arrival_s", "jupiter_to_solo_light_time_s"], views)
    spice.kclear()


if __name__ == "__main__":
    main()
