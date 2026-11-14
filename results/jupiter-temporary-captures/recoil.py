"""Sensitivity of the OE31 forecast to assumed cometary recoil from the fit epoch."""
import argparse
import csv
from math import sqrt
from pathlib import Path

import numpy as np

from integrate import AU, DAY, FIELDS, WINDOWS, initialize, summarize


def propagate(law, radial, transverse):
    sim, epoch, gm, solar_gm = initialize("oe31", "jpl", 0)

    def apply_recoil(pointer):
        particles = pointer.contents.particles
        sun, target = particles[0], particles[9]
        x, y, z = target.x-sun.x, target.y-sun.y, target.z-sun.z
        r = sqrt(x*x + y*y + z*z)
        x, y, z = x/r, y/r, z/r
        vx, vy, vz = target.vx-sun.vx, target.vy-sun.vy, target.vz-sun.vz
        projection = x*vx + y*vy + z*vz
        tx, ty, tz = vx-projection*x, vy-projection*y, vz-projection*z
        vt = sqrt(tx*tx + ty*ty + tz*tz)
        if law == "water":
            u = r / 2.808
            activity = 0.111262 * u**-2.15 * (1 + u**5.093)**-4.6142
        else:
            activity = 1 / r**2
        target.ax += activity * (radial*x + transverse*tx/vt)
        target.ay += activity * (radial*y + transverse*ty/vt)
        target.az += activity * (radial*z + transverse*tz/vt)

    sim.additional_forces = apply_recoil
    sim.force_is_velocity_dependent = 1
    start, stop = WINDOWS["oe31"][1]
    times = np.arange(start, stop+1, 5.)
    radii, energies, hills, longitudes, planet_longitudes = [], [], [], [], []
    for t in times:
        sim.integrate(t-epoch)
        sun, jupiter, target = sim.particles[0], sim.particles[5], sim.particles[9]
        dr = np.array(target.xyz)-jupiter.xyz
        rj = np.array(jupiter.xyz)-sun.xyz
        dv = (np.array(target.vxyz)-jupiter.vxyz)*AU/DAY
        r = np.linalg.norm(dr)
        radii.append(r)
        energies.append(dv@dv/2 - gm/(r*AU))
        hills.append(np.linalg.norm(rj)*(gm/(3*solar_gm))**(1/3))
        longitudes.append(np.arctan2(dr[1], dr[0]))
        planet_longitudes.append(np.arctan2(rj[1], rj[0]))
    angle = np.unwrap(longitudes)
    return summarize(times, np.array(radii), np.array(energies), np.array(hills),
                     angle, angle-np.unwrap(planet_longitudes))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with args.output.open("w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["distance_law","radial_at_1au","transverse_at_1au",*FIELDS])
        for law in ("inverse_square","water"):
            cases = [(0., 0.)]
            for strength in (1e-9, 1e-8, 1e-7):
                cases.extend([(strength,0.), (0.,strength), (0.,-strength)])
            for radial, transverse in cases:
                values = propagate(law, radial, transverse)
                writer.writerow([law, radial, transverse, *[f"{v:.7f}" for v in values]])
                out.flush()


if __name__ == "__main__":
    main()
