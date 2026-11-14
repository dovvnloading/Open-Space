"""Propagate published JPL and MPC solutions and correlated orbit samples."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import rebound

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
AU = 149597870.7
DAY = 86400
SAMPLES = 128
FIELDS = ("epoch_jd", "years_negative_energy", "years_bound_inside_hill",
          "closest_au", "min_energy_km2_s2", "net_turns_inertial",
          "net_turns_rotating", "negative_energy_episodes",
          "max_inertial_excursion_turns_any_episode")
# Windows deliberately extend beyond both ends of each nominal episode.
WINDOWS = {
    "uq146": [(2444970.5, 2447527.5), (2476017.5, 2479304.5)],
    "oe31": [(2453736.5, 2457754.5), (2485158.5, 2493329.5)],
}


def read_json(name):
    return json.loads((DATA / name).read_text())


def initialize(target, source, count):
    fit = read_json(f"{target}-{source}.json")
    if source == "jpl":
        fit = fit["orbit"]["covariance"]
        assert fit["labels"] == ["e", "q", "tp", "node", "peri", "i"]
        epoch = float(fit["epoch"])
        mean = [float(v["value"]) for v in fit["elements"]]
        covariance = np.array(fit["data"], dtype=float)
    else:
        fit = fit[0]["mpc_orb"][0]
        epoch = float(fit["epoch_data"]["epoch"]) + 2400000.5
        cart = fit["CAR"]
        mean = cart["coefficient_values"]
        covariance = np.array([[cart["covariance"][f"cov{min(i,j)}{max(i,j)}"]
                                for j in range(6)] for i in range(6)])
    rng = np.random.default_rng(20261114)
    samples = np.vstack([mean, mean + rng.standard_normal((count, 6)) @
                         np.linalg.cholesky(covariance).T])
    states = np.array(read_json("planets.json")[str(epoch)])
    sim = rebound.Simulation()
    sim.G = states[0, 0] * DAY**2 / AU**3
    sim.integrator = "ias15"
    for gm, x, y, z, vx, vy, vz in states:
        sim.add(m=gm / states[0, 0], x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    sim.N_active = sim.N
    for row in samples:
        if source == "jpl":
            e, q, tp, node, peri, inc = row
            a = q / (1 - e)
            sim.add(primary=sim.particles[0], a=a, e=e,
                    M=np.sqrt(sim.G / a**3) * (epoch - tp),
                    inc=np.radians(inc), Omega=np.radians(node), omega=np.radians(peri))
        else:
            x, y, z, vx, vy, vz = row + states[0, 1:]
            sim.add(x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    sim.move_to_com()
    return sim, epoch, states[5, 0], states[0, 0]


def zero_crossing(t0, t1, value0, value1):
    """Interpolate a bracketed zero in the same time and value convention."""
    return t0 - value0 * (t1-t0) / (value1-value0)


def negative_time(t, values):
    """Integrate time spent below zero on a piecewise-linear trajectory."""
    negative = values < 0
    crossing = negative[1:] != negative[:-1]
    dt = np.diff(t)
    total = dt[negative[:-1] & negative[1:]].sum()
    starts, ends = t[:-1][crossing], t[1:][crossing]
    roots = zero_crossing(starts, ends, values[:-1][crossing], values[1:][crossing])
    return total + np.where(negative[:-1][crossing], roots-starts, ends-roots).sum()


def summarize(t, radius, energy, hill, longitude, rotating):
    """Measure the deepest episode and the largest angular excursion of any episode."""
    negative = energy < 0
    if negative[0] or negative[-1]:
        raise ValueError("A negative-energy episode reaches the sampling boundary; "
                         "extend the time window to bracket entry and exit")
    if not np.any(negative):
        return (np.nan, 0, 0, min(radius), min(energy), 0, 0, 0, 0)
    entries = np.flatnonzero(~negative[:-1] & negative[1:]) + 1
    exits = np.flatnonzero(negative[:-1] & ~negative[1:])
    deepest = int(np.argmin(energy))
    max_excursion = 0.0
    for lo, hi in zip(entries, exits):
        start = zero_crossing(t[lo-1], t[lo], energy[lo-1], energy[lo])
        stop = zero_crossing(t[hi], t[hi+1], energy[hi], energy[hi+1])
        # Include every sampled vertex as well as the interpolated boundaries.
        times = np.r_[start, t[lo:hi+1], stop]
        angles = np.interp(times, t, longitude)
        max_excursion = max(max_excursion, np.ptp(angles) / (2*np.pi))
        if lo <= deepest <= hi:
            gap = np.interp(times, t, radius-hill)
            hill_days = negative_time(times, gap)
            imin = lo + int(np.argmin(radius[lo:hi+1]))
            lon_rot = np.interp([start, stop], t, rotating)
            result = (t[imin], (stop-start)/365.25, hill_days/365.25,
                      radius[imin], energy[deepest], (angles[-1]-angles[0])/(2*np.pi),
                      (lon_rot[1]-lon_rot[0])/(2*np.pi))
    return (*result, len(entries), max_excursion)


def propagate(target, source, count, window):
    sim, epoch, gm, solar_gm = initialize(target, source, count)
    times = np.arange(window[0], window[1]+1, 5.0)
    sim.integrate(times[0]-epoch)
    pos = np.empty((sim.N, 3))
    vel = np.empty_like(pos)
    radius, energy, hill, theta, jlong = [], [], [], [], []
    for jd in times:
        sim.integrate(jd-epoch)
        sim.serialize_particle_data(xyz=pos, vxvyvz=vel)
        dr, dv = pos[9:]-pos[5], (vel[9:]-vel[5])*AU/DAY
        r = np.linalg.norm(dr, axis=1)
        rj = pos[5]-pos[0]
        radius.append(r)
        energy.append(np.sum(dv*dv, axis=1)/2 - gm/(AU*r))
        hill.append(np.linalg.norm(rj)*(gm/(3*solar_gm))**(1/3))
        theta.append(np.arctan2(dr[:, 1], dr[:, 0]))
        jlong.append(np.arctan2(rj[1],rj[0]))
    theta = np.unwrap(theta, axis=0)
    rotating = theta - np.unwrap(jlong)[:, None]
    return [summarize(times, r, e, np.array(hill), lon, rot)
            for r, e, lon, rot in zip(np.array(radius).T, np.array(energy).T,
                                      theta.T, rotating.T)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--samples", type=int, default=SAMPLES)
    args = parser.parse_args()
    with args.output.open("w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["target","fit","window_start_jd","sample",*FIELDS])
        for target, windows in WINDOWS.items():
            for source in ("jpl","mpc"):
                for bounds in windows:
                    for ix, row in enumerate(propagate(target,source,args.samples,bounds)):
                        writer.writerow([target,source,bounds[0],ix,
                                         *[f"{v:.7f}" for v in row]])
                    out.flush()


if __name__ == "__main__":
    main()
