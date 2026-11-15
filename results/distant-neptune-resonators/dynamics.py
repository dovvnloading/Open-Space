"""Integrate an SBDB orbital solution against the four giant-planet systems."""

from pathlib import Path
import json
import math

import numpy as np
import rebound
from jplephem.spk import SPK

ROOT = Path(__file__).resolve().parent
WORK = ROOT / "work"
AU = 149597870.7  # km
YEAR = 365.25  # days
DAY = 86400.0
GM_SUN = 132712440041.279419  # km^3 s^-2
GM_INNER = [22031.868551, 324858.592000, 398600.435507 + 4902.800118, 42828.375816]
GM_GIANTS = [126712764.100000, 37940584.841800, 5794556.400000, 6836527.100580]
EPS = np.deg2rad(84381.448 / 3600)
ROT = np.array(
    [[1.0, 0.0, 0.0], [0.0, np.cos(EPS), np.sin(EPS)], [0.0, -np.sin(EPS), np.cos(EPS)]]
)


def states(epoch):
    """Barycentric mean ecliptic of J2000, au and au/Julian year."""
    with SPK.open(WORK / "de440s.bsp") as kernel:

        def vector(body):
            pos, vel = kernel[0, body].compute_and_differentiate(epoch)
            return np.r_[ROT @ pos / AU, ROT @ vel * YEAR / AU]

        sun = vector(10)
        inner = np.average(
            [sun] + [vector(i) for i in range(1, 5)],
            axis=0,
            weights=[GM_SUN] + GM_INNER,
        )
        return sun, inner, [vector(i) for i in range(5, 9)]


def init(epoch, step):
    sim = rebound.Simulation()
    sim.G = GM_SUN * (YEAR * DAY) ** 2 / AU**3
    sun, inner, giants = states(epoch)
    for mass, v in zip(
        [1 + sum(GM_INNER) / GM_SUN] + [g / GM_SUN for g in GM_GIANTS], [inner] + giants
    ):
        sim.add(m=mass, x=v[0], y=v[1], z=v[2], vx=v[3], vy=v[4], vz=v[5])
    primary = rebound.Particle(
        m=1.0, x=sun[0], y=sun[1], z=sun[2], vx=sun[3], vy=sun[4], vz=sun[5]
    )
    sim.N_active = 5
    sim.integrator = "whfast"
    sim.dt = step
    return sim, primary


def add(sim, sun, a, e, inc, node, peri, anomaly):
    sim.add(
        primary=sun,
        a=a,
        e=e,
        inc=np.deg2rad(inc),
        Omega=np.deg2rad(node),
        omega=np.deg2rad(peri),
        M=np.deg2rad(anomaly),
    )


def integrate(sim, years, interval):
    """Jacobi elements; massless orbits are referred to the massive-system barycenter."""
    times = np.arange(0, abs(years) + 1, interval) * np.sign(years)
    counts = sim.N - sim.N_active
    data = np.empty((len(times), counts, 4))
    for k, t in enumerate(times):
        sim.integrate(t, exact_finish_time=int(sim.integrator != "whfast"))
        orbits = sim.orbits()
        neptune = orbits[3]
        data[k] = [(o.a, o.l, o.pomega, neptune.l) for o in orbits[4:]]
    return times, data


def angle(data, p, q):
    return (p * data[:, 1] - q * data[:, 3] - (p - q) * data[:, 2]) % (2 * np.pi)


def arc(phi):
    ordered = np.sort(np.rad2deg(phi) % 360.0)
    gaps = np.diff(ordered, append=ordered[0] + 360.0)
    largest = np.argmax(gaps)
    width = 360 - gaps[largest]
    center = (ordered[(largest + 1) % len(ordered)] + width / 2) % 360
    return width, center


def ratios(period_ratio):
    for q in (1, 2, 3):
        for p in range(
            max(q + 1, math.floor(q * period_ratio) - 1),
            math.ceil(q * period_ratio) + 2,
        ):
            if math.gcd(p, q) == 1:
                yield p, q


def screen(years=400000, step=0.5, interval=100):
    catalog = json.loads((WORK / "catalog.json").read_text())
    rows = [dict(zip(catalog["fields"], r)) for r in catalog["data"]]
    (epoch,) = {float(r["epoch"]) for r in rows}
    sim, sun = init(epoch, step)
    for r in rows:
        add(sim, sun, *(float(r[k]) for k in ("a", "e", "i", "om", "w", "ma")))
    times, data = integrate(sim, years, interval)
    np.savez_compressed(WORK / "nominal.npz", times=times, data=data)
    return summarize(rows, data)


def summarize(rows, data):
    import csv

    with (WORK / "screen.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "designation",
                "spkid",
                "median_a_au",
                "p",
                "q",
                "arc_deg",
                "center_deg",
                "arc_first_half",
                "arc_second_half",
            ]
        )
        for i, r in enumerate(rows):
            series = data[:, i, :]
            a = np.median(series[:, 0])
            for p, q in ratios((a / 30.07) ** 1.5):
                phi = angle(series, p, q)
                width, center = arc(phi)
                h = len(phi) // 2
                writer.writerow(
                    [
                        r["full_name"].strip(),
                        r["spkid"],
                        round(a, 4),
                        p,
                        q,
                        round(width, 3),
                        round(center, 3),
                        round(arc(phi[:h])[0], 3),
                        round(arc(phi[h:])[0], 3),
                    ]
                )


if __name__ == "__main__":
    screen()
