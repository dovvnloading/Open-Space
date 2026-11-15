"""Replay nominal Earth-launch and planetary-encounter checks for VE138/VV15."""
import argparse
import csv
from datetime import datetime, timedelta
import hashlib
import json
from pathlib import Path

from jplephem.spk import SPK
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize_scalar

AU_KM = 149597870.7
DAY = 86400.
C = 299792.458 * DAY / AU_KM
GM_CONVERSION = DAY**2 / AU_KM**3
GM_SUN = 132712440041.279 * GM_CONVERSION
# JPL DE440 mass parameters in km^3/s^2. Barycenters include satellites.
BODIES = (
    ("Mercury", 1, 22031.86855), ("Venus", 2, 324858.592),
    ("Earth", 399, 398600.43551), ("Moon", 301, 4902.800118),
    ("Mars", 4, 42828.37582), ("Jupiter", 5, 126712764.1),
    ("Saturn", 6, 37940584.84), ("Uranus", 7, 5794556.4),
    ("Neptune", 8, 6836527.1),
)
GM_PLANETS = np.array([item[2] for item in BODIES]) * GM_CONVERSION
FIRST_JD, LAST_JD = 2436115.5, 2461359.5  # 1957 Oct 4 through 2026 Nov 15 TDB
KERNEL_HASH = "c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2"


def planet_ephemeris(path, epoch, step):
    with path.open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != KERNEL_HASH:
            raise ValueError("Expected the specified JPL DE440s kernel")
    offsets = np.arange(FIRST_JD-epoch-step, LAST_JD-epoch+2*step, step)
    jd = epoch + offsets
    eps = np.deg2rad(23.439291111)
    rotation = np.array([[1, 0, 0], [0, np.cos(eps), np.sin(eps)],
                         [0, -np.sin(eps), np.cos(eps)]])
    with SPK.open(str(path)) as kernel:
        sun = np.array(kernel[0, 10].compute_and_differentiate(jd))
        planets = []
        for _, naif, _ in BODIES:
            if naif in (399, 301):
                pv = (np.array(kernel[0, 3].compute_and_differentiate(jd))
                      + np.array(kernel[3, naif].compute_and_differentiate(jd)))
            else:
                pv = np.array(kernel[0, naif].compute_and_differentiate(jd))
            # Order is position/velocity, sample time, Cartesian component.
            planets.append(np.einsum("ij,sjt->sti", rotation, (pv-sun)/AU_KM))
    pv = np.stack(planets, axis=2)
    return CubicHermiteSpline(offsets, pv[0], pv[1], axis=0)


def orbit_trajectory(state, radial, transverse, planets, epoch):
    def acceleration(time, y):
        positions = planets(time)
        r, v = y[:3], y[3:]
        radius = np.linalg.norm(r)
        delta = positions-r
        force = -(GM_SUN-radial)*r/radius**3
        force += delta.T @ (GM_PLANETS/np.linalg.norm(delta, axis=1)**3)
        force -= positions.T @ (GM_PLANETS/np.linalg.norm(positions, axis=1)**3)
        force += GM_SUN/C**2/radius**3*((4*GM_SUN/radius-v@v)*r+4*(r@v)*v)
        if transverse:
            tangent = v-r*(v@r)/radius**2
            force += transverse*tangent/(radius**2*np.linalg.norm(tangent))
        return np.concatenate((v, force))

    def solution(stop):
        result = solve_ivp(acceleration, (0., stop), state, method="DOP853",
                           rtol=3e-12, atol=5e-15, dense_output=True, max_step=7.)
        if not result.success:
            raise RuntimeError(result.message)
        return result.sol

    return ((FIRST_JD-epoch, 0., solution(FIRST_JD-epoch)),
            (0., LAST_JD-epoch, solution(LAST_JD-epoch)))


def closest_approaches(arcs, planets, epoch):
    minima = [[] for _ in BODIES]
    for start, end, trajectory in arcs:
        grid = np.linspace(start, end, int(np.ceil((end-start)*4))+1)
        relative = trajectory(grid)[:3].T[:, None, :]-planets(grid)
        distances = np.linalg.norm(relative, axis=2)
        for idx, (body, _, _) in enumerate(BODIES):
            d = distances[:, idx]
            local = np.flatnonzero((d[1:-1]<d[:-2]) & (d[1:-1]<d[2:]))+1
            # Only interior minima can be genuine encounters. Include endpoints
            # when identifying the minimum separation over the full interval.
            def separation(t):
                return np.linalg.norm(trajectory(t)[:3]-planets(t)[idx])
            candidates = [(grid[0], False), (grid[-1], False)]
            for i in local:
                if d[i] < .16 or d[i] <= np.min(d)+.01:
                    fit = minimize_scalar(separation, bounds=(grid[i-1], grid[i+1]),
                                          method="bounded")
                    candidates.append((fit.x, True))
            for t, interior in candidates:
                speed = np.linalg.norm(trajectory(t)[3:]-planets(t, 1)[idx])*AU_KM/DAY
                date = datetime(2000,1,1,12)+timedelta(days=float(epoch+t-2451545))
                minima[idx].append((body, date.isoformat(timespec="minutes"),
                                    float(separation(t)), float(speed), interior))
    return minima


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("kernel", type=Path)
    cli.add_argument("gravity_fit", type=Path, help="VE138 joint six-state JSON")
    cli.add_argument("acceleration_fits", type=Path, help="VE138 A1/A12/A12_site JSON")
    cli.add_argument("--step",type=float,default=.5,help="ephemeris tabulation step in days")
    cli.add_argument("--out",type=Path, default=Path(__file__).resolve().parent)
    args = cli.parse_args()
    if not 0 < args.step <= 1.:
        cli.error("--step must be in (0, 1] days")
    base=json.loads(args.gravity_fit.read_text())
    accel=json.loads(args.acceleration_fits.read_text())
    epoch=float(base["epoch_jd_tdb"])
    models={"gravity": [*base["full"]["state"],0.,0.]}
    if epoch != float(accel["epoch_jd_tdb"]):
        raise ValueError("Orbit epochs differ")
    for key, fit in (("A1", "all/A1"), ("A12", "all/A12"), ("A12_site", "sites/A12")):
        models[key] = accel[fit]["state"]
    planets=planet_ephemeris(args.kernel, epoch, args.step)
    args.out.mkdir(parents=True,exist_ok=True)
    with (args.out/"encounters.csv").open("w", newline="") as stream:
        writer=csv.writer(stream,lineterminator="\n")
        writer.writerow(["model","selection","body","date_tdb","distance_au","relative_km_s"])
        for model, values in models.items():
            arcs=orbit_trajectory(values[:6],values[6],values[7],planets,epoch)
            minima=closest_approaches(arcs,planets,epoch)
            for entries in minima:
                body,date,distance,speed,_=min(entries,key=lambda x:x[2])
                writer.writerow([model,"minimum",body,date,f"{distance:.7f}",f"{speed:.4f}"])
            for body,date,distance,speed,interior in sorted(minima[2],key=lambda x:x[1]):
                if interior and distance < .15:
                    writer.writerow([model,"earth_lt_0.15_au",body,date,f"{distance:.7f}",f"{speed:.4f}"])

if __name__=="__main__":
    main()
