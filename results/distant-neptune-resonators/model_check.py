"""Independent eight-planet perturbation model and matched old/new orbit comparison."""

from concurrent.futures import ProcessPoolExecutor
import csv
import numpy as np
import rebound
from jplephem.spk import SPK
from clones import setup as jpl_setup
from current import TARGETS, setup as mpc_setup
from dynamics import WORK, YEAR, DAY, AU, GM_SUN, GM_INNER, GM_GIANTS, ROT, angle, arc


def compare_old(task):
    obj, sign = task
    name, seed, p, q = obj
    from dynamics import integrate

    sim = jpl_setup({"spkid": str(seed)}, 128, sign * 0.5)
    _, data = integrate(sim, sign * 400000, 100)
    ranges = []
    for i in range(129):
        phi = np.rad2deg(np.unwrap(angle(data[:, i], p, q)))
        ranges.append([np.min(phi) - phi[0], np.max(phi) - phi[0]])
    return np.array(ranges)


def eight_planets(task):
    (sign,) = task
    epoch = 2461200.5
    sim = rebound.Simulation()
    sim.G = GM_SUN * (YEAR * DAY) ** 2 / AU**3
    with SPK.open(WORK / "de440s.bsp") as kernel:
        for idx, mass in zip(
            [10, 1, 2, 3, 4, 5, 6, 7, 8], [GM_SUN] + GM_INNER + GM_GIANTS
        ):
            pos, vel = kernel[0, idx].compute_and_differentiate(epoch)
            pos = ROT @ pos / AU
            vel = ROT @ vel * YEAR / AU
            sim.add(
                m=mass / GM_SUN,
                x=pos[0],
                y=pos[1],
                z=pos[2],
                vx=vel[0],
                vy=vel[1],
                vz=vel[2],
            )
    sim.N_active = 9
    for obj in TARGETS:
        model = mpc_setup(obj, 0, 0.5)
        particle = model.particles[5]
        sim.add(
            x=particle.x,
            y=particle.y,
            z=particle.z,
            vx=particle.vx,
            vy=particle.vy,
            vz=particle.vz,
        )
    sim.integrator = "whfast"
    sim.dt = sign * 0.02
    phases = [[], []]
    for t in np.arange(0, 400001, 100):
        sim.integrate(sign * t, exact_finish_time=0)
        orbits = sim.orbits()
        for i, (_, _, p, q) in enumerate(TARGETS):
            orbit = orbits[8 + i]
            phases[i].append(
                (p * orbit.l - q * orbits[7].l - (p - q) * orbit.pomega) % (2 * np.pi)
            )
    return [(arc(np.array(phi)), np.ptp(np.rad2deg(np.unwrap(phi)))) for phi in phases]


if __name__ == "__main__":
    with ProcessPoolExecutor(6) as pool:
        old = [pool.submit(compare_old, (obj, s)) for obj in TARGETS for s in (-1, 1)]
        full = [pool.submit(eight_planets, (s,)) for s in (-1, 1)]
        with (WORK / "old-joint.csv").open("w") as f:
            wr = csv.writer(f)
            wr.writerow(["object", "clone", "joint_800kyr_range_deg"])
            for i, obj in enumerate(TARGETS):
                back = old[2 * i].result()
                fwd = old[2 * i + 1].result()
                spans = np.maximum(back[:, 1], fwd[:, 1]) - np.minimum(
                    back[:, 0], fwd[:, 0]
                )
                wr.writerows((obj[0], k, round(v, 4)) for k, v in enumerate(spans))
        with (WORK / "eight-planet.csv").open("w") as f:
            wr = csv.writer(f)
            wr.writerow(
                [
                    "object",
                    "years",
                    "step_days",
                    "arc_deg",
                    "center_deg",
                    "unwrapped_range_deg",
                ]
            )
            for s, job in zip((-1, 1), full):
                for obj, ((width, center), span) in zip(TARGETS, job.result()):
                    wr.writerow(
                        [
                            obj[0],
                            s * 400000,
                            0.02 * YEAR,
                            round(width, 4),
                            round(center, 4),
                            round(span, 4),
                        ]
                    )
