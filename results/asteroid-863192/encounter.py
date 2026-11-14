"""Propagate the astrometric fits and sensitivity controls backward in time."""

import numpy as np
import rebound

from orbits import AU_KM, DAY, GM, IDS, EPOCH, EQUATOR_TO_ECLIPTIC, DATA
import json


def solar_system(gravity, massive_asteroids=False):
    sim = rebound.Simulation()
    sim.G = GM[0] * DAY**2 / AU_KM**3
    for gm, body in zip(GM, IDS):
        pos, vel = gravity.state(0, body, EPOCH) @ EQUATOR_TO_ECLIPTIC.T
        sim.add(
            m=gm / GM[0], x=pos[0], y=pos[1], z=pos[2], vx=vel[0], vy=vel[1], vz=vel[2]
        )
    if massive_asteroids:
        for body in json.loads((DATA / "perturbers.json").read_text()):
            state = np.array(body["barycentric_icrf_au_au_per_day"]).reshape(2, 3)
            pos, vel = state @ EQUATOR_TO_ECLIPTIC.T
            sim.add(
                m=float(body["gm"]["value"]) / GM[0],
                x=pos[0],
                y=pos[1],
                z=pos[2],
                vx=vel[0],
                vy=vel[1],
                vz=vel[2],
            )
    sim.N_active = sim.N
    return sim


def mpc_state(line, gravity):
    sim = solar_system(gravity)
    mean, peri, node, inc = np.deg2rad(
        [float(line[a:b]) for a, b in [(26, 35), (37, 46), (48, 57), (59, 68)]]
    )
    sim.add(
        primary=sim.particles[0],
        a=float(line[92:103]),
        e=float(line[70:79]),
        inc=inc,
        Omega=node,
        omega=peri,
        M=mean,
    )
    particle = sim.particles[-1]
    return (np.array([particle.xyz, particle.vxyz]) @ EQUATOR_TO_ECLIPTIC).ravel()


def uncertainty_root(rows, jacobian):
    site_date = rows.site + " " + rows.utc.str[:10]
    weights = site_date.map(site_date.value_counts()).to_numpy(dtype=float) ** -0.5
    _, singular, right = np.linalg.svd(
        jacobian * np.repeat(weights, 2)[:, None], full_matrices=False
    )
    return right.T / singular


def encounter_grid(
    gravity,
    first,
    second,
    step_days,
    drift_au_myr=None,
    lookbacks=None,
    massive_asteroids=False,
):
    sim = solar_system(gravity, massive_asteroids)
    active = sim.N_active
    for state in np.vstack((first, second)):
        pos, vel = state.reshape(2, 3) @ EQUATOR_TO_ECLIPTIC.T
        sim.add(x=pos[0], y=pos[1], z=pos[2], vx=vel[0], vy=vel[1], vz=vel[2])
    sim.move_to_com()
    sim.integrator = "whfast"
    sim.dt = -step_days
    sim.ri_whfast.safe_mode = 0
    pos = np.empty((sim.N, 3))
    vel = np.empty_like(pos)
    count = len(first)
    best = np.tile([np.inf, 0.0, 0.0], (count, 1))
    track = []
    previous_time = 0.0
    if drift_au_myr is not None:
        # At each output apply a velocity-parallel impulse. Gauss's energy
        # equation da/dt = 2 a^2 (v dot A)/mu sets the imposed drift exactly
        # to first order. Use the actual signed elapsed time for backward kicks.
        rates = np.r_[np.zeros(count), drift_au_myr] / (1e6 * 365.25)
    if lookbacks is None:
        lookbacks = np.arange(0.0, 4000.0, 0.1)
    for year in lookbacks:
        sim.integrate(-year * 365.25, exact_finish_time=0)
        sim.synchronize()
        sim.serialize_particle_data(xyz=pos, vxvyvz=vel)
        distance = (
            np.linalg.norm(pos[active : active + count] - pos[active + count :], axis=1)
            * AU_KM
        )
        update = distance < best[:, 0]
        if update.any():
            speed = (
                np.linalg.norm(
                    vel[active : active + count] - vel[active + count :], axis=1
                )
                * AU_KM
                * 1000
                / DAY
            )
            best[update] = np.column_stack(
                (distance, np.full(count, -sim.t / 365.25), speed)
            )[update]
        if drift_au_myr is not None and previous_time != sim.t:
            relative_pos = pos[active:] - pos[0]
            relative_vel = vel[active:] - vel[0]
            speed_sq = (relative_vel**2).sum(axis=1)
            mu = sim.G * sim.particles[0].m
            semimajor = 1 / (2 / np.linalg.norm(relative_pos, axis=1) - speed_sq / mu)
            scale = mu * rates * (sim.t - previous_time) / (2 * semimajor**2 * speed_sq)
            vel[active:] += relative_vel * scale[:, None]
            sim.set_serialized_particle_data(vxvyvz=vel)
            sim.ri_whfast.recalculate_coordinates_this_timestep = 1
        previous_time = sim.t
        if len(track) == 0 or (-sim.t / 365.25 - track[-1][0]) >= 5:
            track.append((-sim.t / 365.25, distance[0]))
    return best, track
