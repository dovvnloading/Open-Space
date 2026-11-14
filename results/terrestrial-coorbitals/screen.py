"""Replay the three discovery screens from the frozen, eligible SBDB rows."""
import csv
import json
import math
from pathlib import Path

import numpy as np
import rebound

from replay import G, MASS

ROOT = Path(__file__).resolve().parent
BODIES = (10, 1, 2, 3, 4, 5, 6, 7, 8)
with (ROOT / 'screen.csv').open() as fp:
    rows = list(csv.DictReader(fp))
states = json.loads((ROOT / 'orbits.json').read_text())['planet_vectors']['2461200.5']

for planet, idx in [('earth', 3), ('mars', 4), ('venus', 2)]:
    eligible = [row for row in rows if row['planet'] == planet]
    sim = rebound.Simulation()
    sim.G = G
    for body in BODIES:
        x, y, z, vx, vy, vz = states[str(body)]
        mass = 1 / 328900.5597 if body == 3 else MASS[body]
        sim.add(m=mass, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    sim.N_active = sim.N
    sim.integrator = 'mercurius'
    sim.dt = 4
    for row in eligible:
        inc, node, peri, mean = [math.radians(float(row[k])) for k in ('i','om','w','ma')]
        sim.add(primary=sim.particles[0], a=float(row['a']), e=float(row['e']),
                inc=inc, Omega=node, omega=peri, M=mean)
    sim.move_to_com()
    angles = []
    for year in range(301):
        sim.integrate(year * 365.25)
        sun = sim.particles[0]
        longitude = sim.particles[idx].orbit(primary=sun).l
        angles.append([(math.degrees(p.orbit(primary=sun).l - longitude) + 180) % 360 - 180
                       for p in list(sim.particles)[sim.N_active:]])
    spans = np.ptp(np.degrees(np.unwrap(np.radians(angles), axis=0)), axis=0)
    # Long integrations of close-encounter objects can diverge across platforms.
    saved = np.array([float(r['unwrapped_span_300yr_deg']) for r in eligible])
    print(planet, len(eligible), 'median absolute span difference (deg)',
          round(float(np.median(np.abs(saved - spans))), 6))
