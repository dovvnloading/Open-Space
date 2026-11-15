"""Two-body, aligned-spin apsidal structure constant, using Table 2 of Liakos (2025)."""
import json
from pathlib import Path
import numpy as np

def structure(eccentricity, apsidal_rate, spin_ratio=1., radius=(.177,.162), q=1.037):
    """Return the c2-weighted stellar apsidal constant; rates are degrees/year."""
    q_companion = np.array([q, 1/q])
    radius = np.array(radius)
    rotation = spin_ratio**2 * (1+q_companion)/(1-eccentricity**2)**2
    tidal_shape = (1+1.5*eccentricity**2+eccentricity**4/8)/(1-eccentricity**2)**5
    coefficients = radius**5*(rotation+15*q_companion*tidal_shape)
    # GR scales with (1-e^2)^-1 for fixed period and total mass.
    relativistic = .1552497004015139*(1-.13**2)/(1-eccentricity**2)
    deg_per_cycle = (apsidal_rate-relativistic) * 1.8434336/365.25
    k2 = deg_per_cycle/(360*coefficients.sum())
    return dict(e=eccentricity, rate_deg_yr=apsidal_rate, spin_ratio=spin_ratio,
                gr_deg_yr=relativistic, c21=coefficients[0], c22=coefficients[1],
                mean_k2=k2, log10_mean_k2=np.log10(k2))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('orbit_json', type=Path)
    args=parser.parse_args()
    orbit=json.loads(args.orbit_json.read_text())
    for e in (.12,.13,.14):
        r=orbit[str(e)]['rate_deg_per_year']
        print(json.dumps(structure(e,r)))
