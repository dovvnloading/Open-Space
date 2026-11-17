"""Synthetic-only check of ixpeobssim's event-wise Stokes convention."""
import numpy as np
from ixpeobssim.evt.kislat2015 import xStokesAnalysis


def run(seed, polarization, angle, n=500000, mu=0.4):
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0.0, 2*np.pi, n * 2)
    threshold = (1 + mu * polarization * np.cos(2*(phi-angle))) / (1+mu*polarization)
    phi = phi[rng.random(phi.size) < threshold][:n]
    if phi.size != n:
        raise RuntimeError("Not enough synthetic events")
    q = xStokesAnalysis.stokes_q(phi)
    u = xStokesAnalysis.stokes_u(phi, None)
    response = lambda e: np.full_like(e, mu)
    area = lambda e: np.full_like(e, 1.0)
    analysis = xStokesAnalysis(q, u, np.full(n, 4.0), response, area, 1.0, acceptcorr=False)
    i, qsum, usum = analysis.sum_stokes_parameters(2.0, 8.0)
    measured = np.array([qsum/i, usum/i])
    expected = polarization * np.array([np.cos(2*angle), np.sin(2*angle)])
    print(f"synthetic n={n} measured={measured} expected={expected}")
    if np.max(np.abs(measured-expected)) > 0.03:
        raise AssertionError("Synthetic polarization recovery failed")


if __name__ == "__main__":
    run(239, 0.0, 0.0)
    run(240, 0.3, np.deg2rad(30.0))
