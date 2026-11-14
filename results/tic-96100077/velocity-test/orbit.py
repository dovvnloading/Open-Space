"""Frozen outer orbit of the RV tracer. Input epochs are full barycentric JD(TCB)."""

import sys
from pathlib import Path

import numpy as np
from astropy.time import Time

DOSSIER = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DOSSIER / "validation/spectroscopy"))
from orbit_audit import AU_KM, PARAMETERS, read_one

EPOCH = Time(2016.0, format="jyear", scale="tcb").jd
FIELDS = (
    "c_thiele_innes",
    "h_thiele_innes",
    "center_of_mass_velocity",
    "eccentricity",
    "period",
    "t_periastron",
)
ROW = read_one("data/nss-orbit.csv")
NOMINAL = np.array([float(ROW[name]) for name in FIELDS])
# Eclipse clock scale from the merged TESS quadratic-ephemeris fit.
CLOCK_SCALE = -0.8016409589619811


def sample_orbits(n=4096, seed=96100077):
    """Use the marginal six-dimensional block of Gaia's 15-parameter covariance."""
    packed = np.fromstring(ROW["corr_vec"].strip("()"), sep=",")
    lower = np.zeros((len(PARAMETERS), len(PARAMETERS)))
    lower[np.tril_indices(len(PARAMETERS), -1)] = packed
    indices = [PARAMETERS.index(name) for name in FIELDS]
    corr = (np.eye(len(PARAMETERS)) + lower + lower.T)[np.ix_(indices, indices)]
    errors = np.array([float(ROW[f"{name}_error"]) for name in FIELDS])
    covariance = corr * np.outer(errors, errors)
    rng = np.random.default_rng(seed)
    return (
        NOMINAL
        + rng.standard_normal((n, len(FIELDS))) @ np.linalg.cholesky(covariance).T
    )


def velocity(jd_tcb, elements):
    """Differentiate Gaia's line-of-sight displacement, positive away from observer."""
    c, h, gamma, e, period, t_peri = np.asarray(elements).T
    m = ((np.asarray(jd_tcb) - EPOCH - t_peri) / period % 1) * (2 * np.pi)
    anomaly = m.copy()
    for _ in range(10):
        anomaly -= (anomaly - e * np.sin(anomaly) - m) / (1 - e * np.cos(anomaly))
    dz = (-c * np.sin(anomaly) + h * np.sqrt(1 - e**2) * np.cos(anomaly)) / (
        1 - e * np.cos(anomaly)
    )
    return gamma + AU_KM * (2 * np.pi / 86400) * dz / period


def predict_jd(tcb_jd):
    """Nominal published single-lined spectroscopic velocity in km/s, including gamma."""
    return velocity(tcb_jd, NOMINAL)
