"""Point-source point-lens microlensing models with annual parallax.

Conventions
-----------
The sign conventions of u0 and the parallax vector match MulensModel 3.12.
Times are HJD - 2450000. Fluxes use a zero point of I = 18 mag.
The parallax model is geocentric (Gould 2004): the Sun's position relative to
Earth, projected on the sky at the event, has its value and velocity at a
fixed reference time `t0par` removed. `pi_n` and `pi_e` are the north and east
components of the microlens parallax vector in that frame.
"""
import numpy as np
import erfa
from scipy.optimize import least_squares

ZP = 18.0


def read_phot(path):
    data = np.loadtxt(path, usecols=(0, 1, 2, 3, 4), ndmin=2)
    t, mag, err = data[:, 0] - 2450000.0, data[:, 1], data[:, 2]
    good = (err > 0) & (err < 1.0) & (mag > 5) & (mag < 25) & np.isfinite(mag)
    return t[good], mag[good], err[good], data[good, 3], data[good, 4]


def to_flux(mag, err, floor=0.004):
    flux = 10 ** (-0.4 * (mag - ZP))
    sigma = flux * 0.4 * np.log(10) * np.hypot(err, floor)
    return flux, sigma


def sexagesimal(ra, dec):
    h, m, s = (float(x) for x in ra.split(":"))
    sign = -1.0 if dec.strip().startswith("-") else 1.0
    d, dm, ds = (abs(float(x)) for x in dec.split(":"))
    return np.radians(15 * (h + m / 60 + s / 3600)), np.radians(sign * (d + dm / 60 + ds / 3600))


def sun_projected(t, ra, dec):
    """North and east components (au) of the Sun's position seen from Earth,
    with velocity components (au/day). `t` is HJD - 2450000."""
    jd = np.atleast_1d(t) + 2450000.0
    pvh, _ = erfa.epv00(jd, np.zeros_like(jd))
    position, velocity = -pvh["p"], -pvh["v"]
    east = np.array([-np.sin(ra), np.cos(ra), 0.0])
    north = np.array([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)])
    return position @ north, position @ east, velocity @ north, velocity @ east


class Parallax:
    """Pre-computed offsets of the Sun for one event and one reference time."""

    def __init__(self, t, ra, dec, t0par):
        sn, se, _, _ = sun_projected(t, ra, dec)
        sn0, se0, vn0, ve0 = sun_projected(np.array([t0par]), ra, dec)
        self.dn = sn - sn0[0] - (t - t0par) * vn0[0]
        self.de = se - se0[0] - (t - t0par) * ve0[0]

    def shift(self, pi_n, pi_e):
        return pi_n * self.dn + pi_e * self.de, pi_n * self.de - pi_e * self.dn


def magnification(u2):
    u2 = np.maximum(u2, 1e-12)
    return (u2 + 2) / np.sqrt(u2 * (u2 + 4))


def trajectory(theta, t, par=None):
    t0, u0, te = theta[:3]
    tau = (t - t0) / te
    beta = np.full_like(t, u0)
    if par is not None:
        dtau, dbeta = par.shift(theta[3], theta[4])
        tau = tau + dtau
        beta = beta + dbeta
    return tau, beta


def linear_fluxes(amp, flux, sigma):
    """Weighted least-squares source and blend flux for a magnification curve."""
    w = 1.0 / sigma ** 2
    s_aa, s_a, s_1 = np.sum(w * amp * amp), np.sum(w * amp), np.sum(w)
    s_af, s_f = np.sum(w * amp * flux), np.sum(w * flux)
    det = s_aa * s_1 - s_a * s_a
    if det <= 0:
        return 0.0, s_f / s_1
    return (s_af * s_1 - s_a * s_f) / det, (s_aa * s_f - s_a * s_af) / det


def model(theta, t, flux, sigma, par=None):
    tau, beta = trajectory(theta, t, par)
    amp = magnification(tau * tau + beta * beta)
    fs, fb = linear_fluxes(amp, flux, sigma)
    return fs * amp + fb, fs, fb


BLEND_FLOOR = -0.158     # blend flux of a "negative" I = 20 star; fainter blends are penalised
BLEND_SOFT = 0.03


def residuals(theta, t, flux, sigma, par=None):
    """Scaled residuals plus one term that rejects strongly negative blending."""
    if theta[2] <= 0:
        return np.full(len(t) + 1, 1e6)
    curve, _, fb = model(theta, t, flux, sigma, par)
    penalty = min(0.0, fb - BLEND_FLOOR) / BLEND_SOFT
    return np.append((flux - curve) / sigma, penalty)


def fit(theta0, t, flux, sigma, par=None):
    scale = np.array([1.0, 0.05, 5.0, 0.05, 0.05])[: len(theta0)]
    result = least_squares(residuals, theta0, args=(t, flux, sigma, par), x_scale=scale,
                           method="lm", max_nfev=600)
    chi2 = float(np.sum(result.fun ** 2))
    try:
        cov = np.linalg.inv(result.jac.T @ result.jac)
        err = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        cov, err = None, np.full(len(theta0), np.nan)
    return result.x, err, chi2, cov
