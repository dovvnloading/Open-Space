"""Fit barycentric Cartesian orbits to the archived MPC sky positions."""

import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.utils import iers
from jplephem.spk import SPK
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import least_squares

DATA = Path(__file__).resolve().parent / "data"
AU_KM = 149597870.7
DAY = 86400.0
EPOCH = 2461200.5
ARCSEC = np.rad2deg(1.0) * 3600
LIGHT_AU_DAY = 173.144632674
# Sun then planetary systems Mercury through Neptune, km^3/s^2.
GM = np.array(
    [
        132712440041.27942,
        22031.86855,
        324858.592,
        403503.235625,
        42828.375816,
        126712764.1,
        37940584.8418,
        5794556.4,
        6836527.10058,
    ]
)
IDS = [10, 1, 2, 3, 4, 5, 6, 7, 8]
angle = np.deg2rad(23.4392911111)
EQUATOR_TO_ECLIPTIC = np.array(
    [[1, 0, 0], [0, np.cos(angle), np.sin(angle)], [0, -np.sin(angle), np.cos(angle)]]
)
iers.conf.auto_download = False


class Gravity:
    def __init__(self, kernel_path):
        archive = json.loads((DATA / "sources.json").read_text())
        if (
            hashlib.sha256(kernel_path.read_bytes()).hexdigest()
            != archive["DE440s"]["sha256"]
        ):
            raise ValueError("Kernel hash differs from the frozen DE440s input.")
        self.kernel = SPK.open(str(kernel_path))
        dates = np.arange(2448500.0, 2461803.0, 3.0)
        planets = np.array([self.state(0, body, dates) for body in IDS])
        earth_moon = planets[3].copy()
        planets[3] = earth_moon + self.state(3, 399, dates)
        planets = np.concatenate(
            (planets, (earth_moon + self.state(3, 301, dates))[None, :])
        )
        self.gm = np.r_[GM, 4902.800118]
        self.gm[3] -= self.gm[-1]
        self.gm *= DAY**2 / AU_KM**3
        self.planet_positions = CubicHermiteSpline(
            dates - EPOCH,
            planets[:, 0].transpose(1, 0, 2),
            planets[:, 1].transpose(1, 0, 2),
        )

    def state(self, center, body, times):
        position, velocity = self.kernel[center, body].compute_and_differentiate(times)
        return np.array([position.T, velocity.T]) / AU_KM

    def acceleration(self, time, augmented_state):
        state = augmented_state[:6]
        transition = augmented_state[6:].reshape(6, 6)
        offset = self.planet_positions(time) - state[:3]
        distance_sq = np.sum(offset**2, axis=1)
        coefficients = self.gm / distance_sq**1.5
        acceleration = coefficients @ offset
        gradient = (offset.T * (3 * coefficients / distance_sq)) @ offset - np.eye(
            3
        ) * coefficients.sum()
        return np.r_[
            state[3:],
            acceleration,
            transition[3:].ravel(),
            (gradient @ transition[:3]).ravel(),
        ]


def observer_coordinates(rows, gravity):
    times = Time(rows.utc.str.rstrip("Z").tolist(), format="isot", scale="utc")
    jd = times.tdb.jd
    positions = gravity.state(0, 3, jd)[0] + gravity.state(3, 399, jd)[0]
    sites = {
        line[:3]: (float(line[4:13]), float(line[13:21]), float(line[21:30]))
        for line in (DATA / "stations.txt").read_text().splitlines()
    }
    for station, indices in rows.groupby("site").indices.items():
        longitude, cosine, sine = sites[station]
        phi = np.deg2rad(longitude)
        site = EarthLocation.from_geocentric(
            6378.137 * cosine * np.cos(phi),
            6378.137 * cosine * np.sin(phi),
            6378.137 * sine,
            unit="km",
        )
        positions[indices] += (
            site.get_gcrs_posvel(times[indices])[0].xyz.to_value("au").T
        )
    return jd - EPOCH, positions


class SkyFit:
    def __init__(self, rows, gravity):
        self.t, self.observers = observer_coordinates(rows, gravity)
        self.gravity = gravity
        ra, dec = np.deg2rad(rows[["ra", "dec"]].values).T
        self.sightline = np.column_stack(
            (np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec))
        )
        east = np.column_stack((-np.sin(ra), np.cos(ra), np.zeros(len(ra))))
        self.axes = np.stack((east, np.cross(self.sightline, east)), axis=1)
        self.last = None

    def prediction(self, state):
        if self.last is not None and np.array_equal(self.last, state):
            return self.residual, self.jacobian
        initial = np.r_[state, np.eye(6).ravel()]
        limits = [min(self.t.min() - 0.06, -0.1), max(self.t.max() + 0.01, 0.1)]
        paths = [
            solve_ivp(
                self.gravity.acceleration,
                (0, limit),
                initial,
                method="DOP853",
                dense_output=True,
                atol=5e-13,
                rtol=5e-12,
            )
            for limit in limits
        ]
        if not all(path.success for path in paths):
            raise RuntimeError("Orbit integration failed.")

        def sample(times):
            result = np.empty((len(times), 42))
            for path, mask in zip(paths, [times <= 0, times > 0]):
                if mask.any():
                    result[mask] = path.sol(times[mask]).T
            return result

        emission = self.t.copy()
        for _ in range(3):
            pos = sample(emission)[:, :3] - self.observers
            emission = self.t - np.linalg.norm(pos, axis=1) / LIGHT_AU_DAY
        states = sample(emission)
        pos = states[:, :3] - self.observers
        distance = np.linalg.norm(pos, axis=1)
        direction = pos / distance[:, None]
        residual = (
            np.einsum("nai,ni->na", self.axes, direction - self.sightline) * ARCSEC
        )
        derivatives = states[:, 6:].reshape(-1, 6, 6)[:, :3, :]
        velocity = states[:, 3:6]
        radial = np.einsum("ni,nij->nj", direction, derivatives)
        derivatives = (
            derivatives
            - velocity[:, :, None]
            * radial[:, None, :]
            / (LIGHT_AU_DAY + np.sum(direction * velocity, axis=1))[:, None, None]
        )
        projected = (
            self.axes
            - np.einsum("nai,ni->na", self.axes, direction)[:, :, None]
            * direction[:, None, :]
        )
        jacobian = (
            np.einsum("nai,nij->naj", projected, derivatives)
            * (ARCSEC / distance)[:, None, None]
        )
        self.last, self.residual, self.jacobian = (
            state.copy(),
            residual.ravel(),
            jacobian.reshape(-1, 6),
        )
        return self.residual, self.jacobian

    def fit(self, initial):
        result = least_squares(
            lambda x: self.prediction(x)[0],
            initial,
            jac=lambda x: self.prediction(x)[1],
            x_scale="jac",
            loss="soft_l1",
            f_scale=0.6,
            max_nfev=80,
            ftol=1e-8,
            gtol=1e-8,
            xtol=2e-13,
        )
        if not result.success:
            raise RuntimeError(result.message)
        return result.x, self.prediction(result.x)
