"""Compare SAT455 to discovery astrometry and refit both line-of-sight solutions."""

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
DAY = 86400.0
C = 299792.458
EPOCH = 2460200.5
# GM (km^3/s^2) of the Saturn system, Sun and Jupiter system, respectively.
GM_SATURN = 37940584.8418
GM_SUN = 132712440041.2794
GM_JUPITER = 126712764.1
LENGTH = 1e7  # Scale the optimizer's position parameters, in km.


def table(name):
    return np.loadtxt(ROOT / "inputs" / (name + ".csv"), skiprows=1, delimiter=",", ndmin=2)


def trajectory(name):
    rows = table(name)
    return CubicHermiteSpline(rows[:, 0] - EPOCH, rows[:, 1:4], rows[:, 4:] * DAY)


def inverse_square_vector(vector):
    return vector / np.linalg.norm(vector) ** 3


def elements(state):
    pos, vel = state[:3] * LENGTH, state[3:]
    norm = np.linalg.norm(pos)
    axis = 1 / (2 / norm - (vel @ vel) / GM_SATURN)
    ecc = np.linalg.norm(np.cross(vel, np.cross(pos, vel)) / GM_SATURN - pos / norm)
    return dict(a_km=round(axis), e=round(ecc, 5), q_km=round(axis * (1 - ecc)))


def state_km_kms(state):
    return np.round(np.r_[state[:3] * LENGTH, state[3:]], 10).tolist()


def offsets(sky, observed):
    residual = (sky - observed) * 3600
    residual[:, 0] *= np.cos(np.deg2rad(observed[:, 1]))
    return residual


def rms(residual):
    return float(np.sqrt(np.mean(np.sum(residual**2, axis=1))))


class ShortArc:
    def __init__(self, observations):
        self.observed = observations[:, 2:4]
        site = table("ssb_from_t14")
        self.epochs = site[:, 0] - EPOCH
        if not np.allclose(site[:, 0], observations[:, 1] + 69.182 / DAY, rtol=0, atol=1e-8):
            raise ValueError("Observer vectors must match the observation epochs")
        self.observer = -site[:, 1:4]
        self.saturn = trajectory("saturn_barycenter")
        self.perturbers = [(trajectory("sun"), GM_SUN), (trajectory("jupiter_barycenter"), GM_JUPITER)]

    def derivative(self, t, state):
        position = state[:3] * LENGTH
        acc = -GM_SATURN * inverse_square_vector(position)
        for trajectory_fn, gm in self.perturbers:
            body = trajectory_fn(t)
            acc += gm * (inverse_square_vector(body - position) - inverse_square_vector(body))
        return np.r_[state[3:] * DAY / LENGTH, acc * DAY]

    def sky(self, state, *, rtol=3e-11):
        arcs = [solve_ivp(self.derivative, (0, end), state, method="DOP853",
                          rtol=rtol, atol=rtol * .03, dense_output=True).sol
                for end in (min(-1, self.epochs.min() - 1), max(1, self.epochs.max() + 1))]
        emit = self.epochs.copy()
        for _ in range(4):
            position = np.array([arcs[int(t > 0)](t)[:3] * LENGTH for t in emit])
            relative = self.saturn(emit) + position - self.observer
            emit = self.epochs - np.linalg.norm(relative, axis=1) / (C * DAY)
        x, y, z = relative.T
        ra = np.rad2deg(np.arctan2(y, x)) % 360
        dec = np.rad2deg(np.arctan2(z, np.hypot(x, y)))
        return np.column_stack((ra, dec))

    def residuals(self, state):
        return offsets(self.sky(state), self.observed)

    def fit(self, seed, use):
        result = least_squares(lambda state: self.residuals(state)[use].ravel(), seed,
                               x_scale="jac", diff_step=1e-6, max_nfev=110,
                               ftol=1e-10, xtol=1e-10, gtol=1e-8)
        if not result.success:
            raise RuntimeError(result.message)
        return result.x


def analyze():
    astrometry = table("astrometry")
    target = astrometry[astrometry[:, 0] == 38]
    model = ShortArc(target)
    seed = table("sat455_state")[0, 1:].copy()
    seed[:3] /= LENGTH
    public = offsets(target[:, 4:6], target[:, 2:4])
    reproduction = offsets(model.sky(seed), target[:, 4:6])
    output = {
        "target_astrometry_count": len(target),
        "arc_days": round(np.ptp(target[:, 1]), 5),
        "sat455_rms_arcsec": round(rms(public), 4),
        "model_sat455_max_sky_difference_arcsec": round(np.linalg.norm(reproduction, axis=1).max(), 6),
        "control_sat455_rms_arcsec": {},
        "solutions": [],
        "night_exclusions": [],
    }
    for moon in np.unique(astrometry[:, 0]):
        if moon != 38:
            rows = astrometry[astrometry[:, 0] == moon]
            output["control_sat455_rms_arcsec"][int(moon)] = round(rms(offsets(rows[:, 4:6], rows[:, 2:4])), 4)
    sightline = model.saturn(0) - model.observer[2]
    sightline /= np.linalg.norm(sightline)
    use = np.ones(len(target), bool)
    # Switch near/far line-of-sight position at the same measured sky direction.
    # Two distinct velocity starts at each signed depth test local convergence.
    starts = []
    for depth in (-1, 1):
        results = []
        for speed in (-1, 1):
            trial = seed.copy()
            trial[:3] += (depth - trial[:3] @ sightline) * sightline
            trial[3:] += (speed - trial[3:] @ sightline) * sightline
            solution = model.fit(trial, use)
            results.append((rms(model.residuals(solution)), solution))
        error, best = min(results, key=lambda result: result[0])
        starts.append(best)
        output["solutions"].append({
            "depth_seed_10m_km": depth, **elements(best),
            "icrf_saturn_system_state_km_kms": state_km_kms(best),
            "rms_arcsec": round(error, 4),
            "speed_seed_rms_arcsec": [round(result[0], 4) for result in results],
            "residuals_predicted_minus_observed_arcsec": np.round(model.residuals(best), 3).tolist(),
        })
    # CFHT observing nights at Hawaii span two neighboring UTC dates in August.
    night = np.floor(target[:, 1] - 0.25).astype(int)
    low_variants = []
    for date in np.unique(night):
        include = night != date
        fitted = [model.fit(state, include) for state in starts]
        low_variants.append(fitted[0])
        output["night_exclusions"].append({
            "moderate_icrf_saturn_system_state_km_kms": state_km_kms(fitted[0]),
            "jd_night": int(date), "withheld_count": int((~include).sum()),
            "withheld_rms_arcsec": [round(rms(model.residuals(state)[~include]), 3) for state in fitted],
        })
    # Forecast illustration. Re-evaluate even the public state with this model
    # to expose the effect of simplifying the force model over longer baselines.
    observer = table("prediction_observer")
    model.epochs = observer[:, 0] - EPOCH
    model.observer = -observer[:, 1:4]
    model.saturn = trajectory("prediction_saturn")
    model.perturbers = [(trajectory("prediction_sun"), GM_SUN),
                        (trajectory("prediction_jupiter"), GM_JUPITER)]
    moderate, extreme, public_check = [model.sky(state, rtol=2e-10) for state in [*starts, seed]]
    public_prediction = table("prediction_sat455")
    variants = [model.sky(state, rtol=2e-10) for state in low_variants]
    def separation_deg(first, second):
        ra1, dec1 = np.deg2rad(first).T
        ra2, dec2 = np.deg2rad(second).T
        hav = np.sin((dec2 - dec1) / 2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin((ra2 - ra1) / 2)**2
        return np.rad2deg(2 * np.arcsin(np.sqrt(np.clip(hav, 0, 1))))
    output["illustrative_predictions"] = []
    for i, jd in enumerate(public_prediction[:, 0]):
        output["illustrative_predictions"].append({
            "jd_utc": round(jd, 7),
            "moderate_ra_dec_deg": np.round(moderate[i], 5).tolist(),
            "extreme_ra_dec_deg": np.round(extreme[i], 5).tolist(),
            "official_sat455_ra_dec_deg": np.round(public_prediction[i, 1:], 5).tolist(),
            "moderate_to_extreme_arcmin": round(separation_deg(moderate, extreme)[i] * 60, 2),
            "moderate_to_sat455_arcmin": round(separation_deg(moderate, public_prediction[:, 1:])[i] * 60, 2),
            "moderate_leave_night_out_max_shift_arcmin": round(max(
                separation_deg(moderate, variant)[i] * 60 for variant in variants), 2),
            "public_state_reproduction_arcsec": round(
                separation_deg(public_check, public_prediction[:, 1:])[i] * 3600, 3),
        })
    return output


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2))
