"""Construct preliminary heliocentric orbits from training sky directions."""

import numpy as np

from skyfit import EPOCH, LIGHT_AU_DAY


def kepler_drift(position, velocity, elapsed, mu):
    """Propagate an array of bound two-body states with elliptic f and g."""
    radius = np.linalg.norm(position, axis=-1)
    semimajor = 1 / (2 / radius - np.sum(velocity**2, axis=-1) / mu)
    motion = np.sqrt(mu / semimajor**3)
    # These are e cos(E0) and e sin(E0) at the initial epoch.
    e_cos = 1 - radius / semimajor
    e_sin = np.sum(position * velocity, axis=-1) / np.sqrt(mu * semimajor)
    mean = elapsed * motion
    mean = mean - 2 * np.pi * np.rint(mean / (2 * np.pi))
    anomaly = mean.copy()
    for _ in range(16):
        cosine, sine = np.cos(anomaly), np.sin(anomaly)
        residual = anomaly - e_cos * sine + e_sin * (1 - cosine) - mean
        derivative = 1 - e_cos * cosine + e_sin * sine
        anomaly -= np.clip(residual / derivative, -1.0, 1.0)

    cosine, sine = np.cos(anomaly), np.sin(anomaly)
    residual = anomaly - e_cos * sine + e_sin * (1 - cosine) - mean
    # A non-converged trial is excluded from the grid search.
    anomaly = np.where(np.abs(residual) < 1e-10, anomaly, np.nan)
    cosine, sine = np.cos(anomaly), np.sin(anomaly)
    f = 1 - semimajor / radius * (1 - cosine)
    g = (mean + sine - anomaly) / motion
    new_position = f[..., None] * position + g[..., None] * velocity
    new_radius = np.linalg.norm(new_position, axis=-1)
    f_dot = -np.sqrt(mu * semimajor) * sine / (radius * new_radius)
    g_dot = 1 - semimajor / new_radius * (1 - cosine)
    new_velocity = f_dot[..., None] * position + g_dot[..., None] * velocity
    return new_position, new_velocity


def training_seed(model, training_dates):
    """Estimate a starting state using only the supplied model and its UTC dates."""
    dates = np.unique(training_dates)
    recent = np.isin(training_dates, dates[-2:])
    epoch = model.t[recent].mean()
    times = model.t[recent] - epoch
    sightline = model.sightline[recent].mean(axis=0)
    sightline /= np.linalg.norm(sightline)
    # Linear slopes across the latest two WFST nights in the training set.
    weights = times / (times @ times)
    angular_velocity = weights @ model.sightline[recent]
    observer = model.observers[recent].mean(axis=0)
    observer_velocity = weights @ model.observers[recent]

    distances, speeds = np.meshgrid(
        np.linspace(0.5, 9, 300), np.linspace(-0.015, 0.015, 211)
    )
    distances, speeds = distances.ravel(), speeds.ravel()
    sun, sun_velocity = model.gravity.state(0, 10, EPOCH + epoch)
    position = observer - sun + distances[:, None] * sightline
    velocity = (
        observer_velocity - sun_velocity
        + distances[:, None] * angular_velocity
        + speeds[:, None] * sightline
    )
    mu = model.gravity.gm[0]
    # Restrict the grid to bound orbits below sqrt(1.93 / 2) of escape speed.
    bound = np.sum(velocity**2, axis=1) < 1.93 * mu / np.linalg.norm(position, axis=1)
    position, velocity = position[bound], velocity[bound]
    distances, speeds = distances[bound], speeds[bound]

    score = np.zeros(len(distances))
    for date in dates[:-2]:
        i = np.flatnonzero(training_dates == date)[0]
        predicted, _ = kepler_drift(position, velocity, model.t[i] - epoch, mu)
        sun_then, _ = model.gravity.state(0, 10, EPOCH + model.t[i])
        relative = predicted + sun_then - model.observers[i]
        cosine = (relative @ model.sightline[i]) / np.linalg.norm(relative, axis=1)
        score += 1 - np.clip(cosine, -1, 1)
    if not np.isfinite(score).any():
        raise RuntimeError("No finite preliminary orbit on the training grid.")
    i = np.nanargmin(score)
    # The observed direction locates the body at emission, before reception.
    elapsed = -epoch + distances[i] / LIGHT_AU_DAY
    initial = kepler_drift(position[i], velocity[i], elapsed, mu)
    state = (np.array(initial) + model.gravity.state(0, 10, EPOCH)).ravel()
    diagnostic = {
        "source": "training sky positions and DE440s",
        "training_nights_utc": dates.tolist(),
        "initial_range_au": distances[i],
        "initial_range_rate_au_day": speeds[i],
    }
    return state, diagnostic
