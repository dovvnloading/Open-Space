"""Reconstruct Judith's along-scan residuals and the two ATLAS model predictions."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
AU_KM = 149597870.7
LIGHT_DAYS_PER_AU = 0.005775518331
MAS_PER_RAD = 180 / np.pi * 3600000
MODELS = ("5664", "5665")


def vectors(filename):
    response = json.loads((DATA / filename).read_text())["result"]
    lines = response.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    records = [line.split(",") for line in lines]
    return np.array([[float(r[0]), *map(float, r[2:8])] for r in records])


def tdb(epoch):
    jd = epoch + 2455197.5
    return jd - (jd - 2443144.5003725) * 1.550519768e-8 - 6.55e-5 / 86400


def unit(vectors_):
    length = np.linalg.norm(vectors_, axis=1)
    return vectors_ / length[:, None], length


def geometry(rows, states, indices, sun_states):
    """Retard each target state to emission time using the receiving Gaia position."""
    t = tdb(rows.epoch.to_numpy())
    origin = rows[["x_gaia", "y_gaia", "z_gaia"]].to_numpy()
    base = states[indices]
    pos, vel = base[:, 1:4], base[:, 4:7]
    solar_position = np.column_stack([
        np.interp(t, sun_states[:, 0], sun_states[:, k]) for k in range(1, 4)
    ])
    solar_direction, radius = unit(solar_position - pos)
    acceleration = 0.0002959122083 * solar_direction / radius[:, None] ** 2
    emission = t.copy()
    for _ in range(4):
        dt = (emission - base[:, 0])[:, None]
        position = pos + dt * vel + 0.5 * dt**2 * acceleration
        view, distance = unit(origin - position)
        emission = t - LIGHT_DAYS_PER_AU * distance
    sun, heliocentric = unit(solar_position - position)
    ra, dec, angle = np.deg2rad(rows[["ra", "dec", "position_angle_scan"]].to_numpy().T)
    east = np.column_stack([-np.sin(ra), np.cos(ra), np.zeros(len(ra))])
    north = np.column_stack([-np.cos(ra) * np.sin(dec), -np.sin(ra) * np.sin(dec), np.cos(dec)])
    scan = east * np.sin(angle)[:, None] + north * np.cos(angle)[:, None]
    return view, sun, scan, emission, distance, heliocentric


def observed_transits():
    ccd = pd.read_csv(DATA / "gaia_fpr.csv")
    groups = ccd.groupby("transit_id")
    visits = groups.mean(numeric_only=True)
    states = np.vstack([vectors("horizons-0.json"), vectors("horizons-50.json")])
    if not np.allclose(states[:, 0], tdb(visits.epoch.to_numpy()), atol=3e-8, rtol=0):
        raise ValueError("Horizons epochs must match Gaia transit means")
    solar = vectors("sun.json")
    geometries = geometry(visits, states, np.arange(len(visits)), solar)
    accepted = ccd.loc[~ccd.is_rejected].copy()
    view = geometry(accepted, states, visits.index.get_indexer(accepted.transit_id), solar)[0]
    ra, dec, scan = np.deg2rad(accepted[["ra", "dec", "position_angle_scan"]].to_numpy().T)
    predicted_ra = np.arctan2(-view[:, 1], -view[:, 0])
    predicted_dec = np.arcsin(-view[:, 2])
    east = ((ra - predicted_ra + np.pi) % (2 * np.pi) - np.pi) * np.cos(dec)
    north = dec - predicted_dec
    accepted["al"] = MAS_PER_RAD * (east * np.sin(scan) + north * np.cos(scan))
    visits["fpr_oc"] = accepted.groupby("transit_id").al.median()
    visits["accepted_ccds"] = accepted.groupby("transit_id").size().reindex(visits.index, fill_value=0)
    res2 = pd.read_csv(DATA / "gaia_dr2.csv")
    visits["dr2_oc"] = res2.loc[res2.selected].groupby("transit_id").residual_al.median()
    phot = pd.read_csv(DATA / "gaia_dr3_photometry.csv")
    visits["reduced_g"] = phot.groupby("transit_id").g_mag.median().reindex(visits.index) - 5 * np.log10(
        geometries[4] * geometries[5]
    )
    return visits, geometries


def facets(model, diameter):
    lines = (DATA / f"{model}-shape.txt").read_text().splitlines()
    count, _ = map(int, lines[0].split())
    vertices = np.loadtxt(lines[1:count + 1])
    triangles = vertices[np.loadtxt(lines[count + 1:], dtype=int) - 1]
    centroids = triangles.mean(axis=1)
    area_vectors = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]) / 2
    volumes = np.sum(area_vectors * centroids, axis=1) / 3
    volume = volumes.sum()
    center = (centroids * volumes[:, None]).sum(axis=0) * 0.75 / volume
    scale = diameter / (2 * (3 * abs(volume) / (4 * np.pi))**(1 / 3))
    normals, areas = unit(area_vectors * np.sign(volume))
    return (centroids - center) * scale, normals, areas


def illuminate(model, geom, diameter, law, phase_shift=0):
    alpha, delta, rate, epoch, phase = np.fromstring(
        (DATA / f"{model}-IAUspin").read_text(), sep=" "
    )
    a, d = np.deg2rad([alpha, delta])
    x_axis = np.array([-np.sin(a), np.cos(a), 0])
    pole = np.array([np.cos(a) * np.cos(d), np.sin(a) * np.cos(d), np.sin(d)])
    y_axis = np.cross(pole, x_axis)
    view, sun, scan, emission, distance, _ = geom
    theta = np.deg2rad((phase + phase_shift + rate * (emission - epoch)) % 360)
    c, s = np.cos(theta)[:, None], np.sin(theta)[:, None]
    body_x, body_y = c * x_axis + s * y_axis, -s * x_axis + c * y_axis

    def body(vector):
        return np.column_stack([(vector * body_x).sum(axis=1), (vector * body_y).sum(axis=1), vector @ pole])

    centers, normals, areas = facets(model, diameter)
    mu = np.maximum(body(view) @ normals.T, 0)
    mu0 = np.maximum(body(sun) @ normals.T, 0)
    weight = mu * mu0 * areas
    if law != "Lambert":
        weight *= 1 / np.maximum(mu + mu0, 1e-15) + (0.1 if law == "mixed" else 0)
    flux = weight.sum(axis=1)
    pc = (weight @ centers) / flux[:, None]
    offset = np.sum(pc * body(scan), axis=1) * MAS_PER_RAD / (distance * AU_KM)
    return -2.5 * np.log10(flux), offset


def detrend(values, design):
    return values - design @ np.linalg.lstsq(design, values, rcond=None)[0]


def rms(y):
    return float(np.sqrt(np.mean(y**2)))


def statistics(y, p):
    scale = np.dot(p, y) / np.dot(p, p)
    return {
        "n": len(y), "correlation": float(np.corrcoef(y, p)[0, 1]),
        "observed_rms": rms(y), "predicted_rms": rms(p), "fixed_residual_rms": rms(y - p),
        "fitted_scale": float(scale), "scaled_residual_rms": rms(y - scale * p)
    }


def clustered(visits, observed, predicted):
    sample = visits[["epoch", observed, predicted]].dropna().copy()
    # A mean per group removes geometry and orbital offsets on timescales longer than a visit sequence.
    group = (sample.epoch.diff() > 5).cumsum()
    sample["group"] = group
    sample = sample.loc[sample.groupby("group").epoch.transform("size") >= 2]
    cols = [observed, predicted]
    residual = sample[cols] - sample.groupby("group")[cols].transform("mean")
    return sample, residual[observed].to_numpy(), residual[predicted].to_numpy()


def main():
    for path, receipt in json.loads((DATA / "sources.json").read_text()).items():
        if hashlib.sha256((DATA / path).read_bytes()).hexdigest() != receipt["sha256"]:
            raise ValueError(f"Input hash mismatch: {path}")
    visits, geom = observed_transits()
    candidate = next(csv.DictReader((DATA / "gaia_candidate.csv").read_text().splitlines()))
    diameter = float(candidate["Diameter(km)"])
    lo, hi = float(candidate["time0"]), float(candidate["time1"])
    selected = visits.epoch.between(lo - 0.002, hi + 0.002) & visits.fpr_oc.notna()
    t = visits.epoch.to_numpy()[selected]
    x = t - t.mean()
    designs = {
        "constant": np.ones((len(t), 1)),
        "linear": np.column_stack([np.ones(len(t)), x]),
        "quadratic": np.column_stack([np.ones(len(t)), x, x*x])
    }
    design = designs["quadratic"]
    wave = 2 * np.pi * x * 24 / float(candidate["Wobble_Period(h)"])
    periodic_design = np.column_stack([design, np.cos(wave), np.sin(wave)])
    def amplitude(values):
        return float(np.linalg.norm(np.linalg.lstsq(periodic_design, values, rcond=None)[0][-2:]))
    observed = visits.fpr_oc.to_numpy()[selected]
    earlier = visits.dr2_oc.notna()
    results = {
        "keep_fpr_ccds": int(visits.accepted_ccds[selected].sum()),
        "report_transits": int(selected.sum()), "dr2_total": int(earlier.sum()),
        "shared_transit_ids": int((earlier & selected).sum()),
        "earlier_to_report_gap_days": float(t.min() - visits.epoch[earlier].max()),
        "report_observed_semiamplitude_mas": amplitude(observed), "model": {}
    }
    alt_diameter = float(next(p["value"] for p in json.loads((DATA / "jpl_sbdb.json").read_text())["phys_par"] if p["name"] == "diameter"))
    visits["catalog_window"] = selected
    for model in MODELS:
        mag, pc = illuminate(model, geom, diameter, "mixed")
        visits[f"mag_{model}"] = mag
        visits[f"al_{model}"] = pc
        sampled = pc[selected]
        stats = statistics(detrend(observed, design), detrend(sampled, design))
        old, old_y, old_p = clustered(visits, "dr2_oc", f"al_{model}")
        old_stats = statistics(old_y, old_p)
        phot, phot_y, phot_p = clustered(visits, "reduced_g", f"mag_{model}")
        model_results = {
            "catalog_window": stats, "earlier_dr2": old_stats, "earlier_dr2_groups": int(old.group.nunique()),
            "gaia_photometry": statistics(phot_y, phot_p),
            "semiamplitude_mas": amplitude(sampled),
            "residual_semiamplitude_fixed_mas": amplitude(observed - sampled),
            "residual_semiamplitude_scaled_mas": amplitude(observed - stats["fitted_scale"] * sampled),
            "earlier_scaled_from_late_rms_mas": rms(old_y - stats["fitted_scale"] * old_p),
            "late_scaled_from_earlier_rms_mas": rms(detrend(observed - old_stats["fitted_scale"] * sampled, design)),
            "alternative_diameter_km": alt_diameter,
            "alternative_diameter_residual_rms_mas": rms(detrend(observed - sampled * alt_diameter / diameter, design)),
            "sensitivity": {}
        }
        for law in ["LS", "mixed", "Lambert"]:
            _, offsets = illuminate(model, geom, diameter, law)
            model_results["sensitivity"][law] = {
                name: statistics(detrend(observed, d), detrend(offsets[selected], d))
                for name, d in designs.items()
            }
        results["model"][model] = model_results
    evidence = ROOT / "measurements"
    evidence.mkdir(exist_ok=True)
    visits[["epoch", "accepted_ccds", "catalog_window", "fpr_oc", "dr2_oc", "reduced_g",
            "al_5664", "al_5665", "mag_5664", "mag_5665"]].to_csv(evidence / "transits.csv", float_format="%.10f")
    (evidence / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    return results


if __name__ == "__main__":
    report = main()
    print(json.dumps({m: report["model"][m]["catalog_window"] for m in MODELS}, indent=2))
