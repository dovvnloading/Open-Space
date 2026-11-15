"""Replay public MOA, KMTNet and Kepler photometry. See README for units and provenance."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from scipy.linalg import svd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
RA, DEC = 271.27245833, -27.82268
OBS_WINDOWS = [(4330, 4350), (4557, 4577), (5317, 5337),
               (6420, 6441), (6920, 6941)]
# These windows enclose candidates found by the full nine-year scan below.
BASE = 2450000


def write_table(path, header, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def fetch_data(cache):
    manifest = list(csv.DictReader((HERE / "inputs.csv").open()))
    cache.mkdir(parents=True, exist_ok=True)
    for row in manifest:
        path = cache / row["file"]
        if not path.exists():
            with urlopen(row["url"], timeout=180) as response, path.with_suffix(".part").open("wb") as out:
                shutil.copyfileobj(response, out)
            path.with_suffix(".part").replace(path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            raise ValueError(f"Snapshot differs for {path}")


def night_groups(time):
    night = np.floor(time - 0.65)
    for value in np.unique(night):
        yield night == value


def moa_nights(path):
    data = np.loadtxt(path, skiprows=4, usecols=(0, 1, 2, 3, 6))
    # cor_flux is NaN for excluded observations. Both flux columns are retained.
    data = data[np.isfinite(data[:, 2]) & (data[:, 3] < 800) & (data[:, 4] < 3)]
    rows = []
    for mask in night_groups(data[:, 0]):
        t, raw, flux, error, _seeing = data[mask].T
        center = np.median(flux)
        scatter = 1.4826 * np.median(abs(flux - center))
        err = max(15., max(scatter, np.median(error)) / np.sqrt(len(t)))
        rows.append([np.mean(t), center, err, len(t), np.median(raw)])
    if not rows:
        return np.empty((0, 5)), 0
    rows = np.array(rows)
    rows[:, 1] -= np.median(rows[:, 1])
    rows[:, 4] -= np.median(rows[:, 4])
    return rows[rows[:, 3] >= 3], len(data)


def moa_results(cache, out):
    rows, count = moa_nights(cache / "moa13data")
    scatter = 1.4826 * np.median(abs(rows[:, 1] - np.median(rows[:, 1])))
    score = rows[:, 1] / np.hypot(rows[:, 2], scatter)
    strong = rows[score >= 6]
    groups = np.split(strong, np.where(np.diff(strong[:, 0]) > 21)[0] + 1)
    details = []
    for group in groups:
        if not len(group):
            continue
        peak = group[np.argmax(group[:, 1])]
        details.append([float(peak[0]), int(len(group))])
    results = []
    with (HERE / "inputs.csv").open() as manifest:
        names = [r["file"] for r in csv.DictReader(manifest) if r["file"].startswith("control_")]
    controls = {name: moa_nights(cache / name)[0] for name in names}
    for lo, hi in OBS_WINDOWS:
        interval = rows[(rows[:, 0] >= lo) & (rows[:, 0] <= hi)]
        peak = interval[np.argmax(interval[:, 1])]
        active = interval[interval[:, 1] >= 1000]
        ctrl_scores = []
        for series in controls.values():
            if not len(series):
                continue
            width = 1.4826 * np.median(abs(series[:, 1] - np.median(series[:, 1])))
            z = series[:, 1] / np.hypot(series[:, 2], width)
            matched = [np.min(abs(active[:, 0] - epoch)) < 0.15 for epoch in series[:, 0]]
            if np.any(matched):
                ctrl_scores.append(max(z[matched]))
        results.append([round(peak[0], 5), round(peak[1], 1), round(peak[4], 1),
                        round(peak[1] / np.hypot(peak[2], scatter), 2),
                        len(active), int(active[:, 3].sum()), len(ctrl_scores),
                        round(max(ctrl_scores), 2)])
    write_table(out / "moa-episodes.csv",
                ["HJD_minus_2450000_at_brightest_sampled_night", "corrected_flux",
                 "raw_flux", "robust_peak_score", "nights_above_1000",
                 "exposures_on_those_nights", "controls_with_data", "max_control_score"],
                results)
    write_table(out / "moa-nights.csv", ["HJD_minus_2450000", "corrected_flux", "error",
                                        "exposures", "raw_flux"],
                [[f"{t:.5f}", f"{y:.2f}", f"{e:.2f}", int(n), f"{raw:.2f}"]
                 for t, y, e, n, raw in rows])
    return {"good_exposures": count, "nights_with_three_exposures": len(rows),
            "full_scan_groups_peak_and_nights_above_6": details,
            "night_scatter": round(scatter, 2)}


def kmt_series(cache):
    result = []
    with tarfile.open(cache / "kmtobs") as archive:
        for member in archive:
            if not member.name.endswith("_I.pysis"):
                continue
            data = np.loadtxt(archive.extractfile(member))
            data = data[(data[:, 5] < 8) & (data[:, 2] > 0) & (data[:, 2] < 2000)]
            baseline = data[(data[:, 0] < 7525) | (data[:, 0] > 7548)]
            # Brightening has negative difference flux in pySIS; use the event
            # median to normalize each overlap field and instrument separately.
            base = np.median(baseline[:, 1])
            scale = base - np.median(data[(data[:, 0] > 7533) & (data[:, 0] < 7539), 1])
            for t, flux in data[:, :2]:
                if 7530 < t < 7552:
                    result.append((member.name, t, (base - flux) / scale))
    return result


def k2_results(cache, out):
    with fits.open(cache / "k2.fits.gz") as hdus:
        table = hdus[1].data
        times = table["TIME"].copy() + 4833
        cadence = table["CADENCENO"].copy()
        clean = table["QUALITY"] == 0
        images = np.array(table["FLUX"], dtype=float)
        origin = np.array([hdus[1].header["1CRV5P"], hdus[1].header["2CRV5P"]])
    grid = np.genfromtxt(cache / "g92_52")
    basis = [1, RA, DEC, RA**2, RA*DEC, DEC**2]
    pos = np.full((len(times), 2), np.nan)
    idx = np.searchsorted(cadence, grid[:, 0])
    hit = idx < len(cadence)
    hit[hit] &= cadence[idx[hit]] == grid[hit, 0]
    pos[idx[hit]] = np.column_stack((grid[hit, 4:10] @ basis,
                                      grid[hit, 10:16] @ basis)) - origin
    use = clean & np.isfinite(pos).all(axis=1)
    t, pos, images = times[use], pos[use], images[use]
    x, y = np.meshgrid(np.arange(images.shape[2]), np.arange(images.shape[1]))
    distance2 = (x - 42)**2 + (y - 1.5)**2
    reference = (distance2 > 8**2) & (x >= 2) & (x <= 47) & (y >= 2) & (y <= 55)
    ref = images[:, reference]
    ref = ref[:, np.isfinite(ref).all(axis=0)]
    ref -= np.median(ref, axis=0)
    ref /= np.maximum(np.std(ref, axis=0), 1)
    components = svd(ref, full_matrices=False, check_finite=False)[0]
    inside = distance2 < 8**2
    pix = images[:, inside]
    available = np.isfinite(pix).all(axis=0)
    pix = pix[:, available]
    xx, yy = x[inside][available], y[inside][available]
    weights = np.exp(-((pos[:, 0, None] - xx)**2 + (pos[:, 1, None] - yy)**2) / (2 * .85**2))
    # Protect the full eruption. The last 4.4 days remain an unused check.
    train = ((t < 7532) | (t > 7544)) & (t < 7568)
    val = t >= 7568
    event = (t > 7533) & (t < 7541)
    checks = []
    for size in (3, 5, 10, 20, 40, 80, 150):
        design = np.column_stack((np.ones(len(t)), pos - [42, 1.5], components[:, :size]))
        coeff = np.linalg.lstsq(design[train], pix[train], rcond=None)[0]
        residual = pix - design @ coeff
        f = np.sum(residual * weights, axis=1)
        checks.append([size, round(np.median(f[event]), 2),
                       round(np.std(f[val]), 2)])
        if size == 40:
            main = f.copy()
            diff = np.median(residual[event], axis=0)
            largest = np.argmax(diff)
            pixel = (int(xx[largest] + origin[0]), int(yy[largest] + origin[1]))
            separation = np.linalg.norm(np.array(pixel) - np.median(pos[event], axis=0) - origin)
            write_table(out / "k2-pixels.csv", ["detector_column", "detector_row",
                        "median_event_residual_e_per_s"],
                        [[int(u+origin[0]), int(v+origin[1]), f"{z:.3f}"]
                         for u, v, z in zip(xx, yy, diff)])
    kmt = kmt_series(cache)
    daily = []
    for day in range(7531, 7573):
        use = np.floor(t) == day
        if not use.any():
            continue
        tbin, fbin = np.median(t[use]), np.median(main[use])
        nearby = {}
        for series, time, flux in kmt:
            if int(time) == day:
                nearby.setdefault(series, []).append(flux)
        fields = [np.median(a) for a in nearby.values() if len(a) >= 3]
        ground = np.median(fields) if fields else np.nan
        daily.append([day, tbin, fbin, ground, int(use.sum()), len(fields)])
    common = np.array([r for r in daily if r[0] <= 7545 and np.isfinite(r[3])])
    corr = spearmanr(common[:, 2], common[:, 3]).statistic
    write_table(out / "daily.csv", ["integer_JD_minus_2450000", "median_K2_BJD_minus_2450000",
                "K2_weighted_residual_e_per_s", "KMT_normalized_ground_flux", "K2_cadences",
                "KMT_field_site_series"], [[d, f"{tbin:.5f}", f"{flux:.3f}",
                f"{ground:.4f}" if np.isfinite(ground) else "", n, ng]
                for d, tbin, flux, ground, n, ng in daily])
    return {"good_spacecraft_cadences": len(t), "reference_pixels": ref.shape[1],
            "event_brightest_detector_pixel": pixel, "predicted_median_detector_position":
            np.round(np.median(pos[event], axis=0) + origin, 3).tolist(),
            "brightest_pixel_offset_pix": round(separation, 3),
            "K2_KMT_daily_spearman": round(corr, 4), "paired_days": len(common),
            "PC_count_event_median_and_heldout_cadence_rms": checks}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=Path("cache"))
    parser.add_argument("--out", type=Path, default=Path("replay-output"))
    args = parser.parse_args()
    fetch_data(args.cache)
    args.out.mkdir(parents=True, exist_ok=True)
    summary = {"MOA": moa_results(args.cache, args.out), "Kepler": k2_results(args.cache, args.out)}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))
