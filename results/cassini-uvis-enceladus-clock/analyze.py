"""Replay temporal alignment and a held-out spectral check of the UVIS event."""

import argparse
from hashlib import file_digest
from itertools import product
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import xlogy

from check_clocks import HERE, calculate


FULL_FUV = "fuv2017_086_14_30.dat"
PHOTOMETER = "HSP2017_086_14_30.DAT"
WINDOW_FUV = "fuv-active-rows.u16be"
DT_FUV = 2.0
DT_HSP = 0.002
NULL = 65535
BASELINE = (510, 550)
ENCOUNTER = (550, 645)
QUIET = (60, 360)
SPECTRAL = (555, 605)
NEGATIVE_SPECTRAL = (300, 350)


def sha256(path):
    with path.open("rb") as stream:
        return file_digest(stream, "sha256").hexdigest()


def archive_entries():
    path = HERE / "archive" / "stage-receipt.json"
    return {entry["file"]: entry for entry in json.loads(path.read_text(encoding="utf-8"))["files"]}


def verify_file(path, entry):
    if path.stat().st_size != entry["bytes"] or sha256(path) != entry["sha256"]:
        raise ValueError(f"Input source disagrees with receipt: {path}")


def read_detector(root):
    receipt = archive_entries()
    if (root / FULL_FUV).is_file():
        source = root / FULL_FUV
        for filename in (FULL_FUV, "FUV2017_086_14_30_PDS3.DAT"):
            verify_file(root / filename, receipt[filename])
        if receipt[FULL_FUV]["sha256"] != receipt["FUV2017_086_14_30_PDS3.DAT"]["sha256"]:
            raise ValueError("Expected byte identity for the two FUV versions")
        cube = np.memmap(source, dtype=">u2", mode="r", shape=(465, 64, 1024))
        if np.any(cube[:, :30, :] != NULL) or np.any(cube[:, 35:, :] != NULL):
            raise ValueError("Unexpected samples outside the declared spatial window")
        active = np.asarray(cube[:, 30:35, :]).copy()
    else:
        manifest = json.loads((HERE / "detector-subset-receipt.json").read_text(encoding="utf-8"))
        source = root / WINDOW_FUV
        verify_file(source, manifest["fuv_subset"])
        active = np.fromfile(source, dtype=">u2").reshape(465, 5, 1024)
    hsp_source = root / PHOTOMETER
    verify_file(hsp_source, receipt[PHOTOMETER])
    hsp = np.fromfile(hsp_source, dtype=">u2")
    if active.shape != (465, 5, 1024) or hsp.shape != (466032,):
        raise ValueError("Unexpected detector array shape")
    return active, hsp


def hsp_for_fuv(hsp, fuv_end_first, hsp_end_first, frames):
    # HSP samples cover half-open intervals from one end-of-integration to the next.
    edges = hsp_end_first - DT_HSP + DT_HSP * np.arange(len(hsp) + 1)
    cumulative = np.hstack((0, np.cumsum(hsp, dtype=np.float64)))
    end = fuv_end_first + DT_FUV * frames
    start = end - DT_FUV
    return (np.interp(end, edges, cumulative, left=np.nan, right=np.nan)
            - np.interp(start, edges, cumulative, left=np.nan, right=np.nan))


def synthetic():
    fake = np.zeros(3000)
    fake[1500] = 12
    observed = hsp_for_fuv(fake, 4.001, DT_HSP, np.array([0]))[0]
    overlap = hsp_for_fuv(fake, 3.001, DT_HSP, np.array([0]))[0]
    absent = hsp_for_fuv(fake, 5.5, DT_HSP, np.array([0]))[0]
    if not np.allclose([observed, overlap, absent], [12, 6, 0]):
        raise RuntimeError("Integration partition did not recover a planted impulse")
    return {"contained_counts": round(float(observed), 8), "half_counts": round(float(overlap), 8),
            "outside_counts": round(float(absent), 8)}


def window(t, interval):
    return (t >= interval[0]) & (t < interval[1])


def diff_score(t, fuv, hsp, interval):
    use = window(t, interval) & np.isfinite(fuv) & np.isfinite(hsp)
    groups = np.floor((t - interval[0]) / 8).astype(int)
    g = groups[1:]
    pairs = use[:-1] & use[1:] & (groups[:-1] == g)
    train = pairs & (g % 2 == 0)
    test = pairs & (g % 2 == 1)
    x, y = np.diff(hsp), np.diff(fuv)
    denominator = np.dot(x[train], x[train])
    if not train.any() or not test.any() or denominator <= 0:
        return None
    slope = max(0, np.dot(x[train], y[train]) / denominator)
    residual = (y - slope * x) ** 2
    sums = {str(block): float(residual[test & (g == block)].sum()) for block in np.unique(g[test])}
    return {"training_pairs": int(train.sum()), "heldout_pairs": int(test.sum()),
            "slope": float(slope), "heldout_rms_fraction": float(np.sqrt(residual[test].mean())),
            "heldout_loss_by_block": sums}


def flips(a, b):
    blocks = set(a["heldout_loss_by_block"]) & set(b["heldout_loss_by_block"])
    delta = np.array([b["heldout_loss_by_block"][k] - a["heldout_loss_by_block"][k] for k in sorted(blocks)])
    if not len(delta):
        return None
    if len(delta) <= 14:
        signs = np.array(list(product((-1, 1), repeat=len(delta))))
        test_type = "exact sign enumeration"
    else:
        signs = np.random.default_rng(19860227).choice((-1, 1), size=(32768, len(delta)))
        test_type = "fixed-seed sign draws"
    sims = signs @ delta
    total = np.sum(delta)
    return {"method": test_type, "blocks": len(delta), "sclk_wins": int((delta > 1e-12).sum()),
            "mixed_wins": int((delta < -1e-12).sum()), "sclk_minus_mixed_loss": float(-total),
            "one_sided_rank_for_sclk": float(np.mean(sims >= total - 1e-12))}


def compare(t, fuv, counts, interval):
    reference = window(t, BASELINE)
    curve = fuv / np.nanmean(fuv[reference])
    norms = {key: values / np.nanmean(values[reference]) for key, values in counts.items()}
    scores = {key: diff_score(t, curve, values, interval) for key, values in norms.items()}
    use = window(t, interval)
    change = np.abs(counts["mixed_archive"][use] - counts["sclk_pair"][use])
    normalized_change = 100 * np.abs(norms["mixed_archive"][use] - norms["sclk_pair"][use])
    idx = np.flatnonzero(use)[np.nanargmax(change)]
    return {
        "frames": int(use.sum()),
        "median_abs_HSP_count_change": float(np.nanmedian(change)), "max_abs_HSP_count_change": float(np.nanmax(change)),
        "max_abs_HSP_normalized_change_percent": float(np.nanmax(normalized_change)),
        "max_count_change_frame_zero_based": int(idx),
        "max_count_change_sclk_midpoint_s_from_14_30": float(t[idx]),
        "max_count_frame_counts": {key: float(a[idx]) for key, a in counts.items()},
        "fast_change_scores": scores,
        "sclk_vs_mixed_blocks": flips(scores["sclk_pair"], scores["mixed_archive"])
        if scores["sclk_pair"] and scores["mixed_archive"] else None,
    }


def fit_poisson(train, design, offset, counts, bounds):
    def loss(params):
        eta = offset[train] + design[train] @ params
        a = np.exp(eta)
        val = np.sum(a - counts[train] * eta)
        return val, design[train].T @ (a - counts[train])

    fit = minimize(loss, np.zeros(design.shape[1]), jac=True, method="L-BFGS-B", bounds=bounds)
    return fit.x, bool(fit.success)


def check_spectrum(spectra, wavelength, sigma, t, interval):
    ref = window(t, BASELINE)
    obs = window(t, interval)
    use = np.flatnonzero((wavelength >= 1200) & (wavelength <= 1800))
    reference = np.sum(spectra[ref], axis=0)
    measured = np.sum(spectra[obs], axis=0)
    rows = []
    for block in range(len(use) // 16):
        at = use[16 * block:16 * (block + 1)]
        weight = reference[at]
        total = float(weight.sum())
        if total >= 40 and np.isfinite(sigma[at]).all():
            rows.append((block, float(np.average(wavelength[at], weights=weight)), total, float(measured[at].sum()),
                         float(np.average(sigma[at], weights=weight) * 1e16)))
    g, wav, n_ref, n_meas, optical = np.asarray(rows).T
    train = g % 2 == 0
    chrom = (wav - 1500) / 300
    offset = np.log(n_ref * obs.sum() / ref.sum())
    linear = np.column_stack((np.ones(len(g)), chrom))
    models = {}
    for name, design, bounds in (
        ("gray_color", linear, [(-5, 5), (-3, 3)]),
        ("gray_color_water", np.column_stack((linear, -optical)), [(-5, 5), (-3, 3), (0, 20)]),
    ):
        parameters, success = fit_poisson(train, design, offset, n_meas, bounds)
        prediction = np.exp(offset + design @ parameters)
        deviance = np.sum(2 * (xlogy(n_meas[~train], n_meas[~train] / prediction[~train])
                               - n_meas[~train] + prediction[~train]))
        models[name] = {"parameters": list(map(float, parameters)), "converged": success,
                        "heldout_deviance": float(deviance)}
    return {
        "event_frames": int(obs.sum()), "reference_frames": int(ref.sum()), "training_blocks": int(train.sum()),
        "heldout_blocks": int((~train).sum()), "models": models,
        "heldout_improvement_with_water": models["gray_color"]["heldout_deviance"]
        - models["gray_color_water"]["heldout_deviance"],
    }


def analysis(root):
    impulse = synthetic()
    active, hsp = read_detector(root)
    clock = calculate()
    ends = clock["first_end_seconds_since_utc_143000"]
    frames = np.arange(active.shape[0])
    t = ends["sclk_pair"][0] + (frames - 0.5) * DT_FUV
    band = np.loadtxt(HERE / "calibration" / "fuv-wavelength-bin-centers.csv", delimiter=",").ravel()
    water = np.loadtxt(HERE / "calibration" / "H2O-Mota-2005-298K.txt")
    sigma = np.interp(band / 10, water[:, 0], water[:, 1], left=np.nan, right=np.nan)
    if band.shape != (1024,):
        raise ValueError("Unexpected FUV wavelength table shape")
    broadband = (band >= 1200) & (band <= 1800)
    available = active != NULL
    counts_array = np.where(available, active, 0).astype(np.int64)
    summed_rows = counts_array[:, :, broadband].sum(axis=2)
    brightest = int(np.argmax(summed_rows[window(t, QUIET)].sum(axis=0)))
    selected_rows = slice(max(0, brightest - 1), min(5, brightest + 2))
    fuv = summed_rows.sum(axis=1)
    # Ignore intervals with less than 95% of the usually measured number of window cells.
    coverage = available[:, :, broadband].sum(axis=(1, 2))
    fuv_checked = np.where(coverage >= 0.95 * np.median(coverage), fuv, np.nan)
    core = np.where(np.isfinite(fuv_checked), summed_rows[:, selected_rows].sum(axis=1), np.nan)
    masks = {}
    for key, (fuv_end, hsp_end) in ends.items():
        masks[key] = hsp_for_fuv(hsp, fuv_end, hsp_end, frames)
    spec = counts_array.sum(axis=1)
    tested = {
        "array_support": {
            "fuv_frames": len(active), "hsp_samples": len(hsp),
            "median_valid_columns_by_spatial_row": list(map(int, np.median(available.sum(axis=2), axis=0))),
            "core_spatial_rows_zero_based": list(range(30 + selected_rows.start, 30 + selected_rows.stop)),
            "aperture_rows": [30, 31, 32, 33, 34],
            "hsp_peak_count": int(hsp.max()),
            "fuv_reference_row_counts": summed_rows[window(t, BASELINE)].sum(axis=0).tolist(),
            "fuv_encounter_row_counts": summed_rows[window(t, ENCOUNTER)].sum(axis=0).tolist(),
            "reference_frames": int(window(t, BASELINE).sum()),
        },
        "synthetic_integration": impulse,
        "encounter": compare(t, fuv_checked, masks, ENCOUNTER),
        "encounter_core": compare(t, core, masks, ENCOUNTER),
        "quiet": compare(t, fuv_checked, masks, QUIET),
        "spectral_encounter": check_spectrum(spec, band, sigma, t, SPECTRAL),
        "spectral_quiet": check_spectrum(spec, band, sigma, t, NEGATIVE_SPECTRAL),
    }
    curves = np.column_stack((t, fuv_checked, *(v for v in masks.values())))
    return tested, curves, active


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--curve", type=Path)
    parser.add_argument("--subset", type=Path)
    args = parser.parse_args()
    result, curves, active = analysis(args.raw)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    if args.curve:
        np.savetxt(args.curve, curves, delimiter=",", fmt="%.8f", comments="",
                   header="sclk_midpoint_s_from_14_30,fuv,both_pds3,mixed_archive,both_sclk")
    if args.subset:
        active.astype(">u2").tofile(args.subset)


if __name__ == "__main__":
    main()
