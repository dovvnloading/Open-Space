"""Census of sustained factor-ten switches in the Fermi-LAT LCR.

The quantity measured is a median flux level in each of four nonoverlapping
half-years. Per-bin intervals express the published uncertainty convention.
Their median envelopes are descriptive bounds, not simultaneous confidence
intervals. The two native cadence products in each check reuse photons.
"""
import numpy as np


HALF_YEAR = 365.25 / 2
MONTH = 365.25 / 12


def _prepare(curve):
    if curve is None or not len(curve.get("t", [])):
        return np.empty((0, 3))
    t = np.asarray(curve["t"], dtype=float)
    flux = np.asarray(curve["flux"], dtype=float)
    err = np.asarray(curve["err"], dtype=float)
    ul = np.asarray(curve["ul"], dtype=float)
    fit = np.asarray(curve.get("fit", np.full(t.size, np.nan)), dtype=float)
    good = np.isfinite(t) & ((~np.isfinite(fit)) | (fit == 0))
    detection = good & np.isfinite(flux) & (flux > 0) & np.isfinite(err) & (err > 0)
    limit = good & ~detection & np.isfinite(ul) & (ul > 0)
    keep = detection | limit
    lower = np.where(detection, np.maximum(0, flux - 1.96 * err), 0)
    upper = np.where(detection, flux + 1.96 * err, ul)
    rows = np.column_stack([t[keep], lower[keep], upper[keep]])
    rows = rows[np.argsort(rows[:, 0])]
    if len(rows) > 1 and np.any(np.diff(rows[:, 0]) == 0):
        raise ValueError("duplicate time bins in normalized state input")
    return rows


def _blocks(rows, split, count):
    """Bounds valid for any values assigned to missing or censored bins."""
    edges = split + np.arange(-2, 3) * HALF_YEAR
    indices = np.searchsorted(rows[:, 0], edges)
    blocks = []
    for i, j in zip(indices[:-1], indices[1:]):
        valid = int(j - i)
        if valid < int(np.ceil(0.80 * count)):
            return None
        pad = max(0, count - valid)
        lower = np.concatenate((rows[i:j, 1], np.zeros(pad)))
        upper = np.concatenate((rows[i:j, 2], np.full(pad, np.inf)))
        blocks.append([float(np.median(lower)), float(np.median(upper)), valid])
    return blocks


def _ratios(blocks):
    low, high, _ = np.asarray(blocks).T
    return {
        "on": min(low[2:]) / max(high[:2]),
        "off": min(low[:2]) / max(high[2:]),
    }


def scan_state(monthly, weekly, three_day=None):
    """Prefer monthly/weekly; use weekly/three-day only if monthly is unavailable."""
    has_monthly = monthly is not None and len(monthly.get("t", [])) > 0
    week = _prepare(weekly)
    month = _prepare(monthly)
    three = _prepare(three_day)
    if has_monthly:
        first_name, second_name = "monthly", "weekly"
        first_rows, second_rows = month, week
        first_count, second_count = 6, 26
        first_width, second_width = MONTH, 7.
    else:
        first_name, second_name = "weekly", "3day"
        first_rows, second_rows = week, three
        first_count, second_count = 26, 61
        first_width, second_width = 7., 3.
    out = {
        "status": "insufficient_coverage", "transitions": [], "n_splits": 0,
        "primary_cadence": first_name, "confirmation_cadence": second_name,
        "n_monthly_valid": len(month), "n_weekly_valid": len(week), "n_3day_valid": len(three),
        "n_joint_splits": 0, "best_ratio": None, "best_split_mjd": None,
    }
    if (len(first_rows) < 4 * int(np.ceil(0.80 * first_count))
            or len(second_rows) < 4 * int(np.ceil(0.80 * second_count))):
        return out
    first = max(first_rows[0, 0] - first_width / 2,
                second_rows[0, 0] - second_width / 2) + 2 * HALF_YEAR
    last = min(first_rows[-1, 0] + first_width / 2,
               second_rows[-1, 0] + second_width / 2) - 2 * HALF_YEAR
    # Fixed monthly trial grid in MJD, independent of individual fluxes.
    start = 54682.0 + np.ceil((first - 54682.0) / MONTH) * MONTH
    qualified = []
    for split in np.arange(start, last + 1e-6, MONTH):
        out["n_splits"] += 1
        first_bounds = _blocks(first_rows, split, first_count)
        second_bounds = _blocks(second_rows, split, second_count)
        if first_bounds is None or second_bounds is None:
            continue
        out["n_joint_splits"] += 1
        first_ratio, second_ratio = _ratios(first_bounds), _ratios(second_bounds)
        for direction in ("on", "off"):
            ratio = float(min(first_ratio[direction], second_ratio[direction]))
            if out["best_ratio"] is None or ratio > out["best_ratio"]:
                out.update(best_ratio=ratio, best_split_mjd=float(split))
            if ratio >= 10:
                row = dict(direction=direction, split_mjd=float(split), ratio_floor=ratio,
                           primary_cadence=first_name, confirmation_cadence=second_name,
                           monthly_ratio=None, weekly_ratio=None, three_day_ratio=None,
                           monthly_blocks=None, weekly_blocks=None, three_day_blocks=None)
                for cadence, bounds, by_direction in (
                        (first_name, first_bounds, first_ratio), (second_name, second_bounds, second_ratio)):
                    field = "three_day" if cadence == "3day" else cadence
                    row[field + "_ratio"] = float(by_direction[direction])
                    # Each row contains median low, median high, usable bins.
                    row[field + "_blocks"] = bounds
                qualified.append(row)
    if not out["n_joint_splits"]:
        return out
    out["status"] = "no_qualified_transition"
    # Adjacent passing trial dates identify one event; retain positional range
    # and the strongest joint flux bound in each contiguous group.
    for row in qualified:
        transitions = out["transitions"]
        if (transitions and row["direction"] == transitions[-1]["direction"]
                and row["split_mjd"] - transitions[-1]["split_last_mjd"] < 1.01 * MONTH):
            prev = transitions[-1]
            if row["ratio_floor"] > prev["ratio_floor"]:
                prev.update(row)
            prev["split_last_mjd"] = row["split_mjd"]
            prev["n_passing_splits"] += 1
        else:
            transitions.append(dict(
                row, split_first_mjd=row["split_mjd"], split_last_mjd=row["split_mjd"],
                n_passing_splits=1,
            ))
    if out["transitions"]:
        out["status"] = "qualified"
    return out
