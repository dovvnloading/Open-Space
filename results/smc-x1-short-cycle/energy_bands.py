"""Measure three MAXI/GSC bands in four cycles defined by separate Swift data."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path

import numpy as np

from replay import ROOT

BANDS = ("2-4", "4-10", "10-20")
PINNED = ROOT / "data/maxi-orbit-bands.csv.gz"
HEADER = "mjd,full,full_err,soft,soft_err,medium,medium_err,hard,hard_err"
SOURCE = json.loads((ROOT / "data/sources.json").read_text())["products"]["smc_maxi_orb"]
MINIMA = np.array(json.loads((ROOT / "measurements.json").read_text())["minima"]["Swift"]["mjd"])


def pack(raw_file, output):
    """Keep complete original rows across the four cycles, before any quality cuts."""
    raw = raw_file.read_bytes()
    if hashlib.sha256(raw).hexdigest() != SOURCE["sha256"]:
        raise ValueError("Mission product differs from the original manifest")
    rows = [line.split() for line in raw.splitlines() if line.strip()]
    selected = [b",".join(row) for row in rows if MINIMA[0] <= float(row[0]) < MINIMA[-1]]
    payload = HEADER.encode() + b"\n" + b"\n".join(selected) + b"\n"
    with output.open("wb") as stream, gzip.GzipFile(filename="", fileobj=stream, mode="wb", mtime=0) as zipped:
        zipped.write(payload)


def daily_flux(rows, selected):
    # The same weights for all three bands preserve their relative contributions.
    # Equal weight per ISS orbit bin within a day, then equal weight per sampled day.
    t = rows[selected, 0]
    flux = rows[selected, 3::2]
    return np.stack([flux[t.astype(int) == day].mean(axis=0) for day in np.unique(t.astype(int))])


def states(rows, phase_mask, half_width):
    dates = rows[:, 0]
    report = []
    daily = []
    for left, right in zip(MINIMA[:-1], MINIMA[1:]):
        phase = (dates - left) / (right - left)
        inside = phase_mask & (phase >= 0) & (phase < 1)
        middle = inside & (np.abs(phase - 0.5) <= half_width)
        edge = inside & (np.minimum(phase, 1 - phase) <= half_width)
        mid, low = (daily_flux(rows, sel) for sel in (middle, edge))
        daily.append((mid, low))
        report.append({
            "start_mjd": left, "stop_mjd": right,
            "middle_edge_days": [len(mid), len(low)],
            "middle_edge_orbits": [int(middle.sum()), int(edge.sum())],
            "middle_minus_edge": rounded(mid.mean(axis=0) - low.mean(axis=0))
        })
    return daily, report


def rounded(values):
    return np.asarray(values).round(6).tolist()


def metrics(mid, low, hardness):
    diff = mid - low
    # Fractions of the sum of the three subbands, with no use of the overlapping full band.
    if not hardness:
        return np.concatenate([mid, low, diff], axis=-1)
    mid_fraction = mid[..., 2] / mid.sum(axis=-1)
    low_fraction = low[..., 2] / low.sum(axis=-1)
    return np.concatenate([mid, low, diff, (mid_fraction - low_fraction)[..., None]], axis=-1)


def day_resamples(daily, rng, draws=10000):
    # Keep band covariance within each day and resample within phase and cycle.
    samples = [[sample[rng.integers(len(sample), size=(draws, len(sample)))].mean(axis=1)
                for sample in pair] for pair in daily]
    return np.array(samples).mean(axis=0)


def summarize(daily, hardness):
    means = np.array([[sample.mean(axis=0) for sample in pair] for pair in daily])
    values = metrics(*means.mean(axis=0), hardness)
    samples = day_resamples(daily, np.random.default_rng(32564))
    bootstrap = metrics(*samples, hardness)
    lows, highs = np.quantile(bootstrap, [0.025, 0.975], axis=0)
    drop_one = np.array([metrics(*np.delete(means, i, axis=0).mean(axis=0), hardness)
                         for i in range(len(daily))])
    result = {
        "middle": rounded(values[:3]), "edge": rounded(values[3:6]),
        "middle_minus_edge": rounded(values[6:9]),
        "day_bootstrap_95_contrast": rounded(np.c_[lows[6:9], highs[6:9]]),
        "leave_one_cycle_out_contrast_range": rounded(np.c_[drop_one[:, 6:9].min(0), drop_one[:, 6:9].max(0)]),
    }
    if hardness:
        result.update({
            "middle_minus_edge_hard_fraction": round(float(values[9]), 6),
            "day_bootstrap_95_hard_fraction_difference": rounded([lows[9], highs[9]]),
            "leave_one_cycle_out_hard_fraction_difference": rounded(drop_one[:, 9])
        })
    return result


def measure(path):
    payload = gzip.decompress(path.read_bytes())
    rows = np.loadtxt(io.BytesIO(payload), delimiter=",", skiprows=1)
    errors = rows[:, 2::2]
    keep = np.isfinite(rows).all(axis=1) & (errors > 0).all(axis=1) & (rows[:, 2] < SOURCE["error_ceiling"])
    dates = rows[:, 0]
    if not (np.all(np.diff(dates) > 0) and dates.min() >= MINIMA[0] and dates.max() < MINIMA[-1]):
        raise ValueError("Dates must be unique, ordered and inside the four cycles")
    dt = dates - 60505.7361918
    orbital_phase = np.abs((.25696017 * dt + .5 * 2.379e-9 * dt**2 + .5) % 1 - .5)
    result = {
        "source_sha256": SOURCE["sha256"], "csv_sha256": hashlib.sha256(payload).hexdigest(),
        "bands_keV": BANDS, "orbit_bins": len(rows), "eligible_orbit_bins": int(keep.sum()),
        "units": "photons cm^-2 s^-1; hard fraction is dimensionless"
    }
    selections = {}
    for name, mask, half_width in [
        ("out_of_eclipse", orbital_phase >= .17, 1/6),
        ("out_of_eclipse_narrower_windows", orbital_phase >= .17, 1/8),
        ("eclipse_core", orbital_phase <= .055, 1/6),
    ]:
        daily, cycles = states(rows, keep & mask, half_width)
        selections[name] = daily
        result[name] = {"half_width_cycles": half_width, "cycles": cycles, **summarize(daily, hardness=(name != "eclipse_core"))}
    # Subtract the same-phase eclipse core; it can contain real emission from the binary.
    outside, core = selections["out_of_eclipse"], selections["eclipse_core"]
    def state_means(daily):
        return np.mean([[sample.mean(0) for sample in pair] for pair in daily], axis=0)

    mid, edge = state_means(outside) - state_means(core)
    rng = np.random.default_rng(60700)
    boot_mid, boot_edge = day_resamples(outside, rng) - day_resamples(core, rng)
    physical = (boot_mid.sum(axis=1) > 0) & (boot_edge.sum(axis=1) > 0)
    corrected = metrics(mid, edge, True)
    boot = metrics(boot_mid[physical], boot_edge[physical], True)
    result["out_of_eclipse_minus_eclipse_core"] = {
        "middle": rounded(mid), "edge": rounded(edge), "middle_minus_edge": rounded(corrected[6:9]),
        "middle_minus_edge_hard_fraction": round(float(corrected[-1]), 6),
        "day_bootstrap_95_hard_fraction_difference": rounded(np.quantile(boot[:, -1], [0.025, 0.975])),
        "bootstrap_draws_with_positive_summed_flux_in_both_states": int(physical.sum()),
        "bootstrap_draws_total": len(physical)
    }
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-orbit", type=Path, help="Rebuild the pinned CSV using these mission bytes")
    parser.add_argument("--data", type=Path, default=PINNED, help="Read or rebuild this gzip CSV")
    parser.add_argument("--output", type=Path, help="Write the JSON result to this path")
    args = parser.parse_args()
    if args.raw_orbit:
        pack(args.raw_orbit, args.data)
    output = json.dumps(measure(args.data), indent=2) + "\n"
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end="")
