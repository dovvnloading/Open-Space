"""Check a spatially resolved Gaia DR3 light curve against a TESS counterpart."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.timeseries import LombScargle


def analyze(path):
    raw = path.read_bytes()
    rows = list(csv.DictReader(raw.decode().splitlines()))
    if len({r["source_id"] for r in rows}) != 1:
        raise ValueError("An epoch photometry file must contain exactly one Gaia source")
    good = [r for r in rows if r["g_transit_flux"] and r["g_transit_flux_error"]
            and r["variability_flag_g_reject"] == "false" and r["rejected_by_photometry"] == "false"]
    t, f, err = (np.array([float(row[k]) for row in good])
                 for k in ("g_transit_time","g_transit_flux","g_transit_flux_error"))
    ok = np.isfinite(t) & np.isfinite(f) & np.isfinite(err) & (err > 0)
    t, f, err = t[ok], f[ok], err[ok]
    if len(t) < 10:
        raise ValueError("Too few accepted Gaia G epochs")
    ls = LombScargle(t, f, err)
    freq, power = ls.autopower(minimum_frequency=1/8, maximum_frequency=1/.15, samples_per_peak=16)
    k = int(np.argmax(power))
    report = {"gaia_dr3_source_id": rows[0]["source_id"], "epoch_csv_sha256": hashlib.sha256(raw).hexdigest(),
              "rows_total": len(rows), "accepted_g_transits": len(t), "baseline_days": float(np.ptp(t)),
              "frequency_bounds_per_day": [1/8, 1/.15], "frequencies_tested": len(freq), "samples_per_peak": 16,
              "peak_period_days": float(1/freq[k]), "peak_normalized_power": float(power[k]),
              "observed_flux_unit": "electrons/s", "phase_zero_gaia_bjd_minus_2455197p5": float(np.min(t))}
    return report, (t, f, err)


def plot(report, curve, output, name):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    t, f, err = curve
    p = report["peak_period_days"]
    phase = ((t-np.min(t))/p) % 1
    plt.rcParams.update({"svg.fonttype": "none", "svg.hashsalt": "resolved-gaia-period"})
    fig, ax = plt.subplots(figsize=(7.6, 3.6), layout="constrained")
    ax.errorbar(phase, f/1000, err/1000, fmt="o", ms=5, lw=1, color="#00766f")
    ax.set(xlabel="Phase (cycles)", ylabel="Gaia G flux (10³ electrons/s)", xlim=(-.03,1.03))
    ax.set_title(f"{name}  ·  {len(t)} Gaia G measurements at {p:.6f} days", loc="left", weight="bold")
    ax.grid(axis="y", color="#e0e3e2")
    ax.spines[["top","right"]].set_visible(False)
    fig.savefig(output, dpi=160, metadata={"Date":None} if output.suffix.lower()==".svg" else None)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("epochs", type=Path)
    p.add_argument("--figure", type=Path)
    p.add_argument("--name", help="Counterpart name in figure title")
    args = p.parse_args()
    report, curve = analyze(args.epochs)
    if args.figure:
        plot(report, curve, args.figure, args.name or "Gaia DR3 " + report["gaia_dr3_source_id"])
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
