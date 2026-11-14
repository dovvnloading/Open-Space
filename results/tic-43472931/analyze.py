"""Test the published TESS periods on the two separately resolved ZTF sources."""
from pathlib import Path
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
# Kostov et al. (2024) epochs and durations, with the independently published
# 2025 TESS Ten Thousand period for A. Epochs are BJD minus 2,457,000.
EPHEMERIDES = {"A": (1.361339, 1440.2465, 1.8), "B": (1.4851, 1443.9580, 2.3)}
BRIGHT = (74.6819714436, -10.7558382952)
FAINT = (74.6808637906, -10.7570526063)


def distance(ra, dec, position):
    """Small-angle separation in arcseconds, valid over this nine-arcsecond cone."""
    ra0, dec0 = position
    return 3600 * np.hypot((ra - ra0) * np.cos(np.deg2rad(dec0)), dec - dec0)


def phase_distance(time, period, epoch):
    return abs(((time - epoch) / period + 0.5) % 1.0 - 0.5)


def masks(time, clock, offset=0):
    period, epoch, duration = EPHEMERIDES[clock]
    # Select the central 80% of the reported eclipse duration. Keep the
    # between-eclipse reference nearby in phase, and clear of the primary.
    phase = phase_distance(time, period, epoch + period * offset)
    halfwidth = 0.4 * duration / (24 * period)
    inside = phase <= halfwidth
    outside = (phase > 0.12) & (phase < 0.32)
    other = "A" if clock == "B" else "B"
    op, oe, _ = EPHEMERIDES[other]
    keep = phase_distance(time, op, oe) > 0.085
    if other == "A":
        keep &= phase_distance(time, op, oe + op * 0.5038) > 0.085
    return inside & keep, outside & keep


def depth(group, clock, rng):
    # Each UTC night contributes once per mask, source, and band. Bootstrap the
    # union of whole nights, preserving simultaneous baseline and eclipse points.
    time = group.hjd.to_numpy() - 2457000
    inside, outside = masks(time, clock)
    table = group.assign(night=np.floor(group.mjd).astype(int))
    def per_night(mask):
        return table[mask].groupby("night").mag.median()
    events, baseline = per_night(inside), per_night(outside)
    both = pd.concat([events.rename("in"), baseline.rename("out")], axis=1)
    delta = events.median() - baseline.median()
    lo, hi = np.nan, np.nan
    if min(len(events), len(baseline)) >= 4:
        take = rng.integers(0, len(both), (2500, len(both)))
        in_draws = both["in"].to_numpy()[take]
        out_draws = both["out"].to_numpy()[take]
        ok = np.isfinite(in_draws).any(axis=1) & np.isfinite(out_draws).any(axis=1)
        sims = np.nanmedian(in_draws[ok], axis=1) - np.nanmedian(out_draws[ok], axis=1)
        lo, hi = np.percentile(sims, [2.5, 97.5])
    return {"in_nights": len(events), "out_nights": len(baseline),
            "delta_mag": round(delta, 6), "bootstrap_low": round(lo, 6), "bootstrap_high": round(hi, 6)}


def main():
    raw = pd.read_csv(DATA / "ztf.csv.gz")
    clean = raw[raw.catflags == 0].copy()
    clean["source"] = np.where(distance(clean.ra, clean.dec, BRIGHT) <
                               distance(clean.ra, clean.dec, FAINT), "bright", "faint")
    rng = np.random.default_rng(43472931)
    rows = []
    diagnostics = []
    for (source, band), group in clean.groupby(["source", "filtercode"]):
        center = BRIGHT if source == "bright" else FAINT
        diagnostics.append(dict(source=source, band=band, epochs=len(group),
                                median_mag=round(group.mag.median(), 5),
                                median_sep_arcsec=round(float(np.median(distance(group.ra, group.dec, center))), 5)))
        for clock in EPHEMERIDES:
            for subset, data in [("full", group), ("early", group[group.hjd < 2459600]),
                                 ("late", group[group.hjd >= 2459600])]:
                if data.empty:
                    continue
                rows.append(dict(source=source, band=band, clock=clock, subset=subset,
                                 **depth(data, clock, rng)))
    pd.DataFrame(rows).to_csv(HERE / "depths.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(HERE / "photometry.csv", index=False)
    # Test the maximum possible contribution of the neighboring star in each
    # ZTF filter, allowing complete extinction of its median measured flux.
    diag = pd.DataFrame(diagnostics).pivot(index="band", columns="source", values="median_mag")
    contrast = diag["faint"] - diag["bright"]
    out = {"raw_epochs": len(raw), "clean_epochs": len(clean),
           "hjd_range": [raw.hjd.min(), raw.hjd.max()],
           "neighbor_max_combined_depth": dict((1 / (1 + 10**(0.4 * contrast))).round(6))}
    (HERE / "summary.json").write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
