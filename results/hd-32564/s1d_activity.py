"""Measure H-alpha and calcium line-core indices in the public ESPRESSO spectra."""
import argparse
import csv
import io
import json
import tempfile
from pathlib import Path

import numpy as np

from replay import EPOCH, ROOT, design, remove
from s1d_velocities import C, spectrum
from verify_sources import read_csv


def mean_flux(wave, flux, error, quality, center, width, triangle=False):
    # Rectangles are specified by full width and triangles by their FWHM.
    half_extent = width if triangle else width / 2
    offset = (wave - center) / half_extent
    in_band = np.abs(offset) < 1
    weights = np.gradient(wave[in_band])
    if triangle:
        weights *= 1 - np.abs(offset[in_band])
    keep = (quality[in_band] == 0) & np.isfinite(flux[in_band]) & np.isfinite(error[in_band]) & (error[in_band] > 0)
    if weights[keep].sum() < 0.95 * weights.sum():
        raise ValueError(f"Insufficient usable coverage at {center} Angstrom")
    return np.average(flux[in_band][keep], weights=weights[keep])


def indices(row, cache):
    wave, flux, err, quality, digest = spectrum(row, cache)
    # Fixed system velocity in m/s, common to all exposures.
    beta = -40_695.0 / C
    rest = wave / np.sqrt((1 + beta) / (1 - beta))
    air = rest / (1 + 2.735182e-4 + 131.4182 / rest**2 + 2.76249e8 / rest**4)

    def band(center, width, triangle=False):
        return mean_flux(air, flux, err, quality, center, width, triangle)

    ha_reference = band(6550.87, 10.75) + band(6580.31, 8.75)
    ca_reference = band(3901.07, 20) + band(4001.07, 20)
    ca_cores = band(3933.664, 1.09, True) + band(3968.47, 1.09, True)
    values = (band(6562.808, 1.6) / ha_reference, band(6562.808, 0.6) / ha_reference,
              ca_cores / ca_reference, band(6572.795, 0.8) / ha_reference)
    profiles = []
    for name, center, reference in [("CaII_K", 3933.664, ca_reference),
                                    ("CaII_H", 3968.47, ca_reference),
                                    ("Halpha", 6562.808, ha_reference)]:
        # Non-overlapping 0.05 Angstrom bins, normalized to the mean continuum.
        for delta in np.arange(-40, 41) * 0.05:
            profiles.append((name, center + delta, 2 * band(center + delta, 0.05) / reference))
    return digest, values, profiles


def extract(rows, cache):
    buffer = io.StringIO()
    out = csv.writer(buffer, lineterminator="\n")
    out.writerow(["product", "bjd", "prefix_sha256", "halpha16", "halpha06", "cahk", "cai_control"])
    profiles = []
    for row in rows:
        digest, measurements, profile = indices(row, cache)
        profiles.append(profile)
        out.writerow([row["product"], row["bjd"], digest, *[f"{x:.9f}" for x in measurements]])
    profile_csv = io.StringIO()
    writer = csv.writer(profile_csv, lineterminator="\n")
    writer.writerow(["feature", "rest_air_angstrom", "p16", "median", "p84",
                     *[row["product"] for row in rows]])
    for epoch_bins in zip(*profiles):
        label, center, _ = epoch_bins[0]
        fluxes = [flux for _, _, flux in epoch_bins]
        quantiles = np.percentile(fluxes, [16, 50, 84])
        writer.writerow([label, f"{center:.3f}", *[f"{v:.7f}" for v in (*quantiles, *fluxes)]])
    return buffer.getvalue(), profile_csv.getvalue()


def correlations(rows):
    t = np.array([float(r["bjd"]) for r in rows]) - EPOCH
    base = np.column_stack([np.ones(len(t)), t - t.mean()])
    models = json.loads((ROOT / "harps_predictions.json").read_text())
    best = min(models, key=lambda row: row["rms"])
    terms = (design(t, best["P"]) * best["params"][:8]).reshape(-1, 4, 2).sum(axis=2)
    signals = {"short": terms[:, 0], "second": terms[:, 1], "full": terms.sum(axis=1)}
    for key in ("halpha16", "halpha06", "cahk", "cai_control"):
        values = np.array([float(r[key]) for r in rows])
        print(key, "median", round(float(np.median(values)), 6), "r",
              {name: round(float(np.corrcoef(remove(values, base), remove(wave, base))[0, 1]), 4)
               for name, wave in signals.items()})
    activity = np.array([float(r["halpha16"]) for r in rows])
    measured_rv = np.array([float(r["rv"]) * 1000 for r in read_csv(ROOT / "espresso.csv")])
    season = np.floor((t + EPOCH - 2460175) / 365.25)
    season_base = np.column_stack([season == value for value in np.unique(season)])
    seasonal_r = np.corrcoef(remove(activity, season_base), remove(signals["full"], season_base))[0, 1]
    print("H-alpha versus full waveform after season offsets, r", round(float(seasonal_r), 4))
    cai = np.array([float(r["cai_control"]) for r in rows])
    print("H-alpha versus Ca I, r", round(float(np.corrcoef(remove(activity, base), remove(cai, base))[0, 1]), 4))
    control = np.column_stack([base, cai])
    r_partial = np.corrcoef(remove(activity, control), remove(signals["full"], control))[0, 1]
    print("H-alpha versus full waveform after Ca I correction, r", round(float(r_partial), 4))
    base = np.column_stack([base, activity])
    for key, wave in [("null", 0), ("short", signals["short"]), ("full", signals["full"])]:
        print(key, "CCF RMS with H-alpha, m/s", round(float(np.std(remove(measured_rv - wave, base))), 4))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", action="store_true")
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--output-dir", type=Path, help="write remeasured tables in this directory")
    args = parser.parse_args()
    if args.output_dir and not args.extract:
        parser.error("--output-dir requires --extract")
    rows = list(csv.DictReader((ROOT / "espresso.csv").read_text().splitlines()))
    released = ROOT / "s1d-activity.csv"
    if args.extract:
        with tempfile.TemporaryDirectory() as tmp:
            cache = args.cache or Path(tmp)
            cache.mkdir(exist_ok=True, parents=True)
            data, profile_data = extract(rows, cache)
        if args.output_dir:
            args.output_dir.mkdir(exist_ok=True, parents=True)
            (args.output_dir / released.name).write_text(data)
            (args.output_dir / "s1d-profiles.csv").write_text(profile_data)
        elif data != released.read_text() or profile_data != (ROOT / "s1d-profiles.csv").read_text():
            raise ValueError("Measured indices or profiles differ from release")
    else:
        data = released.read_text()
    measured = list(csv.DictReader(data.splitlines()))
    manifest = list(csv.DictReader((ROOT / "s1d-velocities.csv").read_text().splitlines()))
    for key in ("product", "bjd", "prefix_sha256"):
        if [r[key] for r in measured] != [r[key] for r in manifest]:
            raise ValueError(f"Activity table differs from velocity manifest for {key}")
    correlations(measured)


if __name__ == "__main__":
    main()
