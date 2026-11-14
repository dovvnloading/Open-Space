"""Localize the known eclipse clock using SkyMapper DR4 catalogue measurements."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

ROOT = Path(__file__).resolve().parent
TARGET = 363533777
SCALE = 3600.0
KEYS = ["east", "north", "mag_psf", "mag_apc02", "mag_apc05", "mag_apc08"]


def expected_positions(table):
    """Position on a tangent plane, relative to proper-motion/parallax-propagated Gaia."""
    time = Time(table.date + table.exp_time / 172800, format="mjd", scale="utc")
    earth = get_body_barycentric("earth", time).xyz.to_value(u.au)
    ra, dec = np.radians(table.ra), np.radians(table.dec)
    east_parallax = earth[0] * np.sin(ra) - earth[1] * np.cos(ra)
    north_parallax = (earth[0] * np.cos(ra) + earth[1] * np.sin(ra)) * np.sin(
        dec
    ) - earth[2] * np.cos(dec)
    dt = time.tcb.jyear - 2016.0
    table["east"] = (table.ra_img - table.ra) * SCALE * np.cos(dec) - (
        table.pmra * dt + table.parallax * east_parallax
    ) / 1000
    table["north"] = (table.decl_img - table.dec) * SCALE - (
        table.pmdec * dt + table.parallax * north_parallax
    ) / 1000


def gaia_waveform(images, row, target):
    time = Time(images.date + images.exp_time / 172800, format="mjd", scale="utc")
    ra, dec = np.radians([target.ra, target.dec])
    direction = [np.cos(ra) * np.cos(dec), np.sin(ra) * np.cos(dec), np.sin(dec)]
    light_days = (
        np.array(direction)
        @ get_body_barycentric("earth", time).xyz.to_value(u.au)
        * 499.00478384
        / 86400
    )
    phase = (
        (time.tcb.jd + light_days - 2455197.5 - row.reference_time) * row.frequency
    ) % 1
    amp = row.geom_model_cosine_half_period_amplitude
    model = amp * (
        1 + np.cos(4 * np.pi * (phase - row.geom_model_cosine_half_period_phase))
    )
    for i in [1, 2]:
        field = f"geom_model_gaussian{i}_"
        delta = (phase - row[field + "phase"] + 0.5) % 1 - 0.5
        model += row[field + "depth"] * np.exp(
            -0.5 * (delta / row[field + "sigma"]) ** 2
        )
    images["phase"] = phase
    images["inverse_flux_excess"] = 10 ** (0.4 * model) - 1


def register(epochs, source):
    good_ref = epochs[
        (epochs.object_id != source)
        & (epochs.object_id != TARGET)
        & (epochs["flags"] == 0)
        & (epochs.imaflags == 0)
        & (epochs.e_mag_psf < 0.08)
        & (epochs.ruwe < 1.4)
        & (epochs.ipd_frac_multi_peak <= 1)
    ]
    refs = good_ref.groupby(["image_id", "ccd"])
    zero = refs[["east", "north"]].median()
    zero.columns = ["east_zero", "north_zero"]
    zero["references"] = refs.size()
    result = epochs[epochs.object_id == source].join(zero, on=["image_id", "ccd"])
    result = result[
        (result.references >= 3)
        & (result["flags"] == 0)
        & (result.imaflags == 0)
        & (result.e_mag_psf < 0.03)
    ].copy()
    for k in ("east", "north"):
        result[k] -= result[k + "_zero"]
    return result


def coefficients(table, columns):
    # Filter-specific zero points permit fixed chromatic offsets.
    filters = pd.get_dummies(table["filter"], dtype=float).to_numpy()
    x = np.column_stack((filters, table.inverse_flux_excess))
    y = table[columns].to_numpy()
    return np.linalg.lstsq(x, y, rcond=None)[0][-1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dossier", type=Path, help="Open-Space results/tic-96100077 directory"
    )
    args = parser.parse_args()
    sources = pd.read_csv(ROOT / "sources.csv")
    objects = pd.read_csv(ROOT / "detections.csv.gz").merge(
        sources, on="object_id", validate="many_to_one"
    )
    images = pd.read_csv(ROOT / "images.csv")
    gaia = pd.read_csv(args.dossier / "data/gaia-model.csv").iloc[0]
    gaia_waveform(images, gaia, sources.loc[sources.object_id == TARGET].iloc[0])
    epochs = objects.merge(images, on="image_id", validate="many_to_one")
    expected_positions(epochs)
    target = register(epochs, TARGET)
    fitted = coefficients(target, KEYS)
    years = sorted(target.night_mjd.unique())
    jack = np.array(
        [coefficients(target[target.night_mjd != night], KEYS) for night in years]
    )
    centered = jack - jack.mean(axis=0)
    cov = (len(years) - 1) / len(years) * (centered.T @ centered)
    unc = np.sqrt(np.diag(cov))
    # Normalize the removed-light location to the measured PSF depth scale.
    # Magnitudes are locally 2.5 / ln(10) times inverse flux excess. Finite
    # depth is accounted for by fitting both brightness curves below.
    depth_reference = coefficients(
        target.assign(model_mag=2.5 * np.log10(1 + target.inverse_flux_excess)),
        ["model_mag"],
    )[0]
    measured_depth_scale = fitted[2] / depth_reference
    reference_mags = target.assign(
        model_mag=2.5 * np.log10(1 + target.inverse_flux_excess)
    )
    jack_depth_reference = np.array(
        [
            coefficients(
                reference_mags[reference_mags.night_mjd != night], ["model_mag"]
            )[0]
            for night in years
        ]
    )
    jack_offsets = -jack[:, :2] / (jack[:, 2] / jack_depth_reference)[:, None]
    offset_residuals = jack_offsets - jack_offsets.mean(axis=0)
    offset_cov = (len(years) - 1) / len(years) * (offset_residuals.T @ offset_residuals)
    summary = {
        "object_id": TARGET,
        "total_target_detections": int(sum(objects.object_id == TARGET)),
        "accepted_exposures": len(target),
        "nights": len(years),
        "filters": sorted(target["filter"].unique()),
        "reference_count_range": [
            int(target.references.min()),
            int(target.references.max()),
        ],
        "seeing_fwhm_mean_arcsec_range": target.fwhm_mean.agg(["min", "max"]).to_list(),
        "photometric_scale_relative_to_Gaia": dict(
            zip(KEYS[2:], fitted[2:] / depth_reference)
        ),
        "night_jackknife_se": dict(zip(KEYS, unc)),
        "fitted_slopes": dict(zip(KEYS, fitted)),
        "offset_east_north_arcsec": (-fitted[:2] / measured_depth_scale).tolist(),
        "offset_jackknife_east_north_arcsec": np.sqrt(np.diag(offset_cov)).tolist(),
        "leave_night_out_slope_min": dict(zip(KEYS, jack.min(axis=0))),
        "leave_night_out_slope_max": dict(zip(KEYS, jack.max(axis=0))),
    }
    cols = [
        "image_id",
        "ccd",
        "filter",
        "night_mjd",
        "phase",
        "inverse_flux_excess",
        "references",
        "fwhm_mean",
        *KEYS,
    ]
    target[cols].sort_values("image_id").to_csv(
        ROOT / "localized-epochs.csv", index=False, float_format="%.8f"
    )
    # Field controls are matched to the exact accepted image list, including the
    # target's temporal sampling. Their fitted flux and centroid slopes test
    # responses to shared frame errors at the target's ephemeris.
    rows = []
    for star in sources.object_id:
        if star == TARGET:
            continue
        phot = register(epochs, star)
        phot = phot[phot.image_id.isin(target.image_id)]
        if (
            len(phot) < 15
            or phot.night_mjd.nunique() < 5
            or np.ptp(phot.inverse_flux_excess) < 0.07
        ):
            continue
        vals = coefficients(phot, KEYS)
        rows.append([star, len(phot), phot.night_mjd.nunique(), *vals])
    pd.DataFrame(rows, columns=["object_id", "exposures", "nights", *KEYS]).to_csv(
        ROOT / "comparison-controls.csv", index=False, float_format="%.7f"
    )
    (ROOT / "localization.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
