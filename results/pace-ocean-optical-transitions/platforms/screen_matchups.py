"""Screen the NASA PACE Rrs matchup release, including fixed-platform time series."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE = "https://pace.oceansciences.org/data/matchup_Rrs.csv"
DIGEST = "72f3d3187c9f39117fa3d926e652963839d830aae61f871c5c0774d6d6733b1d"
BANDS = range(400, 701, 5)
METRICS = {"blue": (443, 490, 560), "red": (560, 620, 665)}
DATES = ["2024-03-11", "2024-03-15", "2024-04-10", "2024-04-15"]
SITE = "USC_SEAPRISM"


def curvature(frame, instrument, a, b, c):
    bands = frame[[f"{instrument}_Rrs_{w}" for w in (a, b, c)]]
    logs = np.log(bands.where(bands > 0))
    return logs.iloc[:, 1] - ((c - b) * logs.iloc[:, 0] + (b - a) * logs.iloc[:, 2]) / (c - a)


def write_csv(frame, target):
    frame.to_csv(target, index=False, float_format="%.9g", na_rep="")


def add_metrics(frame):
    for label, triplet in METRICS.items():
        for instrument in ("SAT", "FIELD"):
            frame[f"{instrument}_{label}"] = curvature(frame, instrument, *triplet)
    # All eight interior bands must be present for this shape check.
    residuals = [curvature(frame, "SAT", 445, middle, 560) for middle in range(475, 511, 5)]
    frame["SAT_broad"] = pd.concat(residuals, axis=1).median(axis=1, skipna=False)
    return frame


def screen(site_days):
    records = []
    metric_cols = ["SAT_blue", "FIELD_blue", "SAT_red", "FIELD_red", "SAT_broad"]
    for site, series in site_days.groupby("site"):
        series = series.sort_values("date")
        for previous, current in zip(series.iloc[:-1].itertuples(), series.iloc[1:].itertuples()):
            gap = (current.date - previous.date).days
            if not 1 <= gap <= 8:
                continue
            reference_phase = current.date.replace(year=2000).dayofyear
            diff = (series["phase"] - reference_phase).abs()
            same_season = np.minimum(diff, 366 - diff) <= 30
            controls = series.loc[(series.date.dt.year != current.date.year) & same_season]
            row = {"site": site, "start": previous.date.date(), "end": current.date.date(), "days": gap}
            row["min_valid_pixels"] = min(previous.SAT_valid_pixels, current.SAT_valid_pixels)
            row["max_abs_qwip"] = max(abs(previous.SAT_QWIP_Score), abs(current.SAT_QWIP_Score))
            for key in metric_cols:
                values = controls[key].dropna()
                endpoint = getattr(current, key)
                row[f"{key}_delta"] = endpoint - getattr(previous, key)
                row[f"{key}_n_control"] = len(values)
                row[f"{key}_endpoint_z"] = np.nan
                if len(values) >= 6:
                    center = values.median()
                    scale = 1.4826 * (values - center).abs().median()
                    if scale > 0:
                        row[f"{key}_endpoint_z"] = (endpoint - center) / scale
            records.append(row)
    return pd.DataFrame(records)


def main(source, out, allow_new_release):
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != DIGEST and not allow_new_release:
        raise ValueError(f"Input release changed: {digest}. Use --allow-new-release to analyze it.")
    out.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(source, comment="#", na_values=-999, low_memory=False)
    data["date"] = pd.to_datetime(data.SAT_datetime, utc=True, format="ISO8601").dt.normalize()
    data["season"] = data.date.dt.month.map(lambda m: ("DJF", "MAM", "JJA", "SON")[(m % 12) // 3])
    data["site"] = data.FIELD_cruise.str.removeprefix("aoc_v4l15_")
    data = add_metrics(data)
    # The fixed platforms provide true revisits. Do not treat a cruise name as a location.
    fixed = data.loc[data.FIELD_experiment == "AERONET_OC_V4"].copy()
    fixed["abs_time_difference"] = pd.to_timedelta(fixed.SATFIELD_time_diff).abs()
    fixed = fixed.sort_values(["date", "abs_time_difference", "SAT_datetime"]).drop_duplicates(["site", "date"])
    fixed["phase"] = fixed.date.map(lambda date: date.replace(year=2000).dayofyear)
    intervals = screen(fixed)
    write_csv(intervals, out / "fixed_platform_intervals.csv")
    counts = fixed.groupby(["site", "season"]).agg(
        days=("date", "size"), latitude=("FIELD_latitude", "median"),
        longitude=("FIELD_longitude", "median"), first_day=("date", "min"), last_day=("date", "max"),
    ).reset_index()
    write_csv(counts, out / "site_season_coverage.csv")
    sat_cols = [f"SAT_Rrs_{w}" for w in BANDS]
    field_cols = [f"FIELD_Rrs_{w}" for w in BANDS]
    data["complete_oci"] = data[sat_cols].notna().all(axis=1)
    data["complete_field"] = data[field_cols].notna().all(axis=1)
    coverage = data.groupby(["FIELD_experiment", "season"]).agg(
        matchups=("date", "size"), complete_oci_61=("complete_oci", "sum"),
        complete_field_61=("complete_field", "sum"), min_lat=("FIELD_latitude", "min"),
        max_lat=("FIELD_latitude", "max"), min_lon=("FIELD_longitude", "min"),
        max_lon=("FIELD_longitude", "max"),
    ).reset_index()
    write_csv(coverage, out / "hyperspectral_coverage.csv")
    # The complete spectral pairs also cover moving ships and different water types.
    # A spectral angle measures direction after normalization, retaining signed retrievals.
    paired = data.loc[data.complete_oci & data.complete_field].copy()
    sat = paired[sat_cols].to_numpy(dtype=float)
    field = paired[field_cols].to_numpy(dtype=float)
    denom = np.linalg.norm(sat, axis=1) * np.linalg.norm(field, axis=1)
    cosine = np.divide((sat * field).sum(axis=1), denom, out=np.full(len(paired), np.nan), where=denom > 0)
    paired["spectral_angle_deg"] = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
    paired["median_abs_band_error_sr_inv"] = np.median(np.abs(sat - field), axis=1)
    pair_keys = ["FIELD_experiment", "FIELD_cruise", "FIELD_index", "FIELD_datetime", "FIELD_latitude",
                 "FIELD_longitude", "SAT_parent_file", "SAT_datetime", "SAT_cv", "SAT_valid_pixels",
                 "spectral_angle_deg", "median_abs_band_error_sr_inv"]
    write_csv(paired[pair_keys], out / "complete_field_spectral_pairs.csv")
    # The release exports 5-nm samples and the heritage wavelengths. Retain both.
    interesting = fixed.loc[(fixed.site == SITE) & fixed.date.dt.strftime("%Y-%m-%d").isin(DATES)].sort_values("date")
    metrics = ["SAT_blue", "FIELD_blue", "SAT_red", "FIELD_red", "SAT_broad"]
    keys = ["site", "SAT_parent_file", "SAT_datetime", "SAT_latitude", "SAT_longitude",
            "SAT_center_line", "SAT_center_pixel", "FIELD_datetime", "FIELD_latitude", "FIELD_longitude",
            "FIELD_index", "SAT_valid_pixels", "SAT_cv", "SAT_QWIP_Score", "SAT_senz", "SAT_solz",
            "SAT_elev", "SAT_chisqr_mbac", "SATFIELD_time_diff"]
    write_csv(interesting[keys + metrics], out / "california_matchups.csv")
    rows = []
    wavelengths = sorted(set(BANDS) | {412, 443, 532})
    for _, obs in interesting.iterrows():
        for wave in wavelengths:
            row = {"date": obs.date.date(), "wavelength_nm": wave}
            for field, key in (("SAT_Rrs", "oci_rrs"), ("FIELD_Rrs", "field_rrs"),
                               ("SAT_Rrs_unc", "oci_unc"), ("FIELD_Rrs_unc", "field_unc"),
                               ("SAT_Rrs_spatstd", "oci_spatial_sd")):
                row[key] = obs[f"{field}_{wave}"]
            rows.append(row)
    write_csv(pd.DataFrame(rows), out / "california_spectra.csv")
    print(f"Source {digest}; {len(data)} matchups; {len(fixed)} platform-days; {len(intervals)} intervals.")
    for metric in METRICS:
        same_sign = intervals[f"SAT_{metric}_delta"] * intervals[f"FIELD_{metric}_delta"] > 0
        coherent = intervals.loc[same_sign].copy()
        coherent["minimum_abs_delta"] = coherent[[f"SAT_{metric}_delta", f"FIELD_{metric}_delta"]].abs().min(axis=1)
        print(metric, coherent.nlargest(3, "minimum_abs_delta")[["site", "start", "end", "minimum_abs_delta"]].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help=f"Local copy of {SOURCE}")
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--allow-new-release", action="store_true")
    args = parser.parse_args()
    main(args.source, args.out, args.allow_new_release)
