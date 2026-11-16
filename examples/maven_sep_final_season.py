#!/usr/bin/env python3
"""Independent 2025 Aug--Dec MAVEN SEP dual-species, four-look event audit.

Reanalysis of LASP MAVEN KP via a public, dated Parquet mirror. The mirror's
SEP '*_quality' values are Poisson 1-sigma flux uncertainties (NOT 0/1
quality flags); see the PDS KP SIS and SPASE parameter definition.

This is a *screen*, not an assertion of a previously unknown solar eruption.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

MIRROR = "https://huggingface.co/datasets/juliensimon/nasa-maven-kp-insitu/resolve/main/data/{0}/{0}.parquet"
NAMES = ("1f", "1r", "2f", "2r")
ION = ["sep_ion_fov" + name for name in NAMES]
ELECTRON = ["sep_electron_fov" + name for name in NAMES]
FLUX_UNITS = "(cm^2 s sr)^-1"
SOURCE_SIS = "https://pds-ppi.igpp.ucla.edu/data/maven-insitu-calibrated/document/maven_sis_kp_insitu_r05.2.pdf"
SOURCE_KP = "https://lasp.colorado.edu/maven/sdc/public/data/sci/kp/insitu/"
MAX_VALUE = 1e10


def download_parquet(instrument: str, destination: Path) -> dict:
    url = MIRROR.format(instrument)
    sha = hashlib.sha256()
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "Open-Space-research/1.0 (public reproducibility)"})
            with urllib.request.urlopen(request, timeout=180) as response, destination.open("wb") as handle:
                for block in iter(lambda: response.read(1 << 20), b""):
                    sha.update(block)
                    handle.write(block)
            return {"url": url, "sha256": sha.hexdigest(), "bytes": destination.stat().st_size}
        except (OSError, urllib.error.URLError) as exc:
            if attempt == 3:
                raise RuntimeError("Download unavailable for " + url) from exc
            sha = hashlib.sha256()
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def read_table(path: Path, names: list[str]) -> pd.DataFrame:
    frame = pd.read_parquet(path, columns=names)
    frame["time"] = pd.to_datetime(frame["time"], utc=True, errors="coerce")
    frame = frame.dropna(subset=["time"]).sort_values("time")
    return frame.drop_duplicates(subset=["time"]).reset_index(drop=True)


def clean_flux(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    return x.where(x.ge(0) & x.lt(MAX_VALUE))


def sep_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Predeclared all-look and 3-of-4 SNR checks, rejecting invalid values."""
    data = df[["time"]].copy()
    for kind, columns in (("ion", ION), ("electron", ELECTRON)):
        values = np.column_stack([clean_flux(df[c]).to_numpy(dtype=float) for c in columns])
        errors = np.column_stack([clean_flux(df[c + "_quality"]).to_numpy(dtype=float) for c in columns])
        errors[errors <= 0] = np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            snr = values / errors
        # A median is meaningful only if >=3 independent FOVs have valid fluxes.
        valid = np.isfinite(values).sum(axis=1) >= 3
        med = np.nanmedian(values, axis=1)
        med[~valid] = np.nan
        data[kind + "_median"] = med
        data[kind + "_snr3"] = (np.sum(snr >= 8, axis=1) >= 3) & valid
        data[kind + "_valid"] = valid
    # The predeclared numeric cuts select materially strong fluxes; four
    # fields are directional views (not 4 statistically independent events).
    data["joint"] = (
        data["ion_snr3"] & data["electron_snr3"]
        & data["ion_median"].ge(400) & data["electron_median"].ge(1000)
    )
    return data


def summarize_hours(data: pd.DataFrame) -> pd.DataFrame:
    d = data.copy()
    d["hour"] = d["time"].dt.floor("h")
    by = d.groupby("hour", sort=True)
    h = by.agg(
        samples=("joint", "size"),
        joint_samples=("joint", "sum"),
        ion_median=("ion_median", "median"),
        electron_median=("electron_median", "median"),
        ion_snr3_samples=("ion_snr3", "sum"),
        electron_snr3_samples=("electron_snr3", "sum"),
        valid_ion=("ion_valid", "sum"),
        valid_electron=("electron_valid", "sum"),
        highalt_samples=("highalt", "sum"),
        joint_highalt=("joint_highalt", "sum"),
        distinct_orbits=("spice_orbit_number", "nunique"),
    ).reset_index()
    h["joint_fraction"] = h["joint_samples"] / h["samples"]
    h["event_hour"] = (
        h["samples"].ge(100)
        & h["joint_samples"].ge(75)
        & h["joint_fraction"].ge(0.6)
        & h["joint_highalt"].ge(20)
    )
    return h


def find_events(hourly: pd.DataFrame) -> list[dict]:
    """Merge selected hours if <=5 h apart (Mars-orbit visibility gaps)."""
    selected = hourly.loc[hourly["event_hour"]].sort_values("hour")
    if selected.empty:
        return []
    groups: list[list] = []
    for row in selected.itertuples(index=False):
        if not groups or (row.hour - groups[-1][-1].hour).total_seconds() > 5 * 3600:
            groups.append([row])
        else:
            groups[-1].append(row)
    events = []
    for group in groups:
        if len(group) < 2:  # Require sustained, not a single selected clock hour.
            continue
        sub = pd.DataFrame([r._asdict() for r in group])
        events.append({
            "start_utc": group[0].hour.isoformat(),
            "last_selected_hour_utc": group[-1].hour.isoformat(),
            "selected_hours": len(group),
            "joint_8s_samples": int(sub["joint_samples"].sum()),
            "joint_highalt_8s_samples": int(sub["joint_highalt"].sum()),
            "peak_hourly_median_ion": float(sub["ion_median"].max()),
            "peak_hourly_median_electron": float(sub["electron_median"].max()),
            "median_hourly_ion": float(sub["ion_median"].median()),
            "median_hourly_electron": float(sub["electron_median"].median()),
        })
    events.sort(key=lambda e: (e["joint_highalt_8s_samples"], e["peak_hourly_median_ion"]), reverse=True)
    return events


def summarize_independent_hourly(df: pd.DataFrame, value: str) -> dict:
    tmp = df[["time", value]].copy()
    tmp[value] = pd.to_numeric(tmp[value], errors="coerce")
    tmp = tmp[np.isfinite(tmp[value])]
    if value == "swia_dynamic_pressure":
        tmp = tmp[tmp[value].ge(0)]
    by = tmp.groupby(tmp["time"].dt.floor("h"))[value].median()
    return {x.isoformat(): float(y) for x, y in by.items() if math.isfinite(y)}


def run(output: Path, daily_output: Path) -> dict:
    input_records: dict[str, dict] = {}
    with tempfile.TemporaryDirectory(prefix="openspace-maven-") as root:
        rootpath = Path(root)

        def load(instrument: str, cols: list[str]) -> pd.DataFrame:
            path = rootpath / (instrument + ".parquet")
            input_records[instrument] = download_parquet(instrument, path)
            return read_table(path, cols)

        sep_columns = ["time"] + [v for c in ION + ELECTRON for v in (c, c + "_quality")]
        sep = load("sep", sep_columns)
        s = sep_metrics(sep)
        del sep

        # Must use the spacecraft geometry to control orbital effects.
        spice = load("spice", ["time", "spice_altitude", "spice_orbit_number", "spice_mso_x"])
        s = pd.merge_asof(s, spice, on="time", direction="nearest", tolerance=pd.Timedelta("9s"))
        s["highalt"] = s["spice_altitude"].ge(3000)
        s["joint_highalt"] = s["highalt"] & s["joint"]
        h = summarize_hours(s)
        events = find_events(h)
        first = s["time"].iloc[0].isoformat()
        last = s["time"].iloc[-1].isoformat()
        day = h.copy()
        day["date"] = day["hour"].dt.strftime("%Y-%m-%d")
        daily = day.groupby("date", sort=True).agg(
            observed_hours=("hour", "size"),
            event_hours=("event_hour", "sum"),
            n_8s_samples=("samples", "sum"),
            joint_8s_samples=("joint_samples", "sum"),
            joint_highalt_8s_samples=("joint_highalt", "sum"),
            max_hourly_median_ion=("ion_median", "max"),
            max_hourly_median_electron=("electron_median", "max"),
        ).reset_index()
        daily_output.parent.mkdir(parents=True, exist_ok=True)
        daily.to_csv(daily_output, index=False, float_format="%.6g")

        # Independent environmental context, not automatic confirmation:
        # magnetic field and solar-wind pressure are geometry dependent.
        mag = load("mag", ["time", "mag_mso_x", "mag_mso_y", "mag_mso_z"])
        magvec = mag[["mag_mso_x", "mag_mso_y", "mag_mso_z"]].apply(pd.to_numeric, errors="coerce")
        mag["field_magnitude_nt"] = np.sqrt((magvec ** 2).sum(axis=1, min_count=3))
        mag_context = summarize_independent_hourly(mag, "field_magnitude_nt")
        del mag

        swia = load("swia", ["time", "swia_dynamic_pressure", "swia_dynamic_pressure_quality"])
        # SWIA quality flag is nonzero in some valid samples; keep both
        # analyses descriptive without pretending it means a SEP detection.
        swia_context = summarize_independent_hourly(swia, "swia_dynamic_pressure")
        del swia
        hr = h.loc[h["samples"].ge(100)]
        def median(series: pd.Series) -> float | None:
            x = series.dropna()
            return float(x.median()) if len(x) else None
        context = []
        for event in events[:20]:
            start = pd.Timestamp(event["start_utc"])
            end = pd.Timestamp(event["last_selected_hour_utc"])
            window = pd.date_range(start, end, freq="h", tz="UTC")
            mag_vals = [mag_context.get(t.isoformat()) for t in window if t.isoformat() in mag_context]
            swia_vals = [swia_context.get(t.isoformat()) for t in window if t.isoformat() in swia_context]
            selected_s = s[(s["time"] >= start) & (s["time"] < end + pd.Timedelta(hours=1)) & s["joint"]]
            event["distinct_joint_orbits"] = int(selected_s["spice_orbit_number"].nunique())
            event["mag_hourly_median_nt"] = float(np.median(mag_vals)) if mag_vals else None
            event["swia_hourly_median_dynamic_pressure_npa"] = float(np.median(swia_vals)) if swia_vals else None
            event["complete_independent_mag_hours"] = len(mag_vals)
            event["complete_independent_swia_hours"] = len(swia_vals)
            context.append(event)

        output_data = {
            "status": "instrument-space screening; no claim of novel solar eruption, radiation dose or MAVEN failure cause",
            "input_product": "2026-updated third-party Parquet mirror of the original NASA LASP MAVEN KP in-situ product; NOT the separately released 2026-10-30 PDS Level-2 SEP product",
            "independent_followup_archive": "urn:nasa:pds:maven.sep.calibrated / Release 45 / 2026-10-30; not used in this computation",
            "input": input_records,
            "original_kp_archive": SOURCE_KP,
            "original_instrument_sis": SOURCE_SIS,
            "flux_units": FLUX_UNITS,
            "sensor_energy_windows": {"ions": "30-1000 keV", "electrons": "30-300 keV"},
            "field_quality_semantics": "SEP *_quality stores standard Poisson flux uncertainties, NOT binary quality flags; see PDS SIS/SPASE",
            "time_convention": "UTC; hourly aggregation by floor UTC, no implicit interpolation across missing hours",
            "sample_cadence": "time-indexed original KP at nominal ~4-8s; no assumption of complete coverage",
            "selection_frozen_before_full_run": {
                "valid_views": ">=3/4 for both species",
                "signal_to_uncertainty_per_view": ">=8, at least 3/4 views of each species",
                "row_flux": {"ion_median_min": 400, "electron_median_min": 1000},
                "event_hour": ">=100 source samples; >=75 joint; joint fraction>=0.60; >=20 joint at altitude>=3000km",
                "event": ">=2 selected UTC hours linked by <=5h gaps; at least 2 distinct orbits to qualify for independent-orbit evidence",
            },
            "coverage": {
                "sep_distinct_times": int(len(s)),
                "start_utc": first,
                "end_utc": last,
                "distinct_spice_matched": int(s["spice_altitude"].notna().sum()),
                "hours_with_100_samples": int(len(hr)),
                "selected_hours": int(h["event_hour"].sum()),
                "days_observed": int(len(daily)),
                "dual_species_joint_rows": int(s["joint"].sum()),
                "dual_species_joint_rows_highalt": int(s["joint_highalt"].sum()),
            },
            "survey_hourly_reference": {
                "median_ion_flux_at_100plus_sample_hours": median(hr["ion_median"]),
                "median_electron_flux_at_100plus_sample_hours": median(hr["electron_median"]),
                "median_mag_magnitude_nt_across_valid_hours": float(np.median(list(mag_context.values()))) if mag_context else None,
                "median_swia_dynamic_pressure_npa_across_valid_hours": float(np.median(list(swia_context.values()))) if swia_context else None,
            },
            "ranked_events": context,
            "top_observed_hourly_fluxes": [
                {
                    "hour_utc": v.hour.isoformat(),
                    "ion_median": float(v.ion_median),
                    "electron_median": float(v.electron_median) if math.isfinite(v.electron_median) else None,
                    "joint_samples": int(v.joint_samples),
                    "joint_highalt": int(v.joint_highalt),
                    "samples": int(v.samples),
                }
                for v in h.dropna(subset=["ion_median"]).nlargest(16, "ion_median").itertuples(index=False)
            ],
            "essential_limits": [
                "Four directional look channels are not four independently calibrated instruments.",
                "The file is a third-party Parquet repackaging of NASA LASP KP. Verify candidate timings against original daily KP and independently calibrated SEP L2 files.",
                "MAG and SWIA are environmental context, not independent confirmation of SEP flux.",
                "High-altitude does not ensure upstream solar wind. Analyze spacecraft viewing geometry and solar magnetic connection before claiming interplanetary impact.",
                "Multiple-hypothesis top-hour selection is exploratory; no survey false-alarm rate calibrated yet.",
                "MAVEN spacecraft loss in Dec 2025 cannot be connected causally to earlier SEP events without independent evidence.",
                "SEP *_quality suffix is incorrectly called a binary quality flag in the third-party mirror README; the PDS specifies a one-sigma Poisson statistical uncertainty.",
            ],
        }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(output_data, indent=2, allow_nan=False) + "\n")
    return output_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="results/maven-2025-sep-events/summary.json", type=Path)
    parser.add_argument("--daily", default="results/maven-2025-sep-events/daily.csv", type=Path)
    args = parser.parse_args()
    result = run(args.output, args.daily)
    print(json.dumps({"coverage": result["coverage"], "top_events": result["ranked_events"][:3]}, indent=2))
