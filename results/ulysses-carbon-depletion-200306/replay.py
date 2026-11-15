"""Screen the full Ulysses SWICS ratio series, then test the selected episodes.

From this directory: python replay.py [--fetch DIR]. The fetch option reconstructs
the compact reductions from the manifest, verifying the original response bytes.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tempfile
from urllib.request import urlopen
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
DATE = ["year", "day", "hour", "minute", "second"]
RATIO = ["speed", "C6_C5", "O7_O6", "Fe_O", "QFe", "cycles"]
SPECIES = ["He", "C", "O", "Ne", "Mg", "Si", "S", "Fe"]
MOMENT = [f"{species}_{key}" for species in SPECIES for key in ["den", "v", "temp"]]


def get_members(fetch_dir):
    manifest = json.loads((ROOT / "sources.json").read_text())
    fetch_dir.mkdir(parents=True, exist_ok=True)

    def read_entry(entry):
        name = entry["member"]
        file = fetch_dir / name
        if file.exists():
            raw = file.read_bytes()
        else:
            with urlopen(entry["url"], timeout=120) as response:
                raw = response.read()
        if hashlib.sha256(raw).hexdigest() != entry["sha256"]:
            raise ValueError(f"Checksum mismatch: {name}")
        if not file.exists():
            file.write_bytes(raw)
        return name, raw

    with ThreadPoolExecutor(max_workers=6) as pool:
        return dict(pool.map(read_entry, manifest))


def read_swics(raw, columns):
    f = pd.read_csv(io.BytesIO(raw), sep=r"\s+", names=DATE + columns)
    f["utc"] = pd.to_datetime(f.year.astype(str) + "-01-01", utc=True)
    for key, unit in [("day", "D"), ("hour", "h"), ("minute", "m"), ("second", "s")]:
        f["utc"] += pd.to_timedelta(f[key] - int(key == "day"), unit=unit)
    return f.drop(columns=DATE)


def contiguous_runs(series, cutoff=1, cycles=9):
    good = ((series.cycles >= cycles) & (series.C6_C5 > 0)
            & (series.O7_O6 > 0) & (series.QFe > 0))
    flag = good & (series.O7_O6 > 0.15) & (series.contrast < cutoff)
    breaks = (flag != flag.shift()) | (series.utc.diff() > pd.Timedelta(hours=4.5))
    runs = series.loc[flag].groupby(breaks.cumsum()[flag]).agg(
        first=("utc", "first"), last=("utc", "last"), samples=("utc", "size"),
        median_contrast=("contrast", "median"))
    return good, flag, runs.sort_values(["samples", "first"], ascending=[False, True])


def medians(frame, columns):
    return {k: round(float(frame[k].median()), 5) for k in columns}


def matrix_rates(members):
    import cdflib

    # The native CDF contains uint32 count vectors. The HAPI conversion of
    # these variables can return empty fields even when the CDF has counts.
    frames = []
    with tempfile.TemporaryDirectory() as scratch:
        for name, raw in members.items():
            if not name.endswith(".cdf"):
                continue
            path = Path(scratch) / name
            path.write_bytes(raw)
            cdf = cdflib.CDF(path)
            stamp = pd.to_datetime(cdflib.cdfepoch.to_datetime(cdf.varget("Epoch")), utc=True)
            frame = pd.DataFrame({"utc": stamp, "nominal": cdf.varget("NON_NOMINAL_FLAG") == 0,
                                  "voltage": cdf.varget("PAPS_VOLTAGE")})
            for number in (5, 6, 8):
                variable = f"MATRIX_RATE_{number}"
                values = cdf.varget(variable).astype(float)
                limits = cdf.varattsget(variable)
                valid = ((values >= limits["VALIDMIN"]) & (values <= limits["VALIDMAX"])).all(axis=1)
                frame[f"ch{number}"] = np.where(valid, values.sum(axis=1), np.nan)
            frames.append(frame)
    scans = pd.concat(frames)
    scans = scans[scans.nominal & (scans.voltage > 0)].dropna()
    # The ratio product uses a three-hour cadence and labels each bin with
    # its first actual scan, which may fall after the UTC boundary.
    scans["window"] = scans.utc.dt.floor("3h")
    return scans.groupby("window").agg(
        scans=("utc", "size"), ch5=("ch5", "sum"), ch6=("ch6", "sum"), ch8=("ch8", "sum"),
        voltage_min=("voltage", "min"), voltage_max=("voltage", "max"))


def summarize_matrix(rates, stamps):
    windows = rates.reindex(stamps.dt.floor("3h"))
    windows = windows[(windows.scans >= 9) & (windows.ch8 > 0)]
    return {"bins": len(windows), "scans": int(windows.scans.sum()),
            "nominal_C6_over_O6": round(float((windows.ch6/windows.ch8).median()), 5),
            "nominal_C45_over_O6": round(float((windows.ch5/windows.ch8).median()), 5),
            "ch5_counts": int(windows.ch5.sum()), "ch6_counts": int(windows.ch6.sum()),
            "ch8_counts": int(windows.ch8.sum()), "voltage_min_kV": round(float(windows.voltage_min.min()), 2),
            "voltage_max_kV": round(float(windows.voltage_max.max()), 2)}


def write_reductions(members):
    # Keep only the samples and channels needed by the two tests. Integer
    # encoding is lossless at the precision of the published charge ratios.
    ratios = pd.concat([read_swics(raw, RATIO) for name, raw in members.items()
                        if name.startswith("uswichst")], ignore_index=True)
    reduced = pd.DataFrame({
        "utc_seconds": (ratios.utc - pd.Timestamp("1970-01-01", tz="UTC")).dt.total_seconds().astype("int64"),
        "carbon_x1000": (ratios.C6_C5 * 1000).round().astype("int32"),
        "oxygen_x1000": (ratios.O7_O6 * 1000).round().astype("int32"),
        "iron_valid": (ratios.QFe > 0).astype("int8"), "cycles": ratios.cycles})
    reduced.to_csv(ROOT / "mission_ratio_inputs.csv.gz", index=False,
                   compression={"method":"gzip", "mtime":0})
    month = ratios[ratios.utc.dt.strftime("%Y-%m") == "2003-06"].merge(
        read_swics(members["swi03152.asc"], MOMENT)[["utc", "C_den"]], on="utc", validate="one_to_one")
    month.to_csv(ROOT / "june_swics.csv", index=False)
    plasma = pd.read_csv(io.BytesIO(members["june_coho.csv"]),
                         names="utc BR BT BN B speed density alpha_density temperature".split())
    plasma[["utc", "B", "speed", "density", "alpha_density"]].to_csv(ROOT / "june_plasma.csv", index=False)
    catalog = pd.read_csv(io.BytesIO(members["icmecat.csv"]))
    catalog.loc[catalog.sc_insitu.eq("ULYSSES") & catalog.mo_start_time.str.startswith("2003-06"),
                ["icmecat_id", "mo_start_time", "mo_end_time", "mo_sc_heliodistance"]].to_csv(
                    ROOT / "june_ejecta.csv", index=False)
    matrix_rates(members).to_csv(ROOT / "june_matrix_windows.csv", index_label="window", float_format="%.4f")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", metavar="DIR", type=Path)
    args = parser.parse_args()
    if args.fetch:
        write_reductions(get_members(args.fetch))
    compact = pd.read_csv(ROOT / "mission_ratio_inputs.csv.gz")
    ratios = pd.DataFrame({"utc": pd.to_datetime(compact.utc_seconds, unit="s", utc=True),
                           "C6_C5": compact.carbon_x1000 / 1000,
                           "O7_O6": compact.oxygen_x1000 / 1000,
                           "QFe": compact.iron_valid, "cycles": compact.cycles})
    # Multiple versions of 11 timestamps occur outside the selected intervals.
    # Exclude ambiguous timestamps from the comparison instead of choosing a version.
    duplicates = ratios.utc.duplicated(keep=False)
    series = ratios.loc[~duplicates].sort_values("utc").reset_index(drop=True)
    series["contrast"] = series.C6_C5 / series.O7_O6.where(series.O7_O6 > 0)
    good, flag, runs = contiguous_runs(series)
    peaks = runs.iloc[:2]
    summary = {
        "source_rows": len(ratios), "duplicate_timestamp_rows_excluded": int(duplicates.sum()),
        "qualifying_samples": int(good.sum()), "anomalous_samples": int(flag.sum()),
        "mission_median_contrast": round(float(series.loc[good, "contrast"].median()), 4),
        "run_counts_by_length": {str(k): int(v) for k, v in runs.samples.value_counts().sort_index().items()},
        "episodes": [], "sensitivity": []
    }
    all_runs = runs.copy()
    all_runs.to_csv(ROOT / "screened_runs.csv", index=False, float_format="%.6f")
    catalog = pd.read_csv(ROOT / "june_ejecta.csv")
    catalog["start"] = pd.to_datetime(catalog.mo_start_time, utc=True)
    catalog["end"] = pd.to_datetime(catalog.mo_end_time, utc=True)
    june = pd.read_csv(ROOT / "june_swics.csv", parse_dates=["utc"])
    june["contrast"] = june.C6_C5 / june.O7_O6.where(june.O7_O6 > 0)
    # Both products use overlapping SWICS counts. Zero indicates missing.
    june["C6_O6"] = june.C_den.where(june.C_den > 0)
    june["C5_O6_inferred"] = june.C6_O6 / june.C6_C5.where(june.C6_C5 > 0)
    plasma = pd.read_csv(ROOT / "june_plasma.csv", parse_dates=["utc"])
    plasma = plasma.set_index("utc").mask(lambda x: abs(x) > 1e10)
    plasma["alpha_proton"] = plasma.alpha_density / plasma.density.where(plasma.density > 0)
    inside = pd.Series(False, index=june.index)
    for row in peaks.itertuples():
        inside |= june.utc.between(row.first, row.last)
    control = june[~inside & (june.O7_O6 >= 0.5) & (june.QFe >= 12)
                   & (june.speed >= 500) & (june.cycles >= 9)
                   & june.C6_O6.notna() & june.C5_O6_inferred.notna()]
    summary["june_high_charge_fast_control"] = {
        "count": len(control), **medians(control, ["C6_O6", "C5_O6_inferred"]),
        "min_C6_O6": round(float(control.C6_O6.min()), 5)}
    # Use a 3-hour window for each ratio timestamp. These are approximate
    # support intervals since the original accumulation lengths vary.
    control_hours = pd.Series(False, index=plasma.index)
    for stamp in control.utc:
        control_hours |= (plasma.index >= stamp) & (plasma.index < stamp + pd.Timedelta(hours=3))
    summary["alpha_proton_controls"] = {
        "june_valid_hours": int(plasma.alpha_proton.count()),
        "june_median": round(float(plasma.alpha_proton.median()), 5),
        "high_charge_fast_control_valid_hours": int(plasma.loc[control_hours, "alpha_proton"].count()),
        "high_charge_fast_control_median": round(float(plasma.loc[control_hours, "alpha_proton"].median()), 5)}
    rates = pd.read_csv(ROOT / "june_matrix_windows.csv", parse_dates=["window"]).set_index("window")
    summary["june_matrix_control"] = summarize_matrix(rates, control.utc)
    measurement_rows = []
    for row in peaks.itertuples():
        part = june[june.utc.between(row.first, row.last)]
        cat_rows = catalog[(catalog.start <= row.first) & (catalog.end >= row.last)]
        if len(cat_rows) != 1:
            raise ValueError("No unique enclosing catalog interval")
        event = cat_rows.iloc[0]
        stop = row.last + pd.Timedelta(hours=3)
        flow = plasma[(plasma.index >= row.first) & (plasma.index < stop)]
        same_icme = june[(june.utc >= event.start) & (june.utc < event.end)
                         & ~inside & june.C6_O6.notna()]
        summary["episodes"].append({
            "first_utc": row.first.isoformat(), "last_utc": row.last.isoformat(),
            "samples": row.samples, "catalog_id": event.icmecat_id,
            "catalog_start": event.mo_start_time, "catalog_end": event.mo_end_time,
            "distance_AU": float(event.mo_sc_heliodistance),
            "matrix_counts": summarize_matrix(rates, part.utc),
            **medians(part, ["C6_C5", "O7_O6", "QFe", "contrast", "C6_O6", "C5_O6_inferred"]),
            "control_over_episode_C6_O6": round(float(control.C6_O6.median()/part.C6_O6.median()), 2),
            "same_icme_other_valid_C6_O6_count": len(same_icme),
            "same_icme_other_median_C6_O6": round(float(same_icme.C6_O6.median()), 5),
            "coho_hours": len(flow), "alpha_proton_valid_hours": int(flow.alpha_proton.count()),
            "alpha_proton": round(float(flow.alpha_proton.median()), 5),
            "median_proton_speed_km_s": round(float(flow.speed.median()), 2),
            "median_field_nT": round(float(flow.B.median()), 3)})
        measurement_rows.extend(part.index)
    june.loc[measurement_rows, ["utc", *RATIO, "contrast", "C6_O6", "C5_O6_inferred"]].to_csv(
        ROOT / "episode_samples.csv", index=False, float_format="%.6f")
    for cutoff in (0.3, 0.5, 0.75, 1, 1.5, 2):
        for cycles in (9, 11):
            _, _, tests = contiguous_runs(series, cutoff, cycles)
            e_lengths = []
            for peak in peaks.itertuples():
                vicinity = tests[(tests["last"] >= peak.first) & (tests["first"] <= peak.last)]
                e_lengths.append(int(vicinity.samples.max()) if len(vicinity) else 0)
            summary["sensitivity"].append({
                "cutoff": cutoff, "min_cycles": cycles, "mission_max": int(tests.samples.max()),
                "june_run_lengths": e_lengths})
    (ROOT / "measurements.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
