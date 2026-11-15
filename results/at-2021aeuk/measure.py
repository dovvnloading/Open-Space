"""Replay the Swift measurements with Python 3.10+, NumPy and SciPy."""
import csv
import json
from pathlib import Path
from statistics import median

from scipy.stats import chi2, gamma, poisson

ROOT = Path(__file__).resolve().parent


def summarize_xrt(rows):
    n = sum(int(row["source_counts"]) for row in rows)
    background = sum(float(row["expected_background"]) for row in rows)
    effective = sum(float(row["exposure_s"]) / float(row["correction"]) for row in rows)
    return {
        "visits": len(rows), "counts": n, "background": background, "effective_s": effective,
        "net_rate": (n - background) / effective,
        "source_rate_lower_97p5": max(0, chi2.ppf(.025, 2*n)/2 - background) / effective,
        "all_counts_rate_upper_97p5": chi2.ppf(.975, 2*(n+1))/2 / effective,
        "flat_prior_source_upper_95": (gamma.isf(.05*gamma.sf(background, n+1), n+1) - background) / effective,
        "background_only_p_ge_counts": poisson.sf(n-1, background),
    }


def summarize_uv(products):
    result = {}
    for band in ("UVW1", "UVM2", "UVW2"):
        selected = [row for item in products if item["mode"] == "UVOT"
                    for row in item["rows"] if row["band"] == band and row["flags"] == 0
                    and row["sep"] <= 2 and row["magerr"] < .2]
        if selected:
            mags = [row["mag"] for row in selected]
            result[band] = dict(visits=len(mags), median_mag=median(mags), range_mag=[min(mags), max(mags)])
    return result


def aperture_check(products, radius, require_clean_center):
    exposures = {item["obs"]: item["rows"][0] for item in products if item["mode"] == "EXPO"}
    selected = []
    for item in products:
        if item["mode"] != "XRT":
            continue
        row = item["rows"][0]
        if require_clean_center and exposures[item["obs"]]["center"] < .8*row["exp"]:
            continue
        selected.append(row)
    on = sum(row[f"n{radius}"] for row in selected)
    off = sum(row["annulus"] for row in selected)
    t = sum(row["exp"] for row in selected)
    return dict(visits=len(selected), counts=on, annulus=off, exposure_s=t,
                raw_net_rate=(on - off*radius**2/(60**2 - 30**2))/t)


def measure(rows, packet):
    total = [r for r in rows if r["band"] == "Total"]
    if len({r["obsid"] for r in total}) != len(total):
        raise ValueError("An observation appears twice")
    windows = {"2025": (60800, 61000), "2025_last_six": (60898, 61000),
               "2026": (61200, 61400), "published_2023_low_state": (59924,60045)}
    result = {"xrt": {}, "uvot": {}, "apertures": {}}
    for name, (start, end) in windows.items():
        chosen = [r for r in total if start <= float(r["mjd"]) < end]
        result["xrt"][name] = summarize_xrt(chosen)
        if name not in ("2025", "2026"):
            continue
        obsids = {r["obsid"] for r in chosen}
        products = [p for p in packet["products"] if p["obs"] in obsids]
        if {p["obs"] for p in products if p["mode"] == "XRT"} != obsids:
            raise ValueError("Missing XRT event lists")
        result["uvot"][name] = summarize_uv(products)
        result["apertures"][name] = {f"{radius}px_{label}": aperture_check(products,radius,clean)
                                     for radius in (10,15,20) for label,clean in
                                     (("all",False), ("clean_center",True))}
    late = result["xrt"]["2026"]
    result["rate_ratio_joint_lower_95"] = {
        name: result["xrt"][name]["source_rate_lower_97p5"] / late["all_counts_rate_upper_97p5"]
        for name in ("2025", "2025_last_six")
    }
    ids25 = {r["obsid"] for r in total if 60800 <= float(r["mjd"]) < 61000}
    shifts = []
    for star in packet["controls"]:
        early = [m for obs,m in star["visits"] if obs in ids25]
        late_mags = [m for obs,m in star["visits"] if obs not in ids25]
        shifts.append(median(early) - median(late_mags))
    result["control_uvw1_brightening_mag"] = dict(count=len(shifts), median=median(shifts),
                                                   minimum=min(shifts), maximum=max(shifts))
    result["uv_flux_factors"] = {
        band: 10**(.4*(a["median_mag"] - result["uvot"]["2026"][band]["median_mag"]))
        for band,a in result["uvot"]["2025"].items()
    }
    return result


if __name__ == "__main__":
    with (ROOT / "lsxps.csv").open() as stream:
        measurements = measure(list(csv.DictReader(stream)), json.loads((ROOT / "swift-products.json").read_text()))
    expected = json.loads((ROOT / "result.json").read_text())
    def check(actual, reference):
        if isinstance(reference, dict):
            assert actual.keys() == reference.keys()
            for key in reference:
                check(actual[key], reference[key])
        elif isinstance(reference, list):
            assert len(actual) == len(reference)
            for value, test in zip(actual, reference):
                check(value, test)
        else:
            from math import isclose
            assert isclose(actual, reference, rel_tol=1e-8, abs_tol=1e-12), (actual, reference)
    check(measurements, expected)
    print(json.dumps(measurements, indent=2))
