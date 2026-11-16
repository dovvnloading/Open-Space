#!/usr/bin/env python3
"""Reproducible, frame-aware comparison of NASA PDS Mars gravity models.

The input files are archived JPL SHADR tables, not new spacecraft measurements.
No third-party packages. This tool does not assign p-values to correlated fits.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
import re
import time
import urllib.request

ARCHIVE = "https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/"
FILES = {"MRO120F": "jgmro_120f_sha.tab", "MRO120G": "jgmro_120g_sha.tab"}
LABELS = {"MRO120F": "jgmro_120f_sha.lbl", "MRO120G": "jgmro_120g_sha.lbl"}
ARCSEC_PER_RAD = 180.0 * 3600.0 / math.pi
RADIUS_M = 3396000.0
BANDS = [(2, 5), (6, 10), (11, 20), (21, 40), (41, 60), (61, 79), (80, 120)]
K = ["C", "S", "sigmaC", "sigmaS"]


def download(url: str, path: pathlib.Path) -> bytes:
    if path.is_file():
        return path.read_bytes()
    error = None
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Open-Space-Mars-Gravity-Audit/1.0"})
            with urllib.request.urlopen(req, timeout=55) as response:
                content = response.read()
            if not content:
                raise RuntimeError("Empty archive response: " + url)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
            return content
        except Exception as exc:
            error = exc
            time.sleep(attempt + 1)
    raise RuntimeError("Could not fetch %s: %s" % (url, error))


def parse_table(content: bytes) -> dict:
    text = content.decode("ascii")
    header = None
    coeff = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = [x.strip().replace("D", "E") for x in line.split(",")]
        if len(parts) == 8 and header is None:
            header = [float(parts[0]), float(parts[1]), float(parts[2]), *[int(x) for x in parts[3:6]]]
            continue
        if len(parts) != 6:
            raise ValueError("Bad SHADR row: " + repr(line[:150]))
        n, m = int(parts[0]), int(parts[1])
        assert 1 <= n <= 120 and 0 <= m <= n, (n, m)
        assert (n, m) not in coeff, ("Duplicate SHADR row", n, m)
        coeff[n, m] = tuple(map(float, parts[2:]))
    assert header is not None and header[0] == 3396.0, header
    assert header[3:6] == [120, 120, 1], header
    expected = {(n, m) for n in range(1, 121) for m in range(n + 1)}
    assert len(coeff) == 7380 and set(coeff) == expected, (len(coeff), len(expected - set(coeff)))
    assert all(coeff[n, 0][1] == 0 for n in range(1, 121))
    return {"reference_radius_km": header[0], "GM_km3_s2": header[1],
            "GM_sigma_km3_s2": header[2], "coefficients": coeff}


def validate_label(data: bytes, model: str) -> None:
    label = data.decode("ascii", errors="replace").upper()
    assert ("ORIGINAL_PRODUCT_ID" in label and model in label), model
    assert re.search(r"ROWS\s*=\s*7380\b", label), model
    assert "FULLY NORMALIZED" in label, model
    assert "REFERENCE RADIUS = 3396.0 KM" in label, model


def fit_yaw(a: dict, b: dict, lo: int, hi: int, equal_degree: bool = False) -> float:
    """Gauss-Newton least-squares yaw in radians; avoid falsely treating model errors as independent."""
    weights = {}
    for n in range(lo, hi + 1):
        mag = sum((m*m) * (a[n, m][0]**2 + a[n, m][1]**2) for m in range(1, n + 1))
        weights[n] = (1.0 / mag if equal_degree and mag else 1.0)
    yaw = 0.0
    for _ in range(8):
        numerator = denominator = 0.0
        for n in range(lo, hi + 1):
            for m in range(1, n + 1):
                ca, sa = a[n, m][:2]
                cb, sb = b[n, m][:2]
                z = complex(ca, sa) * complex(math.cos(m*yaw), math.sin(m*yaw))
                diff = complex(cb, sb) - z
                deriv = 1j * m * z
                w = weights[n]
                numerator += w * (diff * deriv.conjugate()).real
                denominator += w * (abs(deriv)**2)
        assert denominator > 0
        step = numerator / denominator
        yaw += step
        if abs(step) < 1e-14:
            break
    return yaw


def band_stats(a: dict, b: dict, lo: int, hi: int, yaw: float, gm: float, altitude_km: float) -> dict:
    raw_power = aligned_power = a_power = b_power = 0.0
    graw2 = galigned2 = 0.0
    r_km = 3396.0 + altitude_km
    base = gm / (r_km*r_km) * 1000.0  # km/s^2 -> m/s^2
    n_rows = 0
    for n in range(lo, hi + 1):
        factor = base * (n + 1) * (3396.0/r_km)**n
        for m in range(n + 1):
            ca, sa = a[n, m][:2]
            cb, sb = b[n, m][:2]
            old = complex(ca, sa)
            new = complex(cb, sb)
            rotated = old * complex(math.cos(m*yaw), math.sin(m*yaw))
            raw = abs(new - old)**2
            residual = abs(new - rotated)**2
            raw_power += raw
            aligned_power += residual
            a_power += abs(old)**2
            b_power += abs(new)**2
            graw2 += raw * factor * factor
            galigned2 += residual * factor * factor
            n_rows += 1
    return {
        "degrees": [lo, hi], "coefficient_rows": n_rows,
        "raw_coefficient_l2": math.sqrt(raw_power),
        "aligned_coefficient_l2": math.sqrt(aligned_power),
        "fraction_raw_power_explained_by_rotation": 1.0 - aligned_power / raw_power if raw_power else None,
        "degree_band_norm_change_fraction": math.sqrt(b_power/a_power) - 1.0 if a_power else None,
        "radial_acceleration_rms_at_%skm_m_s2" % int(altitude_km): {
            "raw": math.sqrt(graw2), "aligned": math.sqrt(galigned2)
        },
    }


def run(cache: pathlib.Path, out: pathlib.Path, altitude_km: float = 400.0) -> dict:
    tables = {}
    receipts = {}
    for model in FILES:
        path = cache / FILES[model]
        content = download(ARCHIVE + FILES[model], path)
        label_bytes = download(ARCHIVE + LABELS[model], cache / LABELS[model])
        validate_label(label_bytes, model)
        tables[model] = parse_table(content)
        receipts[model] = {
            "table_url": ARCHIVE + FILES[model], "label_url": ARCHIVE + LABELS[model],
            "sha256_table": hashlib.sha256(content).hexdigest(),
            "sha256_label": hashlib.sha256(label_bytes).hexdigest(),
            "table_bytes": len(content), "label_bytes": len(label_bytes),
            "rows": len(tables[model]["coefficients"]),
            "GM_km3_s2": tables[model]["GM_km3_s2"],
            "GM_sigma_km3_s2": tables[model]["GM_sigma_km3_s2"]
        }
    f = tables["MRO120F"]["coefficients"]
    g = tables["MRO120G"]["coefficients"]
    fit = fit_yaw(f, g, 2, 10)
    ex2 = fit_yaw(f, g, 3, 10)
    equal = fit_yaw(f, g, 2, 10, True)
    controls = {}
    for name, (lo, hi) in {"train": (2, 10), "holdout_low": (11, 20),
                           "holdout_mid": (21, 40)}.items():
        controls[name] = {
            "raw": band_stats(f, g, lo, hi, 0.0, tables["MRO120G"]["GM_km3_s2"], altitude_km),
            "fit_heldout_arcsec": fit_yaw(f, g, lo, hi) * ARCSEC_PER_RAD,
            "rotated": band_stats(f, g, lo, hi, fit, tables["MRO120G"]["GM_km3_s2"], altitude_km)
        }
    model_delta_gm = tables["MRO120G"]["GM_km3_s2"] - tables["MRO120F"]["GM_km3_s2"]
    g_mono = abs(model_delta_gm) / ((3396.0 + altitude_km)**2) * 1000.0
    report = {
        "description": "Two correlated model estimates, not independent new Mars observations",
        "data_release": "2026-09-24",
        "source_publication": "Konopliv et al. 2026, doi:10.3847/PSJ/ae9ba7",
        "model_hashes_and_metadata": receipts,
        "source_coefficients": 7380,
        "fit_training_degrees": [2, 10],
        "effective_yaw": {"radians": fit, "arcseconds": fit * ARCSEC_PER_RAD,
                          "equatorial_distance_m": fit * RADIUS_M,
                          "fit_excluding_degree_2_arcsec": ex2 * ARCSEC_PER_RAD,
                          "fit_equal_degree_weights_arcsec": equal * ARCSEC_PER_RAD},
        "controls": controls,
        "degree_bands": [band_stats(f, g, lo, hi, fit, tables["MRO120G"]["GM_km3_s2"], altitude_km)
                         for lo, hi in BANDS],
        "GM_difference": {"G_minus_F_km3_s2": model_delta_gm,
                          "monopole_radial_accel_difference_at_%skm_m_s2" % int(altitude_km): g_mono},
        "limitations": [
            "The fitted yaw is an effective longitude-registration angle, not proof of a physical rotation of Mars.",
            "Both fields fit mostly overlapping tracking measurements; quoted coefficient sigmas are correlated and cannot establish a detection significance for their difference.",
            "Comparison only accounts for a z-axis frame rotation, not polar-axis, epoch, tide, or force-model changes.",
            "Low-degree fit is frozen before inspecting held-out 11-40 and 41-120 degrees.",
            "This is a reanalysis of published model products; do not claim independent discovery of decadal spin changes documented in the primary paper.",
            "Monopole GM changes are reported separately and are not removed by longitude alignment.",
            "No raw Doppler residuals, independent holdout tracking, or mission-specific orbit propagation were tested."
        ]
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"angle_arcsec": report["effective_yaw"]["arcseconds"],
        "train_fraction_explained": controls["train"]["rotated"]["fraction_raw_power_explained_by_rotation"],
        "holdout_fraction_explained": controls["holdout_low"]["rotated"]["fraction_raw_power_explained_by_rotation"],
        "output": str(out)}, indent=2))
    return report


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--cache", type=pathlib.Path, default=pathlib.Path("/tmp/open-space-mro120g-pds"))
    p.add_argument("--output", type=pathlib.Path, default=pathlib.Path("/tmp/mars-gravity-mro120g-summary.json"))
    p.add_argument("--altitude-km", type=float, default=400.)
    a = p.parse_args()
    run(a.cache, a.output, a.altitude_km)


if __name__ == "__main__":
    main()
