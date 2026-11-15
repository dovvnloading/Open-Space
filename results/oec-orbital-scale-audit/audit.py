#!/usr/bin/env python3
"""Reproduce the independent, frozen 2023 Open Exoplanet Catalogue orbit audit.

No third-party dependencies. This flags inconsistencies; it does not identify a
cause or override original orbital solutions. The separately curated 2026
source-XML checks and primary-paper values live in primary_checks.csv.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

OEC_TABLE_COMMIT = "bf1efad9e4bdc9f0887fe392b77308c3585e263f"
OEC_TABLE_BLOB = "db832c9d70e953ca6e1851122e1d46578c3b5d58"
OEC_SOURCE_XML_COMMIT = "f4228675cfde38cdadc2f3dc93688edb11d89f77"
SOURCE_URL = (
    "https://raw.githubusercontent.com/OpenExoplanetCatalogue/oec_tables/"
    f"{OEC_TABLE_COMMIT}/comma_separated/open_exoplanet_catalogue.txt"
)
DAY_PER_YEAR = 365.25  # Julian year, not calendar-year iteration
SOLAR_RADIUS_IN_AU = 695700000.0 / 149597870700.0  # IAU nominal values
SOLAR_EFFECTIVE_TEMPERATURE_K = 5772.0
PLANET_MASS_SOLAR_PER_JUPITER = 0.000954588
REQUIRED = (
    "name", "period", "semimajoraxis", "hoststar_mass", "binaryflag",
    "hoststar_radius", "discoverymethod", "lastupdate", "list",
)
RESULT_FIELDS = (
    "name", "period_d", "a_reported_au", "host_mass_solar",
    "host_radius_solar", "mass_inferred_solar", "mass_ratio_q",
    "disagreement_factor", "kepler_a_au", "discoverymethod", "catalogue_lastupdate",
)
SOURCE_FIELDS = (
    "planet", "period_d", "paper_a_au", "paper_a_over_rstar",
    "paper_host_mass_solar", "paper_host_radius_solar",
    "mass_implied_by_p_and_a_solar", "mass_ratio_to_paper_host",
    "axis_predicted_by_p_and_host_mass_au", "axis_from_paper_a_over_rstar_au",
    "paper_axis_to_ratio_axis_percent",
)
PRIMARY_FIELDS = (
    "planet", "oec_reported_a_au", "paper_a_over_rstar",
    "oec_host_radius_solar", "paper_flux_earth", "temperature_k",
    "paper_a_over_rstar_times_one_solar_radius_au", "corrected_a_au",
    "oec_to_paper_scale_error_percent", "paper_flux_reconstructed_earth",
    "using_oec_a_flux_earth", "paper_flux_error_percent",
    "apparent_mass_ratio_q",
)


def positive_float(value: str) -> float | None:
    """Reject empty/NaN/inf/nonpositive inputs without changing the denominator."""
    try:
        parsed = float(value)
    except (ValueError, TypeError):
        return None
    return parsed if math.isfinite(parsed) and parsed > 0 else None


def implied_mass_solar(period_days: float, semi_major_au: float) -> float:
    """Keplerian total dynamical mass in solar units, a^3 / P^2."""
    return semi_major_au**3 / (period_days / DAY_PER_YEAR) ** 2


def predicted_axis_au(period_days: float, host_mass_solar: float) -> float:
    """Massless-planet approximation, for a conservative screening statistic."""
    return (host_mass_solar * (period_days / DAY_PER_YEAR) ** 2) ** (1 / 3)


def audit_stream(stream, threshold: float = 1.5):
    """Return sorted flags and transparent eligibility/exclusion denominators.

    One catalogue row = one observing/analysis unit. No cross-record imputation.
    Uncertainties and massive companions are not modeled; every row is a *flag*.
    """
    if not math.isfinite(threshold) or threshold <= 1:
        raise ValueError("threshold must be finite and greater than one")
    reader = csv.DictReader(stream)
    missing = set(REQUIRED) - set(reader.fieldnames or ())
    if missing:
        raise ValueError(f"missing columns: {', '.join(sorted(missing))}")
    counts = Counter()
    seen_names = set()
    flags = []
    for row in reader:
        counts["source_rows"] += 1
        name = (row.get("name") or "").strip()
        if not name:
            counts["excluded_blank_name"] += 1
            continue
        if name in seen_names:
            counts["duplicate_planet_names"] += 1
        seen_names.add(name)
        # A list entry can contain quoted commas, e.g. "Confirmed planets, ...".
        if "Confirmed planets" not in (row.get("list") or "").split(", "):
            counts["excluded_not_confirmed"] += 1
            continue
        if (row.get("binaryflag") or "").strip() != "0":
            counts["excluded_nonzero_binaryflag"] += 1
            continue
        period = positive_float(row.get("period"))
        a = positive_float(row.get("semimajoraxis"))
        host_mass = positive_float(row.get("hoststar_mass"))
        if period is None or a is None or host_mass is None:
            counts["excluded_incomplete_or_invalid"] += 1
            continue
        counts["eligible"] += 1
        mass_dyn = implied_mass_solar(period, a)
        q = mass_dyn / host_mass
        factor = max(q, 1 / q)
        if factor <= threshold:
            continue
        counts["flags_over_threshold"] += 1
        for bound in (2, 3, 5, 10):
            if factor > bound:
                counts[f"flags_over_{bound}x"] += 1
        flags.append({
            "name": name,
            "period_d": f"{period:.10g}",
            "a_reported_au": f"{a:.10g}",
            "host_mass_solar": f"{host_mass:.10g}",
            "host_radius_solar": (row.get("hoststar_radius") or "").strip(),
            "mass_inferred_solar": f"{mass_dyn:.9g}",
            "mass_ratio_q": f"{q:.9g}",
            "disagreement_factor": f"{factor:.8g}",
            "kepler_a_au": f"{predicted_axis_au(period, host_mass):.9g}",
            "discoverymethod": (row.get("discoverymethod") or "").strip(),
            "catalogue_lastupdate": (row.get("lastupdate") or "").strip(),
        })
    flags.sort(key=lambda x: (-float(x["disagreement_factor"]), x["name"]))
    for key in (
        "source_rows", "excluded_blank_name", "excluded_not_confirmed",
        "excluded_nonzero_binaryflag", "excluded_incomplete_or_invalid",
        "duplicate_planet_names", "eligible", "flags_over_threshold",
        "flags_over_2x", "flags_over_3x", "flags_over_5x", "flags_over_10x",
    ):
        counts.setdefault(key, 0)
    return flags, dict(counts)


def primary_paper_diagnostics(stream):
    """Independent test using published *a/R_star* and irradiation, not P/a."""
    reader = csv.DictReader(stream)
    required = {
        "planet", "period_d", "oec_reported_a_au", "oec_host_mass_solar",
        "oec_host_radius_solar", "temperature_k", "paper_a_over_rstar",
        "paper_flux_earth",
    }
    if required - set(reader.fieldnames or ()):
        raise ValueError("missing required original-paper check fields")
    result = []
    for row in reader:
        def field(key):
            v = positive_float(row[key])
            if v is None:
                raise ValueError(f"{row['planet']}: invalid {key}")
            return v
        a_oec = field("oec_reported_a_au")
        radius = field("oec_host_radius_solar")
        a_over_r = field("paper_a_over_rstar")
        temp = field("temperature_k")
        published_irradiation = field("paper_flux_earth")
        p = field("period_d")
        mass = field("oec_host_mass_solar")
        a_one_solar_radius = a_over_r * SOLAR_RADIUS_IN_AU
        corrected_axis = a_one_solar_radius * radius
        lum_solar = radius**2 * (temp / SOLAR_EFFECTIVE_TEMPERATURE_K) ** 4
        reconstructed_flux = lum_solar / corrected_axis**2
        wrong_flux = lum_solar / a_oec**2
        result.append({
            "planet": row["planet"],
            "oec_reported_a_au": f"{a_oec:.9g}",
            "paper_a_over_rstar": f"{a_over_r:.7g}",
            "oec_host_radius_solar": f"{radius:.7g}",
            "paper_flux_earth": f"{published_irradiation:.7g}",
            "temperature_k": f"{temp:.7g}",
            "paper_a_over_rstar_times_one_solar_radius_au": f"{a_one_solar_radius:.9g}",
            "corrected_a_au": f"{corrected_axis:.9g}",
            "oec_to_paper_scale_error_percent": f"{100*(a_oec/a_one_solar_radius-1):.6g}",
            "paper_flux_reconstructed_earth": f"{reconstructed_flux:.8g}",
            "using_oec_a_flux_earth": f"{wrong_flux:.8g}",
            "paper_flux_error_percent": f"{100*(reconstructed_flux/published_irradiation-1):.6g}",
            "apparent_mass_ratio_q": f"{implied_mass_solar(p,a_oec)/mass:.8g}",
        })
    return result



def source_paper_consistency(stream):
    """Disentangle source-publication conflicts from OEC copying mistakes.

    This independent check evaluates K2-285 reported transit a/Rstar vs its
    quoted semimajor axis, and separately checks the shared stellar mass.
    """
    reader = csv.DictReader(stream)
    required = {
        "planet", "period_d", "paper_a_au", "paper_a_over_rstar",
        "paper_host_mass_solar", "paper_host_radius_solar",
    }
    if required - set(reader.fieldnames or ()):
        raise ValueError("missing source-paper check columns")
    out = []
    for row in reader:
        def field(key):
            val = positive_float(row[key])
            if val is None:
                raise ValueError(f"{row['planet']}: invalid {key}")
            return val
        a = field("paper_a_au")
        p = field("period_d")
        mass = field("paper_host_mass_solar")
        radius = field("paper_host_radius_solar")
        ar = field("paper_a_over_rstar")
        a_ratio = ar * radius * SOLAR_RADIUS_IN_AU
        m_dyn = implied_mass_solar(p, a)
        out.append({
            "planet": row["planet"],
            "period_d": f"{p:.10g}",
            "paper_a_au": f"{a:.9g}",
            "paper_a_over_rstar": f"{ar:.7g}",
            "paper_host_mass_solar": f"{mass:.7g}",
            "paper_host_radius_solar": f"{radius:.7g}",
            "mass_implied_by_p_and_a_solar": f"{m_dyn:.9g}",
            "mass_ratio_to_paper_host": f"{m_dyn/mass:.9g}",
            "axis_predicted_by_p_and_host_mass_au": f"{predicted_axis_au(p,mass):.9g}",
            "axis_from_paper_a_over_rstar_au": f"{a_ratio:.9g}",
            "paper_axis_to_ratio_axis_percent": f"{100*(a/a_ratio-1):.7g}",
        })
    return out


def write_csv(path: Path, rows: list[dict], fieldnames):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path,
                        help="local original 2023 CSV snapshot; no auto-download")
    parser.add_argument("--out", required=True, type=Path, help="new output directory")
    parser.add_argument("--threshold", type=float, default=1.5)
    parser.add_argument("--primary-checks", type=Path, default=None,
                        help="curated Kepler-445/446 primary publication and 2026 XML values")
    parser.add_argument("--source-model-checks", type=Path, default=None,
                        help="four-planet K2-285 original-paper model consistency")
    parser.add_argument("--expect-sha256", default=None,
                        help="optional fail-closed SHA-256 of original CSV bytes")
    args = parser.parse_args(argv)
    blob = args.input.read_bytes()
    checksum = hashlib.sha256(blob).hexdigest()
    if args.expect_sha256 and checksum.lower() != args.expect_sha256.lower():
        parser.error(f"SHA256 mismatch for {args.input}: got {checksum}")
    with args.input.open(encoding="utf-8-sig", newline="") as f:
        flags, counts = audit_stream(f, args.threshold)
    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(args.out / "flags.csv", flags, RESULT_FIELDS)
    summary = {
        "description": "screening only; flagged rows are NOT confirmed catalog errors",
        "units": {"period": "day", "a": "au", "mass": "solar"},
        "source_url": SOURCE_URL,
        "source_github_commit": OEC_TABLE_COMMIT,
        "source_github_blob_sha1": OEC_TABLE_BLOB,
        "source_sha256": checksum,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "threshold_factor_strict_gt": args.threshold,
        "exclusions": "confirmed, binaryflag==0, positive finite P,a,Mstar",
        "counts": counts,
    }
    if args.primary_checks:
        with args.primary_checks.open(encoding="utf-8", newline="") as f:
            primary = primary_paper_diagnostics(f)
        write_csv(args.out / "unit_conversion_checks.csv", primary, PRIMARY_FIELDS)
        summary["independently_checked_primary_paper_rows"] = len(primary)
    if args.source_model_checks:
        with args.source_model_checks.open(encoding="utf-8", newline="") as f:
            source_diagnostics = source_paper_consistency(f)
        write_csv(args.out / "source_model_checks.csv", source_diagnostics, SOURCE_FIELDS)
        summary["source_paper_model_rows_checked"] = len(source_diagnostics)
    (args.out / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
