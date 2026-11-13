"""W1 single-exposure astrometry across NEOWISE mission visits."""

import argparse
import hashlib
from collections import Counter
from dataclasses import dataclass
import json
import math
from pathlib import Path
from statistics import median
import sys

from .wise_archive import acquire, field_params

REF_MJD = 57388.5  # J2016.0, Gaia DR3's reference epoch
DAYS_YEAR = 365.25
SOLAR_EPOCH = 51544.5


def number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def tangent(ra, dec, field):
    dra = (ra - field["ra_deg"] + 180.0) % 360.0 - 180.0
    return (3600 * dra * math.cos(math.radians(field["dec_deg"])), 3600 * (dec - field["dec_deg"]))


def sky(x, y, field):
    return (
        (field["ra_deg"] + x / (3600 * math.cos(math.radians(field["dec_deg"])))) % 360,
        field["dec_deg"] + y / 3600,
    )


@dataclass
class Point:
    mjd: float
    x: float
    y: float
    mag: float
    snr: float
    frame: str
    detection: str
    allwise: str


@dataclass
class Epoch:
    mjd: float
    x: float
    y: float
    mag: float
    points: list


def quality_reason(row):
    n = {
        k: number(row[k])
        for k in (
            "w1mpro",
            "w1sigmpro",
            "w1snr",
            "w1rchi2",
            "w1sat",
            "qual_frame",
            "qi_fact",
            "saa_sep",
            "sso_flg",
            "nb",
            "na",
            "sigra",
            "sigdec",
            "ra",
            "dec",
            "mjd",
        )
    }
    if any(v is None for v in n.values()):
        return "missing_measurement"
    if not (n["w1snr"] >= 8 and n["w1mpro"] < 15 and n["w1sigmpro"] > 0):
        return "faint_or_uncertain"
    if not (
        n["qual_frame"] > 0
        and n["qi_fact"] > 0
        and n["saa_sep"] > 0
        and row["moon_masked"][:1] == "0"
        and row["cc_flags"][:1] == "0"
    ):
        return "exposure_or_artifact_flag"
    if (
        n["w1sat"] != 0
        or not 0 < n["w1rchi2"] < 10
        or n["nb"] != 1
        or n["na"] != 0
        or n["sso_flg"] != 0
        or min(n["sigra"], n["sigdec"]) <= 0
    ):
        return "saturation_blend_or_fit"
    return None


def _epoch_clusters(points):
    # Highest-S/N positions seed within-visit associations. Exclude ambiguous
    # assignments and keep just the best measurement per scan/frame and source.
    groups = []
    ambiguous = 0
    for p in sorted(points, key=lambda p: (-p.snr, p.detection)):
        close = [g for g in groups if math.hypot(p.x - g[0].x, p.y - g[0].y) <= 1.6]
        if not close:
            groups.append([p])
        elif len(close) == 1:
            close[0].append(p)
        else:
            ambiguous += 1
    out = []
    for group in groups:
        unique = {}
        for p in group:
            unique.setdefault(p.frame, p)
        ps = list(unique.values())
        if len(ps) >= 4:
            out.append(
                Epoch(*(median(getattr(p, k) for p in ps) for k in ("mjd", "x", "y", "mag")), ps)
            )
    return sorted(out, key=lambda e: (e.x, e.y)), ambiguous


def seasonal_epochs(points):
    visits = []
    for p in sorted(points, key=lambda p: (p.mjd, p.detection)):
        if not visits or p.mjd - visits[-1][-1].mjd > 20:
            visits.append([])
        visits[-1].append(p)
    result, confused = [], 0
    for visit in visits:
        entries, ambiguous = _epoch_clusters(visit)
        result.append(entries)
        confused += ambiguous
    return result, confused


def _predict(track, date):
    last = track[-1]
    if len(track) < 2:
        return last.x, last.y
    first = track[max(0, len(track) - 6)]
    dt = (date - last.mjd) / (last.mjd - first.mjd)
    return last.x + (last.x - first.x) * dt, last.y + (last.y - first.y) * dt


def link_visits(visits):
    tracks = []
    ambiguous = 0
    for visit in visits:
        pairs = []
        for ti, track in enumerate(tracks):
            for vi, v in enumerate(visit):
                dt = v.mjd - track[-1].mjd
                if not 0 < dt <= 760:
                    continue
                x, y = _predict(track, v.mjd)
                sep = math.hypot(v.x - x, v.y - y)
                if sep <= min(5.0, 1.1 + 2.4 * dt / DAYS_YEAR):
                    pairs.append((sep, ti, vi))
        # If either endpoint can plausibly join two objects, discard both
        # ambiguous comparisons; a faint blend should not push a track forward.
        tc, vc = Counter(ti for _, ti, _ in pairs), Counter(vi for _, _, vi in pairs)
        used = set()
        for _, ti, vi in pairs:
            if tc[ti] == vc[vi] == 1:
                tracks[ti].append(visit[vi])
                used.add(vi)
            else:
                ambiguous += 1
        for vi, v in enumerate(visit):
            if vi not in used and vc[vi] == 0:
                tracks.append([v])
    return tracks, ambiguous


def fit(track):
    times = [(v.mjd - REF_MJD) / DAYS_YEAR for v in track]
    tx = sum(times) / len(times)
    xx = sum((t - tx) ** 2 for t in times)
    if not xx or len(times) < 3:
        raise ValueError("At least three separated epochs needed")
    fits = []
    for axis in ("x", "y"):
        values = [getattr(v, axis) for v in track]
        avg = sum(values) / len(values)
        b = sum((t - tx) * (y - avg) for t, y in zip(times, values)) / xx
        a = avg - b * tx
        res = [y - a - b * t for t, y in zip(times, values)]
        # Floor the epoch scatter at 100 mas; errors are approximate and omit
        # calibration covariance and parallax. Each visit has equal weight.
        sigma = max(0.1, math.sqrt(sum(y * y for y in res) / (len(res) - 2)))
        fits.append((a, b, sigma / math.sqrt(xx)))
    return fits


def _separation(track, fitted):
    return [
        math.hypot(
            *[
                getattr(v, axis) - a - b * (v.mjd - REF_MJD) / DAYS_YEAR
                for axis, (a, b, _) in zip(("x", "y"), fitted)
            ]
        )
        for v in track
    ]


def _scan_phase(mjd):
    # Half-year observing-geometry proxy. Explicitly not a flight attitude.
    return ((mjd - SOLAR_EPOCH) / DAYS_YEAR) % 1


def _geometry_subsets(track):
    # Split by annual phase relative to the first observed season.
    reference = _scan_phase(track[0].mjd)
    return [
        [v for v in track if (math.cos(2 * math.pi * (_scan_phase(v.mjd) - reference)) > 0) == side]
        for side in (True, False)
    ]


def gaia_match(x, y, field, rows):
    matches = []
    for row in rows:
        ra, dec, ref = (number(row[k]) for k in ("ra", "dec", "ref_epoch"))
        if ref is None or abs(ref - 2016.0) > 0.01 or ra is None or dec is None:
            continue
        xx, yy = tangent(ra, dec, field)
        d = math.hypot(xx - x, yy - y)
        if d <= 2.0:
            matches.append((d, row))
    matches.sort(key=lambda pair: (pair[0], pair[1]["source_id"]))
    if not matches:
        return None
    sep, row = matches[0]
    return {
        "source_id": row["source_id"],
        "separation_arcsec": round(sep, 4),
        "matches_within_2arcsec": len(matches),
        "pmra_mas_yr": number(row["pmra"]),
        "pmdec_mas_yr": number(row["pmdec"]),
        "ruwe": number(row["ruwe"]),
    }


def assess(track, field, gaia):
    span = track[-1].mjd - track[0].mjd
    if len(track) < 8 or span < 4 * DAYS_YEAR:
        return None
    chronology = max(4, int(len(track) * 0.6))
    train, test = track[:chronology], track[chronology:]
    early = fit(train)
    pred = median(_separation(test, early))
    x0, y0 = median(v.x for v in train), median(v.y for v in train)
    stationary_pred = median(math.hypot(v.x - x0, v.y - y0) for v in test)
    initial = fit(track)
    residuals = _separation(track, initial)
    med = median(residuals)
    mad = median(abs(r - med) for r in residuals)
    keep = [v for v, r in zip(track, residuals) if r <= max(0.55, med + 4.5 * 1.4826 * mad)]
    if len(keep) < 8 or keep[-1].mjd - keep[0].mjd < 4 * DAYS_YEAR:
        return None
    (x, px, ex), (y, py, ey) = best = fit(keep)
    ra, dec = sky(x, y, field)
    # x uses a fixed tangent reference; report mu_alpha*cos(delta_source).
    alpha_factor = math.cos(math.radians(dec)) / math.cos(math.radians(field["dec_deg"]))
    px, ex = px * alpha_factor, ex * alpha_factor
    speed = math.hypot(px, py)
    significance = speed / max(ex, ey)
    phase = []
    for sub in _geometry_subsets(keep):
        if len(sub) < 4 or sub[-1].mjd - sub[0].mjd < 3 * DAYS_YEAR:
            phase.append(None)
        else:
            model = fit(sub)
            phase.append([round(model[0][1] * alpha_factor * 1000, 3), round(model[1][1] * 1000, 3)])
    phase_ok = all(
        p is not None
        and math.hypot(p[0] - px * 1000, p[1] - py * 1000) <= max(80.0, 0.35 * speed * 1000)
        for p in phase
    )
    eligible = (
        speed >= 0.15
        and significance >= 6
        and pred <= 0.8
        and pred < 0.6 * stationary_pred
        and phase_ok
        and len(keep) >= 0.8 * len(track)
    )
    association = gaia_match(x, y, field, gaia)
    allwise = Counter(p.allwise for v in keep for p in v.points if p.allwise not in ("null", "", "0"))
    ids = [{"id": k, "exposures": n} for k, n in allwise.most_common(3)]
    disposition = "below_motion_threshold"
    agrees = None
    vector_difference = None
    if (
        association
        and association["matches_within_2arcsec"] == 1
        and all(association[k] is not None for k in ("pmra_mas_yr", "pmdec_mas_yr"))
    ):
        vector_difference = math.hypot(
            1000 * px - association["pmra_mas_yr"], 1000 * py - association["pmdec_mas_yr"]
        )
        agrees = vector_difference <= max(80, 4 * max(ex, ey) * 1000)
    if eligible:
        if association and association["source_id"] == field.get("control_gaia_id") and agrees:
            disposition = "known_positive_control"
        elif agrees:
            disposition = "catalog_motion_recovered"
        elif association is not None:
            disposition = "catalog_comparison_followup"
        else:
            disposition = "requires_catalog_and_image_vetting"
    anchor = min(p.detection for v in keep for p in v.points)
    return {
        "source_id": f"wise-{anchor}",
        "ra_deg": round(ra, 8),
        "dec_deg": round(dec, 8),
        "reference_mjd": REF_MJD,
        "pmra_mas_yr": round(px * 1000, 3),
        "pmdec_mas_yr": round(py * 1000, 3),
        "pmra_error_mas_yr": round(ex * 1000, 3),
        "pmdec_error_mas_yr": round(ey * 1000, 3),
        "motion_mas_yr": round(speed * 1000, 3),
        "motion_significance_proxy": round(significance, 3),
        "disposition": disposition,
        "gaia": association,
        "gaia_motion_consistent": agrees,
        "gaia_vector_difference_mas_yr": (
            round(vector_difference, 3) if vector_difference is not None else None
        ),
        "gaia_precision_check": (
            "difference_exceeds_4x_wise_axis_error"
            if vector_difference is not None and vector_difference > 4 * 1000 * max(ex, ey)
            else None
        ),
        "allwise_associations": ids,
        "quality": {
            "epochs": len(keep),
            "rejected_epochs": len(track) - len(keep),
            "first_mjd": keep[0].mjd,
            "last_mjd": keep[-1].mjd,
            "exposures": sum(len(v.points) for v in keep),
            "median_w1_mag": median(v.mag for v in keep),
            "training_epochs": len(train),
            "late_epochs": len(test),
            "forecast_error_arcsec": round(pred, 4),
            "stationary_error_arcsec": round(stationary_pred, 4),
            "annual_phase_motions_mas_yr": phase,
            "compatible_phases": phase_ok,
            "median_residual_arcsec": round(median(_separation(keep, best)), 4),
        },
        "epoch_measurements": [
            {
                "mjd": round(v.mjd, 7),
                "ra_deg": round(sky(v.x, v.y, field)[0], 8),
                "dec_deg": round(sky(v.x, v.y, field)[1], 8),
                "w1_mag": round(v.mag, 4),
                "n": len(v.points),
                "example_detection": min(p.detection for p in v.points),
            }
            for v in keep
        ],
    }


def prepare(rows, field):
    field_params(field)
    points, excluded = [], Counter()
    for r in rows:
        reason = quality_reason(r)
        if reason:
            excluded[reason] += 1
            continue
        x, y = tangent(float(r["ra"]), float(r["dec"]), field)
        points.append(
            Point(
                float(r["mjd"]),
                x,
                y,
                float(r["w1mpro"]),
                float(r["w1snr"]),
                r["scan_id"] + "-" + r["frame_num"],
                r["source_id"],
                r["allwise_cntr"],
            )
        )
    return points, excluded


def detect(rows, gaia, field):
    points, excluded = prepare(rows, field)
    visits, confused = seasonal_epochs(points)
    tracks, ambiguous = link_visits(visits)
    outputs = [candidate for t in tracks if (candidate := assess(t, field, gaia)) is not None]
    outputs.sort(key=lambda s: (-s["motion_mas_yr"], s["source_id"]))
    counts = {
        "input_exposures_rows": len(rows),
        "quality_accepted": len(points),
        "exclusions": dict(excluded),
        "seasonal_visits": len(visits),
        "visit_sources": sum(map(len, visits)),
        "tracks_started": len(tracks),
        "qualified_tracks": len(outputs),
        "ambiguous_within_visit_assignments": confused,
        "ambiguous_cross_visit_links": ambiguous,
    }
    return outputs, counts


def run_field(field, cache_dir):
    measurements, neowise = acquire(field, cache_dir, "neowise")
    gaia_rows, gaia = acquire(field, cache_dir, "gaia")
    sources, counts = detect(measurements, gaia_rows, field)
    for s in sources:
        if s["disposition"] == "below_motion_threshold":
            del s["epoch_measurements"]
    candidates = [
        s
        for s in sources
        if s["disposition"] in ("catalog_comparison_followup", "requires_catalog_and_image_vetting")
    ]
    return {
        "schema_version": 1,
        "field": field,
        "inputs": {"neowise": neowise, "gaia": gaia},
        "counts": counts,
        "sources": sources,
        "candidates": candidates,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan", type=Path, default=Path("data/wise-fields.json"))
    p.add_argument("--field", action="append", help="Process one field ID; repeat to choose several")
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--compact", action="store_true", help="Emit compact JSON for archival storage")
    a = p.parse_args(argv)
    try:
        plan = json.loads(a.plan.read_bytes())
        fields = plan["fields"]
        if len({f["id"] for f in fields}) != len(fields):
            raise ValueError("Field IDs must be unique in a plan")
        if a.field:
            chosen = set(a.field)
            fields = [f for f in fields if f["id"] in chosen]
            if len(chosen) != len(fields):
                raise ValueError("Unknown field ID")
        result = {
            "schema_version": 1,
            "plan_sha256": hashlib.sha256(a.plan.read_bytes()).hexdigest(),
            "fields": [run_field(field, a.cache) for field in fields],
        }
        if a.compact:
            print(json.dumps(result, separators=(",", ":"), sort_keys=True, allow_nan=False))
        else:
            print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        print(f"wise: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
