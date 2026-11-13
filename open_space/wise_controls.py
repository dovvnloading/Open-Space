"""Inject known angular drifts into actual field detections and replay blind association."""

import argparse
from dataclasses import replace
import json
import math
from pathlib import Path
import random
from statistics import median

from .wise import DAYS_YEAR, REF_MJD, assess, fit, link_visits, prepare, seasonal_epochs, sky, tangent
from .wise_archive import acquire


def run(fields, cache, seed=30417):
    rng = random.Random(seed)
    records = []
    speeds = (0, 100, 250, 500, 1000)
    for field in fields:
        if field["role"] != "exploration":
            continue
        rows, meta = acquire(field, cache, "neowise")
        points, _ = prepare(rows, field)
        visits, _ = seasonal_epochs(points)
        tracks, _ = link_visits(visits)
        # Source positions and slopes remain blinded to the detector on each
        # trial. Conditioning on a well-sampled baseline limits applicability.
        usable = [t for t in tracks if len(t) >= 14 and t[-1].mjd - t[0].mjd >= 6 * DAYS_YEAR]
        if not usable:
            records.append({"field": field["id"], "result": "no_qualifying_control_track"})
            continue
        track = min(
            usable,
            key=lambda t: (
                -len(t),
                -sum(len(v.points) for v in t),
                min(p.detection for v in t for p in v.points),
            ),
        )
        all_times = [v.mjd for v in track]
        baseline = fit(track)
        x0, y0 = (v[0] for v in baseline)
        cos_factor = math.cos(math.radians(sky(x0, y0, field)[1])) / math.cos(
            math.radians(field["dec_deg"])
        )
        # Include isolated/excluded visit samples in the same 2.0 arcsec tube.
        # All original quality-passing field neighbors remain in the blind search.
        membership = set()
        for p in points:
            t = (p.mjd - REF_MJD) / DAYS_YEAR
            if math.hypot(p.x - x0 - baseline[0][1] * t, p.y - y0 - baseline[1][1] * t) <= 2.0:
                membership.add(p.detection)
        tube_points = [p for p in points if p.detection in membership]
        original_visits = []
        for p in sorted(tube_points, key=lambda p: p.mjd):
            if not original_visits or p.mjd - original_visits[-1][-1].mjd > 20:
                original_visits.append([])
            original_visits[-1].append(p)
        visit_residuals = []
        for group in original_visits:
            dt = (median(p.mjd for p in group) - REF_MJD) / DAYS_YEAR
            visit_residuals.append(
                (
                    median(p.x for p in group) - x0 - baseline[0][1] * dt,
                    median(p.y for p in group) - y0 - baseline[1][1] * dt,
                )
            )
        for speed in speeds:
            # Resample the measured visit offsets to avoid cancelling the
            # original best-fit slope algebraically when injecting another.
            visit_offsets = {}
            for group, original_residual in zip(original_visits, visit_residuals):
                draw = rng.choice(visit_residuals)
                for p in group:
                    visit_offsets[p.detection] = (
                        draw[0] - original_residual[0],
                        draw[1] - original_residual[1],
                    )
            angle = rng.random() * 2 * math.pi
            px = speed * math.cos(angle)
            py = speed * math.sin(angle)
            perturbed = []
            for p in points:
                if p.detection in membership:
                    t = (p.mjd - REF_MJD) / DAYS_YEAR
                    p = replace(
                        p,
                        x=p.x
                        + (px / 1000 / cos_factor - baseline[0][1]) * t
                        + visit_offsets[p.detection][0],
                        y=p.y + (py / 1000 - baseline[1][1]) * t + visit_offsets[p.detection][1],
                    )
                # Enforce the query boundary after shifting; missing flux outside
                # the original table cannot be simulated.
                if math.hypot(p.x, p.y) <= field["radius_arcmin"] * 60:
                    perturbed.append(p)
            visits, _ = seasonal_epochs(perturbed)
            estimates = []
            for inferred in link_visits(visits)[0]:
                s = assess(inferred, field, [])
                if s is None:
                    continue
                x, y = tangent(s["ra_deg"], s["dec_deg"], field)
                if math.hypot(x - x0, y - y0) <= 1.0:
                    estimates.append(s)
            matched = len(estimates) == 1
            passed = matched and estimates[0]["disposition"] != "below_motion_threshold"
            error = (
                math.hypot(estimates[0]["pmra_mas_yr"] - px, estimates[0]["pmdec_mas_yr"] - py)
                if matched
                else None
            )
            recovered = passed and error <= max(80, speed * 0.25)
            records.append(
                {
                    "field": field["id"],
                    "source_table_sha256": meta["sha256"],
                    "original_anchor": min(p.detection for v in track for p in v.points),
                    "baseline_epochs": len(all_times),
                    "modified_input_rows": len(membership),
                    "injected_speed_mas_yr": speed,
                    "injected_pmra_mas_yr": round(px, 4),
                    "injected_pmdec_mas_yr": round(py, 4),
                    "matched": matched,
                    "passed_motion_gate": passed,
                    "recovered_vector": recovered,
                    "vector_error_mas_yr": round(error, 4) if error is not None else None,
                }
            )
    summaries = []
    for speed in speeds:
        tests = [x for x in records if x.get("injected_speed_mas_yr") == speed]
        errors = [x["vector_error_mas_yr"] for x in tests if x["matched"]]
        summaries.append(
            {
                "injected_speed_mas_yr": speed,
                "trials": len(tests),
                "tracks_reassociated": len(errors),
                "passed_motion_gate": sum(x["passed_motion_gate"] for x in tests),
                "recovered_vectors": sum(x["recovered_vector"] for x in tests),
                "median_vector_error_mas_yr": median(errors) if errors else None,
            }
        )
    return {
        "experiment": "real_field_astrometric_drift_injection_v1",
        "seed": seed,
        "summary": summaries,
        "trials": records,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=Path("data/wise-fields.json"))
    parser.add_argument("--cache", required=True, type=Path)
    a = parser.parse_args()
    print(
        json.dumps(run(json.loads(a.plan.read_bytes())["fields"], a.cache), indent=2, allow_nan=False)
    )
