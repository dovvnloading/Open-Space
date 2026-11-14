"""Score authenticated post-DR3 Gaia single-lined epoch velocities against the frozen orbit."""

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from orbit import EPOCH, NOMINAL, predict_jd, sample_orbits, velocity

# End of the DR3 observation interval in native Gaia six-hour revolutions.
DR3_END_OBMT = 5230.0880
TIME_ZERO = 2455197.5
DAY_NS = 86400e9
PRODUCT = "EPOCH_PARAMETERS_RVS_SINGLE"
REQUIRED = {
    "source_id",
    "transit_id",
    "rv_validity_flag",
    "radial_velocity",
    "radial_velocity_error",
    "obs_time_rv",
    "radial_velocity_status",
}


def transit_revolutions(transit_id):
    return (((transit_id >> 17) & 0x000003FFFFFFFFFF) * 204800) / 21600e9


def read_rows(path, source_id, clock):
    seen = set()
    rows = []
    counts = Counter()
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        needed = REQUIRED | ({"obs_time_bary_corr"} if clock == "spacecraft" else set())
        missing = needed - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                "Missing documented DR4 single-lined CSV fields: "
                + ", ".join(sorted(missing))
            )
        for record in reader:
            counts["archive_rows"] += 1
            if record["source_id"] != source_id:
                raise ValueError("Unexpected release source ID")
            transit = int(record["transit_id"])
            if transit in seen or transit <= 0:
                raise ValueError("Duplicate or invalid Gaia transit_id")
            seen.add(transit)
            rev = transit_revolutions(transit)
            if rev <= DR3_END_OBMT:
                counts["overlapping_dr3"] += 1
                continue
            valid = record["rv_validity_flag"].strip().lower()
            if valid not in {"true", "false", "1", "0", "t", "f"}:
                raise ValueError(f"Indeterminate validity flag for transit {transit}")
            if valid in {"false", "0", "f"}:
                counts["invalidated_by_archive"] += 1
                continue
            fields = ["obs_time_rv", "radial_velocity", "radial_velocity_error"]
            if clock == "spacecraft":
                fields.append("obs_time_bary_corr")
            try:
                numbers = [float(record[k]) for k in fields]
            except (ValueError, TypeError):
                numbers = []
            if not numbers or not np.isfinite(numbers).all() or numbers[2] <= 0:
                counts["missing_or_invalid_measurement"] += 1
                continue
            day, rv, error = numbers[:3]
            jd = TIME_ZERO + day
            # The transit ID encodes OBMT, independently of the RV product clock.
            # This loose (1-day) sanity check catches a JD/MJD/epoch origin mismatch.
            if abs(jd - (2457023.75 + (rev - 1717.6256) / 4)) > 1:
                raise ValueError(
                    f"Observing epoch disagrees with transit ID for {transit}"
                )
            if clock == "spacecraft":
                if abs(numbers[3]) > 600e9:
                    raise ValueError(
                        "Light-time correction outside spacecraft range (expected ns)"
                    )
                jd += numbers[3] / DAY_NS
            status = record["radial_velocity_status"].strip()
            flagged = bool(int(status) & ((1 << 55) | (1 << 56))) if status else None
            rows.append((transit, jd, rv, error, flagged))
    counts["scored"] = len(rows)
    return sorted(rows, key=lambda row: row[1]), counts


def compare(rows):
    if not rows:
        return {"n": 0, "status": "no eligible velocities"}
    jd, observed, error = np.asarray([row[1:4] for row in rows]).T
    prediction = predict_jd(jd)
    weights = 1 / error**2
    constant = np.average(observed, weights=weights)
    residual = observed - prediction
    chisq_fixed = float(np.sum(weights * residual**2))
    chisq_flat = float(np.sum(weights * (observed - constant) ** 2))
    angles = np.sort(((jd - EPOCH - NOMINAL[5]) / NOMINAL[4]) % 1)
    # Correlated orbit draws evaluate the same orbital realization at all epochs.
    draws = velocity(jd[:, None], sample_orbits())
    model_cov = np.atleast_2d(np.cov(draws))
    total_cov = np.diag(error**2) + model_cov
    whitened = np.linalg.solve(
        np.linalg.cholesky(total_cov), observed - np.mean(draws, axis=1)
    )
    result = {
        "n": len(rows),
        "span_days": float(np.ptp(jd)),
        "largest_phase_gap_cycles": float(
            np.max(np.diff(angles, append=angles[0] + 1))
        ),
        "predicted_range_km_s": float(np.ptp(prediction)),
        "epoch_velocity_error_median_km_s": float(np.median(error)),
        "weighted_rms_fixed_km_s": float(np.sqrt(chisq_fixed / weights.sum())),
        "weighted_rms_fitted_constant_km_s": float(np.sqrt(chisq_flat / weights.sum())),
        "fitted_constant_km_s": float(constant),
        "chisq_fixed_measurement_only": chisq_fixed,
        "chisq_fitted_constant_measurement_only": chisq_flat,
        "chisq_flat_minus_fixed": chisq_flat - chisq_fixed,
        "quadratic_form_with_joint_orbit_and_measurement_covariance": float(
            whitened @ whitened
        ),
    }
    x = prediction - np.average(prediction, weights=weights)
    leverage = np.sum(weights * x**2)
    if len(rows) >= 2 and np.ptp(prediction) > 1e-10 and leverage > 0:
        result["diagnostic_free_slope"] = float(
            np.sum(weights * x * (observed - constant)) / leverage
        )
        result["free_slope_sigma_measurement_only"] = float(1 / np.sqrt(leverage))
    return result


def run(receipt_path, clock):
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("dr3_id") != "5946882743944130304":
        raise ValueError("Wrong DR3 target in retrieval receipt")
    product = receipt.get("retrieved", {}).get(PRODUCT)
    if not product:
        return {
            "status": "no single-lined product retrieved",
            "access_status": receipt["status"],
            "n": 0,
        }
    path = receipt_path.parent / product["file"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != product["sha256"]:
        raise ValueError("Product differs from retrieval receipt")
    rows, counts = read_rows(path, receipt["dr4_id"], clock)
    if counts["archive_rows"] != product["rows"]:
        raise ValueError("Source row count differs from receipt")
    return {
        "dr4_id": receipt["dr4_id"],
        "clock": clock,
        "counts": dict(counts),
        "all_valid": compare(rows),
        "known_unblended": compare([r for r in rows if r[4] is False]),
        "epochs": [
            {
                "transit_id": str(t),
                "bjd_tcb": jd,
                "rv_km_s": rv,
                "error_km_s": err,
                "predicted_tracer_km_s": float(pred),
                "deblended_or_reblended": flag,
            }
            for (t, jd, rv, err, flag), pred in zip(
                rows, predict_jd(np.array([r[1] for r in rows]))
            )
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    parser.add_argument(
        "--clock",
        required=True,
        choices=["barycentric", "spacecraft"],
        help="frame of obs_time_rv according to final release metadata",
    )
    args = parser.parse_args()
    print(json.dumps(run(args.receipt, args.clock), indent=2, allow_nan=False))
