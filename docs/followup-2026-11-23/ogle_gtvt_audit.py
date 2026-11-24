"""Audit the captured GTVT 1.2 result against its dated Horizons vectors.

No network request is made and no installed package is patched. The unmodified
tool CSV is preserved. A separate object receives vectors matched by JDTDB,
then runs GTVT's existing fixed-target calculations. Requires jwst_gtvt 1.2
and its dependencies; run this script from any directory.
"""
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time
import jwst_gtvt.jwst_tvt as gtvt

HERE = Path(__file__).resolve().parent
OUT = HERE / "ogle-evidence"
EPHEMERIS = OUT / "gtvt-horizons-ephemeris.txt"
ORIGINAL = OUT / "gtvt-original.csv"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def windows(frame, mask):
    indices = np.flatnonzero(mask)
    if not len(indices):
        return []
    pieces = np.split(indices, np.flatnonzero(np.diff(indices) > 1) + 1)
    result = []
    for piece in pieces:
        first, last = frame.iloc[piece[0]], frame.iloc[piece[-1]]
        result.append({
            "first_daily_sample_tdb_date": str(first["display_date"]),
            "last_daily_sample_tdb_date": str(last["display_date"]),
            "daily_samples": len(piece),
            "first_jd_tdb": float(first["JDTDB"]),
            "last_jd_tdb": float(last["JDTDB"]),
            "first_sample_utc": Time(first["JDTDB"], format="jd", scale="tdb").utc.isot,
            "last_sample_utc": Time(last["JDTDB"], format="jd", scale="tdb").utc.isot,
        })
    return result


def main():
    original = pd.read_csv(ORIGINAL)
    lines = EPHEMERIS.read_text().splitlines()
    # __init__ calls Horizons and shifts vectors. Bypass only initialization;
    # the existing parser and downstream geometry/angle methods remain intact.
    obj = object.__new__(gtvt.Ephemeris)
    obj.ephemeris = np.array(lines)
    raw = obj.convert_ephemeris_to_df(obj.ephemeris)
    assert len(raw) == len(original) == 365
    assert raw["JDTDB"].is_unique and original["JDTDB"].is_unique
    assert np.all(np.diff(raw["JDTDB"]) == 1)
    matched = raw.set_index("JDTDB").loc[original["JDTDB"].to_numpy()]
    xyz = ["X", "Y", "Z"]
    assert original[xyz].iloc[0].isna().all()
    lag_error = np.max(np.abs(original[xyz].to_numpy()[1:] - matched[xyz].to_numpy()[:-1]))
    assert lag_error < 1e-6

    obj.dataframe = original.copy(deep=True)
    obj.dataframe[xyz] = matched[xyz].to_numpy()
    assert np.isfinite(obj.dataframe[xyz].to_numpy()).all()
    ra = float(original["ra"].iloc[0])
    dec = float(original["dec"].iloc[0])
    aligned = obj.get_fixed_target_positions(str(ra), str(dec))
    assert np.array_equal(aligned["JDTDB"].to_numpy(), original["JDTDB"].to_numpy())
    assert np.allclose(aligned[xyz].to_numpy(), matched[xyz].to_numpy(), rtol=0, atol=0)

    # Independent vector-dot-product check of the recomputed separation.
    r, d = np.deg2rad([ra, dec])
    target = np.array([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])
    sun = -matched[xyz].to_numpy()
    sun /= np.linalg.norm(sun, axis=1)[:, None]
    direct = np.arccos(np.clip(sun @ target, -1, 1))
    angle_error = float(np.max(np.abs(direct - aligned["dist"].to_numpy())))
    assert angle_error < 1e-10
    direct_mask = (direct > gtvt.MIN_SUN_ANGLE) & (direct < gtvt.MAX_SUN_ANGLE)
    assert np.array_equal(direct_mask, aligned["in_FOR"].to_numpy())
    july = aligned["display_date"].astype(str).str.startswith("2027-07-").to_numpy()
    assert july.sum() == 31 and not direct_mask[july].any()

    aligned_path = OUT / "gtvt-date-aligned.csv"
    aligned.to_csv(aligned_path, index=False)
    header = lines[:lines.index("$$SOE")]
    metadata_keys = ("Target body", "Center body", "Start time", "Stop  time", "Step-size",
                     "Output units", "Output type", "Reference frame")
    report = {
        "audit_date": "2026-11-23",
        "position_deg_icrf": {"ra": ra, "dec": dec},
        "software": {p: importlib.metadata.version(p) for p in
                     ("jwst_gtvt", "numpy", "pandas", "astropy", "pysiaf")},
        "module": "jwst_gtvt/jwst_tvt.py",
        "module_sha256": sha256(Path(gtvt.__file__)),
        "script_sha256": sha256(Path(__file__)),
        "input_sha256": {EPHEMERIS.name: sha256(EPHEMERIS), ORIGINAL.name: sha256(ORIGINAL)},
        "aligned_csv_sha256": sha256(aligned_path),
        "horizons_header_metadata": [line for line in header if line.startswith(metadata_keys)],
        "horizons_request_template": gtvt.URL,
        "time_frame": "JDTDB and Calendar Date (TDB); daily midnight TDB samples. Display dates must not be labeled UTC.",
        "coverage": {"raw_vectors": len(raw), "original_csv_rows": len(original),
                     "original_missing_position_rows": int(original[xyz].isna().any(axis=1).sum()),
                     "aligned_missing_position_rows": int(aligned[xyz].isna().any(axis=1).sum())},
        "finding": "Installed GTVT 1.2 __init__ replaces each coordinate by (coordinate.shift()-coordinate)*1+coordinate. This equals the previous row while retaining the current JDTDB/date.",
        "previous_row_position_matches": {"rows": len(original) - 1,
                                          "maximum_absolute_difference_km": float(lag_error)},
        "aligned_method": "Deep copy of original table; replace only X/Y/Z from raw Horizons rows joined exactly on JDTDB; rerun existing get_fixed_target_positions/build_dataframe. No package file changed.",
        "field_of_regard_criterion": {"minimum_degrees_strict": float(np.rad2deg(gtvt.MIN_SUN_ANGLE)),
                                      "maximum_degrees_strict": float(np.rad2deg(gtvt.MAX_SUN_ANGLE))},
        "unmodified_gtvt_windows": windows(original, original["in_FOR"].to_numpy()),
        "date_aligned_gtvt_windows": windows(aligned, direct_mask),
        "date_aligned_windows_with_85_degree_lower_bound": windows(aligned, (direct >= np.deg2rad(85)) & (direct <= np.deg2rad(135))),
        "independent_vector_check_maximum_angle_difference_radians": angle_error,
        "july_2027": {"daily_samples": 31, "samples_in_field_of_regard": int(direct_mask[july].sum()),
                      "minimum_elongation_degrees": float(np.rad2deg(direct[july]).min()),
                      "maximum_elongation_degrees": float(np.rad2deg(direct[july]).max())},
        "limits": [
            "First/last accepted daily samples are not exact visibility boundary instants or full-day scheduling guarantees.",
            "Horizons states are the supplied LT-corrected JWST_merged vectors centered on the Sun (DE441), in ICRF and KM-S; this audit does not replace or improve that ephemeris.",
            "The date-aligned output is an explicitly modified calculation using the tool's methods, not an unmodified official GTVT result.",
            "Fixed input target coordinates have no propagated stellar proper motion; this is a feasibility estimate, not a pointing or astrometric reduction solution.",
            "Instrument angle columns are tool outputs only; no target-specific roll/detector/guide-star optimization, ETC, APT, allocation or scheduling check was completed.",
            "Daily TDB labels convert to approximately 69 seconds earlier in UTC under the supplied time-scale tables; retain explicit time-scale labels.",
        ],
    }
    (OUT / "gtvt-audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({k: report[k] for k in
                      ("software", "previous_row_position_matches", "unmodified_gtvt_windows",
                       "date_aligned_gtvt_windows", "july_2027")}, indent=2))


if __name__ == "__main__":
    main()
