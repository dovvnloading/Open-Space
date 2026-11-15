"""Test the 2024 Venus-flyby radio contrast within a single spacecraft attitude."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import cdflib
import numpy as np

from analyze import ANTENNAS, BANDS, contrast, load_minutes, rounded

HERE = Path(__file__).resolve().parent
DATE = "2024-11-06"
URL = ("https://spdf.gsfc.nasa.gov/pub/data/psp/sweap/spi/l3/spi_sf00_l3_mom/"
       "2024/psp_swp_spi_sf00_l3_mom_20241106_v04.cdf")
SHA = "9b435d624e0f5cc84af997e966bbd6090aef7e8f7cd5b1ba7fd933b21926fc8e"


def attitude_summary(cdf_path, center):
    cdf = cdflib.CDF(cdf_path)
    times = cdflib.cdfepoch.to_datetime(cdf.varget("Epoch"))
    q = cdf.varget("QUAT_SC_TO_RTN").astype(float)
    valid = np.isfinite(q).all(axis=1) & (np.abs(np.linalg.norm(q, axis=1) - 1) < 1e-4)
    q, times = q[valid], times[valid]
    q /= np.linalg.norm(q, axis=1)[:, None]
    minutes = (times - np.datetime64(DATE)) / np.timedelta64(1, "m")
    offset = minutes - center
    q0 = q[np.argmin(abs(offset))]
    # The absolute quaternion dot product gives relative rotation independent of sign,
    # scalar ordering, and active/passive convention.
    angle = np.rad2deg(2 * np.arccos(np.clip(abs(q @ q0), 0, 1)))
    windows = {}
    for half_width in (3, 10, 50):
        use = abs(offset) <= half_width
        dot = np.min(np.abs(q[use] @ q[use].T))
        windows[str(half_width)] = {
            "records": int(use.sum()),
            "max_pairwise_angle_deg": round(float(np.rad2deg(2*np.arccos(np.clip(dot, 0, 1)))), 6),
            "max_sampling_gap_s": round(float(np.max(np.diff(minutes[use]))*60), 3),
        }
    fixed = (abs(offset) <= 30) & (angle <= 0.03)
    fringe = (abs(offset) >= 20) & (abs(offset) < 50)
    return {
        "cdf_url": URL, "cdf_sha256": SHA, "frame": "spacecraft to RTN",
        "half_window_minutes": windows, "central_attitude_threshold_deg": 0.03,
        "first_and_last_central_attitude_within_half_hour_utc": [
            str(times[fixed][0]), str(times[fixed][-1])],
        "max_departure_across_this_entire_interval_deg":
            round(float(angle[(times>=times[fixed][0]) & (times<=times[fixed][-1])].max()), 6),
        "median_attitude_difference_from_center_in_20_to_50_min_flanks_deg":
            round(float(np.median(angle[fringe])), 6),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cdf", type=Path, required=True,
                        help="cache path for the verified SPAN-I daily CDF")
    parser.add_argument("--minutes", type=Path, default=HERE / "minutes.csv.gz")
    args = parser.parse_args()
    if not args.cdf.exists():
        with urlopen(URL, timeout=120) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != SHA:
            raise ValueError("Upstream CDF differs from the recorded version")
        args.cdf.write_bytes(data)
    if hashlib.sha256(args.cdf.read_bytes()).hexdigest() != SHA:
        raise ValueError("Cached CDF does not match its recorded digest")
    manifest = json.loads((HERE / "manifest.json").read_text())
    center = next(e["maximum_radius_utc_minute"] for e in manifest if e["date"] == DATE)
    res = attitude_summary(args.cdf, center)
    res.update({"date": DATE, "center_utc_minute": center, "bands": BANDS,
                "receiver_nsum": 40, "core_half_width_min": 3,
                "two_flanks_abs_offset_min": [4, 10]})
    series = load_minutes(args.minutes)
    radio = {}
    for antenna in ANTENNAS:
        effects = contrast(series[DATE, antenna, 40], center, 3, 4, 10)
        radio[antenna] = dict(zip(("vs_average_flanks_db", "vs_before_db", "vs_after_db"),
                                  map(rounded, effects)))
    res["radio"] = radio
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
