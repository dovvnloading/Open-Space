"""Reconstruct the spherical Enceladus impact height for the epsilon CMa ray."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
import spiceypy as spice

from check_clocks import HERE, ZERO
from analyze import sha256


def replay(root):
    receipt = json.loads((HERE / "geometry-sources.json").read_text(encoding="utf-8"))
    spice.kclear()
    try:
        for entry in receipt["files"]:
            path = root / entry["file"]
            if sha256(path) != entry["sha256"]:
                raise ValueError(f"SPICE kernel mismatch: {path}")
            spice.furnsh(str(path))
        star = receipt["star"]
        delta_year = star["propagated_epoch"] - 2000
        ra, dec = star["icrs_j2000_ra_deg"], star["icrs_j2000_dec_deg"]
        ra += delta_year * star["pmra_cosdec_mas_year"] / (3600000 * math.cos(math.radians(dec)))
        dec += delta_year * star["pmdec_mas_year"] / 3600000
        ray = spice.radrec(1, *np.deg2rad((ra, dec)))
        origin = spice.str2et(ZERO)

        def estimate(et):
            position, light_time = spice.spkpos("ENCELADUS", et, "J2000", "CN", "CASSINI")
            along = np.dot(position, ray)
            if along <= 0:
                raise ValueError("Target is behind the stellar ray")
            off = position - along * ray
            return float(np.linalg.norm(off) - 252.1), float(light_time)

        rows = []
        for second in range(540, 902):
            alt, light_time = estimate(origin + second)
            shifted, _ = estimate(origin + second + 1)
            rows.append({
                "seconds_from_14_30_utc": second,
                "utc": spice.et2utc(origin + second, "ISOC", 0),
                "spherical_altitude_km": round(alt, 6), "target_light_time_s": round(light_time, 9),
                "altitude_difference_one_second_later_km": round(shifted-alt, 6),
            })
        inside = [row for row in rows if 550 <= row["seconds_from_14_30_utc"] < 645]
        max_difference = max(abs(row["altitude_difference_one_second_later_km"]) for row in inside)
        return {"spherical_radius_km": 252.1, "propagated_star_ra_dec_deg": [ra, dec],
                "minimum_sample": min(rows, key=lambda row: row["spherical_altitude_km"]),
                "encounter_max_one_second_altitude_change_km": max_difference, "samples": rows}
    finally:
        spice.kclear()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kernel_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = replay(args.kernel_root)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(result["minimum_sample"], result["encounter_max_one_second_altitude_change_km"])


if __name__ == "__main__":
    main()
