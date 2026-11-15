"""Screen all SPARC radial pairs. Run with Python 3 and two pinned archive files."""

import argparse
import csv
import hashlib
import itertools
import math
import sys
from pathlib import Path
import zipfile

HASHES = {
    "Rotmod_LTG.zip": "0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588",
    "SPARC_Lelli2016c.mrt": "5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243",
}


def scan_galaxy(points, distance, distance_err, inc, inc_err):
    scale_low = max(0.01, 1 - 2 * distance_err / distance)
    scale_high = 1 + 2 * distance_err / distance
    sin_inc = math.sin(math.radians(inc))
    scale_low *= (math.sin(math.radians(max(1, inc - 2 * inc_err))) / sin_inc) ** 2
    scale_high *= (math.sin(math.radians(min(90, inc + 2 * inc_err))) / sin_inc) ** 2
    nominal, robust = [], []
    gas_only = 0
    for radius, velocity, err, gas, disk, bulge in points:
        gas_only += gas * abs(gas) > velocity**2
    for inner, outer in itertools.combinations(points, 2):
        r1, v1, e1, *terms1 = inner
        r2, v2, e2, *terms2 = outer
        error = 2 * math.hypot(r1 * v1 * e1, r2 * v2 * e2)
        total = r1 * v1**2 - r2 * v2**2
        # r1 V_component^2 - r2 V_component^2; a negative V denotes an outward force.
        components = [r1 * a * abs(a) - r2 * b * abs(b) for a, b in zip(terms1, terms2)]
        fixed = sum(x * y for x, y in zip(components, (1, 0.5, 0.7)))
        maximum = sum(
            max(x * low, x * high)
            for x, (low, high) in zip(components, ((0.8, 1.2), (0.2, 0.8), (0.3, 1.0)))
        )
        nominal.append(((total - fixed) / error, r1, r2))
        robust.append(
            ((total - max(scale_low * maximum, scale_high * maximum)) / error, r1, r2)
        )
    return (len(points), len(robust), gas_only, *max(nominal), *max(robust))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path)
    args = parser.parse_args()
    for filename, expected in HASHES.items():
        if (
            hashlib.sha256((args.inputs / filename).read_bytes()).hexdigest()
            != expected
        ):
            raise ValueError(f"Input changed: {filename}")
    meta = {}
    for line in (args.inputs / "SPARC_Lelli2016c.mrt").read_text().splitlines():
        f = line.split()
        if len(f) >= 19 and f[1].isdigit():
            meta[f[0]] = [float(f[i]) for i in (2, 3, 5, 6)] + [int(f[17])]
    with zipfile.ZipFile(args.inputs / "Rotmod_LTG.zip") as z:
        records = []
        for path in sorted(z.namelist()):
            if "/" in path or not path.endswith("_rotmod.dat"):
                continue
            name = path.removesuffix("_rotmod.dat")
            distance, distance_err, inc, inc_err, quality = meta[name]
            points = [
                list(map(float, line.split()[:6]))
                for line in z.read(path).splitlines()
                if line.strip() and not line.startswith(b"#")
            ]
            records.append(
                [
                    name,
                    quality,
                    *scan_galaxy(points, distance, distance_err, inc, inc_err),
                ]
            )
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(
        "galaxy quality radii pairs gas_exceeds_obs nominal_score r1_kpc r2_kpc worst_case_score worst_r1_kpc worst_r2_kpc".split()
    )
    for row in records:
        writer.writerow([round(a, 6) if isinstance(a, float) else a for a in row])


if __name__ == "__main__":
    main()
