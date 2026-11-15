"""Fixed-epoch seasonal sensitivity and a separate Hipparcos aperture check."""

from collections import defaultdict
from statistics import median

import numpy as np
from replay import HERE, load_rows


def seasonal(maxima):
    for name, low, high in [
        ("V Del", -1, 25),
        ("S Vir", -6, 18),
        ("S Cyg", 19, 40),
        ("S Vir", 20, 35),
    ]:
        subset = [
            r for r in maxima if r["star"] == name and low <= int(r["cycle"]) <= high
        ]
        t = np.array([int(r["jd"]) for r in subset])
        mag = np.array([float(r["magnitude"]) for r in subset])
        odd = np.array([int(r["cycle"]) % 2 for r in subset])
        phase = 2 * np.pi * t / 365.2422
        for harmonics in (0, 1, 2):
            columns = [np.ones(len(t)), odd, (t - t.mean()) / 365.2422]
            for k in range(1, harmonics + 1):
                columns.extend((np.cos(k * phase), np.sin(k * phase)))
            fit = np.linalg.lstsq(np.column_stack(columns), mag, rcond=None)[0]
            print(
                name,
                low,
                high,
                "seasonal harmonics",
                harmonics,
                "odd minus even",
                round(fit[1], 4),
            )


def hipparcos(maxima):
    cycles = {int(r["cycle"]): int(r["jd"]) for r in maxima if r["star"] == "S Vir"}
    data = np.loadtxt(HERE / "hip-66100.dat", delimiter="|", skiprows=1)
    groups = defaultdict(list)
    for reduced_jd, mag, error, flag in data:
        if flag > 1:
            continue
        jd = reduced_jd + 2440000
        cycle = min(cycles, key=lambda k: abs(jd - cycles[k]))
        offset = jd - cycles[cycle]
        if abs(offset) <= 25:
            groups[cycle, round(offset)].append((offset, mag))
    for (cycle, key), points in sorted(groups.items()):
        print(
            "S Vir Hipparcos",
            cycle,
            "offset",
            round(median(p[0] for p in points), 2),
            "Hp",
            round(median(p[1] for p in points), 4),
            "transits",
            len(points),
        )


if __name__ == "__main__":
    maxima = load_rows("comparison-maxima.csv")
    seasonal(maxima)
    hipparcos(maxima)
