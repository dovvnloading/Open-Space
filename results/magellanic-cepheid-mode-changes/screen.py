"""Screen catalogued OGLE F/1O and 1O/2O Cepheids for a reversal of mode dominance."""

import csv
import io
import tarfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs"


def catalog():
    stars = {}
    for cloud in ["lmc", "smc"]:
        for pair in ["F1O", "1O2O"]:
            for line in (INPUT / f"{cloud}4_cep{pair}.dat").read_text().splitlines():
                if line.strip():
                    stars[line[:17]] = (pair, float(line[33:44]), float(line[100:110]))
    return stars


def read_curves(stars):
    for cloud in ["lmc", "smc"]:
        for phase in [3, 4]:
            with tarfile.open(INPUT / f"{cloud}{phase}_phot.tar.gz", "r|gz") as archive:
                for member in archive:
                    name = Path(member.name)
                    if name.stem not in stars or name.suffix != ".dat":
                        continue
                    band = name.parent.name
                    if band not in ["I", "V"]:
                        continue
                    raw = archive.extractfile(member).read()
                    data = np.loadtxt(io.BytesIO(raw), ndmin=2)
                    valid = (
                        np.isfinite(data[:, :3]).all(axis=1)
                        & (data[:, 2] > 0)
                        & (data[:, 2] < 0.2)
                    )
                    yield (name.stem, phase, band), data[valid, :3]


def basis(t, p, q):
    a = 2 * np.pi * (t - np.median(t)) / p
    b = 2 * np.pi * (t - np.median(t)) / q
    args = [a, b, 2 * a, 2 * b, a + b, a - b, 3 * a]
    return np.column_stack(
        [np.ones(len(t))] + [fun(x) for x in args for fun in [np.sin, np.cos]]
    )


def estimate(data, p, q):
    t, y, e = data.T
    X = basis(t, p, q)
    w = 1 / np.hypot(e, 0.007)
    use = np.full(len(t), True)
    for _ in range(3):
        Xw = X[use] * w[use, None]
        beta, _, rank, _ = np.linalg.lstsq(Xw, y[use] * w[use], rcond=None)
        if rank < X.shape[1] or np.linalg.cond(Xw) > 40:
            return None
        r = (y - X @ beta) * w
        cut = 5 * max(1, 1.4826 * np.median(abs(r[use] - np.median(r[use]))))
        new = abs(r - np.median(r[use])) < cut
        if new.sum() < X.shape[1] + 7:
            return None
        if np.array_equal(new, use):
            break
        use = new
    return np.hypot(beta[1], beta[2]), np.hypot(beta[3], beta[4])


def window(data, which):
    t = data[:, 0]
    return data[t < t.min() + 900] if which == "first" else data[t > t.max() - 900]


def main():
    stars = catalog()
    curves = dict(read_curves(stars))
    out = []
    for name, (pair, p, q) in sorted(stars.items()):
        a = curves.get((name, 3, "I"))
        b = curves.get((name, 4, "I"))
        if a is None or b is None:
            continue
        a = window(a, "first")
        b = window(b, "last")
        if min(len(a), len(b)) < 55 or min(np.ptp(a[:, 0]), np.ptp(b[:, 0])) < 350:
            continue
        x = estimate(a, p, q)
        y = estimate(b, p, q)
        if x is None or y is None:
            continue
        d1 = x[0] - x[1]
        d2 = y[0] - y[1]
        # Positive only if the dominant mode reverses, with both endpoints away from equality.
        score = min(abs(d1), abs(d2)) * np.sign(-d1 * d2)
        out.append((name, pair, p, q, len(a), len(b), *x, *y, score))
    out.sort(key=lambda r: -r[-1])
    with (ROOT / "screen.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "star",
                "mode_pair",
                "period_long_d",
                "period_short_d",
                "n_early",
                "n_late",
                "early_A_long",
                "early_A_short",
                "late_A_long",
                "late_A_short",
                "rank_score_mag",
            ]
        )
        w.writerows(out)
    print(
        "catalog",
        len(stars),
        "screened",
        len(out),
        "positive",
        sum(r[-1] > 0 for r in out),
    )
    for row in out[:12]:
        print(row[0], row[1], *[round(float(x), 4) for x in row[6:]])


if __name__ == "__main__":
    main()
