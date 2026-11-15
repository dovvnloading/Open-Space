"""Rebuild the MeerKAT flux–delay experiment from the original TPA release."""
import argparse
import csv
from collections import Counter, defaultdict
from decimal import Decimal
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.stats import rankdata, spearmanr


DIGEST = "3c4d3993983299337d9653b728eec13c2354b056aae6deda8434b24e8d7092ef"
TARGETS = ("J1621-5039", "J1840-0815", "J1825-0935", "J1840-0809")
PAIRS = [(0, 7), (1, 6), (2, 5), (0, 3), (4, 7)]
K = 4148.808  # s MHz² / (pc cm⁻³)


def save(path, rows):
    with path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def table(archive, name, kind):
    return np.loadtxt(BytesIO(archive.read(f"tpa_data/{name}.{kind}")), ndmin=2)


def match(archive, name):
    flux = table(archive, name, "fluxtable")
    dm = table(archive, name, "dmtable")
    # Both published tables report the epoch to the same hundredth of a day.
    dex = {t: (value, err) for t, value, err in dm}
    counts_dm = Counter(dm[:, 0])
    counts_flux = Counter(flux[:, 0])
    return sorted((t, f, fe, *dex[t]) for t, f, fe in flux
                  if counts_dm[t] == 1 and counts_flux[t] == 1)


def residual(values, time):
    return values - np.polyval(np.polyfit(time, values, 1), time)


def subbands(archive, name, epochs):
    content = archive.read(f"tpa_data/{name}.tim").decode()
    timing = defaultdict(dict)
    for line in content.splitlines():
        fields = line.split()
        if len(fields) < 8 or not fields[0].startswith("/"):
            continue
        ch = int(fields[fields.index("-chan") + 1])
        assert ch not in timing[fields[0]]
        timing[fields[0]][ch] = (float(fields[1]), Decimal(fields[2]),
                                  float(fields[fields.index("-tobs") + 1]))
    pars = archive.read(f"tpa_data/{name}.par").decode().splitlines()
    period = 1 / next(float(p.split()[1]) for p in pars if p.split()[:1] == ["F0"])
    times, _, _, dms, _ = np.asarray(epochs).T
    selected = {}
    for observation, channels in timing.items():
        if len(channels) != 8:
            continue
        freq, mjd, duration = zip(*(channels[k] for k in range(8)))
        idx = np.argmin(abs(times - float(mjd[0])))
        if abs(times[idx] - float(mjd[0])) > 0.008:
            continue
        assert idx not in selected
        freq = np.array(freq)
        assert np.all(np.diff(freq) > 0)
        # Difference Decimal timestamps before conversion to float. Remove the
        # reference delay and wrap offsets to the nearest pulse rotation.
        delay = np.array([float(z - mjd[-1]) * 86400 for z in mjd])
        model = K * np.median(dms) * (freq**-2 - freq[-1]**-2)
        rem = (delay - model + period / 2) % period - period / 2
        pairs = {(i, j): (rem[i] - rem[j]) / (K * (freq[i]**-2 - freq[j]**-2))
                 for i, j in PAIRS}
        selected[idx] = (observation.rsplit("/", 1)[-1], duration[0], freq, pairs)
    return selected


def design(t, f):
    return np.column_stack((np.ones(len(t)), (t - t.mean()) / 365.25, f))


def confirm(rows):
    t, f, _, d, _ = np.array(rows).T
    x = design(t, f)
    split = len(t) // 2
    bright = d[f >= np.quantile(f, .75)]
    faint = d[f <= np.quantile(f, .25)]
    def test_rms(cols):
        fit = np.linalg.lstsq(x[:split, :cols], d[:split], rcond=None)[0]
        return np.sqrt(np.mean((d[split:] - x[split:, :cols] @ fit)**2))
    return {
        "quartile_difference": np.median(bright) - np.median(faint),
        "raw_rank_r": spearmanr(f, d).statistic,
        "early_rank_r": spearmanr(f[:split], d[:split]).statistic,
        "late_rank_r": spearmanr(f[split:], d[split:]).statistic,
        "later_rms_time": test_rms(2),
        "later_rms_time_flux": test_rms(3),
    }


def slope(x, y):
    return np.linalg.lstsq(x, y, rcond=None)[0][2]


def curvature(rows, groups):
    ii = sorted(groups)
    t, f, *_ = np.array(rows)[ii].T
    x = design(t, f)
    low = np.array([groups[k][3][(0, 3)] for k in ii])
    high = np.array([groups[k][3][(4, 7)] for k in ii])
    y = high - low
    # Resample complete contiguous epoch blocks, preserving the dependence
    # of the two receiver halves measured in the same exposure.
    rng = np.random.default_rng(1715)
    blocks = rng.integers(0, len(t), size=(4000, (len(t) + 3)//4, 1))
    ids = (blocks + np.arange(4)).reshape(4000, -1)[:, :len(t)] % len(t)
    ci = np.quantile([slope(x[q], y[q]) for q in ids], [.025, .975])
    phase = (t - t.mean()) * 2 * np.pi / 365.25
    extra = np.column_stack((x, [groups[k][1] for k in ii], np.sin(phase),
                             np.cos(phase), [groups[k][2] for k in ii]))
    return {"n_complete_bands": len(ii), "low_slope": slope(x, low),
            "high_slope": slope(x, high), "high_minus_low_slope": slope(x, y),
            "block_ci_low": ci[0], "block_ci_high": ci[1],
            "adjusted_slope": slope(extra, y)}


def null_ranks(t, f, d, rng, count=4000):
    """Permute complete four-visit blocks of the already detrended DM ranks."""
    a = rankdata(residual(f, t))
    b = rankdata(residual(d, t))
    a = a - a.mean()
    b = b - b.mean()
    norm = np.linalg.norm(a) * np.linalg.norm(b)
    # A short trailing block must travel intact, visiting every rank exactly once.
    blocks = np.arange((len(t) + 3) // 4)
    order = np.argsort(rng.random((count, len(blocks))), axis=1)
    ids = (order[:, :, None] * 4 + np.arange(4)).reshape(count, -1)
    ids = ids[ids < len(t)].reshape(count, len(t))
    return abs(np.sum(a * b[ids], axis=1) / norm)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("release", type=Path, help="Keith et al. (2024) tpa_data.zip")
    parser.add_argument("outdir", type=Path)
    args = parser.parse_args()
    if sha256(args.release.read_bytes()).hexdigest() != DIGEST:
        raise ValueError("Unexpected input SHA-256")
    args.outdir.mkdir(parents=True, exist_ok=True)
    screen, visits, results = [], [], []
    rng = np.random.default_rng(264918)
    maxima = np.zeros(4000)
    with ZipFile(args.release) as archive:
        names = sorted(p.split("/")[-1].split(".")[0] for p in archive.namelist()
                       if p.startswith("tpa_data/J") and p.endswith(".fluxtable"))
        for name in names:
            rows = match(archive, name)
            item = {"pulsar": name, "matched_epochs": len(rows), "detrended_rank_r": ""}
            if len(rows) >= 20:
                t, f, _, d, _ = np.array(rows).T
                item["detrended_rank_r"] = spearmanr(residual(f, t), residual(d, t)).statistic
                maxima = np.maximum(maxima, null_ranks(t, f, d, rng))
            screen.append(item)
            if name not in TARGETS:
                continue
            groups = subbands(archive, name, rows)
            for k, (t, f, fe, d, de) in enumerate(rows):
                item = dict(pulsar=name, mjd=t, flux_mjy=f, flux_err=fe, dm=d, dm_err=de)
                obs, duration, freq, pairs = groups[k]
                item.update(observation=obs, duration_s=duration)
                for ch, nu in enumerate(freq):
                    item[f"f{ch}_mhz"] = nu
                for i, j in PAIRS:
                    item[f"relative_dm_{i}_{j}"] = pairs[i, j]
                visits.append(item)
            results.append({"pulsar": name, **confirm(rows), **curvature(rows, groups)})
    for item in results:
        value = next(r["detrended_rank_r"] for r in screen if r["pulsar"] == item["pulsar"])
        item["block_null_surveys_ge_rank"] = int(np.sum(maxima >= abs(value)))
    save(args.outdir / "screen.csv", screen)
    save(args.outdir / "visits.csv", visits)
    save(args.outdir / "results.csv", results)


if __name__ == "__main__":
    main()
