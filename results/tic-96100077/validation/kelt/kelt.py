"""Predict held-out KELT nights with the TESS amplitude and Gaia orbit frozen."""

import argparse
import hashlib
import json
import sys
import re
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.stats import t as student

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))
from measure import PERIOD, REF, gaia_delay

SCALE = -0.80164096
PSEC = PERIOD * 86400


def load(folder, fetch):
    catalog = pd.read_csv(ROOT / "catalog.csv")
    records = []
    for name, source in json.loads((ROOT / "sources.json").read_text()).items():
        path = folder / name
        if not path.exists() and fetch:
            folder.mkdir(parents=True, exist_ok=True)
            data = urlopen(source["url"], timeout=90).read()
            if hashlib.sha256(data).hexdigest() != source["sha256"]:
                raise ValueError(f"Archive has changed: {name}")
            path.write_bytes(data)
        if hashlib.sha256(path.read_bytes()).hexdigest() != source["sha256"]:
            raise ValueError(f"Input mismatch: {name}")
        orientation, reduction = name.split("_")[5:7]
        meta = catalog.query(
            "kelt_orientation == @orientation and proc_type == @reduction"
        ).iloc[0]
        d = pd.read_csv(path, sep=r"\s+", skiprows=15, names=["bjd", "mag", "err"])
        if reduction == "tfa":
            # IPAC TFA files have the same timestamps as raw, plus a constant.
            d.bjd -= 2450000
        header_count = int(re.search(r"NUMRECORDS = (\d+)", path.read_text()).group(1))
        if (
            len(d) != header_count
            or header_count != meta.npts
            or abs(d.mag.mean() - meta["mean"]) > 2e-6
            or abs(d.mag.min() - meta.minvalue) > 2e-6
            or abs(d.mag.max() - meta.maxvalue) > 2e-6
        ):
            raise ValueError(f"Source photometry differs from catalog metadata: {name}")
        if (
            abs(d.bjd.min() - meta.bjdstart) > 2e-6
            or abs(d.bjd.max() - meta.bjdstop) > 2e-6
        ):
            raise ValueError(f"Time span/count differs from archive metadata: {name}")
        d["orientation"], d["reduction"] = orientation, reduction
        # Fixed reference magnitudes from the raw archive metadata.
        ref = (
            catalog.query("kelt_orientation == @orientation and proc_type == 'raw'")
            .iloc[0]
            .kelt_mag
        )
        d["flux"] = 10 ** (-0.4 * (d.mag - ref))
        records.append(d)
    table = pd.concat(records, ignore_index=True)
    for _, group in table.groupby("orientation"):
        raw = np.sort(group.loc[group.reduction == "raw", "bjd"])
        tfa = group.loc[group.reduction == "tfa", "bjd"].to_numpy()
        idx = np.searchsorted(raw, tfa)
        distance = np.minimum(
            abs(raw[np.clip(idx, 0, len(raw) - 1)] - tfa),
            abs(raw[np.clip(idx - 1, 0, len(raw) - 1)] - tfa),
        )
        if distance.max() > 2e-6:
            raise ValueError("A corrected TFA timestamp has no matching raw exposure")
    return table


def basis(phase):
    harmonic = 2 * np.pi * np.arange(1, 5)
    angle = phase[:, None] * harmonic
    c, s = np.cos(angle), np.sin(angle)
    return (
        np.column_stack((np.ones(len(phase)), c, s)),
        np.column_stack((np.zeros(len(phase)), -s * harmonic, c * harmonic)),
    )


def fit(data, train, degree, delay, scale):
    # At a given time shift each orientation has independent linear Fourier terms.
    t = data.bjd.to_numpy() - 2457000
    phase = ((t - REF) / PERIOD) % 1
    time = (t - np.mean(t[train])) / 1000
    drift = [time**k for k in range(1, degree + 1)]
    if scale is None:
        drift.append(delay / PSEC)
    else:
        phase -= scale * delay / PSEC
    x = np.column_stack(drift)
    y = data.flux.to_numpy()
    groups = [
        data.orientation.to_numpy() == orientation for orientation in ("east", "west")
    ]
    beta = np.zeros(x.shape[1])
    for step in range(30):
        prediction, slope = np.empty_like(y), np.empty_like(y)
        shape, derivative = basis(phase - x @ beta)
        templates = []
        for group in groups:
            use = train & group
            coef = np.linalg.lstsq(shape[use], y[use], rcond=None)[0]
            prediction[group] = shape[group] @ coef
            slope[group] = derivative[group] @ coef
            templates.append(shape * group[:, None])
        jac = np.column_stack([*templates, -slope[:, None] * x])
        residual = y - prediction
        change = np.linalg.lstsq(jac[train], residual[train], rcond=None)[0][
            -len(beta) :
        ]
        if np.max(abs(x @ change)) < 1e-8:
            break
        beta += change
    else:
        raise RuntimeError("Ephemeris fit failed to converge")
    return residual, beta, jac


def paired_interval(null, rival, blocks):
    # Positive values favor the rival. Paired block resampling conditions on the CV fits.
    loss = null**2 - rival**2
    groups = (
        pd.DataFrame({"loss": loss, "id": blocks})
        .groupby("id")
        .loss.agg(["sum", "size"])
    )
    rng = np.random.default_rng(96100077)
    draws = rng.integers(0, len(groups), (3000, len(groups)))
    means = np.sum(groups["sum"].to_numpy()[draws], axis=1) / np.sum(
        groups["size"].to_numpy()[draws], axis=1
    )
    return [float(v * 1e6) for v in np.quantile(means, [0.025, 0.975])]


def run(data):
    nights = np.floor(data.bjd).astype(int).to_numpy()
    # Interleave complete observing nights across the full baseline, keeping orientations together.
    fold = np.searchsorted(np.unique(nights), nights) % 5
    instant = data.bjd.to_numpy() - 2457000
    blocks = (nights - np.min(nights)) // 21
    shifts = [i / 100 for i in range(100)]
    delays = [gaia_delay(instant, s) for s in shifts]
    summary, controls = [], []
    for degree in (1, 2):

        def predict(delay, scale):
            residual = np.empty(len(data))
            for f in range(5):
                test = fold == f
                residual[test] = fit(data, ~test, degree, delay, scale)[0][test]
            return residual

        zero = predict(delays[0], 0.0)
        fixed = predict(delays[0], SCALE)
        wrong = predict(delays[0], -SCALE)
        rs, beta, jac = fit(
            data, np.ones(len(data), dtype=bool), degree, delays[0], None
        )
        # Cluster sandwich with 21-day blocks protects against intra-visit correlations.
        scores = pd.DataFrame(jac * rs[:, None]).groupby(blocks).sum().to_numpy()
        bread = np.linalg.pinv(jac.T @ jac)
        nb = len(scores)
        cov = (
            bread
            @ (scores.T @ scores)
            @ bread
            * nb
            / (nb - 1)
            * (len(rs) - 1)
            / (len(rs) - jac.shape[1])
        )
        width = student.ppf(0.975, nb - 1) * np.sqrt(cov[-1, -1])

        def mse(v):
            return float(np.mean(v * v) * 1e6)

        loss = [mse(fixed)] + [mse(predict(d, SCALE)) for d in delays[1:]]
        controls.extend(
            dict(degree=degree, shift=s, held_out_mse_ppt2=l)
            for s, l in zip(shifts, loss)
        )
        summary.append(
            dict(
                degree=degree,
                n=len(data),
                nights=len(np.unique(nights)),
                blocks=nb,
                start_bjd=min(data.bjd),
                end_bjd=max(data.bjd),
                fixed_scale=SCALE,
                free_scale=beta[-1],
                free_scale_cluster95=[beta[-1] - width, beta[-1] + width],
                zero_rms_ppt=np.sqrt(mse(zero)),
                fixed_rms_ppt=np.sqrt(mse(fixed)),
                opposite_rms_ppt=np.sqrt(mse(wrong)),
                zero_minus_fixed_mse_ppt2=mse(zero) - mse(fixed),
                night_ci95_ppt2=paired_interval(zero, fixed, nights),
                block21_ci95_ppt2=paired_interval(zero, fixed, blocks),
                other_phases_as_good=int(np.sum(np.array(loss[1:]) <= loss[0])),
            )
        )
    return summary, controls


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--raw-dir", type=Path, default=ROOT / "raw")
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--output", type=Path, default=ROOT)
    args = p.parse_args()
    data = load(args.raw_dir, args.fetch)
    results, phases = {}, []
    for reduction in ("tfa", "raw"):
        results[reduction], rows = run(
            data[data.reduction == reduction].reset_index(drop=True)
        )
        phases.extend(dict(reduction=reduction, **row) for row in rows)
    (args.output / "fits.json").write_text(json.dumps(results, indent=2) + "\n")
    pd.DataFrame(phases).to_csv(
        args.output / "phase-controls.csv", index=False, float_format="%.9f"
    )
