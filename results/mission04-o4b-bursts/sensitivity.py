"""Measure incident burst sensitivity by injecting into frozen survey strain.

Requires completed survey foreground ledgers. From this directory:
OPENBLAS_NUM_THREADS=1 python sensitivity.py CACHE --out results
"""
import argparse
import csv
import gzip
import io
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta, chi2

import acquire
import core
import vetting

HERE = Path(__file__).resolve().parent
SEED = 406150
CASES = (
    ("primary", 150., 1e-21, 2048),
    ("lower", 150., 2e-22, 512),
    ("middle", 150., 5e-22, 512),
    ("higher_frequency", 350., 1e-21, 512),
)
THRESHOLD = 20.0
Q = 9.
# Earth-fixed arm unit vectors and site coordinates (metres) from
# https://lscsoft.docs.ligo.org/lalsuite/7.26/lal/_l_a_l_detectors_8h_source.html
X = np.array([[-.22389266154, .79983062746, .55690487831],
              [-.95457412153, -.14158077340, -.26218911324]])
Y = np.array([[-.91397818574, .02609403989, -.40492342125],
              [.29774156894, -.48791033647, -.82054461286]])
SITES = np.array([[-2161414.92636, -3834695.17889, 4600350.22664],
                  [-74276.0447238, -5496283.71971, 3224257.01744]])


def response(z, azimuth, psi):
    """Sky direction points toward source. A nearer site sees the wave first."""
    xy = np.sqrt(1 - z*z)
    ca, sa = np.cos(azimuth), np.sin(azimuth)
    direction = np.array([xy*ca, xy*sa, z])
    theta, phi = np.array([z*ca, z*sa, -xy]), np.array([-sa, ca, 0.])
    cp, sp = np.cos(psi), np.sin(psi)
    p, q = cp*theta + sp*phi, -sp*theta + cp*phi
    fplus = ((X@p)**2 - (X@q)**2 - (Y@p)**2 + (Y@q)**2) / 2
    delay = -(SITES[1] - SITES[0]) @ direction / 299792458.
    return float(fplus[0]), float(fplus[1]), float(delay)


def plan(selection, out):
    blocks = [b for b in selection["blocks"] if b["role"] == "search"]
    parent, time, score = [], [], []
    for b in blocks:
        with gzip.open(out / "epochs" / f"{b['gps']}.csv.gz", "rt") as s:
            for row in csv.DictReader(s):
                parent.append(b["gps"])
                time.append(int(row["gps_second"]) - b["gps"])
                score.append(float(row["score"]))
    rng = np.random.default_rng(SEED)
    draws = {b["gps"]: [] for b in blocks}
    for case, f, amp, count in CASES:
        for i in range(count):
            idx = int(rng.integers(len(parent)))
            subsecond = rng.random()
            z = rng.uniform(-1, 1)
            az, psi, phase = rng.random(3) * [2*np.pi, np.pi, 2*np.pi]
            fh, fl, lag = response(z, az, psi)
            draws[parent[idx]].append(dict(case=case, trial=i, block=parent[idx], second=time[idx],
                offset=subsecond, frequency=f, hrss=amp, cos_theta=z, azimuth=az,
                polarization=psi, phase=phase, F_H=fh, F_L=fl, delay_s=lag, baseline=score[idx]))
    return blocks, draws, len(time), sum(v >= THRESHOLD for v in score)


def run_block(cache, block, draws):
    starts = block["gps"]
    inputs = []
    for detector in ("H1", "L1"):
        gps, strain, good = acquire.load_detector(cache / block["files"][detector])
        if gps != starts:
            raise ValueError("Wrong block")
        good = vetting.mask_seconds(gps, good)
        inputs.append((strain, good, core.prepare(strain, acquire.FS, good)))

    def scan_second(arrays, start, stop, sec):
        prepared = [core.prepare(a, acquire.FS, good[start:stop],
                    reference=ref, first_second=start)
                    for a, (_, good, ref) in zip(arrays, inputs)]
        if not all(p.live[sec - start] for p in prepared):
            raise ValueError("Draw outside survey exposure")
        result = core.scan(*prepared)[0]
        return max((t.score for t in result.triggers if int(t.time) == sec), default=0.)

    for draw in draws:
        sec = draw["second"]
        t = sec + draw["offset"]
        start, stop = max(0, sec - 2), min(acquire.DURATION, sec + 3)
        base = [raw[start*acquire.FS:stop*acquire.FS] for raw, good, ref in inputs]
        if draw["trial"] % 32 == 0:
            orig = scan_second(base, start, stop, sec)
            if not np.isclose(orig, draw["baseline"], rtol=1e-4, atol=1e-4):
                raise ValueError(f"Crop replay disagrees with ledger: {orig}, {draw}")
        arrays = []
        for raw, antenna, delta in zip(base, (draw["F_H"], draw["F_L"]), (0., draw["delay_s"])):
            h = core.sine_gaussian(acquire.FS, stop-start, t-start+delta,
                                   draw["frequency"], Q, draw["hrss"], draw["phase"])
            arrays.append(raw + antenna*h)
        score = scan_second(arrays, start, stop, sec)
        draw["injected_score"] = score
        draw["recovered"] = int(draw["baseline"] < THRESHOLD <= score)
    print("injections", starts, len(draws), "completed", flush=True)


def summarize(rows, live, count):
    out = {"seed": SEED, "threshold": THRESHOLD, "Q": Q, "HL_live_s": live,
           "foreground_seconds_ge_threshold": count, "confidence_budget": {
               "efficiency_one_sided_alpha": .05, "count_one_sided_alpha": .05},
           "cases": {}}
    upper = chi2.ppf(.95, 2*(count+1)) / 2
    for case, freq, hrss, n in CASES:
        subset = [d for d in rows if d["case"] == case]
        k = sum(d["recovered"] for d in subset)
        assert len(subset) == n
        eps = beta.ppf(.05, k, n+1-k) if k else 0.
        out["cases"][case] = {"f_Hz": freq, "hrss_Hz_mhalf": hrss, "n": n, "recovered": k,
            "efficiency": k/n, "lower_efficiency_95": eps, "poisson95_count_upper": upper,
            "effective_live_s_lower": live*eps,
            "rate_upper_per_day_pointwise90": upper*86400/(live*eps) if eps>0 else None}
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache", type=Path)
    p.add_argument("--out", type=Path, default=HERE / "results")
    args = p.parse_args()
    blocks, draws, live, count = plan(json.loads((args.cache/"search_selection.json").read_text()), args.out)
    for b in blocks:
        run_block(args.cache, b, draws[b["gps"]])
    rows = [r for b in blocks for r in draws[b["gps"]]]
    with (args.out / "injections.csv.gz").open("wb") as raw, gzip.GzipFile(filename="", fileobj=raw,
          mode="wb", mtime=0) as compressed, io.TextIOWrapper(compressed, newline="") as out:
        w = csv.DictWriter(out, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    result = summarize(rows, live, count)
    (args.out / "sensitivity.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
