"""Fixed incident 150-Hz burst ensemble for the prospective calendar extension."""
import numpy as np
from scipy.stats import beta, chi2

import acquire
import core
from sensitivity import response, Q, THRESHOLD

SEED = 404150
N = 8192
FREQUENCY = 150.
HRSS = 1e-21


def draw_plan(blocks):
    gps = sorted(b["gps"] for b in blocks)
    if len(gps) != len(set(gps)):
        raise ValueError("Duplicate observing interval")
    rng = np.random.default_rng(SEED)
    output = {g: [] for g in gps}
    for trial in range(N):
        block, sec = divmod(int(rng.integers(len(gps)*acquire.DURATION)), acquire.DURATION)
        offset = rng.random()
        z = rng.uniform(-1, 1)
        az, psi, phase = rng.random(3) * [2*np.pi, np.pi, 2*np.pi]
        fh, fl, delay = response(z, az, psi)
        output[gps[block]].append(dict(trial=trial, second=sec, offset=offset,
            cos_theta=z, azimuth=az, psi=psi, phase=phase, F_H=fh, F_L=fl, delay_s=delay))
    return output


def run(block, inputs, prepared, live, baseline_by_second, draws_for_block):
    rows = []
    for draw in draws_for_block:
        sec = draw["second"]
        record = {**draw, "gps_start": block["gps"], "live": bool(live[sec]),
                  "baseline": float(baseline_by_second[sec]) if live[sec] else None,
                  "injected_score": None, "recovered": False}
        rows.append(record)
        if not live[sec]:
            continue
        start, stop = max(0, sec-2), min(len(live), sec+3)
        span = slice(start*acquire.FS, stop*acquire.FS)
        peaks = []
        for (raw, mask), ref, factor, delta in zip(inputs, prepared,
                  (draw["F_H"], draw["F_L"]), (0., draw["delay_s"])):
            center = sec + draw["offset"] + delta - start
            wave = core.sine_gaussian(acquire.FS, stop-start, center, FREQUENCY, Q, HRSS, draw["phase"])
            peaks.append(core.prepare(raw[span] + factor*wave, acquire.FS, mask[start:stop],
                                      reference=ref, first_second=start))
        if not all(p.live[sec-start] for p in peaks):
            raise ValueError("Live draw outside local search")
        score = max((tr.score for tr in core.scan(*peaks)[0].triggers if int(tr.time)==sec), default=0.)
        record["injected_score"] = score
        record["recovered"] = bool(record["baseline"] < THRESHOLD <= score)
    return rows


def summarize(total_scheduled, count, recovered, n):
    if n != N:
        raise ValueError("Incomplete scheduled-block injection ensemble")
    result = {"seed": SEED, "scheduled_s": total_scheduled, "draws": n, "recovered": recovered,
              "efficiency_scheduled": recovered/n, "foreground_seconds_ge20": count}
    for label, alpha in (("pointwise90", .05), ("two_campaign_simultaneous90", .025)):
        low = float(beta.ppf(alpha, recovered, n+1-recovered)) if recovered else 0.
        high = float(chi2.ppf(1-alpha, 2*(count+1)) / 2)
        result[label] = {"one_sided_alpha_per_term": alpha, "efficiency_lower": low,
            "effective_time_s_lower": total_scheduled*low, "count_upper": high,
            "rate_per_day_upper": high*86400/(total_scheduled*low) if low>0 else None}
    return result
