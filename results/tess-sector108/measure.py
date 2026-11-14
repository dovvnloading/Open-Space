"""Measure the frozen-period Sector 108 test. Requires astropy and numpy."""

import csv
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parent


def write_rows(name, rows):
    with (ROOT / name).open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def window(t, f, centre, lo, hi, cadence):
    start, end = np.searchsorted(t, [centre + lo, centre + hi])
    if end - start < 0.65 * (hi - lo) / cadence:
        return None
    return t[start:end], f[start:end]


def measure_events(t, f, epoch, period, duration, cadence, outer):
    events = []
    for n in range(
        math.ceil((t[0] - epoch) / period), math.floor((t[-1] - epoch) / period) + 1
    ):
        centre = epoch + n * period
        left = window(t, f, centre, -outer * duration, -0.8 * duration, cadence)
        right = window(t, f, centre, 0.8 * duration, outer * duration, cadence)
        core = window(t, f, centre, -0.4 * duration, 0.4 * duration, cadence)
        if left is None or right is None or core is None:
            events.append((n, centre, None))
            continue
        baseline = np.interp(
            core[0],
            [np.median(left[0]), np.median(right[0])],
            [np.median(left[1]), np.median(right[1])],
        )
        losses = np.sort((1 - core[1] / baseline) * 1e6)
        trim = math.floor(0.1 * len(losses))
        losses = losses[trim : len(losses) - trim]
        events.append((n, centre, float(np.mean(losses))))
    return events


def mean_or_none(values):
    return float(np.mean(values)) if len(values) else None


def compile_result(events, mid):
    good = [(n, t, d) for n, t, d in events if d is not None]
    a = np.array([d for _, _, d in good])
    n = len(a)
    return {
        "opportunities": len(events),
        "covered": n,
        "positive": int((a > 0).sum()),
        "depth_ppm": mean_or_none(a),
        "event_sem_ppm": float(np.std(a, ddof=1) / np.sqrt(n)) if n >= 2 else None,
        "first_half_ppm": mean_or_none([d for _, t, d in good if t < mid]),
        "second_half_ppm": mean_or_none([d for _, t, d in good if t >= mid]),
    }


def analyze(c, product):
    with fits.open(ROOT / "fits" / product["name"]) as hdus:
        assert (
            hdus[0].header["TICID"] == product["tic"]
            and hdus[0].header["SECTOR"] == 108
        )
        hdus.verify("exception")
        table = hdus[1].data
        assert hdus[1].header["TIMESYS"] == "TDB"
        ref = hdus[1].header["BJDREFI"] + hdus[1].header["BJDREFF"]
        times = np.array(table["TIME"], dtype=float)
        p = float(c["period_d"])
        e = float(c["T0_BJD_TDB"]) - ref
        dur = float(c["duration_h"]) / 24
        cadence = float(hdus[1].header["TIMEDEL"])
        summary = []
        all_events = []
        phases = []
        for col, outer, label in [
            ("PDCSAP_FLUX", 2.0, "primary"),
            ("SAP_FLUX", 2.0, "raw"),
            ("PDCSAP_FLUX", 1.5, "narrow"),
        ]:
            f = np.array(table[col], dtype=float)
            mask = (
                (table["QUALITY"] == 0) & np.isfinite(times) & np.isfinite(f) & (f > 0)
            )
            t, f = times[mask], f[mask]
            mid = (t[0] + t[-1]) / 2
            assert np.all(np.diff(t) > 0)
            ev = measure_events(t, f, e, p, dur, cadence, outer)
            result = compile_result(ev, mid)
            prefix = {"tic": product["tic"], "method": label}
            if label == "primary":
                for n, centre, d in ev:
                    all_events.append(
                        {
                            "tic": product["tic"],
                            "cycle": n,
                            "centre_btjd": round(centre, 7),
                            "depth_ppm": d,
                        }
                    )
            for k in range(1, 200):
                offset = k / 200
                if min(offset, 1 - offset) * p < 2.5 * dur:
                    continue
                ctl = measure_events(t, f, e + offset * p, p, dur, cadence, outer)
                s = compile_result(ctl, mid)
                # Coverage varies with the gap pattern. Do not compare a single event to a stack.
                if s["covered"] == 0 or abs(s["covered"] - result["covered"]) > max(
                    1, round(0.2 * result["covered"])
                ):
                    continue
                phases.append(
                    {
                        **prefix,
                        "phase_offset": offset,
                        "covered": s["covered"],
                        "depth_ppm": s["depth_ppm"],
                    }
                )
            depths = np.array([r["depth_ppm"] for r in phases if r["method"] == label])
            count = len(depths)
            primary = result["depth_ppm"]
            greater = int((depths >= primary).sum()) if primary is not None else None
            summary.append(
                {
                    **prefix,
                    "period_d": p,
                    "predicted_depth_ppm": int(c["depth_ppm"]),
                    "duration_h": float(c["duration_h"]),
                    "n_good_cadences": len(t),
                    "crowdsap": hdus[1].header["CROWDSAP"],
                    **result,
                    "offphase_count": count,
                    "offphase_equal_deeper": greater,
                    "offphase_rms_ppm": (
                        float(np.sqrt(np.mean(depths**2))) if count else None
                    ),
                }
            )
    return summary, all_events, phases


def main():
    receipt = json.loads((ROOT / "provenance.json").read_text())
    with (ROOT / "ephemerides.csv").open() as f:
        catalog = list(csv.DictReader(f))
    results = []
    events = []
    phases = []
    for product in receipt["selected"]:
        for c in [r for r in catalog if int(r["TIC"]) == product["tic"]]:
            summary, e, p = analyze(c, product)
            results.extend(summary)
            events.extend(e)
            phases.extend(p)
    write_rows("results.csv", results)
    write_rows("events.csv", events)
    write_rows("phase_controls.csv", phases)


if __name__ == "__main__":
    main()
