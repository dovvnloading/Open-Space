"""Compare the 2037 Sinope eclipse across archived JPL satellite orbit solutions."""

import argparse
from datetime import datetime, timedelta
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import search
from spk_slice import Kernel, evaluate

J2000 = datetime(2000, 1, 1, 12)
EPOCH_SECONDS = (search.EPOCH - J2000).total_seconds()
COMPARISON_TIME = datetime(2037, 10, 17, 20)
NAMES = ("partial_ingress", "total_ingress", "total_egress", "partial_egress")


def days(stamp):
    return (stamp - search.EPOCH).total_seconds() / 86400


def contacts(states):
    mid = days(datetime(2037, 10, 17, 19))
    intervals = [
        (mid - 0.6, mid, 0),
        (mid - 0.6, mid, 1),
        (mid, mid + 0.57, 1),
        (mid, mid + 0.57, 0),
    ]
    out = {}
    for key, (start, end, component) in zip(NAMES, intervals):

        def gap(t):
            return search.metrics(states, 509, 599, t)[component]

        t = brentq(gap, start, end)
        out[key] = search.utc(search.earth_time(states, 509, t))
    return out


def relative_state(kernel, body):
    @lru_cache(maxsize=None)
    def at_epoch(second):
        center, frame, kind, address, record = kernel.coefficients(body, second)
        if (center, frame) != (5, 1) or kind not in (2, 3):
            raise ValueError(
                "Expected ICRF Jupiter-system-barycentric type 2 or type 3 segment"
            )
        n = (len(record) - 2) // (3 if kind == 2 else 6)
        return record[: 2 + 3 * n]

    # Reuse the same position coefficients across the entire event where possible.
    existing = []

    def state(seconds):
        for record in existing:
            if abs(seconds - record[0]) <= record[1]:
                return evaluate(record, seconds)
        record = at_epoch(seconds)
        existing.append(record)
        return evaluate(record, seconds)

    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--horizons-cache", type=Path, required=True, help="Cache used by ../search.py"
    )
    parser.add_argument(
        "--spk-cache",
        type=Path,
        required=True,
        help="Cache for short NAIF HTTP byte ranges",
    )
    parser.add_argument("--output", type=Path, default=HERE / "orbit-comparison.json")
    args = parser.parse_args()
    states = {}
    for body in (10, 399, 509, 599):
        states[body] = search.state_table(args.horizons_cache, body)[0]
    reference = contacts(states)
    out = {
        "reference": "JUP347 (live Horizons Sinope orbit at original search)",
        "reference_contacts_utc": reference,
        "state_comparison_epoch_tdb": COMPARISON_TIME.isoformat(),
        "solutions": {},
    }
    for key in ("jup300", "jup343", "jup344"):
        kernel = Kernel(key, args.spk_cache)
        old_sat = relative_state(kernel, 509)
        old_jupiter = relative_state(kernel, 599)

        def moon(t):
            sec = EPOCH_SECONDS + t * 86400
            return states[599](t) + old_sat(sec) - old_jupiter(sec)

        trial = dict(states)
        trial[509] = moon
        times = contacts(trial)
        shift = {
            k: round(
                (
                    datetime.fromisoformat(times[k])
                    - datetime.fromisoformat(reference[k])
                ).total_seconds()
            )
            for k in NAMES
        }
        diff = moon(days(COMPARISON_TIME)) - states[509](days(COMPARISON_TIME))
        ranges = [
            dict(start=s, length=n, sha256=sha256(kernel.read(s, n)).hexdigest())
            for s, n in sorted(kernel.used)
        ]
        out["solutions"][key] = {
            "kernel": kernel.url,
            "release_notes": kernel.url.replace(".bsp", ".cmt"),
            "contacts_utc": times,
            "contact_minus_reference_seconds": shift,
            "jovicentric_position_minus_reference_km": np.round(diff, 3).tolist(),
            "jovicentric_position_difference_km": round(float(np.linalg.norm(diff)), 3),
            "ranges": ranges,
        }
    args.output.write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
