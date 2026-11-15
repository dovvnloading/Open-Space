"""Repeat the Phoebe prediction using IMCCE's independent PH20 orbit fit."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from jplephem.spk import SPK
from scipy.optimize import brentq

import rings
import search

KERNEL = "https://ftp.imcce.fr/pub/ephem/satel/phoebe/ph20.bsp"


def run(cache, outdir):
    p = cache / "ph20.bsp"
    if not p.exists():
        with urlopen(
            Request(KERNEL, headers={"User-Agent": "Mozilla/5.0"}), timeout=180
        ) as response:
            p.write_bytes(response.read())
    states = {i: search.state_table(cache, i)[0] for i in [6, 10, 399, 609, 699]}
    barycentre, receipt = search.state_table(cache, 6)
    tm = rings.MID
    with SPK.open(p) as kernel:
        segment = next(
            s
            for s in kernel.segments
            if s.target == 609
            and s.center == 6
            and s.start_jd < search.J0 + tm - 2
            and s.end_jd > search.J0 + tm + 2
        )
        alternate = dict(states)
        alternate[609] = (
            lambda t: barycentre(t) + np.asarray(segment.compute(search.J0, t))[:3].T
        )
        diff = alternate[609](tm) - states[609](tm)
        first, second = zip(
            *[
                (
                    brentq(lambda t: search.metrics(v, 609, 699, t)[1], tm - 1, tm),
                    brentq(lambda t: search.metrics(v, 609, 699, t)[1], tm, tm + 1),
                )
                for v in [states, alternate]
            ]
        )
        prof = np.loadtxt(cache / "PS1P01.TAB", delimiter=",")
        profile = (prof[0, 0], 50, np.maximum(prof[:, 3], 0))
        uv = rings.disk(14)
        minima = []
        for times in [
            np.arange(-21 * 60, -16 * 60 + 1, 5),
            np.arange(15 * 60, 20 * 60 + 1, 5),
        ]:
            d = []
            for v in [states, alternate]:
                f = np.array(
                    [rings.sunlight(v, tm + m / 1440, uv, [profile])[1] for m in times]
                )
                i = int(np.argmin(f))
                d.append(
                    {
                        "utc": search.utc(
                            search.earth_time(v, 609, tm + times[i] / 1440)
                        ),
                        "direct_light_fraction": round(float(f[i]), 6),
                    }
                )
            minima.append({"SAT441": d[0], "PH20": d[1]})
        result = {
            "kernel_url": KERNEL,
            "kernel_sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            "segment_start_jd": segment.start_jd,
            "segment_end_jd": segment.end_jd,
            "barycentre_input": receipt,
            "position_difference_km": round(float(np.linalg.norm(diff)), 2),
            "totality_contact_shifts_seconds": [
                (
                    (
                        search.earth_time(alternate, 609, b)
                        - search.earth_time(states, 609, a)
                    )
                    * 86400
                )
                for a, b in [first, second]
            ],
            "ring_minima": minima,
        }
    (outdir / "orbit-check.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--outdir", type=Path, default=search.DIRECTORY)
    args = p.parse_args()
    run(args.cache, args.outdir)
