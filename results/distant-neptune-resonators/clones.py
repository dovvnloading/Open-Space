"""Follow up the smallest nominal angle arcs with full six-element JPL orbit clones."""

import argparse
import csv
import json
import urllib.parse
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from dynamics import WORK, DAY, AU, GM_SUN, init, add, integrate, angle, arc
from get import save


def setup(r, n, dt):
    solution = json.loads((WORK / (r["spkid"] + ".json")).read_text())
    cov = solution["orbit"]["covariance"]
    assert cov["labels"] == ["e", "q", "tp", "node", "peri", "i"]
    mean = np.array([float(v["value"]) for v in cov["elements"]])
    matrix = np.array(cov["data"], dtype=float)
    sigma = np.sqrt(matrix.diagonal())
    chol = np.linalg.cholesky(matrix / np.outer(sigma, sigma))
    rng = np.random.default_rng(int(r["spkid"]))
    samples = mean + (rng.standard_normal((n, 6)) @ chol.T) * sigma
    epoch = float(cov["epoch"])
    sim, sun = init(epoch, dt)
    for e, q, tp, node, peri, inc in np.vstack([mean, samples]):
        a = q / (1 - e)
        mean_motion = np.sqrt(GM_SUN * DAY**2 / AU**3 / a**3)
        anomaly = np.rad2deg((epoch - tp) * mean_motion)
        add(sim, sun, a, e, inc, node, peri, anomaly)
    return sim


def run(params):
    r, n, dt, years, interval = params
    sim = setup(r, n, dt)
    times, data = integrate(sim, years, interval)
    p, q = int(r["p"]), int(r["q"])
    rows = []
    for clone in range(n + 1):
        phi = angle(data[:, clone], p, q)
        width, center = arc(phi)
        unwrapped = np.rad2deg(np.unwrap(phi))
        rows.append(
            [
                r["designation"],
                r["spkid"],
                p,
                q,
                clone,
                years,
                dt,
                interval,
                round(width, 4),
                round(center, 4),
                round(np.ptp(unwrapped), 4),
            ]
        )
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--objects", nargs="*")
    parser.add_argument("--max-arc", type=float, default=340.0)
    parser.add_argument("--years", type=int, default=1000000)
    parser.add_argument("--step", type=float, default=0.5)
    parser.add_argument("--interval", type=int, default=200)
    parser.add_argument("--clones", type=int, default=8)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--out", default="clone-screen.csv")
    args = parser.parse_args()
    rows = [
        r
        for r in csv.DictReader(open(WORK / "screen.csv"))
        if float(r["arc_deg"]) < args.max_arc
        and (not args.objects or r["spkid"] in args.objects)
    ]
    for r in rows:
        url = "https://ssd-api.jpl.nasa.gov/sbdb.api?" + urllib.parse.urlencode(
            {"sstr": r["spkid"], "cov": "mat", "full-prec": "true"}
        )
        save(url, r["spkid"] + ".json")
    with (WORK / args.out).open("w") as out, ProcessPoolExecutor(args.workers) as pool:
        writer = csv.writer(out)
        writer.writerow(
            [
                "designation",
                "spkid",
                "p",
                "q",
                "clone",
                "years",
                "step_years",
                "cadence_years",
                "arc_deg",
                "center_deg",
                "unwrapped_range_deg",
            ]
        )
        for batch in pool.map(
            run, [(r, args.clones, args.step, args.years, args.interval) for r in rows]
        ):
            writer.writerows(batch)
            out.flush()
            print(
                batch[0][0], sum(row[-1] < 355 for row in batch), len(batch), flush=True
            )
