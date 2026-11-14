"""Independent numerical replay of the published pulsar forecast comparison (Python 3.10+)."""

import argparse
import csv
import hashlib
import itertools
import json
import math
import random
from collections import defaultdict
from html.parser import HTMLParser
from pathlib import Path
from statistics import mean, median
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent
URL = "https://www.jb.man.ac.uk/~pulsar/glitches/gTable.html"


class Cells(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = []
        self.current = []
        self.cell = None

    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.current = []
        if tag == "td":
            self.cell = []

    def handle_data(self, text):
        if self.cell is not None:
            self.cell.append(text)

    def handle_endtag(self, tag):
        if tag == "td" and self.cell is not None:
            self.current.append("".join(self.cell).strip())
            self.cell = None
        elif tag == "tr" and self.current and self.current[0].isdigit():
            self.rows.append(self.current)


def training_and_test_pairs(raw):
    parser = Cells()
    parser.feed(raw.decode("utf-8"))
    parser.close()
    stars = defaultdict(list)
    for r in parser.rows:
        source = r[1] if r[2] == "X" else "J" + r[2].removeprefix("J")
        size = None if r[6].upper() == "X" else float(r[6])
        stars[source].append((float(r[4]), int(r[0]), size))
    pairs = {}
    for source, events in stars.items():
        ordered = sorted(events)
        for (a, row_id, size), (b, _, _) in zip(ordered, ordered[1:]):
            pairs[source, row_id] = (size, b - a)
    return pairs


def estimate(points):
    slopes = [(b[1] - a[1]) / (b[0] - a[0]) for a, b in itertools.combinations(points, 2)
              if b[0] != a[0]]
    slope = max(0.0, median(slopes)) if slopes else 0.0
    return median(w - slope * x for x, w in points), slope, median(w for x, w in points)


def predict(intercept, slope, sizes):
    return [max(0.0, intercept + slope * x) for x in sizes]


def absolute_errors(observed, predicted):
    return [abs(a - b) for a, b in zip(observed, predicted)]


def permutation_p(name, waits, pred):
    if len(set(pred)) == 1:
        return 1.0
    loss = sum(absolute_errors(waits, pred))
    threshold = loss + 1e-11 * max(1.0, loss)
    if len(waits) <= 8:
        return sum(sum(absolute_errors(p, pred)) <= threshold
                   for p in itertools.permutations(waits)) / math.factorial(len(waits))
    scrambled = list(waits)
    rng = random.Random("pulsar-glitch-protocol-v1:" + name)
    successes = 0
    for _ in range(19999):
        rng.shuffle(scrambled)
        successes += sum(absolute_errors(scrambled, pred)) <= threshold
    return (successes + 1) / 20000


def check(label, observed, expected):
    if not math.isclose(observed, float(expected), rel_tol=1e-10, abs_tol=1e-9):
        raise ValueError(f"{label}: calculated {observed}, published {expected}")


def score_predictions(obs, pred, baseline):
    losses = absolute_errors(obs, pred)
    baseline_losses = absolute_errors(obs, [baseline] * len(obs))
    gains = [b - a for a, b in zip(losses, baseline_losses)]
    return {"model_mae_days": mean(losses), "baseline_mae_days": mean(baseline_losses),
            "skill": 1 - mean(losses) / mean(baseline_losses), "mean_gain_days": mean(gains),
            "min_loo_gain_days": min(mean(gains[:i] + gains[i+1:]) for i in range(len(gains))),
            "wins": sum(g > 0 for g in gains)}


def verify_screen(rows, pairs):
    training = {}
    p_values = []
    for row in rows:
        name = row["name"]
        train, test = ([pairs[name, int(i)] for i in row[key].split()] for key in
                       ("train_start_ids", "test_start_ids"))
        check(name, len(train), row["n_train"])
        check(name, len(test), row["n_test"])
        intercept, slope, baseline = estimate(train)
        sizes, waits = zip(*test)
        forecast = predict(intercept, slope, sizes)
        metrics = score_predictions(waits, forecast, baseline)
        for key in ("model_mae_days", "baseline_mae_days", "skill", "min_loo_gain_days"):
            check(name + " " + key, metrics[key], row[key])
        p = permutation_p(name, waits, forecast)
        check(name + " permutation", p, row["p_value"])
        if name != "J0537-6910":
            p_values.append((p, row))
        training[name] = train + test
    running = 0.0
    for i, (prob, row) in enumerate(sorted(p_values, key=lambda pr: pr[0])):
        running = max(running, prob * (len(p_values) - i))
        check(row["name"] + " Holm", min(1., running), row["p_adjusted"])
    return training


def gain_bounds(w_lo, w_hi, p_lo, p_hi, baseline):
    # Optima occur on a rectangle edge or a kink of the absolute-value loss.
    gains = []
    for wait in (w_lo, w_hi, baseline, p_lo, p_hi):
        if w_lo <= wait <= w_hi:
            for prediction in (p_lo, p_hi, min(p_hi, max(p_lo, wait))):
                gains.append(abs(wait - baseline) - abs(wait - prediction))
    return min(gains), max(gains)


def verify_external(data, training):
    events = data["events"]
    if [e["ordinal"] for e in events] != list(range(36, 43)):
        raise ValueError("Expected G36 through G42 in order")
    if min(e["mjd"] - e["mjd_err"] for e in events) <= data["jbo_training_cutoff_mjd"]:
        raise ValueError("The external epochs overlap the training period")
    intercept, slope, baseline = estimate(training)
    summary = data["summary"]
    obs, pred, all_bounds = [], [], []
    points = {(p["row_id"], p["next_row_id"]): p for p in data["points"]}
    for left, right in zip(events, events[1:]):
        if left["size"] is None:
            continue
        point = points[left["row_id"], right["row_id"]]
        wait = right["mjd"] - left["mjd"]
        forecast = predict(intercept, slope, [left["size"]])[0]
        radius = left["mjd_err"] + right["mjd_err"]
        p_lo, p_hi = predict(intercept, slope, [left["size"] - left["size_err"],
                                               left["size"] + left["size_err"]])
        bounds = gain_bounds(wait - radius, wait + radius, p_lo, p_hi, baseline)
        obs.append(wait)
        pred.append(forecast)
        all_bounds.append(bounds)
        check("wait", wait, point["wait_days"])
        check("prediction", forecast, point["pred_days"])
        check("lower gain", bounds[0], point["gain_bounds_days"][0])
        check("upper gain", bounds[1], point["gain_bounds_days"][1])
        print(f"{left['row_id']} to {right['row_id']}: {wait:.0f} observed, {forecast:.2f} forecast days")
    metrics = score_predictions(obs, pred, baseline)
    metrics.update(fit_intercept_days=intercept, fit_slope_days=slope, baseline_days=baseline,
                   n_train=len(training), n_test=len(obs), p_value=permutation_p("J1341-6220", obs, pred))
    for key, value in metrics.items():
        check(key, value, summary[key])
    for i in (0, 1):
        check("mean gain bounds", mean(b[i] for b in all_bounds), summary["mean_gain_bounds_days"][i])
    print(f"Parkes MAE: {metrics['model_mae_days']:.2f} vs {metrics['baseline_mae_days']:.2f} days")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog-html", type=Path, help="Local copy of the pinned Jodrell Bank table")
    args = parser.parse_args()
    data = json.loads((HERE / "parkes-intervals.json").read_text(encoding="utf-8"))
    if args.catalog_html:
        raw = args.catalog_html.read_bytes()
    else:
        with urlopen(Request(URL, headers={"User-Agent": "Open-Space catalog research"}), timeout=40) as f:
            raw = f.read()
    if hashlib.sha256(raw).hexdigest() != data["jbo_sha256"]:
        parser.exit(1, "The table differs from the snapshot. Supply the pinned bytes with --catalog-html.\n")
    with (HERE / "screen.csv").open(encoding="utf-8", newline="") as f:
        summary = list(csv.DictReader(f))
    populations = verify_screen(summary, training_and_test_pairs(raw))
    print(f"Verified historical scores and Holm results across {len(summary)} eligible sources")
    verify_external(data, populations["J1341-6220"])


if __name__ == "__main__":
    main()
