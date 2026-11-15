"""Independent numerical replay of the published pulsar forecast comparison (Python 3.10+)."""

import argparse
import csv
import hashlib
import itertools
import json
import math
import random
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from decimal import Decimal
from html.parser import HTMLParser
from pathlib import Path
from statistics import mean, median
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent
URL = "https://www.jb.man.ac.uk/~pulsar/glitches/gTable.html"


ATNF = "https://www.atnf.csiro.au/people/pulsar/psrcat/glitchTbl.html"
JBO_HASH = "e811419f2037080bd2c99c1db30555dd4ee6f462f36b91796908f0fe90aabeb4"
ATNF_HASH = "4e605281a11fb36457a891cf8b0952d7e59f3abc5c91544dbe364f0dbc748c51"
PAPER_HASH = "dad1b82c08c2179da6221fceee8a8b0c946619282c3f96ac36555ea67dcd8a36"
CONTROL = "J0537-6910"
ERA_BOUNDARY = 57000


@dataclass(frozen=True)
class Event:
    row_id: int
    epoch: float
    error: float | None
    size: float | None


class Cells(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = []
        self.current = []
        self.cell = None

    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.current = []
        if tag in ("td", "th"):
            self.cell = []

    def handle_data(self, text):
        if self.cell is not None:
            self.cell.append(text)

    def handle_endtag(self, tag):
        if tag in ("td", "th") and self.cell is not None:
            self.current.append(" ".join("".join(self.cell).split()))
            self.cell = None
        elif tag == "tr" and self.current:
            self.rows.append(self.current)


def table(raw):
    parser = Cells()
    parser.feed(raw.decode("utf-8"))
    parser.close()
    return parser.rows


def number(token):
    return None if token.upper() in ("", "X") else float(token)


def parse_jodrell(raw):
    groups = defaultdict(list)
    ids = []
    for row in table(raw):
        if not row[0].isdigit():
            continue
        if len(row) != 11:
            raise ValueError("Unexpected Jodrell Bank column count")
        source = row[1] if row[2].upper() in ("", "X") else "J" + row[2].removeprefix("J")
        event = Event(int(row[0]), float(row[4]), number(row[5]), number(row[6]))
        ids.append(event.row_id)
        groups[source].append(event)
    if sorted(ids) != list(range(1, 738)) or len(groups) != 226:
        raise ValueError("Source snapshot has missing/duplicate rows or an incorrect pulsar count")
    for events in groups.values():
        events.sort(key=lambda event: (event.epoch, event.row_id))
    return groups


def atnf_intervening_epochs(raw, groups):
    rows = table(raw)
    if rows[0][:3] != ["PSR", "PSRJ", "GLEP (MJD)"] or len(rows[1:]) != 644:
        raise ValueError("ATNF schema or event count differs from the pinned release")
    unmatched = defaultdict(list)
    for row in rows[1:]:
        name, epoch = row[1:3]
        if epoch == "-" or epoch.endswith("[s]") or name not in groups:
            continue
        match = re.fullmatch(r"(\d+(?:\.\d+)?)(?:\(([\d.]+)\))?", epoch)
        if match is None:
            raise ValueError(f"Unexpected ATNF epoch {epoch!r}")
        center, digits = match.groups()
        error = None
        if digits is not None:
            error = (float(digits) if "." in digits else
                     float(Decimal(digits).scaleb(Decimal(center).as_tuple().exponent)))
        center = float(center)
        for event in groups[name]:
            floor = 30 if (error is None or event.error is None) else 3
            tolerance = max(floor, 2 * ((error or 0) + (event.error or 0)))
            if abs(event.epoch - center) <= tolerance:
                break
        else:
            unmatched[name].append(center)
    return unmatched


def documented_gap(name, a, b):
    if name == "J0534+2200":
        return a < 45028 or abs(a - 50489.7) < 1 or abs(b - 50489.7) < 1
    if name == CONTROL:
        return a < ERA_BOUNDARY < b
    if name == "J2229+6114":
        return any(a < epoch < b for epoch in (59868.062, 60063.47))
    return False


def derive_partitions(groups, unmatched):
    partitions = {}
    count = 0
    for name, events in groups.items():
        intervals = []
        for left, right in zip(events, events[1:]):
            if left.size is None or left.size <= 0:
                continue
            wait = right.epoch - left.epoch
            if wait <= (left.error or 0) + (right.error or 0):
                continue
            if documented_gap(name, left.epoch, right.epoch):
                continue
            if any(left.epoch < t < right.epoch for t in unmatched.get(name, ())):
                continue
            intervals.append((left, right, (left.size, wait)))
        count += len(intervals)
        if name == CONTROL:
            train = [i for i in intervals if i[1].epoch < ERA_BOUNDARY]
            test = [i for i in intervals if i[0].epoch > ERA_BOUNDARY]
        else:
            n_test = max(5, math.ceil(len(intervals) / 3))
            train, test = intervals[:-n_test], intervals[-n_test:]
        if len(train) >= 8 and len(test) >= 5:
            partitions[name] = (train, test)
    return partitions, count


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


def verify_screen(rows, groups, partitions):
    listed_names = [r["name"] for r in rows]
    if len(listed_names) != len(set(listed_names)) or set(listed_names) != set(partitions):
        raise ValueError("Published sources disagree with the complete eligible population")
    training = {}
    p_values = []
    for row in rows:
        name = row["name"]
        split = partitions[name]
        for key, intervals in zip(("train_start_ids", "test_start_ids"), split):
            expected = [i[0].row_id for i in intervals]
            if [int(i) for i in row[key].split()] != expected:
                raise ValueError(f"{name} {key} differs from the independently derived chronological split")
        check(name + " events", len(groups[name]), row["n_events"])
        train, test = ([i[2] for i in intervals] for intervals in split)
        check(name + " training", len(train), row["n_train"])
        check(name + " test", len(test), row["n_test"])
        intercept, slope, baseline = estimate(train)
        sizes, waits = zip(*test)
        forecast = predict(intercept, slope, sizes)
        metrics = score_predictions(waits, forecast, baseline)
        for key in ("model_mae_days", "baseline_mae_days", "skill", "min_loo_gain_days"):
            check(name + " " + key, metrics[key], row[key])
        p = permutation_p(name, waits, forecast)
        check(name + " permutation", p, row["p_value"])
        if name != CONTROL:
            p_values.append((p, row, slope))
        elif row["p_adjusted"]:
            raise ValueError("The named control must remain outside the tested family")
        training[name] = train + test
    running = 0.0
    for i, (prob, row, slope) in enumerate(sorted(p_values, key=lambda pr: pr[0])):
        running = max(running, prob * (len(p_values) - i))
        adjusted = min(1., running)
        check(row["name"] + " Holm", adjusted, row["p_adjusted"])
        if slope > 0 and float(row["skill"]) >= .25 and float(row["min_loo_gain_days"]) > 0 and adjusted <= .05:
            raise ValueError(f"Contradicts the published zero promotions: {row['name']} passes")
    return training


def gain_bounds(w_lo, w_hi, p_lo, p_hi, baseline):
    # Optima occur on a rectangle edge or a kink of the absolute-value loss.
    gains = []
    for wait in (w_lo, w_hi, baseline, p_lo, p_hi):
        if w_lo <= wait <= w_hi:
            for prediction in (p_lo, p_hi, min(p_hi, max(p_lo, wait))):
                gains.append(abs(wait - baseline) - abs(wait - prediction))
    return min(gains), max(gains)


def verify_proof(path, events):
    if hashlib.sha256(path.read_bytes()).hexdigest() != PAPER_HASH:
        raise ValueError("The supplied paper is not the pinned PDF proof")
    result = subprocess.run(["pdftotext", "-f", "4", "-l", "4", "-layout", str(path), "-"],
                            check=True, capture_output=True, text=True)
    values = []
    for name, epoch, error, size, size_error in re.findall(
            r"(?m)^\s*(G\d+)\s+(\d+)\((\d+)\)\s+(?:(\d+\.\d)\((\d)\)|-)\s+", result.stdout):
        values.append((name, float(epoch), float(error), float(size) if size else None,
                       float(size_error) / 10 if size_error else None))
    expected = [(e["row_id"], e["mjd"], e["mjd_err"], e["size"], e["size_err"]) for e in events]
    if values != expected:
        raise ValueError("Epoch, size, or uncertainty transcriptions differ from Table 2 of the proof")
    print("Verified all seven event epochs, sizes and uncertainties against the pinned PDF")


def verify_external(data, training, last_catalogue_epoch):
    events = data["events"]
    if [e["ordinal"] for e in events] != list(range(36, 43)):
        raise ValueError("Expected G36 through G42 in order")
    if any(e["row_id"] != f"G{e['ordinal']}" for e in events):
        raise ValueError("External row IDs disagree with the paper's event numbers")
    if data["jbo_training_cutoff_mjd"] != last_catalogue_epoch:
        raise ValueError("Published training cutoff differs from the last catalogue event")
    if not last_catalogue_epoch < data["observing_span_mjd"][0] <= events[0]["mjd"]:
        raise ValueError("The external observing run must start after all catalogue events")
    if events[-1]["mjd"] > data["observing_span_mjd"][1]:
        raise ValueError("The last external event lies outside the observing run")
    for e in events:
        if e["mjd_err"] < 0 or (e["size"] is None) != (e["size_err"] is None):
            raise ValueError("Invalid externally reported event uncertainty")
        if e["size"] is not None and (e["size"] <= 0 or e["size_err"] < 0):
            raise ValueError("Invalid externally reported size or size uncertainty")
    intercept, slope, baseline = estimate(training)
    summary = data["summary"]
    obs, pred, all_bounds = [], [], []
    expected_ids = []
    rejected = []
    for left, right in zip(events, events[1:]):
        ids = (left["row_id"], right["row_id"])
        radius = left["mjd_err"] + right["mjd_err"]
        if right["mjd"] - left["mjd"] <= radius:
            raise ValueError(f"Indeterminate external epoch ordering for {ids}")
        if left["size"] is None or ids == ("G37", "G38"):
            rejected.append(ids)
        else:
            expected_ids.append(ids)
    ids_published = [(p["row_id"], p["next_row_id"]) for p in data["points"]]
    if expected_ids != ids_published:
        raise ValueError("Duplicate, extra, missing or misordered external forecast rows")
    if rejected != [(p["row_id"], p["next_row_id"]) for p in data["excluded"]]:
        raise ValueError("External exclusion list disagrees with source event sequence")
    indexed = {e["row_id"]: e for e in events}
    for row in data["excluded"]:
        left, right = indexed[row["row_id"]], indexed[row["next_row_id"]]
        check("excluded epoch", left["mjd"], row["mjd"])
        check("excluded successor", right["mjd"], row["next_mjd"])
        reasons = ["absent or nonpositive starting size"] if left["size"] is None else []
        if left["row_id"] == "G37" and right["row_id"] == "G38":
            reasons.insert(0, "No continuous timing solution between G37 and G38. Small intervening glitches may have been missed.")
        if reasons != row["reasons"]:
            raise ValueError("External exclusion reasons disagree with the source")
    for point in data["points"]:
        left, right = indexed[point["row_id"]], indexed[point["next_row_id"]]
        wait = right["mjd"] - left["mjd"]
        forecast = predict(intercept, slope, [left["size"]])[0]
        radius = left["mjd_err"] + right["mjd_err"]
        p_lo, p_hi = predict(intercept, slope, [left["size"] - left["size_err"],
                                               left["size"] + left["size_err"]])
        bounds = gain_bounds(wait - radius, wait + radius, p_lo, p_hi, baseline)
        obs.append(wait)
        pred.append(forecast)
        all_bounds.append(bounds)
        expected = dict(wait_days=wait, pred_days=forecast, baseline_days=baseline, mjd=left["mjd"],
                        next_mjd=right["mjd"], mjd_err=left["mjd_err"], next_mjd_err=right["mjd_err"],
                        size=left["size"], size_err=left["size_err"])
        for key, value in expected.items():
            check(key, value, point[key])
        for key, expected_range in (("wait_bounds_days", (wait - radius, wait + radius)),
                                    ("prediction_bounds_days", (p_lo, p_hi)), ("gain_bounds_days", bounds)):
            if len(point[key]) != 2:
                raise ValueError(f"Invalid bounds for {key}")
            for value, supplied in zip(expected_range, point[key]):
                check(key, value, supplied)
        if point["label"] != f"{left['row_id']} → {right['row_id']}":
            raise ValueError("External event label disagrees with the source")
        print(f"{left['row_id']} to {right['row_id']}: {wait:.0f} observed, {forecast:.2f} forecast days")
    metrics = score_predictions(obs, pred, baseline)
    metrics.update(fit_intercept_days=intercept, fit_slope_days=slope, baseline_days=baseline,
                   n_train=len(training), n_test=len(obs), p_value=permutation_p("J1341-6220", obs, pred))
    for key, value in metrics.items():
        check(key, value, summary[key])
    for i in (0, 1):
        check("mean gain bounds", mean(b[i] for b in all_bounds), summary["mean_gain_bounds_days"][i])
    print(f"Parkes MAE: {metrics['model_mae_days']:.2f} vs {metrics['baseline_mae_days']:.2f} days")
    comparison = data["post_hoc_checks"]
    best_constant = median(obs)
    check("retrospective four-outcome median", best_constant, comparison["test_wait_median_days"])
    check("minimum constant loss", mean(absolute_errors(obs, [best_constant] * len(obs))),
          comparison["best_constant_mae_days"])
    g42 = predict(intercept, slope, [events[-1]["size"]])[0]
    check("unscored G42 size forecast", g42, comparison["g42_wait_days"])
    check("unscored G42 forecast epoch", events[-1]["mjd"] + g42, comparison["g42_forecast_mjd"])


def read_snapshot(path, url, digest):
    if path:
        raw = path.read_bytes()
    else:
        with urlopen(Request(url, headers={"User-Agent": "Open-Space catalog research"}), timeout=40) as f:
            raw = f.read()
    if hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError(f"{url} differs from the pinned input. Supply the original snapshot on disk.")
    return raw


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog-html", type=Path, help="Local copy of the pinned Jodrell Bank table")
    parser.add_argument("--atnf-html", type=Path, help="Local copy of the pinned ATNF table")
    parser.add_argument("--paper-pdf", type=Path, help="Also check the Parkes transcription using pdftotext")
    args = parser.parse_args()
    data = json.loads((HERE / "parkes-intervals.json").read_text(encoding="utf-8"))
    if data["jbo_sha256"] != JBO_HASH or data["primary_pdf"]["sha256"] != PAPER_HASH:
        parser.error("The result files and verifier identify different snapshots")
    if args.paper_pdf:
        verify_proof(args.paper_pdf, data["events"])
    else:
        print("Using transcribed Parkes epochs and sizes. Supply --paper-pdf to check the original Table 2.")
    groups = parse_jodrell(read_snapshot(args.catalog_html, URL, JBO_HASH))
    atnf = atnf_intervening_epochs(read_snapshot(args.atnf_html, ATNF, ATNF_HASH), groups)
    splits, kept = derive_partitions(groups, atnf)
    if kept != 483 or len(splits) != 7:
        parser.error(f"The declared coverage differs from the independently derived {kept} intervals and {len(splits)} eligible sources")
    with (HERE / "screen.csv").open(encoding="utf-8", newline="") as f:
        summary = list(csv.DictReader(f))
    populations = verify_screen(summary, groups, splits)
    print(f"Verified {kept} eligible intervals, {len(splits)} source partitions, six Holm results and zero promotions")
    verify_external(data, populations["J1341-6220"], groups["J1341-6220"][-1].epoch)


if __name__ == "__main__":
    main()
