"""Replay the exploratory velocity-channel search on the frozen three middle panels."""
import argparse
import csv
import hashlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from urllib.parse import quote
from urllib.request import urlopen

import numpy as np
from PIL import Image
from scipy.signal import lombscargle

ROOT = Path(__file__).resolve().parent
PALETTE = [(0, 0, 0), (0, 0, 255), (255, 0, 0), (0, 128, 0), (255, 215, 0),
           (0, 180, 175), (176, 0, 176), (255, 155, 0)]
PAGES = ["3rd", "4th", "5th"]


def traces(path):
    pixels = np.asarray(Image.open(path).convert("RGB").resize((1000, 562)), dtype="int16")[82:492, 103:848]
    result = []
    for color in PALETTE:
        mask = np.all(abs(pixels - color) < 40, axis=2)
        x = np.where(mask.sum(axis=0) >= 5)[0]
        if len(x) < 12:
            result.append(None)
            continue
        y = 410 - np.array([np.median(np.where(mask[:, col])[0]) for col in x])
        result.append((4 + x / 745 * 992, y) if np.ptp(np.percentile(y, [5, 95])) >= 5 else None)
    return result


def panel_scores(panels):
    times, values = [], []
    for page, (t, y) in enumerate(panels):
        times.extend(t + page * 1000)
        values.extend((y - np.median(y)) / max(3, np.std(y)))
    t, y = np.array(times), np.array(values)
    bins = np.floor(t / 6).astype(int)
    obs = np.unique(bins)
    t = np.array([np.mean(t[bins == i]) for i in obs])
    y = np.array([np.median(y[bins == i]) for i in obs])
    train = t < 2000
    if min(sum(~train), sum(train) // 2) < 20:
        return None
    y -= np.array([np.median(y[abs(t - center) < 180]) for center in t])
    omega = np.linspace(1 / 340, 1 / 38, 350) * 2 * np.pi
    power = lombscargle(t[train], y[train], omega, normalize=True)
    peak = omega[np.argmax(power)]
    design = np.column_stack((np.sin(peak * t), np.cos(peak * t)))
    model = design @ np.linalg.lstsq(design[train], y[train], rcond=None)[0]
    return 2 * np.pi / peak, np.corrcoef(y[train], model[train])[0, 1], np.corrcoef(y[~train], model[~train])[0, 1]


def inspect(item):
    row, cache = item
    name = row["directory"]
    if not all(row["sha256_" + page] for page in PAGES):
        return []
    panels = []
    for page in PAGES:
        filename = name + "_" + page + ".jpg"
        path = cache / filename
        hasher = row["sha256_" + page]
        if not path.exists():
            url = "https://vlbi.sci.ibaraki.ac.jp/iMet/data/" + quote(name) + "/" + quote(filename)
            with urlopen(url, timeout=45) as response:
                data = response.read()
            if hashlib.sha256(data).hexdigest() != hasher:
                raise ValueError(f"Updated archive panel: {filename}")
            path.write_bytes(data)
        if hashlib.sha256(path.read_bytes()).hexdigest() != hasher:
            raise ValueError(f"Digest mismatch: {filename}")
        panels.append(traces(path))
    result = []
    for index, channel in enumerate(zip(*panels)):
        if any(data is None for data in channel):
            continue
        scores = panel_scores(channel)
        if scores is not None:
            period, train, test = scores
            result.append((min(train, test), name, index, period, train, test))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    with (ROOT / "screen-sources.csv").open() as stream:
        sources = list(csv.DictReader(stream))
    with ProcessPoolExecutor(args.workers) as pool:
        rows = [result for group in pool.map(inspect, [(row, args.cache) for row in sources]) for result in group]
    with args.out.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["min_r", "directory", "color_index", "grid_period_days", "training_r", "selection_r"])
        writer.writerows(sorted(rows, reverse=True))
