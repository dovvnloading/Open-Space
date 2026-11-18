"""Replay the frozen OSI SAF search and publish cell footprints and controls."""

import argparse
from dataclasses import dataclass
from datetime import date, datetime, timedelta
from hashlib import file_digest, sha256
import json
from pathlib import Path
from urllib.request import urlretrieve

from netCDF4 import Dataset, num2date
import numpy as np
from scipy import ndimage


HERE = Path(__file__).resolve().parent
DATES = [date(2026, 11, 1) + timedelta(days=i) for i in range(17)]
PRODUCTS = {"amsr3": "OSI-408-g", "amsr2": "OSI-408-a", "ssmis": "OSI-401-d"}
NEIGHBORS = np.ones((3, 3), dtype=bool)
VARIANTS = {"primary": (70, 5), "pack60": (60, 5), "pack80": (80, 5), "core3": (70, 3)}
COAST_BITS = 32 | 64 | 128


@dataclass
class Map:
    effective: np.ndarray
    valid: np.ndarray
    offshore: np.ndarray
    uncertainty: np.ndarray


@dataclass
class Opening:
    core: set
    ring: set
    body: set


def locate(root, entry, download):
    path = root / entry["sensor"] / entry["file"]
    if not path.exists() and download:
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_suffix(".download")
        try:
            urlretrieve(entry["url"], part)
            part.replace(path)
        finally:
            part.unlink(missing_ok=True)
    with path.open("rb") as stream:
        digest = file_digest(stream, "sha256").hexdigest()
    if digest != entry["sha256"] or path.stat().st_size != entry["bytes"]:
        raise ValueError(f"Source checksum/size mismatch: {path}")
    return path


def field(dataset, name):
    return dataset[name][:].astype("float64").filled(np.nan)


def screened_map(filtered, raw, flags, uncertainty):
    effective = np.where(flags == 2, np.maximum(filtered, raw), filtered)
    valid = ((flags == 0) | (flags == 2)) & np.isfinite(effective)
    valid &= np.isfinite(filtered) & (filtered >= 0) & (filtered <= 100)
    coast = (flags >= 0) & ((flags & COAST_BITS) != 0)
    offshore = ndimage.distance_transform_edt(~coast) > 5 if coast.any() else np.ones_like(coast)
    return Map(effective, valid, offshore, uncertainty)


def load_current(path, entry, grid):
    with Dataset(path) as dataset:
        if dataset.product_id != PRODUCTS[entry["sensor"]] or dataset.product_version != "4.1":
            raise ValueError(f"Unexpected product: {path}")
        time = dataset["time"]
        bounds = num2date(dataset["time_bnds"][0], time.units, time.calendar)
        start = datetime.fromisoformat(entry["date"])
        expected = (start.isoformat(), (start + timedelta(days=1)).isoformat())
        if tuple(t.isoformat() for t in bounds) != expected:
            raise ValueError(f"Unexpected time bounds: {path}")
        x = dataset["xc"][:].filled(np.nan)
        y = dataset["yc"][:].filled(np.nan)
        projection = dataset["Polar_Stereographic_Grid"].proj4_string
        if (len(x), len(y)) != (790, 830) or not (
            np.all(np.diff(x) == 10) and np.all(np.diff(y) == -10)
        ):
            raise ValueError(f"Unexpected grid axes: {path}")
        signature = sha256(x.astype("<f8").tobytes() + y.astype("<f8").tobytes() + projection.encode()).hexdigest()
        if not grid:
            grid.update(signature=signature, lat=field(dataset, "lat"), lon=field(dataset, "lon"))
        if grid["signature"] != signature:
            raise ValueError(f"Grid mismatch: {path}")
        return screened_map(
            field(dataset, "ice_conc")[0],
            field(dataset, "raw_ice_conc_values")[0],
            dataset["status_flag"][0].filled(-1),
            field(dataset, "total_uncertainty")[0],
        )


def cells(mask, row, col):
    return {(int(i) + row, int(j) + col) for i, j in np.argwhere(mask)}


def detect(frame, pack=70, minimum=5):
    low = frame.valid & frame.offshore & (frame.effective < pack)
    labels, count = ndimage.label(low, structure=NEIGHBORS)
    core = low & (frame.effective <= 30)
    core_counts = np.bincount(labels[core], minlength=count + 1)
    rows, cols = labels.shape
    results = []
    failed = {"edge": 0, "ring": 0}
    for number, (ys, xs) in enumerate(ndimage.find_objects(labels), start=1):
        if core_counts[number] < minimum:
            continue
        if ys.start == 0 or xs.start == 0 or ys.stop == rows or xs.stop == cols:
            failed["edge"] += 1
            continue
        y, x = ys.start - 1, xs.start - 1
        window = np.s_[y:ys.stop + 1, x:xs.stop + 1]
        body = labels[window] == number
        ring = ndimage.binary_dilation(body, structure=NEIGHBORS) & ~body
        compact = frame.valid[window] & frame.offshore[window] & (frame.effective[window] >= pack)
        if not compact[ring].all():
            failed["ring"] += 1
            continue
        results.append(Opening(cells(core[window] & body, y, x), cells(ring, y, x), cells(body, y, x)))
    return results, failed


def values(array, points):
    return array[tuple(np.array(sorted(points)).T)]


def compare(opening, frame):
    if frame is None:
        return {"status": "unavailable"}
    core_ok = values(frame.valid & frame.offshore, opening.core)
    ring_ok = values(frame.valid & frame.offshore, opening.ring)
    core = values(frame.effective, opening.core)[core_ok]
    ring = values(frame.effective, opening.ring)[ring_ok]
    low_count = int(np.count_nonzero(core <= 50))
    compact_count = int(np.count_nonzero(ring >= 70))
    if 5 * core.size < 4 * core_ok.size or 5 * ring.size < 4 * ring_ok.size:
        status = "insufficient"
    elif 2 * low_count >= core.size and 5 * compact_count >= 4 * ring.size:
        status = "supporting"
    else:
        status = "contradictory"
    return {
        "status": status, "valid_core": int(core.size), "valid_ring": int(ring.size),
        "low_core": low_count, "compact_ring": compact_count,
        "core_median_pct": round(float(np.median(core)), 3) if core.size else None,
        "ring_median_pct": round(float(np.median(ring)), 3) if ring.size else None,
    }


def overlaps(first, second):
    smaller, bigger = sorted((first, second), key=len)
    expanded = {(y + dy, x + dx) for y, x in bigger for dy in (-1, 0, 1) for dx in (-1, 0, 1)}
    return 10 * len(smaller & expanded) >= 3 * len(smaller)


def record(opening, grid, frame, comparisons):
    lon = np.deg2rad(values(grid["lon"], opening.core))
    latitude = values(grid["lat"], opening.core)
    unc = values(frame.uncertainty, opening.core)
    unc = unc[np.isfinite(unc)]
    return {
        "center_deg": [
            round(float(np.mean(latitude)), 4),
            round(float(np.rad2deg(np.arctan2(np.sin(lon).sum(), np.cos(lon).sum()))), 4),
        ],
        "core": sorted(opening.core), "ring": sorted(opening.ring), "body_cells": len(opening.body),
        "core_median_pct": round(float(np.median(values(frame.effective, opening.core))), 3),
        "quoted_uncertainty": {
            "available_core_cells": int(unc.size),
            "median_pp": round(float(np.median(unc)), 3) if unc.size else None,
            "max_pp": round(float(np.max(unc)), 3) if unc.size else None,
        },
        "cross_sensor": {key: compare(opening, other) for key, other in comparisons.items()},
    }


def synthetic_controls():
    def cases():
        nominal = np.full((35, 35), 95.0)
        nominal[15:17, 16:19] = 10
        raw = np.full_like(nominal, np.nan)
        flags = np.zeros(nominal.shape, dtype=np.int16)
        yield "sealed", nominal, raw, flags, 1
        opened = nominal.copy()
        opened[:15, 16] = 10
        yield "edge_connected", opened, raw, flags, 0
        gap = flags.copy()
        gap[14, 16] = 256
        yield "missing_ring", nominal, raw, gap, 0
        screened = flags.copy()
        screened[15:17, 16:19] = 2
        restored = raw.copy()
        restored[15:17, 16:19] = 95
        yield "screened_false_hole", nominal, restored, screened, 0
        coast = flags.copy()
        coast[11, 16] = 64
        yield "near_coast", nominal, raw, coast, 0

    test_results = {}
    for label, filtered, raw, flags, expected in cases():
        openings, _ = detect(screened_map(filtered, raw, flags, raw))
        result = len(openings)
        if result != expected:
            raise RuntimeError(f"Control {label}: expected {expected}, got {result}")
        test_results[label] = result
    return test_results


def historical_controls(entries, root, download):
    results = []
    for entry in entries:
        with Dataset(locate(root, entry, download)) as dataset:
            concentration = field(dataset, "z")
            valid = np.isfinite(concentration) & (concentration >= 0) & (concentration <= 100)
            if concentration.shape != (1328, 1264):
                raise ValueError("Unexpected historical ASI grid shape")
            x, y = dataset["x"][:], dataset["y"][:]
            if not (np.all(np.diff(x) == 6250) and np.all(np.diff(y) == 6250)):
                raise ValueError("Unexpected historical ASI grid spacing")
            # Treat every masked or otherwise invalid historical cell as a potential coast.
            offshore = ndimage.distance_transform_edt(valid) > 8
            frame = Map(concentration, valid, offshore, np.full_like(concentration, np.nan))
            openings, _ = detect(frame, minimum=13)
            near_maud = []
            # Approximate projected location of 64.4 S, 4.7 E on EPSG:3412.
            for opening in openings:
                r, c = np.array(sorted(opening.core)).T
                center_km = np.array([float(x[c].mean()), float(y[r].mean())]) / 1000
                if np.linalg.norm(center_km - np.array([230, 2780])) <= 300:
                    near_maud.append({
                        "core_cells": len(opening.core), "body_cells": len(opening.body),
                        "center_xy_km": center_km.round(2).tolist(),
                    })
            results.append({"date": entry["date"], "cutoff_70": near_maud})
    return results


def analyze(entries, root, download):
    synthetic = synthetic_controls()
    history = historical_controls([entry for entry in entries if entry["sensor"] == "history"], root, download)
    indexed = {(entry["sensor"], entry["date"]): entry for entry in entries if entry["sensor"] in PRODUCTS}
    grid = {}
    daily = []
    all_links = {key: [] for key in VARIANTS}
    best_length = {key: 0 for key in VARIANTS}
    previous = {key: [] for key in VARIANTS}
    chain_support = []
    for day in map(str, DATES):
        frames = {
            sensor: load_current(locate(root, indexed[sensor, day], download), indexed[sensor, day], grid)
            if (sensor, day) in indexed else None for sensor in PRODUCTS
        }
        primary = frames["amsr3"]
        comparisons = {key: value for key, value in frames.items() if key != "amsr3"}
        row = {
            "date": day, "sensors": [sensor for sensor, frame in frames.items() if frame is not None],
            "eligible_offshore_cells": int((primary.valid & primary.offshore).sum()),
            "compact_offshore_cells": int((primary.valid & primary.offshore & (primary.effective >= 70)).sum()),
            "variants": {},
        }
        for key, (pack, minimum) in VARIANTS.items():
            observations, rejected = detect(primary, pack, minimum)
            observation_records = []
            current = []
            for index, opening in enumerate(observations, start=1):
                name = f"{day}/{index}"
                details = record(opening, grid, primary, comparisons)
                matched = [item for item in previous[key] if overlaps(opening.core, item[1])]
                length = max((item[2] for item in matched), default=0) + 1
                support = [sensor for sensor, item in details["cross_sensor"].items() if item["status"] == "supporting"]
                histories = [(name,)]
                for other, _, _, paths in matched:
                    all_links[key].append([other, name])
                    histories.extend((*path, name) for path in paths if len(path) < 3)
                current.append((name, opening.core, length, histories))
                for path in histories:
                    if key == "primary" and len(path) == 3:
                        earlier = {item["id"]: item["supporting_sensors"] for date_row in daily[-2:]
                                   for item in date_row["variants"]["primary"]["openings"]}
                        triplet = [earlier[path[0]], earlier[path[1]], support]
                        sensors = [sensor for sensor in comparisons if sum(sensor in tags for tags in triplet) >= 2]
                        if sensors:
                            chain_support.append({"path": list(path), "supporting_sensors": sensors})
                best_length[key] = max(best_length[key], length)
                observation_records.append(dict(details, id=name, consecutive_days=length, supporting_sensors=support))
            previous[key] = current
            row["variants"][key] = {"rejected": rejected, "openings": observation_records}
        daily.append(row)

    return {
        "summary": {
            "daily_products": len(indexed), "primary_dates": len(daily),
            "detections": {key: sum(len(item["variants"][key]["openings"]) for item in daily) for key in VARIANTS},
            "longest_chain_days": best_length,
            "cross_sensor_three_day_chains": len(chain_support),
        },
        "controls": {"synthetic_detections": synthetic, "historical_maud_rise": history},
        "grid_signature": grid["signature"], "daily": daily,
        "links": all_links, "supported_primary_chains": chain_support,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--download", action="store_true", help="Fetch missing exact bytes from the source manifest")
    args = parser.parse_args()
    entries = json.loads((HERE / "sources.json").read_text(encoding="utf-8"))
    result = analyze(entries, args.data, args.download)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
