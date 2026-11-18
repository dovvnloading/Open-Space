"""Measure exposed volcanic deposits and test contact in frozen optical crops."""

import argparse
from dataclasses import dataclass
from hashlib import file_digest
import json
from pathlib import Path

from affine import Affine
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from scipy import ndimage


SCENES = (
    "s2c-20261022", "l9-20261026", "s2b-20261027", "s2c-20261101",
    "l8-20261103", "s2b-20261106", "l9-20261111", "s2c-20261111", "s2b-20261116",
)
BASES = SCENES[:2]
ALT_BASES = ("s2b-20261017", "l8-20261018")
PAIR = ("l9-20261111", "s2c-20261111")
VENT = (533869.3452570109, 2738783.8378536175)
BOUNDS = (532440, 2737320, 535320, 2740200)
CROSS = ndimage.generate_binary_structure(2, 1)
REFLECTANCE = ("blue", "green", "red", "nir", "swir1", "swir2")
THRESHOLDS = {"nominal": (0.035, 0.055, 0.020, 0.18), "relaxed": (0.025, 0.045, 0.014, 0.28),
              "strict": (0.048, 0.067, 0.030, 0.10)}
CRS = "EPSG:32654"


@dataclass
class Grid:
    resolution: int

    def __post_init__(self):
        x0, y0, x1, y1 = BOUNDS
        self.shape = ((y1-y0) // self.resolution, (x1-x0) // self.resolution)
        self.transform = Affine(self.resolution, 0, x0, 0, -self.resolution, y1)
        i, j = np.indices(self.shape)
        self.x = x0 + (j + 0.5) * self.resolution
        self.y = y1 - (i + 0.5) * self.resolution
        self.r = np.hypot(self.x - VENT[0], self.y - VENT[1])
        self.region = self.r <= 1100

    def project(self, data, transform, src_nodata=None, dst_nodata=None, method=Resampling.average):
        result = np.zeros(self.shape, dtype=np.float32)
        reproject(data, result, src_transform=transform, dst_transform=self.transform, src_crs=CRS,
                  dst_crs=CRS, src_nodata=src_nodata, dst_nodata=dst_nodata, resampling=method)
        return result


@dataclass
class Frame:
    values: dict
    valid: np.ndarray
    qa_valid: np.ndarray
    timestamp: str
    scene_id: str


def load(root, scene, grid, sentinel_offset):
    directory = root / scene
    receipt = json.loads((directory / "receipt.json").read_text(encoding="utf-8-sig"))
    assets = {item["band"]: item for item in receipt["assets"]}
    sentinel = scene.startswith("s2")
    keys = (*REFLECTANCE, "qa") if sentinel else (*REFLECTANCE, "qa", "qa_radsat")
    values = {}
    support = np.ones(grid.shape, dtype=bool)
    qa_valid = np.ones(grid.shape, dtype=bool)
    for key in keys:
        asset = assets[key]
        filename = directory / asset["output"]
        with filename.open("rb") as stream:
            if file_digest(stream, "sha256").hexdigest() != asset["output_sha256"]:
                raise ValueError(f"Input SHA-256 mismatch for {filename}")
        with rasterio.open(filename) as dataset:
            if str(dataset.crs) != CRS:
                raise ValueError(f"Unexpected source CRS for {filename}")
            raw = dataset.read(1)
            if key == "qa":
                good = ~np.isin(raw, (0, 1, 3, 8, 9, 10, 11)) if sentinel else ((raw & 63) == 0)
                qa_valid &= grid.project(good.astype(np.float32), dataset.transform) >= 0.9999
            elif key == "qa_radsat":
                qa_valid &= grid.project((raw == 0).astype(np.float32), dataset.transform) >= 0.9999
            else:
                nodata = dataset.nodata if dataset.nodata is not None else 0
                native_valid = (raw != nodata) & (raw != 0)
                offset = 0 if sentinel and sentinel_offset == "zero" else float(asset["offset"])
                reflectance = raw.astype(np.float32) * float(asset["scale"]) + offset
                reflectance[~native_valid] = np.nan
                values[key] = grid.project(reflectance, dataset.transform, np.nan, np.nan)
                if key in ("green", "nir", "swir1", "swir2"):
                    support &= grid.project(native_valid.astype(np.float32), dataset.transform) >= 0.9999
    return Frame(values, support & qa_valid, qa_valid, receipt["scene_datetime"], receipt["scene_id"])


def classes(frame, thresholds=THRESHOLDS["nominal"], expand_cloud=False):
    green, swir1, nir, swir2 = (frame.values[key] for key in ("green", "swir1", "nir", "swir2"))
    g = np.clip(green, 0, None)
    s = np.clip(swir1, 0, None)
    ndwi = np.divide(g - s, g + s, out=np.zeros_like(g), where=(g + s) > 0)
    valid = frame.valid & np.isfinite(swir1 + nir + swir2 + green)
    if expand_cloud:
        valid &= ndimage.binary_erosion(frame.qa_valid, structure=np.ones((3, 3), dtype=bool), border_value=1)
    min_swir, min_nir, min_swir2, max_ndwi = thresholds
    dry = valid & (swir1 >= min_swir) & (nir >= min_nir) & (swir2 >= min_swir2) & (ndwi <= max_ndwi)
    water = valid & (swir1 <= 0.020) & ((nir <= 0.055) | (ndwi >= 0.30))
    return valid, dry, water


def components(gain, near_vent, minimum):
    labels, n = ndimage.label(gain, structure=CROSS)
    sizes = np.bincount(labels[gain], minlength=n + 1)
    qualifying = np.flatnonzero(sizes >= minimum)
    close = np.unique(labels[near_vent & gain])
    found = np.intersect1d(qualifying, close)
    keep = np.isin(labels, found)
    return keep, {"components_near_vent": int(len(found)),
                  "distant_component_cells": int(sizes[np.setdiff1d(qualifying, found)].sum()),
                  "small_component_cells": int(sizes[sizes < minimum].sum())}


def connected(dry, ground, inland, erosion):
    mask = ndimage.binary_erosion(dry, structure=CROSS, iterations=erosion) if erosion else dry
    labels, _ = ndimage.label(mask, structure=CROSS)
    seen = np.unique(labels[ground & mask])
    old = np.unique(labels[inland & mask])
    return bool(np.intersect1d(seen, old).size)


def synthetic():
    dry = np.zeros((20, 30), dtype=bool)
    dry[3:17, 2:10] = True
    dry[3:17, 19:27] = True
    ground, inland = dry.copy(), dry.copy()
    ground[:, :15] = False
    inland[:, 15:] = False
    dry[8:13, 10:19] = True
    robust = connected(dry, ground, inland, 1)
    dry[8:13, 10:19] = False
    dry[10, 10:19] = True
    thin_raw = connected(dry, ground, inland, 0)
    thin_eroded = connected(dry, ground, inland, 1)
    dry[10, 14] = False
    broken = connected(dry, ground, inland, 0)
    if not (robust and thin_raw and not thin_eroded and not broken):
        raise RuntimeError("Synthetic controls failed")
    return {"robust_neck": robust, "one_pixel_path": thin_raw,
            "eroded_one_pixel_path": thin_eroded, "one_pixel_gap": broken}


def deep_water(frame, grid):
    sector = (abs(grid.x - VENT[0]) <= 500) & (grid.y >= 2737390) & (grid.y <= 2737650)
    usable = frame.valid & sector
    result = {"valid_cells": int(usable.sum()), "expected_cells": int(sector.sum())}
    for key in ("green", "nir", "swir1", "swir2"):
        samples = frame.values[key][usable]
        result[key] = {"median": float(np.median(samples)), "p10": float(np.quantile(samples, 0.10)),
                       "p90": float(np.quantile(samples, 0.90))} if samples.size else None
    swir = result["swir1"]
    result["passes_plausibility"] = bool(swir and usable.sum() >= 20 and -0.045 <= swir["median"] <= 0.035)
    return result


def register(baseline, frame, grid, inland):
    ref = baseline.values["swir1"]
    post = frame.values["swir1"]
    area = inland & (grid.y > VENT[1] + 700)
    y, x = np.indices(grid.shape)
    training = ((x // 5 + y // 5) % 2) == 0

    def correlation(shifted, available, split):
        use = area & baseline.valid & available & split & np.isfinite(ref) & np.isfinite(shifted)
        if use.sum() < 35:
            return None
        a, b = ref[use], shifted[use]
        if np.std(a) < 1e-6 or np.std(b) < 1e-6:
            return None
        return round(float(np.corrcoef(a, b)[0, 1]), 6)

    trials = []
    for row in range(-2, 3):
        for col in range(-2, 3):
            delta = (row, col)
            shifted = ndimage.shift(post, delta, order=0, cval=np.nan)
            available = ndimage.shift(frame.valid, delta, order=0, cval=False)
            score = correlation(shifted, available, training)
            if score is not None:
                trials.append((score, -(abs(row) + abs(col)), row, col, correlation(shifted, available, ~training)))
    if not trials:
        return {"available": False}
    best = max(trials)
    zero = next((t for t in trials if t[2:4] == (0, 0)), None)
    return {"available": True, "best_shift_row_col": list(best[2:4]), "fit_corr": best[0],
            "heldout_corr": best[4], "normal_shift_heldout_corr": zero[4] if zero else None}


def controls(frame, base_classes, grid):
    valid, dry, water = frame
    base_land, base_water = base_classes
    away_land = ~grid.region & base_land & (ndimage.distance_transform_edt(base_land) * grid.resolution >= 150)
    away_water = ~grid.region & base_water & (ndimage.distance_transform_edt(base_water) * grid.resolution >= 150)
    return {"stable_land_eligible": int((away_land & valid).sum()),
            "stable_land_retained": int((away_land & dry).sum()),
            "open_water_eligible": int((away_water & valid).sum()),
            "false_dry_on_open_water": int((away_water & dry).sum())}


def measure(frame, bases, grid, thresholds, coast_buffer, expand, vent_radius):
    base_a = classes(bases[0], thresholds)
    base_b = classes(bases[1], thresholds)
    old_land = base_a[1] & base_b[1]
    old_water = base_a[2] & base_b[2]
    inland = old_land & (ndimage.distance_transform_edt(old_land) * grid.resolution >= 100)
    valid, dry, water = classes(frame, thresholds, expand)
    separation = ndimage.distance_transform_edt(~old_land) * grid.resolution
    delta = frame.values["swir1"] - np.maximum(bases[0].values["swir1"], bases[1].values["swir1"])
    gain = old_water & dry & grid.region & (separation > coast_buffer) & (delta >= 0.02)
    new, extra = components(gain, grid.r <= vent_radius, int(np.ceil(6000 / grid.resolution ** 2)))
    close = grid.r <= 450
    nominal_contact = connected(dry, new, inland, 0)
    robust_contact = connected(dry, new, inland, 1)
    outside = controls((valid, dry, water), (old_land, old_water), grid)
    return {
        "new_cells": int(new.sum()), "optical_area_km2": round(float(new.sum() * grid.resolution ** 2 / 1e6), 6),
        "region_cells": int(grid.region.sum()), "region_valid": int((grid.region & valid).sum()),
        "vent_cells": int(close.sum()), "vent_valid": int((close & valid).sum()),
        "vent_water_cells": int((close & water).sum()), "vent_dry_cells": int((close & dry).sum()),
        "nearest_dry_m": round(float(grid.r[dry].min()), 3) if dry.any() else None,
        "baseline_water_cells": int((old_water & grid.region).sum()),
        "nominal_contact": nominal_contact, "eroded_contact": robust_contact,
        "controls": outside, **extra,
    }, {"new": new, "old_land": old_land, "old_water": old_water, "valid": valid, "dry": dry, "water": water, "inland": inland}


def compare_pair(pairs, common):
    new, known = {}, {}
    for name, (masks, grid) in pairs.items():
        valid_fraction = common.project(masks["valid"].astype(np.float32), grid.transform)
        fraction = common.project(masks["new"].astype(np.float32), grid.transform)
        known[name] = valid_fraction >= 2 / 3 - 1e-6
        new[name] = (fraction >= 2 / 3 - 1e-6) & known[name]
    both = known[PAIR[0]] & known[PAIR[1]] & common.region
    both_land = new[PAIR[0]] & new[PAIR[1]] & both
    either_land = (new[PAIR[0]] | new[PAIR[1]]) & both
    union = int(either_land.sum())
    overlap = int(both_land.sum())
    return {"both_valid_comparison_cells": int(both.sum()), "both_new_cells": overlap,
            "either_new_cells": union, "union_km2": union * 0.0036, "agreement_km2": overlap * 0.0036,
            "jaccard": round(overlap / union, 6) if union else None}


def run(root, sentinel_offset="provider", baseline_names=BASES, vent_radius=450):
    checks = synthetic()
    grids = {res: Grid(res) for res in (20, 30, 60)}
    bases = {res: [load(root, name, grids[res], sentinel_offset) for name in baseline_names] for res in (20, 30)}
    baseline_checks = {str(res): {name: deep_water(frame, grids[res]) for name, frame in zip(baseline_names, inputs)}
                       for res, inputs in bases.items()}
    comparisons = {}
    rows = []
    plot_data = {}
    for name in dict.fromkeys((*baseline_names, *SCENES)):
        res = 20 if name.startswith("s2") else 30
        grid = grids[res]
        if not (root / name / "receipt.json").exists():
            rows.append({"id": name, "status": "source unavailable"})
            continue
        frame = load(root, name, grid, sentinel_offset)
        spl = deep_water(frame, grid)
        values, masks = measure(frame, bases[res], grid, THRESHOLDS["nominal"], 40, False, vent_radius)
        sensitivity = {
            case: measure(frame, bases[res], grid, thresh, buffer, expand, vent_radius)[0]
            for case, thresh, buffer, expand in (
                ("relaxed", THRESHOLDS["relaxed"], 40, False),
                ("strict", THRESHOLDS["strict"], 40, False),
                ("coast20", THRESHOLDS["nominal"], 20, False),
                ("coast60", THRESHOLDS["nominal"], 60, False),
                ("cloud_buffer", THRESHOLDS["nominal"], 40, True),
            )
        }
        alignment = register(bases[res][0 if res == 20 else 1], frame, grid, masks["inland"])
        valid_baselines = all(item["passes_plausibility"] for item in baseline_checks[str(res)].values())
        evidence = "usable" if valid_baselines and spl["passes_plausibility"] else "limited by deep-water control"
        rows.append(dict(values, id=name, source=frame.scene_id, utc=frame.timestamp, analysis_pixel_m=res,
                         status=evidence, deep_water=spl, registration=alignment, sensitivity=sensitivity))
        plot_data[name] = (frame, grid, masks, evidence)
        if name in PAIR and evidence == "usable":
            comparisons[name] = masks, grid
    pair = compare_pair(comparisons, grids[60]) if all(name in comparisons for name in PAIR) else None
    return {"sentinel_offset": sentinel_offset, "vent_radius": vent_radius, "baseline_names": list(baseline_names),
            "controls": checks, "baselines": baseline_checks, "scenes": rows, "nov11_60m_comparison": pair}, plot_data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--sentinel-offset", choices=("provider", "zero"), default="provider")
    parser.add_argument("--alternate-baselines", action="store_true")
    parser.add_argument("--vent-radius", type=int, choices=(450, 800), default=450)
    args = parser.parse_args()
    base = ALT_BASES if args.alternate_baselines else BASES
    findings, _ = run(args.data, args.sentinel_offset, base, args.vent_radius)
    args.output.write_text(json.dumps(findings, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    for row in findings["scenes"]:
        print(row["id"], row["status"], row.get("optical_area_km2"), row.get("eroded_contact"))
    print("60 m check", findings["nov11_60m_comparison"])


if __name__ == "__main__":
    main()
