"""Verify the native THEMIS VIS-EQR files and aggregate common map cells."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np

ROOT = Path(__file__).resolve().parent
BANDS = (0.42, 0.75)
SCALES = (360, 720, 1440)
PAIRS = {"primary": ("VA5997004", "VA6883007"), "control": ("VA5423003", "VA6596006")}


def xml_value(tree, tag):
    return next(element.text.strip() for element in tree.iter()
                if element.tag.rsplit("}", 1)[-1] == tag)


def checked(path, expected_size, expected_sha):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    if path.stat().st_size != expected_size or digest.hexdigest() != expected_sha:
        raise ValueError(f"Unexpected bytes at {path}")


def pvl_location(header, keyword):
    return float(re.search(rf"^\s*{keyword}\s*=\s*([-0-9.e+]+)", header, re.M).group(1))


def read_cube(directory, record):
    cube_path, label_path = directory / record["image"], directory / record["label"]
    checked(cube_path, record["bytes"], record["sha256"])
    checked(label_path, record["label_bytes"], record["label_sha256"])
    label = ET.parse(label_path)
    array = next(node for node in label.iter() if node.tag.endswith("}Array_3D_Spectrum"))
    dims = tuple(int(node.text) for node in array.iter() if node.tag.endswith("}elements"))
    offset = int(xml_value(array, "offset"))
    if (xml_value(array, "axis_index_order") != "Last Index Fastest"
            or xml_value(array, "data_type") != "IEEE754LSBSingle"
            or xml_value(array, "scaling_factor") != "1.0"
            or xml_value(array, "value_offset") != "0.0"):
        raise ValueError("Unexpected EQR recording convention")
    with cube_path.open("rb") as raw:
        header = raw.read(offset).split(b"\0", 1)[0].decode("ascii")
    if not re.search(r"^\s*Format\s*=\s*BandSequential\b", header, re.M):
        raise ValueError("Expected a band-sequential EQR cube")
    wavelengths = [float(node.text) for node in label.iter()
                   if node.tag.endswith("}center_wavelength")]
    band_group = re.search(r"\bGroup\s*=\s*BandBin\b.*?\bEnd_Group", header, re.S).group(0)
    declared = re.search(r"^\s*Center\s*=\s*\(([^)]+)\)", band_group, re.M).group(1)
    if wavelengths != [float(value) for value in declared.split(",")] or len(wavelengths) != dims[0]:
        raise ValueError("Wavelength list disagrees with embedded BandBin")
    geog = {key: float(xml_value(label, key)) for key in (
        "a_axis_radius", "latitude_of_projection_origin", "longitude_of_central_meridian",
        "standard_parallel_1", "upperleft_corner_x", "upperleft_corner_y", "pixel_resolution_x",
        "pixel_resolution_y", "minimum_latitude", "maximum_latitude", "minimum_longitude",
        "maximum_longitude", "incidence_angle", "emission_angle", "subsolar_azimuth",
    )}
    if (xml_value(label, "map_projection_name") != "Equirectangular"
            or xml_value(label, "latitude_type") != "Planetocentric"
            or xml_value(label, "longitude_direction") != "Positive East"
            or geog["latitude_of_projection_origin"] != 0):
        raise ValueError("Unexpected geographic projection")
    for header_key, xml_key, tolerance in (
        ("UpperLeftCornerX", "upperleft_corner_x", 8), ("UpperLeftCornerY", "upperleft_corner_y", 1),
        ("CenterLongitude", "longitude_of_central_meridian", 0.001),
        ("CenterLatitude", "standard_parallel_1", 0.001),
    ):
        if abs(pvl_location(header, header_key) - geog[xml_key]) > tolerance:
            raise ValueError(f"Header and XML disagree for {header_key}")
    radius = geog["a_axis_radius"] * 1000
    east = geog["upperleft_corner_x"] + (0.5 + np.arange(dims[2])) * geog["pixel_resolution_x"]
    north = geog["upperleft_corner_y"] - (0.5 + np.arange(dims[1])) * geog["pixel_resolution_y"]
    lat = np.rad2deg(north / radius)
    lon = geog["longitude_of_central_meridian"] + np.rad2deg(
        east / (radius * np.cos(np.deg2rad(geog["standard_parallel_1"])))
    )
    bounds_error = np.max(np.abs([
        min(lat) - geog["minimum_latitude"], max(lat) - geog["maximum_latitude"],
        min(lon) - geog["minimum_longitude"], max(lon) - geog["maximum_longitude"],
    ]))
    if bounds_error > 0.025:
        raise ValueError("Analytical projection cannot reproduce geographic label limits")
    geog["rectangle_to_label_max_abs_difference_deg"] = float(bounds_error)
    geog["bands_um"] = wavelengths
    profile = np.memmap(cube_path, dtype="<f4", offset=offset, shape=dims, mode="r")
    return {"east": east, "north": north, "geog": geog, "bands": [profile[wavelengths.index(b)] for b in BANDS]}


def aggregate(cube, x_edges, y_edges, scale):
    x = np.floor((cube["east"] - x_edges[0]) / scale).astype(int)
    y = np.floor((cube["north"] - y_edges[0]) / scale).astype(int)
    rows, cols = len(y_edges) - 1, len(x_edges) - 1
    valid_cols = (0 <= x) & (x < cols)
    valid_rows = (0 <= y) & (y < rows)
    keys = (y[valid_rows, None] * cols + x[None, valid_cols]).ravel()
    minimum = np.ceil(0.70 * scale**2 / (cube["geog"]["pixel_resolution_x"]
                                          * cube["geog"]["pixel_resolution_y"]))
    mapped = []
    for source in cube["bands"]:
        values = source[np.ix_(valid_rows, valid_cols)].ravel()
        usable = np.isfinite(values) & (values > 0)
        count = np.bincount(keys[usable], minlength=rows * cols)
        total = np.bincount(keys[usable], weights=values[usable], minlength=rows * cols)
        means = np.full(rows * cols, np.nan)
        np.divide(total, count, out=means, where=count >= minimum)
        mapped.append(means.reshape(rows, cols))
    return np.stack(mapped)


def crop_as_pair(earlier, later, scale):
    radius = later["geog"]["a_axis_radius"] * 1000
    keys = ("a_axis_radius", "standard_parallel_1", "longitude_of_central_meridian")
    if any(earlier["geog"][key] != later["geog"][key] for key in keys):
        raise ValueError("The scenes must use a common map projection")
    xmin = max(min(earlier["east"]), min(later["east"]))
    xmax = min(max(earlier["east"]), max(later["east"]))
    ymin = np.deg2rad(-10.65) * radius
    ymax = np.deg2rad(-7.85) * radius
    x_edges = np.arange(np.floor(xmin / scale) - 2, np.ceil(xmax / scale) + 3) * scale
    y_edges = np.arange(np.floor(ymin / scale), np.ceil(ymax / scale) + 1) * scale
    return {
        "east_m": x_edges[:-1] + scale / 2,
        "north_m": y_edges[:-1] + scale / 2,
        "lat_deg": np.rad2deg((y_edges[:-1] + scale / 2) / radius),
        "earlier": aggregate(earlier, x_edges, y_edges, scale),
        "later": aggregate(later, x_edges, y_edges, scale),
    }


def freeze_inputs(data_root, bundle, receipt_path):
    receipts = json.loads((ROOT / "sources.json").read_text(encoding="utf-8"))
    raw = {row["image"][:9]: read_cube(data_root, row) for row in receipts}
    compact = {}
    for name, (early, late) in PAIRS.items():
        for scale in SCALES:
            for key, array in crop_as_pair(raw[early], raw[late], scale).items():
                compact[f"{name}_{scale}_{key}"] = array
    np.savez_compressed(bundle, **compact)
    evidence = {"bands_um": BANDS, "pairs": PAIRS,
                "products": {name: cube["geog"] for name, cube in raw.items()}}
    receipt_path.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--inputs-out", type=Path, required=True)
    parser.add_argument("--geometry-out", type=Path, required=True)
    args = parser.parse_args()
    freeze_inputs(args.data_root, args.inputs_out, args.geometry_out)
