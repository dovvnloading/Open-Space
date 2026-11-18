"""Fit no spatial model. Compare calibrated VIS-RDR detector-plane summaries."""

import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

from decode import BANDS, checked, xml_value

ROOT = Path(__file__).resolve().parent


def summarize_one(data_root, source):
    raw, label = data_root / source["image"], data_root / source["label"]
    checked(raw, source["bytes"], source["sha256"])
    checked(label, source["label_bytes"], source["label_sha256"])
    tree = ET.parse(label)
    item = next(node for node in tree.iter() if node.tag.endswith("}Array_3D_Spectrum"))
    wavelengths = [float(node.text) for node in tree.iter() if node.tag.endswith("}center_wavelength")]
    dims = tuple(int(node.text) for node in item.iter() if node.tag.endswith("}elements"))
    if (xml_value(item, "data_type") != "SignedMSB2"
            or xml_value(item, "axis_index_order") != "Last Index Fastest"
            or len(wavelengths) != dims[0]):
        raise ValueError("Unexpected RDR array layout")
    header_bytes = int(xml_value(item, "offset"))
    factor = float(xml_value(item, "scaling_factor"))
    addend = float(xml_value(item, "value_offset"))
    minimum = int(xml_value(item, "valid_minimum"))
    data = np.memmap(raw, dtype=">i2", mode="r", offset=header_bytes, shape=dims)
    height, width = dims[1:]
    filtered = data[[wavelengths.index(band) for band in BANDS]]
    median, used = [], []
    for subset in (filtered, filtered[:, height//4:3*height//4, width//4:3*width//4]):
        all_values, counts = [], []
        for image in subset:
            radiance = image[image >= minimum].astype(float) * factor + addend
            valid = radiance[np.isfinite(radiance) & (radiance > 0)]
            all_values.append(float(np.median(valid)))
            counts.append(len(valid))
        median.append(all_values)
        used.append(counts)
    return {
        "file": raw.name,
        "median_all_and_central": median,
        "samples_all_and_central": used,
        "solar_distance_au": float(xml_value(tree, "target_heliocentric_distance")),
        "incidence_deg": float(xml_value(tree, "incidence_angle")),
        "subsolar_azimuth_deg": float(xml_value(tree, "subsolar_azimuth")),
    }


def run(source_root):
    records = json.loads((ROOT / "native-rdr-sources.json").read_text(encoding="utf-8"))
    early, late = [summarize_one(source_root, source) for source in records]
    distance_factor = (early["solar_distance_au"] / late["solar_distance_au"])**2
    cosine_factor = np.cos(np.deg2rad(late["incidence_deg"])) / np.cos(
        np.deg2rad(early["incidence_deg"])
    )
    combined = float(distance_factor * cosine_factor)
    uncorrected = np.asarray(late["median_all_and_central"]) / np.asarray(early["median_all_and_central"])
    return {
        "bands_um": BANDS, "products": [early, late],
        "center_inverse_square_factor": distance_factor,
        "center_lambert_cosine_factor": float(cosine_factor),
        "center_geometry_factor": combined,
        "ratio_all_and_central": uncorrected.tolist(),
        "ratio_after_center_geometry_factor": (uncorrected / combined).tolist(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.dumps(run(args.data_root), indent=2, allow_nan=False) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
