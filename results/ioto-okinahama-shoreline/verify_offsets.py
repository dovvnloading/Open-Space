"""Replay the original Sentinel-2 versus Earth Search COG encoding check."""

import argparse
from collections import Counter
from hashlib import file_digest
import json
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import rasterio


HERE = Path(__file__).resolve().parent


def sha256(path):
    with path.open("rb") as stream:
        return file_digest(stream, "sha256").hexdigest()


def replay(root):
    manifest = json.loads((HERE / "native-sources.json").read_text(encoding="utf-8"))
    results = []
    for item in manifest["comparisons"]:
        native = root / item["native_crop"]
        converted = root / item["cog_crop"]
        meta = root / item["product_xml"]
        for path, key in ((native, "native_crop_sha256"), (converted, "cog_crop_sha256"), (meta, "product_xml_sha256")):
            if sha256(path) != item[key]:
                raise ValueError(f"SHA-256 mismatch for {path}")
        xml = ElementTree.parse(meta)
        offset_tag = [int(elem.text) for elem in xml.iter() if elem.tag.endswith("BOA_ADD_OFFSET")]
        quant_tag = [int(elem.text) for elem in xml.iter() if elem.tag.endswith("BOA_QUANTIFICATION_VALUE")]
        if set(offset_tag) != {-1000} or quant_tag != [10000]:
            raise ValueError(f"Unexpected per-product native calibration in {meta}")
        with rasterio.open(native) as original, rasterio.open(converted) as cog:
            if original.crs != cog.crs or original.transform != cog.transform or original.shape != cog.shape:
                raise ValueError(f"Nonmatching comparison grids for {cog}")
            a, b = original.read(1).astype(np.int32), cog.read(1).astype(np.int32)
        common = (a != 0) & (b != 0)
        changes, counts = np.unique(a[common] - b[common], return_counts=True)
        differences = {str(int(key)): int(value) for key, value in zip(changes, counts)}
        results.append({
            "scene": item["scene_id"], "band": item["band"], "pixels": int(a.size), "nonzero_pairs": int(common.sum()),
            "source_nonzero_cog_zero": int(((a != 0) & (b == 0)).sum()),
            "difference_counts_dn": differences,
        })
    by_difference = Counter()
    for row in results:
        by_difference.update(row["difference_counts_dn"])
    return {
        "summary": {
            "scene_band_crops": len(results), "scene_count": len({row["scene"] for row in results}),
            "pixels_with_unique_difference": dict(by_difference),
            "total_paired_pixels": sum(row["nonzero_pairs"] for row in results),
            "source_nonzero_cog_zero": sum(row["source_nonzero_cog_zero"] for row in results),
            "native_quantification": 10000, "native_offset_DN": -1000,
            "extra_STAC_reflectance_subtraction": 0.1,
        },
        "comparisons": results,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("extracted_bundle", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = replay(args.extracted_bundle)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(result["summary"]))


if __name__ == "__main__":
    main()
