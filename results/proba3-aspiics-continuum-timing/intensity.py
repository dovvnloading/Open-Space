"""Post-analysis diagnostic of absolute line-intensity sensitivity."""

from argparse import ArgumentParser
import json
from pathlib import Path

import numpy as np

from measure import ANNULI, WCS, change, flux_per_msb, load, wideband_parents


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    selection = json.loads(Path(__file__).with_name("selection.json").read_text(encoding="utf-8"))
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    axis = (np.arange(2048) + 1 - WCS["CRPIX1"]) * WCS["CDELT1"]
    distance = np.hypot(axis[:, None], axis[None, :])
    result = []
    for orbit, group in selection.items():
        for sample in group["sample"]:
            header, line = load(args.raw, sample["line"], manifest, "photon/s/cm^2/sr")
            earlier, donor = load(args.raw, sample["donor"], manifest, "MSB")
            _, following = load(args.raw, sample["next_brightness"], manifest, "MSB")
            if wideband_parents(header) != wideband_parents(earlier):
                raise ValueError("Wideband lineage differs")
            before = -sample["donor_offset_s"]
            weight = before / (before + sample["next_offset_s"])
            delta = change(donor, following, weight) * flux_per_msb(header)
            keep = (np.isfinite(line) & np.isfinite(donor) & np.isfinite(following)
                    & (line > 0) & (donor > 0) & (following > 0))
            radius = distance / header["RSUN_OBS"]
            for lower, upper in ANNULI:
                mask = keep & (lower <= radius) & (radius < upper)
                ratio = np.abs(delta[mask]) / line[mask]
                result.append({
                    "orbit": int(orbit), "line": sample["line"]["name"], "radius_Rsun": [lower, upper],
                    "positive_pixel_count": len(ratio),
                    "median_ratio": float(np.median(ratio)) if len(ratio) else None,
                    "p95_ratio": float(np.quantile(ratio, 0.95)) if len(ratio) else None,
                    "count_at_least_10_percent": int((ratio >= 0.10).sum()),
                    "count_at_least_50_percent": int((ratio >= 0.50).sum()),
                })
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
