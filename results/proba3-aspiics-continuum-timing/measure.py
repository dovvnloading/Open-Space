"""Evaluate a bracketing continuum substitution in frozen ASPIICS L3 FITS."""

from argparse import ArgumentParser
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re

from astropy.io import fits
import numpy as np
from scipy.ndimage import shift


ANNULI = ((1.2, 1.4), (1.4, 1.7))
JITTER = ((-0.5, 0), (0.5, 0), (0, -0.5), (0, 0.5))
WCS = {
    "CTYPE1": "HPLN-TAN", "CTYPE2": "HPLT-TAN",
    "CUNIT1": "arcsec", "CUNIT2": "arcsec",
    "PC1_1": 1, "PC1_2": 0, "PC2_1": 0, "PC2_2": 1,
    "CRVAL1": 0, "CRVAL2": 0, "CRPIX1": 1024.5, "CRPIX2": 1024.5,
    "CDELT1": 2.8172, "CDELT2": 2.8172,
}


def load(raw, product, manifest, unit):
    path = raw / product["name"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != manifest[path.name]["sha256"]:
        raise ValueError(f"Input digest mismatch for {path.name}")
    with fits.open(path, memmap=False, checksum=True) as hdus:
        hdu = hdus[0]
        header = hdu.header
        if (header["QUALITY"] != product["quality"]
                or datetime.fromisoformat(header["DATE-AVG"]) != datetime.fromisoformat(product["date_avg"])
                or header["CYCLE_ID"] != product["cycle"] or header["BUNIT"] != unit
                or header["N_PARENT"] != product["n_parent"]):
            raise ValueError(f"Metadata disagreement for {path.name}")
        if any(header.get(name) != value for name, value in WCS.items()):
            raise ValueError(f"Unexpected delivered WCS in {path.name}")
        array = hdu.data.astype(np.float64)
        # The v03 product documents infinite saturation markers.
        array[~np.isfinite(array)] = np.nan
    if array.shape != (2048, 2048):
        raise ValueError(f"Unexpected array shape in {path.name}")
    return header, array


def wideband_parents(header):
    cards = (header[f"PARENT{i}"] for i in range(1, header["N_PARENT"] + 1))
    return {value for value in cards if value.startswith("aspiics_wb_l2_")}


def flux_per_msb(header):
    pattern = r"MSB equals ([0-9.]+e[+-]?\d+) photon s-1 cm-2 sr-1"
    found = [
        re.fullmatch(pattern, str(history).strip())
        for history in header["HISTORY"]
    ]
    conversions = {float(match[1]) for match in found if match}
    if len(conversions) != 1:
        raise ValueError("Expected a unique Fe XIV FITS conversion from MSB")
    return conversions.pop()


def change(donor, following, weight):
    return weight * (donor - following)


def metrics(region, k_mask, positive_line, perturbation, flips, robust):
    continuum = region & k_mask
    line = region & positive_line
    coverage = int(continuum.sum())
    count = int(line.sum())
    fractions = np.abs(perturbation[continuum])
    flip_count = int((line & flips).sum())
    robust_count = int((line & robust).sum())
    return {
        "geometric_pixels": int(region.sum()),
        "continuum_valid_pixels": coverage,
        "line_positive_pixels": count,
        "abs_fractional_K_median": float(np.median(fractions)) if coverage else None,
        "abs_fractional_K_p95": float(np.quantile(fractions, .95)) if coverage else None,
        "flips": flip_count,
        "robust_flips": robust_count,
        "flip_fraction": flip_count / count if count else None,
        "robust_fraction": robust_count / count if count else None,
    }


def measure(raw, example, manifest):
    line_header, line = load(raw, example["line"], manifest, "photon/s/cm^2/sr")
    donor_header, donor = load(raw, example["donor"], manifest, "MSB")
    next_header, following = load(raw, example["next_brightness"], manifest, "MSB")
    donor_parents = wideband_parents(donor_header)
    if not donor_parents or donor_parents != wideband_parents(line_header):
        raise ValueError(f"Wideband lineage differs for {example['line']['name']}")
    if any(h["ORBIT_ID"] != line_header["ORBIT_ID"] for h in (donor_header, next_header)):
        raise ValueError("Brightness images must come from the same orbit as the line")

    elapsed_before = -example["donor_offset_s"]
    weight = elapsed_before / (elapsed_before + example["next_offset_s"])
    scale = flux_per_msb(line_header)
    signed = change(donor, following, weight)
    k_mask = (donor > 0) & (following > 0) & np.isfinite(donor) & np.isfinite(following)
    positive_line = k_mask & np.isfinite(line) & (line > 0)
    fraction = np.divide(signed, donor, out=np.full_like(signed, np.nan), where=k_mask)
    flips = positive_line & (line + scale * signed <= 0)
    robust = flips.copy()
    for offset in JITTER:
        moved = shift(following, offset, order=1, mode="constant", cval=np.nan, prefilter=False)
        alternative = change(donor, moved, weight)
        robust &= (moved > 0) & np.isfinite(moved) & (line + scale * alternative <= 0)

    axis = (np.arange(2048) + 1 - WCS["CRPIX1"]) * WCS["CDELT1"]
    radius = np.hypot(axis[:, None], axis[None, :]) / line_header["RSUN_OBS"]
    angle = np.degrees(np.arctan2(axis[:, None], axis[None, :])) % 360
    annuli = []
    for lower, upper in ANNULI:
        region = (lower <= radius) & (radius < upper)
        annuli.append({
            "radius_Rsun": [lower, upper],
            **metrics(region, k_mask, positive_line, fraction, flips, robust),
            "sectors": [
                {"start_deg": start,
                 **metrics(region & (start <= angle) & (angle < start + 30),
                           k_mask, positive_line, fraction, flips, robust)}
                for start in range(0, 360, 30)
            ],
        })
    return {
        "orbit": line_header["ORBIT_ID"], "line": example["line"]["name"],
        "date_avg": example["line"]["date_avg"], "quality": line_header["QUALITY"],
        "donor": example["donor"]["name"], "following": example["next_brightness"]["name"],
        "wideband_L2_parents": sorted(donor_parents),
        "donor_offset_s": example["donor_offset_s"], "next_offset_s": example["next_offset_s"],
        "weight": weight, "photons_per_msb": scale,
        "solar_radius_arcsec": [h["RSUN_OBS"] for h in (line_header, donor_header, next_header)],
        "annuli": annuli,
    }


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    selection = json.loads(Path(__file__).with_name("selection.json").read_text(encoding="utf-8"))
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    result = [
        measure(args.raw, sample, manifest)
        for group in selection.values()
        for sample in group["sample"]
    ]
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
