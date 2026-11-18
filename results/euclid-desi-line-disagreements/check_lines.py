"""Check companion lines at fixed redshifts in native Euclid or DESI arrays."""

from argparse import ArgumentParser
import csv
import hashlib
import json
from pathlib import Path
import warnings

from astropy.io.votable import parse_single_table
from astropy.units import UnitsWarning
import numpy as np


LINES = {
    "C_IV": 1549.5, "He_II_UV": 1640.4, "O_III_UV": 1666.2,
    "C_III]": 1908.7, "C_II]": 2326.0, "Mg_II": 2798.8,
    "[O_II]": 3728.5, "H_delta": 4102.9, "H_gamma": 4341.7,
    "H_beta": 4862.72, "[O_III]_4959": 4960.295, "[O_III]_5007": 5008.239,
    "H_alpha": 6564.614, "[N_II]": 6585.27, "[S_II]_6716": 6718.29,
    "[S_II]_6731": 6732.67, "[S_III]_9069": 9071.1, "[S_III]_9531": 9533.2,
    "He_I": 10833.3, "Pa_gamma": 10941.1, "Pa_beta": 12821.6,
}


def fixed_fit(wave, flux, var, valid, centers, line, fwhm):
    target = centers[line]
    relative = (wave - target) / fwhm
    core = np.abs(relative) <= 0.8
    before = (relative >= -4.0) & (relative <= -1.6)
    after = (relative >= 1.6) & (relative <= 4.0)
    clean = valid.copy()
    for neighbor, center in centers.items():
        if neighbor != line and abs(target - center) > fwhm:
            clean &= np.abs(wave - center) > fwhm

    sizes = [int((mask & clean).sum()) for mask in (core, before, after)]
    record = {"wavelength": target, "fwhm": fwhm, "core_grid": int(core.sum()),
              "core_usable": sizes[0], "continuum_usable": sizes[1:],
              "peak": None, "peak_error": None, "formal_sigma": None}
    if target < wave[0] or target > wave[-1]:
        return dict(record, status="outside")
    if sizes[0] < 3 or sizes[0] < 0.75 * core.sum() or min(sizes[1:]) < 4:
        return dict(record, status="inadequate-native-pixels")

    select = (core | before | after) & clean
    x = relative[select]
    design = np.column_stack((np.ones(len(x)), x, np.exp(-4 * np.log(2) * x*x)))
    error = np.sqrt(var[select])
    lhs = design / error[:, None]
    rhs = flux[select] / error
    coef = np.linalg.lstsq(lhs, rhs, rcond=None)[0]
    formal = float(np.sqrt(np.linalg.inv(lhs.T @ lhs)[2, 2]))
    return dict(record, status="fit", peak=float(coef[2]), peak_error=formal,
                formal_sigma=float(coef[2] / formal))


def decode(path, instrument):
    if instrument == "euclid":
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UnitsWarning)
            table = parse_single_table(path).array
        wave, flux, var = [np.asarray(table[key], dtype=float) for key in ("WAVELENGTH", "SIGNAL", "VAR")]
        mask, ndith = [np.asarray(table[key], dtype=int) for key in ("MASK", "NDITH")]
        allowed = (mask >= 0) & (mask < 64) & (mask % 2 == 0) & (ndith > 0)
    else:
        with path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
        wave, flux, inverse = (np.array([float(row[key]) for row in rows])
                               for key in ("wavelength_A", "flux", "ivar"))
        var = np.divide(1, inverse, out=np.full_like(inverse, np.nan), where=inverse > 0)
        mask = np.array([int(row["mask"]) for row in rows])
        allowed = mask == 0
    if len(wave) < 5 or not np.all(np.isfinite(wave)) or np.any(np.diff(wave) <= 0):
        raise ValueError(f"Invalid wavelength grid: {path}")
    valid = allowed & np.isfinite(var) & np.isfinite(flux) & (var > 0)
    return wave, flux, var, valid, mask


def inspect(item, input_root):
    path = input_root / item["name"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
        raise ValueError(f"Spectrum digest mismatch: {path.name}")
    wave, flux, var, valid, mask = decode(path, item["instrument"])
    step = float(np.median(np.diff(wave)))
    result = {
        "name": item["name"], "object": item["object"], "instrument": item["instrument"],
        "wavelength_range": [float(wave[0]), float(wave[-1])],
        "bins": len(wave), "usable_bins": int(valid.sum()), "median_step": step, "hypotheses": {},
    }
    for name, z in item["redshifts"].items():
        centers = {line: rest * (1+z) for line, rest in LINES.items()}
        checks = {}
        for line, center in centers.items():
            width = 2.5 * step if item["instrument"] == "euclid" else max(2.5*step, center*200/299792.458)
            primary = fixed_fit(wave, flux, var, valid, centers, line, width)
            stress = fixed_fit(wave, flux, var, valid, centers, line, 1.4*width)
            nearest = int(np.argmin(np.abs(wave-center)))
            checks[line] = dict(primary, stress=stress, nearest_grid_wavelength=float(wave[nearest]),
                                nearest_grid_mask=int(mask[nearest]), nearest_usable=bool(valid[nearest]))
        result["hypotheses"][name] = {"z": z, "lines": checks}
    return result


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_bytes())
    result = [inspect(item, args.manifest.parent) for item in manifest]
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
