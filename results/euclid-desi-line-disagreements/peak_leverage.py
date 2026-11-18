"""Remove observed high-leverage bins from fixed fits as a post-selection diagnostic."""

from argparse import ArgumentParser
import hashlib
import json
from pathlib import Path

import numpy as np

from check_lines import LINES, decode, fixed_fit


CASES = (
    ("2720393417678497186-9998.vot", "euclid_rank0", "[O_II]", (16201.4,)),
    ("2720534960677393578-10554.vot", "euclid_rank0", "H_beta", (17434.2, 17420.8)),
    ("2716462137676529071-3261.vot", "independent_desi", "Pa_beta", ()),
)


def check(case, inputs, source):
    name, hypothesis, line, wavelengths = case
    record = inputs[name]
    path = source / name
    if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
        raise ValueError(f"Spectrum digest mismatch: {name}")
    wave, flux, var, valid, mask = decode(path, "euclid")
    redshift = record["redshifts"][hypothesis]
    centers = {key: rest * (1 + redshift) for key, rest in LINES.items()}
    width = 2.5 * float(np.median(np.diff(wave)))
    original = fixed_fit(wave, flux, var, valid, centers, line, width)
    core = np.flatnonzero((np.abs(wave - centers[line]) <= 0.8 * width) & valid)
    # For the agreement control, remove its brightest eligible line bin.
    indices = [int(np.argmin(np.abs(wave - target))) for target in wavelengths]
    if not indices:
        indices = [int(core[np.argmax(flux[core])])]
    if any(index not in core for index in indices):
        raise ValueError(f"Requested bin lies outside the usable line core: {name}")
    remaining = valid.copy()
    omissions = []
    for index in indices:
        remaining[index] = False
        omissions.append({
            "removed_wavelength": float(wave[index]), "removed_flux": float(flux[index]),
            "removed_error": float(np.sqrt(var[index])), "removed_mask": int(mask[index]),
            "standard_fit": fixed_fit(wave, flux, var, remaining, centers, line, width),
            "stress_fit": fixed_fit(wave, flux, var, remaining, centers, line, width,
                                    stress_core_support=True),
        })
    return {"name": name, "hypothesis": hypothesis, "z": redshift, "line": line,
            "original": original, "cumulative_omissions": omissions}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    records = {item["name"]: item for item in json.loads(args.manifest.read_bytes())}
    result = [check(case, records, args.manifest.parent) for case in CASES]
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
