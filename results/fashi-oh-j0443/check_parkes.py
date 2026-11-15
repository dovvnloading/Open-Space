"""Measure two fixed line frequencies in an archived Parkes HIPASS spectrum."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
FREQUENCIES = (1407.29410, 1405.64218)


def amplitude(freq, flux, center, sigma):
    indices = abs(freq - center) < 1.0
    offsets = freq[indices] - center
    design = np.column_stack([np.ones(len(offsets)), offsets,
                              np.exp(-0.5 * (offsets / sigma)**2)])
    return np.linalg.lstsq(design, flux[indices], rcond=None)[0][2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=ROOT / "hipass_spectrum.txt")
    parser.add_argument("--centers", nargs=2, type=float, default=FREQUENCIES, metavar=("MAIN_MHZ", "COMPANION_MHZ"))
    parser.add_argument("--width-kms", type=float, default=44.91)
    args = parser.parse_args()
    path = args.spectrum
    centers = args.centers
    fwhm = centers[0] * args.width_kms / 299792.458
    profile_sigma = fwhm / np.sqrt(8 * np.log(2))
    tokens = [line.split() for line in path.read_text().splitlines()]
    data = np.array([[float(t) for t in row] for row in tokens
                     if len(row) == 3 and row[0][:1].isdigit()])
    _, freq, flux = data.T
    flux *= 1000  # Jy/beam to mJy/beam
    # Omit both test locations and their fit windows, and both band edges.
    controls = freq[(freq < max(freq) - 1) & (freq > min(freq) + 1)]
    for center in centers:
        controls = controls[abs(controls - center) > 1.0]
    vals = np.array([amplitude(freq, flux, center, profile_sigma) for center in controls])
    sigma = 1.4826 * np.median(abs(vals - np.median(vals)))
    out = {"spectrum_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
           "channels": len(freq), "fit_fwhm_mhz": fwhm,
           "comparison_centers": len(controls), "off_line_median_mjy": np.median(vals),
           "off_line_mad_sigma_mjy": sigma, "lines": []}
    for center in centers:
        amp = amplitude(freq, flux, center, profile_sigma)
        out["lines"].append({"center_mhz": center, "fit_peak_mjy": amp,
                             "empirical_signal_to_noise": (amp - np.median(vals)) / sigma})
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
