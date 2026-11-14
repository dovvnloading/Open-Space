"""Measure unmodeled intensity in CHIME's dedispersed data and intersect its HEALPix regions."""

import argparse
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlopen

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent


def inspect_spectrum(path):
    with h5py.File(path, "r") as f:
        good = f["good_freq"][:].astype(bool)
        freqs = f["index_map/freqs"][:][good]
        data = f["data"][:][good].astype(float)
        valid = f["flag"][:][good].astype(bool)
        times = (f["index_map/times"][:] - f.attrs["center_time"]) * 1000
        beam = int(f.attrs["beam_number"])
    off = np.abs(times) > 20
    on = np.abs(times) <= 6
    masked = np.where(valid, data, np.nan)
    baseline = np.nanmean(masked[:, off], axis=1)
    std = np.nanstd(masked[:, off], axis=1, ddof=1)
    keep = np.isfinite(std) & (std > 0)
    whitened = (masked[keep] - baseline[keep, None]) / std[keep, None]
    freq = freqs[keep]
    results = []
    for low, high in [(400, 600), (600, 700), (700, 800.21), (600, 800.21)]:
        in_band = (freq >= low) & (freq < high)
        profile = np.nanmean(whitened[in_band], axis=0)
        fluence = profile[on].sum()
        # Equal-width disjoint windows, all at least 20 ms from the fitted center.
        edges = np.flatnonzero(off)
        width = int(on.sum())
        controls = []
        for segment in np.split(edges, np.flatnonzero(np.diff(edges) != 1) + 1):
            for start in range(0, len(segment) - width + 1, width):
                controls.append(profile[segment[start:start + width]].sum())
        results.append({"band_mhz": [low, high], "channels": int(in_band.sum()),
                        "on_window_sum": float(fluence), "off_window_count": len(controls),
                        "off_window_mean": float(np.mean(controls)), "off_window_sd": float(np.std(controls, ddof=1)),
                        "on_excess_in_off_window_sd": float((fluence-np.mean(controls))/np.std(controls,ddof=1))})
    return {"beam": beam, "on_pulse_samples": int(on.sum()), "bands": results}


def inspect_maps(paths):
    data = []
    for path in paths:
        with h5py.File(path, "r") as handle:
            h = handle["healpix"]
            # The sparse file repeats native HEALPix indices with the same CL.
            # Count each equal-area sky pixel once, including in the individual maps.
            pixel, first, inverse = np.unique(h["ipix"][:], return_index=True, return_inverse=True)
            confidence = h["CL"][:]
            if not np.array_equal(confidence, confidence[first][inverse]):
                raise ValueError("Repeated pixels have conflicting confidence levels")
            data.append((pixel, confidence[first], int(h.attrs["nside"]), h.attrs["ordering"]))
    (p1, cl1, n1, ord1), (p2, cl2, n2, ord2) = data
    if n1 != n2 or ord1 != ord2:
        raise ValueError("Maps require resampling before intersection")
    _, ia, ib = np.intersect1d(p1, p2, return_indices=True)
    area = 4 * np.pi * (180/np.pi)**2 / (12 * n1**2)
    return {str(level): {"first_deg2": float((cl1 < level).sum()*area),
                         "second_deg2": float((cl2 < level).sum()*area),
                         "overlap_deg2": float(((cl1[ia] < level) & (cl2[ib] < level)).sum()*area)}
            for level in (.5,.68,.9,.95)}


def check(directory, fetch):
    sources = json.loads((ROOT / "sources.json").read_text())["archive_products"]
    maps = []
    observations = {}
    for record in sources:
        path = directory / record["filename"]
        if fetch and not path.exists():
            with urlopen(record["url"], timeout=120) as response:
                contents = response.read()
            if hashlib.sha256(contents).hexdigest() != record["sha256"]:
                raise ValueError(f"Unexpected bytes for {path.name}")
            path.write_bytes(contents)
        if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError(f"Hash mismatch for {path.name}")
        if record["kind"] == "localization":
            maps.append(path)
        else:
            observations[record["burst"]] = inspect_spectrum(path)
    return {"intensity": observations, "localization_regions": inspect_maps(maps)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, help="Reuse locally downloaded original files")
    args = parser.parse_args()
    if args.input_dir:
        result = check(args.input_dir, False)
    else:
        with TemporaryDirectory() as tmp:
            result = check(Path(tmp), True)
    print(json.dumps(result, indent=2))
