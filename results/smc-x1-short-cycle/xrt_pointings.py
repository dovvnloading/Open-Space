"""Extract pulse-invariant energy channels at SMC X-1 from pointed Swift/XRT PC files."""
import argparse
import csv
import hashlib
import gzip
import io
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from energy_bands import MINIMA, ROOT

SKY = (19.27144054, -73.44333744)
INPUTS = ROOT / "data/xrt-pointings.csv"
EVENTS = ROOT / "data/xrt-events.csv"
PI_BINS = (30, 100, 200, 1001)
REGIONS = {"circle": (0, 20), "wing6": (6, 20), "wing10": (10, 20), "background": (30, 60)}


def mission_bytes(url, digest, raw_dir):
    if raw_dir:
        raw = (raw_dir / url.rsplit("/", 1)[1]).read_bytes()
    else:
        with urlopen(url, timeout=60) as response:
            raw = response.read()
    if hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError(f"Mission input digest mismatch: {url}")
    return raw


def exposure_areas(raw):
    with fits.open(io.BytesIO(gzip.decompress(raw))) as hdus:
        image = hdus[0].data
        x, y = WCS(hdus[0].header, naxis=2).all_world2pix(*SKY, 0)
        yy, xx = np.indices(image.shape)
        radius = np.hypot(xx - x, yy - y)
        return {name: float(image[(radius >= low) & (radius < high)].sum(dtype=float))
                for name, (low, high) in REGIONS.items()}


def extract(manifest, output, raw_dir):
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["obsid", "radius_px", "pi", "grade"])
        for row in manifest:
            raw = mission_bytes(row["url"], row["sha256"], raw_dir)
            map_raw = mission_bytes(row["exposure_map_url"], row["exposure_map_sha256"], raw_dir)
            for region, value in exposure_areas(map_raw).items():
                if not np.isclose(value, float(row[region + "_px_s"]), rtol=1e-8, atol=1e-4):
                    raise ValueError("Pinned exposure-area differs from the mission image")
            with fits.open(io.BytesIO(gzip.decompress(raw))) as product:
                hdu = product["EVENTS"]
                if hdu.header["OBS_ID"] != row["obsid"] or hdu.header["DATAMODE"] != "PHOTON":
                    raise ValueError("Expected the listed photon-counting observation")
                if not np.isclose(hdu.header["EXPOSURE"], float(row["exposure_s"]), rtol=0, atol=1e-6):
                    raise ValueError("Manifest exposure differs from FITS exposure")
                gti = product["GTI"].data
                duration = gti["STOP"] - gti["START"]
                mid_met = np.average((gti["START"] + gti["STOP"]) / 2, weights=duration)
                mid_mjd = hdu.header["MJDREFI"] + hdu.header["MJDREFF"] + mid_met / 86400
                if abs(mid_mjd - float(row["mjd_mid"])) > 1e-7:
                    raise ValueError("Pinned midpoint differs from the mission GTIs")
                sky = WCS(hdu.header, keysel=["pixel"], colsel=[2, 3])
                x, y = sky.all_world2pix(*SKY, 1)
                ev = hdu.data
                radii = np.hypot(ev["X"]-x, ev["Y"]-y)
                for radius, pi, grade in zip(radii, ev["PI"], ev["GRADE"]):
                    if radius < 60:
                        writer.writerow([row["obsid"], f"{radius:.9f}", int(pi), int(grade)])


def phase_labels(midpoint):
    dt = midpoint - 60505.7361918
    orbit = float(abs((.25696017 * dt + .5 * 2.379e-9 * dt**2 + .5) % 1 - .5))
    cycle = int(np.searchsorted(MINIMA, midpoint, side="right") - 1)
    tag = "outside_four_cycles"
    if 0 <= cycle < len(MINIMA) - 1:
        phase = (midpoint - MINIMA[cycle]) / (MINIMA[cycle+1] - MINIMA[cycle])
        if abs(phase - .5) <= 1/6:
            tag = "middle"
        elif min(phase, 1 - phase) <= 1/6:
            tag = "edge"
        else:
            tag = "intermediate"
    return orbit, tag


def measure(manifest, event_path):
    with event_path.open() as stream:
        rows = list(csv.DictReader(stream))
    ids = {row["obsid"] for row in manifest}
    if any(row["obsid"] not in ids for row in rows):
        raise ValueError("Unrecognized observation in event subset")
    visits = []
    for row in manifest:
        mid = float(row["mjd_mid"])
        orbit, label = phase_labels(mid)
        selected = [r for r in rows if r["obsid"] == row["obsid"] and 0 <= int(r["grade"]) <= 12]
        pi = np.array([int(r["pi"]) for r in selected])
        radius = np.array([float(r["radius_px"]) for r in selected])
        counts = {k: np.histogram(pi[(radius >= lo) & (radius < hi) & (pi <= 1000)], bins=PI_BINS)[0].tolist()
                  for k, (lo, hi) in REGIONS.items()}
        area = {k: float(row[k + "_px_s"]) for k in REGIONS}
        visits.append(dict(obsid=row["obsid"], exposure_s=float(row["exposure_s"]),
                           orbital_distance=orbit, state=label, counts=counts, exposure_area_px_s=area))
    summary = {}
    for label in ("middle", "edge"):
        good = [v for v in visits if v["orbital_distance"] >= .17 and v["state"] == label]
        exposure = sum(v["exposure_s"] for v in good)
        counts = {k: np.sum([v["counts"][k] for v in good], axis=0) for k in REGIONS}
        band_rates = {}
        poisson_se = {}
        mean_region_exposure = {}
        for region in ("wing6", "wing10"):
            lo, hi = REGIONS[region]
            areas = np.array([v["exposure_area_px_s"][region] for v in good])
            bkg_areas = np.array([v["exposure_area_px_s"]["background"] for v in good])
            source_counts = np.array([v["counts"][region] for v in good])
            bkg_counts = np.array([v["counts"]["background"] for v in good])
            scale = (areas / bkg_areas)[:, None]
            region_seconds = areas.sum() / (np.pi * (hi**2 - lo**2))
            net = (source_counts - scale * bkg_counts).sum(axis=0)
            variance = (source_counts + scale**2 * bkg_counts).sum(axis=0)
            band_rates[region] = np.round(net / region_seconds, 6).tolist()
            poisson_se[region] = np.round(np.sqrt(variance) / region_seconds, 6).tolist()
            mean_region_exposure[region] = round(float(region_seconds), 5)
        summary[label] = dict(obsids=[v["obsid"] for v in good], exposure_s=round(exposure, 5),
                              counts={k: v.tolist() for k, v in counts.items()},
                              net_annular_rates_per_s=band_rates, poisson_se_per_s=poisson_se,
                              aperture_mean_exposure_s=mean_region_exposure)
    return dict(event_subset_sha256=hashlib.sha256(event_path.read_bytes()).hexdigest(),
                pi_bands_keV=["[0.3,1)", "[1,2)", "[2,10]"], outside_eclipse=summary, visits=visits)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--extract", action="store_true", help="Regenerate the compact event subset")
    p.add_argument("--raw-dir", type=Path, help="Read saved mission files when extracting")
    p.add_argument("--events", type=Path, default=EVENTS)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    with INPUTS.open() as f:
        manifest = list(csv.DictReader(f))
    if args.extract:
        extract(manifest, args.events, args.raw_dir)
    output = json.dumps(measure(manifest, args.events), indent=2) + "\n"
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end="")
