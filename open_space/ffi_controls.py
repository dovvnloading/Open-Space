"""Calibrate blind image localization by injecting synthetic point sources into real full frames."""

import argparse
import json
from pathlib import Path

import numpy as np

from .ffi import detect, reconcile
from .ffi_archive import acquire
from .tess_archive import atomic_bytes
from .tess_eb_catalog import separation


def calibration(field, cache):
    # Coordinates are fixed on the sky; detector pixels can differ between visits.
    ra, dec = field["ra_deg"]+0.02, field["dec_deg"]-0.006
    injections = [(0.411, 1.0), (0.832, 2.5), (1.69, 5.0)]
    cubes = [acquire(field, product, cache)[0] for product in field["products"]]
    result = []
    for period, peak_rate in injections:
        recovered = []
        for time, flux, error, wcs in cubes:
            x,y = wcs.world_to_pixel_values(ra,dec)
            yy,xx = np.indices(flux.shape[1:])
            psf = np.exp(-((xx-x)**2+(yy-y)**2)/(2*0.85**2))
            wave = np.sin(2*np.pi*(time-1325.0)/period)
            injected = flux + wave[:,None,None]*peak_rate*psf
            peaks, _ = detect((time, injected, error, wcs))
            recovered.append(peaks)
        linked, _ = reconcile(*recovered)
        relevant = sorted(linked, key=lambda s: separation(s["ra_deg"],s["dec_deg"],ra,dec))
        match = relevant[0] if relevant else None
        distance = separation(match["ra_deg"],match["dec_deg"],ra,dec) if match else None
        success = bool(match and match["candidate"] and distance < 21 and
                       abs(match["discovery"]["period_days"]/period-1)<0.025)
        result.append({"period_days": period, "peak_semiamplitude_e_per_s": peak_rate, "injected": True,
                       "position_match_arcsec": distance, "recovered": success, "nearest_detection": match})
    # Permute time labels for complete images while retaining each frame's spatial structure.
    null_peaks = []
    for j,(time,flux,error,wcs) in enumerate(cubes):
        perm = np.random.default_rng(91831+j).permutation(len(time))
        peaks, _ = detect((time,flux[perm],error[perm],wcs))
        null_peaks.append(peaks)
    linked, _ = reconcile(*null_peaks)
    return {"schema_version": 1, "field_id": field["field_id"], "source_science_sha256":
            [p["science_sha256"] for p in field["products"]], "ra_deg": ra, "dec_deg": dec,
            "gaussian_psf_sigma_pixels": 0.85, "injections": result, "permuted_image_null":
            {"seeds": [91831,91832], "peaks_per_sector": list(map(len,null_peaks)),
             "recurring_peaks": sum(s["candidate"] for s in linked)}}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan", type=Path, default=Path(__file__).resolve().parents[1]/"data/ffi-fields.json")
    p.add_argument("--field", default="blind")
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    field = next(f for f in json.loads(args.plan.read_bytes())["fields"] if f["field_id"] == args.field)
    atomic_bytes(args.output, (json.dumps(calibration(field,args.cache),indent=2,allow_nan=False)+"\n").encode())


if __name__ == "__main__":
    main()
