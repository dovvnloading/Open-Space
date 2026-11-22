"""Estimate observed aperture energy flux at the detected count-rate peak.

https://cxc.cfa.harvard.edu/ciao/threads/eff2evt/
https://cxc.cfa.harvard.edu/ciao/ahelp/eff2evt.html
"""
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
from astropy.io import fits

from archive import dump

EV_IN_ERG = 1.602176634e-12


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def peak_flux(event_path, trigger, output_dir):
    """Return net observed aperture flux. Exceptions leave diagnostics on disk."""
    event_path = Path(event_path).resolve()
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    x, y, radius = (float(trigger[k]) for k in ("aperture_x", "aperture_y", "aperture_radius_px"))
    ccd = int(trigger["ccd_id"])
    t0 = float(trigger["peak_start_met_s"])
    t1 = t0 + float(trigger["peak_window_s"])
    live = float(trigger["peak_live_s"])
    inner, outer = 2 * radius, 4 * radius
    if not all(map(math.isfinite, (x, y, radius, t0, t1, live))):
        raise ValueError("non-finite peak aperture or time")
    if not 0 < live <= t1 - t0 + 1e-5 or radius <= 0 or not 0 <= ccd <= 9:
        raise ValueError("invalid peak aperture or live exposure")
    area_ratio = radius**2 / (outer**2 - inner**2)
    # Keep region and pixel selection explicit in the saved command. Reapply
    # the exact half-open interval in Python at the inclusive CIAO filter edge.
    source = (
        f"{event_path}[EVENTS][ccd_id={ccd}][energy=500:7000]"
        f"[grade=0,2,3,4,6][status=0][time={t0:.17g}:{t1:.17g}]"
        f"[sky=circle({x:.17g},{y:.17g},{outer:.17g})]"
    )
    ciao_env = os.environ.get("CIAO_PREFIX")
    prefix_cmd = ["conda", "run", "--no-capture-output", "-p", ciao_env] if ciao_env else []
    commands = []
    def run(tool, *args):
        command = prefix_cmd + [tool, *args]
        commands.append(command)
        dump(out / "commands.json", commands)
        with (out / (tool + ".log")).open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                           check=True, cwd=out, timeout=900, stdin=subprocess.DEVNULL)

    run("ciaover", "-v")
    if ciao_env:
        package_metadata = []
        for path in sorted((Path(ciao_env) / "conda-meta").glob("*.json")):
            record = json.loads(path.read_text())
            package_metadata.append({key: record.get(key) for key in
                                     ("name", "version", "build", "url", "sha256", "md5")})
        dump(out / "environment.json", package_metadata)
    # Materialize the filtered events before calibration so the per-CCD
    # exposure metadata passed to eff2evt reflects the time filter.
    filtered = out / "peak_events.fits"
    weighted = out / "peak_weights.fits"
    run("dmcopy", source, str(filtered), "clobber=yes")
    run("eff2evt", str(filtered), str(weighted), "energy=INDEF", "dafile=CALDB",
        "rmfimg=", "verbose=2", "clobber=yes")
    filtered.unlink()
    with fits.open(weighted, memmap=False) as hdus:
        events = hdus["EVENTS"].data
        names = {name.upper() for name in events.names}
        if not {"DACORR", "EA", "QE", "FLUX"}.issubset(names):
            raise ValueError("CIAO output lacks required response weights")
        radial2 = (events["X"].astype(float) - x)**2 + (events["Y"].astype(float) - y)**2
        in_time = (events["TIME"] >= t0) & (events["TIME"] < t1)
        on = in_time & (radial2 <= radius**2)
        off = in_time & (radial2 >= inner**2) & (radial2 <= outer**2)
        if int(on.sum()) != trigger["peak_counts"] or int(off.sum()) != trigger["peak_annulus_counts"]:
            raise ValueError("CIAO and the blind search selected different peak photons")
        # Reconstruct the fluence per event from the calibrated response, then
        # divide by the measured peak GTI live seconds for this CCD.
        selected = on | off
        area = np.ones(len(events))
        for key in ("EA", "QE", "DACORR"):
            factor = events[key].astype(float)
            if np.any(~np.isfinite(factor[selected]) | (factor[selected] <= 0)):
                raise ValueError(f"undefined {key} for a selected photon")
            area[selected] *= factor[selected]
        flux_weight = events["ENERGY"].astype(float) * EV_IN_ERG / (area * live)
        if not on.any() or np.any(~np.isfinite(flux_weight[selected])):
            raise ValueError("empty source aperture or invalid calibrated photon")
        on_flux, off_flux = flux_weight[on], flux_weight[off]
        net = float(on_flux.sum() - area_ratio * off_flux.sum())
        # Mark each exact calibrated source/background row for audit.
        dump(out / "peak_selection.json", {
            "source_rows_zero_based": np.flatnonzero(on).tolist(),
            "background_rows_zero_based": np.flatnonzero(off).tolist(),
        })
    with weighted.open("rb") as data, (out / "peak_weights.fits.gz").open("wb") as packed:
        with gzip.GzipFile(filename="", mode="wb", fileobj=packed, mtime=0) as compressed:
            shutil.copyfileobj(data, compressed)
    weighted.unlink()
    result = {
        "peak_flux_erg_cm2_s": net,
        "flux_counting_sigma_erg_cm2_s": float(np.sqrt(np.dot(on_flux, on_flux)
                                                       + area_ratio**2 * np.dot(off_flux, off_flux))),
        "gross_flux_erg_cm2_s": float(on_flux.sum()),
        "scaled_background_flux_erg_cm2_s": float(area_ratio * off_flux.sum()),
        "aperture_counts": int(on.sum()), "annulus_counts": int(off.sum()),
        "band_keV": [0.5, 7.0], "peak_start_met": t0, "peak_stop_met": t1,
        "peak_live_s": live, "ccd_id": ccd, "sky_pixels": [x, y],
        "source_radius_pixels": radius, "annulus_radii_pixels": [inner, outer],
        "area_ratio": area_ratio, "event_filename": event_path.name,
        "event_sha256": digest(event_path), "calibrated_fits_sha256": digest(out / "peak_weights.fits.gz"),
        "interpretation": (
            "Observed aperture flux in the count-rate-selected peak interval. The annulus "
            "uses geometric areas. PSF losses, chip-gap coverage, pile-up, spectral "
            "redistribution, absorption and peak-selection bias have no correction. "
            "Counting sigma is a plug-in weighted Poisson standard error conditional on the window."
        ),
    }
    dump(out / "flux.json", result)
    return result
