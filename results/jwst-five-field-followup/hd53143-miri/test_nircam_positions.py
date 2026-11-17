"""Measure F2550W response at each April F444W peak under both motion hypotheses."""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from astropy.time import Time
from scipy import ndimage

from freeze_inputs import HERE, PREVIOUS, main as verify
sys.path.insert(0, str(PREVIOUS))
import reduce as rdi
sys.path.insert(0, str(PREVIOUS.parent / "vetting"))
from other_tracks import parpos, offset

NIRCAM = HERE.parent / "nircam-hd53143-hd38206/outputs"


def main(cache):
    rows = verify(cache, False)
    with (NIRCAM / "HD53143-gaia.csv").open() as stream:
        star = next(r for r in csv.DictReader(stream) if r["source_id"] == "5479222240596469632")
    std = [rdi.prepare(cache / r["filename"]) for r in rows if r["role"] == "standard"]
    result = []
    for s in (r for r in rows if r["role"] in ("roll", "combined")):
        sci = rdi.prepare(cache / s["filename"])
        flux = rdi.core(sci["image"]).sum()
        kernel = rdi.core(std[0]["image"])
        kernel *= flux / kernel.sum()
        use = rdi.R[rdi.PATCH] <= rdi.CUT
        kernel[use] -= kernel[use].mean()
        amplitude = rdi.matched(sci, rdi.radial_reference(std), kernel)
        support = ndimage.minimum_filter(sci["valid"].astype(int), 2*rdi.CUT+5, mode="constant").astype(bool)
        visit = s["obs_id"].split("_")[0]
        curves = list(csv.DictReader((HERE / "measurements" / visit / "contrast_recovery.csv").open()))
        for sourcefile in sorted(NIRCAM.glob("hd53143*_peaks.csv")):
            meta = json.loads(sourcefile.with_name(sourcefile.name.replace("_peaks.csv", "_measurement.json")).read_text())
            pos1 = parpos(star, Time(meta["expstart_mjd"], format="mjd"))
            pos2 = parpos(star, Time(sci["mjd"], format="mjd"))
            for i, candidate in enumerate(csv.DictReader(sourcefile.open()), 1):
                rel = np.array([float(candidate[k]) for k in ("east_arcsec", "north_arcsec")])
                for track, predicted in [("comoving", rel), ("fixed", rel - offset(pos2, pos1))]:
                    x, y = rdi.HALF + np.linalg.solve(sci["mat"], predicted)
                    radius = np.linalg.norm(predicted)
                    ann = next((k for k in curves if float(k["inner_arcsec"]) <= radius < float(k["outer_arcsec"])), None)
                    record = dict(visit=visit, sourcefile=sourcefile.name, peak=i, track=track,
                                  east_arcsec=predicted[0], north_arcsec=predicted[1], radius_arcsec=radius)
                    valid = support[round(y), round(x)] and ann is not None
                    record.update(inside_search=int(valid), forced_contrast="", local_sigma="", forced_sigma="", threshold_sigma="")
                    if valid:
                        reg = support & (rdi.R*sci["scale"] >= float(ann["inner_arcsec"])) & (rdi.R*sci["scale"] < float(ann["outer_arcsec"]))
                        med = np.median(amplitude[reg])
                        sigma = 1.4826*np.median(abs(amplitude[reg]-med))
                        amp = ndimage.map_coordinates(amplitude, [[y], [x]], order=1)[0] - med
                        record.update(forced_contrast=amp, local_sigma=sigma, forced_sigma=amp/sigma,
                                      threshold_sigma=ann["threshold_sigma"])
                    result.append(record)
    rdi.save_table(HERE / "nircam_forced_tracks.csv", result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache", type=Path)
    main(p.parse_args().cache)
