"""Regenerate the narrow spectral windows from the checksum-pinned FITS files."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from measure import C_KMS, HERE, LINES


def samples(path, entry):
    with fits.open(path, memmap=False) as hdus:
        if entry["source"] == "XSHOOTER":
            row = hdus[1].data[0]
            wave = row["WAVE"].astype(float) * 10
            flux = row["FLUX"].astype(float) * 1e17
            error = row["ERR"].astype(float) * 1e17
            good = (row["QUAL"] == 0) & np.isfinite(flux) & np.isfinite(error) & (error > 0)
            weight = np.divide(1, error**2, out=np.zeros_like(error), where=good)
            yield wave, np.nan_to_num(flux), weight
        else:
            for exp_id in entry["exposures"]:
                hdu = next(ext for ext in hdus if ext.header.get("EXPOSURE") == exp_id)
                data = hdu.data
                vacuum = 10**data["LOGLAM"].astype(float)
                sigma_sq = (10000/vacuum)**2
                wave = vacuum/(1.0000834254 + 0.02406147/(130-sigma_sq) +
                               0.00015998/(38.9-sigma_sq))
                flux = data["FLUX"]
                ivar = data["IVAR"]
                good = (data["AND_MASK"] == 0) & np.isfinite(flux) & np.isfinite(ivar) & (ivar > 0)
                yield wave, np.nan_to_num(flux), np.where(good, ivar, 0)


def main(directory):
    files = {}
    for path in sorted(directory.glob("*.fits")):
        files[hashlib.sha256(path.read_bytes()).hexdigest()] = path
    entries = json.loads((HERE / "observations.json").read_text())
    missing = [e["url"] for e in entries if e["sha256"] not in files]
    if missing:
        raise SystemExit("Missing original FITS files or digest mismatch:\n" + "\n".join(missing))
    counters = {"XSHOOTER": 0, "BOSS": 0}
    with (HERE / "line-pixels.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["source", "epoch", "line_air_A", "wavelength_air_A", "flux", "ivar"])
        for entry in entries:
            src = entry["source"]
            for wave, flux, weight in samples(files[entry["sha256"]], entry):
                for line in LINES:
                    use = np.flatnonzero(np.abs(C_KMS*np.log(wave/line)) < 410)
                    for j in use:
                        writer.writerow([src, counters[src], f"{line:.3f}", f"{wave[j]:.7f}",
                                         f"{flux[j]:.7g}", f"{weight[j]:.7g}"])
                counters[src] += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    main(parser.parse_args().directory)
