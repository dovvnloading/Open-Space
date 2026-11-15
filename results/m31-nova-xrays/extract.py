"""Extract photon positions around fixed optical nova positions from Chandra HRC."""
import argparse
import csv
import gzip
import hashlib
from pathlib import Path
import urllib.request

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache", type=Path, help="Directory for the official compressed FITS files")
    p.add_argument("output", type=Path, help="Output event cutout CSV.gz")
    args = p.parse_args()
    root = Path(__file__).resolve().parent
    targets = list(csv.DictReader((root / "targets.csv").open()))
    files = list(csv.DictReader((root / "observations.csv").open()))
    args.cache.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, "wt", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["obsid", "target", "x_offset_arcsec", "y_offset_arcsec", "time_met", "pi"])
        for obs in files:
            filename = obs["url"].rsplit("/", 1)[1]
            file = args.cache / filename
            if not file.exists():
                urllib.request.urlretrieve(obs["url"], file)
            digest = hashlib.sha256()
            with file.open("rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
            if digest.hexdigest() != obs["sha256"]:
                raise ValueError(f"Archive file digest differs for {filename}")
            with fits.open(file) as product:
                events = product["EVENTS"]
                wcs = WCS(naxis=2)
                xy = [events.columns.names.index(axis) + 1 for axis in ("x", "y")]
                for prop, keyword in (("ctype","TCTYP"), ("crval","TCRVL"), ("crpix","TCRPX"), ("cdelt","TCDLT")):
                    setattr(wcs.wcs, prop, [events.header[f"{keyword}{col}"] for col in xy])
                scale = abs(wcs.wcs.cdelt[0]) * 3600
                for target in targets:
                    x, y = wcs.all_world2pix([[float(target["ra_deg"]), float(target["dec_deg"])]], 1)[0]
                    dx = (events.data["x"].astype("float64") - x) * scale
                    dy = (events.data["y"].astype("float64") - y) * scale
                    rows = np.flatnonzero(dx * dx + dy * dy <= 100)
                    for i in rows:
                        writer.writerow([obs["obsid"], target["name"], f"{dx[i]:.6f}", f"{dy[i]:.6f}",
                                         f"{events.data['time'][i]:.6f}", events.data["pi"][i]])


if __name__ == "__main__":
    main()
