"""Check the reversal on the additional 20-second TESS light-curve products."""
import argparse
import csv
from pathlib import Path

from astropy.io import fits
import numpy as np

from measure import ROOT, amplitudes, prepare, verified_product
from slow_signal import write_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=ROOT/"fast_check.csv")
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    with (ROOT/"fast_manifest.csv").open() as handle:
        entries = list(csv.DictReader(handle))
    rows = []
    for row in entries:
        path = verified_product(row, args.cache)
        data = fits.getdata(path, 1)
        for key in ("SAP_FLUX", "PDCSAP_FLUX"):
            t, y = prepare(data, key)
            # Full coverage of each two-minute bin, built from independent SPOC
            # 20-second products. Both cadence streams record the same star and time.
            ix = np.floor((t - t[0]) * 720 + 1e-7).astype(int)
            n = np.bincount(ix)
            good = n >= 5
            tx = np.bincount(ix, weights=t)[good] / n[good]
            yx = np.bincount(ix, weights=y)[good] / n[good]
            _, amp = amplitudes(tx, yx)
            rows.append(dict(sector=row["sector"], flux=key, fast_cadences=len(t),
                             binned_cadences=len(tx), low_ppm=amp[0], high_ppm=amp[1],
                             ratio=amp[0]/amp[1]))
    write_csv(args.out, rows)


if __name__ == "__main__":
    main()
