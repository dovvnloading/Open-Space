"""Fixed-redshift Mg II aperture check on neighboring quasars in the discovery healpixel."""
import argparse
from pathlib import Path
import sys

from astropy.io import fits
import numpy as np

from environment import angular_distance
from search import directions, write_csv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "spectra"))
from measure import measure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("redrock", type=Path)
    parser.add_argument("coadd", type=Path)
    parser.add_argument("--out", type=Path, default=HERE)
    args = parser.parse_args()
    absorbers = np.genfromtxt(args.out / "j0108-absorbers.csv", delimiter=",", names=True, dtype=None)
    center = directions(absorbers["RA"].mean(), absorbers["DEC"].mean())
    reference = np.median(absorbers["Z_MGII"])
    rows = []
    with fits.open(args.redrock) as redrock, fits.open(args.coadd) as coadd:
        r = redrock["REDSHIFTS"].data
        f = redrock["FIBERMAP"].data
        idx = {int(x): i for i, x in enumerate(coadd["FIBERMAP"].data["TARGETID"])}
        if not np.array_equal(r["TARGETID"], f["TARGETID"]):
            raise ValueError("Misaligned redrock arrays")
        radius = angular_distance(directions(f["TARGET_RA"], f["TARGET_DEC"]), center)
        include = ((radius <= 12) & (r["ZWARN"] == 0) & (r["SPECTYPE"] == "QSO")
                   & (299792.458*np.log((1+r["Z"])/(1+reference)) > 5000))
        wave = coadd["B_WAVELENGTH"].data
        win = abs(wave - 2800*(1+reference)) < 70
        for j in np.where(include)[0]:
            targetid = r["TARGETID"][j]
            i = idx[int(targetid)]
            values = [coadd[key].data[i, win] for key in ["B_FLUX","B_IVAR","B_MASK"]]
            measurements = {rest: (ew, err) for rest, ew, err in
                            measure(wave[win], *values, reference, "MgII", 325, "linear")}
            row = [targetid, f["TARGET_RA"][j], f["TARGET_DEC"][j], r["Z"][j], radius[j]]
            for rest in (2796.3543, 2803.5315):
                ew, err = measurements.get(rest, (np.nan, np.nan))
                row.extend((ew, err))
            rows.append(row)
    write_csv(args.out / "neighboring-sightlines.csv",
              ["TARGETID","ra_deg","dec_deg","z_qso","radius_arcmin","EW2796_A","err2796_A",
               "EW2803_A","err2803_A"], sorted(rows, key=lambda row: row[4]))


if __name__ == "__main__":
    main()
