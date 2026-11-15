"""Verify the archived STIS exposure lengths and image diagnostic from MAST FITS files."""
import argparse
import csv
import hashlib
from itertools import groupby
from pathlib import Path
from urllib.request import urlopen

from astropy.io import fits
from scipy.ndimage import median_filter

ROOT = Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    cache=parser.parse_args().cache
    cache.mkdir(parents=True,exist_ok=True)
    with (ROOT/"stis-exposures.csv").open() as stream:
        rows=list(csv.DictReader(stream))
    for name, group in groupby(rows,key=lambda row:row["product"]):
        entries=list(group)
        path=cache/name
        if not path.exists():
            url="https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/"+name
            path.write_bytes(urlopen(url,timeout=90).read())
        assert hashlib.sha256(path.read_bytes()).hexdigest()==entries[0]["sha256_first_row"]
        with fits.open(path) as product:
            for row in entries:
                hdu=product["SCI",int(row["SCI_EXTVER"])]
                assert hdu.header["BUNIT"]=="COUNTS"
                assert hdu.header["EXPTIME"]==float(row["exposure_s"])
                assert hdu.header["EXPSTART"]==float(row["EXPSTART_MJD"])
                peak=median_filter(hdu.data[35:155,455:565],size=3).max()
                assert abs(float(peak)-float(row["median3_peak_counts"]))<0.0051
    print(f"Verified {len(rows)} calibrated SCI extensions.")

if __name__=="__main__":
    main()
