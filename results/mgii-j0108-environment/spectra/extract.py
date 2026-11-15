"""Extract local vacuum-wavelength windows from identified public DR1/DR17 spectra."""
import argparse
import csv
import hashlib
from pathlib import Path
from urllib.request import urlretrieve

from astropy.io import fits
import numpy as np

HERE = Path(__file__).resolve().parent
DESI_URL = ("https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix/"
            "main/dark/179/17997/coadd-main-dark-17997.fits")
SDSS_URL = "https://data.sdss.org/sas/dr17"
SIGHTLINES = {
    "A": (39627754810908189, 0.9601519467, [(694, 52209, 122), (3736, 55214, 291), (4226, 55475, 228)]),
    "B": (39627754815095242, 0.9604029861, [(694, 52209, 128), (3736, 55214, 284), (4226, 55475, 230)]),
    "C": (39627760854893276, 0.9591435925, [(4226, 55475, 232)]),
}
# Boundaries selected to contain Fe II 2344,2374,2382,2586,2600 and Mg II 2796,2803, Mg I 2852.
WINDOWS = [(4535, 4725), (5000, 5160), (5415, 5545), (5530, 5650)]


def write_csv(name, header, rows):
    with (HERE / name).open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def pixels(wave, flux, ivar, mask):
    selected = np.logical_or.reduce([(wave >= lo) & (wave <= hi) for lo, hi in WINDOWS])
    return zip(wave[selected], flux[selected], ivar[selected], mask[selected])


def get_file(cache, name, url):
    file = cache / name
    if not file.exists():
        temporary = file.with_suffix(".partial")
        urlretrieve(url, temporary)
        temporary.replace(file)
    return file


def digest(file):
    with file.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def extract(cache):
    cache.mkdir(exist_ok=True, parents=True)
    coadd = get_file(cache, "coadd.fits", DESI_URL)
    sources = [["desi", coadd.name, digest(coadd), DESI_URL]]
    data, visits = [], []
    with fits.open(coadd) as hdus:
        fibers = hdus["FIBERMAP"].data
        exposures = hdus["EXP_FIBERMAP"].data
        for beam, (target, redshift, observations) in SIGHTLINES.items():
            row, = np.flatnonzero(fibers["TARGETID"] == target)
            specid = f"DESI-{target}"
            # All windows lie in the blue arm below 5650 Angstrom.
            values = pixels(hdus["B_WAVELENGTH"].data, hdus["B_FLUX"].data[row],
                            hdus["B_IVAR"].data[row], hdus["B_MASK"].data[row])
            data.extend((beam, specid, *v) for v in values)
            ra, dec = fibers[row]["TARGET_RA"], fibers[row]["TARGET_DEC"]
            for exposure in exposures[exposures["TARGETID"] == target]:
                visits.append([beam, specid, redshift, ra, dec, exposure["MJD"], exposure["EXPID"],
                               exposure["FIBER"], exposure["FIBERSTATUS"]])
            for plate, mjd, fiber in observations:
                name = f"spec-{plate:04d}-{mjd}-{fiber:04d}.fits"
                run = "sdss/spectro/redux/26/spectra" if plate == 694 else "eboss/spectro/redux/v5_13_2/spectra/lite"
                url = f"{SDSS_URL}/{run}/{plate:04d}/{name}"
                file = get_file(cache, name, url)
                specid = name[:-5]
                sources.append([specid, name, digest(file), url])
                with fits.open(file) as sdss:
                    table = sdss["COADD"].data
                    header = sdss[0].header
                    sdss_ra, sdss_dec = header["PLUG_RA"], header["PLUG_DEC"]
                    sep = np.hypot((ra - sdss_ra) * np.cos(np.deg2rad(dec)), dec - sdss_dec) * 3600
                    assert sep < 1.0 and header["MJD"] == mjd and header["FIBERID"] == fiber
                    visits.append([beam, specid, redshift, sdss_ra, sdss_dec, mjd, "", fiber, ""])
                    values = pixels(10**table["loglam"].astype(float), table["flux"], table["ivar"], table["and_mask"])
                    data.extend((beam, specid, *v) for v in values)
    (HERE / "pixels").mkdir(exist_ok=True)
    for spectrum in sorted({row[1] for row in data}):
        write_csv(f"pixels/{spectrum}.csv", ["beam","spectrum","wave_A","flux","ivar","mask"],
                  (row for row in data if row[1] == spectrum))
    write_csv("visits.csv", ["beam","spectrum","z_abs_desi","ra_deg","dec_deg","mjd","expid","fiber","fiberstatus"], visits)
    write_csv("sources.csv", ["product","cache_filename","sha256","url"], sources)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    extract(parser.parse_args().cache)
