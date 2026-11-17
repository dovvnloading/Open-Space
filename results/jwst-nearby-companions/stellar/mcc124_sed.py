"""Estimate blended host F444W flux from 2MASS Ks and AllWISE W2.

Run: python mcc124_sed.py /path/to/model_cache > mcc124_sed.csv
Requires numpy>=2 and astropy. Inputs and digests in photometry_inputs.json.
"""
import csv
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import urlopen
import xml.etree.ElementTree as ET

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent


def get_inputs(cache):
    cache.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((HERE / "photometry_inputs.json").read_text())
    for name, source in manifest.items():
        dest = cache / name
        if not dest.is_file():
            with urlopen(source["url"], timeout=90) as resp:
                data = resp.read()
            if hashlib.sha256(data).hexdigest() != source["sha256"]:
                raise ValueError("Input changed: " + name)
            dest.write_bytes(data)
        elif hashlib.sha256(dest.read_bytes()).hexdigest() != source["sha256"]:
            raise ValueError("Input changed: " + name)


def get_band(path):
    document = ET.parse(path)
    params = {n.get("name"): n.get("value") for n in document.findall(".//{*}PARAM")}
    wav, response = np.array([[float(k.text) for k in row.findall("{*}TD")]
                              for row in document.findall(".//{*}TR")]).T
    return wav, response, int(params["DetectorType"]), float(params["ZeroPoint"])


def band_integral(wavelength, flam, band):
    wav, response, photon, _ = band
    # Integrate on the union so that both the spectrum and filter features are retained.
    use = wavelength[(wavelength > wav[0]) & (wavelength < wav[-1])]
    grid = np.union1d(use, wav)
    weights = np.interp(grid, wav, response)
    if photon == 1:
        weights *= grid
    return np.trapezoid(np.interp(grid, wavelength, flam) * weights, grid)


def main(cache):
    get_inputs(cache)
    bands = [get_band(cache / (name + ".xml")) for name in ("ks", "w2", "f444w")]
    with fits.open(cache / "vega.fits") as f:
        w, flux = f[1].data["WAVELENGTH"], f[1].data["FLUX"]
        vega = np.array([band_integral(w, flux, b) for b in bands])
    with (HERE / "catalogue_receipts" / "twomass.csv").open() as f:
        ks = float(next(csv.DictReader(f))["k_m"])
    with (HERE / "catalogue_receipts" / "allwise.csv").open() as f:
        w2 = float(next(csv.DictReader(f))["w2mpro"])
    out = csv.writer(sys.stdout, lineterminator="\n")
    out.writerow(["teff_K", "logg", "f444w_minus_ks_vega", "f444w_minus_w2_vega",
                  "ks_anchored_f444w_mJy", "w2_anchored_f444w_mJy"])
    for name in sorted(cache.glob("ph[0-9][0-9][0-9][0-9].fits")):
        with fits.open(name) as f:
            table = f[1].data
            for logg, col in ((4.5, "g45"), (5., "g50")):
                observed = [band_integral(table["WAVELENGTH"], table[col], b) for b in bands]
                mags = -2.5 * np.log10(observed / vega)
                dks, dw2 = mags[2] - mags[0], mags[2] - mags[1]
                jy0 = bands[2][3]
                fluxes = 1000 * jy0 * 10 ** (-0.4 * np.array([ks + dks, w2 + dw2]))
                out.writerow([name.stem[2:], logg, f"{dks:.4f}", f"{dw2:.4f}",
                              *(f"{v:.2f}" for v in fluxes)])


if __name__ == "__main__":
    main(Path(sys.argv[1]))
