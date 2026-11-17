"""Photon-counting, band-averaged flux density from an archived CALSPEC SED."""

from functools import lru_cache
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table


@lru_cache(maxsize=8)
def band_jy(model: Path, bandpass: Path):
    with fits.open(model) as h:
        sed = h[1].data
        wavelength = np.asarray(sed["WAVELENGTH"], dtype=float)
        flux = np.asarray(sed["FLUX"], dtype=float)
        assert h[1].columns["WAVELENGTH"].unit.lower() in ("angstrom", "angstroms")
    band = Table.read(bandpass, format="votable")
    wave = band["Wavelength"].quantity.to_value("Angstrom")
    transmission = np.asarray(band["Transmission"], dtype=float)
    if wavelength[0] > wave[0] or wavelength[-1] < wave[-1]:
        raise ValueError("Spectrum does not cover the entire bandpass.")
    grid = np.unique(
        np.concatenate(
            [wave, wavelength[(wavelength > wave[0]) & (wavelength < wave[-1])]]
        )
    )
    response = np.interp(grid, wave, transmission)
    flam = np.interp(grid, wavelength, flux)
    return float(
        np.trapezoid(flam * grid * response, grid)
        / np.trapezoid(response / grid, grid)
        / 2.99792458e18
        * 1e23
    )
