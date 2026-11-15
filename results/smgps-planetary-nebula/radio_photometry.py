"""Check the integrated radio flux against independently dereddened optical photometry."""

import csv
import hashlib
import json
import math
from pathlib import Path
from urllib.request import urlopen

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

HERE = Path(__file__).resolve().parent
RADIO = HERE / "cache" / "radio.fits"
RADIO_URL = "https://archive-gw-1.kat.ac.za/public/repository/10.48479/0n8c-5q84/data/SMGPS_fits/SMGPS_Candidate_G276.229%2B00.419.fits"
RADIO_SHA = "6e332859389ac55d02603e9b5001b6bfa6ff1b7c7afe7365ce75bad381bd68e4"
POSITION = SkyCoord(145.505 * u.deg, -52.345 * u.deg)
# Goedhart et al. 2024 report a common 8 arcsec circular beam.
BEAM_FWHM_AS = 8.0
RADIUS_AS = 90
HALPHA_LOG, HALPHA_LOG_ERR = -11.58, 0.09
EBV, EBV_ERR = 0.75, 0.14
K_HALPHA = 2.53  # Cardelli R_V=3.1 at 656.3 nm
TE, NU_GHZ, HELIUM_FACTOR, BALMER_RATIO = 10000.0, 1.3, 1.10, 2.86


def model(log_halpha, ebv):
    hbeta = 10 ** (log_halpha + 0.4 * K_HALPHA * ebv) / BALMER_RATIO
    # Pottasch optically thin expression (Ruffle et al. 2004).
    # mW m^-2 and erg cm^-2 s^-1 are numerically equivalent.
    return 1000 * 2.51e7 * TE**0.53 * NU_GHZ**-0.1 * HELIUM_FACTOR * hbeta


def main():
    RADIO.parent.mkdir(exist_ok=True)
    if not RADIO.exists():
        with urlopen(RADIO_URL, timeout=90) as remote:
            RADIO.write_bytes(remote.read())
    if hashlib.sha256(RADIO.read_bytes()).hexdigest() != RADIO_SHA:
        raise ValueError("Radio image differs from the reported release")
    with fits.open(RADIO) as hdus:
        im, hdr = hdus[0].data, hdus[0].header
    if hdr["BUNIT"].strip().upper() != "JY/BEAM":
        raise ValueError("Unrecognized intensity unit")
    cx, cy = WCS(hdr, naxis=2).world_to_pixel(POSITION)
    y, x = np.indices(im.shape)
    dx, dy = (x - cx) * abs(hdr["CDELT1"]) * 3600, (y - cy) * abs(hdr["CDELT2"]) * 3600
    rad = np.hypot(dx, dy)
    az = np.arctan2(dy, dx)
    pix_area = abs(hdr["CDELT1"] * hdr["CDELT2"]) * 3600**2
    beam_area = math.pi / (4 * math.log(2)) * BEAM_FWHM_AS**2
    jybeam_to_mjy = 1000 * pix_area / beam_area

    def measure(radius, sky):
        mask = rad <= radius
        bg = float(np.median(im[sky]))
        return float(np.sum(im[mask] - bg) * jybeam_to_mjy), bg * 1e6

    data = []
    for aperture in [70, 80, 85, 90, 95, 100, 110]:
        for inner, outer in [(110, 135), (120, 150), (135, 165), (150, 175)]:
            sky = (rad > inner) & (rad < outer)
            flux, bg = measure(aperture, sky)
            data.append((aperture, inner, outer, "all", flux, bg))
    standard_sky = (rad > 120) & (rad < 150)
    sector_fluxes = []
    for i in range(8):
        sec = (
            standard_sky
            & ((az + math.pi) >= i * math.pi / 4)
            & ((az + math.pi) < (i + 1) * math.pi / 4)
        )
        flux, bg = measure(RADIUS_AS, sec)
        sector_fluxes.append(flux)
        data.append((RADIUS_AS, 120, 150, i, flux, bg))
    with (HERE / "radio-apertures.csv").open("w", newline="") as out:
        w = csv.writer(out, lineterminator="\n")
        w.writerow(
            [
                "aperture_radius_arcsec",
                "sky_inner_arcsec",
                "sky_outer_arcsec",
                "sky_octant",
                "flux_mjy",
                "sky_ujy_per_beam",
            ]
        )
        w.writerows(data)
    observed, _ = measure(RADIUS_AS, standard_sky)
    predicted = model(HALPHA_LOG, EBV)
    siglog = math.hypot(HALPHA_LOG_ERR, 0.4 * K_HALPHA * EBV_ERR)
    aperture_fluxes = [r[4] for r in data if r[3] == "all" and r[0] >= 85]
    result = {
        "radio_sha256": RADIO_SHA,
        "aperture_radius_arcsec": RADIUS_AS,
        "sky_annulus_arcsec": [120, 150],
        "common_beam_fwhm_arcsec": BEAM_FWHM_AS,
        "flux_mjy": observed,
        "flux_range_85_to_110_arcsec_across_four_annuli_mjy": [
            min(aperture_fluxes),
            max(aperture_fluxes),
        ],
        "flux_range_90_arcsec_aperture_eight_sky_octants_mjy": [
            min(sector_fluxes),
            max(sector_fluxes),
        ],
        "annular_flux_from_90_to_110_arcsec_mjy": measure(110, standard_sky)[0]
        - observed,
        "optical_halpha_log10_flux_cgs": HALPHA_LOG,
        "optical_halpha_error_dex": HALPHA_LOG_ERR,
        "cs_photometry_E_BV": EBV,
        "cs_photometry_E_BV_error": EBV_ERR,
        "thermal_model_fixed_assumptions": {
            "Te_kelvin": TE,
            "nu_GHz": NU_GHZ,
            "Halpha_Hbeta": BALMER_RATIO,
            "helium_factor": HELIUM_FACTOR,
            "A_Halpha_over_E_BV": K_HALPHA,
        },
        "predicted_freefree_mjy": predicted,
        "predicted_1sigma_from_flux_and_reddening_mjy": [
            predicted / 10**siglog,
            predicted * 10**siglog,
        ],
        "observed_over_predicted": observed / predicted,
        "implied_E_BV_for_all_thermal": math.log10(observed / model(HALPHA_LOG, 0))
        / (0.4 * K_HALPHA),
    }
    (HERE / "radio-photometry.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
