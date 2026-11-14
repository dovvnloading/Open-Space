"""Measure GX 1+4 pulses from the six NuSTAR event files listed in sources.json."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import urlretrieve

import numpy as np
from astropy import units as u
from astropy.coordinates import get_body_barycentric_posvel
from astropy.io import fits
from astropy.time import Time
from scipy.optimize import minimize_scalar


def spectrum(counts, exposure):
    """Sum Fourier power at the fundamental and second harmonic on a broad grid."""
    residual = counts - exposure * counts.sum() / exposure.sum()
    sample_time = 2.0
    fft_size = 2**20
    transform = np.fft.rfft(residual, n=fft_size)
    frequencies = np.fft.rfftfreq(fft_size, sample_time)
    band = np.flatnonzero((frequencies >= 0.0039) & (frequencies <= 0.005))
    power = abs(transform[band])**2 + abs(transform[2*band])**2
    peak = band[np.argmax(power)]
    seconds = (np.arange(counts.size)+0.5)*sample_time

    def negative_power(frequency):
        phasor = np.exp(2j * np.pi * seconds * frequency)
        return -2 * (abs(residual @ phasor)**2 + abs(residual @ (phasor**2))**2) / counts.sum()

    result = minimize_scalar(negative_power, bounds=frequencies[[peak-1,peak+1]],
                             method="bounded", options={"xatol":1e-12})
    return result.x, -result.fun


def measure(path):
    with fits.open(path) as hdus:
        events, header = hdus[1].data, hdus[1].header
        # The targets of these pointings coincide with the tangent-plane reference
        # pixel to less than one arcsecond. Choose sky radius in arcseconds.
        pixscale = 3600 * abs(header["TCDLT39"])
        radius = pixscale*np.hypot(events["X"]-header["TCRPX38"], events["Y"]-header["TCRPX39"])
        pi = events["PI"]
        good = (pi >= 35) & (pi < 835)  # E = 1.6 + 0.04 PI, in keV
        gti = hdus["GTI"].data
        origin = gti["START"][0]
        span = gti["STOP"][-1] - origin
        bins = int(np.ceil(span/2)) + 1
        indices = ((events["TIME"]-origin)/2).astype(int)
        source = np.bincount(indices[good & (radius<60)], minlength=bins)
        annulus = np.bincount(indices[good & (radius>100) & (radius<180)], minlength=bins)
        edges = np.arange(bins+1)*2.0
        exposure = np.zeros(bins)
        for start,stop in zip(gti["START"]-origin, gti["STOP"]-origin):
            exposure += np.maximum(0, np.minimum(edges[1:],stop)-np.maximum(edges[:-1],start))
        epoch = header["MJDREFI"]+header["MJDREFF"]+(origin+span/2)/86400
        ra,dec = np.deg2rad([header["RA_OBJ"],header["DEC_OBJ"]])
    unit = np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])
    _,velocity = get_body_barycentric_posvel("earth", Time(epoch,format="mjd",scale="tt"))
    projected = velocity.xyz.to_value(u.km/u.s) @ unit
    raw_frequency, power = spectrum(source, exposure)
    annulus_frequency, annulus_power = spectrum(annulus, exposure)
    conversion = 1e6/(1+projected/299792.458)
    return [path.name, f"{epoch:.6f}",source.sum(),f"{raw_frequency*conversion:.5f}",
            f"{power:.2f}",annulus.sum(),f"{annulus_frequency*conversion:.5f}",f"{annulus_power:.2f}"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--download", action="store_true", help="retrieve missing mission event files")
    args = parser.parse_args()
    sources = json.loads(Path(__file__).with_name("sources.json").read_text())
    output = csv.writer(sys.stdout, lineterminator="\n")
    output.writerow(["file","mjd_tt","source_photons","frequency_uHz","source_power",
                     "annulus_photons","annulus_frequency_uHz","annulus_power"])
    for ref in sources["nustar"]:
        path = args.cache/ref["url"].split("/")[-1]
        if not path.exists() and args.download:
            args.cache.mkdir(parents=True,exist_ok=True)
            urlretrieve(ref["url"],path)
        with path.open("rb") as stream:
            sha = hashlib.sha256()
            for block in iter(lambda:stream.read(1048576), b""):
                sha.update(block)
        if sha.hexdigest() != ref["sha256"]:
            raise ValueError(f"NuSTAR input differs from recorded source bytes: {path}")
        output.writerow(measure(path))
