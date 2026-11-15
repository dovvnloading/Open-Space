"""Measure a restoring-beam template at known nova coordinates in native VLASS pixels."""
import argparse
import csv
import hashlib
import io
import sys
from datetime import datetime
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
CUTOUT = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/sync"


def gaussian(x, y, angle, major, minor):
    # x increases east and y north. FITS BPA is measured east of north.
    u = x * np.sin(angle) + y * np.cos(angle)
    v = x * np.cos(angle) - y * np.sin(angle)
    return np.exp(-4 * np.log(2) * ((u / major)**2 + (v / minor)**2))


def estimate(data, east, north, angle, major, minor, dx=0, dy=0):
    x, y = east - dx, north - dy
    region = (x*x + y*y < 6**2) & np.isfinite(data)
    template = gaussian(x[region], y[region], angle, major, minor)
    amplitude, background = np.linalg.lstsq(
        np.column_stack((template, np.ones(template.size))), data[region], rcond=None
    )[0]
    return float(amplitude)


def measure(record, path):
    with fits.open(path) as hdus:
        header = hdus[0].header
        assert header["BUNIT"].lower() == "jy/beam"
        data = hdus[0].data.squeeze().astype(float) * 1000
        wcs = WCS(header).celestial
    ra, dec = float(record["ra"]), float(record["dec"])
    # Gnomonic tangent offsets relative to target using the local WCS Jacobian.
    x, y = wcs.world_to_pixel_values(ra, dec)
    ys, xs = np.indices(data.shape)
    east, north = np.einsum('ij,jkl->ikl', wcs.pixel_scale_matrix * 3600,
                            np.stack((xs-x, ys-y)))
    major, minor = header["BMAJ"]*3600, header["BMIN"]*3600
    angle = np.deg2rad(header["BPA"])
    amplitude = estimate(data, east, north, angle, major, minor)
    # Disjoint centers 8 arcsec apart, retaining a complete fitting aperture.
    controls = [
        estimate(data, east, north, angle, major, minor, dx, dy)
        for dx in range(-24, 25, 8) for dy in range(-24, 25, 8)
        if 64 <= dx*dx + dy*dy <= 28**2
    ]
    median = np.median(controls)
    sigma = 1.482602218 * np.median(np.abs(np.array(controls)-median))
    radius = np.hypot(east, north)
    sky = data[(radius>10) & (radius<35) & np.isfinite(data)]
    image_rms = 1.482602218*np.median(np.abs(sky-np.median(sky)))
    mask = (radius <= 6) & np.isfinite(data)
    def residual(p):
        a, dx, dy, offset = p
        return a*gaussian(east[mask]-dx,north[mask]-dy,angle,major,minor)+offset-data[mask]
    fit = least_squares(residual,[max(amplitude,0),0,0,np.median(sky)],
                        bounds=([0,-2,-2,-np.inf],[np.inf,2,2,np.inf]))
    a,dx,dy,bg = fit.x
    when = datetime.fromisoformat(header['DATE-OBS'])
    delay = ((when-datetime.fromisoformat(record['optical_date'])).total_seconds()/86400
             if record['optical_date'] else None)
    return dict(name=record["name"],epoch=record["cutout"].split('VLASS')[-1][:-5],
                utc=when.isoformat(timespec='seconds'),days=round(delay,2) if delay is not None else '',
                forced_mJy_beam=round(amplitude,3),image_rms_mJy_beam=round(image_rms,3),
                control_sigma_mJy_beam=round(float(sigma),3),control_count=len(controls),
                max_abs_control_mJy_beam=round(max(abs(z) for z in controls),3),
                peak_pixel_mJy_beam=round(np.nanmax(data[radius<2]),3),
                free_centroid_mJy_beam=round(float(a),3) if amplitude > 5*image_rms else '',
                centroid_separation_as=round(float(np.hypot(dx,dy)),3) if amplitude > 5*image_rms else '',
                major_as=round(major,3),minor_as=round(minor,3),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch",action="store_true",help="Download a fresh CADC cutout when the local file is missing")
    args = parser.parse_args()
    with (ROOT/"observations.csv").open() as stream:
        records = list(csv.DictReader(stream))
    output = []
    for record in records:
        path = ROOT/"cutouts"/record["cutout"]
        if args.fetch and not path.exists():
            query = urlencode({"ID":record["uri"],"CIRCLE":f"{record['ra']} {record['dec']} 0.0125"})
            with urlopen(CUTOUT+"?"+query,timeout=90) as response:
                data=response.read()
            with fits.open(io.BytesIO(data)) as hdus:
                hdus.verify("exception")
            path.parent.mkdir(exist_ok=True)
            path.write_bytes(data)
        output.append(measure(record,path))
    writer=csv.DictWriter(sys.stdout,fieldnames=list(output[0]),lineterminator="\n")
    writer.writeheader()
    writer.writerows(output)


if __name__ == "__main__":
    main()
