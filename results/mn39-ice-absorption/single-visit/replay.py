"""Check the 4.26-micron dip in one four-exposure SPHEREx observing sequence."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy import log
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time

log.setLevel("WARNING")
ROOT = Path(__file__).resolve().parent
SOURCES = [('MN39', 242.61055777, -51.35703930),
           ('AllWISE J161036.07-512153.8', 242.650318, -51.364950)]
SOURCE_BIT = np.uint32(0x200000)


def extract(source, dest):
    with (ROOT / 'manifest.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    dest.mkdir(parents=True, exist_ok=True)
    source.mkdir(parents=True, exist_ok=True)
    for row in rows:
        path = source / row['filename']
        if not path.exists():
            with urlopen(row['url'], timeout=120) as response:
                data = response.read()
            if hashlib.sha256(data).hexdigest() != row['sha256']:
                raise ValueError(path.name)
            path.write_bytes(data)
        if hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError(path.name)
        with fits.open(path) as hdul:
            header = hdul['IMAGE'].header
            celestial = WCS(header, naxis=2)
            positions = np.stack(celestial.world_to_pixel_values(
                [s[1] for s in SOURCES], [s[2] for s in SOURCES]), axis=1)
            wave_wcs = WCS(header, hdul, key='W')
            # The alternate WCS is tabulated in detector coordinates. wcs_pix2world
            # bypasses the celestial SIP transformation, already used for positions.
            if [str(unit) for unit in wave_wcs.wcs.cunit] != ['um', 'um']:
                raise ValueError('Unexpected spectral WCS units')
            wave = wave_wcs.wcs_pix2world(positions, 0)
            low = np.rint(positions).astype(int).min(axis=0) - 10
            high = np.rint(positions).astype(int).max(axis=0) + 11
            cut = np.s_[low[1]:high[1], low[0]:high[0]]
            np.savez_compressed(dest / (path.stem + '.npz'), image=hdul['IMAGE'].data[cut],
                                flags=hdul['FLAGS'].data[cut], xy=positions-low,
                                wave=wave, mjd=header['MJD-BEG'])


def compare(cutouts):
    paths = sorted(cutouts.glob('*.npz'))
    if len(paths) != 4:
        raise ValueError('This comparison requires all four detector-5 exposures of pointing 0225')
    results = {}
    for number, (name, _, _) in enumerate(SOURCES):
        pix, used, all_flux, wavelengths, times, records = [], [], [], [], [], []
        for path in paths:
            with np.load(path) as data:
                im, flags = data['image'], data['flags'].astype(np.uint32)
                if im.shape != flags.shape:
                    raise ValueError('Image and quality mask must match')
                cx, cy = data['xy'][number]
                y, x = np.indices(im.shape)
                dist = np.hypot(x-cx, y-cy)
                good = (flags & ~SOURCE_BIT) == 0
                annulus = (dist >= 4) & (dist <= 7) & good & np.isfinite(im)
                bg = np.median(im[annulus])
                aperture = dist <= 2.3
                flux = float(np.sum(im[aperture] - bg))
                x0, y0 = np.rint([cx, cy]).astype(int)
                iy, ix = np.mgrid[-2:3, -2:3]
                iy, ix = iy + y0, ix + x0
                pix.append(im[iy, ix] - bg)
                used.append(good[iy, ix] & aperture[iy, ix] & np.isfinite(im[iy, ix]))
                all_flux.append(flux)
                wavelengths.append(float(data['wave'][number, 0]))
                times.append(float(data['mjd']))
                records.append({'fits': path.stem + '.fits', 'utc': Time(times[-1], format='mjd').isot,
                                'wavelength_um': wavelengths[-1], 'bandpass_um': float(data['wave'][number, 1]),
                                'sum_mjy_sr': flux, 'source_pixels': int(aperture.sum()),
                                'warning_pixels': int((aperture & ~good).sum()),
                                'nearest_pixel_quality_hex': hex(int(flags[y0, x0]))})
        common = np.all(used, axis=0)
        flux_clean = np.sum(np.array(pix)[:, common], axis=1)
        for row, cleaned in zip(records, flux_clean):
            row['fixed_mask_sum_mjy_sr'] = float(cleaned)
        results[name] = {'minutes_first_to_last': (times[-1]-times[0])*1440,
                         'ratio_to_mean_first_two': all_flux[-1] / np.mean(all_flux[:2]),
                         'fixed_mask_ratio_to_mean_first_two': float(flux_clean[-1] / np.mean(flux_clean[:2])),
                         'fixed_mask_pixels': int(common.sum()), 'exposures': records}
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fits-cache', type=Path, help='Verify/download the native files and regenerate cutouts')
    parser.add_argument('--cutouts', type=Path, default=ROOT/'cutouts')
    args = parser.parse_args()
    if args.fits_cache:
        extract(args.fits_cache, args.cutouts)
    print(json.dumps(compare(args.cutouts), indent=2))
