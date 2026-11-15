"""Reproduce local-background aperture checks and compact evidence tables."""
import csv
from pathlib import Path
import numpy as np
from scipy.stats import median_abs_deviation

Y, X = np.indices((97, 97))
TARGET = (49, 48)  # visit 52 peak within two pixels of the published coordinate


def aperture(image, x, y, radius=2.7):
    r2 = (X - x)**2 + (Y - y)**2
    core = r2 <= radius**2
    sky = (r2 > 16) & (r2 <= 36)
    return image[core].sum() - core.sum() * np.median(image[sky])


def z(value, sample):
    sigma = median_abs_deviation(sample, scale='normal', nan_policy='omit')
    return (value - np.nanmedian(sample)) / sigma


def measure(directory):
    rows = []
    for seq, band in [(1, 'F275W'), (2, 'F606W')]:
        files = sorted(directory.glob(f'id5o??0{seq}0_drc.fits.npz'))
        visits = np.array([int(p.name[4:6]) for p in files])
        if not np.array_equal(visits, np.arange(1, 54)):
            raise ValueError(f'Missing visits for {band}')
        images = np.array([np.load(p)['image'] for p in files])
        baseline = visits <= 41
        template = np.nanmedian(images[baseline], axis=0)
        differences = images - template
        fluxes = np.array([aperture(im, *TARGET) for im in differences])
        if band == "F606W":
            centers = []
            for im in differences[[42, 51]]:
                r2 = (X - TARGET[0])**2 + (Y - TARGET[1])**2
                sky = np.median(im[(r2 > 16) & (r2 <= 36)])
                weights = np.maximum(im - sky, 0) * (r2 < 9)
                centers.append([(weights * axis).sum() / weights.sum() for axis in (X, Y)])
            print(f"June/July centroid offset: {0.04 * np.linalg.norm(np.subtract(*centers)):.4f} arcsec")
        for visit, file, im, flux in zip(visits, files, differences, fluxes):
            controls = [aperture(im, x, y) for x in range(7, 91, 6) for y in range(7, 91, 6)
                        if (x-49)**2 + (y-48)**2 > 100]
            rows.append([file.name.removesuffix('.npz'), band, visit, float(np.load(file)['mjd']),
                         flux, z(flux, controls), z(flux, fluxes[baseline]), aperture(im, 49, 40)])
    with open('visits.csv', 'w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['product', 'band', 'visit', 'mjd_start', 'target_e_s', 'spatial_mad_z',
                         'baseline_mad_z', 'neighbor_e_s'])
        writer.writerows([[round(s, 5) if isinstance(s, (float, np.floating)) else s for s in row] for row in rows])
    rows = []
    for path in sorted(directory.glob('*_flc.npz')):
        with np.load(path) as data:
            image = data['image'] / float(data['exptime'])
            rows.append([path.stem+'.fits', str(data['filter']), round(float(data['mjd']), 5),
                         round(aperture(image, *TARGET), 5), round(aperture(image, 49, 40), 5)])
    with open('exposures.csv', 'w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['product', 'band', 'mjd_start', 'target_e_s', 'neighbor_e_s'])
        writer.writerows(rows)


if __name__ == '__main__':
    measure(Path('stamps'))
