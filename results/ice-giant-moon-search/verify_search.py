"""Deterministic end-to-end software check. All photons and orbits here are synthetic."""
import csv
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np

import search


def bank(offsets, x0, y0, time, exposure_id, label):
    # Two visits see the same curved path at different phases.
    t = np.linspace(-1.5, 1.5, 13)
    x = x0[None, :] + offsets[:, None]
    y = np.broadcast_to(y0[None, :], x.shape).copy()
    return dict(x=x, y=y, x_path=x[:, :, None]+t,
                y_path=y[:, :, None]+0.12*np.sin(time[None, :, None]+t),
                exposure_id=exposure_id, orbit_id=np.array([f'{label}_{i}' for i in range(len(x))]))


def plant(image, tracks, row, kernel, fluxes):
    for i, flux in enumerate(fluxes):
        y, x, stamp = search.trail_kernel(tracks['x_path'][row, i], tracks['y_path'][row, i], kernel)
        cut, cropped = search.clipped_stamp(image[i].shape, y, x, stamp)
        image[i][cut] += cropped*flux


def main():
    rng = np.random.default_rng(1471)
    n, side = 6, 100
    yy, xx = np.indices((side, side))
    radius = np.hypot(yy-50, xx-50)
    image = rng.normal(0, 1, (n, side, side)) + 12*np.exp(-(radius/20)**2)
    visit = np.repeat(['synthetic_A', 'synthetic_B'], 3)
    ids = np.array([f'SYNTH{i}' for i in range(n)])
    t = np.arange(n)
    xbase = np.array([45, 46, 47, 43, 44, 45])
    ybase = np.array([23, 24, 25, 74, 75, 76])
    search_bank = bank(np.array([-19., -10., 0., 10., 22.]), xbase, ybase, t, ids, 'grid')
    plants = bank(np.array([-18.8, -9.8]), xbase, ybase, t, ids, 'plant')
    kernel = search.gaussian(1.2)
    plant(image, search_bank, 2, kernel, [90]*n)      # consistent moving source
    plant(image, search_bank, 3, kernel, [160]*3+[0]*3)  # single-visit contaminant
    plant(image, search_bank, 4, kernel, [240, 0, 0, 240, 0, 0])  # cosmic rays
    valid = np.broadcast_to(radius > 13, image.shape).copy()
    cube = dict(image=image, variance=np.ones_like(image), valid=valid, visit=visit, center_xy=[50, 50],
                mjd=61000+np.array([0, .001, .002, 7, 7.001, 7.002]), exptime_s=np.full(n, 60),
                exposure_id=ids, psf=kernel, flux_unit='synthetic counts/pixel', band='synthetic')
    with TemporaryDirectory() as tmp:
        path = Path(tmp)
        np.savez(path/'cube.npz', **cube)
        np.savez(path/'tracks.npz', **search_bank)
        np.savez(path/'plants.npz', **plants)
        subprocess.run([sys.executable, str(Path(__file__).with_name('search.py')),
                        str(path/'cube.npz'), str(path/'tracks.npz'), str(path/'out'),
                        '--plants', str(path/'plants.npz'), '--flux', '0.001', '100',
                        '--n-plants', '2'], check=True)
        report = json.loads((path/'out/report.json').read_text())
        with open(path/'out/passing_templates.csv') as f:
            hits = list(csv.DictReader(f))
        assert [h['orbit_id'] for h in hits] == ['grid_2'], hits
        assert np.isclose(float(hits[0]['flux']), 90, rtol=0.20), hits
        assert report['negative_passing'] == 0, report
        assert [g['recovered'] for g in report['plant_recovery']] == [0, 2], report
        assert report['timing_model'] == 'exposure-integrated'
        # Check a fully masked second visit cannot produce two-visit evidence.
        cube['valid'][3:] = False
        maps = search.prepare(cube, [kernel]*n, search.Cuts())
        masked = search.scan(maps, search_bank, visit, search.Cuts())
        assert not masked['eligible'].any()
        print('Synthetic only: 1/1 two-visit source, 0/2 contaminant tracks, '
              '0/2 faint plants and 2/2 bright off-grid plants; masked visit rejected.')


if __name__ == '__main__':
    main()
