"""Prepare an encircled-energy PSF and independent circular orbit banks for pilot search."""
import argparse
import hashlib
from pathlib import Path
import sys
from urllib.request import urlopen

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import geometry

EE_URL = ('https://jwst-docs.stsci.edu/files/216457060/216457083/1/'
          '1762453590926/Encircled_Energy_SW_ETCv2.txt')


def radial_psf(ee_file, scale, filter_name):
    table = np.genfromtxt(ee_file, names=True)
    radius = np.r_[0, table['aper_radius']]
    ee = np.r_[0, table[filter_name]]
    density = np.diff(ee) / (np.pi * np.diff(radius**2))
    oversample, half = 16, 6
    subpixels = (np.arange(oversample) + .5)/oversample-.5
    pix = (np.arange(-half, half+1)[:, None]+subpixels).ravel()*scale
    radial = np.hypot(pix[:, None], pix[None, :])
    intensity = density[np.searchsorted(radius[1:], radial)]*(scale/oversample)**2
    kernel = intensity.reshape(2*half+1, oversample, 2*half+1, oversample).sum(axis=(1, 3))
    return kernel


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cube', type=Path)
    p.add_argument('out', type=Path)
    p.add_argument('--ee', type=Path, help='local copy of STScI table')
    p.add_argument('--seed', type=int, default=42031)
    p.add_argument('--n', type=int, default=120)
    p.add_argument('--grid', type=float, nargs=4, default=[79500, 80500, 250, .75])
    args = p.parse_args()
    out = args.out
    out.mkdir(exist_ok=True, parents=True)
    ee = out/'Encircled_Energy_SW_ETCv2.txt'
    ee.write_bytes(args.ee.read_bytes() if args.ee else urlopen(EE_URL, timeout=60).read())
    with np.load(args.cube) as f:
        cube = {k: f[k] for k in f.files}
    band = str(cube['band']).split('/')[-1]
    if band == 'CLEAR':
        band = str(cube['band']).split('/')[-2]
    psf = radial_psf(ee, float(cube['scale_arcsec']), band)
    np.savez_compressed(out/'science.npz', **{**cube, 'psf': psf, 'psf_fraction': psf.sum()})
    r0, r1, dr, dphi = args.grid
    geometry.generate(args.cube, out/'tracks.npz', grid=[r0, r1, dr, dphi])
    rng = np.random.default_rng(args.seed)
    rows = np.zeros((args.n, 6))
    rows[:, 0] = rng.uniform(r0, r1, args.n)
    rows[:, -1] = rng.uniform(0, 360, args.n)
    np.savetxt(out/'injected_elements.csv', rows, delimiter=',', comments='',
               header='a_km,e,i_deg,node_deg,peri_deg,mean_deg')
    with np.load(out/'tracks.npz') as f:
        epoch = float(f['epoch_jd_tdb'])
    geometry.generate(args.cube, out/'plants.npz', elements=out/'injected_elements.csv',
                      model='circular', epoch=epoch)
    with np.load(out/'plants.npz') as f:
        bank = {k: f[k] for k in f.files}
    bank['orbit_id'] = np.array([f'plant_{i}' for i in range(args.n)])
    np.savez_compressed(out/'plants.npz', **bank)
    print(f'EE source SHA-256 {hashlib.sha256(ee.read_bytes()).hexdigest()}, '
          f'{band} model stamp fraction {psf.sum():.6f}')


if __name__ == '__main__':
    main()
