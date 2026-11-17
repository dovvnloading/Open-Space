"""Measure F182M uptake of off-grid artificial satellites in the three real wide mosaics."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree

from wide_survey import EPOCH, projected, rank, visit_inputs
from wide_mosaic import block_stats, robust_map
from search import file_hash
from geometry import ARCSEC


FLUX = np.array([0., .08, .16, .20, .24, .28, .32, .48, .64])
PHASE = np.deg2rad(np.arange(0., 360., .25))
RADII = np.arange(200000., 650001., 500.)
STEP_S = 107.37
SEED = 27680465


def middle_paths(visit, radius, phase):
    return projected(visit, radius, np.atleast_1d(phase)) @ visit['M'] + visit['center']


def flux_template(visit, a, phase, bounds, model):
    """Conserve admitted encircled energy at arbitrary subpixel centers and trail phases."""
    x0, x1, y0, y1 = bounds
    result = np.zeros((y1-y0, x1-x0), 'f4')
    radius, density = model
    sub = (np.arange(12) + .5)/12 - .5
    for dt in (np.arange(5) - 2)*STEP_S/5:
        trail = middle_paths({**visit, 'time': visit['time'] + dt}, a, phase)[0]
        for x, y in trail:
            xx, yy = int(np.rint(x)), int(np.rint(y))
            xgrid = (np.arange(xx-4, xx+5) - x)[:, None] + sub
            ygrid = (np.arange(yy-4, yy+5) - y)[:, None] + sub
            pix = float(visit['pixel_scale'])
            r = np.hypot(ygrid.ravel()[:, None], xgrid.ravel()[None, :])*pix
            v = density[np.searchsorted(radius[1:], r)]*(pix/len(sub))**2
            response = v.reshape(9, 12, 9, 12).sum(axis=(1, 3))/25
            result[yy-y0-4:yy-y0+5, xx-x0-4:xx-x0+5] += response
    # SCI is MJy/sr. An input microjansky corresponds to 1e-12 MJy.
    return result * (1e-12/(float(visit['pixel_scale'])/ARCSEC)**2), float(result.sum())


def science_inputs(maps, folder):
    visits = []
    for m in sorted(maps):
        v = visit_inputs(m)
        path = folder / str(v['source_filename'])
        with fits.open(path) as h:
            w = h['WHT'].data
            if h[0].header['RESWHT'].strip() != 'exptime' or h[0].header['NDRIZ'] != 5:
                raise ValueError('This experiment requires five exposure-time-weighted dithers')
            v['sci'], v['valid'] = block_stats(h['SCI'].data, w, 3)
            v['con'] = np.memmap(path.with_suffix('.con.bin'), mode='r', dtype='>i4', shape=w.shape)
        visits.append(v)
    return visits


def context_support(v, path):
    good = map_coordinates(v['good'].astype(np.uint8), path.T[::-1], order=0, mode='constant')
    if not good.all():
        return False
    for x, y in path:
        # Full five-dither coverage in the native PSF core, and at least 95%
        # of the 0.83-arcsecond square containing each planted stamp.
        x, y = np.rint(np.array([x,y])*3+1).astype(int)
        c = v['con'][y-13:y+14,x-13:x+14]
        if c.shape != (27,27):
            return False
        full = np.ones_like(c, dtype=bool)
        for group in range(5):
            full &= ((c >> (4*group)) & 15) != 0
        if not full[11:16,11:16].all() or full.mean() < .95:
            return False
    return True


def draw_plants(n_per_bin):
    rng = np.random.default_rng(SEED)
    a = np.concatenate([rng.uniform(lo,lo+50000,n_per_bin) for lo in range(200000,650000,50000)])
    return list(zip(a, rng.uniform(0,2*np.pi,len(a))))


def nearby_templates(v, plants):
    # Complete subset: success below requires the trial and input midpoints to
    # agree within 2 binned pixels in all three visits, including visit one.
    mids = np.vstack([middle_paths(v,a,PHASE)[:,2,:] for a in RADII])
    tree = cKDTree(mids)
    return [tree.query_ball_point(middle_paths(v,a,ph)[0,2], 2.) for a,ph in plants]


def interval_bounds(v, true):
    x,y = np.mean(true, axis=0)
    # Tile edges align with the original 48x48 local MAD grid. A three-tile
    # margin supplies all Gaussian and per-tile neighbors for the center.
    tx, ty = int(x)//48, int(y)//48
    return (max(0,(tx-3)*48), min(v['z'].shape[1],(tx+4)*48),
            max(0,(ty-3)*48), min(v['z'].shape[0],(ty+4)*48))


def sample(image, paths):
    return map_coordinates(image, paths.transpose(2,0,1)[::-1], order=1,
                           prefilter=False, mode='constant', cval=0.)


def run(args):
    visits = science_inputs(args.maps, args.mosaics)
    ee = np.genfromtxt(args.ee, names=True)
    radius = np.r_[0,ee['aper_radius']]
    model = (radius, np.diff(np.r_[0,ee['F182M']])/(np.pi*np.diff(radius**2)))
    plants = draw_plants(args.per_bin)
    neighbors = nearby_templates(visits[0], plants)
    args.out.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for i, ((a,ph), idx) in enumerate(zip(plants, neighbors)):
        truths = [middle_paths(v,a,ph)[0] for v in visits]
        ok = np.array([context_support(v,t) for v,t in zip(visits,truths)])
        row = dict(plant=i, a_km=a, phase_deg=float(np.rad2deg(ph)),
                   full_context_visits=int(ok.sum()), eligible=False, compatible_count=0)
        if ok.all() and idx:
            ap, pb = np.divmod(idx,len(PHASE))
            preds = [np.stack([middle_paths(v,RADII[j],PHASE[k])[0] for j,k in zip(ap,pb)]) for v in visits]
            match = np.logical_and.reduce([np.linalg.norm(p-t,axis=2).max(1)<=2. for p,t in zip(preds,truths)])
            # Require complete bilinear support for the original search and all
            # 15 predicted positions, including the third, held-out visit.
            for v,p in zip(visits,preds):
                match &= sample(v['good'].astype('f4'),p).min(1)>.9999
            idx = np.array(idx)[match]
            row['compatible_count'] = len(idx)
            if len(idx):
                scores = []
                energy = []
                worst = 0.
                for v,t,p in zip(visits,truths,preds):
                    p = p[match]
                    bounds = interval_bounds(v,t)
                    x0,x1,y0,y1 = bounds
                    image = v['sci'][y0:y1,x0:x1]
                    valid = v['valid'][y0:y1,x0:x1]
                    added, enc = flux_template(v,a,ph,bounds,model)
                    energy.append(enc)
                    relative = p-[x0,y0]
                    values = []
                    original = sample(v['z'],p)
                    for flux in FLUX:
                        z,_ = robust_map(image+added*flux,valid)
                        points = sample(z,relative)
                        if flux == 0:
                            worst = max(worst, np.max(np.abs(original-points)))
                        values.append(points.mean(1))
                    scores.append(values)
                row['zero_flux_map_max_abs_error'] = worst
                row['mean_stamp_fraction'] = np.mean(energy)
                if worst < .05:
                    row['eligible'] = True
                    for f, z in zip(FLUX,np.stack(scores,axis=-1)):
                        statistic,_ = rank(z,np.ones_like(z,dtype=bool))
                        row[f'rank_{f:g}uJy'] = float(statistic.max())
        all_rows.append(row)
        if i%30==0:
            print(i,'eligible',sum(r['eligible'] for r in all_rows),flush=True)
    keys = ['plant','a_km','phase_deg','full_context_visits','eligible','compatible_count',
            'zero_flux_map_max_abs_error','mean_stamp_fraction',*[f'rank_{f:g}uJy' for f in FLUX]]
    with (args.out/'plants.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(all_rows)
    rows = [r for r in all_rows if r['eligible']]
    summary = dict(plant_count=len(plants),eligible=len(rows),seed=SEED,epoch_tdb_jd=EPOCH,
                   on_sky_unit='microJy total F182M flux density', psf='STScI EE radial density, 12x subpixel integration',
                   source_ee_sha256=file_hash(args.ee),maps_sha256=[file_hash(m) for m in sorted(args.maps)],
                   success='same on-sky grid, all-visit matching within 2 binned pixels, two visit scores >= 5, pooled >= 8')
    summary['counts']={f'{f:g}':sum(r[f'rank_{f:g}uJy']>=1 for r in rows) for f in FLUX}
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary['eligible'],summary['counts'],flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mosaics',type=Path);p.add_argument('ee',type=Path);p.add_argument('out',type=Path)
    p.add_argument('maps',nargs='+');p.add_argument('--per-bin',type=int,default=64)
    run(p.parse_args())
