"""Orbit-matched filtering and blind plant recovery on registered single-band cubes."""
import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import scipy
from scipy import ndimage, signal
from scipy.stats import binomtest


@dataclass
class Cuts:
    visit_snr: float = 5.0
    total_snr: float = 7.0
    dropout_snr: float = 2.0
    min_exposures: int = 2
    glare_sigma: float = 12.0
    psf_support: float = 0.9


def file_hash(path):
    with open(path, 'rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load_npz(path):
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def strvec(values):
    return np.asarray(values).astype(str)


def gaussian(sigma):
    if not np.isfinite(sigma) or sigma < 0.3:
        raise ValueError('PSF sigma must be finite and at least 0.3 pixels')
    radius = int(np.ceil(5 * sigma))
    dy, dx = np.mgrid[-radius:radius + 1, -radius:radius + 1]
    kernel = np.exp(-(dx**2 + dy**2) / (2 * sigma**2))
    return kernel / kernel.sum()


def validate_cube(cube):
    image = cube['image']
    if image.ndim != 3 or min(image.shape[1:]) < 12:
        raise ValueError('image must have axes exposure,y,x with spatial coverage')
    n = len(image)
    if any(cube[key].shape != image.shape for key in ['valid', 'variance']):
        raise ValueError('variance and valid must have the image shape')
    for key in ['mjd', 'exposure_id', 'visit']:
        if cube[key].shape != (n,):
            raise ValueError(f'{key} must contain one entry per exposure')
    if not np.isfinite(cube['mjd']).all() or np.unique(strvec(cube['exposure_id'])).size != n:
        raise ValueError('finite epochs and unique exposure IDs required')
    if np.any(strvec(cube['visit']) == '') or np.unique(strvec(cube['visit'])).size < 2:
        raise ValueError('two distinct archive visit IDs required')
    for key in ['band', 'flux_unit']:
        if cube[key].size == 0 or np.unique(strvec(cube[key])).size != 1 or not strvec(cube[key]).flat[0]:
            raise ValueError(f'supply one nonempty {key} per run')
    if cube['valid'].dtype.kind != 'b':
        raise ValueError('valid must be boolean with all geometry and quality masks applied')
    if 'ring_mask' in cube and np.any(cube['valid'] & cube['ring_mask']):
        raise ValueError('ring_mask is not excluded in valid')
    centers = np.broadcast_to(cube['center_xy'], (n, 2))
    if not np.isfinite(centers).all():
        raise ValueError('finite registered planet centers required')
    if 'psf' in cube:
        psfs = np.asarray(cube['psf'], dtype=float)
        if psfs.ndim == 2:
            psfs = np.broadcast_to(psfs, (n, *psfs.shape))
        if psfs.ndim != 3 or len(psfs) != n or any(s % 2 == 0 for s in psfs.shape[1:]):
            raise ValueError('supply odd-sized pixel-centered PSF stamps')
        if not np.isfinite(psfs).all() or np.any(psfs < 0) or np.any(psfs.sum(axis=(1, 2)) <= 0):
            raise ValueError('PSF must contain finite nonnegative source flux')
        return [p / p.sum() for p in psfs]
    return [gaussian(s) for s in np.broadcast_to(cube['psf_sigma_px'], (n,))]


def validate_tracks(tracks, cube):
    x, y = tracks['x'], tracks['y']
    if x.ndim != 2 or x.shape != y.shape or x.shape[1] != len(cube['image']) or len(x) == 0:
        raise ValueError('x,y must have shape (nonempty orbit bank, all exposures)')
    if not np.array_equal(strvec(tracks['exposure_id']), strvec(cube['exposure_id'])):
        raise ValueError('exposure order differs between cube and orbit bank')
    if ('x_path' in tracks) != ('y_path' in tracks):
        raise ValueError('both x_path and y_path required for exposure-integrated templates')
    if 'x_path' in tracks:
        paths = tracks['x_path'], tracks['y_path']
        if (paths[0].shape != paths[1].shape or paths[0].ndim != 3
                or paths[0].shape[:2] != x.shape or paths[0].shape[-1] < 2):
            raise ValueError('exposure paths require shape (orbit, exposure, >=2 time samples)')
        for start in range(0, len(x), 128):
            px, py = (arr[start:start + 128] for arr in paths)
            spacing = np.hypot(np.diff(px), np.diff(py))
            if np.any(spacing > 0.75):
                raise ValueError('increase orbit time samples: maximum adjacent displacement is 0.75 pixel')
    elif 'exptime_s' in cube:
        raise ValueError('finite exposures require x_path and y_path; integrate all search and plant tracks')
    ids = strvec(tracks['orbit_id'])
    if ids.shape != (len(x),) or np.unique(ids).size != len(x):
        raise ValueError('orbit_id must uniquely identify each template')


def wide_background(data, mask, width):
    def smooth(use):
        den = ndimage.gaussian_filter(use.astype(float), width, mode='constant')
        num = ndimage.gaussian_filter(np.where(use, data, 0.0), width, mode='constant')
        return np.divide(num, den, out=np.zeros_like(num), where=den > 0.01), den

    rough, den = smooth(mask)
    residual = (data - rough)[mask & (den > 0.1)]
    if not residual.size:
        return rough, den
    median = np.median(residual)
    mad = 1.4826 * np.median(np.abs(residual - median))
    if mad <= 0:
        return rough, den
    return smooth(mask & (np.abs(data - rough - median) <= 4.5 * mad))


def filtered_frame(data, variance, valid, kernel, center, cuts):
    good = valid & np.isfinite(data) & np.isfinite(variance) & (variance > 0)
    bg, support = wide_background(data, good, cuts.glare_sigma)
    good &= support > 0.1
    residual = np.where(good, data - bg, 0.0)
    invvar = np.divide(1.0, variance, where=good, out=np.zeros(data.shape))
    # Reverse the stamp to compute cross-correlation for asymmetric PSFs.
    def corr(array, stamp):
        return signal.fftconvolve(array, stamp[::-1, ::-1], mode='same')

    numerator = corr(residual * invvar, kernel)
    information = np.maximum(corr(invvar, kernel**2), 0.0)
    pixel_good = good.copy()
    good &= (corr((invvar > 0).astype(float), kernel) >= cuts.psf_support) & (information > 0)
    snr = np.divide(numerator, np.sqrt(information), where=good, out=np.zeros(data.shape))
    yy, xx = np.indices(data.shape)
    annulus = (np.hypot(xx - center[0], yy - center[1]) / 12).astype(int)
    # Broad radial correction for underestimated variance/correlated noise.
    # Resampling still correlates neighboring trial orbits. Scores are diagnostics.
    raw_num = np.zeros_like(invvar)
    raw_info = np.zeros_like(invvar)
    for zone in np.unique(annulus[good]):
        use = (annulus == zone) & good
        sky = snr[use]
        if sky.size < 64:
            sky = snr[good]
        median = np.median(sky)
        scale = max(1.0, 1.4826 * np.median(np.abs(sky - median)))
        numerator[use] = (numerator[use] - median * np.sqrt(information[use])) / scale**2
        information[use] /= scale**2
        raw = (annulus == zone) & pixel_good
        raw_info[raw] = invvar[raw] / scale**2
        raw_num[raw] = (residual[raw] - np.median(residual[raw])) * raw_info[raw]
    return (np.where(good, numerator, 0), np.where(good, information, 0), good.astype(float),
            raw_num, raw_info, kernel)


def prepare(cube, kernels, cuts):
    centers = np.broadcast_to(cube['center_xy'], (len(kernels), 2))
    return [filtered_frame(*args, cuts) for args in zip(
        cube['image'], cube['variance'], cube['valid'], kernels, centers, strict=True)]


def sample(maps, x, y):
    ntrial, nexp = x.shape
    nums, info = np.zeros((2, ntrial, nexp))
    for t, (num, weight, valid, *_) in enumerate(maps):
        finite = np.isfinite(x[:, t]) & np.isfinite(y[:, t])
        coords = np.where(finite, [y[:, t], x[:, t]], -1)
        def interpolate(image):
            return ndimage.map_coordinates(image, coords, order=1, mode='constant', cval=0)
        supported = interpolate(valid) >= 1 - 1e-6
        nums[:, t] = np.where(supported, interpolate(num), 0)
        info[:, t] = np.where(supported, interpolate(weight), 0)
    return nums, info


def trail_kernel(px, py, psf):
    """Integrate a unit PSF along uniformly timed path samples (endpoint half-weights)."""
    left, top = int(np.floor(px.min())), int(np.floor(py.min()))
    x, y = px-left, py-top
    ix, iy = np.floor(x).astype(int), np.floor(y).astype(int)
    occupancy = np.zeros((iy.max() + 2, ix.max() + 2))
    time_weights = np.ones(len(px))
    time_weights[[0, -1]] = 0.5
    time_weights /= time_weights.sum()
    for dx, wx in [(0, 1 - (x-ix)), (1, x-ix)]:
        for dy, wy in [(0, 1 - (y-iy)), (1, y-iy)]:
            np.add.at(occupancy, (iy+dy, ix+dx), wx * wy * time_weights)
    kernel = signal.convolve(occupancy, psf, mode='full', method='auto')
    return top-psf.shape[0]//2, left-psf.shape[1]//2, kernel


def clipped_stamp(shape, top, left, stamp):
    ny, nx = stamp.shape
    y0, x0 = max(0, top), max(0, left)
    y1, x1 = min(shape[0], top+ny), min(shape[1], left+nx)
    if y1 <= y0 or x1 <= x0:
        return None
    return (np.s_[y0:y1, x0:x1], stamp[y0-top:y1-top, x0-left:x1-left])


def sample_paths(maps, px, py, cuts):
    num, info = np.zeros((2, *px.shape[:2]))
    for t, (_, _, center_ok, image_weighted, weight, psf) in enumerate(maps):
        ny, nx = weight.shape
        for k in range(len(px)):
            x, y = px[k, t], py[k, t]
            if not (np.isfinite(x+y).all() and np.all((x >= 0) & (y >= 0) & (x <= nx-1) & (y <= ny-1))):
                continue
            centers = ndimage.map_coordinates(center_ok, [y, x], order=1)
            if centers.min() < 1-1e-6:
                continue
            cut, stamp = clipped_stamp(weight.shape, *trail_kernel(x, y, psf))
            w = weight[cut]
            if stamp[w > 0].sum() < cuts.psf_support:
                continue
            num[k, t] = np.sum(stamp * image_weighted[cut])
            info[k, t] = np.sum(stamp**2 * w)
    return num, info


def snr(num, info):
    return np.divide(num, np.sqrt(info), out=np.zeros_like(num, dtype=float), where=info > 0)


def score(nums, info, visits, cuts):
    visits = strvec(visits)
    per_visit, dropout, counts = [], [], []
    pooled_num, pooled_info = np.zeros((2, len(nums)))
    for visit in np.unique(visits):
        take = visits == visit
        num, weight = nums[:, take], info[:, take]
        summed_num, summed_weight = num.sum(1), weight.sum(1)
        count = (weight > 0).sum(1)
        usable = count >= cuts.min_exposures
        view = snr(summed_num, summed_weight)
        # Worst supported leave-one-exposure-out score in this visit.
        remaining = snr(summed_num[:, None] - num,
                        np.maximum(summed_weight[:, None] - weight, 0))
        worst = np.min(np.where(weight > 0, remaining, np.inf), axis=1)
        per_visit.append(np.where(usable, view, 0))
        dropout.append(np.where(usable, worst, 0))
        counts.append(count)
        pooled_num += summed_num * usable
        pooled_info += summed_weight * usable
    per_visit, dropout, counts = map(np.asarray, (per_visit, dropout, counts))
    visit_score = np.minimum(per_visit / cuts.visit_snr, dropout / cuts.dropout_snr)
    second = np.partition(visit_score, -2, axis=0)[-2]
    total = snr(pooled_num, pooled_info)
    rank = np.minimum(total / cuts.total_snr, second)
    flux = np.divide(pooled_num, pooled_info, out=np.zeros_like(pooled_num), where=pooled_info > 0)
    return dict(rank=rank, snr=total, flux=flux, visit_snr=per_visit.T,
                dropout=dropout.T, n=counts.T, eligible=(counts >= cuts.min_exposures).sum(0) >= 2)


def scan(maps, tracks, visits, cuts, negative=False, chunk=4096):
    fields = {}
    for start in range(0, len(tracks['x']), chunk):
        stop = start + chunk
        if 'x_path' in tracks:
            num, weight = sample_paths(maps, tracks['x_path'][start:stop],
                                       tracks['y_path'][start:stop], cuts)
        else:
            num, weight = sample(maps, tracks['x'][start:stop], tracks['y'][start:stop])
        for key, value in score(-num if negative else num, weight, visits, cuts).items():
            fields.setdefault(key, []).append(value)
    return {key: np.concatenate(parts) for key, parts in fields.items()}


def add_stamp(image, x, y, psf, flux):
    if not np.isfinite(x + y):
        return
    ix, iy = int(np.rint(x)), int(np.rint(y))
    hy, hx = np.asarray(psf.shape) // 2
    y0, x0 = iy - hy, ix - hx
    y1, x1 = max(y0, 0), max(x0, 0)
    y2, x2 = min(y0 + psf.shape[0], image.shape[0]), min(x0 + psf.shape[1], image.shape[1])
    if y1 >= y2 or x1 >= x2:
        return
    shifted = ndimage.shift(psf, (y - iy, x - ix), order=1, mode='constant')
    shifted /= shifted.sum()
    image[y1:y2, x1:x2] += flux * shifted[y1-y0:y2-y0, x1-x0:x2-x0]


def matched(index, plant_x, plant_y, tracks, samples, visits, cuts, radius):
    # Require positional agreement in two visits which independently passed.
    for k in index:
        passed = ((samples['visit_snr'][k] >= cuts.visit_snr)
                  & (samples['dropout'][k] >= cuts.dropout_snr))
        close = 0
        for j, visit in enumerate(np.unique(visits)):
            if not passed[j]:
                continue
            use = visits == visit
            delta = np.hypot(plant_x[use] - tracks['x'][k, use], plant_y[use] - tracks['y'][k, use])
            # Conservative when an injected ephemeris has gaps.
            if np.all(delta <= radius):
                close += 1
        if close >= 2:
            return int(k)
    return None


def plants(cube, kernels, tracks, injections, cuts, baseline, levels, n_per_level, seed, radius):
    if len(injections['x']) < n_per_level:
        raise ValueError('off-grid bank must contain at least --n-plants distinct orbits')
    if np.intersect1d(strvec(injections['orbit_id']), strvec(tracks['orbit_id'])).size:
        raise ValueError('plant-bank IDs must be distinct from the search bank IDs')
    rng = np.random.default_rng(seed)
    visits = strvec(cube['visit'])
    base_hits = np.flatnonzero(baseline['rank'] >= 1)
    for flux in levels:
        for i in rng.choice(len(injections['x']), n_per_level, replace=False):
            x, y = injections['x'][i], injections['y'][i]
            if np.any(np.all((tracks['x'] == x) & (tracks['y'] == y), axis=1)):
                raise ValueError('plant path lies exactly on search grid: generate off-grid trials')
            row = dict(flux=float(flux), orbit_id=str(injections['orbit_id'][i]), recovered=False,
                       status='missed', recovered_orbit_id='')
            if 'a_km' in injections:
                row['a_km'] = float(injections['a_km'][i])
            if matched(base_hits, x, y, tracks, baseline, visits, cuts, radius) is not None:
                row['status'] = 'preexisting_peak'
                yield row
                continue
            image = cube['image'].astype(float, copy=True)
            if 'x_path' in injections:
                for j, (frame, psf) in enumerate(zip(image, kernels, strict=True)):
                    px, py = injections['x_path'][i, j], injections['y_path'][i, j]
                    if not np.isfinite(px+py).all():
                        continue
                    # Out-of-field tracks cannot add light beyond the PSF support.
                    if (px.max() < -psf.shape[1] or py.max() < -psf.shape[0]
                            or px.min() > frame.shape[1]+psf.shape[1]
                            or py.min() > frame.shape[0]+psf.shape[0]):
                        continue
                    patch = clipped_stamp(frame.shape, *trail_kernel(px, py, psf))
                    if patch is not None:
                        cut, stamp = patch
                        frame[cut] += flux * stamp
            else:
                for frame, px, py, psf in zip(image, x, y, kernels, strict=True):
                    add_stamp(frame, px, py, psf, flux)
            maps = prepare({**cube, 'image': image}, kernels, cuts)
            result = scan(maps, tracks, visits, cuts)
            peaks = np.flatnonzero(result['rank'] >= 1)
            peaks = peaks[np.argsort(-result['rank'][peaks])]
            k = matched(peaks, x, y, tracks, result, visits, cuts, radius)
            if k is not None:
                row.update(recovered=True, status='recovered', recovered_orbit_id=str(tracks['orbit_id'][k]))
            yield row


def summarize_plants(rows):
    groups = []
    for flux in sorted({r['flux'] for r in rows}):
        use = [r for r in rows if r['flux'] == flux]
        n = len(use)
        recovered = sum(r['recovered'] for r in use)
        ci = binomtest(recovered, n).proportion_ci()
        groups.append(dict(flux=flux, recovered=recovered, planted=n,
                           efficiency=recovered/n, exact_95_lo=float(ci.low),
                           exact_95_hi=float(ci.high),
                           preexisting=sum(r['status'] == 'preexisting_peak' for r in use)))
    return groups


def write_rows(path, rows, fieldnames):
    with open(path, 'x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cube', type=Path)
    p.add_argument('tracks', type=Path)
    p.add_argument('out', type=Path, help='new output directory')
    p.add_argument('--plants', type=Path, help='independently sampled off-grid orbit bank')
    p.add_argument('--flux', nargs='+', type=float, help='integrated PSF flux in cube units')
    p.add_argument('--n-plants', type=int, default=100, help='per flux per supplied orbital population')
    p.add_argument('--seed', type=int, default=23017)
    p.add_argument('--match-px', type=float, default=1.0)
    p.add_argument('--visit-snr', type=float, default=5.0)
    p.add_argument('--total-snr', type=float, default=7.0)
    p.add_argument('--dropout-snr', type=float, default=2.0)
    p.add_argument('--glare-sigma', type=float, default=12.0)
    args = p.parse_args()
    if bool(args.flux) != bool(args.plants):
        p.error('--plants and --flux must be supplied together')
    positive = [args.match_px, args.visit_snr, args.total_snr, args.dropout_snr,
                args.glare_sigma, *(args.flux or [])]
    if not np.isfinite(positive).all() or min(positive) <= 0 or args.n_plants < 1:
        p.error('thresholds, scales, fluxes and plant count must be positive and finite')
    if args.flux and len(set(args.flux)) != len(args.flux):
        p.error('duplicate injected flux level')
    cube, tracks = load_npz(args.cube), load_npz(args.tracks)
    kernels = validate_cube(cube)
    validate_tracks(tracks, cube)
    injections = load_npz(args.plants) if args.plants else None
    if injections is not None:
        validate_tracks(injections, cube)
        if ('x_path' in tracks) != ('x_path' in injections):
            p.error('plant and search bank must share the same exposure timing model')
    cuts = Cuts(visit_snr=args.visit_snr, total_snr=args.total_snr, dropout_snr=args.dropout_snr,
                glare_sigma=args.glare_sigma)
    maps = prepare(cube, kernels, cuts)
    visits = strvec(cube['visit'])
    measured = scan(maps, tracks, visits, cuts)
    negative = scan(maps, tracks, visits, cuts, negative=True)
    out = args.out
    out.mkdir(parents=True, exist_ok=False)
    report = dict(status='image screening only', cube_sha256=file_hash(args.cube),
                  tracks_sha256=file_hash(args.tracks), search_sha256=file_hash(Path(__file__)),
                  cuts=vars(cuts), visits=np.unique(visits).tolist(), n_exposures=len(visits),
                  band=str(strvec(cube['band']).flat[0]), flux_unit=str(strvec(cube['flux_unit']).flat[0]),
                  psf_source='supplied stamp' if 'psf' in cube else 'Gaussian approximation',
                  timing_model='exposure-integrated' if 'x_path' in tracks else 'midpoint approximation',
                  templates=len(tracks['x']), eligible_templates=int(measured['eligible'].sum()),
                  negative_passing=int((negative['rank'] >= 1).sum()),
                  negative_max_rank=float(negative['rank'].max()),
                  positive_passing=int((measured['rank'] >= 1).sum()),
                  positive_max_rank=float(measured['rank'].max()),
                  numpy=np.__version__, scipy=scipy.__version__, python=sys.version.split()[0])
    cols = ['orbit_id', 'rank', 'snr', 'flux', 'visit_snr', 'dropout_snr', 'n_exposures']
    def peak_rows():
        for k in np.flatnonzero(measured['rank'] >= 1):
            yield dict(zip(cols, [tracks['orbit_id'][k], measured['rank'][k], measured['snr'][k],
                                  measured['flux'][k], *[json.dumps(measured[f][k].tolist())
                                  for f in ['visit_snr', 'dropout', 'n']]], strict=True))
    write_rows(out/'passing_templates.csv', peak_rows(), cols)
    if injections is not None:
        report.update(plants_sha256=file_hash(args.plants), match_px=args.match_px,
                      seed=args.seed, flux_levels=args.flux, n_plants_per_level=args.n_plants)
        rows = []
        fields = ['flux', 'orbit_id', 'recovered', 'status', 'recovered_orbit_id']
        if 'a_km' in injections:
            fields.append('a_km')
        with open(out/'plants.csv', 'x', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for row in plants(cube, kernels, tracks, injections, cuts, measured, args.flux,
                              args.n_plants, args.seed, args.match_px):
                rows.append(row)
                writer.writerow(row)
                stream.flush()
        report['plant_recovery'] = summarize_plants(rows)
    (out/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
