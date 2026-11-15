"""Compare stationary and catalog-moving positions using the four frozen VLASS images."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates, spline_filter
from scipy.optimize import minimize
from scipy.signal import correlate

from replay import Image, predicted


def amplitude_map(image, origin, half_width):
    y, x = np.mgrid[-half_width:half_width + 1, -half_width:half_width + 1]
    x, y = x.ravel(), y.ravel()
    # The first coefficient gives the restoring-beam amplitude allowing a local plane.
    design = np.stack([image.kernel(x, y), np.ones_like(x), x, y], axis=1)
    weights = np.linalg.pinv(design)[0].reshape((2 * half_width + 1,) * 2)
    amplitude = correlate(np.nan_to_num(image.image), weights, mode='same')
    center = np.asarray(image.wcs.world_to_pixel_values(*origin))
    offsets = np.array([(a, b) for a in np.arange(-27, 28, 3.1)
                        for b in np.arange(-27, 28, 3.1) if 8 < np.hypot(a, b) < 27]).T
    sample = map_coordinates(amplitude, (center[:, None] + offsets)[::-1], order=3)
    background = np.median(sample)
    rms = 1.4826 * np.median(abs(sample - background))
    return spline_filter(amplitude, order=3), center, np.linalg.inv(image.jac), background, rms


def best_position(maps, path, center, signed_2020, radius=3):
    """Fit a common sky offset and independent amplitudes. Score uses one value per epoch."""
    def evaluate(offsets, amplitudes_only=False):
        offsets = np.atleast_2d(offsets)
        voltages, scores = [], []
        for epoch, (pixels, origin, inverse_jac, bg, noise) in enumerate(maps):
            sky = center + path[epoch] + offsets
            xy = origin[:, None] + inverse_jac @ sky.T
            amplitude = map_coordinates(pixels, xy[::-1], order=3, prefilter=False)
            voltages.append(amplitude)
            z = (amplitude - bg) / noise
            if epoch == 1 and signed_2020 == 'omit':
                z = np.zeros_like(z)
            elif epoch != 1 or signed_2020 == 'positive':
                z = np.maximum(z, 0)
            scores.append(z ** 2)
        return np.array(voltages) if amplitudes_only else np.sum(scores, axis=0)

    grid = np.arange(-radius, radius + .01, .3)
    dx, dy = np.meshgrid(grid, grid)
    offsets = np.column_stack([dx.ravel(), dy.ravel()])
    start = offsets[np.argmax(evaluate(offsets))]
    opt = minimize(lambda v: -evaluate(v)[0], start, method='Nelder-Mead',
                   bounds=[(-radius, radius)] * 2, options={'xatol': .003})
    p = opt.x
    return {'score': float(-opt.fun), 'offset_east_north_arcsec': np.round(p, 3).tolist(),
            'amplitude_mjy': np.round(evaluate(p, True).ravel(), 3).tolist()}


def compare(maps, delta, center, treatment, radius=3):
    moving = best_position(maps, delta, center, treatment, radius)
    fixed = best_position(maps, np.zeros_like(delta), center, treatment, radius)
    return {'moving': moving, 'stationary': fixed,
            'moving_minus_stationary_score': round(moving['score'] - fixed['score'], 3)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent
    parser.add_argument('--inputs', type=Path, default=root / 'inputs')
    parser.add_argument('--output', type=Path, default=root / 'track-fit.json')
    args = parser.parse_args()
    images = [Image(args.inputs / f'VLASS{k}.fits') for k in range(1, 5)]
    origin = np.asarray(predicted(images[0].date))
    path = (np.array([predicted(im.date) for im in images]) - origin) * 3600
    path[:, 0] *= np.cos(np.deg2rad(origin[1]))
    midpoint = path.mean(axis=0)
    delta = path - midpoint
    out = {'anchor': 'stationary offsets are relative to the mean predicted epoch position; '
                     'moving offsets are relative to the predicted track',
           'score': 'sum of squared background-subtracted beam amplitudes divided by '
                    'empirical per-epoch amplitude-map noise squared; descriptive only'}
    for n in (7, 9, 11):
        maps = [amplitude_map(im, origin, n) for im in images]
        out[f'{2*n+1}_pixel_window'] = compare(maps, delta, midpoint, 'signed')
        if n != 9:
            continue
        for treatment in ('omit', 'positive'):
            out[f'2020_{treatment}'] = compare(maps, delta, midpoint, treatment)
        out['wider_5_arcsec_astrometric_search'] = compare(maps, delta, midpoint, 'signed', radius=5)
        offsets = [(e, north) for e in np.arange(-16., 16.1, 3.8)
                   for north in np.arange(-16., 16.1, 3.8) if 7 < np.hypot(e, north) < 18]
        controls = [compare(maps, delta, midpoint + v, 'signed')['moving_minus_stationary_score']
                    for v in offsets]
        out['shifted_controls'] = {'count': len(controls), 'median': float(np.median(controls)),
                                   'min': min(controls), 'max': max(controls),
                                   '95th_percentile': float(np.percentile(controls, 95))}
    args.output.write_text(json.dumps(out, indent=2) + '\n')


if __name__ == '__main__':
    main()
