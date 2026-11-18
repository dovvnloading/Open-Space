"""Post-hoc sensitivity test frozen in protocol_off_control.md."""
import json
import numpy as np
import pandas as pd

from search import CONTROLS, DM, FIELDS, OUT, RADII, STEP, FieldGrid, catalogue, peaks, weights

N_PERMUTATIONS = 199
SEED = 20261118
EXCLUSION_ARCMIN = 26.0


def secondary_grid(field):
    grid = FieldGrid(field)
    ra, dec = grid.coords(grid.y, grid.x)
    for cra, cdec, _ in CONTROLS.values():
        grid.valid &= (grid.dist(ra, dec, cra, cdec) > EXCLUSION_ARCMIN)[None, :, :]
    return grid


def evaluate_field(field, rng):
    grid = secondary_grid(field)
    stars = catalogue(field)
    pixels, inside = grid.pixels(stars.ra.to_numpy(), stars.dec.to_numpy())
    raw_weights = weights(*(stars[c].to_numpy()[inside] for c in ['g', 'r', 'gerr', 'rerr']))
    observed = grid.scores(*grid.moments(pixels, raw_weights))
    observed_peaks = [dict(field=field, **row) for row in peaks(grid, observed)]
    geom = [
        dict(field=field, radius_arcmin=radius, n_centres=valid.sum(),
             grid_area_sqdeg=valid.sum()*STEP**2)
        for radius, valid in zip(RADII, grid.valid)
    ]
    cell = (pixels // grid.nx // 10)*100 + pixels % grid.nx // 10
    strata = [np.where(cell == key)[0] for key in np.unique(cell)]
    perm = np.arange(len(pixels))
    maxima = []
    for _ in range(N_PERMUTATIONS):
        for group in strata:
            perm[group] = rng.permutation(group)
        maxima.append(grid.scores(*grid.moments(pixels, raw_weights[:, perm])).max())
    return observed_peaks, geom, maxima


def main():
    rng = np.random.default_rng(SEED)
    observed, geometry, maxima = [], [], {}
    for field in FIELDS:
        field_peaks, field_geom, maxima[field] = evaluate_field(field, rng)
        observed.extend(field_peaks)
        geometry.extend(field_geom)
    null = pd.DataFrame(maxima)
    null['global'] = null.max(axis=1)
    scores = null['global'].to_numpy()
    allowed_exceedances = int(0.05*(N_PERMUTATIONS+1)) - 1
    threshold = np.sort(scores)[-(allowed_exceedances+1)]
    frame = pd.DataFrame(observed)
    frame['global_p'] = [(1+(scores >= score).sum()) / (N_PERMUTATIONS+1) for score in frame.score]
    frame['retained'] = frame.global_p <= 0.05
    for field in FIELDS:
        subset = frame.field == field
        rows = frame.loc[subset]
        frame.loc[subset, 'known_control_distance_arcmin'] = FieldGrid(field).dist(
            rows.ra, rows.dec, *CONTROLS[field][:2])
    frame.to_csv(OUT/'off_control_peaks.csv', index=False, float_format='%.8f')
    frame[frame.retained].to_csv(OUT/'off_control_retained.csv', index=False, float_format='%.8f')
    pd.DataFrame(geometry).to_csv(OUT/'off_control_geometry.csv', index=False, float_format='%.8f')
    null.to_csv(OUT/'off_control_field_maxima_null.csv', index=False, float_format='%.12g')
    winning_fields = null[list(FIELDS)].idxmax(axis=1).value_counts().to_dict()
    (OUT/'off_control_calibration.json').write_text(json.dumps(dict(
        seed=SEED, permutations=N_PERMUTATIONS, excluded_centres_within_arcmin=EXCLUSION_ARCMIN,
        strict_threshold=float(threshold), global_null_maxima_by_winning_field=winning_fields,
        n_retained=int(frame.retained.sum()), n_moduli=len(DM),
        p_method='(1 + count of paired global maxima >= observed)/200'), indent=2)+'\n')


if __name__ == '__main__':
    main()
