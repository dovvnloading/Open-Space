"""Predict NGC 6517 pulsar dispersion measures from angular position.

Run with NumPy and SciPy from any working directory. Input data and resulting JSON
are colocated with this script. All positions are east/north offsets in arcsec.
"""
import csv
import json
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parent
RA0 = (18 + 1 / 60 + 50.52 / 3600) * 15
DEC0 = -(8 + 57 / 60 + 31.6 / 3600)


def read(name):
    with (ROOT / name).open() as stream:
        return list(csv.DictReader(stream))


def sexagesimal(value):
    h, m, s = map(float, value.split(':'))
    return (-1 if value.startswith('-') else 1) * (abs(h) + m/60 + s/3600)


def offsets(rows):
    ra = np.array([15 * sexagesimal(row['ra_j2000']) for row in rows])
    dec = np.array([sexagesimal(row['dec_j2000']) for row in rows])
    return np.column_stack([(ra - RA0) * np.cos(np.radians(DEC0)), dec - DEC0]) * 3600


def design(x):
    return np.column_stack([np.ones(len(x)), x])


def rms(x):
    return float(np.sqrt(np.mean(np.square(x))))


def fit(x, y):
    matrix = design(x)
    inverse = np.linalg.pinv(matrix)
    coeff = inverse @ y
    residual = y - matrix @ coeff
    leverage = (matrix @ inverse).diagonal()
    return coeff, residual, leverage


def summary(x, y):
    coeff, residual, leverage = fit(x, y)
    return dict(coeff=coeff.tolist(), rms=rms(residual),
                r2=float(1 - residual.var()/y.var()),
                loo_rms=rms(residual/(1-leverage)))


def gradient_angle(coeff):
    return float(np.degrees(np.arctan2(coeff[1], coeff[2])) % 360)


def galactic_gradient(coeff):
    # IAU J2000 ICRS->Galactic rotation (Hipparcos catalogue, vol. 1, §1.5.3).
    rot = np.array([[-.0548755604, -.8734370902, -.4838350155],
                    [ .4941094279, -.4448296300,  .7469822445],
                    [-.8676661490, -.1980763734,  .4559837762]])
    a, d = np.radians([RA0, DEC0])
    unit = [np.cos(a)*np.cos(d), np.sin(a)*np.cos(d), np.sin(d)]
    ge = rot @ np.array([-np.sin(a), np.cos(a), 0])
    gn = rot @ np.array([-np.sin(d)*np.cos(a), -np.sin(d)*np.sin(a), np.cos(d)])
    gx, gy, gz = rot @ unit
    lon, lat = np.arctan2(gy, gx), np.arcsin(gz)
    gal_e = [-np.sin(lon), np.cos(lon), 0]
    gal_n = [-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)]
    tangent = coeff[1]*ge + coeff[2]*gn
    return [float(np.dot(tangent, gal_e)), float(np.dot(tangent, gal_n))]


def spatial_null(x, y, n=1000000):
    # Is the amount of variance captured appreciable for generic *correlated*
    # columns? The two-direction fitted statistic is used throughout. Simulations
    # condition on real sampling and do not tune the three structure exponents.
    rng = np.random.default_rng(65172026)
    centered = np.eye(len(y)) - 1/len(y)
    projection = design(x) @ np.linalg.pinv(design(x)) - np.ones((len(y), len(y)))/len(y)
    observed = summary(x, y)['r2']
    results = {}
    for kind in ['exchangeable', 1, 5/3, 1.95]:
        if kind != 'exchangeable':
            variogram = cdist(x, x)**kind
            eigen, vec = np.linalg.eigh(-0.5 * centered @ variogram @ centered)
            transform = vec * np.sqrt(np.maximum(eigen, 0))
        count = 0
        for first in range(0, n, 10000):
            batch = min(n-first,10000)
            if kind == 'exchangeable':
                # Permute entire pulsars; pairs of pulsars aren't independent.
                draws = rng.permuted(np.broadcast_to(y-y.mean(), (batch, len(y))), axis=1)
            else:
                draws = rng.normal(size=(batch, len(y))) @ transform.T
            r2 = ((draws @ projection) ** 2).sum(axis=1) / (draws ** 2).sum(axis=1)
            count += int(np.count_nonzero(r2 >= observed))
        results[str(kind)] = dict(exceedances=count, draws=n, p=(count+1)/(n+1))
    return results


def holdout(x, y, xy_gbt, dm_gbt):
    # GBT D is remote from the core group A-C in the original discovery paper.
    # Four nearly collinear GBT positions do not constrain a full 2D plane.
    axis = xy_gbt[3] - xy_gbt[:3].mean(axis=0)
    axis /= np.linalg.norm(axis)
    coeff = np.polyfit(xy_gbt @ axis, dm_gbt, 1)
    prediction = np.polyval(coeff, x @ axis)
    # An exact label-permutation test for the ten FAST-held-out DMs.
    # Scale ranks to integers and use the squared rank differences to avoid
    # floating point tie handling at the exact boundary.
    from itertools import islice, permutations
    rank_order = np.argsort(prediction)
    obs = stats.rankdata(y[rank_order]).astype(int)
    target = np.arange(1,len(y)+1)
    observed = abs(int(np.dot(2*target-(len(y)+1),2*obs-(len(y)+1))))
    null = permutations(target.tolist())
    extreme, total = 0, 0
    while batch := list(islice(null, 25000)):
        products = (2 * np.array(batch) - (len(y)+1)) @ (2*target-(len(y)+1))
        extreme += int(np.count_nonzero(np.abs(products) >= observed))
        total += len(batch)
    return dict(axis_e_n=axis.tolist(), slope=float(coeff[0]), intercept=float(coeff[1]),
                predicted=prediction.tolist(), measured=y.tolist(),
                pearson=float(stats.pearsonr(prediction, y).statistic),
                spearman=float(stats.spearmanr(prediction, y).statistic),
                exact_rank_perm=dict(extreme=extreme, permutations=total, p=extreme/total),
                rmse=rms(prediction-y), gbt_mean_rmse=rms(dm_gbt.mean()-y),
                gbt_core_mean_rmse=rms(dm_gbt[:3].mean()-y))


def analyze():
    fast, gbt = read('fast.csv'), read('green-bank.csv')
    x, xgbt = offsets(fast), offsets(gbt)
    names = [row['pulsar'][-1] for row in fast]
    y = np.array([float(row['dm_yin']) for row in fast])
    stack = np.array([float(row['dm_dai']) for row in fast])
    dm_gbt = np.array([float(row['dm_lynch']) for row in gbt])
    newer = np.array([p not in 'ABCD' for p in names])
    first = np.array([p in 'ABCDEFGHI' for p in names])
    coeff, residual, h = fit(x,y)
    inverse = np.linalg.pinv(design(x))
    hc3 = (inverse * (residual / (1-h))) @ (inverse * (residual / (1-h))).T
    se = np.sqrt(np.diag(hc3))
    early_coeff = fit(x[first], y[first])[0]
    early_pred = design(x[~first]) @ early_coeff
    baseline_r = np.linalg.norm(x,axis=1)
    jack = []
    for i in range(len(y)):
        keep = np.arange(len(y))!=i
        c = fit(x[keep], y[keep])[0]
        jack.append(dict(omit=names[i], position_angle=gradient_angle(c),
                         **summary(x[keep], y[keep])))
    proxy = np.array([-float(row['spin_frequency_derivative_hz_s'])/float(row['spin_frequency_hz'])
                      for row in fast])
    d_leaf = jack[names.index('D')]
    radial = summary(baseline_r[:,None],y)
    stack_fit=summary(x,stack)
    result = dict(n=len(y), origin=[RA0,DEC0], origin_radec=['18:01:50.52','-08:57:31.6'],
                  data_order=names, full_fit=summary(x,y), east_north_hc3_se=se[1:].tolist(),
                  pa_degrees=gradient_angle(coeff), magnitude=float(np.hypot(*coeff[1:])),
                  d_dm_per_arcsec_l_cosb_b=galactic_gradient(coeff),
                  constant_loo_rms=rms((y-y.mean())*len(y)/(len(y)-1)),
                  angular_radius_model=radial, east_only=summary(x[:,:1],y),
                  north_only=summary(x[:,1:],y), leave_one_out=jack,
                  omit_D_prediction=float((design(x[names.index('D'):names.index('D')+1]) @ d_leaf['coeff'])[0]),
                  positions_arcsec=x.tolist(), fitted_residuals=residual.tolist(),
                  nulls=spatial_null(x,y), historical_gbt=holdout(x[newer],y[newer],xgbt,dm_gbt),
                  first_nine_to_last_five=dict(predicted=early_pred.tolist(), measured=y[~first].tolist(),
                      rmse=rms(early_pred-y[~first]), mean_rmse=rms(y[first].mean()-y[~first])),
                  same_fast_stack=dict(fit=stack_fit, dm_rms_difference=rms(stack-y)),
                  residual_spin_fractional_derivative=dict(
                      spearman=float(stats.spearmanr(residual,proxy).statistic),
                      p=float(stats.spearmanr(residual,proxy).pvalue)))
    return result


if __name__ == '__main__':
    (ROOT / 'radio-analysis.json').write_text(json.dumps(analyze(),indent=2,allow_nan=False)+'\n')
