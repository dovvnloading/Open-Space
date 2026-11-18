"""Predict three other survey curves from already fitted KMT-2024-BLG-1516 shapes.

Each alternative uses the same epochs, ±32-day peak window, 0.18-day bins,
off-peak baseline constraint, formal uncertainties and independently fitted
positive source flux/zero-point. No held-out survey adjusts the lens geometry.
"""
import json
from pathlib import Path

import numpy as np

from long_fs import HERE, KMT, fit, jobs, read_curve, write_csv


def predict(job, geometry):
    event = job[0]
    t, f, e, g, diff = read_curve(job)
    t, f, e, g = fit.clean(t, f, e, g)
    if diff:
        e = fit.reported_error_floor(e, g)
    t = t-geometry['reference_jd']
    window = geometry['window_halfwidth_d']
    far = abs(t) > 3*window
    if far.sum()<4:
        far = abs(t)>window
    z = (f[far]-np.median(f[far]))/e[far]
    scale = max(1., np.median(abs(z))/.67448975)
    good = abs(z)<6*scale
    w_far = 1/(e[far][good]*scale)**2
    baseline = np.average(f[far][good], weights=w_far)
    within = abs(t)<=window
    _, group = np.unique(np.floor(t[within]/geometry['bin_width_d']), return_inverse=True)
    w = 1/(e[within]*scale)**2
    w_bin = np.bincount(group, weights=w)
    tb = np.bincount(group, weights=w*t[within])/w_bin
    fb = np.append(np.bincount(group, weights=w*f[within])/w_bin, baseline)
    w = np.append(w_bin, sum(w_far))
    cap = np.inf if diff else geometry['cap_factor']*np.median(f)
    vb = fit.engine()
    vb.Tol = .001
    score = {}
    fluxes = {}
    for model in ('ps', 'fs', 'bin', 'exp'):
        v = geometry['models'][model]['vector']
        dt, lte, lu = v[:3]
        te, u = np.exp([lte, lu])
        if model == 'ps':
            a = fit.lens.magnification(((tb-dt)/te)**2 + u*u)-1
        elif model == 'fs':
            a = np.array(vb.ESPLLightCurve([lu,lte,dt,v[3]],tb.tolist())[0])-1
        elif model == 'bin':
            a = np.array(vb.BinaryLightCurve([v[4],v[5],u,v[6],v[3],lte,dt],tb.tolist())[0])-1
        else:
            a = np.exp(np.clip(np.where(tb<dt,(tb-dt)/te,-(tb-dt)/u),-745,0))
        a = np.append(a,0.)
        a -= np.average(a,weights=w)
        centered = fb-np.average(fb,weights=w)
        source = np.clip(np.sum(w*a*centered)/np.sum(w*a*a),0.,np.inf if model=='exp' else cap)
        fluxes[model] = source
        score[model] = np.sum(w*(centered-source*a)**2)
    return dict(event=event,heldout_epochs=int(sum(within)),heldout_bins=len(tb),baseline_pseudopoints=1,
                error_inflation=round(scale,4),**{k+'_chi2':round(v,3) for k,v in score.items()},
                fs_minus_bin=round(score['fs']-score['bin'],3),ps_minus_fs=round(score['ps']-score['fs'],3),
                fs_source_counts_or_flux=round(fluxes['fs'],5),bin_source_counts_or_flux=round(fluxes['bin'],5))


if __name__ == '__main__':
    geometry = json.loads((HERE.parent/'peaks'/(KMT+'.json')).read_text())
    write_csv(HERE/'long_fs_heldout.csv',[predict(j,geometry) for j,_,_ in jobs()])
