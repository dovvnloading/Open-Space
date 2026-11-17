"""Diagnostic for a point source dominating the putative extended emission."""
import numpy as np
from scipy import ndimage


def extended_profile(image, valid, row):
    dx, dy = row["x1"]-row["x0"], row["y1"]-row["y0"]
    length = np.hypot(dx,dy)
    width = row["fwhm_px"]
    u = np.array([dx/length,dy/length])
    v = np.array([-u[1],u[0]])
    along = (np.arange(max(24, round(length)))+.5) / max(24,round(length))*length
    cross = width*np.array([-3.5,-3,-2.5, -.45,0,.45, 2.5,3,3.5])
    sx = row["x0"] + u[0]*along + v[0]*cross[:,None]
    sy = row["y0"] + u[1]*along + v[1]*cross[:,None]
    coverage = ndimage.map_coordinates(valid.astype("u1",copy=False), [sy,sx], order=0,
                                       mode="constant", cval=0)
    if np.min(coverage) == 0:
        return dict(profile_gate=False, profile_reason="mask or boundary")
    vals = ndimage.map_coordinates(image, [sy,sx], order=1, mode="constant", cval=0)
    low, center, high = [np.mean(z,axis=0) for z in np.split(vals, 3, axis=0)]
    contrasts = center - .5*(low+high)
    flux = np.array([z.sum() for z in np.array_split(contrasts,8)])
    pos = np.maximum(0,flux)
    dominance = float(pos.max()/pos.sum()) if pos.sum() > 0 else 1.
    both = [np.array_split(z,8) for z in [center,low,high]]
    n = sum(c.mean() > max(l.mean(),h.mean()) for c,l,h in zip(*both))
    return dict(profile_gate=bool(n >= 6 and dominance <= .40),
                center_above_both_flanks=int(n), max_of_eight_flux_fraction=dominance)
