"""Measure fixed-center, signed LSF fluxes in observed spectra and RVSpecFit models."""
import argparse
import csv
from pathlib import Path
import numpy as np
from desi_shard import targets

LI = 6709.7  # vacuum Angstrom, blended resonance doublet


def profile(e, center, half):
    obs = center * e.doppler
    j = int(np.searchsorted(e.wave, obs)) - 1
    m = len(e.resolution)
    radius = m // 2
    a, b = j-half, j+half+2
    if a < 1 or b >= len(e.wave) or half < radius:
        return None
    w = e.wave[a:b]
    d = np.gradient(e.wave)
    t = (obs-e.wave[j]) / (e.wave[j+1]-e.wave[j])
    p = np.zeros(len(w))
    for column, weight in [(j, 1-t), (j+1, t)]:
        k = column-radius
        r = np.asarray(e.resolution[:,column], dtype=float)
        norm = r @ d[k:k+m]
        if not norm > 0:
            return None
        p[k-a:k-a+m] += weight*r/norm
    weights = e.ivar[a:b]
    flux = e.flux[a:b]
    model = e.template[a:b]
    ok = (weights>0) & (e.mask[a:b]==0) & np.isfinite(weights+flux+model+p)
    if not ok[np.abs(w-obs)<2.5].all() or ok.sum() < len(w)//2:
        return None
    X = np.column_stack((p, np.ones(len(w)), w-obs))
    W = np.sqrt(weights[ok])
    Xw = X[ok]*W[:,None]
    xy = np.column_stack((flux[ok],model[ok]))*W[:,None]
    coef, _, rank, _ = np.linalg.lstsq(Xw, xy, rcond=None)
    if rank != 3 or np.any(coef[1] <= 0):
        return None
    sigma = np.sqrt(np.linalg.inv(Xw.T @ Xw)[0,0])
    return dict(observed_flux17=coef[0,0], template_flux17=coef[0,1], flux_sigma17=sigma,
                observed_cont17=coef[1,0], template_cont17=coef[1,1],
                observed_ew_A=coef[0,0]/(coef[1,0]*e.doppler),
                template_ew_A=coef[0,1]/(coef[1,1]*e.doppler),
                residual_snr=(coef[0,0]-coef[0,1])/sigma,
                observed_snr=coef[0,0]/sigma, n_good=int(ok.sum()))


def measure(cache, survey, program, hp, output):
    rows = []
    for meta, arms in targets(cache, survey, program, hp):
        for exp in arms["r"]:
            for center in [LI,LI-12,LI+12]:
                for half in (15,25):
                    measurement = profile(exp, center, half)
                    rows.append(dict(targetid=meta["targetid"], expid=exp.expid, rest_A=center,
                                     half_window_px=half, **(measurement or {})))
    path = output / f"lithium-{survey}-{program}-{hp}.csv"
    output.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    parser.add_argument("survey")
    parser.add_argument("program")
    parser.add_argument("healpix",type=int)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    measure(args.cache,args.survey,args.program,args.healpix,args.output)
